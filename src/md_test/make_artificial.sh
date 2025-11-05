#!/usr/bin/env bash
# Create one artificial cfDNA BAM by mixing a Mom sample and a Fetus sample,
# applying a heterozygous microdeletion (loss) to the fetal component only,
# then downsampling to the target paired-read count.
#
# Requirements: samtools, bedtools
#
# Inputs (as arguments):
#   --mom_bam   PATH  (e.g., Mom/GNCIxxxx.sorted.bam)
#   --fetus_bam PATH  (e.g., Female/GNCIxxxx.sorted.bam  or  Male/GNCIxxxx.sorted.bam)
#   --ff_map    PATH  (TSV: sample_id<TAB>ff_percent; sample_id is BAM basename without .sorted.bam)
#   --md_bed    PATH  (BED of target microdeletion region(s); 1-based half-open OK)
#   --ff_target FLOAT (e.g., 5 or 10 or 15)  # desired fetal fraction (%) in final mixture
#   --reads     INT   (target read count in final BAM, e.g., 10000000, 15000000, 20000000)
#                    Note: Uses MAPQ>=30 reads (eff_reads) for NIPT analysis
#   --outdir    PATH
#   [--gender   M|F]  (fetal gender; if not provided, will be detected from fetus BAM)
#   [--mom_idx  INT]  (mom sample index in TSV file, e.g., 1, 2, 3...)
#   [--fetus_idx INT] (fetus sample index in TSV file, e.g., 1, 2, 3...)
#   [--sample_id STR] (output sample ID; if not provided, will be generated from parameters)
#   [--seed     INT]  (default 42)
#   [--C_reads  INT]  (pre-merge total reads before final downsample; default = max(reads, 21000000))
#
# Notes:
# - We assume input BAMs are coordinate-sorted, duplicate-removed, and standard primary alignments.
# - Mixing is done by probabilistic downsampling on each source with -s <seed.frac>.
# - Fetal deletion: ~50% drop of fetal reads overlapping md region (heterozygous loss).
#   This is approximate at read-level; for exact pair-level dropping you'd need name-sorted BAM + pair-wise filter.
# - Output BAM is sorted + indexed.
#
set -euo pipefail

# ------------- parse args -------------
MOM_BAM=""; FETUS_BAM=""; FF_MAP=""; MD_BED=""; FF_TGT=""; READS=""; OUTDIR=""
GENDER=""; SAMPLE_ID=""; MOM_IDX=""; FETUS_IDX=""
SEED=42; C_READS=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mom_bam) MOM_BAM="$2"; shift 2;;
    --fetus_bam) FETUS_BAM="$2"; shift 2;;
    --ff_map) FF_MAP="$2"; shift 2;;
    --md_bed) MD_BED="$2"; shift 2;;
    --ff_target) FF_TGT="$2"; shift 2;;
    --reads) READS="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --gender) GENDER="$2"; shift 2;;
    --mom_idx) MOM_IDX="$2"; shift 2;;
    --fetus_idx) FETUS_IDX="$2"; shift 2;;
    --sample_id) SAMPLE_ID="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    --C_reads) C_READS="$2"; shift 2;;
    *) echo "Unknown arg: $1"; exit 1;;
  esac
done

for x in MOM_BAM FETUS_BAM FF_MAP MD_BED FF_TGT READS OUTDIR; do
  [[ -n "${!x}" ]] || { echo "Missing --${x,,}"; exit 1; }
done
[[ -s "$MOM_BAM" && -s "$FETUS_BAM" && -s "$FF_MAP" && -s "$MD_BED" ]] || { echo "Input file missing"; exit 1; }

mkdir -p "$OUTDIR"

# ------------- helpers -------------
bname() { 
  bn="$(basename "$1")"
  # Handle both .proper_paired.bam and .sorted.bam
  bn="${bn%%.proper_paired.bam}"
  bn="${bn%%.sorted.bam}"
  echo "$bn"
}
sam_count_reads() { # count reads (using properly paired primary alignments with MAPQ>=30)
  local bam="$1"
  samtools view -c -f 2 -F 256 -F 2048 -q 30 "$bam" 2>/dev/null || echo "0"
}
ff_lookup() { # from ff_map.tsv: sample_id \t ff_percent
  local sid="$1"
  awk -v s="$sid" -F'\t' '$1==s{print $2}' "$FF_MAP" | head -n1
}
detect_gender_from_bam() { # detect gender from BAM (simple Y chromosome check)
  local bam="$1"
  local y_reads=$(samtools view -c -f 2 -F 256 -F 2048 -q 30 "$bam" chrY 2>/dev/null || echo "0")
  local a_reads=$(samtools view -c -f 2 -F 256 -F 2048 -q 30 "$bam" chr1 2>/dev/null || echo "0")
  if [[ "$a_reads" -eq 0 ]]; then
    echo "F"
    return
  fi
  local ratio=$(awk -v y="$y_reads" -v a="$a_reads" 'BEGIN{printf "%.6f", y/a}')
  # Threshold: > 0.0001 for male
  if (( $(echo "$ratio > 0.0001" | bc -l) )); then
    echo "M"
  else
    echo "F"
  fi
}

# ------------- identify IDs & FFs -------------
MOM_ID="$(bname "$MOM_BAM")"
FETUS_ID="$(bname "$FETUS_BAM")"
FF_A="$(ff_lookup "$MOM_ID")"
FF_B="$(ff_lookup "$FETUS_ID")"

[[ -n "$FF_A" && -n "$FF_B" ]] || { echo "FF not found in ff_map.tsv for $MOM_ID or $FETUS_ID"; exit 1; }

# ------------- detect gender if not provided -------------
if [[ -z "$GENDER" ]]; then
  echo "[INFO] Detecting gender from fetus BAM..."
  GENDER=$(detect_gender_from_bam "$FETUS_BAM")
  echo "[INFO] Detected gender: $GENDER"
else
  GENDER=$(echo "$GENDER" | tr '[:lower:]' '[:upper:]')
  GENDER="${GENDER:0:1}"  # Take first character (M/F)
  echo "[INFO] Using provided gender: $GENDER"
fi

# ------------- generate sample_id if not provided -------------
if [[ -z "$SAMPLE_ID" ]]; then
  # Extract disease name from BED file (4th column) or use filename
  DISEASE_NAME=$(awk 'NR==1 {print $4}' "$MD_BED" 2>/dev/null || echo "")
  if [[ -z "$DISEASE_NAME" ]]; then
    DISEASE_NAME=$(basename "$MD_BED" .bed)
  fi
  # Clean disease name (remove spaces, special chars)
  DISEASE_NAME=$(echo "$DISEASE_NAME" | sed 's/[^a-zA-Z0-9_]//g' | tr '[:upper:]' '[:lower:]')
  
  # Calculate deletion size
  DEL_START=$(awk 'NR==1 {print $2}' "$MD_BED")
  DEL_END=$(awk 'NR==1 {print $3}' "$MD_BED")
  DEL_SIZE_MB=$(( (DEL_END - DEL_START) / 1000000 ))
  
  # Format: {mom_idx}_{fetus_idx}_{disease}_{ff}FF_{reads}M_{del}Mb_{gender}
  READS_M=$((READS / 1000000))
  
  # Include indices if provided
  if [[ -n "$MOM_IDX" && -n "$FETUS_IDX" ]]; then
    SAMPLE_ID="${MOM_IDX}_${FETUS_IDX}_${DISEASE_NAME}_FF${FF_TGT}_${READS_M}M_${DEL_SIZE_MB}Mb_${GENDER}"
  else
    # Fallback to old format if indices not provided
    SAMPLE_ID="${DISEASE_NAME}_FF${FF_TGT}_${READS_M}M_${DEL_SIZE_MB}Mb_${GENDER}"
  fi
  echo "[INFO] Generated sample_id: $SAMPLE_ID"
fi

# convert to fractions
fA=$(awk -v x="$FF_A" 'BEGIN{print x/100.0}')
fB=$(awk -v x="$FF_B" 'BEGIN{print x/100.0}')
fT=$(awk -v x="$FF_TGT" 'BEGIN{print x/100.0}')

# ------------- mixing weights -------------
# beta = (fT - fA) / (fB - fA), alpha = 1 - beta
beta=$(awk -v fA="$fA" -v fB="$fB" -v fT="$fT" 'BEGIN{
  if (fB==fA) {print -1; exit}
  print (fT - fA) / (fB - fA)
}')

# Clamp beta to [0.0, 1.0] if target FF is not achievable
beta_result=$(python3 <<PYEOF
import sys
b = float($beta)
if b < 0.0:
    print("0.0:1")
elif b > 1.0:
    print("1.0:2")
else:
    print(f"{b:.10f}:0")
PYEOF
)

# Parse result: format is "value:exit_code"
beta_clamped=$(echo "$beta_result" | cut -d: -f1)
beta_exit_code=$(echo "$beta_result" | cut -d: -f2)
beta="$beta_clamped"

# Warn if beta was clamped
if [[ "$beta_exit_code" == "1" ]]; then
  echo "[WARNING] Target FF=${FF_TGT}% is lower than achievable (min=${FF_A}%). Using beta=0.0 (pure mom)."
elif [[ "$beta_exit_code" == "2" ]]; then
  echo "[WARNING] Target FF=${FF_TGT}% is higher than achievable (max=${FF_B}%). Using beta=1.0 (pure fetus)."
fi

alpha=$(awk -v b="$beta" 'BEGIN{print 1.0-b}')

# ------------- choose pre-merge C (in reads) -------------
# Input READS is already in reads
TARGET_READS="$READS"
if [[ "$C_READS" -le 0 ]]; then
  # small overhead above target to avoid rounding loss
  if   [[ "$READS" -ge 40000000 ]]; then C_READS=41000000
  elif [[ "$READS" -ge 30000000 ]]; then C_READS=31000000
  else C_READS=21000000
  fi
fi

# compute read quotas for mom/fetus before deletion
A_reads=$(awk -v a="$alpha" -v C="$C_READS" 'BEGIN{printf "%.0f", a*C}')
B_reads=$(awk -v b="$beta"  -v C="$C_READS" 'BEGIN{printf "%.0f", b*C}')

echo "[INFO] Mom=$MOM_ID (FF_A=${FF_A}%), Fetus=$FETUS_ID (FF_B=${FF_B}%), target FF=${FF_TGT}%"
echo "[INFO] alpha=${alpha}, beta=${beta}, target reads=${TARGET_READS}, C_reads=${C_READS}"
echo "[INFO] quotas: A=${A_reads} reads ($((A_reads/2)) pairs), B=${B_reads} reads ($((B_reads/2)) pairs)"

# ------------- compute subsample fractions relative to current sizes -------------
A_total=$(sam_count_reads "$MOM_BAM")
B_total=$(sam_count_reads "$FETUS_BAM")
[[ "$A_total" -gt 0 && "$B_total" -gt 0 ]] || { echo "Zero reads in input"; exit 1; }

pA=$(awk -v need="$A_reads" -v tot="$A_total" 'BEGIN{printf "%.6f", (tot>0 ? need/tot : 0)}')
pB=$(awk -v need="$B_reads" -v tot="$B_total" 'BEGIN{printf "%.6f", (tot>0 ? need/tot : 0)}')
# clamp 0..1
for v in pA pB; do
  eval "x=\$$v"; eval "$v=$(awk -v x="$x" 'BEGIN{if(x<0) x=0; if(x>1) x=1; print x}')"
done

# ------------- temp files -------------
TMPD="$(mktemp -d)"
trap 'rm -rf "$TMPD"' EXIT

A_sub="$TMPD/A.sub.bam"
B_sub="$TMPD/B.sub.bam"
B_in="$TMPD/B.in_md.bam"
B_out="$TMPD/B.out_md.bam"
B_in_half="$TMPD/B.in_md.half.bam"
B_del="$TMPD/B.with_del.bam"
MIX_BAM="$TMPD/mix.pre.bam"

# Output BAM name based on sample_id
OUT_BAM="$OUTDIR/${SAMPLE_ID}.proper_paired.bam"

# ------------- performance settings -------------
THREADS="${THREADS:-8}"  # Use more threads for faster processing
COMPRESS_LEVEL="${COMPRESS_LEVEL:-6}"  # Compression level (1-9, higher = better compression but slower)

# ------------- subsample each source by quota -------------
# use samtools -s SEED.frac ; mates keep/skip consistently based on QNAME hash
samtools view -b -@ "$THREADS" -s "${SEED}.$(printf "%06d" "$(awk -v p="$pA" 'BEGIN{print int(p*1000000)}')")" "$MOM_BAM" -o "$A_sub"
# Skip indexing for intermediate files - not needed for merging
samtools view -b -@ "$THREADS" -s "$((SEED+1)).$(printf "%06d" "$(awk -v p="$pB" 'BEGIN{print int(p*1000000)}')")" "$FETUS_BAM" -o "$B_sub"

# ------------- apply fetal microdeletion (heterozygous loss → ~50% drop in-region) -------------
# Split B into reads overlapping md8.bed vs outside; then random-half keep on the in-region set; merge back.
# (Approximate at read-level; acceptable for shallow WGS)

# Check BAM chromosome format and convert BED if needed
FIRST_CHR=$(samtools view -H "$B_sub" | grep "^@SQ" | head -n1 | awk '{print $2}' | cut -d: -f2)
MD_BED_CONVERTED="$TMPD/md_bed_converted.bed"

if [[ "$FIRST_CHR" =~ ^chr ]]; then
  # BAM uses chr prefix, ensure BED does too
  awk -v OFS='\t' '{
    chr=$1
    if (chr !~ /^chr/) {
      if (chr == "M") chr="chrM"
      else chr="chr"chr
    }
    print chr, $2, $3, $4, $5, $6, $7, $8, $9, $10
  }' "$MD_BED" > "$MD_BED_CONVERTED"
else
  # BAM doesn't use chr prefix, remove from BED if present
  awk -v OFS='\t' '{
    chr=$1
    sub(/^chr/, "", chr)
    print chr, $2, $3, $4, $5, $6, $7, $8, $9, $10
  }' "$MD_BED" > "$MD_BED_CONVERTED"
fi

bedtools intersect -abam "$B_sub" -b "$MD_BED_CONVERTED"            > "$B_in"
bedtools intersect -abam "$B_sub" -b "$MD_BED_CONVERTED" -v         > "$B_out"
samtools view -b -@ "$THREADS" -s "$((SEED+2)).500000" "$B_in"          > "$B_in_half"
samtools merge -@ "$THREADS" -f "$B_del" "$B_out" "$B_in_half"
samtools sort -@ "$THREADS" -l "$COMPRESS_LEVEL" -o "$B_del" "$B_del"

# ------------- merge A + B -------------
samtools merge -@ "$THREADS" -f "$MIX_BAM" "$A_sub" "$B_del"
samtools sort -@ "$THREADS" -l "$COMPRESS_LEVEL" -o "$MIX_BAM" "$MIX_BAM"

# ------------- final downsample to exact target reads -------------
CUR=$(sam_count_reads "$MIX_BAM")
frac=$(awk -v want="$TARGET_READS" -v cur="$CUR" 'BEGIN{printf "%.6f", (cur>0 ? want/cur : 0)}')
# Downsample and sort with high compression level
samtools view -b -@ "$THREADS" -s "$((SEED+3)).$(printf "%06d" "$(awk -v p="$frac" 'BEGIN{print int(p*1000000)}')")" "$MIX_BAM" | \
  samtools sort -@ "$THREADS" -l "$COMPRESS_LEVEL" -o "$OUT_BAM" -
samtools index -@ "$THREADS" "$OUT_BAM"

# ------------- deletion region validation -------------
# Extract region from BED file and add "chr" prefix if needed
CHR=$(awk 'NR==1 {print $1}' "$MD_BED")
START=$(awk 'NR==1 {print $2}' "$MD_BED")
END=$(awk 'NR==1 {print $3}' "$MD_BED")

# Use the same chromosome format as BAM (should match MD_BED_CONVERTED)
FIRST_CHR=$(samtools view -H "$MOM_BAM" | grep "^@SQ" | head -n1 | awk '{print $2}' | cut -d: -f2)
if [[ "$FIRST_CHR" =~ ^chr ]]; then
  # BAM uses chr prefix
  if [[ ! "$CHR" =~ ^chr ]]; then
    if [[ "$CHR" == "M" ]]; then
      CHR="chrM"
    else
      CHR="chr${CHR}"
    fi
  fi
else
  # BAM doesn't use chr prefix, remove if present
  CHR="${CHR#chr}"
fi

# Define deletion region and flanking regions (same size as deletion)
REGION_SIZE=$((END - START))
DEL_REGION="${CHR}:${START}-${END}"

# Upstream flanking region (before deletion)
UPSTREAM_START=$(awk -v s="$START" -v r="$REGION_SIZE" 'BEGIN{x=s-r; if(x<1) x=1; printf "%.0f", x}')
UPSTREAM_END=$((START - 1))
if [[ "$UPSTREAM_END" -lt "$UPSTREAM_START" ]]; then
  UPSTREAM_END=$((UPSTREAM_START + 1))
fi
UPSTREAM_REGION="${CHR}:${UPSTREAM_START}-${UPSTREAM_END}"

# Downstream flanking region (after deletion)
DOWNSTREAM_START=$((END + 1))
DOWNSTREAM_END=$((END + REGION_SIZE))
DOWNSTREAM_REGION="${CHR}:${DOWNSTREAM_START}-${DOWNSTREAM_END}"

# Count reads in all three regions for both Mom and Output
# Deletion region
MOM_DEL_READS=$(samtools view -c "$MOM_BAM" "$DEL_REGION" 2>/dev/null || echo "0")
OUT_DEL_READS=$(samtools view -c "$OUT_BAM" "$DEL_REGION" 2>/dev/null || echo "0")

# Upstream flanking
MOM_UP_READS=$(samtools view -c "$MOM_BAM" "$UPSTREAM_REGION" 2>/dev/null || echo "0")
OUT_UP_READS=$(samtools view -c "$OUT_BAM" "$UPSTREAM_REGION" 2>/dev/null || echo "0")

# Downstream flanking
MOM_DOWN_READS=$(samtools view -c "$MOM_BAM" "$DOWNSTREAM_REGION" 2>/dev/null || echo "0")
OUT_DOWN_READS=$(samtools view -c "$OUT_BAM" "$DOWNSTREAM_REGION" 2>/dev/null || echo "0")

# Calculate ratios (Output/Mom for each region)
if [[ "$MOM_DEL_READS" -gt 0 ]]; then
  DEL_RATIO=$(awk -v out="$OUT_DEL_READS" -v mom="$MOM_DEL_READS" 'BEGIN{printf "%.4f", out/mom}')
else
  DEL_RATIO="N/A"
fi

if [[ "$MOM_UP_READS" -gt 0 ]]; then
  UP_RATIO=$(awk -v out="$OUT_UP_READS" -v mom="$MOM_UP_READS" 'BEGIN{printf "%.4f", out/mom}')
else
  UP_RATIO="N/A"
fi

if [[ "$MOM_DOWN_READS" -gt 0 ]]; then
  DOWN_RATIO=$(awk -v out="$OUT_DOWN_READS" -v mom="$MOM_DOWN_READS" 'BEGIN{printf "%.4f", out/mom}')
else
  DOWN_RATIO="N/A"
fi

# ------------- report -------------
FINAL_READS=$(sam_count_reads "$OUT_BAM")
FINAL_PAIRS=$((FINAL_READS / 2))
echo "[DONE] Output: $OUT_BAM (reads=${FINAL_READS}, pairs=${FINAL_PAIRS})"
echo "[DELETION CHECK]"
echo "  Upstream   (${UPSTREAM_REGION}): Mom=$MOM_UP_READS, Output=$OUT_UP_READS, Ratio=$UP_RATIO"
echo "  Deletion   (${DEL_REGION}): Mom=$MOM_DEL_READS, Output=$OUT_DEL_READS, Ratio=$DEL_RATIO"
echo "  Downstream (${DOWNSTREAM_REGION}): Mom=$MOM_DOWN_READS, Output=$OUT_DOWN_READS, Ratio=$DOWN_RATIO"

# ------------- create JSON metadata file -------------
JSON_FILE="$OUTDIR/${SAMPLE_ID}.json"
DEL_START=$(awk 'NR==1 {print $2}' "$MD_BED")
DEL_END=$(awk 'NR==1 {print $3}' "$MD_BED")
DEL_CHR=$(awk 'NR==1 {print $1}' "$MD_BED")
# Remove chr prefix if present for consistency
DEL_CHR="${DEL_CHR#chr}"
DEL_SIZE=$((DEL_END - DEL_START))
DEL_SIZE_MB=$((DEL_SIZE / 1000000))
DISEASE_NAME=$(awk 'NR==1 {print $4}' "$MD_BED" 2>/dev/null || echo "")

# Create JSON with metadata
python3 <<PYEOF > "$JSON_FILE"
import json
import datetime

metadata = {
    "sample_id": "$SAMPLE_ID",
    "creation_date": datetime.datetime.now().isoformat(),
    "pipeline": "make_artificial.sh",
    "source_samples": {
        "mom": {
            "sample_id": "$MOM_ID",
            "bam_path": "$MOM_BAM",
            "ff_percent": float("$FF_A"),
            "index": int("$MOM_IDX") if "$MOM_IDX" else None
        },
        "fetus": {
            "sample_id": "$FETUS_ID",
            "bam_path": "$FETUS_BAM",
            "ff_percent": float("$FF_B"),
            "index": int("$FETUS_IDX") if "$FETUS_IDX" else None
        }
    },
    "target_parameters": {
        "ff_target_percent": float("$FF_TGT"),
        "target_reads": int("$READS"),
        "target_pairs": int("$READS") // 2
    },
    "actual_output": {
        "reads": int("$FINAL_READS"),
        "pairs": int("$FINAL_PAIRS"),
        "bam_path": "$OUT_BAM"
    },
    "deletion": {
        "chromosome": "$DEL_CHR",
        "start": int("$DEL_START"),
        "end": int("$DEL_END"),
        "size_bp": int("$DEL_SIZE"),
        "size_mb": int("$DEL_SIZE_MB"),
        "bed_file": "$MD_BED",
        "disease_name": "$DISEASE_NAME"
    },
    "gender": "$GENDER",
    "mixing_parameters": {
        "alpha": float("$alpha"),
        "beta": float("$beta"),
        "c_reads": int("$C_READS")
    },
    "deletion_validation": {
        "upstream_ratio": "$UP_RATIO" if "$UP_RATIO" != "N/A" else None,
        "deletion_ratio": "$DEL_RATIO" if "$DEL_RATIO" != "N/A" else None,
        "downstream_ratio": "$DOWN_RATIO" if "$DOWN_RATIO" != "N/A" else None
    },
    "calculated_ff": {
        "yff": None,
        "seqff": None,
        "fragment_ff": None,
        "final_ff": None,
        "final_method": None
    }
}

with open("$JSON_FILE", "w") as f:
    json.dump(metadata, f, indent=2)

PYEOF

echo "[INFO] Metadata saved to: $JSON_FILE"
