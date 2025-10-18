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
#   --pairs     INT   (target paired-read count in final BAM, e.g., 10000000, 15000000, 20000000)
#   --outdir    PATH
#   [--seed     INT]  (default 42)
#   [--C_pairs  INT]  (pre-merge total pairs before final downsample; default = max(pairs, 20500000))
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
MOM_BAM=""; FETUS_BAM=""; FF_MAP=""; MD_BED=""; FF_TGT=""; PAIRS=""; OUTDIR=""
SEED=42; C_PAIRS=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mom_bam) MOM_BAM="$2"; shift 2;;
    --fetus_bam) FETUS_BAM="$2"; shift 2;;
    --ff_map) FF_MAP="$2"; shift 2;;
    --md_bed) MD_BED="$2"; shift 2;;
    --ff_target) FF_TGT="$2"; shift 2;;
    --pairs) PAIRS="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    --C_pairs) C_PAIRS="$2"; shift 2;;
    *) echo "Unknown arg: $1"; exit 1;;
  esac
done

for x in MOM_BAM FETUS_BAM FF_MAP MD_BED FF_TGT PAIRS OUTDIR; do
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
sam_count_pairs() { # count paired reads as pairs (using properly paired primary alignments; /2 since -c counts reads)
  local bam="$1"; local reads
  reads=$(samtools view -c -f 2 -F 256 -F 2048 "$bam")
  echo $(( reads / 2 ))
}
ff_lookup() { # from ff_map.tsv: sample_id \t ff_percent
  local sid="$1"
  awk -v s="$sid" -F'\t' '$1==s{print $2}' "$FF_MAP" | head -n1
}

# ------------- identify IDs & FFs -------------
MOM_ID="$(bname "$MOM_BAM")"
FETUS_ID="$(bname "$FETUS_BAM")"
FF_A="$(ff_lookup "$MOM_ID")"
FF_B="$(ff_lookup "$FETUS_ID")"

[[ -n "$FF_A" && -n "$FF_B" ]] || { echo "FF not found in ff_map.tsv for $MOM_ID or $FETUS_ID"; exit 1; }

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
# sanity
python3 - "$beta" <<'PY' || { echo "Invalid beta"; exit 1; }
import sys
b=float(sys.argv[1])
sys.exit(0 if 0.0 <= b <= 1.0 else 1)
PY
alpha=$(awk -v b="$beta" 'BEGIN{print 1.0-b}')

# ------------- choose pre-merge C -------------
if [[ "$C_PAIRS" -le 0 ]]; then
  # small overhead above target to avoid rounding loss
  if   [[ "$PAIRS" -ge 20000000 ]]; then C_PAIRS=20500000
  elif [[ "$PAIRS" -ge 15000000 ]]; then C_PAIRS=15500000
  else C_PAIRS=10500000
  fi
fi

# compute pair quotas for mom/fetus before deletion
A_pairs=$(awk -v a="$alpha" -v C="$C_PAIRS" 'BEGIN{printf "%.0f", a*C}')
B_pairs=$(awk -v b="$beta"  -v C="$C_PAIRS" 'BEGIN{printf "%.0f", b*C}')

echo "[INFO] Mom=$MOM_ID (FF_A=${FF_A}%), Fetus=$FETUS_ID (FF_B=${FF_B}%), target FF=${FF_TGT}%"
echo "[INFO] alpha=${alpha}, beta=${beta}, C_pairs=${C_PAIRS}, quotas: A=${A_pairs}, B=${B_pairs}"

# ------------- compute subsample fractions relative to current sizes -------------
A_total=$(sam_count_pairs "$MOM_BAM")
B_total=$(sam_count_pairs "$FETUS_BAM")
[[ "$A_total" -gt 0 && "$B_total" -gt 0 ]] || { echo "Zero pairs in input"; exit 1; }

pA=$(awk -v need="$A_pairs" -v tot="$A_total" 'BEGIN{printf "%.6f", (tot>0 ? need/tot : 0)}')
pB=$(awk -v need="$B_pairs" -v tot="$B_total" 'BEGIN{printf "%.6f", (tot>0 ? need/tot : 0)}')
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

# Extract disease label from BED filename for output naming
MD_LABEL=$(basename "$MD_BED" .bed)
# Simplified output name (no sample IDs, just sequential)
OUT_BAM="$OUTDIR/output.bam"  # Will be renamed by caller with proper index

# ------------- subsample each source by quota -------------
# use samtools -s SEED.frac ; mates keep/skip consistently based on QNAME hash
samtools view -b -s "${SEED}.$(printf "%06d" "$(awk -v p="$pA" 'BEGIN{print int(p*1000000)}')")" "$MOM_BAM" -o "$A_sub"
samtools index -@ 2 "$A_sub"
samtools view -b -s "$((SEED+1)).$(printf "%06d" "$(awk -v p="$pB" 'BEGIN{print int(p*1000000)}')")" "$FETUS_BAM" -o "$B_sub"
samtools index -@ 2 "$B_sub"

# ------------- apply fetal microdeletion (heterozygous loss → ~50% drop in-region) -------------
# Split B into reads overlapping md8.bed vs outside; then random-half keep on the in-region set; merge back.
# (Approximate at read-level; acceptable for shallow WGS)
bedtools intersect -abam "$B_sub" -b "$MD_BED"            > "$B_in"
bedtools intersect -abam "$B_sub" -b "$MD_BED" -v         > "$B_out"
samtools view -b -s "$((SEED+2)).500000" "$B_in"          > "$B_in_half"
samtools merge -@ 2 -f "$B_del" "$B_out" "$B_in_half"
samtools sort -@ 2 -o "$B_del" "$B_del"
samtools index -@ 2 "$B_del"

# ------------- merge A + B -------------
samtools merge -@ 2 -f "$MIX_BAM" "$A_sub" "$B_del"
samtools sort -@ 4 -o "$MIX_BAM" "$MIX_BAM"
samtools index -@ 2 "$MIX_BAM"

# ------------- final downsample to exact target pairs -------------
CUR=$(sam_count_pairs "$MIX_BAM")
frac=$(awk -v want="$PAIRS" -v cur="$CUR" 'BEGIN{printf "%.6f", (cur>0 ? want/cur : 0)}')
samtools view -b -s "$((SEED+3)).$(printf "%06d" "$(awk -v p="$frac" 'BEGIN{print int(p*1000000)}')")" "$MIX_BAM" -o "$OUT_BAM"
samtools sort -@ 4 -o "$OUT_BAM" "$OUT_BAM"
samtools index -@ 2 "$OUT_BAM"

# ------------- deletion region validation -------------
# Extract region from BED file and add "chr" prefix if needed
CHR=$(awk 'NR==1 {print $1}' "$MD_BED")
START=$(awk 'NR==1 {print $2}' "$MD_BED")
END=$(awk 'NR==1 {print $3}' "$MD_BED")

# Check if BAM uses "chr" prefix
FIRST_CHR=$(samtools view -H "$MOM_BAM" | grep "^@SQ" | head -n1 | awk '{print $2}' | cut -d: -f2)
if [[ "$FIRST_CHR" =~ ^chr ]]; then
  # Add chr prefix if not present
  if [[ ! "$CHR" =~ ^chr ]]; then
    CHR="chr${CHR}"
  fi
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
FINAL=$(sam_count_pairs "$OUT_BAM")
echo "[DONE] Output: $OUT_BAM (pairs=${FINAL})"
echo "[DELETION CHECK]"
echo "  Upstream   (${UPSTREAM_REGION}): Mom=$MOM_UP_READS, Output=$OUT_UP_READS, Ratio=$UP_RATIO"
echo "  Deletion   (${DEL_REGION}): Mom=$MOM_DEL_READS, Output=$OUT_DEL_READS, Ratio=$DEL_RATIO"
echo "  Downstream (${DOWNSTREAM_REGION}): Mom=$MOM_DOWN_READS, Output=$OUT_DOWN_READS, Ratio=$DOWN_RATIO"
