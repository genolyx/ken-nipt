#!/usr/bin/env bash
# Summarize read counts of 4 BAMs (sorted/dedup/uniq/proper) per sample and compute rates.
# Targets: sample dirs whose name starts with GNCI* or OPC* under $ROOT/*/
# stdout: TSV (sample, sorted_read, dedup_read, uniq_read, proper_read, dedup_rate, uniq_rate, proper_rate, proper_rate_sorted)
# stderr: progress & warnings
# Usage:
#   chmod +x summarize_read_counts.sh
#   ./summarize_read_counts.sh > read_summary.tsv 2>read_summary.log
#
# Env:
#   ROOT=/home/ken/ken-nipt/analysis   (default)
#   THREADS=4                          (default)

set -euo pipefail
export LC_ALL=C

ROOT="${ROOT:-/home/ken/ken-nipt/analysis}"
THREADS="${THREADS:-4}"

shopt -s nullglob

# -------- collect sample dirs (GNCI* or OPC*) --------
declare -A DIRSET=()
for d in "$ROOT"/*/GNCI*/; do DIRSET["${d%/}"]=1; done
for d in "$ROOT"/*/OPC*/;  do DIRSET["${d%/}"]=1; done

SAMPLE_DIRS=("${!DIRSET[@]}")
IFS=$'\n' SAMPLE_DIRS=($(printf "%s\n" "${SAMPLE_DIRS[@]}" | sort))  # stable order
TOTAL=${#SAMPLE_DIRS[@]}
if (( TOTAL == 0 )); then
  echo "No GNCI* or OPC* sample directories found under $ROOT" >&2
  exit 1
fi

echo "Found ${TOTAL} sample directories under $ROOT (prefix GNCI* or OPC*)." >&2
printf "sample\tsorted_read\tdedup_read\tuniq_read\tproper_read\tdedup_rate\tuniq_rate\tproper_rate\tproper_rate_sorted\n"

i=0
for d in "${SAMPLE_DIRS[@]}"; do
  ((i++)) || true
  sample="$(basename "$d")"

  sorted_bam="$d/${sample}.sorted.bam"
  dedup_bam="$d/${sample}.dedup.bam"
  uniq_bam="$d/${sample}.uniq.bam"
  proper_bam="$d/${sample}.proper_paired.bam"

  # progress (compact, stderr)
  printf "[%d/%d] %s ...\n" "$i" "$TOTAL" "$sample" >&2

  # require all four BAMs; if any is missing, skip (log only)
  if [[ ! -s "$sorted_bam" || ! -s "$dedup_bam" || ! -s "$uniq_bam" || ! -s "$proper_bam" ]]; then
    echo "  ↳ skipped: missing BAM(s) for $sample" >&2
    [[ -s "$sorted_bam" ]] || echo "     - missing: $sorted_bam" >&2
    [[ -s "$dedup_bam"  ]] || echo "     - missing: $dedup_bam"  >&2
    [[ -s "$uniq_bam"   ]] || echo "     - missing: $uniq_bam"   >&2
    [[ -s "$proper_bam" ]] || echo "     - missing: $proper_bam" >&2
    continue
  fi

  # read counts (individual reads)
  sorted_read=$(samtools view -@ "$THREADS" -c "$sorted_bam" 2>/dev/null || echo 0)
  dedup_read=$(samtools view -@ "$THREADS" -c "$dedup_bam"  2>/dev/null || echo 0)
  uniq_read=$(samtools view -@ "$THREADS" -c "$uniq_bam"   2>/dev/null || echo 0)
  proper_read=$(samtools view -@ "$THREADS" -c "$proper_bam" 2>/dev/null || echo 0)

  # compute rates safely (NA if denominator is 0)
  awk -v s="$sample" -v a="$sorted_read" -v b="$dedup_read" -v c="$uniq_read" -v p="$proper_read" '
    function rate(num, den) { return (den==0)?"NA":num/den; }
    BEGIN{
      dr  = rate(b,a);  ur  = rate(c,b);
      pr  = rate(p,c);  prs = rate(p,a);
      if (dr  != "NA") dr  = sprintf("%.6f", dr);
      if (ur  != "NA") ur  = sprintf("%.6f", ur);
      if (pr  != "NA") pr  = sprintf("%.6f", pr);
      if (prs != "NA") prs = sprintf("%.6f", prs);
      printf("%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n", s, a, b, c, p, dr, ur, pr, prs);
    }'
done

echo "Done. Processed ${TOTAL} sample directories." >&2
