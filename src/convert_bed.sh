#!/usr/bin/env bash
# Usage: bash convert_bed.sh input.bed > output.bed
# Input : chr<TAB>start<TAB>end<TAB>disease_name<TAB>gain_loss<TAB>summary
# Output: chrom(without 'chr')<TAB>start<TAB>end<TAB>disease_name<TAB>gain_loss<TAB>overlap<TAB>-<TAB>-<TAB>-<TAB>-

set -euo pipefail
in="${1:-/dev/stdin}"

awk -v FS='\t' -v OFS='\t' '
# skip empty/comment
/^$/ || /^#/ { next }

{
  # 1) chr 컬럼에서 "chr" 제거 (chr1→1, chrX→X, chrY→Y, chrM→M)
  chr = $1
  sub(/^chr/, "", chr)

  start   = $2
  end     = $3
  disease = $4
  gl      = $5    # gain/loss 그대로 사용

  # summary(6번째)는 버림, overlap 추가, '-' 4개 추가
  print chr, start, end, disease, gl, "overlap", "-", "-", "-", "-"
}' "$in"
