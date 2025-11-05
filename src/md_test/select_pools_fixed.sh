#!/usr/bin/env bash

# Select mother (A-pool) and fetus (B-pool) samples from candidate TSV.

# Input TSV columns (from select_mother_fetus_candidates.sh):
#   batch  sample  gender  ff_type  ff_value  eff_reads  eff_pairs  ok_10M  ok_15M  ok_20M  bam_path
#
# Default criteria (override by env vars):
#   A(모체) 풀: FF <= 5.0 (percent), gender=XX 우선, eff_pairs >= 15,000,000
#   B(태아) 풀: FF >= 17.0 (percent), gender mix (XX/XY 모두 허용), eff_pairs >= 15,000,000
#   목표 개수: A_N=10, B_N=12  (부하 낮추려면 줄이세요)
#   배치 편중 제한: BATCH_MAX=999 (제한 끄기). 예: 4로 설정하면 배치별 최대 4개.
#
# Usage:
#   ./select_pools.sh sample_candidates.tsv
# Output:
#   mother_pool.tsv
#   fetus_pool.tsv
#
# Notes:
# - ff_value 가 퍼센트 기호 없이 실수라고 가정(예: 15.21). 퍼센트 기호가 있어도 제거 처리함.
# - eff_pairs는 properly paired 기준 "pair 수".
# - 진행/요약은 stderr로, 결과 TSV는 stdout이 아닌 개별 파일로 저장.

set -euo pipefail
export LC_ALL=C

IN="${1:-}"
if [[ -z "${IN}" || ! -s "${IN}" ]]; then
  echo "Usage: $0 sample_candidates.tsv" >&2
  exit 1
fi

# ---- configurable thresholds ----
A_FF_MAX="${A_FF_MAX:-5.0}"         # mother FF <= 5%
B_FF_MIN="${B_FF_MIN:-17.0}"        # fetus FF >= 17%
A_MIN_PAIRS="${A_MIN_PAIRS:-15000000}"
B_MIN_PAIRS="${B_MIN_PAIRS:-15000000}"
A_N="${A_N:-10}"
B_N="${B_N:-12}"
BATCH_MAX="${BATCH_MAX:-999}"       # per-batch cap (999 = no cap)

echo "Criteria:" >&2
echo "  A(mother): FF <= ${A_FF_MAX} %, eff_pairs >= ${A_MIN_PAIRS}, prefer gender=XX" >&2
echo "  B(fetus) : FF >= ${B_FF_MIN} %, eff_pairs >= ${B_MIN_PAIRS}, gender=XX/XY" >&2
echo "  Targets  : A_N=${A_N}, B_N=${B_N}, BATCH_MAX=${BATCH_MAX}" >&2

# strip % if present, and guard non-numeric
cleanup_tsv="$(mktemp)"
awk -F'\t' 'BEGIN{OFS="\t"}
NR==1{print; next}
{
  g=$3; ff=$5; gsub(/%/,"",ff);
  if(ff=="" || ff=="NA"){ff="NA"}
  print $1,$2,g,$4,ff,$6,$7,$8,$9,$10,$11
}' "$IN" > "$cleanup_tsv"

# headers
echo -e "sample\tbatch\tgender\tff\tff_type\teff_pairs\tbam_path\trank" > mother_pool.tsv
echo -e "sample\tbatch\tgender\tff\tff_type\teff_pairs\tbam_path\trank" > fetus_pool.tsv

# ---- select A-pool (mother) ----
# rule: gender=XX first, FF<=A_FF_MAX, eff_pairs>=A_MIN_PAIRS
# sort key: (1) gender priority XX>XY, (2) FF asc (낮을수록 좋음), (3) eff_pairs desc
# batch cap applied

echo "Selecting A(mother) pool ..." >&2
awk -v OFS="\t" -v amax="$A_FF_MAX" -v minp="$A_MIN_PAIRS" '
NR==1{next}
{
  batch=$1; sample=$2; gender=$3; fftype=$4; ff=$5; eff_pairs=$7; bam=$11;
  if(ff=="NA" || eff_pairs=="" || eff_pairs=="NA"){next}
  if(ff+0 <= amax+0 && eff_pairs+0 >= minp+0){
    prio=(gender=="XX"?0:1); # XX 우선
    printf "%s\t%s\t%s\t%f\t%s\t%d\t%s\t%d\n",
      sample,batch,gender,ff,fftype,eff_pairs,bam,prio
  }
}' "$cleanup_tsv" \
| sort -t $'\t' -k8,8n -k4,4n -k6,6nr \
> .A_candidates.tsv

# batch cap & top N
awk -v OFS="\t" -v N="$A_N" -v CAP="$BATCH_MAX" '
BEGIN{print_count=0}
{
  sample=$1; batch=$2;
  if(cap[batch] >= CAP) next;
  print; cap[batch]++; print_count++;
  if(print_count>=N) exit
}' .A_candidates.tsv \
| awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,"A"NR}' \
>> mother_pool.tsv

# ---- select B-pool (fetus) ----
# rule: FF>=B_FF_MIN, eff_pairs>=B_MIN_PAIRS, gender mix (both XX/XY allowed)
# sort key: (1) FF desc (높을수록 좋음), (2) eff_pairs desc

echo "Selecting B(fetus) pool ..." >&2
awk -v OFS="\t" -v bmin="$B_FF_MIN" -v minp="$B_MIN_PAIRS" '
NR==1{next}
{
  batch=$1; sample=$2; gender=$3; fftype=$4; ff=$5; eff_pairs=$7; bam=$11;
  if(ff=="NA" || eff_pairs=="" || eff_pairs=="NA"){next}
  if(ff+0 >= bmin+0 && eff_pairs+0 >= minp+0){
    printf "%s\t%s\t%s\t%f\t%s\t%d\t%s\n",
      sample,batch,gender,ff,fftype,eff_pairs,bam
  }
}' "$cleanup_tsv" \
| sort -t $'\t' -k4,4nr -k6,6nr \
> .B_candidates.tsv

# batch cap & top N
awk -v OFS="\t" -v N="$B_N" -v CAP="$BATCH_MAX" '
BEGIN{print_count=0}
{
  sample=$1; batch=$2;
  if(cap[batch] >= CAP) next;
  print; cap[batch]++; print_count++;
  if(print_count>=N) exit
}' .B_candidates.tsv \
| awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,"B"NR}' \
>> fetus_pool.tsv

rm -f "$cleanup_tsv" .A_candidates.tsv .B_candidates.tsv

# summary
A_CNT=$(( $(wc -l < mother_pool.tsv) - 1 ))
B_CNT=$(( $(wc -l < fetus_pool.tsv) - 1 ))
echo "Done. Selected A(mother)=${A_CNT}, B(fetus)=${B_CNT}" >&2
echo "  -> mother_pool.tsv / fetus_pool.tsv" >&2

