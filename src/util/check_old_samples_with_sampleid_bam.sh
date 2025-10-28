#!/usr/bin/env bash
set -euo pipefail

FASTQ_DIR="/home/ken/ken-nipt/fastq"
ANALYSIS_DIR="/home/ken/ken-nipt/analysis"
OUTPUT_TSV="result.tsv"
INPUT_FILE=""

usage() {
  cat <<'USAGE'
Usage: check_order_sample_bam_depth2.sh -i <two-col.txt|-|omit> -o <output.tsv>

입력(two-col): <Order_id><TAB><Sample_id>
출력 헤더:
Order_id  Sample_id  Order_id_fastq  Order_id_analysis  Sample_id_fastq  Sample_id_analysis  sorted_bam  proper_paired_bam

규칙:
- fastq/analysis 아래의 work_dir(예: 2411, 250703_01 등)을 모두 훑고,
  그 바로 하위(depth=2)에 존재하는 디렉터리명으로 판단(정확 일치).
- sorted_bam / proper_paired_bam 은 analysis 쪽 해당 Order_id 디렉터리 내부(depth<=2)에서
  ORDER_ID.sorted.bam / ORDER_ID.proper_paired.bam 이 존재하면 Y.
USAGE
}

while getopts ":i:o:h" opt; do
  case "$opt" in
    i) INPUT_FILE="$OPTARG" ;;
    o) OUTPUT_TSV="$OPTARG" ;;
    h) usage; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 2 ;;
  esac
done
shift $((OPTIND-1))

# 입력 스트림
if [[ -n "${INPUT_FILE}" && "${INPUT_FILE}" != "-" ]]; then
  [[ -f "${INPUT_FILE}" ]] || { echo "입력 파일 없음: ${INPUT_FILE}" >&2; exit 1; }
  INPUT_STREAM=("cat" "${INPUT_FILE}")
else
  INPUT_STREAM=("cat" "-")
fi

# ---- depth=2의 모든 디렉터리명 세트 & 경로 매핑 ----
# 모든 depth=2 디렉터리명(샘플/오더 공통)
declare -A SET_FASTQ_D2 SET_ANALYSIS_D2
# 각 디렉터리명(보통 Order_id)에 해당하는 실제 경로들(::로 연결)
declare -A FASTQ_PATHS ANALYSIS_PATHS

# fastq: work_dir/* (depth=2) 디렉터리 수집
# -printf '%f\t%p\n' : basename<TAB>fullpath
while IFS=$'\t' read -r name full; do
  [[ -z "$name" ]] && continue
  SET_FASTQ_D2["$name"]=1
  if [[ -n "${FASTQ_PATHS[$name]:-}" ]]; then
    FASTQ_PATHS["$name"]+="${FASTQ_PATHS[$name]:+::}$full"
  else
    FASTQ_PATHS["$name"]="$full"
  fi
done < <(find "$FASTQ_DIR" -mindepth 2 -maxdepth 2 -type d -printf '%f\t%p\n' 2>/dev/null)

# analysis: work_dir/* (depth=2) 디렉터리 수집
while IFS=$'\t' read -r name full; do
  [[ -z "$name" ]] && continue
  SET_ANALYSIS_D2["$name"]=1
  if [[ -n "${ANALYSIS_PATHS[$name]:-}" ]]; then
    ANALYSIS_PATHS["$name"]+="${ANALYSIS_PATHS[$name]:+::}$full"
  else
    ANALYSIS_PATHS["$name"]="$full"
  fi
done < <(find "$ANALYSIS_DIR" -mindepth 2 -maxdepth 2 -type d -printf '%f\t%p\n' 2>/dev/null)

# ---- analysis 쪽 BAM 존재 여부(ORDER_ID 기준) 미리 계산 ----
declare -A HAS_SORTED_BAM HAS_PROPER_PAIRED_BAM
for id in "${!ANALYSIS_PATHS[@]}"; do
  IFS='::' read -r -a paths <<< "${ANALYSIS_PATHS[$id]}"
  for p in "${paths[@]}"; do
    [[ -z "$p" ]] && continue
    # depth<=2 (즉, 해당 디렉터리 바로 아래와 그 하위 한 단계)만 검색
    if find "$p" -mindepth 1 -maxdepth 2 -type f -name "${id}.sorted.bam" -print -quit | grep -q .; then
      HAS_SORTED_BAM["$id"]=1
    fi
    if find "$p" -mindepth 1 -maxdepth 2 -type f -name "${id}.proper_paired.bam" -print -quit | grep -q .; then
      HAS_PROPER_PAIRED_BAM["$id"]=1
    fi
    [[ ${HAS_SORTED_BAM[$id]+_} && ${HAS_PROPER_PAIRED_BAM[$id]+_} ]] && break
  done
done

# ---- 결과 파일 헤더 ----
printf "Order_id\tSample_id\tOrder_id_fastq\tOrder_id_analysis\tSample_id_fastq\tSample_id_analysis\tsorted_bam\tproper_paired_bam\n" > "$OUTPUT_TSV"

# ---- 입력 처리 ----
"${INPUT_STREAM[@]}" | awk 'NF>=2 {print $1 "\t" $2}' | \
while IFS=$'\t' read -r ORDER_ID SAMPLE_ID; do
  ORDER_ID=$(echo -n "$ORDER_ID" | tr -d '\r' | xargs)
  SAMPLE_ID=$(echo -n "$SAMPLE_ID" | tr -d '\r' | xargs)
  [[ -z "$ORDER_ID" || -z "$SAMPLE_ID" ]] && continue

  # 기본값
  o_fastq="N"; o_analysis="N"; s_fastq="N"; s_analysis="N"; has_sorted="N"; has_proper="N"

  # Order_id는 depth=2에서 정확 일치
  [[ ${SET_FASTQ_D2[$ORDER_ID]+_}    ]] && o_fastq="Y"
  [[ ${SET_ANALYSIS_D2[$ORDER_ID]+_} ]] && o_analysis="Y"

  # Sample_id도 depth=2에서 정확 일치 (샘플 중심 레이아웃 대응)
  [[ ${SET_FASTQ_D2[$SAMPLE_ID]+_}    ]] && s_fastq="Y"
  [[ ${SET_ANALYSIS_D2[$SAMPLE_ID]+_} ]] && s_analysis="Y"

  # BAM 파일(analysis, ORDER_ID 기준)
  [[ ${HAS_SORTED_BAM[$ORDER_ID]+_}        ]] && has_sorted="Y"
  [[ ${HAS_PROPER_PAIRED_BAM[$ORDER_ID]+_} ]] && has_proper="Y"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$ORDER_ID" "$SAMPLE_ID" "$o_fastq" "$o_analysis" "$s_fastq" "$s_analysis" "$has_sorted" "$has_proper" \
    >> "$OUTPUT_TSV"
done

echo "완료: $OUTPUT_TSV"
