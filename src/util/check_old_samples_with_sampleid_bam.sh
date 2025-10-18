#!/usr/bin/env bash
set -euo pipefail

FASTQ_DIR="/home/ken/ken-nipt/fastq"
ANALYSIS_DIR="/home/ken/ken-nipt/analysis"
OUTPUT_TSV="result.tsv"
INPUT_FILE=""

usage() {
  cat <<'USAGE'
Usage: check_old_samples_with_sampleid_bam.sh -i <two-col.txt|-|omit> -o <output.tsv>

입력: <Order_id><TAB><Sample_id>
출력 헤더:
Order_id  Sample_id  Order_id_fastq  Order_id_analysis  Sample_id_fastq  Sample_id_analysis  sorted_bam  proper_paired_bam

판정 규칙:
- Order_id_*: /{fastq,analysis} 트리 depth=2(배치/샘플)에서 디렉터리명 "정확 일치"가 있으면 Y
- Sample_id_*: 해당 Order_id 디렉터리 내부(depth=1)에 'Sample_id*' 로 시작하는 "하위 디렉터리"가 있을 때만 Y (파일명 포함 매칭 X)
- sorted_bam / proper_paired_bam: 해당 Order_id 디렉터리 내부(depth<=2)에
  'Order_id.sorted.bam' / 'Order_id.proper_paired.bam' 파일이 있으면 Y
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

# ---- depth=2: Order_id 디렉터리 목록 & 경로 인덱싱 ----
# 키: Order_id  값: 경로들(::로 구분)
declare -A FASTQ_PATHS ANALYSIS_PATHS SET_FASTQ_DIR SET_ANALYSIS_DIR

# fastq
while IFS=$'\t' read -r order fullpath; do
  [[ -z "$order" ]] && continue
  SET_FASTQ_DIR["$order"]=1
  if [[ -n "${FASTQ_PATHS[$order]:-}" ]]; then
    FASTQ_PATHS["$order"]+="${FASTQ_PATHS[$order]:+::}$fullpath"
  else
    FASTQ_PATHS["$order"]="$fullpath"
  fi
done < <(find "$FASTQ_DIR" -mindepth 2 -maxdepth 2 -type d -printf '%f\t%p\n' 2>/dev/null)

# analysis
while IFS=$'\t' read -r order fullpath; do
  [[ -z "$order" ]] && continue
  SET_ANALYSIS_DIR["$order"]=1
  if [[ -n "${ANALYSIS_PATHS[$order]:-}" ]]; then
    ANALYSIS_PATHS["$order"]+="${ANALYSIS_PATHS[$order]:+::}$fullpath"
  else
    ANALYSIS_PATHS["$order"]="$fullpath"
  fi
done < <(find "$ANALYSIS_DIR" -mindepth 2 -maxdepth 2 -type d -printf '%f\t%p\n' 2>/dev/null)

# ---- 결과 헤더 ----
printf "Order_id\tSample_id\tOrder_id_fastq\tOrder_id_analysis\tSample_id_fastq\tSample_id_analysis\tsorted_bam\tproper_paired_bam\n" > "$OUTPUT_TSV"

# ---- 한 줄씩 처리 ----
"${INPUT_STREAM[@]}" | awk 'NF>=2 {print $1 "\t" $2}' | while IFS=$'\t' read -r ORDER_ID SAMPLE_ID; do
  ORDER_ID=$(echo -n "$ORDER_ID" | tr -d '\r' | xargs)
  SAMPLE_ID=$(echo -n "$SAMPLE_ID" | tr -d '\r' | xargs)
  [[ -z "$ORDER_ID" || -z "$SAMPLE_ID" ]] && continue

  o_fastq="N"; o_analysis="N"; s_fastq="N"; s_analysis="N"; has_sorted="N"; has_proper="N"

  # 1) Order_id 디렉터리 존재 여부
  [[ ${SET_FASTQ_DIR[$ORDER_ID]+_}    ]] && o_fastq="Y"
  [[ ${SET_ANALYSIS_DIR[$ORDER_ID]+_} ]] && o_analysis="Y"

  # 2) Sample_id_*: 해당 Order_id 디렉터리 내부에 'Sample_id*' 디렉터리 존재 확인 (파일은 무시)
  # fastq 쪽
  if [[ -n "${FASTQ_PATHS[$ORDER_ID]:-}" ]]; then
    IFS='::' read -r -a paths <<< "${FASTQ_PATHS[$ORDER_ID]}"
    for p in "${paths[@]}"; do
      if find "$p" -mindepth 1 -maxdepth 1 -type d -name "${SAMPLE_ID}*" -print -quit | grep -q .; then
        s_fastq="Y"; break
      fi
    done
  fi
  # analysis 쪽
  if [[ -n "${ANALYSIS_PATHS[$ORDER_ID]:-}" ]]; then
    IFS='::' read -r -a paths <<< "${ANALYSIS_PATHS[$ORDER_ID]}"
    for p in "${paths[@]}"; do
      if find "$p" -mindepth 1 -maxdepth 1 -type d -name "${SAMPLE_ID}*" -print -quit | grep -q .; then
        s_analysis="Y"; break
      fi
    done
  fi

  # 3) BAM 파일: 해당 Order_id 디렉터리 내부(depth<=2)에 파일 존재 확인
  if [[ -n "${ANALYSIS_PATHS[$ORDER_ID]:-}" ]]; then
    IFS='::' read -r -a paths <<< "${ANALYSIS_PATHS[$ORDER_ID]}"
    for p in "${paths[@]}"; do
      if find "$p" -mindepth 1 -maxdepth 2 -type f -name "${ORDER_ID}.sorted.bam" -print -quit | grep -q .; then
        has_sorted="Y"
      fi
      if find "$p" -mindepth 1 -maxdepth 2 -type f -name "${ORDER_ID}.proper_paired.bam" -print -quit | grep -q .; then
        has_proper="Y"
      fi
      [[ "$has_sorted$has_proper" == "YY" ]] && break
    done
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$ORDER_ID" "$SAMPLE_ID" "$o_fastq" "$o_analysis" "$s_fastq" "$s_analysis" "$has_sorted" "$has_proper" \
    >> "$OUTPUT_TSV"
done

echo "완료: $OUTPUT_TSV"
