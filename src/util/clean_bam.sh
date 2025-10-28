#!/usr/bin/env bash
set -euo pipefail

# 사용법:
#   ./clean_bams.sh [BASE_DIR] [--apply]
# 예:
#   ./clean_bams.sh                          # 드라이런 (기본: ~/ken-nipt/analysis)
#   ./clean_bams.sh --apply                  # 실제 삭제
#   ./clean_bams.sh /data/ken-nipt/analysis  # 드라이런 (커스텀 경로)
#   ./clean_bams.sh /data/ken-nipt/analysis --apply  # 실제 삭제

BASE="${1:-$HOME/ken-nipt/analysis}"
APPLY=0
if [[ "${2-}" == "--apply" ]] || [[ "${1-}" == "--apply" ]]; then
  APPLY=1
  # 인자가 --apply 하나뿐인 경우 BASE 를 기본값으로 재설정
  [[ "${1-}" == "--apply" ]] && BASE="$HOME/ken-nipt/analysis"
fi

if [[ ! -d "$BASE" ]]; then
  echo "에러: 디렉터리가 존재하지 않습니다: $BASE" >&2
  exit 1
fi

echo "베이스 경로: $BASE"
echo "모드: $([[ $APPLY -eq 1 ]] && echo '실제 삭제' || echo '드라이런(목록만 출력)')"
echo

# sample_id 디렉터리(깊이 2)만 순회: ~/ken-nipt/analysis/<work_dir>/<sample_id>/
# 각 sample_id 디렉터리 직속의 .bam / .bai 파일을 대상으로 함
while IFS= read -r sample_dir; do
  echo ">> 처리 중: $sample_dir"

  # 지울 후보: *.bam, *.bai 중에서 아래 4종을 제외
  #   *proper_paired.bam(.bai), *of_orig.bam(.bai), *of_fetus.bam(.bai), *of_mom.bam(.bai)
  mapfile -t TO_DELETE < <(find "$sample_dir" -maxdepth 1 -type f \
    \( -name '*.bam' -o -name '*.bai' \) \
    -not \( \
        -name '*proper_paired.bam' -o -name '*proper_paired.bam.bai' -o \
        -name '*of_orig.bam'       -o -name '*of_orig.bam.bai'       -o \
        -name '*of_fetus.bam'      -o -name '*of_fetus.bam.bai'      -o \
        -name '*of_mom.bam'        -o -name '*of_mom.bam.bai' \
    \) \
    | sort)

  if [[ ${#TO_DELETE[@]} -eq 0 ]]; then
    echo "   삭제할 파일 없음."
    echo
    continue
  fi

  if [[ $APPLY -eq 1 ]]; then
    # 실제 삭제
    printf '%s\0' "${TO_DELETE[@]}" | xargs -0 -r rm -v --
  else
    # 드라이런: 목록만
    for f in "${TO_DELETE[@]}"; do
      echo "   (삭제 예정) $f"
    done
  fi
  echo
done < <(find "$BASE" -mindepth 2 -maxdepth 2 -type d | sort)

if [[ $APPLY -eq 0 ]]; then
  echo "드라이런이었습니다. 실제로 삭제하려면 --apply 옵션을 붙여 다시 실행하세요."
fi
