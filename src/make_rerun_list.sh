#!/usr/bin/env bash
set -euo pipefail

# 1) 설정: 필요 시 수정
ANALYSIS_BASE="$HOME/ken-nipt/analysis"   # analysis root
OUT="sample_sheet.tsv"                    # 출력 TSV
LAB="cordlife"
AGE="30"

# 2) 헤더 작성
echo -e "SAMPLE_NAME\tWORK_DIR\tFQ1\tFQ2\tAGE\tLAB" > "$OUT"

# 3) 인자로 WORK_DIR 주어졌으면 그 하나만, 없으면 모든 서브디렉토리
if [[ $# -gt 0 ]]; then
    WORK_DIRS=("$ANALYSIS_BASE/$1")
else
    WORK_DIRS=("$ANALYSIS_BASE"/*)
fi

# 4) 각 WORK_DIR 순회
for wd_path in "${WORK_DIRS[@]}"; do
    [[ -d "$wd_path" ]] || continue
    WORK_DIR=$(basename "$wd_path")

    # 5) 각 샘플 디렉토리 순회
    for sample_path in "$wd_path"/*; do
        [[ -d "$sample_path" ]] || continue
        SAMPLE=$(basename "$sample_path")

        # 6) FASTQ 이름 패턴
        FQ1="${SAMPLE}_R1.fastq.gz"
        FQ2="${SAMPLE}_R2.fastq.gz"

        # 7) 한 줄 추가
        echo -e "${SAMPLE}\t${WORK_DIR}\t${FQ1}\t${FQ2}\t${AGE}\t${LAB}"
    done
done >> "$OUT"

echo "✅ Generated $OUT ($(wc -l < "$OUT" | tr -d ' ') lines including header)."
