#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 <target_dir>"
  exit 1
}

# 인자 체크
if [ $# -ne 1 ]; then
  usage
fi

TARGET_DIR="$1"

# 디렉터리가 존재하는지 확인
if [ ! -d "$TARGET_DIR" ]; then
  echo "Error: Directory '$TARGET_DIR' not found." >&2
  exit 1
fi

cd "$TARGET_DIR"

# 모든 R1 파일을 순회하며 sample_name 추출 → R2 찾고 이동
for fq1 in *_R1_*.fastq.gz; do
  # 파일이 없으면 루프 탈출
  [[ ! -e $fq1 ]] && break

  sample_name="${fq1%%_*}"     # 언더바(_) 앞부분을 sample_name 으로
  fq2="${fq1/_R1_/_R2_}"       # R1 → R2 치환

  # R2 파일 존재 체크
  if [ ! -e "$fq2" ]; then
    echo "⚠️  Warning: R2 file for '$fq1' not found, skipping." >&2
    continue
  fi

  echo "Processing sample: $sample_name"

  # sample_name 폴더 생성 후 파일 이동
  mkdir -p "$sample_name"
  mv -- "$fq1" "$fq2" "$sample_name"/
done

echo "Done."
