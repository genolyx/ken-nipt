#!/usr/bin/env bash
set -euo pipefail

BASE="/home/ken/ken-nipt/fastq/2504"
TSV="OrderId_SampleId_table.tsv"

echo "베이스 디렉토리 : $BASE"
echo "TSV 파일         : $TSV"
echo

# TSV 파일을 읽어서 두 번째 컬럼(old) → 첫 번째 컬럼(new)으로 이름 변경
while IFS=$'\t' read -r new_id old_id || [[ -n "$new_id" ]]; do
  # 빈줄/주석 무시
  [[ -z "${new_id// }" || "$new_id" =~ ^# ]] && continue

  OLD_DIR="$BASE/$old_id"
  NEW_DIR="$BASE/$new_id"

  if [[ ! -d "$OLD_DIR" ]]; then
    echo "[WARN] 존재하지 않음: $OLD_DIR"
    continue
  fi

  if [[ -e "$NEW_DIR" ]]; then
    echo "[WARN] 이미 존재함: $NEW_DIR (건너뜀)"
    continue
  fi

  echo ">>> mv $OLD_DIR  →  $NEW_DIR"
  mv "$OLD_DIR" "$NEW_DIR"
done < "$TSV"

echo
echo "✅ 완료! 디렉토리 이름이 변경되었습니다."
