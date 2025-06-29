#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "사용법: $0 <대상_디렉토리> <rename.tsv>"
  exit 1
fi

TARGET_DIR=$1
MAP_FILE=$2

# 매핑 파일 존재 확인
if [ ! -f "$MAP_FILE" ]; then
  echo "오류: 매핑 파일 '$MAP_FILE'을 찾을 수 없습니다."
  exit 1
fi

# tsv 파일을 한 줄씩 읽어서 src → target 처리
while IFS=$'\t' read -r SRC TARGET; do
  # 빈 줄 또는 주석(#)으로 시작하면 건너뛰기
  [[ -z "$SRC" || "${SRC:0:1}" == "#" ]] && continue

  echo ">>> 매핑: '$SRC' → '$TARGET'"

  # 깊이 우선(depth-first)으로 하위 항목부터 처리
  find "$TARGET_DIR" -depth -name "${SRC}*" | while read -r ITEM; do
    DIR=$(dirname "$ITEM")
    BASE=$(basename "$ITEM")
    # 파일명 또는 디렉토리명 앞부분 SRC를 TARGET으로 바꿔치기
    NEWBASE="${BASE/#$SRC/$TARGET}"
    if [ "$BASE" != "$NEWBASE" ]; then
      echo "리네임: '$ITEM' → '$DIR/$NEWBASE'"
      mv "$ITEM" "$DIR/$NEWBASE"
    fi
  done

done < "$MAP_FILE"
