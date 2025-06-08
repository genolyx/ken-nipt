#!/bin/bash

# 인자 확인
if [ $# -ne 1 ]; then
  echo "Usage: $0 <target_directory>"
  exit 1
fi

# 전달받은 디렉토리
target_dir="$1"

# 절대 경로로 변환
base_dir=$(realpath "$target_dir")

# 디렉토리 존재 확인
if [ ! -d "$base_dir" ]; then
  echo "Error: Directory '$base_dir' does not exist."
  exit 1
fi

# 하위 디렉토리만 순회
cd "$base_dir" || exit 1
for dir in */ ; do
  # 디렉토리 이름에서 마지막 '/' 제거
  dirname="${dir%/}"

  # tar.gz 파일 생성
  tar --exclude='*.bam' --exclude='*.bai' -czvf "${base_dir}/${dirname}.tar.gz" -C "$base_dir" "$dirname"

  echo "Created archive: ${dirname}.tar.gz"
done
