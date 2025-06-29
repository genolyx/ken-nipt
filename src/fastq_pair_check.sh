#!/usr/bin/env bash

R1="$1"
R2="$2"
SAMPLE_READS=${SAMPLE_READS:-1000}

# 기본 검증
if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "Error: File not found"
    exit 1
fi

# 임시 파일
TMP=$(mktemp -d)
trap "rm -rf $TMP" EXIT

# 샘플 추출
SAMPLE_LINES=$(( SAMPLE_READS * 4 ))
gunzip -c "$R1" | head -n $SAMPLE_LINES > "$TMP/r1.txt" 2>/dev/null
gunzip -c "$R2" | head -n $SAMPLE_LINES > "$TMP/r2.txt" 2>/dev/null

# 완전 동일 파일 검사
if cmp -s "$TMP/r1.txt" "$TMP/r2.txt"; then
    echo "Same"
    exit 2
fi

# 헤더 추출
awk 'NR%4==1 { sub(/^@/, ""); print }' "$TMP/r1.txt" > "$TMP/r1_h.txt"
awk 'NR%4==1 { sub(/^@/, ""); print }' "$TMP/r2.txt" > "$TMP/r2_h.txt"

# ID 매칭 검사
R1_COUNT=$(wc -l < "$TMP/r1_h.txt")
MATCHING=0

for i in $(seq 1 $R1_COUNT); do
    R1_ID=$(sed -n "${i}p" "$TMP/r1_h.txt" | cut -d' ' -f1)
    R2_ID=$(sed -n "${i}p" "$TMP/r2_h.txt" | cut -d' ' -f1)
    
    if [ "$R1_ID" = "$R2_ID" ]; then
        MATCHING=$((MATCHING + 1))
    fi
done

# 결과 판정
MATCH_RATE=$((MATCHING * 100 / R1_COUNT))

if [ $MATCH_RATE -ge 95 ]; then
    echo "Match"
    exit 0
else
    echo "Unmatch"
    exit 1
fi
