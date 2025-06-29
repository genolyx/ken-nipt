#!/usr/bin/env bash
set -euo pipefail

R1="$1"
R2="$2"
SAMPLE_READS=${SAMPLE_READS:-1000}

echo "▶ 단순 페어링 검사"
echo "  • R1: $(basename "$R1")"
echo "  • R2: $(basename "$R2")"

TMP=$(mktemp -d)
cleanup() { rm -rf "$TMP"; }
trap cleanup EXIT

# R1 헤더 추출
echo "  • R1 처리 중..."
gunzip -c "$R1" | head -n $(( SAMPLE_READS * 4 )) | awk 'NR%4==1 { sub(/^@/, ""); print }' > "$TMP/r1_full.txt"
R1_COUNT=$(wc -l < "$TMP/r1_full.txt")
echo "    ✓ R1 완료: $R1_COUNT 개"

# R2 헤더 추출
echo "  • R2 처리 중..."
gunzip -c "$R2" | head -n $(( SAMPLE_READS * 4 )) | awk 'NR%4==1 { sub(/^@/, ""); print }' > "$TMP/r2_full.txt"
R2_COUNT=$(wc -l < "$TMP/r2_full.txt")
echo "    ✓ R2 완료: $R2_COUNT 개"

# 첫 5개 헤더 비교
echo "  • 헤더 비교:"
for i in {1..5}; do
    R1_HEADER=$(sed -n "${i}p" "$TMP/r1_full.txt")
    R2_HEADER=$(sed -n "${i}p" "$TMP/r2_full.txt")
    
    echo "    $i. R1: $R1_HEADER"
    echo "       R2: $R2_HEADER"
    
    # ID 부분 추출 (공백 앞)
    R1_ID="${R1_HEADER%% *}"
    R2_ID="${R2_HEADER%% *}"
    
    if [[ "$R1_ID" == "$R2_ID" ]]; then
        echo "       ✓ ID 일치"
    else
        echo "       ✗ ID 불일치"
    fi
    
    # 메타데이터 부분 (공백 뒤)
    if [[ "$R1_HEADER" == *" "* ]]; then
        R1_META="${R1_HEADER#* }"
        R2_META="${R2_HEADER#* }"
        
        echo "       R1 메타: $R1_META"
        echo "       R2 메타: $R2_META"
        
        if [[ "$R1_META" == "$R2_META" ]]; then
            echo "       ⚠ 메타데이터 동일 (복사 파일 의심)"
        else
            echo "       ✓ 메타데이터 다름"
        fi
    fi
    echo ""
done

# 전체 일치 검사
SAME_COUNT=0
for i in $(seq 1 $R1_COUNT); do
    R1_LINE=$(sed -n "${i}p" "$TMP/r1_full.txt")
    R2_LINE=$(sed -n "${i}p" "$TMP/r2_full.txt")
    
    if [[ "$R1_LINE" == "$R2_LINE" ]]; then
        SAME_COUNT=$((SAME_COUNT + 1))
    fi
done

echo "  • 결과:"
echo "    전체 동일 헤더: $SAME_COUNT / $R1_COUNT"

if [[ $SAME_COUNT -eq $R1_COUNT ]]; then
    echo "    ✗ 모든 헤더가 동일합니다. 복사 파일일 가능성이 높습니다!"
elif [[ $SAME_COUNT -eq 0 ]]; then
    echo "    ✓ 정상적인 페어링입니다."
else
    echo "    ⚠ 일부 헤더가 동일합니다. 확인이 필요합니다."
fi

echo "▶ 검사 완료"
