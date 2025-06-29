#!/usr/bin/env bash

F1="$1"
F2="$2"

echo "▶ 초간단 파일 비교"
echo "파일1: $(basename "$F1")"
echo "파일2: $(basename "$F2")"

# 각 파일에서 첫 100라인만 추출
echo "파일1 처리 중..."
gunzip -c "$F1" | head -1000 > /tmp/f1.txt

echo "파일2 처리 중..."  
gunzip -c "$F2" | head -1000 > /tmp/f2.txt

echo "첫 몇 라인 비교:"
echo "=== 파일1 처음 4라인 ==="
head -4 /tmp/f1.txt

echo "=== 파일2 처음 4라인 ==="
head -4 /tmp/f2.txt

# 동일성 체크
if cmp -s /tmp/f1.txt /tmp/f2.txt; then
    echo "✗ 결론: 두 파일이 완전히 동일합니다!"
else
    echo "✓ 결론: 두 파일이 다릅니다."
    
    # 차이점 보기
    echo "차이점:"
    diff /tmp/f1.txt /tmp/f2.txt | head -10
fi

# 정리
rm -f /tmp/f1.txt /tmp/f2.txt

echo "완료!"
