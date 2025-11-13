# Coverage Check Tools - 사용 가이드

## ⚠️ 중요: 어떤 도구를 사용해야 하나?

이 디렉토리에는 **2개의 커버리지 체크 도구**가 있습니다.

---

## 🎯 NIPT 분석용 (Shallow Depth WGS)

### ✅ `wisecondor_coverage_check.py` 사용!

```bash
# NIPT용 (0.2-0.3x 평균 깊이)
python3 wisecondor_coverage_check.py \
    --bam sample.proper_paired.bam \
    --bed MD_Target_8.bed \
    --output coverable_regions.bed \
    --report
```

**특징:**
- ✅ NIPT shallow depth WGS (0.2~0.3x)를 위해 설계
- ✅ Bin 단위로 read count 계산 (Wisecondor와 동일)
- ✅ 기본값: 200kb bin, 최소 10 reads/bin
- ✅ 통계적으로 의미 있는 분석

**왜 이 도구를 사용해야 하나?**
```
NIPT는 평균 0.2-0.3x 깊이:
- 대부분의 염기는 커버리지 0
- Per-base coverage는 의미 없음
- Bin 단위 read count만 의미 있음

0.2x 깊이에서:
→ 200kb bin당 약 40 reads
→ 통계적 분석 가능
→ Wisecondor/WisecondorX와 호환
```

📖 **상세 문서:** [README_wisecondor_coverage_check.md](README_wisecondor_coverage_check.md)

---

## 🔬 일반 WGS용 (High Depth)

### ❌ `md_coverage_check.py` (NIPT에는 부적합!)

```bash
# 일반 WGS용 (30x 이상)
python3 md_coverage_check.py \
    --bam sample.bam \
    --bed targets.bed \
    --output output.bed \
    --min-coverage 5 \
    --coverage-threshold 0.95
```

**특징:**
- ✅ 일반 WGS (30x 이상)를 위해 설계
- ✅ Per-base coverage 계산
- ✅ 기본값: 1x coverage, 90% threshold
- ❌ **NIPT에는 부적합** (너무 높은 threshold)

**왜 NIPT에 부적합한가?**
```
일반 WGS는 30x 이상의 깊이:
- 모든 염기가 높은 커버리지
- Per-base 분석이 의미 있음
- 1x threshold도 충분히 낮음

NIPT는 0.2-0.3x:
- 1x threshold는 너무 높음
- Per-base는 대부분 0
- 이 도구로는 분석 불가능!
```

📖 **상세 문서:** [README_md_coverage_check.md](README_md_coverage_check.md)

---

## 빠른 비교표

| 항목 | wisecondor_coverage_check.py | md_coverage_check.py |
|------|------------------------------|----------------------|
| **대상 깊이** | 0.2-0.3x (NIPT) | 30x+ (일반 WGS) |
| **분석 단위** | Bin (200kb) | Per-base |
| **측정값** | Bin당 총 read 수 | 각 염기의 깊이 |
| **기본 threshold** | 10 reads/bin | 1x per-base |
| **NIPT 적합성** | ✅ 적합 | ❌ 부적합 |
| **Wisecondor 호환** | ✅ 완전 호환 | ❌ 다른 방식 |

---

## 예시: 동일한 샘플, 다른 결과

### NIPT 샘플 (0.25x 평균 깊이)

**`wisecondor_coverage_check.py` 결과:**
```
✓ PASS MD_22q11.2
  Total reads: 1,000
  Mean reads per bin: 50
  → 분석 가능!
```

**`md_coverage_check.py` 결과:**
```
✗ FAIL MD_22q11.2
  Mean coverage: 0.25x
  Coverage fraction: 25%
  → 분석 불가능! (threshold 1x)
```

**결론:**
- 동일한 샘플이지만 도구에 따라 결과가 완전히 다름
- NIPT에는 **반드시** `wisecondor_coverage_check.py` 사용!

---

## 자주 묻는 질문 (FAQ)

### Q1: 왜 두 개의 도구가 필요한가?

**A:** 시퀀싱 깊이가 다르면 분석 방법도 달라야 합니다.

```
NIPT (0.2-0.3x):
- Shallow depth: 비용 절감
- Bin 단위 분석: CNV 탐지 충분
- Per-base는 의미 없음

일반 WGS (30x+):
- High depth: 정확한 변이 탐지
- Per-base 분석: SNV, InDel 탐지 가능
- Bin 분석은 해상도 낮음
```

### Q2: NIPT에 md_coverage_check.py를 사용하면?

**A:** 거의 모든 영역이 "분석 불가능"으로 나옵니다.

```
0.25x 평균 깊이:
- Per-base로 1x 이상인 염기: ~25%
- Coverage threshold 90% 만족: 불가능
- 결과: ✗ FAIL (실제로는 분석 가능한데도!)
```

### Q3: 일반 WGS에 wisecondor_coverage_check.py를 사용하면?

**A:** 사용 가능하지만, 해상도가 낮습니다.

```
30x 깊이에서 200kb bin:
- Bin당 6,000 reads (충분!)
- 하지만 200kb 해상도는 너무 낮음
- SNV, InDel 탐지는 불가능
→ Per-base 분석이 더 적합
```

### Q4: 중간 깊이 (5-10x)는?

**A:** 목적에 따라 다릅니다.

```
CNV 탐지만 필요:
→ wisecondor_coverage_check.py

정확한 변이 탐지 필요:
→ md_coverage_check.py (threshold 조정)
```

---

## 실전 사용 가이드

### NIPT 샘플 분석 워크플로우

```bash
#!/bin/bash
# NIPT 샘플 커버리지 체크 워크플로우

SAMPLE_ID="OPC241100001"
BAM="/Work/NIPT/analysis/2411/${SAMPLE_ID}/${SAMPLE_ID}.proper_paired.bam"
BED="/Work/NIPT/data/bed/common/MD_Target_8.bed"
OUT_DIR="/Work/NIPT/analysis/2411/${SAMPLE_ID}/coverage_check"

mkdir -p "${OUT_DIR}"

# 1. 깊이 확인
echo "1. Checking sequencing depth..."
TOTAL_READS=$(samtools idxstats "${BAM}" | awk '{sum+=$3} END {print sum}')
DEPTH=$(echo "scale=2; ${TOTAL_READS} * 150 / 3000000000" | bc)
echo "   Total reads: ${TOTAL_READS}"
echo "   Estimated depth: ${DEPTH}x"

# 2. 커버리지 체크 (NIPT용 도구 사용!)
echo "2. Checking coverage with wisecondor_coverage_check.py..."
python3 wisecondor_coverage_check.py \
    --bam "${BAM}" \
    --bed "${BED}" \
    --output "${OUT_DIR}/coverable_regions.bed" \
    --report

# 3. 결과 확인
echo "3. Review results:"
echo "   Report: ${OUT_DIR}/coverable_regions_report.txt"
echo "   BED: ${OUT_DIR}/coverable_regions.bed"

# 4. Wisecondor 분석 (coverable regions 사용)
echo "4. Run Wisecondor analysis with coverable regions..."
# ... (Wisecondor 분석 코드)
```

---

## 요약

### ✅ NIPT 분석 시

```bash
# 이 명령어를 사용하세요!
python3 wisecondor_coverage_check.py \
    --bam sample.proper_paired.bam \
    --bed MD_Target_8.bed \
    --output coverable.bed \
    --report
```

**이유:**
- NIPT는 0.2-0.3x shallow depth
- Bin 단위 read count가 유일한 의미 있는 지표
- Wisecondor/WisecondorX와 호환

### ❌ 잘못된 사용

```bash
# NIPT에 이 명령어 사용하지 마세요!
python3 md_coverage_check.py \
    --bam sample.proper_paired.bam \
    --bed MD_Target_8.bed \
    --output output.bed
# → 결과: 모든 영역 분석 불가능 (잘못된 결과!)
```

---

## 문의 및 지원

문제가 발생하면:
1. 먼저 도구의 README 파일 확인
2. 샘플의 시퀀싱 깊이 확인
3. 올바른 도구 사용 여부 확인

**기억하세요:**
- NIPT = `wisecondor_coverage_check.py` ✅
- 일반 WGS = `md_coverage_check.py` ✅

---

Created by Ken
Version: 1.0

