# Wisecondor Coverage Check Tool (NIPT용)

## ⚠️ 중요: 어떤 도구를 사용해야 하나?

이 디렉토리에는 **2개의 커버리지 체크 도구**가 있습니다:

### 1. `wisecondor_coverage_check.py` ✅ **NIPT용 (이 도구 사용!)**
- **NIPT shallow depth WGS (0.2~0.3x)를 위해 설계됨**
- Bin 단위로 read count 계산
- Wisecondor/WisecondorX와 동일한 방식으로 분석
- 기본값: bin 200kb, 최소 10 reads/bin

### 2. `md_coverage_check.py` ❌ **일반 WGS용 (NIPT에 부적합)**
- 일반 WGS (30x 이상)를 위해 설계됨
- Per-base coverage 계산
- NIPT의 shallow depth에는 적합하지 않음

---

## 개요 (Overview)

`wisecondor_coverage_check.py`는 NIPT에서 사용하는 **shallow depth WGS (0.2~0.3x 평균 깊이)**를 위해 특별히 설계된 도구입니다.

Wisecondor/WisecondorX와 동일한 방식으로 bin 단위의 read count를 분석하여, 각 BED 영역이 분석 가능한지 판단합니다.

## 왜 Bin 단위로 분석하나? (Why Bins?)

NIPT는 shallow depth WGS를 사용하기 때문에:
- 평균 깊이가 0.2~0.3x 밖에 안 됨
- **대부분의 염기 위치는 커버리지가 0임**
- Per-base coverage는 의미가 없음

대신, Wisecondor/WisecondorX는:
- 큰 bin (기본 200kb)으로 게놈을 나눔
- 각 bin 내의 **총 read count**를 계산
- Bin 간의 read count 비율로 CNV 탐지

따라서 이 도구도 동일한 방식으로 분석합니다!

## 분석 기준 (Analysis Criteria)

### Per-base Coverage (❌ 부적합)
```
각 염기 위치의 깊이를 확인
→ 0.2x 평균 깊이에서는 대부분이 0
→ NIPT에 적합하지 않음
```

### Bin-level Read Count (✅ 적합)
```
200kb bin 내의 총 read 수를 계산
→ 0.2x 깊이에서도 bin당 40 reads 예상 (200kb × 0.2x = 40 reads)
→ 통계적으로 의미 있는 분석 가능
→ Wisecondor/WisecondorX와 동일한 방식
```

## 사용법 (Usage)

### 기본 사용 (Basic Usage)

```bash
python3 wisecondor_coverage_check.py \
    --bam sample.proper_paired.bam \
    --bed MD_Target_8.bed \
    --output coverable_regions.bed
```

### 고급 사용 (Advanced Usage)

```bash
python3 wisecondor_coverage_check.py \
    --bam sample.proper_paired.bam \
    --bed MD_Target_8.bed \
    --output coverable_regions.bed \
    --bin-size 200000 \
    --min-reads-per-bin 10 \
    --min-bin-fraction 0.8 \
    --min-total-reads 1000 \
    --report
```

## 파라미터 (Parameters)

### 필수 파라미터 (Required)

- `--bam`: 입력 BAM 파일 (proper_paired.bam 권장)
- `--bed`: 타겟 영역이 정의된 BED 파일
- `--output`, `-o`: 분석 가능한 영역을 출력할 BED 파일

### 선택 파라미터 (Optional)

- `--bin-size`: Wisecondor bin 크기 (기본값: 200000 bp)
  - Wisecondor/WisecondorX의 bin size와 일치시켜야 함
- `--min-reads-per-bin`: Bin당 최소 read 수 (기본값: 10)
  - 0.2x 깊이에서 200kb bin은 약 40 reads
  - 통계적 유의성을 위해 최소 10 reads 권장
- `--min-bin-fraction`: 충분한 reads를 가진 bin의 최소 비율 (기본값: 0.8)
  - 영역의 80% 이상의 bin이 충분한 reads를 가져야 함
- `--min-total-reads`: 영역의 최소 총 read 수 (기본값: 1000)
- `--report`: 상세 리포트 생성
- `--no-partial`: 부분적으로 커버된 영역을 출력하지 않음

## 기본값 설정 근거 (Default Values Rationale)

### Bin Size: 200kb
- Wisecondor/WisecondorX 기본값과 동일
- 0.2x 깊이에서 bin당 40 reads 예상
- CNV 탐지에 적합한 해상도

### Min Reads per Bin: 10
```
0.2x 평균 깊이:
- 200kb bin × 0.2x = 40 reads (평균)
- Poisson 분포를 고려하여 최소 10 reads 권장
- 너무 낮은 read count는 통계적 신뢰도 떨어짐
```

### Min Bin Fraction: 0.8
- 영역의 80% 이상이 충분한 커버리지를 가져야 함
- Wisecondor 분석의 안정성 보장

### Min Total Reads: 1000
```
예: 5 Mb 영역
- 0.2x 깊이: 5,000,000 bp × 0.2x = 1,000 reads (최소)
- 통계적으로 의미 있는 분석을 위한 최소값
```

## 출력 파일 (Output Files)

### 1. Coverable Regions BED File

```
# Format: chrom, start, end, name, score, strand, total_reads, mean_reads, bin_fraction
chr22    18000000    22000000    MD_22q11.2    1000    .    850    42.50    0.9500
chr15    22700000    25200000    MD_15q11.2_partial1    500    .    520    26.00    1.0000
```

**컬럼 설명:**
- `total_reads`: 영역 내 총 read 수
- `mean_reads`: Bin당 평균 read 수
- `bin_fraction`: 충분한 reads를 가진 bin의 비율

### 2. Coverage Report (with --report)

```
================================================================================
Wisecondor Coverage Analysis Report (NIPT Shallow WGS)
================================================================================

BAM File: sample.proper_paired.bam
Total Mapped Reads: 15,234,567
Estimated Average Depth: 0.25x

Analysis Parameters:
  - Bin size: 200,000 bp
  - Minimum reads per bin: 10
  - Minimum bin fraction: 80.00%
  - Minimum total reads: 1,000

Summary:
  - Total regions: 8
  - Fully analyzable: 6 (75.0%)
  - Partially covered: 1 (12.5%)
  - Not analyzable: 1 (12.5%)

Read Count Statistics:
  - Total reads (avg per region): 2,450
  - Mean reads per bin (avg): 45.20
  - Mean reads per bin (median): 42.00

================================================================================
Per-Region Analysis
================================================================================

✓ PASS MD_22q11.2 (chr22:18000000-22000000)
  Size: 4,000,000 bp
  Number of bins: 20
  Total reads: 850
  Mean reads per bin: 42.50
  Median reads per bin: 40.00
  Read range per bin: 30-55
  Bins with sufficient reads: 19/20 (95.00%)
  Bin-level reads: [45, 42, 38, 50, 41, 39, 43, 46, 40, 44, 30, 48, 42, 41, 55, 43, 40, 38, 45, 42]

✗ FAIL MD_1p36 (chr1:1000000-5000000)
  Size: 4,000,000 bp
  Number of bins: 20
  Total reads: 320
  Mean reads per bin: 16.00
  Median reads per bin: 12.00
  Read range per bin: 5-35
  Bins with sufficient reads: 8/20 (40.00%)
  Reason: Low bin fraction (40.00% < 80.00%)
  Coverable bin groups: 1
    1. chr1:2400000-4200000 (9 bins, 1,800,000 bp)
```

## 실제 사용 예시 (Real-world Examples)

### 예시 1: MD 타겟 분석

```bash
SAMPLE_ID="OPC241100001"
BAM="/Work/NIPT/analysis/2411/${SAMPLE_ID}/${SAMPLE_ID}.proper_paired.bam"
BED="/Work/NIPT/data/bed/common/MD_Target_8.bed"
OUT="/Work/NIPT/analysis/2411/${SAMPLE_ID}/coverable_MD_Target_8.bed"

python3 wisecondor_coverage_check.py \
    --bam "${BAM}" \
    --bed "${BED}" \
    --output "${OUT}" \
    --report
```

### 예시 2: 여러 타겟 배치 분석

```bash
#!/bin/bash
SAMPLE_ID="OPC241100001"
BAM="/Work/NIPT/analysis/2411/${SAMPLE_ID}/${SAMPLE_ID}.proper_paired.bam"
OUT_DIR="/Work/NIPT/analysis/2411/${SAMPLE_ID}/coverage_check"

mkdir -p "${OUT_DIR}"

# 모든 MD 타겟 분석
for TARGET in MD_Target_1 MD_Target_4 MD_Target_5 MD_Target_8 MD_Target_15 MD_Target_22; do
    echo "Analyzing ${TARGET}..."
    python3 wisecondor_coverage_check.py \
        --bam "${BAM}" \
        --bed "/Work/NIPT/data/bed/common/${TARGET}.bed" \
        --output "${OUT_DIR}/coverable_${TARGET}.bed" \
        --report
done

echo "All targets analyzed. Results in: ${OUT_DIR}"
```

### 예시 3: 낮은 깊이 샘플을 위한 완화된 기준

```bash
# 낮은 깊이 샘플 (0.15x)을 위한 완화된 파라미터
python3 wisecondor_coverage_check.py \
    --bam low_depth_sample.bam \
    --bed MD_Target_8.bed \
    --output coverable_relaxed.bed \
    --min-reads-per-bin 5 \
    --min-bin-fraction 0.7 \
    --min-total-reads 500 \
    --report
```

## 결과 해석 (Interpreting Results)

### 분석 가능 (Analyzable) ✓
```
영역의 80% 이상의 bin이 충분한 reads (≥10)를 가짐
→ Wisecondor/WisecondorX 분석 가능
→ 신뢰할 수 있는 CNV 탐지 가능
```

### 부분 커버 (Partially Covered) ⚠
```
일부 bin만 충분한 reads를 가짐
→ 커버 가능한 bin group이 출력 BED에 포함됨
→ 해당 영역만 사용하여 분석 가능
→ 전체 타겟 영역은 아니지만 일부 정보 얻을 수 있음
```

### 분석 불가 (Not Analyzable) ✗
```
충분한 reads를 가진 bin이 너무 적음
→ Wisecondor 분석 불가능
→ 가능한 원인:
  - 시퀀싱 깊이가 너무 낮음 (<0.1x)
  - BAM 품질 문제
  - 영역이 mappability가 낮은 곳에 위치
```

## 예상 Read Count (Expected Read Counts)

### 0.2x 평균 깊이 (일반적인 NIPT)
```
200kb bin: 200,000 bp × 0.2x = 40 reads (평균)
1 Mb 영역: 1,000,000 bp × 0.2x = 200 reads
4 Mb 영역: 4,000,000 bp × 0.2x = 800 reads
```

### 0.3x 평균 깊이 (높은 품질 NIPT)
```
200kb bin: 200,000 bp × 0.3x = 60 reads (평균)
1 Mb 영역: 1,000,000 bp × 0.3x = 300 reads
4 Mb 영역: 4,000,000 bp × 0.3x = 1,200 reads
```

### 0.15x 평균 깊이 (낮은 품질)
```
200kb bin: 200,000 bp × 0.15x = 30 reads (평균)
→ 파라미터 조정 필요 (--min-reads-per-bin 5)
```

## 문제 해결 (Troubleshooting)

### 모든 영역이 분석 불가능 ✗

**1. 깊이 확인**
```bash
# Total mapped reads 확인
samtools idxstats sample.bam | awk '{sum+=$3} END {print "Total mapped:", sum}'

# 예상 깊이 계산
# Depth = (Total reads × Read length) / Genome size
# Depth = (15,000,000 × 150) / 3,000,000,000 = 0.75x
```

**2. 파라미터 완화**
```bash
# 낮은 깊이를 위한 완화된 파라미터
python3 wisecondor_coverage_check.py \
    --bam sample.bam \
    --bed targets.bed \
    --output output.bed \
    --min-reads-per-bin 5 \
    --min-bin-fraction 0.6 \
    --min-total-reads 300
```

**3. BAM 품질 확인**
```bash
samtools flagstat sample.bam
samtools stats sample.bam | grep "average length"
```

### 일부 영역만 분석 불가능 ⚠

**특정 영역의 문제일 수 있음:**
- Mappability가 낮은 영역
- Repetitive sequences
- Segmental duplications

**해결 방법:**
1. 부분 커버된 영역 사용 (`--no-partial` 제거)
2. 해당 영역 제외하고 분석
3. 더 나은 시퀀싱 데이터 생성

## Wisecondor와의 호환성

이 도구는 Wisecondor/WisecondorX와 **완전히 호환**됩니다:

1. **동일한 bin 크기**: 기본 200kb
2. **동일한 read counting 방식**: pysam.count()
3. **동일한 분석 단위**: Bin-level read counts

따라서 이 도구로 확인된 영역은 Wisecondor/WisecondorX에서 **신뢰할 수 있게 분석 가능**합니다!

## 다음 단계 (Next Steps)

1. **리포트 확인**
   ```bash
   less coverable_regions_report.txt
   ```

2. **Coverable regions BED 파일 사용**
   ```bash
   # MD 분석 시 coverable regions만 사용
   # (파이프라인에 통합 가능)
   ```

3. **필요시 파라미터 조정**
   - 낮은 깊이: 파라미터 완화
   - 높은 품질 요구: 파라미터 강화

## 기술 세부사항 (Technical Details)

### Read Counting 방식
```python
# pysam을 사용한 bin 단위 read count
for each bin (200kb):
    count = bamfile.count(chrom, bin_start, bin_end)
    # Properly paired reads만 카운트 (proper_paired.bam 사용 시)
```

### 통계적 근거
```
Poisson 분포를 따르는 read count:
- 평균 40 reads/bin (0.2x)
- 표준편차 √40 ≈ 6.3
- 최소 10 reads = 평균의 25% ≈ -4.8σ
- 충분한 통계적 유의성 확보
```

### Bin Group 탐지
```
연속된 sufficient bin들을 그룹화:
- 각 bin이 min_reads_per_bin 이상인지 확인
- 연속된 sufficient bin을 하나의 그룹으로 묶음
- 최소 1개 bin 이상인 그룹만 출력
```

## 비교: 두 도구의 차이

| 특성 | wisecondor_coverage_check.py | md_coverage_check.py |
|------|------------------------------|----------------------|
| 대상 | NIPT (0.2-0.3x) | 일반 WGS (30x+) |
| 분석 단위 | Bin (200kb) | Per-base |
| Read count | Bin당 총 reads | 각 염기의 깊이 |
| 기본 threshold | 10 reads/bin | 1x coverage |
| 적합성 | ✅ NIPT에 적합 | ❌ NIPT에 부적합 |

## 저자 (Author)
Ken

## 버전 (Version)
1.0

## 라이선스 (License)
Internal use only - NIPT Pipeline

