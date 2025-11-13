# MD Coverage Check Tool

## 개요 (Overview)

`md_coverage_check.py`는 proper_paired.bam 파일과 BED 파일을 입력받아 각 BED 영역에 대한 커버리지를 분석하고, Wisecondor/WisecondorX 분석이 가능한지 여부를 판단하는 도구입니다.

This tool analyzes BAM coverage over BED regions to determine if Wisecondor/WisecondorX analysis is feasible, and outputs coverable regions in BED format.

## 주요 기능 (Key Features)

1. **커버리지 분석**: 각 BED 영역에 대한 상세한 커버리지 통계 제공
2. **분석 가능성 판단**: Wisecondor/WisecondorX 분석에 충분한 커버리지 확인
3. **부분 영역 탐지**: 전체 영역을 커버하지 못하는 경우, 커버 가능한 부분 영역 탐지 및 출력
4. **상세 리포트**: 옵션으로 상세한 분석 리포트 생성

## 설치 (Installation)

필수 패키지:
```bash
pip install pysam numpy pandas
```

## 사용법 (Usage)

### 기본 사용 (Basic Usage)

```bash
python3 md_coverage_check.py \
    --bam sample.proper_paired.bam \
    --bed MD_Target_8.bed \
    --output coverable_regions.bed
```

### 고급 사용 (Advanced Usage)

```bash
python3 md_coverage_check.py \
    --bam sample.proper_paired.bam \
    --bed MD_Target_8.bed \
    --output coverable_regions.bed \
    --min-coverage 5 \
    --coverage-threshold 0.95 \
    --min-region-size 5000 \
    --bin-size 200000 \
    --report
```

## 파라미터 (Parameters)

### 필수 파라미터 (Required)

- `--bam`: 입력 BAM 파일 (proper_paired.bam 권장)
- `--bed`: 타겟 영역이 정의된 BED 파일
- `--output`, `-o`: 분석 가능한 영역을 출력할 BED 파일

### 선택 파라미터 (Optional)

- `--min-coverage`: 최소 read depth (기본값: 1)
- `--min-region-size`: 최소 영역 크기 (bp 단위, 기본값: 1000)
- `--coverage-threshold`: 최소 커버된 염기 비율 (0.0-1.0, 기본값: 0.9)
  - 0.9 = 영역의 90% 이상이 최소 커버리지를 만족해야 함
- `--bin-size`: Wisecondor 분석용 bin 크기 (기본값: 200000)
- `--report`: 상세 리포트 생성
- `--no-partial`: 부분적으로 커버된 영역을 출력하지 않음

## 출력 파일 (Output Files)

### 1. Coverable Regions BED File

커버 가능한 영역을 포함하는 BED 파일:

```
# Format: chrom, start, end, name, score, strand, mean_cov, median_cov, coverage_frac
chr1    1000000    2000000    region1    1000    .    15.50    14.00    0.9850
chr2    3000000    3500000    region2_partial1    500    .    10.20    10.00    0.9200
```

**컬럼 설명:**
- `chrom`: 염색체
- `start`: 시작 위치
- `end`: 종료 위치
- `name`: 영역 이름 (부분 영역의 경우 `_partial1`, `_partial2` 등 접미사 추가)
- `score`: 1000 (완전 커버), 500 (부분 커버)
- `strand`: 방향 (항상 `.`)
- `mean_cov`: 평균 커버리지
- `median_cov`: 중간값 커버리지
- `coverage_frac`: 커버된 염기 비율

### 2. Coverage Report (optional, with --report flag)

상세한 분석 리포트:

```
================================================================================
MD Coverage Analysis Report
================================================================================

BAM File: sample.proper_paired.bam
Analysis Parameters:
  - Minimum coverage: 1x
  - Minimum region size: 1000 bp
  - Coverage threshold: 90.00%
  - Bin size (Wisecondor): 200000 bp

Summary:
  - Total regions: 10
  - Fully analyzable: 7 (70.0%)
  - Partially covered: 2 (20.0%)
  - Not analyzable: 1 (10.0%)

Coverage Statistics:
  - Mean coverage (avg): 12.45x
  - Mean coverage (median): 11.20x
  - Mean coverage (min): 0.50x
  - Mean coverage (max): 25.30x

================================================================================
Per-Region Analysis
================================================================================

✓ PASS chr1_region1 (chr1:1000000-2000000)
  Size: 1,000,000 bp
  Mean coverage: 15.50x
  Median coverage: 14.00x
  Coverage range: 8-25x
  Bases covered: 985,000/1,000,000 (98.50%)

✗ FAIL chr2_region2 (chr2:3000000-4000000)
  Size: 1,000,000 bp
  Mean coverage: 2.30x
  Median coverage: 1.50x
  Coverage range: 0-8x
  Bases covered: 450,000/1,000,000 (45.00%)
  Reason: Low coverage fraction (45.00% < 90.00%)
  Coverable subregions: 1
    1. chr2:3000000-3500000 (500,000 bp)

================================================================================
Wisecondor/WisecondorX Bin Analysis
================================================================================
Bin size: 200,000 bp

chr1_region1: 5 bins (1,000,000 bp / 200,000 bp)
chr2_region2: 5 bins (1,000,000 bp / 200,000 bp)
```

## 실제 사용 예시 (Real-world Examples)

### 예시 1: 기본 MD 타겟 분석

```bash
# MD_Target_8 영역에 대한 기본 분석
python3 md_coverage_check.py \
    --bam /Work/NIPT/analysis/2511/SAMPLE001/SAMPLE001.proper_paired.bam \
    --bed /Work/NIPT/data/bed/common/MD_Target_8.bed \
    --output /Work/NIPT/analysis/2511/SAMPLE001/coverable_MD_Target_8.bed \
    --report
```

### 예시 2: 높은 품질 요구사항으로 분석

```bash
# 최소 5x 커버리지, 95% 영역 커버 필요
python3 md_coverage_check.py \
    --bam sample.proper_paired.bam \
    --bed MD_Target_8.bed \
    --output coverable_strict.bed \
    --min-coverage 5 \
    --coverage-threshold 0.95 \
    --min-region-size 5000 \
    --report
```

### 예시 3: 여러 타겟에 대한 배치 분석

```bash
#!/bin/bash
SAMPLE_ID="SAMPLE001"
BAM="/Work/NIPT/analysis/2511/${SAMPLE_ID}/${SAMPLE_ID}.proper_paired.bam"
OUT_DIR="/Work/NIPT/analysis/2511/${SAMPLE_ID}/coverage_analysis"

mkdir -p "${OUT_DIR}"

# 각 MD 타겟에 대해 분석
for TARGET in MD_Target_1 MD_Target_4 MD_Target_5 MD_Target_8 MD_Target_15 MD_Target_22; do
    echo "Analyzing ${TARGET}..."
    python3 md_coverage_check.py \
        --bam "${BAM}" \
        --bed "/Work/NIPT/data/bed/common/${TARGET}.bed" \
        --output "${OUT_DIR}/coverable_${TARGET}.bed" \
        --report
done
```

## 결과 해석 (Interpreting Results)

### 분석 가능 (Analyzable)
- ✓ **Fully analyzable**: 영역 전체가 분석 가능
  - 커버리지 threshold 만족
  - 최소 영역 크기 만족
  - Wisecondor/WisecondorX 분석 가능

### 부분 커버 (Partially Covered)
- ⚠ **Partially covered**: 영역의 일부만 분석 가능
  - 커버 가능한 부분 영역이 출력 BED에 포함됨
  - 해당 부분 영역만 사용하여 분석 가능

### 분석 불가 (Not Analyzable)
- ✗ **Not analyzable**: 분석 불가능
  - 커버리지 부족
  - 영역 크기 부족
  - BAM 품질 확인 필요

## 다음 단계 (Next Steps)

1. **리포트 확인**: `*_report.txt` 파일에서 상세 분석 결과 확인
2. **BED 파일 사용**: 출력된 coverable regions BED 파일을 Wisecondor/WisecondorX 분석에 사용
3. **MD 분석 실행**: 
   ```bash
   # 필터링된 영역으로 MD 분석 실행
   python3 md_pipeline.py \
       --sample_id SAMPLE001 \
       --work_dir 2511 \
       --labcode cordlife
   ```

## 문제 해결 (Troubleshooting)

### BAM 파일 인덱스가 없는 경우
도구가 자동으로 인덱스를 생성합니다:
```
WARNING: BAM file is not indexed. Creating index...
```

수동 인덱스 생성:
```bash
samtools index sample.proper_paired.bam
```

### 모든 영역이 분석 불가능한 경우
- BAM 파일 품질 확인
- 시퀀싱 깊이 확인
- 타겟 BED 파일 확인 (올바른 좌표계 사용 여부)

```bash
# BAM 파일 통계 확인
samtools flagstat sample.proper_paired.bam
samtools stats sample.proper_paired.bam
```

### 메모리 부족
큰 BAM 파일의 경우 메모리가 부족할 수 있습니다. 이 경우:
1. 더 작은 영역으로 나누어 분석
2. 더 큰 메모리가 있는 시스템에서 실행

## 기술 세부사항 (Technical Details)

### 커버리지 계산 방법
- pysam의 pileup을 사용하여 각 염기 위치의 read depth 계산
- 각 영역에 대해 평균, 중간값, 최소, 최대 커버리지 계산

### 분석 가능성 판단 기준
영역이 분석 가능하려면 다음 조건을 모두 만족해야 합니다:
1. `coverage_fraction >= coverage_threshold`
2. `region_size >= min_region_size`
3. `mean_coverage >= min_coverage`

### Wisecondor Bin 분석
- Wisecondor/WisecondorX는 지정된 bin 크기로 게놈을 나누어 분석
- 기본 bin 크기: 200kb
- 각 영역이 몇 개의 bin을 포함하는지 계산하여 리포트에 표시

## 저자 (Author)
Ken

## 버전 (Version)
1.0

## 라이선스 (License)
Internal use only - NIPT Pipeline

