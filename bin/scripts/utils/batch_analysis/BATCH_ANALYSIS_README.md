# NIPT Batch Effect Analysis Pipeline

NIPT 샘플들의 배치 효과(Batch Effect)를 종합적으로 분석하는 통합 파이프라인입니다.

## 📋 목차

- [개요](#개요)
- [분석 내용](#분석-내용)
- [사용 방법](#사용-방법)
- [출력 파일](#출력-파일)
- [결과 해석](#결과-해석)

## 🎯 개요

이 파이프라인은 NIPT 샘플들의 배치 효과를 다각도로 분석합니다:

1. **기본 배치 분석**: 배치별 QC 메트릭 분포 및 통계적 비교
2. **PCA 분석**: 주성분 분석을 통한 샘플 패턴 시각화
3. **클러스터링 분석**: K-means를 이용한 샘플 그룹화 및 특성 비교
4. **타임라인 분석**: 시간에 따른 QC 메트릭 변화 추적

## 📊 분석 내용

### 1. 기본 배치 분석
- 배치별 샘플 수 분포
- QC 메트릭 (mapping rate, duplication rate, GC content 등) 통계
- 염색체별 read 분포 및 coverage uniformity
- Fetal fraction 분포
- Kruskal-Wallis 검정을 통한 배치 간 통계적 유의성 검증

### 2. PCA 분석 (3가지)
- **전체 메트릭 PCA**: 모든 QC 메트릭 + 염색체 비율
- **염색체 전용 PCA**: chr13, chr18, chr21 비율만 사용 (성별 효과 제거)
- **QC 전용 PCA**: GC content, mapping rate, duplication rate, mapping quality, coverage만 사용

### 3. K-means 클러스터링
- 샘플을 2개 그룹으로 자동 분류
- 각 그룹의 특성 통계적 비교 (Mann-Whitney U test)
- 배치별 그룹 구성 분석

### 4. 타임라인 분석
- Duplication rate의 시간에 따른 변화
- Coverage의 시간에 따른 변화
- 배치별 평균값 추적
- 프로토콜 변경 시점 감지

## 🚀 사용 방법

### 기본 사용법

```bash
python3 bin/scripts/utils/batch_analysis/run_batch_analysis.py \
  <sample_list.tsv> \
  <analysis_dir> \
  <output_dir> \
  <lab_name>
```

### 예제: UCL 샘플 분석

```bash
cd ~/ken-nipt

python3 bin/scripts/utils/batch_analysis/run_batch_analysis.py \
  data/refs/ucl/reference_make/reference_sample_list_UCL_all.tsv \
  ~/ken-nipt/analysis \
  analysis/batch_analysis_ucl \
  UCL
```

### 예제: Cordlife 샘플 분석

```bash
cd ~/ken-nipt

python3 bin/scripts/utils/batch_analysis/run_batch_analysis.py \
  data/refs/cordlife/reference_make/reference_sample_list_Cordlife_all.tsv \
  ~/ken-nipt/analysis \
  analysis/batch_analysis_cordlife \
  Cordlife
```

## 📁 출력 파일

### 텍스트 리포트
- `batch_effect_summary.txt`: 전체 배치 효과 종합 리포트
- `pca_clustering_report.txt`: 클러스터링 결과 상세 리포트
- `pca_group_analysis.txt`: PCA 그룹 분석 리포트
- `sample_list_with_clusters.tsv`: 클러스터 정보가 추가된 샘플 리스트

### 시각화 그래프

#### 기본 배치 분석 (5개)
1. `01_batch_sample_counts.png`: 배치별 샘플 수
2. `02_qc_by_batch.png`: 배치별 QC 메트릭 box plot
3. `03_qc_trends.png`: 시간에 따른 QC 메트릭 추세
4. `04_ff_by_batch.png`: 배치별 fetal fraction 분포
5. `05_chromosome_distributions.png`: 배치별 염색체 분포

#### PCA 분석 (3개)
6. `06_pca_batch_and_lab.png`: 전체 메트릭 PCA
7. `07_pca_chromosome_only.png`: 염색체 비율 전용 PCA
8. `08_pca_qc_only.png`: QC 메트릭 전용 PCA

#### 클러스터링 분석 (2개)
9. `pca_clustering_groups.png`: K-means 클러스터링 결과
10. `pca_groups_comparison.png`: 그룹 간 QC 메트릭 비교

#### 타임라인 분석 (4개)
11. `duplication_rate_timeline.png`: Duplication rate 타임라인
12. `duplication_rate_timeline_with_means.png`: 배치 평균 포함
13. `coverage_timeline.png`: Coverage 타임라인
14. `coverage_timeline_with_means.png`: 배치 평균 포함

### HTML 리포트
- `batch_analysis_report.html`: 모든 그래프와 분석 결과를 포함한 종합 HTML 리포트

## 📈 결과 해석

### 1. 배치 효과 확인

**Kruskal-Wallis test p-value 확인:**
- p < 0.001: 배치 간 매우 유의미한 차이 (⚠️ 주의 필요)
- p < 0.05: 배치 간 유의미한 차이
- p ≥ 0.05: 배치 간 차이 없음 (✅ 양호)

### 2. PCA 해석

**PCA 산점도에서:**
- 샘플들이 **무작위로 분산**: 배치 효과 적음 ✅
- 샘플들이 **배치별로 군집**: 배치 효과 있음 ⚠️
- **Feature importance** 확인: 어떤 메트릭이 분리에 가장 영향을 미치는지

### 3. 클러스터링 해석

**두 그룹으로 나뉘는 경우:**
- `pca_groups_comparison.png`에서 통계적으로 유의미한 차이가 있는 메트릭 확인
- 주요 분리 원인:
  - Coverage 차이 → Sequencing depth 변경
  - Duplication rate 차이 → Library prep 변경
  - GC content 차이 → 샘플 특성 또는 prep 방법

### 4. 타임라인 해석

**배치 평균선(빨간선) 패턴:**
- **급격한 변화 시점**: 프로토콜 변경 가능성
- **두 그룹 교차**: 프로토콜 전환기
- **안정적 패턴**: 일관된 프로토콜

### 5. Reference 생성 전략

#### Case 1: 배치 효과가 거의 없는 경우 (p > 0.05)
✅ **모든 샘플을 하나의 reference로 사용 가능**

#### Case 2: 두 그룹으로 명확히 나뉘는 경우
⚠️ **두 개의 별도 reference 생성 권장**
- Reference A: Group 1 샘플들
- Reference B: Group 2 샘플들

#### Case 3: 시간에 따라 프로토콜이 변경된 경우
⚠️ **프로토콜별 reference 생성 권장**
- Reference_old: 변경 이전 샘플들
- Reference_new: 변경 이후 샘플들

#### Case 4: 특정 배치에서만 문제가 있는 경우
⚠️ **문제 배치 제외 권장**
- Outlier 배치 필터링 후 reference 생성

## 🔧 개별 스크립트 사용

필요시 각 단계를 개별적으로 실행할 수 있습니다:

### Step 1: 기본 배치 분석
```bash
python3 bin/scripts/utils/batch_analysis/analyze_batch_effect_v2.py \
  --sample-list <sample_list.tsv> \
  --analysis-dir ~/ken-nipt/analysis \
  --output-dir <output_dir> \
  --lab <lab_name>
```

### Step 2: PCA 그룹 분석
```bash
python3 bin/scripts/utils/batch_analysis/analyze_pca_groups.py \
  <sample_list.tsv> <output_dir>
```

### Step 3: K-means 클러스터링
```bash
python3 bin/scripts/utils/batch_analysis/analyze_pca_clustering.py \
  <sample_list.tsv> <output_dir> <lab_name>
```

### Step 4: 타임라인 분석
```bash
python3 bin/scripts/utils/batch_analysis/plot_duplication_timeline.py \
  <output_dir> <lab_name>
```

## 📝 필요 조건

### Python 패키지
- pandas
- numpy
- matplotlib
- seaborn
- scipy
- scikit-learn

### 설치
```bash
pip install pandas numpy matplotlib seaborn scipy scikit-learn
```

## 🐛 문제 해결

### Q: sklearn import 에러
```bash
pip install scikit-learn
```

### Q: 그래프가 생성되지 않음
- 출력 디렉토리 권한 확인
- matplotlib backend 설정 확인

### Q: Coverage/duplication 데이터가 없음
- `10mb.wig.Normalization.txt` 파일 존재 확인
- `qc.txt` 파일 존재 확인

## 📧 문의

문제가 있거나 개선 사항이 있으면 이슈를 등록해주세요.

---

**Last Updated:** 2026-02-02
