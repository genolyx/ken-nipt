# Dashboard Quick Start Guide

## 현재 데이터 위치

### 샘플 디렉토리
```
/data/md_validation/analysis_result/
├── 1p36/       (2100 samples)
├── 2q33/       (2100 samples)
├── CDC/        (2100 samples)
├── DGS/        (2100 samples)
├── Jacobsen/   (2100 samples)
├── PWS/        (2100 samples)
├── WBS/        (2100 samples)
└── WHS/        (2100 samples)
```

### Z-score TSV 파일
```
/home/ken/ken-nipt/analysis/md_validation/zscore/
├── zscore_extraction_1p36.tsv
├── zscore_extraction_2q33.tsv
├── zscore_extraction_CDC.tsv
├── zscore_extraction_DGS.tsv
├── zscore_extraction_Jacobsen.tsv
├── zscore_extraction_PWS.tsv
├── zscore_extraction_WBS.tsv
└── zscore_extraction_WHS.tsv
```

## 실행 방법

### 1. 패키지 설치

```bash
pip install pandas numpy plotly dash
```

### 2. Dashboard 시작

```bash
cd /home/ken/ken-nipt/src/md_test
python md_dashboard.py --port 8050
```

### 3. 브라우저에서 접속

```
http://localhost:8050
```

## Dashboard 사용법

### Phase 1: 최적 Z-score Threshold 결정

#### Step 1: 데이터 로드
1. **Z-score Extraction TSV** 입력:
   ```
   /home/ken/ken-nipt/analysis/md_validation/zscore/zscore_extraction_1p36.tsv
   ```

2. **Sample Directory** 입력:
   ```
   /data/md_validation/analysis_result/1p36
   ```

3. **Load Data** 버튼 클릭

#### Step 2: ROC 분석
1. **Method** 선택: `WCX_fetus` (권장)
2. **Fetal Fraction** 선택: `10` (%)
3. **Min Detect Length** 선택: `1 Mb` (권장)
4. **Analyze ROC** 버튼 클릭 (10-20초 소요)

#### Step 3: 결과 해석
- **ROC Curve**: TPR vs FPR 곡선
- **Sensitivity vs Threshold**: Z-score에 따른 민감도 변화
- **PPV vs Threshold**: Z-score에 따른 PPV 변화
- **Sensitivity-PPV Trade-off**: 최적 지점 선택용

하단에 추천 threshold가 표시됩니다:
```
PPV ≥ 90%: z = -7.5, Sensitivity = 85.2%, PPV = 90.4%
PPV ≥ 80%: z = -6.2, Sensitivity = 91.3%, PPV = 83.1%
PPV ≥ 70%: z = -5.1, Sensitivity = 95.7%, PPV = 72.8%
```

예: `z = -7.5`를 Phase 2에서 사용

### Phase 2: 성능 시각화

#### Step 1: "Phase 2: Performance Analysis" 탭 클릭

#### Step 2: 파라미터 설정
1. **Method** 선택: `wcx_fetus`
2. **Z-score Threshold** 입력: `-7.5` (Phase 1에서 결정한 값)
3. **Min Length** 입력: `1.0` (Mb)
4. **Calculate Performance** 버튼 클릭 (30-60초 소요)

#### Step 3: 결과 확인
- **Heatmap**: FF × Deletion Length별 Sensitivity
  - 가로축: FF (5%, 10%, 15%)
  - 세로축: Deletion Length (1, 3, 5, 7, 10 Mb)
  - 색상: 초록색(높음) ~ 빨강색(낮음)
  
- **Line Plot**: Deletion Length에 따른 Sensitivity 변화
  - 3개 라인: FF 5%, 10%, 15%
  - X축: Deletion Length
  - Y축: Sensitivity

## 8개 질환 분석하기

각 질환별로 위 과정 반복:

| 질환 | TSV 파일 | 샘플 디렉토리 |
|-----|---------|-------------|
| 1p36 | `/home/ken/ken-nipt/analysis/md_validation/zscore/zscore_extraction_1p36.tsv` | `/data/md_validation/analysis_result/1p36` |
| 2q33 | `/home/ken/ken-nipt/analysis/md_validation/zscore/zscore_extraction_2q33.tsv` | `/data/md_validation/analysis_result/2q33` |
| CDC | `/home/ken/ken-nipt/analysis/md_validation/zscore/zscore_extraction_CDC.tsv` | `/data/md_validation/analysis_result/CDC` |
| DGS | `/home/ken/ken-nipt/analysis/md_validation/zscore/zscore_extraction_DGS.tsv` | `/data/md_validation/analysis_result/DGS` |
| Jacobsen | `/home/ken/ken-nipt/analysis/md_validation/zscore/zscore_extraction_Jacobsen.tsv` | `/data/md_validation/analysis_result/Jacobsen` |
| PWS | `/home/ken/ken-nipt/analysis/md_validation/zscore/zscore_extraction_PWS.tsv` | `/data/md_validation/analysis_result/PWS` |
| WBS | `/home/ken/ken-nipt/analysis/md_validation/zscore/zscore_extraction_WBS.tsv` | `/data/md_validation/analysis_result/WBS` |
| WHS | `/home/ken/ken-nipt/analysis/md_validation/zscore/zscore_extraction_WHS.tsv` | `/data/md_validation/analysis_result/WHS` |

## Detection Mode 설명

### Individual Methods
- **WC_orig** / **wc_orig**: Wisecondor original mode
- **WC_fetus** / **wc_fetus**: Wisecondor fetus mode
- **WCX_orig** / **wcx_orig**: WisecondorX original mode
- **WCX_fetus** / **wcx_fetus**: WisecondorX fetus mode ⭐ 권장

### Aggregated Methods (OR logic)
- **orig**: WC_orig OR WCX_orig (원본 방식 중 하나라도 검출)
- **fetus**: WC_fetus OR WCX_fetus (Fetus 방식 중 하나라도 검출)
- **any**: 4가지 방법 중 하나라도 검출

## 팁

### 권장 설정
- **Method**: `wcx_fetus` (일반적으로 가장 좋은 성능)
- **FF**: `10%` (balanced)
- **Min Detect Length**: `1 Mb` (임상적으로 의미있는 크기)

### PPV 목표
- **PPV ≥ 90%**: 확진 검사용 (False Positive 최소화)
- **PPV ≥ 80%**: 스크리닝 검사용 (균형잡힌 설정)
- **PPV ≥ 70%**: 고민감도 검사용 (False Negative 최소화)

### 예상 성능
- **FF 5%**: 낮은 민감도, 특히 작은 deletion (<3 Mb)
- **FF 10%**: 중간 민감도, 균형잡힌 성능
- **FF 15%**: 높은 민감도, 작은 deletion도 잘 검출

## 문제 해결

### Dashboard가 시작되지 않을 때
```bash
# 패키지 재설치
pip install --upgrade pandas numpy plotly dash

# 포트 변경
python md_dashboard.py --port 8051
```

### 데이터 로딩 실패
- 파일 경로 확인 (절대 경로 사용 권장)
- TSV 파일 형식 확인 (탭 구분자)
- 샘플 디렉토리 내 results/ 폴더 존재 확인

### ROC 분석이 느릴 때
- 정상입니다 (2100개 샘플 × 96개 threshold = ~200,000회 계산)
- 첫 실행: 10-30초
- Min Detect Length를 높이면 약간 빠름

### Performance 계산이 느릴 때
- 정상입니다 (2100개 샘플 분석)
- 예상 시간: 30-60초
- 진행 상황은 터미널 로그에서 확인 가능

## 원격 접속

다른 컴퓨터에서 접속하려면:

```bash
# 서버에서
python md_dashboard.py --host 0.0.0.0 --port 8050

# 클라이언트 브라우저에서
http://192.168.x.x:8050  (서버 IP 주소)
```

## 문의

문제가 있으면 bioinformatics 팀에 문의하세요.


