# run_parallel_v3.py 사용 가이드

## 개요

`run_parallel_v3.py`는 Microdeletion(MD) 분석을 위한 인공 샘플(Artificial Sample)을 생성하는 스크립트입니다.

### 주요 특징

- **슬롯 관리 (Slot Management)**: 최대 동시 실행 개수를 제한하고, 하나가 완료되면 즉시 다음 샘플 시작
- **배치 생성**: 여러 조합(질환 × FF × 커버리지 × Mom/Fetus 조합)을 자동으로 생성
- **실시간 모니터링**: 주기적으로 상태를 확인하고 진행 상황을 출력
- **상세 로그**: 각 샘플별 독립 로그 파일 및 최종 요약 생성

### v2 대비 개선점

| 항목 | v2 (ProcessPoolExecutor) | v3 (Slot Management) |
|------|-------------------------|---------------------|
| 동시 실행 제어 | 고정된 worker pool | 동적 슬롯 관리 |
| 완료 후 처리 | 모든 worker가 끝날 때까지 대기 | 즉시 다음 태스크 시작 |
| 메모리 사용 | 전체 결과를 메모리에 보관 | 스트리밍 방식으로 처리 |
| 실시간 모니터링 | 제한적 | 상태 업데이트 + 로그 출력 |

---

## 사전 준비

### 1. 필수 파일

#### 입력 TSV 파일
- **mom_list.tsv**: 엄마(임산부) 또는 비임신 여성 샘플 목록
- **female_fetus_list.tsv**: 여아 태아 샘플 목록
- **male_fetus_list.tsv**: 남아 태아 샘플 목록

**TSV 파일 형식:**
```tsv
Work_Dir	Sample_ID	Gender	FF_Method	FF_Value	...	BAM_Path
2508	GNCI25080001	F	seqFF	4.5	...	/path/to/GNCI25080001/GNCI25080001.proper_paired.bam
2508	GNCI25080002	F	seqFF	5.2	...	/path/to/GNCI25080002/GNCI25080002.proper_paired.bam
```

**중요 컬럼:**
- Column 2 (index 1): `Sample_ID` - BAM 파일의 샘플 ID
- Column 5 (index 4): `FF_Value` - Fetal Fraction 값 (%)
- 마지막 컬럼: `BAM_Path` - BAM 파일 경로 (절대경로)

#### BED 파일
Microdeletion 영역을 정의하는 BED 파일

**형식:**
```tsv
chr	start	end	disease_name	...
1	10001	12840259	1p36 deletion syndrome	...
5	123000	456000	5p deletion syndrome	...
```

**중요 컬럼:**
- Column 1-3: 영역 정의 (chr, start, end)
- Column 4: 질환 이름 (disease name)

#### make_artificial.sh
실제 인공 샘플을 생성하는 bash 스크립트 (이미 구현됨)

---

## 사용법

### 기본 명령어

```bash
python3 run_parallel_v3.py \
  --mom_tsv <mom_samples.tsv> \
  --female_tsv <female_fetus.tsv> \
  --male_tsv <male_fetus.tsv> \
  --md_bed <diseases.bed> \
  --script <make_artificial.sh> \
  --output <output_directory> \
  --ff_targets <FF1,FF2,...> \
  --coverages <COV1,COV2,...> \
  --max_workers <N> \
  [options]
```

### 필수 인자

| 인자 | 설명 | 예시 |
|------|------|------|
| `--mom_tsv` | 엄마 샘플 TSV 파일 | `mom_list.tsv` |
| `--female_tsv` | 여아 태아 샘플 TSV 파일 | `female_fetus_list.tsv` |
| `--male_tsv` | 남아 태아 샘플 TSV 파일 | `male_fetus_list.tsv` |
| `--md_bed` | Microdeletion BED 파일 | `test_1p36_only.bed` |

### 선택 인자

| 인자 | 기본값 | 설명 |
|------|--------|------|
| `--script` | `make_artificial.sh` | 샘플 생성 스크립트 경로 |
| `--output` | `test_output` | 출력 디렉토리 |
| `--ff_targets` | `5,10,15` | 타겟 FF 값 (콤마 구분) |
| `--coverages` | `10M` | 타겟 커버리지 (예: `10M`, `15M`) |
| `--max_workers` | `4` | 최대 동시 실행 개수 |
| `--poll_interval` | `5` | 상태 확인 주기 (초) |
| `--n_moms` | 모두 | 사용할 엄마 샘플 개수 제한 |
| `--n_fetuses` | 모두 | 사용할 태아 샘플 개수 제한 (성별당) |
| `--limit_samples` | 모두 | 전체 생성 샘플 개수 제한 |

---

## 사용 예시

### 예시 1: 빠른 테스트 (1개 질환, 최소 조합)

```bash
cd /home/ken/ken-nipt/src/md_test

python3 run_parallel_v3.py \
  --mom_tsv mom_list.tsv \
  --female_tsv female_fetus_list.tsv \
  --male_tsv male_fetus_list.tsv \
  --md_bed test_1p36_only.bed \
  --script make_artificial.sh \
  --output test_output \
  --ff_targets 5 \
  --coverages 10M \
  --n_moms 2 \
  --n_fetuses 1 \
  --max_workers 2 \
  --poll_interval 3
```

**생성 샘플 수:** 2 moms × 1 female × 1 FF × 1 coverage = 2 samples

### 예시 2: 중간 규모 테스트 (여러 FF, 슬롯 4개)

```bash
python3 run_parallel_v3.py \
  --mom_tsv mom_list.tsv \
  --female_tsv female_fetus_list.tsv \
  --male_tsv male_fetus_list.tsv \
  --md_bed test_1p36_only.bed \
  --ff_targets 5,10,15 \
  --coverages 10M \
  --n_moms 3 \
  --n_fetuses 2 \
  --max_workers 4
```

**생성 샘플 수:** 3 moms × 2 fetuses (female+male) × 3 FF × 1 coverage = 18 samples

### 예시 3: 대규모 배치 (모든 조합, 슬롯 8개)

```bash
python3 run_parallel_v3.py \
  --mom_tsv mom_list.tsv \
  --female_tsv female_fetus_list.tsv \
  --male_tsv male_fetus_list.tsv \
  --md_bed all_diseases.bed \
  --ff_targets 5,10,15 \
  --coverages 10M,15M,20M \
  --max_workers 8 \
  --poll_interval 10
```

**생성 샘플 수:** N_diseases × N_moms × (N_female + N_male) × 3 FF × 3 coverages

### 예시 4: 샘플 수 제한 (디버깅용)

```bash
python3 run_parallel_v3.py \
  --mom_tsv mom_list.tsv \
  --female_tsv female_fetus_list.tsv \
  --male_tsv male_fetus_list.tsv \
  --md_bed test_1p36_only.bed \
  --ff_targets 5,10 \
  --coverages 10M \
  --limit_samples 5 \
  --max_workers 2
```

**효과:** 조합이 20개 생성되더라도 처음 5개만 실행

---

## 출력 구조

```
md_test/
├── test_output/                          # 생성된 BAM 파일
│   └── 1p36_deletion_syndrome/          # 질환별 디렉토리
│       ├── FF05_10M/                    # FF와 커버리지별 디렉토리
│       │   ├── 0001_1p36_deletion_syndrome_FF05_10M.bam
│       │   ├── 0001_1p36_deletion_syndrome_FF05_10M.bam.bai
│       │   ├── 0002_1p36_deletion_syndrome_FF05_10M.bam
│       │   └── ...
│       ├── FF10_10M/
│       └── FF15_10M/
│
├── test_logs/                            # 로그 파일
│   └── 1p36_deletion_syndrome/
│       ├── FF05_10M/
│       │   ├── 0001_1p36_deletion_syndrome_FF05_10M.log
│       │   └── ...
│       └── ...
│
└── summary/                              # 요약 파일
    ├── all_samples.tsv                   # 전체 샘플 정보 (TSV)
    └── sample_mix.log                    # 샘플 믹싱 로그 (human-readable)
```

### 요약 파일

#### all_samples.tsv
모든 생성된 샘플의 메타데이터 (TSV 형식)

**컬럼:**
- `idx`: 샘플 인덱스
- `bam_name`: 생성된 BAM 파일명
- `ff_target`: 타겟 FF (%)
- `pairs`: 타겟 페어 수
- `disease`: 질환 이름
- `mom_id`, `fetus_id`: 원본 샘플 ID
- `mom_ff`, `fetus_ff`: 원본 샘플의 FF 값
- `alpha`, `beta`: 믹싱 비율
- `upstream_ratio`, `deletion_ratio`, `downstream_ratio`: 영역별 read 비율 (QC)
- `bam_path`: 생성된 BAM 파일 경로
- `elapsed_sec`: 실행 시간

#### sample_mix.log
사람이 읽기 쉬운 형태의 요약 로그

```
0001.	0001_1p36_deletion_syndrome_FF05_10M.bam
	Mom: GNCI25080001 (FF=4.5%)
	Fetus: GNCI25080163 (FF=12.3%)
	Mixing: α=0.654 (Mom) + β=0.346 (Fetus) → Target FF=5%
	Region Check:
	  Upstream   Ratio: 0.6523
	  Deletion   Ratio: 0.5234 ← Should be lower!
	  Downstream Ratio: 0.6498
```

---

## 실행 로그 예시

```
2025-10-29 12:00:00 [   INFO] | Found 1 disease(s): 1p36 deletion syndrome

2025-10-29 12:00:00 [   INFO] | Loading samples...
2025-10-29 12:00:00 [   INFO] |   Moms: 5
2025-10-29 12:00:00 [   INFO] |   Female fetuses: 3
2025-10-29 12:00:00 [   INFO] |   Male fetuses: 3

2025-10-29 12:00:00 [   INFO] | FF targets: [5, 10, 15]
2025-10-29 12:00:00 [   INFO] | Coverages: ['10M']
2025-10-29 12:00:00 [   INFO] | Max workers: 4
2025-10-29 12:00:00 [   INFO] | Poll interval: 5s

2025-10-29 12:00:00 [   INFO] | Total samples to generate: 90

2025-10-29 12:00:01 [   INFO] | Starting generation with slot management...

2025-10-29 12:00:01 [   INFO] | [0001] Started: 0001_1p36_deletion_syndrome_FF05_10M.bam
2025-10-29 12:00:01 [   INFO] | [0002] Started: 0002_1p36_deletion_syndrome_FF05_10M.bam
2025-10-29 12:00:01 [   INFO] | [0003] Started: 0003_1p36_deletion_syndrome_FF05_10M.bam
2025-10-29 12:00:01 [   INFO] | [0004] Started: 0004_1p36_deletion_syndrome_FF05_10M.bam

2025-10-29 12:00:06 [   INFO] | ============================================================
2025-10-29 12:00:06 [   INFO] | Status: Running=4/4, Completed=0, Failed=0, Pending=86, Total=90
2025-10-29 12:00:06 [   INFO] | ============================================================

2025-10-29 12:05:23 [   INFO] | [0001] ✓ Completed: 0001_1p36_deletion_syndrome_FF05_10M.bam (322.5s)
2025-10-29 12:05:23 [   INFO] | [0005] Started: 0005_1p36_deletion_syndrome_FF05_10M.bam

2025-10-29 12:06:15 [   INFO] | [0002] ✓ Completed: 0002_1p36_deletion_syndrome_FF05_10M.bam (374.2s)
2025-10-29 12:06:15 [   INFO] | [0006] Started: 0006_1p36_deletion_syndrome_FF05_10M.bam

...

2025-10-29 15:30:42 [   INFO] | [0090] ✓ Completed: 0090_1p36_deletion_syndrome_FF15_10M.bam (298.7s)

2025-10-29 15:30:47 [   INFO] | ============================================================
2025-10-29 15:30:47 [   INFO] | Generation Complete!
2025-10-29 15:30:47 [   INFO] |   Total time: 12646.3s
2025-10-29 15:30:47 [   INFO] |   Succeeded: 90
2025-10-29 15:30:47 [   INFO] |   Failed: 0
2025-10-29 15:30:47 [   INFO] |   Total: 90
2025-10-29 15:30:47 [   INFO] | ============================================================

2025-10-29 15:30:47 [   INFO] | 
✓ Summary saved to:
2025-10-29 15:30:47 [   INFO] |   - summary/all_samples.tsv
2025-10-29 15:30:47 [   INFO] |   - summary/sample_mix.log
```

---

## 슬롯 관리 동작 방식

### 상태 전이

```
Init → Running → Completed (성공)
                → Failed    (실패)
```

### 슬롯 관리 알고리즘

1. **초기화**: 모든 태스크를 `Init` 상태로 등록
2. **시작 루프**:
   - 사용 가능한 슬롯 확인: `max_workers - running_count`
   - 슬롯이 있고 `Init` 태스크가 있으면 → 시작 → `Running`으로 변경
3. **모니터링 루프** (poll_interval마다):
   - 실행 중인 프로세스 확인 (`poll()`)
   - 완료된 프로세스 → 로그 저장, 결과 파싱, `Completed`/`Failed` 변경
   - 새로 생긴 빈 슬롯 → 다음 `Init` 태스크 시작
4. **종료 조건**: `Completed + Failed == Total`

### 최대 효율 유지

- 슬롯이 빈 즉시 다음 태스크 시작
- 유휴 시간 최소화
- 리소스 사용 제어 (max_workers 제한)

---

## 문제 해결

### 1. TSV 파일 형식 오류

**증상:**
```
ERROR: No moms or fetuses found.
```

**해결:**
- TSV 파일의 컬럼 개수 확인 (최소 11개 이상)
- 마지막 컬럼이 BAM 파일 경로인지 확인
- BAM 파일이 실제로 존재하는지 확인

### 2. BED 파일 형식 오류

**증상:**
```
ERROR: No diseases found in test.bed
```

**해결:**
- BED 파일의 4번째 컬럼에 질환 이름이 있는지 확인
- `#`으로 시작하는 주석 라인이 아닌지 확인

### 3. make_artificial.sh 실행 실패

**증상:**
```
[0001] ✗ Failed: ... (rc=1)
```

**해결:**
- 로그 파일 확인: `test_logs/{disease}/FFxx_xxM/{sample}.log`
- samtools, bedtools 설치 여부 확인
- BED 파일의 chromosome 이름과 BAM의 reference 이름 일치 확인

### 4. 메모리 부족

**증상:**
시스템이 느려지거나 프로세스가 강제 종료됨

**해결:**
- `--max_workers` 값을 줄이기 (예: 8 → 4 → 2)
- `--coverages` 값을 줄이기 (예: 20M → 10M)

### 5. 디스크 공간 부족

**증상:**
```
OSError: [Errno 28] No space left on device
```

**해결:**
- 출력 디렉토리의 디스크 공간 확인
- 불필요한 BAM 파일 정리
- 샘플 수 제한 (`--limit_samples`)

---

## 팁과 권장사항

### 1. 단계적 실행

1. **소규모 테스트** (`--limit_samples 2`)
2. **중간 규모** (실제 조합의 10%)
3. **전체 실행**

### 2. 적절한 max_workers 설정

- **CPU 코어 수**: 서버 코어 수의 50-75%
- **메모리**: 샘플당 약 2-4GB 필요 (coverage에 따라)
- **디스크 I/O**: SSD인 경우 더 많은 worker 가능

**권장:**
- 8 cores, 32GB RAM → `--max_workers 4`
- 16 cores, 64GB RAM → `--max_workers 8`
- 32 cores, 128GB RAM → `--max_workers 12`

### 3. Poll interval 설정

- 짧은 샘플 (< 5분): `--poll_interval 3`
- 중간 샘플 (5-10분): `--poll_interval 5`
- 긴 샘플 (> 10분): `--poll_interval 10`

### 4. Non-pregnancy women 샘플 사용

Non-pregnancy women (FF=0%) 샘플을 mom_list.tsv에 추가하면 더 정확한 FF 제어 가능:

```tsv
NonPreg	NPWXXX001	F	seqFF	0.0	...	/path/to/NPWXXX001.sorted.bam
```

---

## 관련 스크립트

- `make_artificial.sh`: 실제 샘플 생성 (Mom + Fetus 믹싱 + Deletion 적용)
- `run_parallel_v2.py`: 이전 버전 (ProcessPoolExecutor 방식)
- `quick_test.sh`: 빠른 테스트용 래퍼 스크립트

---

## 참고

- 생성된 인공 샘플은 `run_md_pipeline.sh`로 분석 가능
- BAM 파일 크기: 10M pairs ≈ 1-2GB (압축)
- 예상 생성 시간: 샘플당 약 5-10분 (서버 성능에 따라)

