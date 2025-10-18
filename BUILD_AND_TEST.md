# Docker Build & Test Guide for MD Pipeline

## 📦 Docker 이미지 빌드

### 1. 일반 빌드 (캐시 사용)

```bash
cd /home/ken/ken-nipt

# 기본 빌드
./build_docker.sh

# 또는 직접 명령어
docker build -t nipt_docker_v1.1:latest -f ./docker/Dockerfile .
```

### 2. 전체 재빌드 (캐시 없이)

```bash
# Dockerfile에서 CACHE_BREAK 변경
sed -i 's/ARG CACHE_BREAK=1/ARG CACHE_BREAK=2/' docker/Dockerfile

# 빌드
./build_docker.sh

# 원복
sed -i 's/ARG CACHE_BREAK=2/ARG CACHE_BREAK=1/' docker/Dockerfile
```

### 3. 빌드 확인

```bash
# 이미지 확인
docker images | grep nipt_docker

# 이미지 상세 정보
docker inspect nipt_docker_v1.1:latest

# md_pipeline.py 포함 확인
docker run --rm --entrypoint ls nipt_docker_v1.1:latest -la /Work/NIPT/bin/scripts/md_pipeline.py
```

---

## 🧪 테스트

### Test 1: Docker 컨테이너 내부에서 직접 실행

```bash
# 컨테이너 진입
docker run --rm -it \
  -v /home/ken/ken-nipt/analysis:/Work/NIPT/analysis \
  -v /home/ken/ken-nipt/output:/Work/NIPT/output \
  -v /home/ken/ken-nipt/data:/Work/NIPT/data \
  -v /home/ken/ken-nipt/config:/Work/NIPT/config \
  --entrypoint bash \
  nipt_docker_v1.1:latest

# 컨테이너 내부에서
python3 /Work/NIPT/bin/scripts/md_pipeline.py --help
```

### Test 2: run_md_pipeline.sh로 실제 샘플 테스트

```bash
cd /home/ken/ken-nipt/src

# 기존 샘플 테스트 (proper_paired.bam 있는 샘플)
./run_md_pipeline.sh \
  -s GNCI25080163 \
  -l cordlife \
  -root /home/ken/ken-nipt \
  -work 2508

# Force 재실행
./run_md_pipeline.sh \
  -s GNCI25080163 \
  -l cordlife \
  -root /home/ken/ken-nipt \
  -work 2508 \
  -f

# 결과 확인
ls -lh /home/ken/ken-nipt/analysis/2508/GNCI25080163/Output_WC/orig/
ls -lh /home/ken/ken-nipt/analysis/2508/GNCI25080163/Output_WCX/orig/
ls -lh /home/ken/ken-nipt/analysis/2508/GNCI25080163/plots/
```

### Test 3: Artificial Sample 테스트

```bash
cd /home/ken/ken-nipt/src

# 1. Artificial sample을 analysis 디렉토리로 복사
SAMPLE_ID="0001_1p36_FF05_10M"
WORK_DIR="md_test_artificial"

mkdir -p /home/ken/ken-nipt/analysis/$WORK_DIR/$SAMPLE_ID

# artificial sample 복사
cp /home/ken/ken-nipt/src/md_test/test_output/1p36/FF05_10M/0001_1p36_FF05_10M.bam \
   /home/ken/ken-nipt/analysis/$WORK_DIR/$SAMPLE_ID/${SAMPLE_ID}.proper_paired.bam

cp /home/ken/ken-nipt/src/md_test/test_output/1p36/FF05_10M/0001_1p36_FF05_10M.bam.bai \
   /home/ken/ken-nipt/analysis/$WORK_DIR/$SAMPLE_ID/${SAMPLE_ID}.proper_paired.bam.bai

# 2. MD 분석 실행 (성별 지정)
./run_md_pipeline.sh \
  -s $SAMPLE_ID \
  -l cordlife \
  -root /home/ken/ken-nipt \
  -work $WORK_DIR \
  --fetal_gender M

# 3. 결과 확인
cat /home/ken/ken-nipt/analysis/$WORK_DIR/$SAMPLE_ID/Output_WCX/orig/${SAMPLE_ID}.wcx.orig_aberrations.bed
cat /home/ken/ken-nipt/analysis/$WORK_DIR/$SAMPLE_ID/${SAMPLE_ID}_md_analysis.log
```

### Test 4: sorted.bam 자동 변환 테스트

```bash
# sorted.bam만 있는 샘플 테스트 (proper_paired.bam 생성 확인)
SAMPLE_ID="test_sorted_only"
WORK_DIR="md_test_sorted"

mkdir -p /home/ken/ken-nipt/analysis/$WORK_DIR/$SAMPLE_ID

# 기존 샘플의 sorted.bam 복사
cp /path/to/some/sample.sorted.bam \
   /home/ken/ken-nipt/analysis/$WORK_DIR/$SAMPLE_ID/${SAMPLE_ID}.sorted.bam

# MD 분석 실행 (자동으로 proper_paired.bam 생성)
./run_md_pipeline.sh \
  -s $SAMPLE_ID \
  -l cordlife \
  -root /home/ken/ken-nipt \
  -work $WORK_DIR

# proper_paired.bam이 생성되었는지 확인
ls -lh /home/ken/ken-nipt/analysis/$WORK_DIR/$SAMPLE_ID/*.bam
```

---

## ✅ 체크리스트

### 빌드 단계
- [ ] Dockerfile 수정됨 (md_pipeline.py COPY 추가)
- [ ] Docker 이미지 빌드 성공
- [ ] md_pipeline.py가 이미지에 포함됨

### 테스트 단계
- [ ] `--help` 옵션 작동
- [ ] 기존 샘플로 MD 분석 성공
- [ ] WC NPZ 파일 생성 확인
- [ ] WCX NPZ 파일 생성 확인
- [ ] Plot 파일 생성 확인
- [ ] MD detection TSV 파일 생성 확인
- [ ] sorted.bam → proper_paired.bam 자동 변환 작동
- [ ] Artificial sample 분석 성공
- [ ] 성별별 WCX reference 올바르게 사용됨

---

## 🐛 문제 해결

### 1. "md_pipeline.py not found" 에러

**원인**: Docker 이미지가 오래된 버전  
**해결**:
```bash
# 이미지 재빌드
./build_docker.sh

# 확인
docker run --rm --entrypoint ls nipt_docker_v1.1:latest -la /Work/NIPT/bin/scripts/md_pipeline.py
```

### 2. "No module named 'process_md_result'" 에러

**원인**: sys.path 문제  
**해결**: md_pipeline.py의 sys.path.append 확인
```python
sys.path.append("/Work/NIPT/bin")  # 이 줄이 있어야 함
```

### 3. Reference 파일 없음

**원인**: data 디렉토리 마운트 안 됨  
**해결**:
```bash
# run_md_pipeline.sh에서 마운트 확인
-v "$HOST_DATA_DIR:/Work/NIPT/data"

# Reference 파일 확인
ls -lh /home/ken/ken-nipt/data/refs/cordlife/WC/
ls -lh /home/ken/ken-nipt/data/refs/cordlife/WCX/
```

### 4. 성별별 WCX reference 문제

**원인**: 성별 감지 실패 또는 reference 파일 없음  
**해결**:
```bash
# 성별 파일 확인
cat /home/ken/ken-nipt/analysis/2508/SAMPLE_ID/Output_FF/SAMPLE_ID.gender.txt

# Reference 확인 (orig_M_, orig_F_ 모두 있어야 함)
ls -lh /home/ken/ken-nipt/data/refs/cordlife/WCX/orig_*_200k_proper_paired.npz

# 수동으로 성별 지정
./run_md_pipeline.sh -s SAMPLE_ID -l cordlife -root /home/ken/ken-nipt -work 2508 --fetal_gender M
```

### 5. Permission denied

**원인**: User/Group ID 불일치  
**해결**:
```bash
# run_md_pipeline.sh에서 자동 처리됨
--user "${USER_UID}:${USER_GID}"

# 수동으로 권한 수정
sudo chown -R ken:ken /home/ken/ken-nipt/analysis
```

### 6. Output_WC/Output_WC 경로 중복

**원인**: 이전 버전의 md_pipeline.py  
**해결**:
```bash
# 최신 코드로 재빌드 필요
cd /home/ken/ken-nipt
./build_docker.sh
```

---

## 📊 예상 결과

### 성공 시 출력 예시

```
=========================================
MD Pipeline Started for Sample: GNCI25080071
Work Directory: 2508
Lab Code: cordlife
=========================================

=== Step 1: Finding BAM File ===
Found proper_paired.bam: /Work/NIPT/analysis/2508/GNCI25080071/GNCI25080071.proper_paired.bam

=== Step 2: Using existing proper_paired.bam ===

=== Step 2.5: Determining Fetal Gender ===
Reading gender from: /Work/NIPT/analysis/2508/GNCI25080071/Output_FF/GNCI25080071.gender.txt
Gender from gd_2: XX -> F
Using fetal gender from gender.txt: Female

=== Step 3: Running Wisecondor ===
=== Starting Wisecondor Analysis ===
WC Step 1: Converting BAM to NPZ
WC NPZ created: /Work/NIPT/analysis/2508/GNCI25080071/Output_WC/GNCI25080071.wc.proper_paired.npz
WC Step 2: Running prediction
WC Step 3: Generating plot
Wisecondor completed: /Work/NIPT/analysis/2508/GNCI25080071/Output_WC/orig/GNCI25080071.wc.orig.out.npz
Plot saved: /Work/NIPT/analysis/2508/GNCI25080071/Output_WC/orig/GNCI25080071.wc.orig_z.png

=== Step 4: Running WisecondorX ===
=== Starting WisecondorX Analysis ===
WCX Step 1: Converting BAM to NPZ
WCX NPZ created: /Work/NIPT/analysis/2508/GNCI25080071/Output_WCX/GNCI25080071.wcx.proper_paired.npz
Using WCX reference: /Work/NIPT/data/refs/cordlife/WCX/orig_F_200k_proper_paired.npz
WCX Step 2: Running prediction
WisecondorX completed: /Work/NIPT/analysis/2508/GNCI25080071/Output_WCX/orig/GNCI25080071.wcx.orig_aberrations.bed

=== Step 5: Running MD Detection ===
=== Starting MD Detection Pipeline ===
MD detection completed successfully

=========================================
MD Pipeline Completed Successfully!
Results saved in: /Work/NIPT/analysis/2508/GNCI25080071
Plots saved in: /Work/NIPT/analysis/2508/GNCI25080071/plots
=========================================
```

### 생성 파일 목록

```bash
analysis/2508/GNCI25080071/
├── GNCI25080071.proper_paired.bam
├── GNCI25080071.proper_paired.bam.bai
├── Output_WC/
│   ├── GNCI25080071.wc.proper_paired.npz        ← Input NPZ
│   └── orig/
│       ├── GNCI25080071.wc.orig.out.npz         ← Result NPZ
│       └── GNCI25080071.wc.orig_z.png           ← Plot
├── Output_WCX/
│   ├── GNCI25080071.wcx.proper_paired.npz       ← Input NPZ
│   └── orig/
│       ├── GNCI25080071.wcx.orig_aberrations.bed
│       ├── GNCI25080071.wcx.orig_bins.bed
│       ├── GNCI25080071.wcx.orig_segments.bed
│       ├── GNCI25080071.wcx.orig_statistics.txt
│       └── GNCI25080071.wcx.orig.plots/         ← WCX plots
│           └── chr*.png
├── GNCI25080071_md_analysis.log                 ← Analysis log
└── GNCI25080071.md_pipeline_completed.marker    ← Completion marker
```

---

## 🚀 다음 단계

빌드 및 테스트 성공 후:

### 1. 대량 Artificial Sample 분석

```bash
cd /home/ken/ken-nipt/src/md_test

# Artificial samples 생성
python3 run_parallel_v2.py \
  --mom_tsv mom_samples.tsv \
  --female_tsv female_fetus_samples.tsv \
  --male_tsv male_fetus_samples.tsv \
  --md_bed test_1p36_only.bed \
  --ff_targets 5,10,15 \
  --coverages 10M \
  --workers 4 \
  --limit_samples 10

# 생성된 샘플들에 대해 MD 분석 실행
for bam in test_output/*/FF*/*.bam; do
  sample_id=$(basename "$bam" .bam)
  ./run_md_pipeline.sh \
    -s "$sample_id" \
    -l cordlife \
    -root /home/ken/ken-nipt \
    -work md_test_artificial \
    --fetal_gender M
done
```

### 2. 결과 분석

- Detection rate 계산
- FF vs Detection rate plot
- Coverage vs Detection rate plot

### 3. 프로덕션 배포

- 검증 완료 후 이미지 태깅
- Docker Hub push (선택)

---

## 📝 빌드 명령어 요약

```bash
# 1. Docker 이미지 빌드
cd /home/ken/ken-nipt
./build_docker.sh

# 2. 빠른 테스트
docker run --rm --entrypoint python3 nipt_docker_v1.1:latest \
  /Work/NIPT/bin/scripts/md_pipeline.py --help

# 3. 실제 샘플 테스트
cd src
./run_md_pipeline.sh -s GNCI25080071 -l cordlife -root /home/ken/ken-nipt -work 2508 -f
```

---

## 🔄 업데이트 히스토리

### v1.1 (2025-10-16)
- ✅ MD-only pipeline 추가 (`md_pipeline.py`)
- ✅ 성별별 WCX reference 지원
- ✅ Output_FF/gender.txt에서 gd_2 성별 읽기
- ✅ `--fetal_gender` 옵션 추가 (artificial sample용)
- ✅ Output_WC/orig, Output_WCX/orig 구조 사용
- ✅ WC plot 별도 생성

### v1.0 (Initial)
- Full NIPT pipeline

끝! 🎉

