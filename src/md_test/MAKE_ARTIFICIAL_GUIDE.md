# make_artificial.sh 사용 방법 요약

## 검토 결과

### 1. make_artificial.sh 동작 방식
- **입력**: Mom BAM, Fetus BAM, FF map, BED file, FF target, pairs target
- **과정**:
  1. Mom과 Fetus를 FF target에 맞게 mixing
  2. Fetus에서 deletion 영역의 reads를 ~50% 제거 (heterozygous loss)
  3. 최종 pairs 수에 맞게 downsampling
- **출력**: Mixed BAM 파일 + deletion 검증 정보

### 2. 실행 방법

#### 방법 A: run_parallel_v3.py 사용 (권장)
```bash
cd /home/ken/ken-nipt/src/md_test

# 각 deletion 크기별로 실행
python3 run_parallel_v3.py \
  --mom_tsv non_preg_mom_list.tsv \
  --female_tsv female_fetus_list.tsv \
  --male_tsv male_fetus_list.tsv \
  --md_bed test_1p36_10Mb.bed \
  --script make_artificial.sh \
  --output test_output \
  --ff_targets 5,10,15 \
  --coverages 10M,15M \
  --max_workers 4
```

#### 방법 B: 모든 크기 자동 실행
```bash
cd /home/ken/ken-nipt/src/md_test

# 모든 deletion 크기 자동 실행
./run_all_deletion_sizes.sh
```

#### 방법 C: 개별 샘플 직접 생성 (테스트용)
```bash
# FF map 생성 (한 번만)
cat > ff_map.tsv <<EOF
GNCI25100169	0.00
GNCI25080181	27.8
EOF

# 개별 샘플 생성
./make_artificial.sh \
  --mom_bam /home/ken/ken-nipt/analysis/2510/GNCI25100169/GNCI25100169.proper_paired.bam \
  --fetus_bam /home/ken/ken-nipt/analysis/2508/GNCI25080181/GNCI25080181.proper_paired.bam \
  --ff_map ff_map.tsv \
  --md_bed test_1p36_10Mb.bed \
  --ff_target 10 \
  --pairs 10000000 \
  --outdir test_output/test_sample
```

### 3. 생성 조합

**Mom 샘플**: non_preg_mom_list.tsv (5개, FF=0%)
**Fetus 샘플**: female_fetus_list.tsv + male_fetus_list.tsv
**FF**: 5%, 10%, 15% (3개)
**Coverage**: 10M, 15M pairs (2개)
**Deletion 크기**: 1, 3, 5, 7, 10, 15, 20 Mb (7개)

**예상 생성 샘플 수**:
- 5 moms × (N female + M male fetuses) × 3 FF × 2 coverage × 7 sizes
- 예: 5 × 10 × 3 × 2 × 7 = **2,100 샘플** (fetuses에 따라 달라짐)

### 4. 생성된 파일 구조

```
test_output/
├── 1p36_deletion_syndrome/
│   ├── FF05_10M/
│   │   ├── 0001_1p36_deletion_syndrome_FF05_10M.bam
│   │   ├── 0001_1p36_deletion_syndrome_FF05_10M.bam.bai
│   │   └── ...
│   ├── FF05_15M/
│   ├── FF10_10M/
│   └── ...
└── ...

test_logs/
└── 1p36_deletion_syndrome/
    └── FF05_10M/
        └── 0001_1p36_deletion_syndrome_FF05_10M.log

summary/
├── all_samples.tsv
└── sample_mix.log
```

### 5. 주의사항

1. **FF map 파일**: run_parallel_v3.py가 자동 생성 (각 샘플마다 임시 파일)
2. **Non-pregnant 샘플**: FF=0%이므로 더 정확한 FF 제어 가능
3. **Deletion 크기**: 각 크기별 BED 파일 필요 (이미 생성됨)
4. **메모리**: 샘플당 약 2-4GB 필요
5. **실행 시간**: 샘플당 약 5-10분 (서버 성능에 따라)

### 6. 테스트 실행 권장사항

1. **소규모 테스트** (1개 deletion 크기, 최소 조합):
```bash
python3 run_parallel_v3.py \
  --mom_tsv non_preg_mom_list.tsv \
  --female_tsv female_fetus_list.tsv \
  --male_tsv male_fetus_list.tsv \
  --md_bed test_1p36_10Mb.bed \
  --ff_targets 10 \
  --coverages 10M \
  --n_moms 1 \
  --n_fetuses 1 \
  --max_workers 2
```

2. **전체 실행** (모든 조합):
```bash
./run_all_deletion_sizes.sh
```

### 7. 생성된 파일 확인

```bash
# 생성된 샘플 수 확인
find test_output -name "*.bam" | wc -l

# 요약 파일 확인
cat summary/all_samples.tsv | head -20

# 특정 샘플 로그 확인
cat test_logs/1p36_deletion_syndrome/FF10_10M/0001_*.log
```

## 다음 단계

1. ✅ BED 파일 생성 완료 (7개 크기)
2. ✅ 실행 스크립트 작성 완료
3. ⏳ 소규모 테스트 실행
4. ⏳ 전체 배치 실행

