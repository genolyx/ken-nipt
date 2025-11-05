# MD Artificial Sample 생성 가이드

## 현재 상황

### 샘플 파일
- `mom_list.tsv`: 일반 모체 샘플 (FF > 0%)
- `non_preg_mom_list.tsv`: Non-pregnant 샘플 5개 (FF=0%, eff_reads ~15M)
- `female_fetus_list.tsv`: 여아 태아 샘플
- `male_fetus_list.tsv`: 남아 태아 샘플

### 생성 조건
- **Read 수**: 10M, 15M pairs
- **FF**: 5%, 10%, 15%
- **Deletion 크기**: 1, 3, 5, 7, 10, 15, 20 Mb
- **질환**: 1p36 deletion syndrome (test_1p36_only.bed)

## 주요 발견사항

### 1. make_artificial.sh 동작 방식
- BED 파일 하나당 하나의 deletion 영역 사용
- 다양한 크기의 deletion 테스트를 위해서는 **각 크기별 BED 파일**이 필요
- 현재 `test_1p36_only.bed`는 약 12.8Mb 크기

### 2. Non-pregnant 샘플 활용
- FF=0%이므로 더 정확한 FF 제어 가능
- eff_reads ~15M이므로 15M pairs 생성에 적합
- 10M pairs 생성 시에도 충분한 reads 확보 가능

### 3. run_parallel_v3.py 사용 가능
- 이미 병렬 처리 지원
- 여러 BED 파일을 동시에 처리 가능
- 각 조합별로 자동 생성

## 생성 전략

### 방법 1: 각 크기별 BED 파일 생성 후 run_parallel_v3.py 사용

```bash
# 1. 각 크기별 BED 파일 생성
# 1p36 영역의 시작 위치를 기준으로 여러 크기 생성
# 예: chr1:10001을 시작점으로 1Mb, 3Mb, 5Mb, 7Mb, 10Mb, 15Mb, 20Mb

# 2. run_parallel_v3.py 실행
python3 run_parallel_v3.py \
  --mom_tsv non_preg_mom_list.tsv \
  --female_tsv female_fetus_list.tsv \
  --male_tsv male_fetus_list.tsv \
  --md_bed test_1p36_deletions.bed \  # 여러 deletion 크기 포함
  --ff_targets 5,10,15 \
  --coverages 10M,15M \
  --max_workers 4
```

### 방법 2: 스크립트로 각 크기별 BED 파일 생성 후 순차 실행

각 deletion 크기별로 별도 실행:
```bash
# 각 크기별로 BED 파일 생성
for size in 1 3 5 7 10 15 20; do
  # BED 파일 생성
  python3 create_bed_by_size.py --size ${size}Mb --output test_1p36_${size}Mb.bed
  
  # run_parallel_v3.py 실행
  python3 run_parallel_v3.py \
    --mom_tsv non_preg_mom_list.tsv \
    --female_tsv female_fetus_list.tsv \
    --male_tsv male_fetus_list.tsv \
    --md_bed test_1p36_${size}Mb.bed \
    --ff_targets 5,10,15 \
    --coverages 10M,15M \
    --max_workers 4
done
```

## 필요한 작업

### 1. 각 크기별 BED 파일 생성 스크립트 작성
- 1p36 영역 (chr1:10001-12840259) 기준
- 각 크기(1, 3, 5, 7, 10, 15, 20 Mb)별로 BED 파일 생성
- 시작 위치는 동일하게 유지하거나 중앙 정렬

### 2. FF map 파일 생성
- 모든 샘플의 FF 값을 포함한 ff_map.tsv 생성
- non_preg_mom_list.tsv + mom_list.tsv + fetus_list.tsv 통합

### 3. 생성 조합 확인
- Mom: non_preg (5개) + 일반 mom (선택)
- Fetus: female + male
- FF: 5%, 10%, 15% (3개)
- Coverage: 10M, 15M (2개)
- Deletion size: 1, 3, 5, 7, 10, 15, 20 Mb (7개)

**총 조합 수**: (5 moms + 선택) × (N female + M male fetuses) × 3 FF × 2 coverage × 7 sizes

## 다음 단계

1. 각 크기별 BED 파일 생성 스크립트 작성
2. FF map 통합 파일 생성
3. 테스트 실행 (소규모)
4. 전체 배치 실행

