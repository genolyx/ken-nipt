#!/usr/bin/env python3
"""
염색체 파일에서 임계값 기반 샘플 필터링 스크립트
- orig, fetus, mom 디렉토리의 chr1.txt~chr22.txt 파일 처리
- threshold 파일의 UR_min, Z_min 보다 낮은 값을 가진 샘플 삭제
"""

import pandas as pd
import os
import sys
import logging

# 로거 설정
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def filter_chr_samples_by_threshold(chr_dir, threshold_file):
    """
    임계값을 기준으로 염색체 파일의 샘플들을 필터링
    
    중요: type이 'N'인 Normal 샘플은 임계값과 관계없이 모두 유지
    임계값 필터링은 P 샘플이나 기타 타입 샘플에만 적용
    
    Parameters:
    - chr_dir: chr1.txt~chr22.txt 파일들이 있는 디렉토리
    - threshold_file: 임계값 파일 경로 (orig_thresholds.tsv 등)
    
    Returns:
    - dict: 처리 결과 요약
    """
    
    results = {
        'processed_files': 0,
        'total_samples_before': 0,
        'total_samples_after': 0,
        'removed_samples': 0,
        'error_files': []
    }
    
    try:
        print(f"\n=== 염색체 파일 임계값 필터링 시작 ===")
        print(f"디렉토리: {chr_dir}")
        print(f"임계값 파일: {threshold_file}")
        print(f"필터링 정책: N 타입 샘플은 무조건 유지, 다른 타입만 임계값 적용")
        
        # 디렉토리 존재 확인
        if not os.path.exists(chr_dir):
            print(f"ERROR: 디렉토리가 존재하지 않습니다: {chr_dir}")
            return results
        
        # 임계값 파일 로드
        if not os.path.exists(threshold_file):
            print(f"ERROR: 임계값 파일이 존재하지 않습니다: {threshold_file}")
            return results
        
        threshold_df = pd.read_csv(threshold_file, sep='\t')
        print(f"임계값 로드 완료: {len(threshold_df)}개 염색체")
        print(f"임계값 컬럼: {list(threshold_df.columns)}")
        
        # 각 염색체 파일 처리
        for i in range(1, 23):
            chr_name = f'chr{i}'
            chr_file = os.path.join(chr_dir, f'{chr_name}.txt')
            
            if not os.path.exists(chr_file):
                print(f"  {chr_name}: 파일 없음")
                continue
            
            try:
                # 해당 염색체의 임계값 가져오기
                threshold_row = threshold_df[threshold_df['chr'] == chr_name]
                
                if threshold_row.empty:
                    print(f"  {chr_name}: 임계값 정보 없음, 건너뜀")
                    continue
                
                th_data = threshold_row.iloc[0]
                ur_min = th_data.get('UR_min', None)
                z_min = th_data.get('Z_min', None)
                
                if pd.isna(ur_min) or pd.isna(z_min):
                    print(f"  {chr_name}: UR_min 또는 Z_min 값이 없음, 건너뜀")
                    print(f"    UR_min: {ur_min}, Z_min: {z_min}")
                    continue
                
                print(f"  {chr_name}: 임계값 UR_min={ur_min}, Z_min={z_min}")
                
                # 염색체 파일 로드
                print(f"    파일 로드 중: {chr_file}")
                df = pd.read_csv(chr_file, sep='\t')
                original_count = len(df)
                
                print(f"    로드된 데이터: {original_count}행, 컬럼: {list(df.columns)}")
                
                if df.empty:
                    print(f"    빈 파일, 건너뜀")
                    continue
                
                # 첫 몇 행 출력
                if len(df) > 0:
                    print(f"    데이터 샘플:")
                    print(df.head(3).to_string(index=False))
                
                # 필수 컬럼 확인
                if 'UR' not in df.columns or 'Z' not in df.columns:
                    print(f"    ERROR: 필수 컬럼(UR, Z) 없음")
                    print(f"    사용 가능한 컬럼: {list(df.columns)}")
                    results['error_files'].append(f"{chr_name} (missing columns)")
                    continue
                
                # 데이터 타입 확인
                print(f"    UR 컬럼 타입: {df['UR'].dtype}")
                print(f"    Z 컬럼 타입: {df['Z'].dtype}")
                
                # UR, Z 컬럼에 숫자가 아닌 값이 있는지 확인
                ur_numeric = pd.to_numeric(df['UR'], errors='coerce')
                z_numeric = pd.to_numeric(df['Z'], errors='coerce')
                
                ur_na_count = ur_numeric.isna().sum()
                z_na_count = z_numeric.isna().sum()
                
                if ur_na_count > 0:
                    print(f"    WARNING: UR 컬럼에 숫자가 아닌 값 {ur_na_count}개 발견")
                    non_numeric_ur = df[pd.to_numeric(df['UR'], errors='coerce').isna()]
                    print(f"    비숫자 UR 값들: {list(non_numeric_ur['UR'].unique())}")
                
                if z_na_count > 0:
                    print(f"    WARNING: Z 컬럼에 숫자가 아닌 값 {z_na_count}개 발견")
                    non_numeric_z = df[pd.to_numeric(df['Z'], errors='coerce').isna()]
                    print(f"    비숫자 Z 값들: {list(non_numeric_z['Z'].unique())}")
                
                # 숫자로 변환 가능한 데이터만 사용
                df['UR_numeric'] = ur_numeric
                df['Z_numeric'] = z_numeric
                
                # NaN 값 제거
                valid_df = df.dropna(subset=['UR_numeric', 'Z_numeric'])
                
                if len(valid_df) < len(df):
                    print(f"    유효한 데이터: {len(df)} → {len(valid_df)} (숫자 변환 실패 {len(df) - len(valid_df)}개 제거)")
                
                if valid_df.empty:
                    print(f"    ERROR: 유효한 숫자 데이터가 없음")
                    results['error_files'].append(f"{chr_name} (no valid numeric data)")
                    continue
                
                # 필터링 적용: 
                # 1. type이 'N'인 샘플은 무조건 유지
                # 2. type이 'N'이 아닌 샘플만 UR >= UR_min AND Z >= Z_min 조건 적용
                
                # N 타입 샘플 (무조건 유지)
                n_samples = valid_df[valid_df['type'] == 'N']
                
                # N이 아닌 타입 샘플 (임계값 필터링 적용)
                non_n_samples = valid_df[valid_df['type'] != 'N']
                filtered_non_n = non_n_samples[(non_n_samples['UR_numeric'] >= ur_min) & (non_n_samples['Z_numeric'] >= z_min)]
                
                # N 샘플 + 필터링된 non-N 샘플 결합
                filtered_df = pd.concat([n_samples, filtered_non_n], ignore_index=True)
                
                # 원본 컬럼으로 복원 (UR_numeric, Z_numeric 제거)
                filtered_df = filtered_df.drop(['UR_numeric', 'Z_numeric'], axis=1)
                
                filtered_count = len(filtered_df)
                removed_count = original_count - filtered_count
                
                print(f"    샘플 처리:")
                print(f"      N 타입 (유지): {len(n_samples)}개")
                print(f"      non-N 타입 (필터링): {len(non_n_samples)} → {len(filtered_non_n)}개")
                print(f"      전체: {original_count} → {filtered_count} ({removed_count}개 제거)")
                
                # 제거된 샘플들 정보 출력 (처음 5개만)
                if removed_count > 0:
                    # 제거된 샘플들 찾기 (N이 아닌 타입에서만)
                    removed_non_n = non_n_samples[~((non_n_samples['UR_numeric'] >= ur_min) & (non_n_samples['Z_numeric'] >= z_min))]
                    
                    print(f"    제거된 샘플 (N 타입 제외):")
                    for idx, row in removed_non_n.head(5).iterrows():
                        ur_val = row.get('UR', 'N/A')
                        z_val = row.get('Z', 'N/A')
                        sample_name = row.get('sample', 'Unknown')
                        sample_type = row.get('type', 'Unknown')
                        print(f"      {sample_name} ({sample_type}) - UR:{ur_val}, Z:{z_val}")
                    if len(removed_non_n) > 5:
                        print(f"      ... 외 {len(removed_non_n) - 5}개")
                    
                    print(f"    참고: N 타입 샘플 {len(n_samples)}개는 임계값과 관계없이 모두 유지됨")
                
                # 필터링된 데이터 저장 (원본 덮어쓰기)
                filtered_df.to_csv(chr_file, sep='\t', index=False)
                print(f"    파일 저장 완료: {chr_file}")
                
                # 결과 업데이트
                results['processed_files'] += 1
                results['total_samples_before'] += original_count
                results['total_samples_after'] += filtered_count
                results['removed_samples'] += removed_count
                
            except Exception as e:
                print(f"  {chr_name}: 처리 중 상세 오류")
                print(f"    오류 유형: {type(e).__name__}")
                print(f"    오류 메시지: {str(e)}")
                
                # 파일 존재 및 읽기 가능 여부 확인
                try:
                    print(f"    파일 존재: {os.path.exists(chr_file)}")
                    print(f"    파일 크기: {os.path.getsize(chr_file) if os.path.exists(chr_file) else 'N/A'} bytes")
                    
                    if os.path.exists(chr_file):
                        # 파일 첫 몇 줄 읽어보기
                        with open(chr_file, 'r') as f:
                            first_lines = [f.readline().strip() for _ in range(3)]
                        print(f"    파일 첫 3줄:")
                        for i, line in enumerate(first_lines, 1):
                            print(f"      {i}: {line}")
                except Exception as file_error:
                    print(f"    파일 정보 확인 실패: {file_error}")
                
                # 상세 에러 추적
                import traceback
                traceback.print_exc()
                
                results['error_files'].append(f"{chr_name} ({str(e)})")
                logger.error(f"Error processing {chr_name}: {e}")
        
        # 결과 요약
        print(f"\n=== 처리 완료 ===")
        print(f"처리된 파일: {results['processed_files']}/22")
        print(f"전체 샘플: {results['total_samples_before']} → {results['total_samples_after']}")
        print(f"제거된 샘플: {results['removed_samples']}")
        if results['error_files']:
            print(f"오류 파일: {results['error_files']}")
        
        return results
        
    except Exception as e:
        print(f"전체 처리 중 오류: {e}")
        logger.error(f"Overall processing error: {e}")
        return results


def filter_all_groups(base_dir, labcode):
    """
    orig, fetus, mom 모든 그룹에 대해 필터링 수행
    
    Parameters:
    - base_dir: DATA_DIR/refs/<labcode>/EZD 경로
    - labcode: 실험실 코드
    """
    
    groups = ['orig', 'fetus', 'mom']
    total_results = {}
    
    print(f"\n========== 모든 그룹 필터링 시작 ==========")
    print(f"Base 디렉토리: {base_dir}")
    print(f"Lab 코드: {labcode}")
    
    for group in groups:
        print(f"\n{'='*50}")
        print(f"그룹: {group.upper()}")
        print(f"{'='*50}")
        
        # 디렉토리 경로 설정
        chr_dir = os.path.join(base_dir, group)
        threshold_file = os.path.join(base_dir, group, f'{group}_thresholds_new.tsv')
        
        # 대안 경로도 확인
        if not os.path.exists(threshold_file):
            threshold_file = os.path.join(base_dir, f'{group}_thresholds_new.tsv')
        
        # 필터링 실행
        results = filter_chr_samples_by_threshold(chr_dir, threshold_file)
        total_results[group] = results
    
    # 전체 결과 요약
    print(f"\n{'='*60}")
    print(f"전체 처리 결과 요약")
    print(f"{'='*60}")
    
    total_processed = 0
    total_before = 0
    total_after = 0
    total_removed = 0
    
    for group, results in total_results.items():
        print(f"\n{group.upper()}:")
        print(f"  처리된 파일: {results['processed_files']}")
        print(f"  샘플 변화: {results['total_samples_before']} → {results['total_samples_after']}")
        print(f"  제거된 샘플: {results['removed_samples']}")
        
        total_processed += results['processed_files']
        total_before += results['total_samples_before']
        total_after += results['total_samples_after']
        total_removed += results['removed_samples']
    
    print(f"\n전체 합계:")
    print(f"  처리된 파일: {total_processed}")
    print(f"  전체 샘플: {total_before} → {total_after}")
    print(f"  제거된 샘플: {total_removed}")
    print(f"  제거 비율: {(total_removed/total_before*100):.1f}%" if total_before > 0 else "  제거 비율: 0%")
    
    return total_results


def backup_chr_files(chr_dir):
    """
    필터링 전 chr 파일들 백업
    
    Parameters:
    - chr_dir: chr1.txt~chr22.txt 파일들이 있는 디렉토리
    """
    
    print(f"=== 파일 백업 중 ===")
    print(f"디렉토리: {chr_dir}")
    
    backup_count = 0
    for i in range(1, 23):
        chr_file = os.path.join(chr_dir, f'chr{i}.txt')
        if os.path.exists(chr_file):
            backup_file = f'{chr_file}.backup'
            
            if not os.path.exists(backup_file):
                try:
                    with open(chr_file, 'r') as original:
                        with open(backup_file, 'w') as backup:
                            backup.write(original.read())
                    print(f"  백업: chr{i}.txt")
                    backup_count += 1
                except Exception as e:
                    print(f"  백업 실패 chr{i}.txt: {e}")
            else:
                print(f"  이미 존재: chr{i}.txt.backup")
    
    print(f"백업 완료: {backup_count}개 파일")


def main():
    """
    메인 함수 - 명령행 인자 처리
    """
    import sys
    import argparse
    
    parser = argparse.ArgumentParser(description='염색체 파일 임계값 기반 필터링')
    parser.add_argument('--base-dir', required=True, help='Base directory (DATA_DIR/refs/<labcode>/EZD)')
    parser.add_argument('--labcode', required=True, help='Lab code')
    parser.add_argument('--group', help='특정 그룹만 처리 (orig, fetus, mom)')
    parser.add_argument('--backup', action='store_true', help='처리 전 백업 수행')
    parser.add_argument('--dry-run', action='store_true', help='실제 변경 없이 시뮬레이션만')
    
    args = parser.parse_args()
    
    if args.dry_run:
        print("DRY RUN 모드: 실제 파일은 변경되지 않습니다.")
        return
    
    base_dir = args.base_dir
    labcode = args.labcode
    
    if args.backup:
        # 모든 그룹 백업
        groups = [args.group] if args.group else ['orig', 'fetus', 'mom']
        for group in groups:
            chr_dir = os.path.join(base_dir, group)
            if os.path.exists(chr_dir):
                backup_chr_files(chr_dir)
    
    if args.group:
        # 특정 그룹만 처리
        chr_dir = os.path.join(base_dir, args.group)
        threshold_file = os.path.join(base_dir, args.group, f'{args.group}_thresholds_new.tsv')
        
        if not os.path.exists(threshold_file):
            threshold_file = os.path.join(base_dir, f'{args.group}_thresholds_new.tsv')
        
        filter_chr_samples_by_threshold(chr_dir, threshold_file)
    else:
        # 모든 그룹 처리
        filter_all_groups(base_dir, labcode)


if __name__ == "__main__":
    # 사용 예시
    if len(sys.argv) == 1:
        print("염색체 파일 임계값 기반 필터링 스크립트")
        print("\n사용법:")
        print("  python script.py --base-dir /data/refs/LAB001/EZD --labcode LAB001")
        print("  python script.py --base-dir /data/refs/LAB001/EZD --labcode LAB001 --group orig")
        print("  python script.py --base-dir /data/refs/LAB001/EZD --labcode LAB001 --backup")
        print("  python script.py --base-dir /data/refs/LAB001/EZD --labcode LAB001 --dry-run")
        print("\n옵션:")
        print("  --base-dir PATH    # Base directory (DATA_DIR/refs/<labcode>/EZD)")
        print("  --labcode CODE     # Lab code")
        print("  --group GROUP      # 특정 그룹만 처리 (orig, fetus, mom)")
        print("  --backup           # 처리 전 백업 수행")
        print("  --dry-run          # 실제 변경 없이 시뮬레이션")
    else:
        main()
