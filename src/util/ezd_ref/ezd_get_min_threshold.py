#!/usr/bin/env python3
"""
임계값 자동 계산 및 업데이트 스크립트
- chr1.txt~chr22.txt 파일에서 N 타입 샘플의 최대값 기준으로 임계값 계산
- UR_min = max(N samples UR) + 0.01
- Z_min = max(N samples Z) + 0.01
- 기존 threshold.tsv 파일을 직접 수정
"""

import pandas as pd
import os
import logging
from typing import Dict, Tuple, Optional

# 로거 설정
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def calculate_threshold_from_chr_files(chr_dir: str) -> Dict[str, Dict[str, float]]:
    """
    염색체 파일들에서 N 타입 샘플의 최대값을 기준으로 임계값 계산
    
    Parameters:
    - chr_dir: chr1.txt~chr22.txt 파일들이 있는 디렉토리
    
    Returns:
    - dict: {chr_name: {'UR_min': value, 'Z_min': value}}
    """
    
    thresholds = {}
    
    print(f"\n=== 임계값 계산 시작 ===")
    print(f"디렉토리: {chr_dir}")
    
    if not os.path.exists(chr_dir):
        print(f"ERROR: 디렉토리가 존재하지 않습니다: {chr_dir}")
        return thresholds
    
    # 각 염색체 파일 처리
    for i in range(1, 23):
        chr_name = f'chr{i}'
        chr_file = os.path.join(chr_dir, f'{chr_name}.txt')
        
        if not os.path.exists(chr_file):
            print(f"  {chr_name}: 파일 없음")
            continue
        
        try:
            # 염색체 파일 로드
            df = pd.read_csv(chr_file, sep='\t')
            
            if df.empty:
                print(f"  {chr_name}: 빈 파일")
                continue
            
            # 필수 컬럼 확인
            if not all(col in df.columns for col in ['type', 'UR', 'Z']):
                print(f"  {chr_name}: 필수 컬럼 없음 (type, UR, Z)")
                continue
            
            # N 타입 샘플만 필터링
            n_samples = df[df['type'] == 'N']
            
            if n_samples.empty:
                print(f"  {chr_name}: N 타입 샘플 없음")
                continue
            
            # UR, Z 값을 숫자로 변환 (에러 처리)
            ur_values = pd.to_numeric(n_samples['UR'], errors='coerce').dropna()
            z_values = pd.to_numeric(n_samples['Z'], errors='coerce').dropna()
            
            if ur_values.empty or z_values.empty:
                print(f"  {chr_name}: 유효한 숫자 데이터 없음")
                continue
            
            # 최대값 계산
            max_ur = ur_values.max()
            max_z = z_values.max()
            
            # 임계값 계산 (최대값 + 0.01)
            ur_min = max_ur + 0.01
            z_min = max_z + 0.01
            
            thresholds[chr_name] = {
                'UR_min': ur_min,
                'Z_min': z_min
            }
            
            print(f"  {chr_name}: N샘플 {len(n_samples)}개")
            print(f"    UR: max={max_ur:.4f} → UR_min={ur_min:.4f}")
            print(f"    Z:  max={max_z:.4f} → Z_min={z_min:.4f}")
            
        except Exception as e:
            print(f"  {chr_name}: 처리 중 오류 - {e}")
            logger.error(f"Error processing {chr_name}: {e}")
    
    print(f"\n임계값 계산 완료: {len(thresholds)}/22개 염색체")
    return thresholds


def update_threshold_file(threshold_file: str, calculated_thresholds: Dict[str, Dict[str, float]], 
                         backup: bool = True) -> bool:
    """
    기존 threshold 파일의 UR_min, Z_min 값을 업데이트
    
    Parameters:
    - threshold_file: 기존 threshold.tsv 파일 경로
    - calculated_thresholds: 계산된 임계값들
    - backup: 백업 파일 생성 여부
    
    Returns:
    - bool: 성공 여부
    """
    
    try:
        print(f"\n=== 임계값 파일 업데이트 ===")
        print(f"파일: {threshold_file}")
        
        if not os.path.exists(threshold_file):
            print(f"ERROR: 임계값 파일이 존재하지 않습니다: {threshold_file}")
            return False
        
        # 백업 생성
        if backup:
            backup_file = f"{threshold_file}.backup"
            if not os.path.exists(backup_file):
                with open(threshold_file, 'r') as original:
                    with open(backup_file, 'w') as backup_f:
                        backup_f.write(original.read())
                print(f"백업 생성: {backup_file}")
            else:
                print(f"백업 파일 이미 존재: {backup_file}")
        
        # 기존 임계값 파일 로드
        df = pd.read_csv(threshold_file, sep='\t')
        print(f"기존 임계값 로드: {len(df)}개 염색체")
        print(f"기존 컬럼: {list(df.columns)}")
        
        updated_count = 0
        
        # 각 염색체별로 임계값 업데이트
        for chr_name, thresholds in calculated_thresholds.items():
            # 해당 염색체 행 찾기
            mask = df['chr'] == chr_name
            
            if not mask.any():
                print(f"  {chr_name}: 기존 파일에 없음, 건너뜀")
                continue
            
            # 기존 값 가져오기
            old_ur_min = df.loc[mask, 'UR_min'].iloc[0] if 'UR_min' in df.columns else None
            old_z_min = df.loc[mask, 'Z_min'].iloc[0] if 'Z_min' in df.columns else None
            
            # 새 값으로 업데이트
            df.loc[mask, 'UR_min'] = round(thresholds['UR_min'], 2)
            df.loc[mask, 'Z_min'] = round(thresholds['Z_min'], 2)
            
            print(f"  {chr_name}:")
            print(f"    UR_min: {old_ur_min} → {thresholds['UR_min']:.4f}")
            print(f"    Z_min:  {old_z_min} → {thresholds['Z_min']:.4f}")
            
            updated_count += 1
        
        # 업데이트된 파일 저장
        df.to_csv(threshold_file, sep='\t', index=False)
        
        print(f"\n업데이트 완료: {updated_count}개 염색체")
        print(f"파일 저장: {threshold_file}")
        
        return True
        
    except Exception as e:
        print(f"임계값 파일 업데이트 중 오류: {e}")
        logger.error(f"Error updating threshold file: {e}")
        return False


def calculate_and_update_thresholds(chr_dir: str, threshold_file: str, backup: bool = True) -> bool:
    """
    염색체 파일들에서 임계값을 계산하고 threshold 파일을 업데이트하는 메인 함수
    
    Parameters:
    - chr_dir: chr1.txt~chr22.txt 파일들이 있는 디렉토리
    - threshold_file: 업데이트할 threshold.tsv 파일 경로
    - backup: 백업 파일 생성 여부
    
    Returns:
    - bool: 성공 여부
    """
    
    print(f"========== 임계값 자동 계산 및 업데이트 ==========")
    print(f"Chr 디렉토리: {chr_dir}")
    print(f"Threshold 파일: {threshold_file}")
    print(f"백업 생성: {'예' if backup else '아니오'}")
    
    # 1. 임계값 계산
    calculated_thresholds = calculate_threshold_from_chr_files(chr_dir)
    
    if not calculated_thresholds:
        print("ERROR: 계산된 임계값이 없습니다.")
        return False
    
    # 2. 임계값 파일 업데이트
    success = update_threshold_file(threshold_file, calculated_thresholds, backup)
    
    if success:
        print(f"\n{'='*50}")
        print(f"임계값 업데이트 완료!")
        print(f"업데이트된 염색체: {len(calculated_thresholds)}개")
        print(f"{'='*50}")
    else:
        print(f"\n{'='*50}")
        print(f"임계값 업데이트 실패!")
        print(f"{'='*50}")
    
    return success


def preview_threshold_changes(chr_dir: str, threshold_file: str) -> None:
    """
    임계값 변경사항을 미리보기하는 함수
    
    Parameters:
    - chr_dir: chr1.txt~chr22.txt 파일들이 있는 디렉토리
    - threshold_file: 기존 threshold.tsv 파일 경로
    """
    
    print(f"\n=== 임계값 변경사항 미리보기 ===")
    
    # 새로운 임계값 계산
    calculated_thresholds = calculate_threshold_from_chr_files(chr_dir)
    
    if not calculated_thresholds:
        print("계산된 임계값이 없습니다.")
        return
    
    if not os.path.exists(threshold_file):
        print(f"기존 임계값 파일이 없습니다: {threshold_file}")
        return
    
    # 기존 임계값 로드
    existing_df = pd.read_csv(threshold_file, sep='\t')
    
    print(f"\n변경사항 요약:")
    print(f"{'Chromosome':<10} {'기존 UR_min':<12} {'새 UR_min':<12} {'기존 Z_min':<12} {'새 Z_min':<12}")
    print(f"{'-'*60}")
    
    for chr_name in sorted(calculated_thresholds.keys()):
        mask = existing_df['chr'] == chr_name
        
        if mask.any():
            old_ur = existing_df.loc[mask, 'UR_min'].iloc[0] if 'UR_min' in existing_df.columns else 'N/A'
            old_z = existing_df.loc[mask, 'Z_min'].iloc[0] if 'Z_min' in existing_df.columns else 'N/A'
        else:
            old_ur, old_z = 'N/A', 'N/A'
        
        new_ur = calculated_thresholds[chr_name]['UR_min']
        new_z = calculated_thresholds[chr_name]['Z_min']
        
        print(f"{chr_name:<10} {old_ur:<12} {new_ur:<12.4f} {old_z:<12} {new_z:<12.4f}")


def main():
    """
    메인 함수 - 명령행 인자 처리
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='염색체 파일 기반 임계값 자동 계산 및 업데이트')
    parser.add_argument('--chr-dir', required=True, help='chr1.txt~chr22.txt 파일들이 있는 디렉토리')
    parser.add_argument('--threshold-file', required=True, help='업데이트할 threshold.tsv 파일 경로')
    parser.add_argument('--no-backup', action='store_true', help='백업 파일 생성 안함')
    parser.add_argument('--preview', action='store_true', help='변경사항 미리보기만 (실제 변경 안함)')
    
    args = parser.parse_args()
    
    chr_dir = args.chr_dir
    threshold_file = args.threshold_file
    backup = not args.no_backup
    preview = args.preview
    
    if preview:
        # 미리보기만
        preview_threshold_changes(chr_dir, threshold_file)
    else:
        # 실제 업데이트
        calculate_and_update_thresholds(chr_dir, threshold_file, backup)


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) == 1:
        print("임계값 자동 계산 및 업데이트 스크립트")
        print("\n사용법:")
        print("  python script.py --chr-dir /path/to/chr/files --threshold-file /path/to/threshold.tsv")
        print("  python script.py --chr-dir /data/refs/LAB001/EZD/orig --threshold-file /data/refs/LAB001/EZD/orig/orig_thresholds.tsv")
        print("\n옵션:")
        print("  --chr-dir PATH         # chr1.txt~chr22.txt 파일들이 있는 디렉토리")
        print("  --threshold-file PATH  # 업데이트할 threshold.tsv 파일")
        print("  --no-backup           # 백업 파일 생성 안함")
        print("  --preview             # 변경사항 미리보기만 (실제 변경 안함)")
        print("\n예시:")
        print("  # 미리보기")
        print("  python script.py --chr-dir /data/chr --threshold-file /data/threshold.tsv --preview")
        print("  ")
        print("  # 실제 업데이트")
        print("  python script.py --chr-dir /data/chr --threshold-file /data/threshold.tsv")
        print("  ")
        print("  # 백업 없이 업데이트")
        print("  python script.py --chr-dir /data/chr --threshold-file /data/threshold.tsv --no-backup")
    else:
        main()
