#!/usr/bin/env python3
"""
염색체 파일 일괄 처리 스크립트
- chr1.txt ~ chr22.txt 파일들을 처리
- type 컬럼이 'N'이 아닌 경우 'P'로 변경
- sample 컬럼을 P1, P2, P3... 또는 N1, N2, N3... 형태로 변경
"""

import pandas as pd
import os
import glob

def process_chromosome_file(file_path, enable_filtering=False):
    """
    개별 염색체 파일을 처리하는 함수
    
    Parameters:
    - file_path: 처리할 파일 경로
    - enable_filtering: 필터링 기능 활성화 여부
    """
    
    print(f"Processing {file_path}...")
    
    try:
        # 파일 읽기
        df = pd.read_csv(file_path, sep='\t')
        
        # 컬럼명 확인
        if not all(col in df.columns for col in ['sample', 'type', 'UR', 'Z']):
            print(f"Warning: {file_path} does not have expected columns")
            return False
        
        original_count = len(df)
        
        # 필터링 적용 (옵션)
        if enable_filtering:
            print(f"  - Applying filters...")
            
            # 1. sample이 'ON1'로 시작하는 것 제거
            before_sample_filter = len(df)
            df = df[~df['sample'].str.startswith('ON1', na=False)]
            after_sample_filter = len(df)
            print(f"    * Sample filter (remove ON1*): {before_sample_filter} -> {after_sample_filter} rows")
            
            # 2. sample이 'ED'로 시작하는 것 제거
            before_ed_filter = len(df)
            df = df[~df['sample'].str.startswith('ED', na=False)]
            after_ed_filter = len(df)
            print(f"    * Sample filter (remove ED*): {before_ed_filter} -> {after_ed_filter} rows")
            
            if len(df) == 0:
                print(f"    * No rows remaining after sample filters, skipping file")
                return False
            
            # 4. type이 '_5p'로 끝나는 것과 'TP'만 유지 (N, _5p, TP 타입만 유지)
            before_5p_filter = len(df)
            df = df[(df['type'] == 'N') | (df['type'].str.endswith('_5p', na=False)) | (df['type'] == 'TP')]
            after_5p_filter = len(df)
            print(f"    * Keep only N, *5p, and TP filter: {before_5p_filter} -> {after_5p_filter} rows")

            # 3. type에 'FP'가 포함된 것 제거
            #before_fp_filter = len(df)
            #df = df[~df['type'].str.contains('FP', na=False)]
            #after_fp_filter = len(df)
            #print(f"    * FP filter: {before_fp_filter} -> {after_fp_filter} rows")
            
            # 4. type이 '_3p'로 끝나는 것 제거
            #before_3p_filter = len(df)
            #df = df[~df['type'].str.endswith('_3p', na=False)]
            #after_3p_filter = len(df)
            #print(f"    * _3p filter: {before_3p_filter} -> {after_3p_filter} rows")
            
            # 4. type이 '_3p'로 끝나는 것 제거
            #before_15p_filter = len(df)
            #df = df[~df['type'].str.endswith('_15p', na=False)]
            #after_15p_filter = len(df)
            #print(f"    * _15p filter: {before_15p_filter} -> {after_15p_filter} rows")

            print(f"  - Total filtering: {original_count} -> {len(df)} rows")
        
        if len(df) == 0:
            print(f"  - No data remaining after filtering, skipping file")
            return False
        
        # P 그룹과 N 그룹 분리
        p_mask = df['type'] != 'N'
        n_mask = df['type'] == 'N'
        
        # P 그룹 처리
        p_count = 0
        for idx in df[p_mask].index:
            p_count += 1
            df.loc[idx, 'sample'] = f'P{p_count}'
            df.loc[idx, 'type'] = 'P'
        
        # N 그룹 처리
        n_count = 0
        for idx in df[n_mask].index:
            n_count += 1
            df.loc[idx, 'sample'] = f'N{n_count}'
            # type은 이미 'N'이므로 그대로 유지
        
        # 파일 저장 (원본 덮어쓰기)
        df.to_csv(file_path, sep='\t', index=False)
        
        print(f"  - Processed {p_count} P samples and {n_count} N samples")
        if enable_filtering:
            print(f"  - Final data: {len(df)} rows (filtered from {original_count})")
        
        return True
        
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return False


def process_all_chromosome_files(directory=None, enable_filtering=False):
    """
    지정된 디렉토리의 모든 염색체 파일을 처리하는 메인 함수
    
    Parameters:
    - directory: 처리할 파일들이 있는 디렉토리 경로 (None이면 현재 디렉토리)
    - enable_filtering: 필터링 기능 활성화 여부
    """
    
    if directory is None:
        directory = os.getcwd()
    
    # 디렉토리 존재 확인
    if not os.path.exists(directory):
        print(f"Error: Directory '{directory}' does not exist!")
        return False
    
    if not os.path.isdir(directory):
        print(f"Error: '{directory}' is not a directory!")
        return False
    
    print(f"=== 염색체 파일 일괄 처리 시작 ===")
    print(f"처리 디렉토리: {os.path.abspath(directory)}")
    if enable_filtering:
        print("필터링 모드: 활성화")
        print("  - sample이 'ON1'로 시작하는 것 제거")
        print("  - sample이 'ED'로 시작하는 것 제거")
        print("  - type에 'FP'가 포함된 것 제거")
        print("  - type이 '_3p'로 끝나는 것 제거")
    else:
        print("필터링 모드: 비활성화")
    print()
    
    # chr1.txt ~ chr22.txt 파일 목록 생성
    chr_files = []
    for i in range(1, 23):
        file_name = f'chr{i}.txt'
        file_path = os.path.join(directory, file_name)
        if os.path.exists(file_path):
            chr_files.append(file_path)
        else:
            print(f"Warning: {file_name} not found in {directory}")
    
    if not chr_files:
        print(f"No chromosome files found in directory: {directory}")
        return False
    
    print(f"Found {len(chr_files)} chromosome files:")
    for file_path in chr_files:
        print(f"  - {os.path.basename(file_path)}")
    print()
    
    # 각 파일 처리
    success_count = 0
    for file_path in chr_files:
        if process_chromosome_file(file_path, enable_filtering):
            success_count += 1
        print()  # 빈 줄 추가
    
    print("=== 처리 완료 ===")
    print(f"성공: {success_count}/{len(chr_files)} 파일")
    
    if success_count == len(chr_files):
        print("모든 파일이 성공적으로 처리되었습니다!")
        return True
    else:
        print(f"{len(chr_files) - success_count}개 파일에서 오류가 발생했습니다.")
        return False


def preview_changes(file_path, num_rows=5, enable_filtering=False):
    """
    파일 변경사항을 미리보기하는 함수
    
    Parameters:
    - file_path: 미리볼 파일 경로
    - num_rows: 표시할 행 수
    - enable_filtering: 필터링 기능 활성화 여부
    """
    
    try:
        # 원본 파일 읽기
        df_original = pd.read_csv(file_path, sep='\t')
        
        print(f"\n=== {os.path.basename(file_path)} 원본 (상위 {num_rows}행) ===")
        print(df_original.head(num_rows).to_string(index=False))
        print(f"원본 총 행 수: {len(df_original)}")
        
        # 변경 로직 적용 (미리보기용, 파일은 저장하지 않음)
        df_preview = df_original.copy()
        
        # 필터링 적용 (옵션)
        if enable_filtering:
            print(f"\n=== 필터링 적용 ===")
            
            # 1. sample이 'ON1'로 시작하는 것 제거
            before_count = len(df_preview)
            df_preview = df_preview[~df_preview['sample'].str.startswith('ON1', na=False)]
            print(f"Sample 필터 (remove ON1*): {before_count} -> {len(df_preview)} 행")
            
            # 2. sample이 'ED'로 시작하는 것 제거
            before_count = len(df_preview)
            df_preview = df_preview[~df_preview['sample'].str.startswith('ED', na=False)]
            print(f"Sample 필터 (remove ED*): {before_count} -> {len(df_preview)} 행")
            
            # 3. type에 'FP'가 포함된 것 제거
            before_count = len(df_preview)
            df_preview = df_preview[~df_preview['type'].str.contains('FP', na=False)]
            print(f"FP 필터: {before_count} -> {len(df_preview)} 행")
            
            # 4. type이 '_3p'로 끝나는 것 제거
            before_count = len(df_preview)
            df_preview = df_preview[~df_preview['type'].str.endswith('_3p', na=False)]
            print(f"_3p 필터: {before_count} -> {len(df_preview)} 행")
        
        if len(df_preview) == 0:
            print("필터링 후 데이터가 없습니다.")
            return
        
        # P 그룹과 N 그룹 분리
        p_mask = df_preview['type'] != 'N'
        n_mask = df_preview['type'] == 'N'
        
        # P 그룹 처리
        p_count = 0
        for idx in df_preview[p_mask].index:
            p_count += 1
            df_preview.loc[idx, 'sample'] = f'P{p_count}'
            df_preview.loc[idx, 'type'] = 'P'
        
        # N 그룹 처리
        n_count = 0
        for idx in df_preview[n_mask].index:
            n_count += 1
            df_preview.loc[idx, 'sample'] = f'N{n_count}'
        
        print(f"\n=== {os.path.basename(file_path)} 변경 후 (상위 {num_rows}행) ===")
        print(df_preview.head(num_rows).to_string(index=False))
        
        print(f"\n변경 요약:")
        if enable_filtering:
            print(f"  - 필터링: {len(df_original)} -> {len(df_preview)} 행")
        print(f"  - P 샘플: {p_count}개")
        print(f"  - N 샘플: {n_count}개")
        print(f"  - 최종 데이터: {len(df_preview)} 행")
        
    except Exception as e:
        print(f"미리보기 오류: {e}")


def backup_files(directory=None):
    """
    처리 전 원본 파일들을 백업하는 함수
    
    Parameters:
    - directory: 백업할 파일들이 있는 디렉토리 경로 (None이면 현재 디렉토리)
    """
    
    if directory is None:
        directory = os.getcwd()
    
    print(f"=== 원본 파일 백업 중 ===")
    print(f"백업 디렉토리: {os.path.abspath(directory)}")
    
    backup_count = 0
    for i in range(1, 23):
        file_name = f'chr{i}.txt'
        file_path = os.path.join(directory, file_name)
        if os.path.exists(file_path):
            backup_path = f'{file_path}.backup'
            try:
                # 이미 백업 파일이 있으면 건너뛰기
                if not os.path.exists(backup_path):
                    with open(file_path, 'r') as original:
                        with open(backup_path, 'w') as backup:
                            backup.write(original.read())
                    print(f"  - {file_name} -> {file_name}.backup")
                    backup_count += 1
                else:
                    print(f"  - {file_name}.backup already exists, skipping")
            except Exception as e:
                print(f"  - Error backing up {file_name}: {e}")
    
    print(f"백업 완료: {backup_count}개 파일\n")


def get_directory_from_args():
    """
    명령행 인자에서 디렉토리 경로를 추출하는 함수
    
    Returns:
    - str: 디렉토리 경로 또는 None
    """
    import sys
    
    # 명령행 인자에서 디렉토리 찾기
    for i, arg in enumerate(sys.argv):
        if arg in ['--dir', '-d', '--directory']:
            if i + 1 < len(sys.argv):
                return sys.argv[i + 1]
            else:
                print("Error: --dir option requires a directory path")
                return None
        elif arg.startswith('--dir='):
            return arg.split('=', 1)[1]
    
    return None


def has_filter_flag():
    """
    명령행 인자에서 필터링 플래그를 확인하는 함수
    
    Returns:
    - bool: 필터링 활성화 여부
    """
    import sys
    return '--filter' in sys.argv or '-f' in sys.argv


if __name__ == "__main__":
    import sys
    
    # 디렉토리 경로와 필터링 옵션 가져오기
    target_directory = get_directory_from_args()
    enable_filtering = has_filter_flag()
    
    # 명령행 인자 확인
    if len(sys.argv) > 1:
        if '--preview' in sys.argv:
            # 미리보기 모드
            print("=== 미리보기 모드 ===")
            if target_directory:
                print(f"대상 디렉토리: {os.path.abspath(target_directory)}")
            else:
                target_directory = os.getcwd()
                print(f"현재 디렉토리: {os.path.abspath(target_directory)}")
            
            if enable_filtering:
                print("필터링: 활성화")
            else:
                print("필터링: 비활성화")
            
            # 첫 번째 파일만 미리보기
            for i in range(1, 23):
                file_name = f'chr{i}.txt'
                file_path = os.path.join(target_directory, file_name)
                if os.path.exists(file_path):
                    preview_changes(file_path, enable_filtering=enable_filtering)
                    break
        
        elif '--backup' in sys.argv:
            # 백업 후 처리
            backup_files(target_directory)
            process_all_chromosome_files(target_directory, enable_filtering)
        
        elif '--help' in sys.argv or '-h' in sys.argv:
            print("염색체 파일 일괄 처리 스크립트")
            print("\n사용법:")
            print("  python script.py [옵션] [--dir 디렉토리경로] [--filter]")
            print("\n기본 옵션:")
            print("  (없음)           # 바로 처리 (확인 후)")
            print("  --preview        # 미리보기 (파일 변경 없음)")
            print("  --backup         # 백업 후 처리")
            print("  --help, -h       # 도움말")
            print("\n디렉토리 지정:")
            print("  --dir PATH       # 처리할 디렉토리 지정")
            print("  --directory PATH # 처리할 디렉토리 지정")
            print("  -d PATH          # 처리할 디렉토리 지정")
            print("  --dir=PATH       # 처리할 디렉토리 지정")
            print("\n필터링 옵션:")
            print("  --filter, -f     # 필터링 활성화")
            print("    * sample이 'ON1'로 시작하는 것 제거")
            print("    * sample이 'ED'로 시작하는 것 제거")
            print("    * type에 'FP'가 포함된 것 제거")
            print("    * type이 '_3p'로 끝나는 것 제거")
            print("\n예시:")
            print("  # 기본 처리")
            print("  python script.py --dir /path/to/chr/files")
            print("  ")
            print("  # 필터링과 함께 처리")
            print("  python script.py --dir /data/chr --filter")
            print("  ")
            print("  # 필터링 미리보기")
            print("  python script.py --preview --dir /data/chr --filter")
            print("  ")
            print("  # 백업 후 필터링 처리")
            print("  python script.py --backup --dir /data/chr --filter")
        
        else:
            # 기본 모드: 바로 처리
            if target_directory:
                print(f"대상 디렉토리: {os.path.abspath(target_directory)}")
            else:
                target_directory = os.getcwd()
                print(f"현재 디렉토리: {os.path.abspath(target_directory)}")
            
            if enable_filtering:
                print("필터링 모드: 활성화")
                print("  - sample이 'ON1'로 시작하는 것 제거")
                print("  - sample이 'ED'로 시작하는 것 제거")
                print("  - type에 'FP'가 포함된 것 제거")
                print("  - type이 '_3p'로 끝나는 것 제거")
            
            response = input("원본 파일을 덮어쓰게 됩니다. 계속하시겠습니까? (y/N): ")
            if response.lower() in ['y', 'yes']:
                process_all_chromosome_files(target_directory, enable_filtering)
            else:
                print("처리가 취소되었습니다.")
                print("미리보기를 원하면: python script.py --preview --dir <directory> [--filter]")
                print("백업 후 처리를 원하면: python script.py --backup --dir <directory> [--filter]")
    
    else:
        # 인자가 없는 경우 - 현재 디렉토리에서 처리
        print("현재 디렉토리에서 처리합니다.")
        print(f"디렉토리: {os.path.abspath(os.getcwd())}")
        response = input("원본 파일을 덮어쓰게 됩니다. 계속하시겠습니까? (y/N): ")
        if response.lower() in ['y', 'yes']:
            process_all_chromosome_files()
        else:
            print("처리가 취소되었습니다.")
            print("도움말을 보려면: python script.py --help")
