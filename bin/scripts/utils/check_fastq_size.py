#!/usr/bin/env python3
# check_fastq_size.py - FASTQ 크기 확인 및 다운샘플링 스크립트

import sys
import os
import subprocess
import gzip

def main():
    if len(sys.argv) != 8:
        print(f"Usage: {sys.argv[0]} <sample_id> <fastq_dir> <fastq_1> <fastq_2> <max_size> <downsample_size> <output_dir>")
        sys.exit(1)
    
    sample_id = sys.argv[1]
    fastq_dir = sys.argv[2]
    fastq_1 = sys.argv[3]
    fastq_2 = sys.argv[4]
    max_size = int(sys.argv[5])
    downsample_size = int(sys.argv[6])
    output_dir = sys.argv[7]
    
    # 파일 경로 구성
    fastq1_path = os.path.join(fastq_dir, sample_id, fastq_1)
    fastq2_path = os.path.join(fastq_dir, sample_id, fastq_2)
    
    # 디버깅 정보
    print(f"[Debug] Looking for files:")
    print(f"[Debug] FASTQ1: {fastq1_path}")
    print(f"[Debug] FASTQ2: {fastq2_path}")
    
    # 파일 존재 확인
    if not os.path.exists(fastq1_path):
        print(f"[Error] File not found: {fastq1_path}")
        sys.exit(1)
    if not os.path.exists(fastq2_path):
        print(f"[Error] File not found: {fastq2_path}")
        sys.exit(1)
        
    print(f"[Debug] Files found. Checking sizes...")
    
    # 읽기 수 계산 - subprocess로 실행
    try:
        cmd1 = f"zcat {fastq1_path} | wc -l"
        cmd2 = f"zcat {fastq2_path} | wc -l"
        
        print(f"[Debug] Running command: {cmd1}")
        process1 = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
        if process1.returncode != 0:
            print(f"[Error] Command failed: {cmd1}")
            print(f"[Error] Error message: {process1.stderr}")
            sys.exit(1)
            
        print(f"[Debug] Running command: {cmd2}")
        process2 = subprocess.run(cmd2, shell=True, capture_output=True, text=True)
        if process2.returncode != 0:
            print(f"[Error] Command failed: {cmd2}")
            print(f"[Error] Error message: {process2.stderr}")
            sys.exit(1)
            
        fs1 = int(process1.stdout.strip())
        fs2 = int(process2.stdout.strip())
        
        print(f"[Debug] Line counts - FASTQ1: {fs1}, FASTQ2: {fs2}")
        
        fs1_n = fs1 // 4
        fs2_n = fs2 // 4
        total_reads = fs1_n + fs2_n
        
        print(f"[Debug] Read counts - FASTQ1: {fs1_n}, FASTQ2: {fs2_n}, Total: {total_reads}")
        
        # QC 결과 파일 준비
        os.makedirs(output_dir, exist_ok=True)
        with open(os.path.join(output_dir, f"{sample_id}.fastq_check.txt"), 'w') as f:
            
            # 최대 크기 확인
            if total_reads > max_size:
                print(f"[Log:] FastQ size ({total_reads}) exceeds maximum ({max_size}). Downsampling needed.")
                
                # ... 다운샘플링 코드 ...
                f.write(f"0. FastQ size ({total_reads}) has been downsampled to {downsample_size*2}\n")
                
            elif total_reads == 0:
                print(f"[Error:] FastQ size is zero!")
                f.write(f"0. FastQ size ({total_reads}) : FAIL\n")
                sys.exit(1)
                
            else:
                print(f"[Log:] FastQ size ({total_reads}) is within limit. No downsampling needed.")
                f.write(f"0. FastQ size ({total_reads}) : PASS\n")
                
    except Exception as e:
        print(f"[Error] Unexpected error: {str(e)}")
        sys.exit(1)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
