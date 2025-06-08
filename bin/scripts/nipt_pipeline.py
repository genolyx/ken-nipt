#!/usr/bin/env python3
"""
Enhanced NIPT Pipeline with Lab-specific Support and Fetal Fraction Calculation

This pipeline provides comprehensive NIPT analysis including:
- Lab-specific configuration and resource management
- Complete BAM processing pipeline
- Y-based Fetal Fraction (YFF) calculation
- Sequence-based Fetal Fraction (seqFF) calculation
- Quality control and validation
- WiseCondor/WiseCondorX/WiseCondorFF analysis
- HMMcopy analysis

Usage:
    python3 nipt_pipeline.py --sample_name SAMPLE_001 --fastq_r1 R1.fq.gz --fastq_r2 R2.fq.gz --labcode cordlife --age 30

Author: NIPT Pipeline Development Team
Version: 2.0
"""

import os
import sys
import subprocess
import logging
import datetime
import time
import argparse
import json
import numpy as np
import pandas as pd
import shutil
import glob
import hashlib
import pysam
from pathlib import Path

sys.path.append('/Work/NIPT/bin')
try:
    from process_md_result import run_microdeletion_decision_pipeline
except ImportError as e:
    logging.warning(f"Could not import run_microdeletion_decision_pipeline module: {e}")
    process_microdeletion_result = None

try:
    from generate_json_output import build_nipt_json
except ImportError as e:
    logging.warning(f"Could not import build_nipt_json module: {e}")
    fill_json_output = None

try:
    from prizm_runner import run_multiple_prizm_analysis
except ImportError as e:
    logging.warning(f"Could not import run_multiple_prizm_analysis module: {e}")
    run_prizm_pipeline = None

try:
    from ezd_runner import run_ezd_group
except ImportError as e:
    logging.warning(f"Could not import run_ezd_group module: {e}")
    run_ezd_pipeline = None

try:
    from html_review_page import generate_nipt_html_report
except ImportError as e:
    logging.warning(f"Could not import generate_nipt_html_report module: {e}")
    run_html_pipeline = None

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [%(levelname)s] %(filename)s:%(lineno)d %(message)s',
    #format='[%(asctime)s] [%(levelname)s](%(lineno)d) %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Fixed directories for Docker environment
SRC_DIR = "/Work/NIPT/bin"
ANALYSIS_DIR = "/Work/NIPT/analysis"
LOG_DIR = "/Work/NIPT/log"
FASTQ_DIR = "/Work/NIPT/fastq"
OUTPUT_DIR = "/Work/NIPT/output"
DATA_DIR = "/Work/NIPT/data"
VERSION = "v1.0"


# -------------------------------------------
# Progress Tracker
# -------------------------------------------
class ProgressTracker:
    def __init__(self, sample_name, output_dir):
        self.sample_name = sample_name
        self.output_dir = output_dir
        self.progress_file = os.path.join(output_dir, f"{sample_name}_progress.txt")
        self.completed_file = os.path.join(output_dir, f"{sample_name}.completed")
        self.failed_file = os.path.join(output_dir, f"{sample_name}.failed")

        # 시작 시 파일 초기화
        try:
            # 출력 디렉토리가 없으면 생성
            os.makedirs(output_dir, exist_ok=True)

            # Progress 파일 생성/초기화
            with open(self.progress_file, 'w') as f:
                f.write(f"Pipeline started for {sample_name} at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        except Exception as e:
            raise Exception(f"Failed to create progress file: {e}")


    def update_step(self, step_num, step_name, status, details=""):
        """
        단계별 진행상황 업데이트
        step_num: 단계 번호 (0, 1, 2, ...)
        step_name: 단계 이름
        status: PASS 또는 FAIL
        details: 추가 정보 (파일 크기 등)
        """
        with open(self.progress_file, 'a') as f:
            if details:
                f.write(f"{step_num}. {step_name} ({details}) : {status}\n")
            else:
                f.write(f"{step_num}. {step_name} : {status}\n")

        log_and_print(f"Step {step_num} - {step_name}: {status}")

        # FAIL인 경우 즉시 실패 파일 생성
        if status == "FAIL":
            self.mark_failed(f"Failed at step {step_num}: {step_name}")

    def mark_completed(self):
        """파이프라인 성공 완료"""
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open(self.completed_file, 'w') as f:
            f.write(f"Pipeline completed successfully for {self.sample_name} at {timestamp}\n")

        with open(self.progress_file, 'a') as f:
            f.write(f"Pipeline completed successfully at {timestamp}\n")

        log_and_print(f"Pipeline completed successfully for {self.sample_name}")

    def mark_failed(self, reason="Unknown error"):
        """파이프라인 실패"""
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open(self.failed_file, 'w') as f:
            f.write(f"Pipeline failed for {self.sample_name} at {timestamp}\n")
            f.write(f"Reason: {reason}\n")

        with open(self.progress_file, 'a') as f:
            f.write(f"Pipeline failed at {timestamp} - {reason}\n")

        log_and_print(f"Pipeline failed for {self.sample_name}: {reason}", 'ERROR')

def get_lab_paths(labcode):
    """Get lab-specific file paths based on labcode"""
    base_paths = {
        "config": f"/Work/NIPT/config/{labcode}/pipeline_config.json",
        "references": f"/Work/NIPT/references/{labcode}/",
        "bed_files": f"/Work/NIPT/data/bed/{labcode}/",
        "models": f"/Work/NIPT/models/{labcode}/"
    }
    return base_paths

def get_lab_bed_paths(labcode):
    """Get lab-specific BED file paths with fallback to common"""
    lab_bed_dir = Path("/Work/NIPT/data/bed") / labcode
    common_bed_dir = Path("/Work/NIPT/data/bed/common")

    def get_bed_path(filename, lab_specific=False):
        lab_path = lab_bed_dir / filename
        common_path = common_bed_dir / filename

        if lab_specific and lab_path.exists():
            return str(lab_path)
        elif not lab_specific and common_path.exists():
            return str(common_path)
        elif lab_specific:
            return str(lab_path)  # fallback even if missing
        else:
            return str(common_path)  # fallback

    return {
        # Main filter BED files
        "of": get_bed_path("Uniform_2017_allY.bed", lab_specific=True),
        "nf08": get_bed_path("hg19_mappability_0.8_clean_all_36mer.bed", lab_specific=True),
        "nf09": get_bed_path("hg19_mappability_0.9_clean_all_36mer.bed", lab_specific=True),

        # FF calculation BED files
        "y_regions_09": get_bed_path("chrY_target_0.9.bed"),
        "y_regions_noPARs": get_bed_path("chrY_noPARs.bed"),
        "y_regions_target": get_bed_path("chrY_target.bed"),
        "autosome_control": get_bed_path("autosome_control.bed"),
    }

def get_lab_references(labcode):
    """Get lab-specific reference datasets with fallback"""
    lab_ref_dir = f"/Work/NIPT/references/{labcode}"
    common_ref_dir = "/Work/NIPT/references/common"
    
    def get_ref_path(filename):
        """Get reference file path with fallback mechanism"""
        lab_path = f"{lab_ref_dir}/{filename}"
        common_path = f"{common_ref_dir}/{filename}"
        
        if os.path.exists(lab_path):
            return lab_path
        elif os.path.exists(common_path):
            log_and_print(f"Using common reference: {common_path}")
            return common_path
        else:
            log_and_print(f"Reference file not found: {filename}", 'WARNING')
            return lab_path
    
    return {
        "wc_reference": get_ref_path("wc_reference.npz"),
        "wcx_reference": get_ref_path("wcx_reference.npz"),
        "wcff_reference": get_ref_path("wcff_reference.npz"),
        "seqff_model": get_ref_path("seqff_model.pkl"),
        "risk_model": get_ref_path("risk_calculation_model.json"),
        "prizm_mean": get_ref_path("prizm_mean.txt"),
        "prizm_sd": get_ref_path("prizm_sd.txt"),
        "prizm_mean_10mb": get_ref_path("prizm_mean_10mb.txt"),
        "prizm_sd_10mb": get_ref_path("prizm_sd_10mb.txt")
    }

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Enhanced NIPT Analysis Pipeline')
    parser.add_argument('--sample_name', required=True, help='Sample name')
    parser.add_argument('--fastq_r1', required=True, help='R1 FASTQ filename')
    parser.add_argument('--fastq_r2', required=True, help='R2 FASTQ filename')
    parser.add_argument('--labcode', required=True, help='Laboratory code (e.g., cordlife, ucl)')
    parser.add_argument('--age', required=True, help='Age of mom to calculate risk(before)')
    parser.add_argument('--config_dir', default='/Work/NIPT/config', help='Config directory path')
    parser.add_argument('--skip_ff', action='store_true', help='Skip fetal fraction calculation')
    parser.add_argument('--force', action='store_true', help='Force rerun even if completed')
    
    return parser.parse_args()

def get_downsample_size(config):
    """Get downsample parameters from config"""
    qc_config = config.get('QC', {})
    max_fq_size = qc_config.get('max_fq_size', 50000000)
    downsample_size = qc_config.get('downsample_size', 7500000)
    return max_fq_size, downsample_size

def get_qc_thresholds(config):
    """Get QC thresholds from config"""
    qc_config = config.get('QC', {})
    return {
        'orig_biqc': qc_config.get('orig_biqc', 4.0),
        'fetus_biqc': qc_config.get('fetus_biqc', 4.5),
        'mom_biqc': qc_config.get('mom_biqc', 4.5),
        'number_of_reads': qc_config.get('number_of_reads', 5000000),
        'number_of_mapped_reads': qc_config.get('number_of_mapped_reads', 5000000),
        'mapping_rate': qc_config.get('mapping_rate', 50.0),
        'duplication_rate': qc_config.get('duplication_rate', 40.0),
        'GC_content_min': qc_config.get('GC_content_min', 33.0),
        'GC_content_max': qc_config.get('GC_content_max', 55.0),
        'YFF': qc_config.get('YFF', 4.0),
        'seqFF': qc_config.get('seqFF', 4.0),
        'FFGap_1': qc_config.get('FFGap_1', 2.5),
        'FFGap_2': qc_config.get('FFGap_2', 3.0)
    }

def qc_filter(sample_name):
    """
    QC filter function that processes genome_results.txt and creates QC summary

    Args:
        analysis_dir (str): Analysis directory path
        sample_name (str): Sample name

    Returns:
        bool: True if successful, False if error
    """

    # 입력 파일과 출력 파일 경로 설정
    genome_results_file = os.path.join(ANALYSIS_DIR, sample_name, "Output_QC", "genome_results.txt")
    qc_output_file = os.path.join(ANALYSIS_DIR, sample_name, "Output_QC", f"{sample_name}.qc.txt")

    # 입력 파일 존재 확인
    if not os.path.exists(genome_results_file):
        log_and_print(f"Error: genome_results.txt not found at {genome_results_file}")
        return False

    # 출력 디렉토리 생성 (필요한 경우)
    os.makedirs(os.path.dirname(qc_output_file), exist_ok=True)

    log_and_print(f"Processing QC filter for sample: {sample_name}")
    log_and_print(f"Input file: {genome_results_file}")
    log_and_print(f"Output file: {qc_output_file}")

    try:
        # 출력 파일 초기화
        with open(qc_output_file, 'w') as f:
            pass  # 빈 파일 생성

        # 1. number of reads 처리
        cmd1 = f"grep 'number of reads' '{genome_results_file}' | sed -e 's/     //g' | sed -e 's/,//g' | sed -e 's/ /_/g' | sed -e 's/_=_/\\t/g' >> '{qc_output_file}'"
        subprocess.run(cmd1, shell=True, check=True)

        # 2. number of mapped reads 처리 (첫 번째)
        cmd2 = f"grep 'number of mapped reads' '{genome_results_file}' | sed -e 's/     //g' | sed -e 's/= /\\t/g' | cut -f 2 | sed -e 's/ (/\\t/g' | cut -f 1 | sed -e 's/^/number of mapped reads = /g' | sed -e 's/,//g' | sed -e 's/ /_/g' | sed -e 's/_=_/\\t/g' >> '{qc_output_file}'"
        subprocess.run(cmd2, shell=True, check=True)

        # 3. mapping rate 처리 (두 번째)
        cmd3 = f"grep 'number of mapped reads' '{genome_results_file}' | sed -e 's/     //g' | sed -e 's/= /\\t/g' | cut -f 2 | sed -e 's/ (/\\t/g' | cut -f 2 | sed -e 's/)//g' | sed -e 's/^/mapping rate = /g' | sed -e 's/ /_/g' | sed -e 's/_=_/\\t/g' >> '{qc_output_file}'"
        subprocess.run(cmd3, shell=True, check=True)

        # 4. number of duplicated reads 처리
        cmd4 = f"grep 'number of duplicated reads' '{genome_results_file}' | sed -e 's/     //g' | sed -e 's/,//g' | sed -e 's/ /_/g' | sed -e 's/_=_/\\t/g' >> '{qc_output_file}'"
        subprocess.run(cmd4, shell=True, check=True)

        # 5. duplication rate 처리
        cmd5 = f"grep 'duplication rate' '{genome_results_file}' | sed -e 's/     //g' | sed -e 's/ /_/g' | sed -e 's/_=_/\\t/g' >> '{qc_output_file}'"
        subprocess.run(cmd5, shell=True, check=True)

        # 6. mean mapping quality 처리
        cmd6 = f"grep 'mean mapping quality' '{genome_results_file}' | sed -e 's/     //g' | sed -e 's/ /_/g' | sed -e 's/_=_/\\t/g' >> '{qc_output_file}'"
        subprocess.run(cmd6, shell=True, check=True)

        # 7. mean coverageData 처리
        cmd7 = f"grep 'mean coverageData' '{genome_results_file}' | sed -e 's/     //g' | sed -e 's/ /_/g' | sed -e 's/_=_/\\t/g' >> '{qc_output_file}'"
        subprocess.run(cmd7, shell=True, check=True)

        # 8. GC content 처리
        cmd8 = f"grep 'GC' '{genome_results_file}' | sed -e 's/     //g' | sed -e 's/ percentage = /_content\\t/g' >> '{qc_output_file}'"
        subprocess.run(cmd8, shell=True, check=True)

        log_and_print("QC filter completed successfully!")
        log_and_print(f"Output saved to: {qc_output_file}")

        # 결과 파일 미리보기
        if os.path.exists(qc_output_file):
            log_and_print("\nQC results preview:")
            log_and_print("===================")
            with open(qc_output_file, 'r') as f:
                content = f.read().strip()
                if content:
                    log_and_print(content)
                else:
                    log_and_print("Output file is empty - no matching lines found")

        return True

    except subprocess.CalledProcessError as e:
        log_and_print(f"Error executing command: {e}")
        return False
    except Exception as e:
        log_and_print(f"Error during QC filter processing: {e}")
        return False


def qc_pass_fail_filter(qc_file, output_file, config):
    """
    QC 결과를 분석하여 PASS/FAIL을 결정하는 함수

    Args:
        qc_file (str): QC 결과 파일 경로 (*.qc.txt)
        output_file (str): 출력 파일 경로
        config (dict): 이미 로드된 설정 딕셔너리

    Returns:
        bool: True if successful, False if error
    """

    # 기존 get_qc_thresholds 함수 사용
    qc_config = get_qc_thresholds(config)

    # Check if QC file exists
    if not os.path.exists(qc_file):
        log_and_print(f"Error: QC file not found at {qc_file}")
        return False

    log_and_print(f"Processing QC PASS/FAIL filter...")
    log_and_print(f"QC file: {qc_file}")
    log_and_print(f"Output file: {output_file}")
    log_and_print(f"QC thresholds: {qc_config}")

    try:
        # 출력 디렉토리 생성 (필요한 경우)
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        # 출력 파일 초기화
        with open(output_file, 'w') as f:
            pass

        # QC 파일 읽기 및 처리
        with open(qc_file, 'r') as f:
            lines = f.readlines()

        results = []

        for line in lines:
            line_strip = line.strip()
            if not line_strip:
                continue

            line_parts = line_strip.split('\t')
            if len(line_parts) < 2:
                continue

            metric_name = line_parts[0]
            metric_value = line_parts[1]
            status = "PASS"

            try:
                if 'number_of_reads' in metric_name and 'mapped' not in metric_name and 'duplicated' not in metric_name:
                    if int(metric_value) >= qc_config["number_of_reads"]:
                        status = "PASS"
                    else:
                        status = "FAIL"

                elif 'number_of_mapped_reads' in metric_name:
                    if int(metric_value) >= qc_config["number_of_mapped_reads"]:
                        status = "PASS"
                    else:
                        status = "FAIL"

                elif 'mapping_rate' in metric_name:
                    mapping_rate = metric_value.split('%')[0].strip()
                    if float(mapping_rate) >= qc_config["mapping_rate"]:
                        status = "PASS"
                    else:
                        status = "FAIL"

                elif 'number_of_duplicated_reads' in metric_name:
                    status = "PASS"  # Always pass for duplicated reads count

                elif 'duplication_rate' in metric_name:
                    duplication_rate = metric_value.split('%')[0].strip()
                    if float(duplication_rate) < qc_config["duplication_rate"]:
                        status = "PASS"
                    else:
                        status = "FAIL"

                elif 'mean_mapping_quality' in metric_name:
                    status = "PASS"  # Always pass for mapping quality

                elif 'mean_coverageData' in metric_name:
                    status = "PASS"  # Always pass for coverage data

                elif 'GC_content' in metric_name:
                    gc_content = metric_value.split('%')[0].strip()
                    if (float(gc_content) > qc_config["GC_content_min"] and
                        float(gc_content) < qc_config["GC_content_max"]):
                        status = "PASS"
                    else:
                        status = "FAIL"

            except (ValueError, IndexError) as e:
                log_and_print(f"Warning: Could not parse metric {metric_name} with value {metric_value}: {e}")
                status = "FAIL"

            # 결과를 리스트에 추가
            result_line = f"{line_strip}\t{status}"
            results.append(result_line)

        # 결과를 파일에 저장
        with open(output_file, 'w') as f:
            for result in results:
                f.write(result + '\n')

        log_and_print(f"QC PASS/FAIL filter completed successfully!")
        log_and_print(f"Output saved to: {output_file}")

        # 결과 미리보기
        log_and_print("\nQC PASS/FAIL results:")
        log_and_print("=" * 50)
        for result in results:
            log_and_print(result)

        # 통계 출력
        pass_count = sum(1 for result in results if result.endswith('\tPASS'))
        fail_count = sum(1 for result in results if result.endswith('\tFAIL'))
        log_and_print(f"\nSummary: {pass_count} PASS, {fail_count} FAIL")

        return True

    except Exception as e:
        log_and_print(f"Error during QC PASS/FAIL processing: {e}")
        return False


def process_qc_pipeline(sample_name, config):
    """
    QC 파이프라인 전체 처리 함수 (QC filter + PASS/FAIL 판정)

    Args:
        sample_name (str): Sample name
        config (dict): 이미 로드된 설정 딕셔너리

    Returns:
        tuple: (qc_success, pass_fail_success)
    """

    # 1. QC filter 실행
    log_and_print("Running QC filter...")
    qc_success = qc_filter(sample_name)

    if not qc_success:
        log_and_print("QC filter failed!")
        return False, False

    # 2. QC PASS/FAIL 판정
    log_and_print("\nStep 2: Running QC PASS/FAIL filter...")
    qc_file = os.path.join(ANALYSIS_DIR, sample_name, "Output_QC", f"{sample_name}.qc.txt")
    output_file = os.path.join(ANALYSIS_DIR, sample_name, "Output_QC", f"{sample_name}.qc.filter.txt")

    pass_fail_success = qc_pass_fail_filter(qc_file, output_file, config)

    return qc_success, pass_fail_success


def get_ff_config(config):
    """Get fetal fraction configuration"""
    ff_config = config.get('FetalFraction', {})
    return {
        'enable_yff': ff_config.get('enable_yff', True),
        'enable_seqff': ff_config.get('enable_seqff', True),
        'yff_gender_threshold': ff_config.get('yff_gender_threshold', 0.01),
        'ff_min_threshold': ff_config.get('ff_min_threshold', 2.0),
        'ff_max_threshold': ff_config.get('ff_max_threshold', 40.0)
    }


def run_command_simple(description, command):
    """Run a command with timing and logging"""
    log_and_print(f"Running {description}")
    
    start_time = time.time()
    try:
        subprocess.run(command, shell=True, check=True)
        end_time = time.time()
        elapsed_time = end_time - start_time
        log_and_print(f"Completed {description} in {elapsed_time:.2f} seconds")
        return True
    except subprocess.CalledProcessError as e:
        log_and_print(f"Error running {description}: {e}", 'ERROR')
        return False

def run_command_error_log(description, command):
    try:
        log_and_print(f"Running: {description}")
        log_and_print(f"Command: {command}")

        result = subprocess.run(command, shell=True, capture_output=True, text=True, timeout=3600)

        if result.returncode == 0:
            log_and_print(f"Success: {description}")
            return True
        else:
            log_and_print(f"Failed: {description}", 'ERROR')
            log_and_print(f"Return code: {result.returncode}", 'ERROR')
            log_and_print(f"STDOUT: {result.stdout}", 'ERROR')
            log_and_print(f"STDERR: {result.stderr}", 'ERROR')

            # Analysis 로그 파일에도 기록
            if 'analysis_log_file' in globals():
                with open(analysis_log_file, "a") as log_file:
                    log_file.write(f"ERROR: {description}\n")
                    log_file.write(f"Command: {command}\n")
                    log_file.write(f"Return code: {result.returncode}\n")
                    log_file.write(f"STDOUT: {result.stdout}\n")
                    log_file.write(f"STDERR: {result.stderr}\n")
                    log_file.write("-" * 50 + "\n")

            return False

    except Exception as e:
        log_and_print(f"Exception in {description}: {e}", 'ERROR')
        return False

# Write a log for Success or Fail
def run_command(description, command):
    try:
        log_and_print(f"Running: {description}")
        log_and_print(f"Command: {command}")

        result = subprocess.run(command, shell=True, capture_output=True, text=True, timeout=3600)

        # 콘솔 로깅
        if result.returncode == 0:
            log_and_print(f"Success: {description}")

            if result.stdout and len(result.stdout.strip()) > 0:
                log_and_print(f"Output: {result.stdout.strip()}")

            return True
        else:
            log_and_print(f"Failed: {description}", 'ERROR')
            log_and_print(f"Return code: {result.returncode}", 'ERROR')

            if result.stderr:
                log_and_print(f"STDERR: {result.stderr}", 'ERROR')
            return False

    except subprocess.TimeoutExpired:
        log_and_print(f"Timeout in {description}", 'ERROR')
        return False

    except Exception as e:
        log_and_print(f"Exception in {description}: {e}", 'ERROR')
        return False

def run_command_realtime(description, command):
    """긴 작업을 위한 실시간 출력 함수"""
    import subprocess

    log_and_print(f"Running: {description}")
    log_and_print(f"Command: {command}")

    try:
        process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            bufsize=1
        )

        # 실시간 출력
        for line in iter(process.stdout.readline, ''):
            if line.strip():
                log_and_print(f"  {line.strip()}")

        process.wait()

        if process.returncode == 0:
            log_and_print(f"Success: {description}")
            return True
        else:
            log_and_print(f"Failed: {description} (Return code: {process.returncode})", 'ERROR')
            return False

    except Exception as e:
        log_and_print(f"Exception in {description}: {e}", 'ERROR')
        return False

def create_directories(sample_name):
    """Create required analysis and output directories"""
    analysis_subdirs = [
        "Output_QC", "Output_WC", "Output_WCX", "Output_WCFF",
        "Output_Result", "Output_hmmcopy", "Output_FF",
        "Output_EZD", "Output_PRIZM"
    ]
    output_subdirs = [
        "Output_QC", "Output_WC", "Output_WCX", "Output_WCFF",
        "Output_EZD", "Output_PRIZM"
    ]

    directories = [
        f"{ANALYSIS_DIR}/{sample_name}"
    ] + [
        f"{ANALYSIS_DIR}/{sample_name}/{subdir}" for subdir in analysis_subdirs
    ] + [
        f"{OUTPUT_DIR}/{sample_name}/{subdir}" for subdir in output_subdirs
    ] + [
        f"{LOG_DIR}/{sample_name}"
    ]

    for directory in directories:
        try:
            os.makedirs(directory, exist_ok=True)
            log_and_print(f"Created directory: {directory}")
        except Exception as e:
            log_and_print(f"Failed to create directory {directory}: {e}", 'ERROR')
            return False

    return True

def create_symbolic_links(sample_name, fastq_r1, fastq_r2):
    """Create symbolic links to FASTQ files"""
    try:
        fastq_r1_path = f"{FASTQ_DIR}/{sample_name}/{fastq_r1}"
        fastq_r2_path = f"{FASTQ_DIR}/{sample_name}/{fastq_r2}"
        
        # Check if FASTQ files exist
        if not os.path.exists(fastq_r1_path):
            log_and_print(f"FASTQ R1 file not found: {fastq_r1_path}", 'ERROR')
            return False
        if not os.path.exists(fastq_r2_path):
            log_and_print(f"FASTQ R2 file not found: {fastq_r2_path}", 'ERROR')
            return False
        
        # Define target paths for symbolic links
        r1_link_path = f"{ANALYSIS_DIR}/{fastq_r1}"
        r2_link_path = f"{ANALYSIS_DIR}/{fastq_r2}"
        
        # Remove existing symbolic links if they exist
        if os.path.exists(r1_link_path) or os.path.islink(r1_link_path):
            os.unlink(r1_link_path)
            log_and_print(f"Removed existing link: {r1_link_path}")
            
        if os.path.exists(r2_link_path) or os.path.islink(r2_link_path):
            os.unlink(r2_link_path)
            log_and_print(f"Removed existing link: {r2_link_path}")
        
        # Create new symbolic links
        os.symlink(fastq_r1_path, r1_link_path)
        log_and_print(f"Created symbolic link: {r1_link_path}")
        
        os.symlink(fastq_r2_path, r2_link_path)
        log_and_print(f"Created symbolic link: {r2_link_path}")
        
        return True

    except Exception as e:
        log_and_print(f"Error creating symbolic links: {e}", 'ERROR')
        return False

def check_file_exists_advanced(filepath, description, file_type="generic"):
    """Advanced file existence and integrity check based on file type"""
    if not os.path.exists(filepath):
        log_and_print(f"{description} not found: {filepath} - Will create")
        return False
    
    file_size = os.path.getsize(filepath)
    if file_size == 0:
        log_and_print(f"Found empty {description}: {filepath} (0 bytes) - Will recreate", 'WARNING')
        remove_file_safely(filepath)
        return False
    
    # File type specific validation
    if file_type == "bam":
        return validate_bam_file_advanced(filepath, description, file_size)
    elif file_type == "bam_index":
        return validate_bam_index(filepath, description, file_size)
    elif file_type == "npz":
        return validate_npz_file(filepath, description, file_size)
    elif file_type == "wig":
        return validate_wig_file(filepath, description, file_size)
    elif file_type == "zip":
        return validate_zip_file(filepath, description, file_size)
    else:
        # Generic validation - just check it's not empty
        log_and_print(f"Found valid {description}: {filepath} ({file_size} bytes) - Skipping step")
        return True

def validate_bam_file_advanced(filepath, description, file_size):
    """Validate BAM file with size and integrity checks"""
    # BAM files should be at least 1MB for real data
    min_bam_size = 1024 * 1024  # 1MB
    if file_size < min_bam_size:
        log_and_print(f"Found suspiciously small {description}: {filepath} ({file_size} bytes) - Will recreate", 'WARNING')
        remove_file_safely(filepath)
        return False
    
    # Check BAM integrity
    if not validate_bam_file(filepath):
        log_and_print(f"Found corrupted {description}: {filepath} - Will recreate", 'WARNING')
        remove_file_safely(filepath)
        # Also remove corresponding index
        remove_file_safely(filepath + '.bai')
        return False
    
    log_and_print(f"Found valid {description}: {filepath} ({file_size} bytes) - Skipping step")
    return True

def validate_bam_index(filepath, description, file_size):
    """Validate BAM index file"""
    # Index files should be at least 1KB
    if file_size < 1024:
        log_and_print(f"Found suspiciously small {description}: {filepath} ({file_size} bytes) - Will recreate", 'WARNING')
        remove_file_safely(filepath)
        return False
    
    log_and_print(f"Found valid {description}: {filepath} ({file_size} bytes) - Skipping step")
    return True

def validate_npz_file(filepath, description, file_size):
    """Validate NPZ file"""
    # NPZ files should be at least 10KB
    if file_size < 10240:
        log_and_print(f"Found suspiciously small {description}: {filepath} ({file_size} bytes) - Will recreate", 'WARNING')
        remove_file_safely(filepath)
        return False
    
    # Try to load the NPZ file to check integrity
    try:
        with np.load(filepath) as data:
            if len(data.files) == 0:
                log_and_print(f"Found empty NPZ file: {filepath} - Will recreate", 'WARNING')
                remove_file_safely(filepath)
                return False
    except Exception as e:
        log_and_print(f"Found corrupted NPZ file: {filepath} ({e}) - Will recreate", 'WARNING')
        remove_file_safely(filepath)
        return False
    
    log_and_print(f"Found valid {description}: {filepath} ({file_size} bytes) - Skipping step")
    return True

def validate_wig_file(filepath, description, file_size):
    """Validate WIG file"""
    # WIG files should be at least 1KB
    if file_size < 1024:
        log_and_print(f"Found suspiciously small {description}: {filepath} ({file_size} bytes) - Will recreate", 'WARNING')
        remove_file_safely(filepath)
        return False
    
    log_and_print(f"Found valid {description}: {filepath} ({file_size} bytes) - Skipping step")
    return True

def validate_zip_file(filepath, description, file_size):
    """Validate ZIP file"""
    # ZIP files should be at least 1KB
    if file_size < 1024:
        log_and_print(f"Found suspiciously small {description}: {filepath} ({file_size} bytes) - Will recreate", 'WARNING')
        remove_file_safely(filepath)
        return False
    
    # Try to open ZIP file to check integrity
    try:
        import zipfile
        with zipfile.ZipFile(filepath, 'r') as zip_file:
            if len(zip_file.namelist()) == 0:
                log_and_print(f"Found empty ZIP file: {filepath} - Will recreate", 'WARNING')
                remove_file_safely(filepath)
                return False
    except Exception as e:
        log_and_print(f"Found corrupted ZIP file: {filepath} ({e}) - Will recreate", 'WARNING')
        remove_file_safely(filepath)
        return False
    
    log_and_print(f"Found valid {description}: {filepath} ({file_size} bytes) - Skipping step")
    return True

def validate_bam_file(bam_path):
    """Validate BAM file integrity using samtools quickcheck"""
    try:
        sam_tools = os.environ.get('SAMTools', 'samtools')
        result = subprocess.run(
            [sam_tools, 'quickcheck', bam_path], 
            capture_output=True, 
            text=True, 
            timeout=30
        )
        return result.returncode == 0
    except subprocess.TimeoutExpired:
        log_and_print(f"BAM validation timeout for {bam_path}", 'WARNING')
        return False
    except Exception as e:
        log_and_print(f"Could not validate BAM file {bam_path}: {e}", 'WARNING')
        # If we can't validate, assume it's okay to avoid false positives
        return True

def remove_file_safely(filepath):
    """Safely remove a file with error handling"""
    try:
        if os.path.exists(filepath):
            os.remove(filepath)
            log_and_print(f"Removed file: {filepath}")
    except Exception as e:
        log_and_print(f"Could not remove file {filepath}: {e}", 'WARNING')

def count_fastq_reads(fastq_file):
    """Count number of reads in FASTQ file"""
    try:
        if fastq_file.endswith('.gz'):
            cmd = f"zcat {fastq_file} | wc -l"
        else:
            cmd = f"wc -l {fastq_file}"
        
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            line_count = int(result.stdout.strip().split()[0])
            read_count = line_count // 4  # FASTQ has 4 lines per read
            return read_count
        else:
            log_and_print(f"Error counting reads in {fastq_file}", 'ERROR')
            return 0
    except Exception as e:
        log_and_print(f"Error counting reads: {e}", 'ERROR')
        return 0

def check_file_exists(filepath, description):
    """Simple file existence check for directories"""
    if os.path.exists(filepath):
        log_and_print(f"Found existing {description}: {filepath} - Skipping step")
        return True
    else:
        log_and_print(f"{description} not found: {filepath} - Will create")
        return False

def downsample_fastq(sample_name, fastq_r1, fastq_r2, config):
    """Downsample FASTQ files if they exceed maximum size"""
    max_fq_size, downsample_size = get_downsample_size(config)
    
    # Count reads in R1 file
    r1_path = f"{FASTQ_DIR}/{sample_name}/{fastq_r1}"
    r1_read_count = count_fastq_reads(r1_path)
    
    log_and_print(f"R1 read count: {r1_read_count:,}")
    log_and_print(f"Max allowed reads: {max_fq_size:,}")
    
    if r1_read_count > max_fq_size:
        log_and_print(f"Read count exceeds limit. Downsampling to {downsample_size:,} reads")
        
        # Calculate sampling fraction
        sampling_fraction = downsample_size / r1_read_count
        
        # Create backup directory
        backup_dir = f"{FASTQ_DIR}/{sample_name}/backup_original"
        os.makedirs(backup_dir, exist_ok=True)
        log_and_print(f"Created backup directory: {backup_dir}")
        
        # Define file paths
        r2_path = f"{FASTQ_DIR}/{sample_name}/{fastq_r2}"
        r1_backup = f"{backup_dir}/{fastq_r1}"
        r2_backup = f"{backup_dir}/{fastq_r2}"
        r1_downsampled_temp = f"{ANALYSIS_DIR}/{sample_name}_temp_R1.fastq.gz"
        r2_downsampled_temp = f"{ANALYSIS_DIR}/{sample_name}_temp_R2.fastq.gz"
        
        # Backup original files
        log_and_print("Backing up original FASTQ files...")
        if not run_command(
            "Backup R1",
            f"cp {r1_path} {r1_backup}"
        ):
            log_and_print("Failed to backup R1. Pipeline terminated.", 'ERROR')
            return False
            
        if not run_command(
            "Backup R2", 
            f"cp {r2_path} {r2_backup}"
        ):
            log_and_print("Failed to backup R2. Pipeline terminated.", 'ERROR')
            return False
        
        log_and_print(f"Original files backed up to: {backup_dir}")
        
        # Create downsampled files with temporary names
        log_and_print("Creating downsampled FASTQ files...")
        if not run_command(
            "Downsample R1",
            f"seqtk sample -s100 {r1_path} {sampling_fraction} | gzip > {r1_downsampled_temp}"
        ):
            log_and_print("Failed to downsample R1. Pipeline terminated.", 'ERROR')
            return False
        
        if not run_command(
            "Downsample R2",
            f"seqtk sample -s100 {r2_path} {sampling_fraction} | gzip > {r2_downsampled_temp}"
        ):
            log_and_print("Failed to downsample R2. Pipeline terminated.", 'ERROR')
            return False
        
        # Replace original files with downsampled versions
        log_and_print("Replacing original files with downsampled versions...")
        if not run_command(
            "Replace R1 with downsampled",
            f"mv {r1_downsampled_temp} {r1_path}"
        ):
            log_and_print("Failed to replace R1 with downsampled version. Pipeline terminated.", 'ERROR')
            return False
            
        if not run_command(
            "Replace R2 with downsampled", 
            f"mv {r2_downsampled_temp} {r2_path}"
        ):
            log_and_print("Failed to replace R2 with downsampled version. Pipeline terminated.", 'ERROR')
            return False
        
        # Verify downsampled file read counts
        new_r1_count = count_fastq_reads(r1_path)
        log_and_print(f"Downsampling completed successfully!")
        log_and_print(f"Original R1 reads: {r1_read_count:,}")
        log_and_print(f"Downsampled R1 reads: {new_r1_count:,}")
        log_and_print(f"Sampling fraction: {sampling_fraction:.4f}")
        log_and_print(f"Original files backed up to: {backup_dir}")
        
        # Log downsampling to sample log
        sample_log_file = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.log"
        with open(sample_log_file, "a") as log_file:
            log_file.write(f"DOWNSAMPLING PERFORMED:\n")
            log_file.write(f"  Original reads: {r1_read_count:,}\n")
            log_file.write(f"  Downsampled reads: {new_r1_count:,}\n")
            log_file.write(f"  Sampling fraction: {sampling_fraction:.4f}\n")
            log_file.write(f"  Original files backed up to: {backup_dir}\n")
            log_file.write("="*50 + "\n")
        
        return True#fastq_r1, fastq_r2  # Return original filenames (now containing downsampled data)
    else:
        log_and_print("Read count within limits. No downsampling needed.")
        return True

def get_picard_command(picard, picard_memory):
    """Get the appropriate PICARD command based on the environment"""
    if picard.endswith('.jar'):
        # JAR 파일로 직접 실행 - tmp 디렉토리를 먼저 생성
        return f"mkdir -p /tmp && java -Xmx{picard_memory} -Djava.io.tmpdir=/tmp -jar {picard}"
    else:
        # wrapper script 사용하는 경우도 JAR 직접 실행으로 대체
        return f"mkdir -p /tmp && java -Xmx{picard_memory} -Djava.io.tmpdir=/tmp -jar /Work/NIPT/bin/picard/picard.jar"

def generate_proper_paired_bam(sample_name, fastq_r1, fastq_r2, config, progress, base_step):
    """Generate proper_paired.bam file from FASTQ files"""
    # Define command paths from environment variables with defaults
    bwa = os.environ.get('BWA2', 'bwa-mem2')
    sam_tools = os.environ.get('SAMTools', 'samtools')
    picard = os.environ.get('PICARD', 'picard')
    qualimap = os.environ.get('qualimap', 'qualimap')

    # Get config value order : Environment parameter --> Config --> Default value
    try:
        bwa_threads = os.environ.get('QC.bwa_threads') or str(config['QC']['bwa_threads'])
    except (KeyError, TypeError):
        bwa_threads = '16'

    try:
        samtools_threads = os.environ.get('QC.samtools_threads') or str(config['QC']['samtools_threads'])
    except (KeyError, TypeError):
        samtools_threads = '8'

    try:
        samtools_memory = os.environ.get('QC.samtools_memory') or config['QC']['samtools_memory']
    except (KeyError, TypeError):
        samtools_memory = '4G'

    try:
        picard_memory = os.environ.get('QC.picard_memory') or config['QC']['picard_memory']
    except (KeyError, TypeError):
        picard_memory = '12G'

    log_and_print(f"BWA threads: {bwa_threads}")
    log_and_print(f"Samtools threads: {samtools_threads}")
    log_and_print(f"Samtools memory: {samtools_memory}")
    log_and_print(f"Picard memory: {picard_memory}")

    # Define file paths
    sorted_bam = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.sorted.bam"
    sorted_bam_index = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.sorted.bam.bai"
    qualimap_output = f"{ANALYSIS_DIR}/{sample_name}/Output_QC"
    qualimap_zip = f"{ANALYSIS_DIR}/{sample_name}/Output_QC/{sample_name}.Qualimap.zip"
    dedup_bam = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.dedup.bam"
    dedup_bam_index = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.dedup.bam.bai"
    uniq_bam = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.uniq.bam"
    uniq_bam_index = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.uniq.bam.bai"
    proper_paired_bam = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.proper_paired.bam"
    proper_paired_bam_index = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.proper_paired.bam.bai"
    qc_filter_result = f"{ANALYSIS_DIR}/{sample_name}/Output_QC/{sample_name}.qc.filter.txt"

    ref_hg = f"{DATA_DIR}/refs/common/hg19/ucsc.hg19.fasta"

    # 1. BWA MEM alignment and sorting
    if check_file_exists_advanced(sorted_bam, "sorted BAM", "bam"):
        progress.update_step(f"{base_step}.1", "BWA Alignment", "SKIP", "file exists")
    else:
        if run_command(
            "BWA-MEM2 alignment and sorting",
            f"{bwa} mem -t {bwa_threads} {ref_hg} {ANALYSIS_DIR}/{fastq_r1} {ANALYSIS_DIR}/{fastq_r2} | "
            f"{sam_tools} sort -@ {samtools_threads} -m {samtools_memory} -O bam "
            f"-o {sorted_bam} "
            f"-T {ANALYSIS_DIR}/{sample_name}/{sample_name}.sorted"
        ):
            progress.update_step(f"{base_step}.1", "BWA Alignment", "PASS")
        else:
            progress.update_step(f"{base_step}.1", "BWA Alignment", "FAIL")
            log_and_print("BWA-MEM2 alignment and sorting failed. Pipeline terminated.", 'ERROR')
            return False

    # 2. Index sorted BAM
    if check_file_exists_advanced(sorted_bam_index, "sorted BAM index", "bam_index"):
        progress.update_step(f"{base_step}.2", "sorted BAM index", "SKIP", "file exists")
    else:
        if run_command(
            "SAMTools Indexing",
            f"{sam_tools} index {sorted_bam}"
        ):
            progress.update_step(f"{base_step}.2", "sorted BAM index", "PASS")
        else:
            progress.update_step(f"{base_step}.2", "sorted BAM index", "FAIL")
            log_and_print("SAMTools Indexing failed. Pipeline terminated.", 'ERROR')
            return False

    # 3. Run Qualimap
    if check_file_exists_advanced(qualimap_zip, "Qualimap ZIP file", "zip"):
        progress.update_step(f"{base_step}.3", "Qualiimap execution", "SKIP", "file exists")
        progress.update_step(f"{base_step}.4", "Qualiimap zipping", "SKIP", "file exists")
    else:
        if run_command(
            "Qualimap",
            f"{qualimap} bamqc -bam {sorted_bam} "
            f"-nt 10 -outdir {qualimap_output}"
        ):
            progress.update_step(f"{base_step}.3", "Qualimap", "PASS")
        else:
            progress.update_step(f"{base_step}.3", "Qualimap", "FAIL")
            log_and_print("Qualimap failed. Pipeline terminated.", 'ERROR')
            return False

        if run_command(
            "Zip Qualimap results",
            f"zip -r {qualimap_zip} {qualimap_output}"
        ):
            progress.update_step(f"{base_step}.4", "Qualimap ZIP", "PASS")
        else:
            progress.update_step(f"{base_step}.4", "Qualimap ZIP", "FAIL")
            log_and_print("Zip Qualimap results failed. Pipeline terminated.", 'ERROR')
            return False


    # 4.1 QC filter
    if not check_file_exists(qc_filter_result, "QC filter result"):
        try:
            qc_success, pass_fail_success = process_qc_pipeline(
                sample_name,
                config
            )
            log_and_print(f"process_qc_pipeline completed: qc_success={qc_success}, pass_fail_success={pass_fail_success}")
            if qc_success and pass_fail_success:
                progress.update_step("{base_step}.5", "QC filter", "PASS")
                log_and_print("QC filter: PASS")
            else:
                progress.update_step("{base_step}.5.1", "QC filter", "FAIL")
                log_and_print(f"QC filter: FAIL (qc_success={qc_success}, pass_fail_success={pass_fail_success})", "ERROR")

        except Exception as e:
            #logger = logging.getLogger(__name__)
            log_and_print(f"=== QC PIPELINE EXCEPTION ===", "ERROR")
            log_and_print(f"Exception type: {type(e).__name__}", "ERROR")
            log_and_print(f"Exception message: {str(e)}", "ERROR")

            # 전체 traceback 출력
            import traceback
            log_and_print(f"Full traceback:", "ERROR")
            for line in traceback.format_exc().split('\n'):
                if line.strip():
                    log_and_print(f"  {line}", "ERROR")

            progress.update_step(f"{base_step}.5", "QC filter", "FAIL")
            log_and_print("QC filter: FAIL (due to exception)", "ERROR")
    else:
        progress.update_step(f"{base_step}.5", "QC filter", "SKIP", "file exists")

    # 5. PICARD - Remove duplications
    if check_file_exists_advanced(dedup_bam, "deduplicated BAM", "bam"):
        progress.update_step(f"{base_step}.6", "Picard Dedup", "SKIP", "file exists")
    else:
        picard_cmd = get_picard_command(picard, picard_memory)

        if run_command(
            "PICARD MarkDuplicates",
            f"{picard_cmd} MarkDuplicates "
            f"I={sorted_bam} REMOVE_DUPLICATES=true "
            f"VALIDATION_STRINGENCY=LENIENT AS=true "
            f"M=\"{ANALYSIS_DIR}/{sample_name}/{sample_name}\"_dup.metrics "
            f"O=\"{dedup_bam}\""
        ):
            progress.update_step(f"{base_step}.6", "PICARD Dedup", "PASS")
        else:
            progress.update_step(f"{base_step}.6", "PICARD Dedup", "FAIL")
            log_and_print("PICARD MarkDuplicates failed. Pipeline terminated.", 'ERROR')
            return False

    # 6. Index dedup BAM
    if check_file_exists_advanced(dedup_bam_index, "deduplicated BAM index", "bam_index"):
        progress.update_step(f"{base_step}.7", "Picard Dedup index", "SKIP", "file exists")
    else:
        if run_command(
            "PICARD bam indexing",
            f"{sam_tools} index {dedup_bam}"
        ):
            progress.update_step(f"{base_step}.7", "PICARD Dedup index", "PASS")
        else:
            progress.update_step(f"{base_step}.7", "PICARD Dedup index", "FAIL")
            log_and_print("PICARD bam indexing failed. Pipeline terminated.", 'ERROR')
            return False

    # 7. SAMtools - Make unique BAM
    if check_file_exists_advanced(uniq_bam, "unique BAM", "bam"):
        progress.update_step(f"{base_step}.8", "Unique BAM", "SKIP", "file exists")
    else:
        if run_command(
            "unique bam",
            f"{sam_tools} view -bq 1 {dedup_bam} > {uniq_bam}"
        ):
            progress.update_step(f"{base_step}.8", "Unique BAM", "PASS")
        else:
            progress.update_step(f"{base_step}.8", "Unique BAM", "FAIL")
            log_and_print("Unique bam creation failed. Pipeline terminated.", 'ERROR')
            return False

    # 8. Index unique BAM
    if check_file_exists_advanced(uniq_bam_index, "unique BAM index", "bam_index"):
        progress.update_step(f"{base_step}.9", "Unique BAM index", "SKIP", "file exists")
    else:
        if run_command(
            "unique bam indexing",
            f"{sam_tools} index {uniq_bam}"
        ):
            progress.update_step(f"{base_step}.9", "Unique BAM index", "PASS")
        else:
            progress.update_step(f"{base_step}.9", "Unique BAM index", "FAIL")
            log_and_print("Unique bam indexing failed. Pipeline terminated.", 'ERROR')
            return False

    # 9. Extract Proper paired BAM
    if check_file_exists_advanced(proper_paired_bam, "proper paired BAM", "bam"):
        progress.update_step(f"{base_step}.10", "Proper paired BAM", "SKIP", "file exists")
    else:
        if run_command(
            "proper bam",
            f"{sam_tools} view -b -f 0x2 {uniq_bam} > {proper_paired_bam}"
        ):
            progress.update_step(f"{base_step}.10", "Proper paired BAM", "PASS")
        else:
            progress.update_step(f"{base_step}.10", "Proper paired BAM", "FAIL")
            log_and_print("Proper paired bam creation failed. Pipeline terminated.", 'ERROR')
            return False

    # 10. Index proper paired BAM
    if check_file_exists_advanced(proper_paired_bam_index, "proper paired BAM index", "bam_index"):
        progress.update_step(f"{base_step}.11", "Proper paired BAM index", "SKIP", "file exists")
    else:
        if run_command(
            "proper bam indexing",
            f"{sam_tools} index {proper_paired_bam}"
        ):
            progress.update_step(f"{base_step}.11", "Proper paired BAM index", "PASS")
        else:
            progress.update_step(f"{base_step}.11", "Proper paired BAM index", "FAIL")
            log_and_print("Proper paired bam indexing failed. Pipeline terminated.", 'ERROR')
            return False

    log_and_print(f"{proper_paired_bam} generation completed")
    return True

def process_filter(sample_name, filter_type, filter_path):
    """Process a specific filter type and create all size variants"""
    # Define command paths with defaults
    sam_tools = os.environ.get('SAMTools', 'samtools')
    ref_hg = f"{DATA_DIR}/refs/common/hg19/ucsc.hg19.fasta"

    # Define file paths
    proper_paired_bam = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.proper_paired.bam"
    orig_bam = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_type}_orig.bam"
    orig_bam_index = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_type}_orig.bam.bai"
    fetus_bam = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_type}_fetus.bam"
    fetus_bam_index = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_type}_fetus.bam.bai"
    mom_bam = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_type}_mom.bam"
    mom_bam_index = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_type}_mom.bam.bai"

    # Create orig BAM from proper_paired.bam
    if not check_file_exists_advanced(orig_bam, f"{filter_type}_orig BAM", "bam"):
        log_and_print(f"SAMTools {filter_type}_orig bam start")

        # Debug: Check if files exist before running command
        log_and_print(f"Checking BED file: {filter_path}")
        if not os.path.exists(filter_path):
            log_and_print(f"BED file not found: {filter_path}", 'ERROR')
            return False

        if not os.path.exists(proper_paired_bam):
            log_and_print(f"BAM file not found: {proper_paired_bam}", 'ERROR')
            return False

        if not run_command(
            f"{filter_type}_orig bam",
            f"cd /Work/NIPT && {sam_tools} view -b -L {filter_path} {proper_paired_bam} > {orig_bam}"
        ):
            log_and_print(f"{filter_type}_orig bam creation failed. Pipeline terminated.", 'ERROR')
            return False

    # Index orig BAM
    if not check_file_exists_advanced(orig_bam_index, f"{filter_type}_orig BAM index", "bam_index"):
        if not run_command(
            f"{filter_type}_orig bam indexing",
            f"{sam_tools} index {orig_bam}"
        ):
            log_and_print(f"{filter_type}_orig bam indexing failed. Pipeline terminated.", 'ERROR')
            return False

    log_and_print(f"SAMTools {filter_type}_orig bam end")

    # Create fetus BAM (size filter using fetus.awk)
    if not check_file_exists_advanced(fetus_bam, f"{filter_type}_fetus BAM", "bam"):
        log_and_print(f"SAMTools {filter_type}_fetus bam start")
        if not run_command(
            f"{filter_type}_fetus bam",
            f"{sam_tools} view {orig_bam} | "
            f"awk -f {SRC_DIR}/fetus.awk | {sam_tools} view -bt {ref_hg}.fai -o {fetus_bam} -"
        ):
            log_and_print(f"{filter_type}_fetus bam creation failed. Pipeline terminated.", 'ERROR')
            return False

    # Index fetus BAM
    if not check_file_exists_advanced(fetus_bam_index, f"{filter_type}_fetus BAM index", "bam_index"):
        if not run_command(
            f"{filter_type}_fetus bam indexing",
            f"{sam_tools} index {fetus_bam}"
        ):
            log_and_print(f"{filter_type}_fetus bam indexing failed. Pipeline terminated.", 'ERROR')
            return False

    log_and_print(f"SAMTools {filter_type}_fetus bam end")

    # Create mom BAM (size filter using mom.awk)
    if not check_file_exists_advanced(mom_bam, f"{filter_type}_mom BAM", "bam"):
        log_and_print(f"SAMTools {filter_type}_mom bam start")
        if not run_command(
            f"{filter_type}_mom bam",
            f"{sam_tools} view {orig_bam} | "
            f"awk -f {SRC_DIR}/mom.awk | {sam_tools} view -bt {ref_hg}.fai -o {mom_bam} -"
        ):
            log_and_print(f"{filter_type}_mom bam creation failed. Pipeline terminated.", 'ERROR')
            return False

    # Index mom BAM
    if not check_file_exists_advanced(mom_bam_index, f"{filter_type}_mom BAM index", "bam_index"):
        if not run_command(
            f"{filter_type}_mom bam indexing",
            f"{sam_tools} index {mom_bam}"
        ):
            log_and_print(f"{filter_type}_mom bam indexing failed. Pipeline terminated.", 'ERROR')
            return False

    log_and_print(f"SAMTools {filter_type}_mom bam end")

def create_npz_files(sample_name, bam_file, bam_suffix, config):
    """Create WC, WCX, and WFF NPZ files for a given BAM file using lab-specific references"""
    # Define command paths with defaults
    python2 = os.environ.get('PYTHON2', 'python2.7')
    wc_path = os.environ.get('WC', '/opt/wisecondor/wisecondor.py')
    wcx_path = os.environ.get('WCX', 'wisecondorx')
    wcff_path = os.environ.get('WCFF', 'wisecondor-ff')

    # Get lab-specific references
    lab_references = config.get('lab_references', {})

    # Define output file paths
    wc_output = f"{ANALYSIS_DIR}/{sample_name}/Output_WC/{sample_name}.wc.{bam_suffix}.npz"
    wcx_output = f"{ANALYSIS_DIR}/{sample_name}/Output_WCX/{sample_name}.wcx.{bam_suffix}.npz"
    wcff_output = f"{ANALYSIS_DIR}/{sample_name}/Output_WCFF/{sample_name}.wcff.{bam_suffix}.npz"

    # Get binsize from config
    wc_binsize = config.get('Wisecondor', {}).get('binsize', 200000)
    wcx_binsize = config.get('WisecondorX', {}).get('binsize', 200000)
    wcff_binsize = config.get('WisecondorFF', {}).get('binsize', 200000)

    # WC convert
    if not check_file_exists_advanced(wc_output, f"WC NPZ file for {bam_suffix}", "npz"):
        if not run_command(
            f"WC convert {bam_suffix}",
            f"{python2} {wc_path} convert {bam_file} {wc_output} -binsize {wc_binsize}"
        ):
            log_and_print(f"WC convert {bam_suffix} failed. Pipeline terminated.", 'ERROR')
            return False

    # WCX convert
    if not check_file_exists_advanced(wcx_output, f"WCX NPZ file for {bam_suffix}", "npz"):
        if not run_command(
            f"WCX convert {bam_suffix}",
            f"{wcx_path} convert {bam_file} {wcx_output} --binsize {wcx_binsize}"
        ):
            log_and_print(f"WCX convert {bam_suffix} failed. Pipeline terminated.", 'ERROR')
            return False

    # WCFF convert
    if not check_file_exists_advanced(wcff_output, f"WCFF NPZ file for {bam_suffix}", "npz"):
        if not run_command(
            f"WCFF convert {bam_suffix}",
            f"{wcff_path} convert -i {bam_file} -o {wcff_output} -b {wcff_binsize}"
        ):
            log_and_print(f"WCFF convert {bam_suffix} failed. Pipeline terminated.", 'ERROR')
            return False

def create_hmmcopy_files(sample_name, bam_file, bam_suffix):
    """Create 50kb and 10mb wig files and run HMMcopy for a given BAM file"""
    # Define command paths with defaults
    hmmcopy_path = os.environ.get('HMMcopy', '/opt/conda/envs/nipt')
    rscript = os.environ.get('Rscript', 'Rscript')

    # Common chromosomes list
    chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

    # Ensure Output_hmmcopy directory exists
    os.makedirs(f"{ANALYSIS_DIR}/{sample_name}/Output_hmmcopy", exist_ok=True)

    # 50kb wig file
    wig_50kb = f"{ANALYSIS_DIR}/{sample_name}/Output_hmmcopy/{sample_name}.{bam_suffix}.50kb.wig"
    if not check_file_exists_advanced(wig_50kb, f"50kb WIG file for {bam_suffix}", "wig"):
        if not run_command(
            f"readCounter 50kb {bam_suffix}",
            f"{hmmcopy_path}/bin/readCounter -w 50000 -c {chromosomes} {bam_file} > {wig_50kb}"
        ):
            log_and_print(f"readCounter 50kb {bam_suffix} failed. Pipeline terminated.", 'ERROR')
            return False

    # Run HMMcopy R script for 50kb and save to Output_hmmcopy
    hmmcopy_50kb_output = f"{ANALYSIS_DIR}/{sample_name}/Output_hmmcopy/{sample_name}.{bam_suffix}.50kb.wig.Normalization.txt"
    if not check_file_exists(hmmcopy_50kb_output, f"HMMcopy 50kb output for {bam_suffix}"):
        if not run_command(
            f"HMMcopy R 50kb {bam_suffix}",
            f"{rscript} --no-save --no-restore {SRC_DIR}/HMMcopy.R {wig_50kb} "
            f"{hmmcopy_path}/share/hmmcopy/hg19.gc.50kb.wig {hmmcopy_path}/share/hmmcopy/hg19.map.50kb.wig 50kb "
            f"{hmmcopy_50kb_output}"
        ):
            log_and_print(f"HMMcopy R 50kb {bam_suffix} failed. Pipeline terminated.", 'ERROR')
            return False

    # 10mb wig file
    wig_10mb = f"{ANALYSIS_DIR}/{sample_name}/Output_hmmcopy/{sample_name}.{bam_suffix}.10mb.wig"
    if not check_file_exists_advanced(wig_10mb, f"10mb WIG file for {bam_suffix}", "wig"):
        if not run_command(
            f"readCounter 10mb {bam_suffix}",
            f"{hmmcopy_path}/bin/readCounter -w 10000000 -c {chromosomes} {bam_file} > {wig_10mb}"
        ):
            log_and_print(f"readCounter 10mb {bam_suffix} failed. Pipeline terminated.", 'ERROR')
            return False

    # Run HMMcopy R script for 10mb and save to Output_hmmcopy
    hmmcopy_10mb_output = f"{ANALYSIS_DIR}/{sample_name}/Output_hmmcopy/{sample_name}.{bam_suffix}.10mb.wig.Normalization.txt"
    if not check_file_exists(hmmcopy_10mb_output, f"HMMcopy 10mb output for {bam_suffix}"):
        if not run_command(
            f"HMMcopy R 10mb {bam_suffix}",
            f"{rscript} --no-save --no-restore {SRC_DIR}/HMMcopy.R {wig_10mb} "
            f"{hmmcopy_path}/share/hmmcopy/hg19.gc.10mb.wig {hmmcopy_path}/share/hmmcopy/hg19.map.10mb.wig 10mb "
            f"{hmmcopy_10mb_output}"
        ):
            log_and_print(f"HMMcopy R 10mb {bam_suffix} failed. Pipeline terminated.", 'ERROR')
            return False

    return True

def count_reads_in_regions(bam_file, bed_file):
    """BED 파일에 정의된 영역에서 reads 카운트"""
    try:
        bamfile = pysam.AlignmentFile(bam_file, "rb")
        regions = []
        counts = []
        total_bases = 0
        
        # BED 파일 읽기
        with open(bed_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                        regions.append((chrom, start, end))
                        total_bases += (end - start)
        
        # 각 영역에서 reads 카운트
        for chrom, start, end in regions:
            count = bamfile.count(chrom, start, end)
            counts.append(count)
        
        bamfile.close()
        return counts, regions, total_bases
        
    except Exception as e:
        log_and_print(f"Error counting reads in regions: {e}", "ERROR")
        return [], [], 0

def calculate_normalized_coverage(counts, regions):
    """Calculate normalized coverage for regions"""
    normalized_counts = []

    try:
        normalized = []
        for i, (chrom, start, end) in enumerate(regions):
            region_length = end - start
            if region_length > 0:
                coverage = counts[i] / region_length
                normalized.append(coverage)
        return normalized
    except Exception as e:
        log_and_print(f"Error calculating normalized coverage: {e}", "ERROR")
        return []

def calculate_yff(sample_name, bam_file, ff_config, lab_bed_paths):
    """Calculate Y-based fetal fraction using lab-specific BED files"""

    log_and_print(f"Calculating Y-based fetal fraction for {sample_name}")

    # Define BED file paths from lab-specific paths
    bed_dir = Path("/Work/NIPT/data/bed/common")
    y_bed = bed_dir / "chrY_target_0.9.bed"
    a_bed = bed_dir / "autosome_control.bed"

    # Check if BED files exist
    if not os.path.exists(y_bed):
        print(f"Y chromosome BED file not found: {y_bed}")
        return {
            'sample_name': sample_name,
            'YFF1': 0,
            'gd_1_value': 0,
            'gd_1_gender': 'UNKNOWN',
            'status': 'FAILED'
        }

    if not os.path.exists(a_bed):
        print(f"Autosome control BED file not found: {a_bed}")
        return {
            'sample_name': sample_name,
            'YFF1': 0,
            'gd_1_value': 0,
            'gd_1_gender': 'UNKNOWN',
            'status': 'FAILED'
        }

    try:
        # Count reads in Y chromosome regions
        log_and_print("Counting reads in Y chromosome regions...")
        y_counts, y_regions, y_total_bases = count_reads_in_regions(bam_file, str(y_bed))

        # Count reads in autosomal regions
        log_and_print("Counting reads in autosomal control regions...")
        a_counts, a_regions, a_total_bases = count_reads_in_regions(bam_file, str(a_bed))

        if not y_counts or not a_counts:
            log_and_print("Failed to count reads in regions", 'WARNING')
            return {
                'sample_name': sample_name,
                'YFF1': 0,
                'gd_1_value': 0,
                'gd_1_gender': 'UNKNOWN',
                'status': 'FAILED'
            }

        # Calculate normalized coverage
        y_normalized = calculate_normalized_coverage(y_counts, y_regions)
        a_normalized = calculate_normalized_coverage(a_counts, a_regions)

        if not y_normalized or not a_normalized:
            log_and_print("Failed to calculate normalized coverage", 'WARNING')
            return {
                'sample_name': sample_name,
                'YFF1': 0,
                'gd_1_value': 0,
                'gd_1_gender': 'UNKNOWN',
                'status': 'FAILED'
            }

        # Calculate median coverage
        y_median_coverage = np.median(y_normalized)
        a_median_coverage = np.median(a_normalized)

        log_and_print(f"Median Y coverage: {y_median_coverage:.6f}")
        log_and_print(f"Median Autosome coverage: {a_median_coverage:.6f}")

        # Check if we have reasonable values to proceed
        if a_median_coverage <= 0:
            log_and_print("Autosomal median coverage is zero or negative. Cannot calculate YFF.", 'WARNING')
            return {
                'sample_name': sample_name,
                'YFF1': 0,
                'gd_1_value': 0,
                'gd_1_gender': 'UNKNOWN',
                'status': 'FAILED',
            }

        # Calculate the ratio of Y to autosome
        y_to_a_ratio = y_median_coverage / a_median_coverage

        # Calculate fetal fraction assuming a male fetus
        # For a male fetus in maternal blood, the formula is:
        # FF = 2 * (Y/A ratio) since only the male fetus contributes Y chromosome
        fetal_fraction = 2 * y_to_a_ratio * 100  # Convert to percentage

        # Determine gender based on Y coverage
        gender_threshold = ff_config.get('gd_1_threshold', 0.01)

        log_and_print(f"gender_threshold = {gender_threshold}")
        if y_to_a_ratio < gender_threshold:
            gd_1_gender = "XX"
            gd_1_value = y_to_a_ratio
            status = "OK"
        else:
            gd_1_gender = "XY"
            gd_1_value = y_to_a_ratio
            status = "OK"

        return {
            'sample_name': sample_name,
            'YFF1': fetal_fraction,
            'gd_1_value': gd_1_value,
            'gd_1_gender': gd_1_gender,
            'status': "OK"
        }

    except Exception as e:
        print(f"Error in calculate_yff: {e}")
        return {
            'sample_name': sample_name,
            'YFF1': 0,
            'gd_1_value': 0,
            'gd_1_gender': 'UNKNOWN',
            'status': f'ERROR: {str(e)}',
        }

def calculate_yff2(wig_norm_file, ff_config, paths):
    """Calculate adjusted Y-based fetal fraction (YFF2) using corrected wig normalization file."""

    if not os.path.exists(wig_norm_file):
        return {
            "YFF2": 0,
            "UAR_X": 0,
            "UAR_Y": 0,
            "gd_2_value": 0.0,
            "gd_2_gender": "XX",
            "gd_3_value": 0.0,
            "gd_3_gender": "XX",
            "status": "FAILED"
        }

    try:
        # Load normalization wig file
        df = pd.read_csv(
            wig_norm_file,
            sep='\t',
            header=0,
            names=["chr", "start", "end", "reads", "gc", "map", "valid", "ideal", "cor.gc", "cor.map", "copy"]
        )

        # Strip spaces and convert types
        df["chr"] = df["chr"].str.strip()
        for col in ["start", "end", "reads"]:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)
        df["cor.gc"] = pd.to_numeric(df["cor.gc"], errors="coerce").fillna(0)

        # Keep only positive cor.gc values
        df = df[df["cor.gc"] > 0]

        # Remove centromere regions
        df_chrX = df[(df["chr"] == "chrX") & ~((df["start"] >= 58100001) & (df["start"] <= 63000000))]
        df_chrY = df[(df["chr"] == "chrY") & ~((df["start"] >= 11600001) & (df["start"] <= 14000000))]

        # ================================
        # GD_2 CALCULATION 
        # ================================

        log_and_print("Calculating gd_2...")

        # 1. unique.Y: Y 염색체에서 특정 영역만 선택 (R script와 동일)
        unique_Y = df_chrY[(df_chrY["start"] > 2650001) & (df_chrY["start"] < 59050000)]

        # 2. unique.Y.specific: 14개 특정 위치의 윈도우만 선택
        specific_positions = [7650001, 7750001, 7800001, 8400001, 8450001, 8500001,
                            8550001, 8600001, 15500001, 18900001, 22250001, 22450001,
                            22900001, 23600001]

        unique_Y_specific = unique_Y[unique_Y["start"].isin(specific_positions)]

        logger.info(len(unique_Y))
        logger.info(unique_Y["cor.gc"].sum())
        # 3. gd_2 값 계산: sum(unique.Y.specific[,4])/sum(unique.Y[,4]) * 100
        if len(unique_Y) == 0 or unique_Y["cor.gc"].sum() == 0:
            log_and_print(f"gd_2 FAIL")
            gd_2_value = 0.0
        else:
            sum_unique_Y_specific = unique_Y_specific["cor.gc"].sum()
            sum_unique_Y = unique_Y["cor.gc"].sum()
            gd_2_value = (sum_unique_Y_specific / sum_unique_Y) * 100
            logger.info(gd_2_value)

        # gd_2 gender detection
        gd_2_threshold = ff_config.get('gd_2_threshold', 0.4)
        gd_2_gender = "XY" if gd_2_value > gd_2_threshold else "XX"

        logger.info(f"gd_2: {gd_2_value:.6f} → {gd_2_gender} (threshold: {gd_2_threshold})")
        logger.info(f"  unique_Y windows: {len(unique_Y)}, unique_Y_specific matched: {len(unique_Y_specific)}")

        # ================================
        # YFF2 Calculation
        # ================================

        log_and_print("Calculating YFF2...")

        # Autosomal sum (chr1 to chr22)
        autosome_sum = 0
        for i in range(1, 23):
            chr_df = df[df["chr"] == f"chr{i}"]
            autosome_sum += chr_df["cor.gc"].sum()

        logger.info(f"autosome_sum : {autosome_sum}")
        if autosome_sum <= 0:
            return {
                "YFF2": 0,
                "UAR_X": 0,
                "UAR_Y": 0,
                "gd_2_value": gd_2_value,
                "gd_2_gender": gd_2_gender,
                "gd_3_value": 0.0,
                "gd_3_gender": "XX",
                "status": "FAILED"
            }

        # Usage ratios
        UAR_Y = df_chrY["cor.gc"].sum() / autosome_sum * 100
        UAR_X = df_chrX["cor.gc"].sum() / autosome_sum * 100

        # Estimate FF
        base = 0.0199917459497307
        max_val = 0.0468805942976191
        FF_chrY = (UAR_Y - base) / (max_val - base)
        FF_chrY_1 = round(FF_chrY, 3) * 10
        FF_chrY_2 = FF_chrY_1 + 0.55

        logger.info(f"UAR_Y : {UAR_Y}, UAR_X : {UAR_X}")

        # ================================
        # GD_3 CALCULATION
        # ================================

        print("Calculating gd_3 (UAR_Y based)...")

        # gd_3 gender detection based on UAR_Y
        gd_3_threshold = ff_config.get('gd_3_threshold', 0.02)

        if UAR_Y > gd_3_threshold:
            gd_3_gender = "XY"  # MALE
            gd_3_value = FF_chrY_2  # Use actual FF value for males
        else:
            gd_3_gender = "XX"  # FEMALE
            gd_3_value = 0.0  # Set FF to 0 for female fetuses

        print(f"gd_4: {gd_3_value:.6f} → {gd_3_gender} (UAR_Y: {UAR_Y:.6f}, threshold: {gd_3_threshold})")

        return {
            # YFF2 원본 결과
            "YFF2": FF_chrY_2,
            "UAR_X": UAR_X,
            "UAR_Y": UAR_Y,
            "status": "OK",
            "gd_2_value": gd_2_value,
            "gd_2_gender": gd_2_gender,
            "gd_3_value": gd_3_value,
            "gd_3_gender": gd_3_gender,
            "status": "OK"
        }


    except Exception as e:
        return {
            "YFF2": 0,
            "UAR_X": 0,
            "UAR_Y": 0,
            "gd_2_value": 0.0,
            "gd_2_gender": "XX",
            "gd_3_value": 0.0,
            "gd_3_gender": "XX",
            "status": f"ERROR: {str(e)}"
        }

def extract_fragment_features(bam_file):
    """Extract fragment size and GC content features from BAM file"""
    try:
        sam_tools = os.environ.get('SAMTools', 'samtools')
        
        # Sample a subset of reads for efficiency
        sample_size = 100000
        
        # Extract fragment sizes
        cmd = f"{sam_tools} view -f 0x2 {bam_file} | head -n {sample_size} | cut -f9 | awk '{{if($1>0) print $1}}'"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            log_and_print(f"Failed to extract fragment sizes: {result.stderr}", 'WARNING')
            return None
        
        fragment_sizes = []
        for line in result.stdout.strip().split('\n'):
            if line.strip():
                try:
                    size = int(line.strip())
                    if 50 <= size <= 1000:  # Filter reasonable fragment sizes
                        fragment_sizes.append(size)
                except ValueError:
                    continue
        
        if not fragment_sizes:
            log_and_print("No valid fragment sizes found", 'WARNING')
            return None
        
        # Calculate fragment size statistics
        fragment_sizes = np.array(fragment_sizes)
        median_fragment_size = np.median(fragment_sizes)
        fragment_size_std = np.std(fragment_sizes)
        
        # Extract GC content from sequences (simplified approach)
        cmd = f"{sam_tools} view -f 0x2 {bam_file} | head -n {sample_size} | cut -f10"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        gc_contents = []
        if result.returncode == 0:
            for line in result.stdout.strip().split('\n'):
                if line.strip():
                    seq = line.strip()
                    if seq and seq != '*':
                        gc_count = seq.count('G') + seq.count('C')
                        gc_content = gc_count / len(seq) if len(seq) > 0 else 0
                        gc_contents.append(gc_content)
        
        if gc_contents:
            median_gc_content = np.median(gc_contents)
            gc_content_std = np.std(gc_contents)
        else:
            median_gc_content = 0.5  # Default value
            gc_content_std = 0.1
        
        return {
            'median_fragment_size': median_fragment_size,
            'fragment_size_std': fragment_size_std,
            'median_gc_content': median_gc_content,
            'gc_content_std': gc_content_std
        }
        
    except Exception as e:
        log_and_print(f"Error extracting fragment features: {e}", 'ERROR')
        return None


def run_fastqc(sample_name, fastq_r1, fastq_r2):
    output_dir = f"{ANALYSIS_DIR}/{sample_name}/Output_QC"
    os.makedirs(output_dir, exist_ok=True)

    fq1_base = os.path.basename(fastq_r1).replace(".gz", "").replace(".fastq", "")
    fq2_base = os.path.basename(fastq_r2).replace(".gz", "").replace(".fastq", "")
    fq1_html = os.path.join(output_dir, f"{fq1_base}_fastqc.html")
    fq2_html = os.path.join(output_dir, f"{fq2_base}_fastqc.html")
    fq1_zip = os.path.join(output_dir, f"{fq1_base}_fastqc.zip")
    fq2_zip = os.path.join(output_dir, f"{fq2_base}_fastqc.zip")

    if all(os.path.isfile(f) for f in [fq1_html, fq1_zip, fq2_html, fq2_zip]):
        log_and_print("FastQC output already exists. Skipping FastQC step.")
        return True

    input_dir = f"{FASTQ_DIR}/{sample_name}"
    fq1_path = os.path.join(input_dir, fastq_r1)
    fq2_path = os.path.join(input_dir, fastq_r2)

    log_and_print("Running FastQC")
    try:
        subprocess.run(f"fastqc -o {output_dir} {fq1_path} {fq2_path}", shell=True, check=True)
    except subprocess.CalledProcessError:
        log_and_print("FastQC execution failed.", 'ERROR')
        return False

    missing_files = [f for f in [fq1_html, fq1_zip, fq2_html, fq2_zip] if not os.path.isfile(f)]
    if missing_files:
        log_and_print(f"FastQC missing output files: {missing_files}", 'ERROR')
        return False

    log_and_print(f"FastQC completed successfully. Results in: {output_dir}")
    return True

def calculate_fragmentff_simplified(sample_name, bam_file, features, ff_config):
    """Simplified seqFF calculation based on fragment size distribution"""
    log_and_print("Calculating simplified FF based on fragment size distribution")
    
    # Use fragment size distribution to estimate fetal fraction
    median_fragment_size = features['median_fragment_size']
    fragment_size_std = features['fragment_size_std']
    
    # Simplified formula based on fragment size deviation from maternal baseline
    maternal_baseline = 195  # Typical maternal fragment size
    fetal_optimal = 167     # Typical fetal fragment size
    
    # Calculate deviation score
    if median_fragment_size < maternal_baseline:
        # Shorter fragments suggest fetal contribution
        deviation_score = (maternal_baseline - median_fragment_size) / (maternal_baseline - fetal_optimal)
        fragmentff_value = min(max(deviation_score * 15, 0), 30)  # Scale to 0-30%
    else:
        # Longer fragments suggest minimal fetal contribution
        fragmentff_value = max(0, 5 - (median_fragment_size - maternal_baseline) / 10)
    
    # Adjust based on fragment size variance
    if fragment_size_std > 40:  # High variance might indicate fetal contribution
        fragmentff_value += min(fragment_size_std / 20, 5)
    
    fragmentff_value = min(fragmentff_value, 40)  # Cap at 40%
    
    log_and_print(f"Simplified seqFF calculation completed: {fragmentff_value:.2f}%")
    
    return {
        'sample_name': sample_name,
        'fragmentff_value': fragmentff_value,
        'median_fragment_size': median_fragment_size,
        'fragment_size_std': fragment_size_std,
        'median_gc_content': features['median_gc_content'],
        'status': 'OK_SIMPLIFIED'
    }

def calculate_fragmentff(sample_name, bam_file, ff_config):
    """Calculate sequence-based fetal fraction"""
    log_and_print(f"Calculating sequence-based fetal fraction for {sample_name}")
    
    # Extract features from BAM file
    features = extract_fragment_features(bam_file)
    if not features:
        log_and_print("Failed to extract features for seqFF calculation", 'WARNING')
        return {}
    
    # Use simplified approach (in production, you would use a trained ML model)
    return calculate_fragmentff_simplified(sample_name, bam_file, features, ff_config)

def determine_final_ff(ff_results, ff_config):
    """Determine the final fetal fraction value from multiple methods"""
    min_threshold = ff_config.get('ff_min_threshold', 2.0)
    max_threshold = ff_config.get('ff_max_threshold', 40.0)
    
    yff_available = 'yff' in ff_results and ff_results['yff'].get('status') not in ['FAILED', 'NA']
    seqff_available = 'seqff' in ff_results and ff_results['seqff'].get('status') != 'FAILED'
    
    # Priority: YFF for males, seqFF for females or when YFF is not available
    if yff_available:
        yff = ff_results['yff']
        if yff['gender'] == 'MALE' and yff['yff_value'] >= min_threshold:
            # Use YFF for male samples
            return {
                'value': yff['yff_value'],
                'method': 'YFF',
                'quality': 'HIGH' if yff['yff_value'] <= max_threshold else 'MODERATE',
                'status': 'OK',
                'note': 'Y-based calculation for male fetus'
            }
    
    if seqff_available:
        seqff = ff_results['seqff']
        if min_threshold <= seqff['seqff_value'] <= max_threshold:
            # Use seqFF
            quality = 'HIGH' if seqff.get('status') == 'OK' else 'MODERATE'
            return {
                'value': seqff['seqff_value'],
                'method': 'seqFF',
                'quality': quality,
                'status': 'OK',
                'note': 'Sequence-based calculation'
            }
    
    # Fallback: try to use any available method with lower quality
    if yff_available and ff_results['yff']['gender'] == 'FEMALE':
        return {
            'value': 0,
            'method': 'YFF',
            'quality': 'LOW',
            'status': 'NA',
            'note': 'Female fetus - YFF not applicable'
        }
    
    if seqff_available:
        return {
            'value': ff_results['seqff']['seqff_value'],
            'method': 'seqFF',
            'quality': 'LOW',
            'status': 'WARNING',
            'note': 'Value outside normal range'
        }
    
    # No reliable method available
    return {
        'value': 0,
        'method': 'NONE',
        'quality': 'FAILED',
        'status': 'FAILED',
        'note': 'No reliable fetal fraction calculation available'
    }

def process_hmmcopy(sample_name, config):
    """Process all BAM files to create NPZ and HMMcopy files using lab-specific settings"""
    
    # List of all BAM files to process with their suffixes
    bam_files = []
    
    # Add proper_paired.bam
    proper_paired_bam = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.proper_paired.bam"
    bam_files.append((proper_paired_bam, "proper_paired"))
    
    # Process each filter type to create BAM files using lab-specific BED files
    # nf09 removes too many reads and it causes WC, WCFF convert error
    #filter_types = ['of', 'nf08', 'nf09']
    filter_types = ['of', 'nf08']
    for filter_type in filter_types:
        # Add orig BAM
        orig_bam = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_type}_orig.bam"
        bam_files.append((orig_bam, f"{filter_type}_orig"))
        
        # Add fetus BAM
        fetus_bam = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_type}_fetus.bam"
        bam_files.append((fetus_bam, f"{filter_type}_fetus"))
        
        # Add mom BAM
        mom_bam = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_type}_mom.bam"
        bam_files.append((mom_bam, f"{filter_type}_mom"))
    
    # Process each BAM file
    for bam_file, bam_suffix in bam_files:
        log_and_print(f"Processing {bam_suffix} BAM file")
        
        # Create NPZ files for all BAM files
        create_npz_files(sample_name, bam_file, bam_suffix, config)
        
        # For all except proper_paired.bam, also create HMMcopy files
        #if bam_suffix != "proper_paired":
        create_hmmcopy_files(sample_name, bam_file, bam_suffix)
    
    log_and_print("All BAM processing and file generation completed")
    return True

def extract_seqff_value(seqff_txt_path):
    """Extract the SeqFF value from output .txt file"""
    if not seqff_txt_path or not os.path.exists(seqff_txt_path):
        return 0
    try:
        with open(seqff_txt_path, 'r') as f:
            for line in f:
                if line.startswith('"SeqFF"'):
                    parts = line.strip().split(',')
                    return float(parts[1])
    except Exception as e:
        log_and_print(f"Failed to extract SeqFF value: {e}", 'WARNING')
    return 0

def calculate_seqff(sample_name, bam_path):
    output_dir = f"{ANALYSIS_DIR}/{sample_name}/Output_FF"
    seqff_txt = os.path.join(output_dir, f"{sample_name}.seqff.txt")

    log_and_print("Running official R-based seqFF script")
    rscript = os.environ.get("Rscript", "Rscript")
    seqff_r_script = "/opt/nipt/bin/scripts/seqFF_R/seqff.r"

    original_cwd = os.getcwd()
    seqff_dir = "/opt/nipt/bin/scripts/seqFF_R"

    cmd = f'{rscript} --vanilla {seqff_r_script} -f {bam_path} -o {seqff_txt}'
    try:
        os.chdir(seqff_dir)
        subprocess.run(cmd, shell=True, check=True)
        log_and_print(f"Official seqFF output written to: {seqff_txt}")
    except subprocess.CalledProcessError:
        log_and_print("R-based seqFF failed to run.", 'ERROR')
        return {"seqff_value": 0, "seqff_file": None}

    # Immediately parse "SeqFF" result from file
    try:
        import pandas as pd
        df = pd.read_csv(seqff_txt, index_col=0)
        seqff_value = float(df.loc["SeqFF", "x"]) * 100
        return {"seqff_value": round(seqff_value, 2), "seqff_file": seqff_txt}
    except Exception as e:
        log_and_print(f"Could not parse seqFF result: {e}", 'WARNING')
        return {"seqff_value": 0, "seqff_file": seqff_txt}

def calculate_fetal_fraction(sample_name, config, paths):

    ff_output_dir = Path(ANALYSIS_DIR) / sample_name / "Output_FF"
    os.makedirs(ff_output_dir, exist_ok=True)
    ff_output_txt = ff_output_dir / f"{sample_name}.fetal_fraction.txt"
    gender_output_txt = ff_output_dir / f"{sample_name}.gender.txt"

    '''
    # 이미 결과가 있으면 생략
    if check_file_exists(ff_json, "Fetal Fraction JSON result"):
        log_and_print(f"Skipping FF calculation for {sample_name} as results already exist.")
        with open(ff_json) as f:
            existing_results = json.load(f)
        log_and_print(f"Fetal fraction results for {sample_name}: {json.dumps(existing_results)}")
        return existing_results
    '''

    bam_path = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.proper_paired.bam"
    ff_results = {}
    gender_results = {}

    # Fragment-based FF
    fragff_result = calculate_fragmentff(sample_name, bam_path, config.get("FF_Gender_Config", {}))
    ff_results['Fragment_FF'] = {
        'value' : round(fragff_result.get('fragmentff_value', 0), 2)
    }

    # YFF (based on autosome vs Y read counts)
    yff_result = calculate_yff(sample_name, bam_path, config.get("FF_Gender_Config", {}), config.get('lab_bed_paths', {}))
    if yff_result['status'] == 'OK':
        ff_results['YFF_1'] = {
                'value' : round(yff_result.get('YFF1', 0), 2)
        }
        gender_results['gd_1'] = {
            'value' : yff_result['gd_1_value'],
            'gender' : yff_result['gd_1_gender']
        }
    else:
        log_and_print(f"calculate_yff error!", 'ERROR')

    # YFF2 (based on Normalization wig)
    wig_norm_file = f"{ANALYSIS_DIR}/{sample_name}/Output_hmmcopy/{sample_name}.proper_paired.50kb.wig.Normalization.txt"
    if os.path.exists(wig_norm_file):
        yff2_result = calculate_yff2(wig_norm_file, config.get("FF_Gender_Config", {}), paths)
        ff_results['YFF_2'] = {
            'value' : round(yff2_result.get('YFF2', 0), 2)
        }
        gender_results['gd_2'] = {
            'value' : yff2_result['gd_2_value'],
            'gender' : yff2_result['gd_2_gender']
        }
        gender_results['gd_3'] = {
            'value' : yff2_result['gd_3_value'],
            'gender' : yff2_result['gd_3_gender']
        }
    else:
        log_and_print(f"Normalization WIG file not found: {wig_norm_file}", 'WARNING')
        yff2_result = {"FF_chrY_adjusted": 0}
        ff_results['YFF_2'] = 0

    # seqFF
    seqff_result = calculate_seqff(sample_name, bam_path)
    ff_results["SeqFF"] = {
        'value' : seqff_result["seqff_value"]
    }

    df_ff = pd.DataFrame.from_dict(ff_results, orient='index')
    df_ff.to_csv(ff_output_txt, sep='\t', index=True, header=True)

    df_gender = pd.DataFrame.from_dict(gender_results, orient='index')
    df_gender.to_csv(gender_output_txt, sep='\t', index=True, header=True)

    #log_and_print(f"Gender inferred from YFF: {gender}")

    return yff_result, yff2_result, ff_results

def run_wcfamily_prediction(sample_name, paths, bam_type, filter_type, gender):
    log_and_print(f"Running CNV prediction for {sample_name} [{bam_type}]")

    # 출력 디렉토리 설정
    wc_out_dir = f"{ANALYSIS_DIR}/{sample_name}/Output_WC/{bam_type}"
    wcx_out_dir = f"{ANALYSIS_DIR}/{sample_name}/Output_WCX/{bam_type}"
    wcff_out_dir = f"{ANALYSIS_DIR}/{sample_name}/Output_WCFF/{bam_type}"

    os.makedirs(wc_out_dir, exist_ok=True)
    os.makedirs(wcx_out_dir, exist_ok=True)
    os.makedirs(wcff_out_dir, exist_ok=True)

    # 완료 파일 경로 정의
    wc_report_txt = os.path.join(wc_out_dir, f"{sample_name}.wc.{bam_type}.report.txt")
    wcx_out_bed = os.path.join(wcx_out_dir, f"{sample_name}.wcx.{bam_type}_bins.bed")
    wcff_output = os.path.join(wcff_out_dir, f"{sample_name}.wcff.{bam_type}.output")  # WisecondorFF 완료 파일

    # 완료 상태 확인
    wc_completed = os.path.exists(wc_report_txt)
    wcx_completed = os.path.exists(wcx_out_bed)

    log_and_print(f"[Prev Completed] WC : {wc_completed}, WCX : {wcx_completed}")
    #wcff_completed = os.path.exists(wcff_output)  # WisecondorFF 완료 파일에 맞게 수정

    # 모든 분석이 완료되었으면 전체 스킵
    #if wc_completed and wcx_completed #and wcff_completed:
    #    log_and_print(f"All CNV predictions already completed for {sample_name} [{bam_type}] - skipping")
    #    return

    # ----------------------------
    # Wisecondor (Python2-based)
    # ----------------------------
    if not wc_completed:
        if filter_type == "proper_paired":
            # proper_paired doesn't have filter_type. It's "orig" as default.
            wc_npz_input = f"{ANALYSIS_DIR}/{sample_name}/Output_WC/{sample_name}.wc.{filter_type}.npz"
        else:
            wc_npz_input = f"{ANALYSIS_DIR}/{sample_name}/Output_WC/{sample_name}.wc.{filter_type}_{bam_type}.npz"

        if not os.path.exists(wc_npz_input):
            log_and_print(f"{wc_npz_input} not found")
            return False

        wc_config = paths["config"].get("WC", {})
        #wc_threshold = wc_config.get(f"{bam_type}_threshold", 6) # WC has not detection z-score
        cyto_file = wc_config.get("cytoband", "/Work/NIPT/data/bed/common/cytoBand.txt")
        ref_wc = paths.get(f"ref_wc_{bam_type}_{filter_type}")

        if not os.path.exists(ref_wc):
            log_and_print(f"{ref_wc} not found")
            return False

        log_and_print(f"{wc_npz_input}")
        log_and_print(f"{cyto_file}")
        log_and_print(f"{ref_wc}")

        try:
            if ref_wc and os.path.exists(ref_wc) and os.path.exists(wc_npz_input):
                log_and_print(f"{wc_npz_input} is being analyzed with {ref_wc}")
                wc_out_dir = os.path.join(ANALYSIS_DIR, sample_name, "Output_WC", bam_type)
                #os.makedirs(wc_out_dir, exist_ok=True)

                wc_out_npz = os.path.join(wc_out_dir, f"{sample_name}.wc.{bam_type}.out.npz")
                # It's defined above
                #wc_report_txt = os.path.join(wc_out_dir, f"{sample_name}.wc.{bam_type}.report.txt")
                wc_plot_png = os.path.join(wc_out_dir, f"{sample_name}.wc.{bam_type}")

                wc_script = os.environ.get("WC", "/opt/wisecondor/wisecondor.py")
                python2 = os.environ.get("PYTHON2", "python2.7")

                run_command("WC test", f"{python2} {wc_script} test {wc_npz_input} {wc_out_npz} {ref_wc}")
                run_command("WC report", f"{python2} {wc_script} report {wc_npz_input} {wc_out_npz} > {wc_report_txt}")
                run_command("WC plot", f"{python2} {wc_script} plot -cytofile {cyto_file} -filetype png {wc_out_npz} {wc_plot_png}")
            else:
                log_and_print(f"Could not find reference or input npz file", 'WARNING')

        except Exception as e:
            log_and_print(f"Could not find reference or input npz file: {e}", 'WARNING')
            return False
    else:
        log_and_print(f"{wc_report_txt} already exist! skip it")

    # ----------------------------
    # WisecondorX
    # ----------------------------
    # Define zscore threshold from config
    if not wcx_completed:
        log_and_print("WisecondorX analysis started.... ")
        log_and_print(f"{bam_type}, {filter_type}")
        if filter_type == "proper_paired":
            # proper_paired doesn't have filter_type. It's "orig" as default.
            wcx_npz_input = f"{ANALYSIS_DIR}/{sample_name}/Output_WCX/{sample_name}.wcx.{filter_type}.npz"
        else:
            wcx_npz_input = f"{ANALYSIS_DIR}/{sample_name}/Output_WCX/{sample_name}.wcx.{filter_type}_{bam_type}.npz"

        if not os.path.exists(wcx_npz_input):
            log_and_print(f"{wcx_npz_input} not found")
            return False

        wcx_config = paths["config"].get("WCX", {})
        wcx_threshold = wcx_config.get(f"{bam_type}_threshold", 6)
        ref_wcx = paths.get(f"ref_wcx_{bam_type}_{filter_type}")

        if not os.path.exists(ref_wcx):
            log_and_print(f"{ref_wcx} not found")
            return False

        try:
            if ref_wcx and os.path.exists(ref_wcx) and os.path.exists(wcx_npz_input):
                log_and_print(f"{wcx_npz_input} is being analyzed with {ref_wcx}")
                wcx_out_dir = os.path.join(ANALYSIS_DIR, sample_name, "Output_WCX", bam_type)
                #os.makedirs(wcx_out_dir, exist_ok=True)

                wcx_out_npz = os.path.join(wcx_out_dir, f"{sample_name}.wcx.{bam_type}")
                # It's defined above
                #wcx_out_bed = os.path.join(wcx_out_dir, f"{sample_name}.wcx.{bam_type}.bed")
                #wcx_out_png = os.path.join(wcx_out_dir, f"{sample_name}.wcx.{bam_type}.png")

                wcx_bin = os.environ.get("WCX", "WisecondorX")
                run_command(
                    f"WCX predict {bam_type}",
                    f"{wcx_bin} predict --zscore {wcx_threshold} {wcx_npz_input} {ref_wcx} {wcx_out_npz} "
                    f"--plot --bed --seed 100"
                )

                # Gender 검증 추가
                wcx_gender_txt = os.path.join(wcx_out_dir, f"{sample_name}.wcx.{bam_type}_gender.txt")
                run_command(
                    f"WCX gender validation {bam_type}",
                    f"{wcx_bin} gender {wcx_npz_input} {ref_wcx} > {wcx_gender_txt}"
                )

        except Exception as e:
            log_and_print(f"Could not find reference or input npz file: {e}", 'WARNING')
            return False

    else:
        log_and_print(f"{wcx_out_bed} already exist! skip it")

    # ----------------------------
    # WisecondorFF
    # ----------------------------
    '''
    wcff_output = f"{ANALYSIS_DIR}/{sample_name}/Output_WCFF/{sample_name}.wcff.{filter_type}.npz"
    ref_wcff = paths.get(f"ref_wcff_{bam_type}")
    if ref_wcff and os.path.exists(ref_wcff) and os.path.exists(wcff_npz_input):
        wcff_out_dir = os.path.join(ANALYSIS_DIR, sample_name, "Output_WCFF", bam_type)
        os.makedirs(wcff_out_dir, exist_ok=True)

        wcff_out_bed = os.path.join(wcff_out_dir, f"{sample_name}.wcff.{bam_type}.bed")

        wcff_bin = os.environ.get("WCFF", "wisecondorff")
        run_command(
            f"WCFF predict {bam_type}",
            f"{wcff_bin} predict {wcff_npz_input} {ref_wcff} {wcff_out_bed}"
        )
    '''
    return True

def run_prizm_pipeline(sample_name, config):
    """Run PRIZM analysis pipeline"""
    log_and_print(f"=== Starting PRIZM Analysis for {sample_name} ===")

    # Check if PRIZM module is available
    if run_prizm_analysis is None:
        log_and_print("PRIZM module not available. Skipping PRIZM analysis.", 'WARNING')
        return False

    # Get PRIZM configuration
    prizm_config = config.get('PRIZM', {})
    qc_cutoff = prizm_config.get('qc_cutoff', 3.0)
    enable_plots = prizm_config.get('enable_plots', True)

    # Get reference files
    lab_references = config.get('lab_references', {})
    mean_file = lab_references.get('prizm_mean')
    sd_file = lab_references.get('prizm_sd')
    mean_10mb_file = lab_references.get('prizm_mean_10mb')
    sd_10mb_file = lab_references.get('prizm_sd_10mb')

    # Check if all reference files exist
    required_files = [mean_file, sd_file, mean_10mb_file, sd_10mb_file]
    missing_files = [f for f in required_files if not f or not os.path.exists(f)]

    if missing_files:
        log_and_print(f"PRIZM reference files missing: {missing_files}", 'WARNING')
        log_and_print("Skipping PRIZM analysis due to missing reference files", 'WARNING')
        return False

    # Create PRIZM output directory
    prizm_output_dir = f"{ANALYSIS_DIR}/{sample_name}/Output_PRIZM"
    os.makedirs(prizm_output_dir, exist_ok=True)

    # Look for 10mb count file (from HMMcopy analysis)
    count_file_10mb = f"{ANALYSIS_DIR}/{sample_name}/Output_hmmcopy/{sample_name}.proper_paired.10mb.wig"

    if not os.path.exists(count_file_10mb):
        log_and_print(f"10mb count file not found: {count_file_10mb}", 'ERROR')
        log_and_print("PRIZM analysis requires HMMcopy 10mb wig file", 'ERROR')
        return False

    try:
        # Convert wig to count format if needed
        count_file_converted = f"{prizm_output_dir}/{sample_name}.10mb.count"
        if not os.path.exists(count_file_converted):
            convert_wig_to_count(count_file_10mb, count_file_converted)

        # Run PRIZM analysis
        log_and_print("Running PRIZM Z-score calculation...")
        results = run_prizm_analysis(
            count_file_10mb=count_file_converted,
            mean_file=mean_file,
            sd_file=sd_file,
            mean_10mb_file=mean_10mb_file,
            sd_10mb_file=sd_10mb_file,
            sample_name=sample_name,
            qc_cutoff=qc_cutoff,
            skip_plots=not enable_plots
        )

        log_and_print("PRIZM analysis completed successfully")

        # Log results summary
        log_and_print(f"Chromosome Z-scores matrix: {results.zscore_chr.shape}")
        log_and_print(f"10mb Z-scores matrix: {results.zscore_10mb.shape}")

        # Create summary report
        create_prizm_summary_report(sample_name, results, prizm_output_dir)

        return True

    except Exception as e:
        log_and_print(f"PRIZM analysis failed: {str(e)}", 'ERROR')
        return False

def copy_outputs_to_final_dir_simple(sample_name: str, analysis_dir: str, output_dir: str):
    """
    Copies relevant output files from analysis_dir to the final output_dir.
    Organizes results under Output_WC, Output_WCX, Output_WCFF, Output_EZD, Output_PRIZM, and Output_QC.
    """
    group_list = ["orig", "fetus", "mom"]

    def safe_copy(src, dst):
        try:
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            if os.path.exists(src):
                shutil.copy2(src, dst)
        except Exception as e:
            log_and_print(f"Failed to copy {src} to {dst}: {e}", 'ERROR')

    def safe_copy_to_dir(src, dst_dir):
        try:
            os.makedirs(dst_dir, exist_ok=True)
            if os.path.exists(src):
                shutil.copy2(src, dst_dir)  # 디렉토리로 복사 (원본 파일명 유지)
        except Exception as e:
            log_and_print(f"Failed to copy {src} to {dst_dir}: {e}", 'ERROR')

    def copy_all_files(src_dir, dst_dir):
        """디렉토리의 모든 파일을 다른 디렉토리로 복사"""
        # 방법 1: glob 사용
        for file_path in glob.glob(os.path.join(src_dir, '*')):
            if os.path.isfile(file_path):
                safe_copy_to_dir(file_path, dst_dir)

    for group in group_list:
        # WC
        wc_src = os.path.join(analysis_dir, sample_name, "Output_WC", group)
        safe_copy_to_dir(os.path.join(wc_src, f"{sample_name}.wc.{group}.report.txt"),
                  os.path.join(output_dir, sample_name, "Output_WC"))

        safe_copy_to_dir(os.path.join(wc_src, f"{sample_name}.wc.{group}.png"),
                  os.path.join(output_dir, sample_name, "Output_WC"))

        # WCX
        wcx_src = os.path.join(analysis_dir, sample_name, "Output_WCX", group)
        safe_copy_to_dir(os.path.join(wcx_src, f"{sample_name}.wcx.{group}_aberrations.bed"),
                  os.path.join(output_dir, sample_name, "Output_WCX"))
        safe_copy(os.path.join(wcx_src, f"{sample_name}.wcx.{group}.plots/genome_wide.png"),
                  os.path.join(output_dir, sample_name, "Output_WCX", f"{sample_name}.wcx.{group}.png"))

        # EZD
        ezd_src = os.path.join(analysis_dir, sample_name, "Output_EZD", group)
        safe_copy_to_dir(os.path.join(ezd_src, f"{group}_EZD_grid.png"), os.path.join(output_dir, sample_name, "Output_EZD"))

        # PRIZM
        prizm_src = os.path.join(analysis_dir, sample_name, "Output_PRIZM", group)
        safe_copy_to_dir(os.path.join(prizm_src, f"{sample_name}_{group}_chromosome_line.png"),
                  os.path.join(output_dir, sample_name, "Output_PRIZM"))
        safe_copy_to_dir(os.path.join(prizm_src, f"{sample_name}_{group}_10mb_line.png"),
                  os.path.join(output_dir, sample_name, "Output_PRIZM"))

    # QC
    qc_src = os.path.join(analysis_dir, sample_name, "Output_QC")

    copy_all_files(os.path.join(qc_src), os.path.join(output_dir, sample_name, "Output_QC"))

    log_and_print(f"Copied final outputs for sample {sample_name} to output directory.")
    return True

# Progress Tracking functions
def run_pipeline_step(step_num, step_name, step_function, progress, *args):
    """공통 단계 실행 패턴"""
    try:
        result = step_function(*args)
        if result:
            progress.update_step(step_num, step_name, "PASS")
            return True
        else:
            progress.update_step(step_num, step_name, "FAIL")
            return False
    except Exception as e:
        progress.update_step(step_num, step_name, "FAIL")
        log_and_print(f"Step {step_num} failed: {e}", 'ERROR')
        return False

def return_with_error(progress, code):
    progress.mark_failed("Pipeline step failed")
    return code

def run_microdeletion_step(sample_name, labcode, config, analysis_dir, output_dir, bed_dir):
    """외부 모듈을 래핑하는 함수"""
    try:
        md_success = run_microdeletion_decision_pipeline(
            sample_name=sample_name,
            labcode=labcode,
            config=config,
            analysis_dir=analysis_dir,
            output_dir=output_dir,
            bed_dir=bed_dir
        )
        return md_success  # True/False 반환

    except Exception as e:
        log_and_print(f"Microdeletion pipeline failed: {e}", 'ERROR')
        return False

def setup_unified_logging(sample_name):
    """통합 로깅 설정 - 외부 모듈까지 포함"""
    import sys
    import logging
    import os

    # Docker 실시간 출력 설정
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)

    # 로그 디렉토리 생성
    os.makedirs(f"{LOG_DIR}/{sample_name}", exist_ok=True)
    os.makedirs(f"{ANALYSIS_DIR}/{sample_name}", exist_ok=True)

    # 로그 파일 경로
    main_log = f"{LOG_DIR}/{sample_name}/{sample_name}.pipeline.log"
    analysis_log = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.analysis.log"

    # Root logger 설정 (모든 모듈에 적용됨)
    root_logger = logging.getLogger()

    # 기존 핸들러 제거
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # 통합 포매터
    formatter = logging.Formatter(
        '[%(asctime)s] [%(levelname)s] %(name)s:%(lineno)d %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # 콘솔 핸들러 (Docker logs)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    console_handler.setLevel(logging.INFO)

    # 파일 핸들러들
    main_file_handler = logging.FileHandler(main_log, mode='a')
    main_file_handler.setFormatter(formatter)
    main_file_handler.setLevel(logging.INFO)

    analysis_file_handler = logging.FileHandler(analysis_log, mode='a')
    analysis_file_handler.setFormatter(formatter)
    analysis_file_handler.setLevel(logging.INFO)

    # Root logger에 핸들러 추가
    root_logger.addHandler(console_handler)
    root_logger.addHandler(main_file_handler)
    root_logger.addHandler(analysis_file_handler)
    root_logger.setLevel(logging.INFO)

    # 초기 메시지
    logger = logging.getLogger(__name__)
    logger.info(f"=== Unified logging initialized for {sample_name} ===")
    logger.info(f"Main log: {main_log}")
    logger.info(f"Analysis log: {analysis_log}")

    return logger

def setup_logging(sample_name):
    """통합 로깅 설정"""
    import sys
    import logging
    import os
    
    # Docker 실시간 출력을 위한 설정
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)
    
    # 로그 디렉토리 생성
    os.makedirs(f"{LOG_DIR}/{sample_name}", exist_ok=True)
    os.makedirs(f"{ANALYSIS_DIR}/{sample_name}", exist_ok=True)
    
    # 로그 파일 경로
    pipeline_log = f"{LOG_DIR}/{sample_name}/{sample_name}.pipeline.log"
    analysis_log = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.analysis.log"
    
    # 기존 핸들러 제거
    logger = logging.getLogger()
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # 통합 로깅 설정
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] [%(levelname)s] %(filename)s:%(lineno)d %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[
            logging.StreamHandler(sys.stdout),  # Docker logs에 표시
            logging.FileHandler(pipeline_log, mode='a'),  # 파이프라인 로그
            logging.FileHandler(analysis_log, mode='a')   # 분석 로그 (기존 호환성)
        ]
    )
    
    logger = logging.getLogger(__name__)
    log_and_print(f"=== Logging initialized for {sample_name} ===")
    log_and_print(f"Pipeline log: {pipeline_log}")
    log_and_print(f"Analysis log: {analysis_log}")
    
    return logger

def log_and_print(message, level='INFO'):
    """로그 파일과 콘솔에 동시 출력"""
    print(f"[{level}] {message}", flush=True)  # flush=True로 즉시 출력
    logger = logging.getLogger(__name__)

    if level == 'INFO':
        logger.info(message)
    elif level == 'ERROR':
        logger.error(message)
    elif level == 'WARNING':
        logger.warning(message)
    elif level == 'DEBUG':
        logger.debug(message)

def get_file_hash(filepath):
    """Calculate hash of a single file"""
    if not os.path.exists(filepath):
        return None

    hash_md5 = hashlib.md5()
    try:
        # Add file path and modification time
        hash_md5.update(filepath.encode())
        mtime = os.path.getmtime(filepath)
        hash_md5.update(str(mtime).encode())

        # Add file size
        file_size = os.path.getsize(filepath)
        hash_md5.update(str(file_size).encode())

        return hash_md5.hexdigest()
    except (IOError, OSError):
        return None

def get_directory_hash(directory_path):
    """Calculate hash of all files in directory"""
    if not os.path.exists(directory_path):
        return None

    hash_md5 = hashlib.md5()

    for root, dirs, files in os.walk(directory_path):
        dirs.sort()
        files.sort()

        for filename in files:
            filepath = os.path.join(root, filename)
            file_hash = get_file_hash(filepath)
            if file_hash:
                hash_md5.update(file_hash.encode())

    return hash_md5.hexdigest()

def get_source_files_hash(analysis_dir, sample_name, group_list):
    """Calculate hash of all source files that will be copied"""
    hash_md5 = hashlib.md5()

    # Hash individual files for each group
    for group in group_list:
        files_to_check = [
            # WC files
            os.path.join(analysis_dir, sample_name, "Output_WC", group, f"{sample_name}.wc.{group}.report.txt"),
            os.path.join(analysis_dir, sample_name, "Output_WC", group, f"{sample_name}.wc.{group}.png"),
            # WCX files
            os.path.join(analysis_dir, sample_name, "Output_WCX", group, f"{sample_name}.wcx.{group}_aberrations.bed"),
            os.path.join(analysis_dir, sample_name, "Output_WCX", group, f"{sample_name}.wcx.{group}.plots", "genome_wide.png"),
            # EZD files
            os.path.join(analysis_dir, sample_name, "Output_EZD", group, f"{group}_EZD_grid.png"),
            # PRIZM files
            os.path.join(analysis_dir, sample_name, "Output_PRIZM", group, f"{sample_name}_{group}_chromosome_line.png"),
            os.path.join(analysis_dir, sample_name, "Output_PRIZM", group, f"{sample_name}_{group}_10mb_line.png"),
        ]

        for filepath in files_to_check:
            file_hash = get_file_hash(filepath)
            if file_hash:
                hash_md5.update(file_hash.encode())

    # Hash QC directory (entire directory)
    qc_dir = os.path.join(analysis_dir, sample_name, "Output_QC")
    qc_hash = get_directory_hash(qc_dir)
    if qc_hash:
        hash_md5.update(qc_hash.encode())

    return hash_md5.hexdigest()

def create_copy_marker(marker_path, source_hash, copy_info):
    """Create marker file with copy information"""
    marker_data = {
        "copy_completed": True,
        "source_hash": source_hash,
        "copy_timestamp": str(pd.Timestamp.now()),
        "copy_info": copy_info
    }

    try:
        with open(marker_path, 'w') as f:
            json.dump(marker_data, f, indent=2)
        log_and_print(f"Copy marker created: {marker_path}")
    except Exception as e:
        log_and_print(f"Failed to create copy marker {marker_path}: {e}", "ERROR")

def check_copy_marker(marker_path, current_source_hash):
    """Check if copy marker exists and is valid"""
    if not os.path.exists(marker_path):
        return False

    try:
        with open(marker_path, 'r') as f:
            marker_data = json.load(f)

        # Check if copy was completed
        if not marker_data.get("copy_completed", False):
            return False

        # Check if source hash matches (detect changes)
        stored_hash = marker_data.get("source_hash", "")
        if stored_hash != current_source_hash:
            log_and_print("Source files have changed, copy needed")
            return False

        log_and_print("Copy marker valid, skipping copy operation", "ERROR")
        return True

    except Exception as e:
        log_and_print(f"Error reading copy marker {marker_path}: {e}")
        return False

def safe_copy(src, dst):
    """Safely copy a file with error handling"""
    try:
        if os.path.exists(src):
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy2(src, dst)
            return True
        else:
            log_and_print(f"Source file not found: {src}", "WARNING")
            return False
    except Exception as e:
        log_and_print(f"Failed to copy {src} to {dst}: {e}", "ERROR")
        return False

def safe_copy_to_dir(src, dst_dir):
    """Safely copy a file to directory with error handling"""
    try:
        if os.path.exists(src):
            os.makedirs(dst_dir, exist_ok=True)
            shutil.copy2(src, dst_dir)
            return True
        else:
            log_and_print(f"Source file not found: {src}", "WARNING")
            return False
    except Exception as e:
        log_and_print(f"Failed to copy {src} to {dst_dir}: {e}", "ERROR")
        return False

def copy_all_files(src_dir, dst_dir):
    """Copy all files from source directory to destination directory"""
    try:
        if os.path.exists(src_dir):
            os.makedirs(dst_dir, exist_ok=True)
            shutil.copytree(src_dir, dst_dir, dirs_exist_ok=True)
            return True
        else:
            log_and_print(f"Source directory not found: {src_dir}", "WARNING")
            return False
    except Exception as e:
        log_and_print(f"Failed to copy directory {src_dir} to {dst_dir}: {e}", "ERROR")
        return False

def copy_outputs_to_final_dir(sample_name, analysis_dir, output_dir, group_list=None, force_copy=False):
    """
    Copy specific output files to final directory with marker-based optimization

    Args:
        sample_name: Sample identifier
        analysis_dir: Source analysis directory
        output_dir: Destination output directory
        group_list: List of groups to process (default: ['orig', 'fetus', 'mom'])
        force_copy: Force copy even if marker exists
    """

    if group_list is None:
        group_list = ['orig', 'fetus', 'mom']

    dest_dir = os.path.join(output_dir, sample_name)
    marker_path = os.path.join(dest_dir, "copied.marker")

    # Calculate current source hash for specific files only
    current_hash = get_source_files_hash(analysis_dir, sample_name, group_list)
    if current_hash is None:
        logger.error(f"Failed to calculate hash for source files: {sample_name}")
        return False

    # Check marker if not forcing copy
    if not force_copy and check_copy_marker(marker_path, current_hash):
        logger.info(f"Copy skipped for {sample_name} - already up to date")
        return True

    logger.info(f"Starting selective copy operation for {sample_name}")

    try:
        # Create destination directory
        os.makedirs(dest_dir, exist_ok=True)

        # Remove existing marker during copy
        if os.path.exists(marker_path):
            os.remove(marker_path)

        # Track copy statistics
        copy_stats = {
            "files_copied": 0,
            "files_failed": 0,
            "total_size": 0
        }

        # Copy files for each group
        for group in group_list:
            logger.info(f"Processing group: {group}")

            # WC files
            wc_src = os.path.join(analysis_dir, sample_name, "Output_WC", group)
            wc_files = [
                (os.path.join(wc_src, f"{sample_name}.wc.{group}.report.txt"),
                 os.path.join(output_dir, sample_name, "Output_WC")),
                (os.path.join(wc_src, f"{sample_name}.wc.{group}_z.png"),
                 os.path.join(output_dir, sample_name, "Output_WC"))
            ]

            for src, dst_dir in wc_files:
                if safe_copy_to_dir(src, dst_dir):
                    copy_stats["files_copied"] += 1
                    if os.path.exists(src):
                        copy_stats["total_size"] += os.path.getsize(src)
                else:
                    copy_stats["files_failed"] += 1

            # WCX files
            wcx_src = os.path.join(analysis_dir, sample_name, "Output_WCX", group)

            # Copy aberrations.bed file
            if safe_copy_to_dir(os.path.join(wcx_src, f"{sample_name}.wcx.{group}_aberrations.bed"),
                              os.path.join(output_dir, sample_name, "Output_WCX")):
                copy_stats["files_copied"] += 1
                src_file = os.path.join(wcx_src, f"{sample_name}.wcx.{group}_aberrations.bed")
                if os.path.exists(src_file):
                    copy_stats["total_size"] += os.path.getsize(src_file)
            else:
                copy_stats["files_failed"] += 1

            # Copy genome_wide.png with rename
            if safe_copy(os.path.join(wcx_src, f"{sample_name}.wcx.{group}.plots", "genome_wide.png"),
                        os.path.join(output_dir, sample_name, "Output_WCX", f"{sample_name}.wcx.{group}.png")):
                copy_stats["files_copied"] += 1
                src_file = os.path.join(wcx_src, f"{sample_name}.wcx.{group}.plots", "genome_wide.png")
                if os.path.exists(src_file):
                    copy_stats["total_size"] += os.path.getsize(src_file)
            else:
                copy_stats["files_failed"] += 1

            # Copy WCX chromosome plots (chr*.png files) for each group
            logger.info(f"Processing WCX chromosome plots for {group}")
            wcx_plots_src = os.path.join(analysis_dir, sample_name, "Output_WCX", group, f"{sample_name}.wcx.{group}.plots")
            wcx_plots_dst = os.path.join(output_dir, sample_name, "Output_WCX", "chr_plots", group)

            if os.path.exists(wcx_plots_src):
                os.makedirs(wcx_plots_dst, exist_ok=True)

                # Find all chr*.png files
                import glob
                chr_plot_files = glob.glob(os.path.join(wcx_plots_src, "chr*.png"))

                if chr_plot_files:
                    logger.info(f"Found {len(chr_plot_files)} chromosome plot files for {group}")
                    for chr_plot_file in chr_plot_files:
                        filename = os.path.basename(chr_plot_file)
                        dst_file = os.path.join(wcx_plots_dst, filename)

                        try:
                            shutil.copy2(chr_plot_file, dst_file)
                            copy_stats["files_copied"] += 1
                            copy_stats["total_size"] += os.path.getsize(chr_plot_file)
                            logger.debug(f"Copied {group} chromosome plot: {filename}")
                        except Exception as e:
                            logger.error(f"Failed to copy {chr_plot_file}: {e}")
                            copy_stats["files_failed"] += 1
                else:
                    logger.warning(f"No chr*.png files found in {wcx_plots_src}")
            else:
                logger.warning(f"WCX plots directory not found: {wcx_plots_src}")
                copy_stats["files_failed"] += 1

            # EZD files
            ezd_src = os.path.join(analysis_dir, sample_name, "Output_EZD", group)
            if safe_copy_to_dir(os.path.join(ezd_src, f"{group}_EZD_grid.png"),
                              os.path.join(output_dir, sample_name, "Output_EZD")):
                copy_stats["files_copied"] += 1
                src_file = os.path.join(ezd_src, f"{group}_EZD_grid.png")
                if os.path.exists(src_file):
                    copy_stats["total_size"] += os.path.getsize(src_file)
            else:
                copy_stats["files_failed"] += 1

            # PRIZM files
            prizm_src = os.path.join(analysis_dir, sample_name, "Output_PRIZM", group)
            prizm_files = [
                (os.path.join(prizm_src, f"{sample_name}_{group}_chromosome_line.png"),
                 os.path.join(output_dir, sample_name, "Output_PRIZM")),
                (os.path.join(prizm_src, f"{sample_name}_{group}_10mb_line.png"),
                 os.path.join(output_dir, sample_name, "Output_PRIZM"))
            ]

            for src, dst_dir in prizm_files:
                if safe_copy_to_dir(src, dst_dir):
                    copy_stats["files_copied"] += 1
                    if os.path.exists(src):
                        copy_stats["total_size"] += os.path.getsize(src)
                else:
                    copy_stats["files_failed"] += 1

        # Copy entire QC directory
        logger.info("Processing QC directory")
        qc_src = os.path.join(analysis_dir, sample_name, "Output_QC")
        if copy_all_files(qc_src, os.path.join(output_dir, sample_name, "Output_QC")):
            # Count QC files
            qc_dst = os.path.join(output_dir, sample_name, "Output_QC")
            if os.path.exists(qc_dst):
                for root, dirs, files in os.walk(qc_dst):
                    copy_stats["files_copied"] += len(files)
                    for file in files:
                        filepath = os.path.join(root, file)
                        if os.path.exists(filepath):
                            copy_stats["total_size"] += os.path.getsize(filepath)
        else:
            copy_stats["files_failed"] += 1

        logger.info(f"Copy completed: {copy_stats['files_copied']} files successful, "
                   f"{copy_stats['files_failed']} files failed, "
                   f"{copy_stats['total_size'] / (1024*1024):.2f} MB total")

        # Create marker file after successful copy
        create_copy_marker(marker_path, current_hash, copy_stats)

        return copy_stats["files_failed"] == 0  # Return True if no failures

    except Exception as e:
        logger.error(f"Copy operation failed for {sample_name}: {e}")

        # Remove incomplete marker if exists
        if os.path.exists(marker_path):
            try:
                os.remove(marker_path)
            except:
                pass

        return False

def copy_outputs_to_final_dir_old(sample_name, analysis_dir, output_dir, group_list=None, force_copy=False):
    """
    Copy specific output files to final directory with marker-based optimization

    Args:
        sample_name: Sample identifier
        analysis_dir: Source analysis directory
        output_dir: Destination output directory
        group_list: List of groups to process (default: ['orig', 'fetus', 'mom'])
        force_copy: Force copy even if marker exists
    """

    if group_list is None:
        group_list = ['orig', 'fetus', 'mom']

    dest_dir = os.path.join(output_dir, sample_name)
    marker_path = os.path.join(dest_dir, "copied.marker")

    # Calculate current source hash for specific files only
    current_hash = get_source_files_hash(analysis_dir, sample_name, group_list)
    if current_hash is None:
        log_and_print(f"Failed to calculate hash for source files: {sample_name}", "ERROR")
        return False

    # Check marker if not forcing copy
    if not force_copy and check_copy_marker(marker_path, current_hash):
        log_and_print(f"Copy skipped for {sample_name} - already up to date")
        return True

    log_and_print(f"Starting selective copy operation for {sample_name}")

    try:
        # Create destination directory
        os.makedirs(dest_dir, exist_ok=True)

        # Remove existing marker during copy
        if os.path.exists(marker_path):
            os.remove(marker_path)

        # Track copy statistics
        copy_stats = {
            "files_copied": 0,
            "files_failed": 0,
            "total_size": 0
        }

        # Copy files for each group
        for group in group_list:
            log_and_print(f"Processing group: {group}")

            # WC files
            wc_src = os.path.join(analysis_dir, sample_name, "Output_WC", group)
            wc_files = [
                (os.path.join(wc_src, f"{sample_name}.wc.{group}.report.txt"),
                 os.path.join(output_dir, sample_name, "Output_WC")),
                (os.path.join(wc_src, f"{sample_name}.wc.{group}.png"),
                 os.path.join(output_dir, sample_name, "Output_WC"))
            ]

            for src, dst_dir in wc_files:
                if safe_copy_to_dir(src, dst_dir):
                    copy_stats["files_copied"] += 1
                    if os.path.exists(src):
                        copy_stats["total_size"] += os.path.getsize(src)
                else:
                    copy_stats["files_failed"] += 1

            # WCX files
            wcx_src = os.path.join(analysis_dir, sample_name, "Output_WCX", group)

            # Copy aberrations.bed file
            if safe_copy_to_dir(os.path.join(wcx_src, f"{sample_name}.wcx.{group}_aberrations.bed"),
                              os.path.join(output_dir, sample_name, "Output_WCX")):
                copy_stats["files_copied"] += 1
                src_file = os.path.join(wcx_src, f"{sample_name}.wcx.{group}_aberrations.bed")
                if os.path.exists(src_file):
                    copy_stats["total_size"] += os.path.getsize(src_file)
            else:
                copy_stats["files_failed"] += 1

            # Copy genome_wide.png with rename
            if safe_copy(os.path.join(wcx_src, f"{sample_name}.wcx.{group}.plots", "genome_wide.png"),
                        os.path.join(output_dir, sample_name, "Output_WCX", f"{sample_name}.wcx.{group}.png")):
                copy_stats["files_copied"] += 1
                src_file = os.path.join(wcx_src, f"{sample_name}.wcx.{group}.plots", "genome_wide.png")
                if os.path.exists(src_file):
                    copy_stats["total_size"] += os.path.getsize(src_file)
            else:
                copy_stats["files_failed"] += 1

            # EZD files
            ezd_src = os.path.join(analysis_dir, sample_name, "Output_EZD", group)
            if safe_copy_to_dir(os.path.join(ezd_src, f"{group}_EZD_grid.png"),
                              os.path.join(output_dir, sample_name, "Output_EZD")):
                copy_stats["files_copied"] += 1
                src_file = os.path.join(ezd_src, f"{group}_EZD_grid.png")
                if os.path.exists(src_file):
                    copy_stats["total_size"] += os.path.getsize(src_file)
            else:
                copy_stats["files_failed"] += 1

            # PRIZM files
            prizm_src = os.path.join(analysis_dir, sample_name, "Output_PRIZM", group)
            prizm_files = [
                (os.path.join(prizm_src, f"{sample_name}_{group}_chromosome_line.png"),
                 os.path.join(output_dir, sample_name, "Output_PRIZM")),
                (os.path.join(prizm_src, f"{sample_name}_{group}_10mb_line.png"),
                 os.path.join(output_dir, sample_name, "Output_PRIZM"))
            ]

            for src, dst_dir in prizm_files:
                if safe_copy_to_dir(src, dst_dir):
                    copy_stats["files_copied"] += 1
                    if os.path.exists(src):
                        copy_stats["total_size"] += os.path.getsize(src)
                else:
                    copy_stats["files_failed"] += 1

        # Copy entire QC directory
        log_and_print("Processing QC directory")
        qc_src = os.path.join(analysis_dir, sample_name, "Output_QC")
        if copy_all_files(qc_src, os.path.join(output_dir, sample_name, "Output_QC")):
            # Count QC files
            qc_dst = os.path.join(output_dir, sample_name, "Output_QC")
            if os.path.exists(qc_dst):
                for root, dirs, files in os.walk(qc_dst):
                    copy_stats["files_copied"] += len(files)
                    for file in files:
                        filepath = os.path.join(root, file)
                        if os.path.exists(filepath):
                            copy_stats["total_size"] += os.path.getsize(filepath)
        else:
            copy_stats["files_failed"] += 1

        log_and_print(f"Copy completed: {copy_stats['files_copied']} files successful, "
                   f"{copy_stats['files_failed']} files failed, "
                   f"{copy_stats['total_size'] / (1024*1024):.2f} MB total")

        # Create marker file after successful copy
        create_copy_marker(marker_path, current_hash, copy_stats)

        return copy_stats["files_failed"] == 0  # Return True if no failures

    except Exception as e:
        log_and_print(f"Copy operation failed for {sample_name}: {e}", "ERROR")

        # Remove incomplete marker if exists
        if os.path.exists(marker_path):
            try:
                os.remove(marker_path)
            except:
                pass

        return False

def cleanup_old_copies(output_dir, max_age_days=30):
    """Clean up old copied outputs based on marker timestamps"""
    import pandas as pd
    from datetime import timedelta

    if not os.path.exists(output_dir):
        return

    cutoff_time = pd.Timestamp.now() - timedelta(days=max_age_days)
    cleaned_count = 0

    for item in os.listdir(output_dir):
        item_path = os.path.join(output_dir, item)
        marker_path = os.path.join(item_path, "copied.marker")

        if not os.path.isdir(item_path) or not os.path.exists(marker_path):
            continue

        try:
            with open(marker_path, 'r') as f:
                marker_data = json.load(f)

            copy_time = pd.Timestamp(marker_data.get("copy_timestamp", "1970-01-01"))

            if copy_time < cutoff_time:
                log_and_print(f"Cleaning up old copy: {item_path}")
                shutil.rmtree(item_path)
                cleaned_count += 1

        except Exception as e:
            logger.warning(f"Error checking marker for cleanup {marker_path}: {e}")
            continue

    if cleaned_count > 0:
        log_and_print(f"Cleaned up {cleaned_count} old copy directories")


def main():
    args = parse_args()

    sample_name = args.sample_name
    fastq_r1 = args.fastq_r1
    fastq_r2 = args.fastq_r2
    #gender = args.gender
    labcode = args.labcode
    age = args.age

    # Base directory structure
    root_dir = Path("/Work/NIPT")
    ref_common_dir = root_dir / "data" / "refs" / "common"
    config_dir = root_dir / "config" / labcode
    data_dir = root_dir / "data"
    fastq_dir = root_dir / "fastq"
    ref_common_dir = data_dir / "refs" / "common"
    ref_lab_dir = data_dir / "refs" / labcode
    bed_dir = data_dir / "bed"
    output_dir = root_dir / "output"

    #logger = setup_logging(sample_name)
    # To make external scripts use the same logging
    logger = setup_unified_logging(sample_name)

    progress = ProgressTracker(sample_name, output_dir)

    # Define filter paths
    filter_paths = {
        "of": bed_dir / "common" / "Uniform_2017_allY.bed",
        "nf08": bed_dir / labcode / "hg19_mappability_0.8_clean_all_36mer.bed"
        #"nf09": bed_dir / labcode / "hg19_mappability_0.9_clean_all_36mer.bed"
    }

    # Load config
    config_file = config_dir / "pipeline_config.json"
    with open(config_file) as cf:
        config = json.load(cf)

    fq1 = os.path.join(fastq_dir, sample_name, fastq_r1)
    fq2 = os.path.join(fastq_dir, sample_name, fastq_r2)
    if not os.path.isfile(fq1) or not os.path.isfile(fq2):
        log_and_print("FASTQ files not found after symbolic linking.", 'ERROR')
        return False

    log_and_print(f"Starting NIPT pipeline for sample: {sample_name}")
    log_and_print(f"Lab code: {labcode}, Age: {age}")
    log_and_print(f"FASTQ files: {fastq_r1}, {fastq_r2}")

    final_marker = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.pipeline_completed.marker"
    # 250529 : Block it temporarily
    if check_file_exists(final_marker, "Pipeline completion marker"):
        log_and_print("Pipeline appears to have been completed previously.")
        return

    log_and_print(f"=== {sample_name}.pipeline_completed.marker was initialized ===")

    log_and_print(f"=== Downsample Fastq check ===")

    # 250529 : Block it temporarily
    if not run_pipeline_step(1, "Downsample", downsample_fastq, progress, sample_name, fastq_r1, fastq_r2, config):
        progress.mark_failed("Downsample failed")
        return 1

    log_and_print(f"=== Running FastQC ===")
    if not run_pipeline_step(2, "FastQC", run_fastqc, progress, sample_name, fastq_r1, fastq_r2):
        progress.mark_failed("FastQC failed")
        return 1

    # Construct paths
    paths = {
        "root": root_dir,
        "config": config,
        "ref_lab": ref_lab_dir,
        "ref_fasta": ref_common_dir / "hg19" / "ucsc.hg19.fasta",
        "wig_gc_50kb": ref_common_dir / "hmmcopy" / "hg19.gc.50kb.wig",
        "wig_map_50kb": ref_common_dir / "hmmcopy" / "hg19.map.50kb.wig",
        "wig_gc_10mb": ref_common_dir / "hmmcopy" / "hg19.gc.10mb.wig",
        "wig_map_10mb": ref_common_dir / "hmmcopy" / "hg19.map.10mb.wig",
        "bed_dir": bed_dir
    }

    #for bam_type in ["orig", "fetus", "mom"]:
    #    paths[f"ref_wc_{bam_type}"] = ref_lab_dir / "WC" / f"{bam_type}_{config['WC']['ref']}.npz"
        # After gender confirmation, allocate belows
        #paths[f"ref_wcx_{bam_type}"] = ref_lab_dir / "WCX" / f"{bam_type}_{gender.upper()}_{config['WCX']['ref']}.npz"
        #paths[f"ref_lomaz_{bam_type}"] = ref_lab_dir / "LoMAz" / bam_type

    # Add MD targets if defined
    paths["md_beds"] = {}
    for md_key in ["MD_Target_8", "MD_Target_87", "MD_Target_320", "MD_Target_116"]:
        if md_key in config:
            paths["md_beds"][md_key] = bed_dir / config[md_key]["bed"]

    log_and_print(f"=== Creating the default directories ===")
    if not run_pipeline_step(3, "Creating directories", create_directories, progress, sample_name):
        progress.mark_failed("Creating directories failed")
        return 1

    # KWON : I need to check if this log files is really needed.
    sample_log_file = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.log"
    with open(sample_log_file, "w") as log_file:
        log_file.write(f"NIPT Pipeline Log for {sample_name}")
        log_file.write(f"Started at: {time.strftime('%Y-%m-%d %H:%M:%S')}")
        log_file.write("="*50 + " ")

    try:
        file_handler = logging.FileHandler(f"{LOG_DIR}/{sample_name}/{sample_name}.pipeline.log")
        file_handler.setFormatter(logging.Formatter('[%(asctime)s] [%(levelname)s] %(message)s'))
        logger.addHandler(file_handler)
        log_and_print("File logging initialized successfully")
    except Exception as e:
        log_and_print(f"Could not initialize file logging: {e}", 'WARNING')

    if not run_pipeline_step(4, "Creating symbolic links of Fastq", create_symbolic_links, progress, sample_name, fastq_r1, fastq_r2):
        progress.mark_failed("Creating directories failed")
        return 1

    log_and_print("=== Starting BAM file generation ===")
    if not generate_proper_paired_bam(sample_name, fastq_r1, fastq_r2, config, progress, base_step = 5):
        progress.mark_failed("BAM generation failed")
        return 1

    log_and_print("=== Starting filter processing ===")
    base_step = 6
    step_counter =0 
    for idx, (filter_type, filter_path) in enumerate(filter_paths.items(), 1):
        step_num = f"{base_step}.{step_counter + 1}"
        log_and_print(f"Processing filter: {filter_type}")

        try:
            process_filter(sample_name, filter_type, filter_path)
            progress.update_step(step_num, f"Process {filter_type} filter", "PASS")
            step_counter += 1

        except Exception as e:
            progress.update_step(step_num, f"Process {filter_type} filter", "FAIL")
            log_and_print(f"Failed to process {filter_type} filter: {e}", 'ERROR')
            return 1

    log_and_print("=== Starting NPZ and HMMcopy file generation ===")
    if not run_pipeline_step(7, "Process HMMcopy", process_hmmcopy, progress, sample_name, paths):
        progress.mark_failed("Process HMMcopy failed")
        return 1

    log_and_print("=== Starting Calculation FF ===")
    try:
        yff_result, yff2_result, ff_results = calculate_fetal_fraction(sample_name, config, paths)
        progress.update_step(8, "Calculate Fetal Fraction", "PASS")

    except Exception as e:
        progress.update_step(8, "Calculate Fetal Fraction", "FAIL")
        log_and_print(f"Failed to calculate fetal fraction: {e}", 'ERROR')
        return 1

    gender_dict = {"XX" : 'F', "XY" : 'M'}
    gender = gender_dict[yff2_result['gd_2_gender']]

    # Gender was determined
    logger.info(f"final gender : {gender}")
    '''
    paths[f"ref_wc_orig_proper_paired"] = ref_lab_dir / "WC" / f"{bam_type}_{config['WC']['ref']}_proper_paired.npz"
    paths[f"ref_wcx_orig_proper_paired"] = ref_lab_dir / "WCX" / f"orig_{confirmed_gender}_{config['WCX']['ref']}_proper_paired.npz"

    for bam_type in ["fetus", "mom"]:
        paths[f"ref_wc_{bam_type}_of"] = ref_lab_dir / "WC" / f"{bam_type}_{config['WC']['ref']}_of.npz"
        paths[f"ref_wcx_{bam_type}_of"] = ref_lab_dir / "WCX" / f"{bam_type}_{confirmed_gender}_{config['WCX']['ref']}_of.npz"
    '''

    # WC/WCX/ WCFF
    '''
    When all references are ready, I'll use the following code.
    # Orig, Fetus ---> Reference is separated by gender
    for bam_type in ["orig", "fetus", "mom"]:
        for filter_type in ['of', 'nf08', 'nf09']:
            paths[f"ref_wc_{bam_type}_{filter_type}"] = paths["ref_lab"] / "WC" / f"{bam_type}_{config['WC']['ref']}_{filter_type}.npz"
            if filter_type == "mom":
                paths[f"ref_wcx_{bam_type}_{filter_type}"] = paths["ref_lab"] / "WCX" / f"{bam_type}_{config['WCX']['ref']}_{filter_type}.npz"
            else:
                paths[f"ref_wcx_{bam_type}_{filter_type}"] = paths["ref_lab"] / "WCX" / f"{bam_type}_{gender}_{config['WCX']['ref']}_{filter_type}.npz"

            paths[f"ref_wcff_{bam_type}_{filter_type}"] = paths["ref_lab"] / "WCFF" / f"{bam_type}_{config['WCFF']['ref']}_{filter_type}.npz"

            run_wcfamily_prediction(sample_name, paths, bam_type, filter_type, gender)
    '''

    md_wcx_orig_result = os.path.join(ANALYSIS_DIR, f"{sample_name}", "Output_WCX", "orig", f"{sample_name}.wcx.orig_bins.bed")
    md_wcx_fetus_result = os.path.join(ANALYSIS_DIR, f"{sample_name}", "Output_WCX", "fetus", f"{sample_name}.wcx.fetus_bins.bed")
    md_wcx_mom_result = os.path.join(ANALYSIS_DIR, f"{sample_name}", "Output_WCX", "mom", f"{sample_name}.wcx.mom_bins.bed")

    # Comment : WC, WCX in run_wcfamily_prediction is performed. If WCFF is added, it would be better to check FF result
    # Orig
    if check_file_exists_advanced(md_wcx_orig_result, "WC/WCX orig"):
        progress.update_step(f"{base_step}.9", "Run WC/WCX/WCFF prediction (orig)", "SKIP", "file exists")
    else:
        filter_type = "proper_paired"
        paths[f"ref_wc_orig_{filter_type}"] = paths["ref_lab"] / "WC" / f"orig_{config['WC']['ref']}_{filter_type}.npz"
        paths[f"ref_wcx_orig_{filter_type}"] = paths["ref_lab"] / "WCX" / f"orig_{gender}_{config['WCX']['ref']}_{filter_type}.npz"

        if not run_pipeline_step(9, "Run WC/WCX/WCFF prediction (orig)", run_wcfamily_prediction, progress, sample_name, paths, "orig", filter_type, gender):
            progress.mark_failed("Run WC/WCX/WCFF prediction (orig) failed")
            return 1

    # Fetus
    if check_file_exists_advanced(md_wcx_fetus_result, "WC/WCX fetus"):
        progress.update_step(f"{base_step}.9", "Run WC/WCX/WCFF prediction (fetus)", "SKIP", "file exists")
    else:
        filter_type = "of"
        paths[f"ref_wc_fetus_{filter_type}"] = paths["ref_lab"] / "WC" / f"fetus_{config['WC']['ref']}_{filter_type}.npz"
        paths[f"ref_wcx_fetus_{filter_type}"] = paths["ref_lab"] / "WCX" / f"fetus_{gender}_{config['WCX']['ref']}_{filter_type}.npz"
        if not run_pipeline_step(10, "Run WC/WCX/WCFF prediction (fetus)", run_wcfamily_prediction, progress, sample_name, paths, "fetus", filter_type, gender):
            progress.mark_failed("Run WC/WCX/WCFF prediction (fetus) failed")
            return 1

    # Mom
    if check_file_exists_advanced(md_wcx_mom_result, "WC/WCX mom"):
        progress.update_step(f"{base_step}.9", "Run WC/WCX/WCFF prediction (mom)", "SKIP", "file exists")
    else:
        filter_type = "of"
        paths[f"ref_wc_mom_{filter_type}"] = paths["ref_lab"] / "WC" / f"mom_{config['WC']['ref']}_{filter_type}.npz"
        paths[f"ref_wcx_mom_{filter_type}"] = paths["ref_lab"] / "WCX" / f"mom_{config['WCX']['ref']}_{filter_type}.npz"
        if not run_pipeline_step(11, "Run WC/WCX/WCFF prediction (mom)", run_wcfamily_prediction, progress, sample_name, paths, "mom", filter_type, gender):
            progress.mark_failed("Run WC/WCX/WCFF prediction (mom) failed")
            return 1

    # -------------------------------------------
    # MD & Trisomy Detection with WC family
    # let the final output file with WCX mom
    # -------------------------------------------
    output_filename_md = f"{sample_name}_WCX_mom_md320.tsv"
    md_result_file = os.path.join(ANALYSIS_DIR, f"{sample_name}", "Output_WCX", "mom", f"{output_filename_md}")
    md_completed = os.path.exists(md_result_file)

    if not md_completed:
        log_and_print("=== Starting MD Analysis ===")
        if not run_pipeline_step(12, "Microdeletion Analysis", run_microdeletion_decision_pipeline, progress, sample_name, labcode, config, ANALYSIS_DIR, OUTPUT_DIR, bed_dir):
            progress.mark_failed("Microdeletion analysis failed")
            return 1
    else:
        progress.update_step(12, "Microdeletion Analysis", "SKIP", "already completed")

    # -------------------------------------------
    # EZD
    # -------------------------------------------
    ezd_mom_result = os.path.join(ANALYSIS_DIR, f"{sample_name}", "Output_EZD", "mom",  "Trisomy_detect_result_mom_with_SCA.tsv")
    ezd_completed = os.path.exists(ezd_mom_result)

    if not ezd_completed:
        log_and_print("=== Starting EZD Analysis ===")
        ezd_input_paths = {
            "orig": f"{ANALYSIS_DIR}/{sample_name}/Output_hmmcopy/{sample_name}.proper_paired.50kb.wig.Normalization.txt",
            "fetus": f"{ANALYSIS_DIR}/{sample_name}/Output_hmmcopy/{sample_name}.of_fetus.50kb.wig.Normalization.txt",
            "mom": f"{ANALYSIS_DIR}/{sample_name}/Output_hmmcopy/{sample_name}.of_mom.50kb.wig.Normalization.txt",
        }

        groups = ['orig', 'fetus', 'mom']
        ezd_success = True

        for idx, group in enumerate(groups, 1):
            step_num = f"13.{idx}"
            wig_path = ezd_input_paths[group]

            try:
                ezd_df, decision_df = run_ezd_group(
                    sample_name=sample_name,
                    group=group,
                    wig_path=wig_path,
                    labcode=labcode,
                    analysis_dir=ANALYSIS_DIR,
                    data_dir=DATA_DIR
                )

                progress.update_step(step_num, f"EZD {group.upper()}", "PASS")
                log_and_print(f"\n=== {group.upper()} 그룹 결과 ===")
                log_and_print("\nDecision Results:")
                log_and_print(decision_df)

            except Exception as e:
                progress.update_step(step_num, f"EZD {group.upper()}", "FAIL")
                log_and_print(f"Error processing {group}: {e}", 'WARNING')
                ezd_success = False
                break

        if not ezd_success:
            progress.mark_failed("EZD analysis failed")
            return 1

        progress.update_step(13, "EZD Analysis", "PASS")
    else:
        progress.update_step(13, "EZD Analysis", "SKIP", "already completed")

    # -------------------------------------------
    # PRIZM
    # -------------------------------------------
    prizm_mom_result = os.path.join(ANALYSIS_DIR, f"{sample_name}", "Output_PRIZM", "mom",  f"{sample_name}_mom.trisomy_detection.tsv")
    prizm_completed = os.path.exists(prizm_mom_result)

    if not prizm_completed:
        log_and_print("=== Starting PRIZM Analysis ===")
        
        if run_multiple_prizm_analysis is None:
            log_and_print("PRIZM module not available. Skipping PRIZM analysis.", 'WARNING')
            progress.update_step(14, "PRIZM Analysis", "SKIP", "module not available")
        else:
            try:
                prizm_success = run_multiple_prizm_analysis(sample_name, gender, labcode, config, ANALYSIS_DIR, DATA_DIR)
                
                if prizm_success:
                    progress.update_step(14, "PRIZM Analysis", "PASS")
                    log_and_print("[PRIZM] Analysis completed.")
                else:
                    progress.update_step(14, "PRIZM Analysis", "FAIL")
                    log_and_print("[PRIZM] Analysis failed.", 'ERROR')
                    # 실패 시 파이프라인 중단할지 결정
                    progress.mark_failed("PRIZM analysis failed")
                    return 1
                    
            except Exception as e:
                progress.update_step(14, "PRIZM Analysis", "FAIL")
                log_and_print(f"[PRIZM] Analysis exception: {e}", 'ERROR')
                # 예외 시 파이프라인 중단할지 결정
                progress.mark_failed(f"PRIZM analysis error: {str(e)}")
                return 1
    else:
        progress.update_step(14, "PRIZM Analysis", "SKIP", "already completed")

    # Copy files to Output directory
    success = copy_outputs_to_final_dir(sample_name, ANALYSIS_DIR, OUTPUT_DIR)
    if success:
        progress.update_step(15, "Copying files to Output directory", "PASS")
    else:
        progress.update_step(15, "Copying files to Output directory", "FAIL")
        log_and_print(f"Failed to copy outputs for {sample_name}", "ERROR")
        return False

    # Previous style
    #copy_outputs_to_final_dir(sample_name, ANALYSIS_DIR, OUTPUT_DIR)

    # Json Output
    json_result = os.path.join(output_dir, f"{sample_name}", f"{sample_name}.json")
    json_completed = os.path.exists(json_result)
    if not json_completed:
        log_and_print("=== Starting Json Generation ===")
        try:
            json_output_path = build_nipt_json(ANALYSIS_DIR, OUTPUT_DIR, f"{DATA_DIR}/refs/{labcode}", sample_name, age, VERSION, f"{bed_dir}/common")
            log_and_print(f"[JSON] Output json file saved to: {json_output_path}")
            progress.update_step(16, "Json Output Generation", "PASS")
        except Exception as e:
            log_and_print(f"[JSON] Failed to generate output json: {e}", "ERROR")
            progress.update_step(16, "Json Output Generation", "FAIL")
    else:
        progress.update_step(16, "Json Output Generation", "SKIP", "already completed")

    # HTML generation - 완전히 독립적으로 실행
    html_result = os.path.join(output_dir, f"{sample_name}", f"{sample_name}.html")
    html_completed = os.path.exists(html_result)
    if not html_completed:
        try:
            html_report_path = generate_nipt_html_report(
                json_file_path=json_result,
                output_dir=OUTPUT_DIR
                )

            if html_report_path:
                log_and_print(f"[HTML] HTML report available at: {html_report_path}")
                progress.update_step(17, "HTML Output Generation", "PASS")
            else:
                log_and_print(f"[HTML] HTML report generation was skipped or failed", 'WARNING')
                progress.update_step(17, "HTML Output Generation", "FAIL")

        except Exception as e:
            log_and_print(f"[HTML] Unexpected error in HTML generation: {e}", 'ERROR')
    else:
        progress.update_step(17, "HTML Output Generation", "SKIP", "already completed")

    with open(final_marker, "w") as marker_file:
        marker_file.write(f"NIPT Pipeline completed successfully for {sample_name} at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    log_and_print(f"NIPT pipeline completed successfully : {sample_name}")
    progress.mark_completed()

if __name__ == "__main__":
    main()
