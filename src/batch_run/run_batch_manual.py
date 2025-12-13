# ------------------------------------------------------------------------------------
#   Docker batch process (Manual version)
#   Author : Hyukjung Kwon
#   Date : 2025. 06. 01
# ------------------------------------------------------------------------------------

"""
---------------------------------------------
Docker batch process (Manual version)

Author: {author}
Contact: {email}
---------------------------------------------
"""

import subprocess
import argparse
import logging
import threading
import time
import os
import sys
import pandas as pd

__author__ = 'Hyukjung Kwon'
__version__ = '0.1'


# Function to parse the output of 'docker ps -a' and update the process status dictionary
def update_process_status(process_status, samples_dic=None):
    try:
        result = subprocess.run(["docker", "ps", "-a", "--format", "{{.Image}}|{{.Names}}|{{.Status}}"], capture_output=True, text=True)
        lines = result.stdout.split('\n')

        containers_to_remove = []  # 제거할 컨테이너 목록

        for line in lines:
            if line:
                image, sample_name, status = line.split('|')

                if sample_name in process_status:
                    # 상태 업데이트 로직 개선
                    current_status = status.strip().split(' ')[0]
                    if current_status == "Up":
                        process_status[sample_name] = "Running"
                    elif current_status == "Exited":
                        # Exited 상태인데 아직 Completed로 표시되지 않은 경우만 처리
                        if process_status[sample_name] != "Completed":
                            process_status[sample_name] = "Completed"
                            containers_to_remove.append(sample_name)
                            logger.info(f"Sample {sample_name} completed, will be removed")
                    else:
                        process_status[sample_name] = current_status

        # 완료된 컨테이너들 정리 및 analysis directory 제거
        for container_name in containers_to_remove:
            logger.info(f"{container_name} will be removed")
            
            # 분석 완료 후 analysis directory 제거 (v2 기능)
            # output에 결과가 정상적으로 생성되었는지 확인 후 제거
            if samples_dic and container_name in samples_dic:
                root_dir = "/home/ken/ken-nipt"
                work_dir = samples_dic[container_name][0]
                output_dir = f"{root_dir}/output/{work_dir}"
                analysis_dir = f"{root_dir}/analysis/{work_dir}"
                json_file = f"{output_dir}/{container_name}/{container_name}.json"
                analysis_sample_dir = f"{analysis_dir}/{container_name}"
                
                # JSON 파일이 존재하면 analysis directory 제거
                if os.path.exists(json_file) and os.path.isdir(analysis_sample_dir):
                    logger.info(f"Removing analysis directory: {analysis_sample_dir}")
                    try:
                        import shutil
                        shutil.rmtree(analysis_sample_dir)
                        logger.info(f"Successfully removed analysis directory for {container_name}")
                    except Exception as e:
                        logger.warning(f"Failed to remove analysis directory for {container_name}: {e}")
            
            try:
                remove_result = subprocess.run(["docker", "rm", container_name],
                                             capture_output=True, text=True, timeout=10)
                if remove_result.returncode == 0:
                    logger.info(f"Successfully removed container: {container_name}")
                else:
                    logger.warning(f"Failed to remove container {container_name}: {remove_result.stderr}")
            except Exception as e:
                logger.error(f"Error removing container {container_name}: {e}")

    except subprocess.CalledProcessError as e:
        logger.error(f"Error running : {e}")

# Function to print the process status dictionary
def print_process_status(process_status):
    logger.info("===================================")
    for sample_name, status in process_status.items():
        logger.info(f"{sample_name}: {status}")
    logger.info("===================================")

# Function to check if a sample is already running
def is_sample_running(sample_name):
    try:
        result = subprocess.run(["docker", "ps", "-a", "--format", "{{.Names}}"], capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running : {e}")

    return sample_name in result.stdout.split()

def all_finished(process_status, sample_number):
    """모든 샘플이 완료되었는지 확인"""
    completed_count = sum(1 for status in process_status.values() if status == "Completed")
    failed_count = sum(1 for status in process_status.values() if status == "Failed")
    running_count = sum(1 for status in process_status.values() if status == "Running")
    init_count = sum(1 for status in process_status.values() if status == "Init")

    logger.info(f"Running: {running_count}, Completed: {completed_count}, Failed: {failed_count}, Init: {init_count}, Total: {sample_number}")

    return (completed_count + failed_count) == sample_number

def get_next_key(process_status):
    """다음에 실행할 샘플 찾기"""
    for sample_name, status in process_status.items():
        if status == "Init":
            return sample_name
    return "Anymore"

# Function to check if additional samples need to be started
def need_additional_samples(max_samples, process_status):
    """추가로 실행할 수 있는 샘플이 있는지 확인"""
    # 현재 실행 중인 샘플 수 계산
    running_count = sum(1 for status in process_status.values() if status == "Running")
    logger.info(f"Current running sample count is {running_count}")

    # 사용 가능한 슬롯 계산
    avail_sample = max_samples - running_count

    if avail_sample > 0:
        # 다음에 실행할 샘플이 있는지 확인
        next_sample = get_next_key(process_status)
        if next_sample != "Anymore":
            return next_sample, avail_sample

    return None, avail_sample

# Function to start additional samples
def start_sample(samples_dic, cmd, labcode, max_samples, process_status):
    """새 샘플 시작"""
    next_sample, avail_sample = need_additional_samples(max_samples, process_status)

    if next_sample is None:
        if get_next_key(process_status) == "Anymore":
            logger.info(f"No more samples to be run")
        else:
            logger.info(f"Maximum capacity reached ({max_samples}), waiting for slots...")
        return False

    if avail_sample > 0:
        logger.info(f"[start_sample] {next_sample}")

        try:
            cmd_to_run = cmd.replace('SAMPLE_NAME', next_sample)
            cmd_to_run = cmd_to_run.replace('WORK_DIR', str(samples_dic[next_sample][0]))
            cmd_to_run = cmd_to_run.replace('FQ1', samples_dic[next_sample][1])
            cmd_to_run = cmd_to_run.replace('FQ2', samples_dic[next_sample][2])
            cmd_to_run = cmd_to_run.replace('AGE', str(samples_dic[next_sample][3]))
            cmd_to_run = cmd_to_run.replace('LABCODE', labcode)

            logger.warning(f"Running : {cmd_to_run}")

            try:
                # Execute the command
                result = subprocess.run(cmd_to_run.split(' '), shell=False, capture_output=True, text=True, timeout=60)

                logger.info(f"Return code: {result.returncode}")
                logger.info(f"STDOUT: {result.stdout}")
                logger.info(f"STDERR: {result.stderr}")

                if result.returncode == 0:
                    # 스킵된 경우 (이미 완료된 샘플) 확인
                    if "SKIPPING" in result.stdout or "Already completed" in result.stdout:
                        logger.info(f"Sample {next_sample} already completed, skipping")
                        # output 디렉토리에 JSON 파일과 output.tar 파일이 있는지 확인
                        root_dir = "/home/ken/ken-nipt"
                        work_dir = samples_dic[next_sample][0]
                        output_dir = f"{root_dir}/output/{work_dir}/{next_sample}"
                        json_file = f"{output_dir}/{next_sample}.json"
                        tar_file = f"{output_dir}/{next_sample}.output.tar"
                        
                        # Both json and output.tar must exist
                        if os.path.exists(json_file) and os.path.exists(tar_file):
                            process_status[next_sample] = "Completed"
                            logger.info(f"Sample {next_sample} marked as Completed (json + output.tar exist)")
                        else:
                            # 파일이 없으면 다시 실행해야 함
                            missing = []
                            if not os.path.exists(json_file):
                                missing.append("json")
                            if not os.path.exists(tar_file):
                                missing.append("output.tar")
                            logger.warning(f"Sample {next_sample} has marker but missing files: {', '.join(missing)}, will retry")
                            process_status[next_sample] = "Init"  # 다시 시도
                    else:
                        logger.info(f"Command completed successfully for {next_sample}")
                        process_status[next_sample] = "Running"  # 실행 중으로 상태 변경
                else:
                    logger.error(f"Command failed with return code {result.returncode}")
                    logger.error(f"Error output: {result.stderr}")
                    process_status[next_sample] = "Failed"  # 실패로 상태 변경

            except subprocess.TimeoutExpired:
                logger.error(f"Command timeout for {next_sample}")
                process_status[next_sample] = "Failed"
            except subprocess.CalledProcessError as e:
                logger.error(f"Error running subprocess: {e}")
                process_status[next_sample] = "Failed"
            except Exception as e:
                logger.error(f"Unexpected error: {e}")
                process_status[next_sample] = "Failed"

            update_process_status(process_status, samples_dic)

        except Exception as e:
            logger.error(f"Error running : {e}")
            process_status[next_sample] = "Failed"

        return True

    return False

# Main function
def run_daemon(cmd, labcode, sample_list, max_samples, intv, test_num=None):

    # Read sample list from file
    samples_dic = {}

    try:
        sample_df = pd.read_csv(sample_list, sep='\t', dtype={"SAMPLE_NAME": str})
        sample_df = sample_df[["SAMPLE_NAME", "WORK_DIR", "FQ1", "FQ2", "AGE"]]

        # Limit to test_num samples if specified
        if test_num is not None and test_num > 0:
            sample_df = sample_df.head(test_num)
            logger.info(f"TEST MODE: Limiting to first {test_num} samples")

        samples_dic = {row["SAMPLE_NAME"]: [row['WORK_DIR'], row["FQ1"], row["FQ2"], row["AGE"]] for _, row in sample_df.iterrows()}

        logger.info(f"Loaded {len(samples_dic)} samples from {sample_list}")

    except Exception as e:
        logger.error(f"Error reading sample list: {e}")
        return

    total_samples = len(samples_dic)
    process_status = {sample: "Init" for sample in samples_dic}
    
    # Check for already completed samples (json + output.tar exist)
    logger.info("Checking for already completed samples...")
    root_dir = "/home/ken/ken-nipt"
    completed_initially = 0
    for sample_name in samples_dic:
        work_dir = samples_dic[sample_name][0]
        output_dir = f"{root_dir}/output/{work_dir}/{sample_name}"
        json_file = f"{output_dir}/{sample_name}.json"
        tar_file = f"{output_dir}/{sample_name}.output.tar"
        
        # Check if both json and output.tar exist
        if os.path.exists(json_file) and os.path.exists(tar_file):
            logger.info(f"Sample {sample_name} already completed (json + output.tar exist)")
            process_status[sample_name] = "Completed"
            completed_initially += 1
    
    logger.info(f"Found {completed_initially} already completed samples")
    print_process_status(process_status)

    # Continuously check if additional samples need to be started
    while True:
        # 새 샘플 시작 시도
        start_sample(samples_dic, cmd, labcode, max_samples, process_status)

        # 프로세스 상태 업데이트
        update_process_status(process_status, samples_dic)
        print_process_status(process_status)

        # 모든 샘플이 완료되었는지 확인
        if all_finished(process_status, total_samples):
            logger.info("All samples have been processed!")
            break

        time.sleep(intv)  # 지정된 간격으로 대기

    # 최종 결과 출력
    completed_count = sum(1 for status in process_status.values() if status == "Completed")
    failed_count = sum(1 for status in process_status.values() if status == "Failed")

    logger.info("=" * 50)
    logger.info("FINAL RESULTS:")
    logger.info(f"Total samples: {total_samples}")
    logger.info(f"Completed: {completed_count}")
    logger.info(f"Failed: {failed_count}")
    logger.info("=" * 50)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="NIPT Sample Processing Daemon (Manual version)")

    parser.add_argument('-l', '--labcode', required=True, help='Lab code')
    parser.add_argument('-s', '--sample_list', required=True, help='Path to sample list TSV file')
    parser.add_argument('-m', '--max_samples', type=int, default=2, help='Maximum concurrent samples (default: 2)')
    parser.add_argument('-i', '--poll_intv', type=int, default=5, help='Polling interval in seconds (default: 5)')
    parser.add_argument('--test-num', type=int, default=None, help='Test mode: process only first N samples (default: process all)')

    parser.add_argument("--log_stdout", default="Y", choices=["Y", "N"], help="Pipeline logging to Standard output (default: NO)")

    args = parser.parse_args()

    # Setup logging
    if args.log_stdout == "N":
        logging.basicConfig(filename=args.sample_list+".log", level=logging.INFO, format="%(asctime)s [%(levelname) 7s] (%(lineno)4d) | %(message)s")
    else:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname) 7s] (%(lineno)4d) | %(message)s")

    logging.addLevelName(logging.ERROR, "\033[91m%7s\033[0m" % logging.getLevelName(logging.ERROR))
    logging.addLevelName(logging.WARNING, "\033[93m%7s\033[0m" % logging.getLevelName(logging.WARNING))
    logger = logging.getLogger(__name__)

    try:
        logger.info("Run daemon (Manual version)...")

        # Use run_nipt_manual_v2.sh instead of run_nipt_manual.sh
        command = "bash /home/ken/ken-nipt/src/batch_run/run_nipt_manual_v2.sh -s SAMPLE_NAME -1 FQ1 -2 FQ2 -a AGE -l LABCODE -root /home/ken/ken-nipt -work WORK_DIR -cf --detached"

        logger.info(f"Starting NIPT daemon with:")
        logger.info(f"  Lab code: {args.labcode}")
        logger.info(f"  Sample list: {args.sample_list}")
        logger.info(f"  Max samples: {args.max_samples}")
        logger.info(f"  Poll interval: {args.poll_intv}s")
        if args.test_num:
            logger.info(f"  Test mode: processing first {args.test_num} samples only")

        run_daemon(command, args.labcode, args.sample_list, args.max_samples, args.poll_intv, args.test_num)

    except Exception as e:
        logger.exception("Error occurred during pipeline execution.")

