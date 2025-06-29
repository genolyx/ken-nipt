# ------------------------------------------------------------------------------------
#   Docker batch process
#   Author : Hyukjung Kwon
#   Date : 2025. 06. 01
# ------------------------------------------------------------------------------------

"""
---------------------------------------------
Docker batch process

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
def update_process_status(process_status):
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
                        process_status[sample_name] = "Completed"
                        containers_to_remove.append(sample_name)
                        logger.info(f"Sample {sample_name} completed, will be removed")
                    else:
                        process_status[sample_name] = current_status

        # 완료된 컨테이너들 정리
        for container_name in containers_to_remove:
            logger.info(f"{container_name} will be removed")
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
    
    logger.info(f"Completed: {completed_count}, Failed: {failed_count}, Total: {sample_number}")
    
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

            update_process_status(process_status)

        except Exception as e:
            logger.error(f"Error running : {e}")
            process_status[next_sample] = "Failed"

        return True

    return False

# Main function
def run_daemon(cmd, labcode, sample_list, max_samples, intv):

    # Read sample list from file
    samples_dic = {}

    try:
        sample_df = pd.read_csv(sample_list, sep='\t', dtype={"SAMPLE_NAME": str})
        sample_df = sample_df[["SAMPLE_NAME", "WORK_DIR", "FQ1", "FQ2", "AGE"]]
        samples_dic = {row["SAMPLE_NAME"]: [row['WORK_DIR'], row["FQ1"], row["FQ2"], row["AGE"]] for _, row in sample_df.iterrows()}

        logger.info(f"Loaded {len(samples_dic)} samples from {sample_list}")

    except Exception as e:
        logger.error(f"Error reading sample list: {e}")
        return

    total_samples = len(samples_dic)
    process_status = {sample: "Init" for sample in samples_dic}
    print_process_status(process_status)

    # Continuously check if additional samples need to be started
    while True:
        # 새 샘플 시작 시도
        start_sample(samples_dic, cmd, labcode, max_samples, process_status)
        
        # 프로세스 상태 업데이트
        update_process_status(process_status)
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

    parser = argparse.ArgumentParser(description="NIPT Sample Processing Daemon")

    parser.add_argument('-l', '--labcode', required=True, help='Lab code')
    parser.add_argument('-s', '--sample_list', required=True, help='Path to sample list TSV file')
    parser.add_argument('-m', '--max_samples', type=int, default=2, help='Maximum concurrent samples (default: 2)')
    parser.add_argument('-i', '--poll_intv', type=int, default=5, help='Polling interval in seconds (default: 5)')

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
        logger.info("Run daemon...")

        #command = "bash /home/ken/ken-nipt/src/run_nipt_dev.sh -s SAMPLE_NAME -1 FQ1 -2 FQ2 -l LABCODE -root /home/ken/ken-nipt -work WORK_DIR --detached"
        #command = "bash /home/ken/ken-nipt/src/run_nipt_v3.sh -s SAMPLE_NAME -1 FQ1 -2 FQ2 -a AGE -l LABCODE -root /home/ken/ken-nipt -work WORK_DIR --detached -cf -ao"
        command = "bash /home/ken/ken-nipt/src/run_nipt_v3.sh -s SAMPLE_NAME -1 FQ1 -2 FQ2 -a AGE -l LABCODE -root /home/ken/ken-nipt -work WORK_DIR --detached -cf"

        logger.info(f"Starting NIPT daemon with:")
        logger.info(f"  Lab code: {args.labcode}")
        logger.info(f"  Sample list: {args.sample_list}")
        logger.info(f"  Max samples: {args.max_samples}")
        logger.info(f"  Poll interval: {args.poll_intv}s")

        run_daemon(command, args.labcode, args.sample_list, args.max_samples, args.poll_intv)

    except Exception as e:
        logger.exception("Error occurred during pipeline execution.")
