# ------------------------------------------------------------------------------------
#   Docker batch process
#   Author : Kenneth Kwon
#   Date : 2024. 05. 13
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

__author__ = 'Kenneth Kwon'
__email__ = "kenneth@ploidyone.com"
__version__ = '0.1'


# Function to parse the output of 'docker ps -a' and update the process status dictionary
def update_process_status(process_status):
    try:
        result = subprocess.run(["docker", "ps", "-a", "--format", "{{.Image}}|{{.Names}}|{{.Status}}"], capture_output=True, text=True)
        lines = result.stdout.split('\n')
        for line in lines:
            if line:
                image, sample_name, status = line.split('|')

                #print(image, sample_name, status)
                #print(process_status)

                if sample_name in process_status:
                    process_status[sample_name] = status.strip().split(' ')[0]

    except subprocess.CalledProcessError as e:
        logger.error(f"Error running : {e}")

    #logger.info(process_status)

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
    exit_count = sum("Exited" in status for status in process_status.values())

    if exit_count == sample_number:
        return True
    else:
        return False

def get_next_key(process_status):
    for sample_name, status in process_status.items():
        if status == "Init":
            return sample_name
    return "Anymore"

# Function to check if additional samples need to be started
def need_additional_samples(max_samples, process_status):
    running_count = sum("Up" in status for status in process_status.values())
    logger.info(f"Current running sample count is {running_count}")

    avail_sample = max_samples - running_count
    if avail_sample > 0:
        return get_next_key(process_status), avail_sample
    else:
        return None, avail_sample 

# Function to start additional samples
def start_sample(samples_dic, cmd, max_samples, process_status):
    next_sample, avail_sample = need_additional_samples(max_samples, process_status)

    if next_sample == "Anymore":
        return False

    if avail_sample > 0:
        logger.info(f"[start_sample] {next_sample}")

        try:
            cmd = cmd.replace('SAMPLE_NAME', next_sample)
            cmd = cmd.replace('FQ1', samples_dic[next_sample][0])
            cmd = cmd.replace('FQ2', samples_dic[next_sample][1])
            #cmd = cmd.replace('LAB', samples_dic[next_sample][2])

            logger.warning(f"Running : {cmd}")

            try:
                subprocess.run(cmd.split(' '))
            except subprocess.CalledProcessError as e:
                logger.error(f"Error running : {e}")

            update_process_status(process_status)

        except subprocess.CalledProcessError as e:
            logger.error(f"Error running : {e}")

        return True


# Main function
def run_daemon(cmd, sample_list, max_samples, intv):

    # Read sample list from file
    samples_dic = {}

    #sample_df = pd.read_excel(sample_list)
    sample_df = pd.read_csv(sample_list, sep='\t', dtype={"SAMPLE_NAME": str})
    #sample_df = sample_df[["SAMPLE_NAME", "FQ1", "FQ2", "LAB"]]
    sample_df = sample_df[["SAMPLE_NAME", "FQ1", "FQ2"]]

    #samples_dic = {row["SAMPLE_NAME"]: [row["FQ1"], row["FQ2"], row["LAB"]] for _, row in sample_df.iterrows()}
    samples_dic = {row["SAMPLE_NAME"]: [row["FQ1"], row["FQ2"]] for _, row in sample_df.iterrows()}

    process_status = {sample: "Init" for sample in samples_dic}
    print_process_status(process_status)

    # Continuously check if additional samples need to be started
    while True:
        retv = start_sample(samples_dic, cmd, max_samples, process_status)
        if retv == False:
            logger.warning("No remained samples to be run!")
            break

        update_process_status(process_status)
        print_process_status(process_status)
        time.sleep(intv)  # Sleep for 5 seconds before checking again


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Docker batch process")
    parser.add_argument("max_samples", type=int, help="Max number of samples running simultaenously")
    parser.add_argument("sample_list", help="Sample list file containing sample name only")
    parser.add_argument("command", help="bash command for docker run")

    # Kenneth : I should make it as optional argument
    parser.add_argument("--poll_intv", type=int, default=5, help="polling interval sec (default: 5)")
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
        run_daemon(args.command, args.sample_list, args.max_samples, args.poll_intv)

    except Exception as e:
        logger.exception("Error occurred during pipeline execution.")

