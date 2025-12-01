#!/usr/bin/env python3
# ------------------------------------------------------------------------------------
#   MD batch process for artificial samples
#   Author : Based on run_batch_dev.py by Hyukjung Kwon
#   Date : 2025. 11. 02
# ------------------------------------------------------------------------------------

"""
---------------------------------------------
MD batch process for artificial samples

Scans batch_output directory and runs MD pipeline for each sample
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
import json
import csv
from pathlib import Path

__author__ = 'Based on run_batch_dev.py'
__version__ = '0.1'


def get_gender_from_json(sample_dir, sample_name):
    """Read gender from sample JSON file"""
    json_file = sample_dir / f"{sample_name}.json"
    if json_file.exists():
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
                gender = data.get('gender', '').upper()
                if gender in ['M', 'MALE', 'XY']:
                    return 'MALE'
                elif gender in ['F', 'FEMALE', 'XX']:
                    return 'FEMALE'
        except Exception as e:
            logger.debug(f"Could not read gender from {json_file}: {e}")
    return 'UNKNOWN'


def check_all_results_exist(root_dir, work_dir, sample_name, filter_type="of"):
    """Check if all required output files exist for a sample
    
    Checks for:
    - Marker file: {sample_name}.md_pipeline_completed.marker (checked first)
    - Output_FF: Based on gender
      - Female: seqFF only required
      - Male: both YFF and seqFF required
      - Unknown: seqFF only required
    - Output_WC/orig: wc.orig.out.npz and wc.orig.report.txt
    - Output_WC/fetus: wc.fetus.out.npz and wc.fetus.report.txt
    - Output_WCX/orig: wcx.proper_paired.npz, wcx.orig_aberrations.bed, wcx.orig.plots
    - Output_WCX/fetus: wcx.of_fetus.npz, wcx.fetus_aberrations.bed, wcx.fetus.plots
    
    If any result file is missing, returns False (sample needs processing).
    The pipeline will automatically create required BAM files (of_orig.bam, of_fetus.bam)
    from proper_paired.bam if they don't exist.
    
    Args:
        root_dir: Root directory
        work_dir: Work directory
        sample_name: Sample name
        filter_type: Filter type (default: "of", reserved for future use)
    
    Returns:
        bool: True if all required outputs exist (skip processing), False otherwise (needs processing)
    """
    sample_dir = Path(root_dir) / "analysis" / work_dir / sample_name
    
    # First, check for marker file (indicates pipeline completed successfully)
    marker_file = sample_dir / f"{sample_name}.md_pipeline_completed.marker"
    if marker_file.exists():
        logger.debug(f"Marker file found for {sample_name}, considering as completed")
        return True
    
    # Get gender from JSON file
    gender = get_gender_from_json(sample_dir, sample_name)
    
    # Check Output_FF - requirements depend on gender
    ff_dir = sample_dir / "Output_FF"
    ff_yff_exists = False
    ff_seqff_exists = False
    
    if ff_dir.exists():
        # Check seqFF file
        seqff_file = ff_dir / f"{sample_name}.seqff.txt"
        if seqff_file.exists():
            try:
                df = pd.read_csv(seqff_file, index_col=0)
                if "SeqFF" in df.index:
                    ff_seqff_exists = True
            except:
                pass
        
        # Check YFF - check JSON file first (preferred method)
        json_file = sample_dir / f"{sample_name}.json"
        if json_file.exists():
            try:
                with open(json_file, 'r') as f:
                    metadata = json.load(f)
                
                # Check calculated_ff section in JSON
                calculated_ff = metadata.get("calculated_ff", {})
                yff_status = calculated_ff.get("yff_status")
                
                # YFF exists if yff_status is "OK" (calculation completed)
                if yff_status == "OK":
                    ff_yff_exists = True
            except:
                pass
        
        # Fallback: also check fetal_fraction.txt file if JSON doesn't have YFF
        if not ff_yff_exists:
            ff_txt = ff_dir / f"{sample_name}.fetal_fraction.txt"
            if ff_txt.exists():
                try:
                    with open(ff_txt, 'r') as f:
                        content = f.read()
                        if "YFF" in content or "yff" in content:
                            try:
                                ff_data = json.loads(content)
                                if ff_data.get("yff", 0) > 0 or ff_data.get("YFF", 0) > 0:
                                    ff_yff_exists = True
                            except:
                                if "YFF" in content and "0.00" not in content.split("YFF")[1][:10]:
                                    ff_yff_exists = True
                except:
                    pass
    
    # Check FF completeness based on gender
    if gender == 'FEMALE':
        # Female: seqFF only required
        ff_complete = ff_seqff_exists
    elif gender == 'MALE':
        # Male: both YFF and seqFF required
        ff_complete = ff_yff_exists and ff_seqff_exists
    else:
        # Unknown: at least seqFF should exist
        ff_complete = ff_seqff_exists
    
    # Check Output_WC/orig
    wc_orig_npz = sample_dir / "Output_WC" / "orig" / f"{sample_name}.wc.orig.out.npz"
    wc_orig_report = sample_dir / "Output_WC" / "orig" / f"{sample_name}.wc.orig.report.txt"
    wc_orig_complete = wc_orig_npz.exists() and wc_orig_report.exists()
    
    # Check Output_WC/fetus
    wc_fetus_npz = sample_dir / "Output_WC" / "fetus" / f"{sample_name}.wc.fetus.out.npz"
    wc_fetus_report = sample_dir / "Output_WC" / "fetus" / f"{sample_name}.wc.fetus.report.txt"
    wc_fetus_complete = wc_fetus_npz.exists() and wc_fetus_report.exists()
    
    # Check Output_WCX/orig
    wcx_orig_npz = sample_dir / "Output_WCX" / f"{sample_name}.wcx.proper_paired.npz"
    wcx_orig_bed = sample_dir / "Output_WCX" / "orig" / f"{sample_name}.wcx.orig_aberrations.bed"
    wcx_orig_plots = sample_dir / "Output_WCX" / "orig" / f"{sample_name}.wcx.orig.plots"
    wcx_orig_complete = wcx_orig_npz.exists() and wcx_orig_bed.exists() and wcx_orig_plots.exists()
    
    # Check Output_WCX/fetus
    wcx_fetus_npz = sample_dir / "Output_WCX" / f"{sample_name}.wcx.of_fetus.npz"
    wcx_fetus_bed = sample_dir / "Output_WCX" / "fetus" / f"{sample_name}.wcx.fetus_aberrations.bed"
    wcx_fetus_plots = sample_dir / "Output_WCX" / "fetus" / f"{sample_name}.wcx.fetus.plots"
    wcx_fetus_complete = wcx_fetus_npz.exists() and wcx_fetus_bed.exists() and wcx_fetus_plots.exists()
    
    # All components must be complete
    # If any result is missing, the sample needs processing (pipeline will create BAM files if needed)
    all_complete = ff_complete and wc_orig_complete and wc_fetus_complete and wcx_orig_complete and wcx_fetus_complete
    
    if all_complete:
        logger.debug(f"All results exist for {sample_name}: FF={ff_complete} (yff={ff_yff_exists}, seqff={ff_seqff_exists}), "
                     f"WC/orig={wc_orig_complete}, WC/fetus={wc_fetus_complete}, "
                     f"WCX/orig={wcx_orig_complete}, WCX/fetus={wcx_fetus_complete}")
    
    return all_complete


# Function to parse the output of 'docker ps -a' and update the process status dictionary
def update_process_status(process_status, samples_dic, root_dir, force=False):
    """Update process status based on Docker containers and result files
    
    Args:
        force: If True, ignore existing result files and only check container status
    """
    try:
        # First, check result files for completion (skip if force=True)
        if not force:
            for sample_name, status in process_status.items():
                if status in ["Completed", "Failed"]:
                    continue  # Skip already completed/failed samples
                
                # Get sample info
                if sample_name not in samples_dic:
                    continue
                
                sample_info = samples_dic[sample_name]
                work_dir = sample_info['work_dir']
                
                # Check if all required result files exist
                if check_all_results_exist(root_dir, work_dir, sample_name):
                    logger.info(f"All result files found for {sample_name}, marking as Completed")
                    process_status[sample_name] = "Completed"
                    continue
        
        # Then check Docker container status
        result = subprocess.run(["docker", "ps", "-a", "--format", "{{.Image}}|{{.Names}}|{{.Status}}"], capture_output=True, text=True)
        lines = result.stdout.split('\n')
        
        # Get list of existing container names for checking removed containers
        existing_containers = set()
        for line in lines:
            if line:
                parts = line.split('|')
                if len(parts) >= 2:
                    existing_containers.add(parts[1])
        
        # Check for samples with containers that were already removed
        # If container was removed but sample is still Running, check marker/result files
        for sample_name, status in process_status.items():
            if status == "Running":
                container_name = f"md_{sample_name}"
                if container_name not in existing_containers:
                    # Container was removed, check for completion markers
                    if sample_name in samples_dic:
                        sample_info = samples_dic[sample_name]
                        work_dir = sample_info['work_dir']
                        
                        # Check all result files (marker file is not sufficient, need actual results)
                        if not force and check_all_results_exist(root_dir, work_dir, sample_name):
                            logger.info(f"Container removed and all results found for {sample_name}, marking as Completed")
                            process_status[sample_name] = "Completed"
                        elif force:
                            # In force mode, if container was removed, assume completed (might be risky)
                            logger.info(f"Container removed for {sample_name} (force mode), marking as Completed")
                            process_status[sample_name] = "Completed"

        containers_to_remove = []  # 제거할 컨테이너 목록

        for line in lines:
            if line:
                parts = line.split('|')
                if len(parts) < 3:
                    continue
                image, container_name, status = parts[0], parts[1], '|'.join(parts[2:])
                
                # Container name format: md_{sample_id}
                if container_name.startswith("md_"):
                    sample_name = container_name[3:]  # Remove "md_" prefix
                    
                    if sample_name in process_status:
                        # Skip if already completed
                        if process_status[sample_name] == "Completed":
                            containers_to_remove.append(container_name)
                            continue
                        
                        # 상태 업데이트 로직
                        current_status = status.strip().split(' ')[0]
                        if current_status == "Up":
                            # Container is running
                            if process_status[sample_name] != "Completed":
                                old_status = process_status[sample_name]
                                process_status[sample_name] = "Running"
                                if old_status != "Running":
                                    logger.info(f"Status updated: {sample_name} {old_status} -> Running (container: {container_name})")
                        elif current_status == "Exited":
                            # Container exited, check exit code
                            exit_code_match = subprocess.run(
                                ["docker", "inspect", "--format", "{{.State.ExitCode}}", container_name],
                                capture_output=True, text=True
                            )
                            if exit_code_match.returncode == 0:
                                exit_code = exit_code_match.stdout.strip()
                                if exit_code == "0":
                                    # Check result file or marker file (skip if force=True)
                                    if sample_name in samples_dic:
                                        sample_info = samples_dic[sample_name]
                                        work_dir = sample_info['work_dir']
                                        if force:
                                            # In force mode, only check container exit code
                                            process_status[sample_name] = "Completed"
                                            containers_to_remove.append(container_name)
                                        else:
                                            # Check all result files (marker file is not sufficient, need actual results)
                                            if check_all_results_exist(root_dir, work_dir, sample_name):
                                                # All result files exist (FF, WC, WCX)
                                                process_status[sample_name] = "Completed"
                                                containers_to_remove.append(container_name)
                                            else:
                                                # Container exited but result files are missing - need to re-run
                                                logger.warning(f"Container {container_name} exited with code 0 but result files are missing for {sample_name}, will re-run")
                                                process_status[sample_name] = "Init"  # Reset to Init so it can be re-run
                                                containers_to_remove.append(container_name)
                                else:
                                    # Container exited with error
                                    process_status[sample_name] = "Failed"
                                    containers_to_remove.append(container_name)
                            else:
                                # Could not check exit code, mark as failed
                                process_status[sample_name] = "Failed"
                                containers_to_remove.append(container_name)
                        else:
                            process_status[sample_name] = current_status

        # 완료된 컨테이너들 정리
        for container_name in containers_to_remove:
            logger.info(f"{container_name} will be removed")
            try:
                # First, check if container is running and stop it if needed
                check_result = subprocess.run(["docker", "inspect", "--format", "{{.State.Running}}", container_name],
                                             capture_output=True, text=True, timeout=5)
                if check_result.returncode == 0:
                    is_running = check_result.stdout.strip() == "true"
                    if is_running:
                        logger.info(f"Stopping running container: {container_name}")
                        stop_result = subprocess.run(["docker", "stop", container_name],
                                                   capture_output=True, text=True, timeout=10)
                        if stop_result.returncode != 0:
                            logger.warning(f"Failed to stop container {container_name}: {stop_result.stderr}")
                
                # Remove container (use -f to force remove if still running)
                remove_result = subprocess.run(["docker", "rm", "-f", container_name],
                                             capture_output=True, text=True, timeout=10)
                if remove_result.returncode == 0:
                    logger.info(f"Successfully removed container: {container_name}")
                else:
                    logger.warning(f"Failed to remove container {container_name}: {remove_result.stderr}")
            except subprocess.TimeoutExpired:
                logger.warning(f"Timeout removing container {container_name}")
            except Exception as e:
                logger.error(f"Error removing container {container_name}: {e}")

    except subprocess.CalledProcessError as e:
        logger.error(f"Error running docker command: {e}")
    except Exception as e:
        logger.error(f"Error updating process status: {e}", exc_info=True)

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
def start_sample(samples_dic, cmd, labcode, max_samples, process_status, root_dir, force=False):
    """새 샘플 시작
    
    Args:
        force: If True, ignore existing result files when updating status
    """
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
            # Extract sample info from dictionary
            sample_info = samples_dic[next_sample]
            sample_id = sample_info['sample_id']
            work_dir = sample_info['work_dir']
            gender = sample_info.get('gender', '')
            
            # Build command
            cmd_to_run = cmd.replace('SAMPLE_NAME', sample_id)
            cmd_to_run = cmd_to_run.replace('WORK_DIR', work_dir)
            cmd_to_run = cmd_to_run.replace('LABCODE', labcode)
            
            # Add gender if available (remove placeholder if not)
            if gender:
                cmd_to_run = cmd_to_run.replace('--fetal_gender GENDER', f'--fetal_gender {gender}')
            else:
                # Remove the gender option entirely
                cmd_to_run = cmd_to_run.replace(' --fetal_gender GENDER', '').strip()

            logger.warning(f"Running : {cmd_to_run}")

            try:
                # Execute the command in background (non-blocking)
                # Since run_md_pipeline.sh uses --detached, it should return quickly after starting the container
                process = subprocess.Popen(
                    cmd_to_run.split(' '),
                    shell=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )
                
                # Give it a moment to start the container
                time.sleep(2)
                
                # Check if process is still running (should have started container and exited)
                # For detached mode, the script should return quickly after starting the container
                if process.poll() is None:
                    # Process still running, wait a bit more (max 10 seconds)
                    try:
                        stdout, stderr = process.communicate(timeout=10)
                        return_code = process.returncode
                    except subprocess.TimeoutExpired:
                        # Process taking too long, kill it
                        process.kill()
                        stdout, stderr = process.communicate()
                        return_code = process.returncode
                        logger.error(f"Command timeout for {next_sample}")
                else:
                    # Process already finished
                    stdout, stderr = process.communicate()
                    return_code = process.returncode
                
                logger.info(f"Return code: {return_code}")
                if stdout:
                    logger.info(f"STDOUT: {stdout[:500]}")  # First 500 chars
                if stderr:
                    logger.info(f"STDERR: {stderr[:500]}")  # First 500 chars

                # Check if container was started successfully
                # In detached mode, the script should return 0 after starting container
                if return_code == 0:
                    logger.info(f"Container started for {next_sample}")
                    process_status[next_sample] = "Running"  # 실행 중으로 상태 변경
                else:
                    logger.error(f"Command failed with return code {return_code}")
                    if stderr:
                        logger.error(f"Error output: {stderr[:500]}")
                    process_status[next_sample] = "Failed"  # 실패로 상태 변경

            except Exception as e:
                logger.error(f"Unexpected error starting {next_sample}: {e}", exc_info=True)
                process_status[next_sample] = "Failed"
            
            # Update status after starting (check if container actually started)
            update_process_status(process_status, samples_dic, root_dir, force=force)

        except Exception as e:
            logger.error(f"Error running : {e}")
            process_status[next_sample] = "Failed"

        return True

    return False

def generate_sample_sheet(scan_dir, output_tsv, root_dir, work_dir):
    """Generate TSV file from specified directory
    
    Scans subdirectories in the specified directory and checks for BAM files.
    """
    rows = []
    summary = {
        'total_dirs': 0,
        'with_bam': 0,
        'without_bam': 0,
        'with_json': 0,
        'without_json': 0
    }
    
    # Determine which directory to scan
    if scan_dir:
        scan_path = Path(scan_dir)
    else:
        # Default to analysis/{work_dir}/ directory
        scan_path = Path(root_dir) / "analysis" / work_dir
    
    logger.info(f"Scanning directory: {scan_path}")
    
    if not scan_path.exists():
        logger.error(f"Directory not found: {scan_path}")
        return False
    
    if not scan_path.is_dir():
        logger.error(f"Path is not a directory: {scan_path}")
        return False
    
    # Scan for sample directories
    for sample_dir in sorted(scan_path.iterdir()):
        if not sample_dir.is_dir():
            continue
        
        summary['total_dirs'] += 1
        sample_id = sample_dir.name
        
        # Check if BAM file exists (try multiple naming patterns)
        bam_file = None
        bam_patterns = [
            f"{sample_id}.proper_paired.bam",
            f"{sample_id}.sorted.bam",
            f"{sample_id}.bam"
        ]
        
        for pattern in bam_patterns:
            candidate = sample_dir / pattern
            if candidate.exists():
                bam_file = candidate
                break
        
        if not bam_file:
            summary['without_bam'] += 1
            logger.warning(f"[{sample_id}] No BAM file found (checked: {', '.join(bam_patterns)})")
            # Still add to TSV but mark BAM path as empty or relative
            bam_path = f"analysis/{work_dir}/{sample_id}/{sample_id}.proper_paired.bam"
        else:
            summary['with_bam'] += 1
            # Always use analysis/{work_dir}/{sample_id}/... format for consistency with run_md_pipeline.sh
            # This format is expected by run_md_pipeline.sh: $ROOT_DIR/analysis/$WORK_DIR/$SAMPLE_ID/...
            bam_filename = bam_file.name
            bam_path = f"analysis/{work_dir}/{sample_id}/{bam_filename}"
        
        # Try to read JSON for gender info
        json_file = sample_dir / f"{sample_id}.json"
        gender = ''
        if json_file.exists():
            summary['with_json'] += 1
            try:
                with open(json_file, 'r') as f:
                    metadata = json.load(f)
                    gender = metadata.get('gender', '')
            except Exception as e:
                logger.warning(f"[{sample_id}] Failed to read JSON: {e}")
        else:
            summary['without_json'] += 1
        
        # Add row: SAMPLE_NAME, WORK_DIR, BAM_PATH, GENDER
        rows.append([sample_id, work_dir, bam_path, gender])
    
    # Write TSV file
    try:
        output_path = Path(output_tsv)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['SAMPLE_NAME', 'WORK_DIR', 'BAM_PATH', 'GENDER'])
            writer.writerows(rows)
        
        # Print summary
        logger.info("=" * 60)
        logger.info("TSV Generation Summary")
        logger.info("=" * 60)
        logger.info(f"Scanned directory: {scan_path}")
        logger.info(f"Total subdirectories found: {summary['total_dirs']}")
        logger.info(f"  ✓ With BAM file: {summary['with_bam']}")
        logger.info(f"  ✗ Without BAM file: {summary['without_bam']}")
        logger.info(f"  ✓ With JSON file: {summary['with_json']}")
        logger.info(f"  ✗ Without JSON file: {summary['without_json']}")
        logger.info(f"Generated TSV file: {output_path}")
        logger.info(f"Total rows written: {len(rows)}")
        logger.info("=" * 60)
        
        return True
    except Exception as e:
        logger.error(f"Failed to write TSV file: {e}")
        return False

def load_samples_from_tsv(sample_list, root_dir):
    """Load samples from TSV file (like run_batch_dev.py)"""
    samples_dic = {}
    
    try:
        sample_df = pd.read_csv(sample_list, sep='\t', dtype={"SAMPLE_NAME": str})
        
        # Support different TSV formats
        if 'WORK_DIR' in sample_df.columns:
            if 'GENDER' in sample_df.columns:
                # Format: SAMPLE_NAME, WORK_DIR, BAM_PATH, GENDER
                samples_dic = {
                    row["SAMPLE_NAME"]: {
                        'sample_id': row["SAMPLE_NAME"],
                        'work_dir': row['WORK_DIR'],
                        'gender': row.get('GENDER', '')
                    }
                    for _, row in sample_df.iterrows()
                }
            else:
                # Format: SAMPLE_NAME, WORK_DIR, BAM_PATH
                samples_dic = {
                    row["SAMPLE_NAME"]: {
                        'sample_id': row["SAMPLE_NAME"],
                        'work_dir': row['WORK_DIR'],
                        'gender': ''
                    }
                    for _, row in sample_df.iterrows()
                }
        else:
            logger.error("TSV file must contain SAMPLE_NAME and WORK_DIR columns")
            return samples_dic
        
        logger.info(f"Loaded {len(samples_dic)} samples from {sample_list}")
        
    except Exception as e:
        logger.error(f"Error reading sample list: {e}")
        return samples_dic
    
    return samples_dic

def scan_batch_output(batch_output_dir, root_dir, work_dir):
    """Scan batch_output directory and create sample dictionary"""
    samples_dic = {}
    batch_path = Path(batch_output_dir)
    
    if not batch_path.exists():
        logger.error(f"Batch output directory not found: {batch_output_dir}")
        return samples_dic
    
    # Scan for sample directories
    for sample_dir in batch_path.iterdir():
        if not sample_dir.is_dir():
            continue
        
        sample_id = sample_dir.name
        
        # Check if BAM file exists
        bam_file = sample_dir / f"{sample_id}.proper_paired.bam"
        if not bam_file.exists():
            logger.warning(f"BAM file not found for {sample_id}: {bam_file}")
            continue
        
        # Try to read JSON for gender info
        json_file = sample_dir / f"{sample_id}.json"
        gender = ''
        if json_file.exists():
            try:
                with open(json_file, 'r') as f:
                    metadata = json.load(f)
                    gender = metadata.get('gender', '')
            except Exception as e:
                logger.warning(f"Failed to read JSON for {sample_id}: {e}")
        
        # Move BAM file to analysis directory if needed
        analysis_dir = Path(root_dir) / "analysis" / work_dir / sample_id
        analysis_dir.mkdir(parents=True, exist_ok=True)
        
        target_bam = analysis_dir / f"{sample_id}.proper_paired.bam"
        
        # Copy BAM if not already in analysis directory
        if not target_bam.exists():
            logger.info(f"Copying BAM file for {sample_id}...")
            try:
                import shutil
                shutil.copy2(bam_file, target_bam)
                # Copy index if exists
                bam_index = bam_file.with_suffix('.bam.bai')
                if bam_index.exists():
                    target_index = target_bam.with_suffix('.bam.bai')
                    shutil.copy2(bam_index, target_index)
            except Exception as e:
                logger.error(f"Failed to copy BAM for {sample_id}: {e}")
                continue
        
        # Add to samples dictionary
        samples_dic[sample_id] = {
            'sample_id': sample_id,
            'work_dir': work_dir,
            'gender': gender
        }
    
    logger.info(f"Found {len(samples_dic)} samples in {batch_output_dir}")
    return samples_dic

# Main function
def run_daemon(cmd, labcode, sample_list, batch_output_dir, root_dir, work_dir, max_samples, intv, force=False):
    """Run MD pipeline daemon for artificial samples
    
    Args:
        force: If True, ignore existing result files and re-run all samples
    """
    
    # Load samples from TSV file or scan directory
    if sample_list:
        # Load from TSV file
        samples_dic = load_samples_from_tsv(sample_list, root_dir)
    elif batch_output_dir:
        # Scan batch output directory
        samples_dic = scan_batch_output(batch_output_dir, root_dir, work_dir)
    else:
        logger.error("Either --sample_list or --batch_output must be provided")
        return
    
    if not samples_dic:
        logger.error("No samples found. Exiting.")
        return
    
    total_samples = len(samples_dic)
    process_status = {sample: "Init" for sample in samples_dic}
    
    # Check for already completed samples (result files exist)
    # Skip this check if force=True
    if force:
        logger.info("Force mode: Ignoring existing result files, will re-run all samples")
    else:
        logger.info("Checking for already completed samples...")
        for sample_name in samples_dic:
            sample_info = samples_dic[sample_name]
            work_dir = sample_info['work_dir']
            
            if check_all_results_exist(root_dir, work_dir, sample_name):
                logger.info(f"Sample {sample_name} already completed (all result files exist)")
                process_status[sample_name] = "Completed"
    
    print_process_status(process_status)

    # Continuously check if additional samples need to be started
    while True:
        try:
            # 새 샘플 시작 시도 (사용 가능한 슬롯만큼 반복)
            started_any = False
            while True:
                # Try to start one sample
                if start_sample(samples_dic, cmd, labcode, max_samples, process_status, root_dir, force=force):
                    started_any = True
                    time.sleep(1)  # Small delay between starts
                else:
                    # No more slots available or no more samples to start
                    break
            
            if started_any:
                logger.info("Finished starting available samples")
            
            # 프로세스 상태 업데이트
            update_process_status(process_status, samples_dic, root_dir, force=force)
            
            # 상태 출력 (30초마다)
            logger.info("=" * 60)
            logger.info(f"Status update at {time.strftime('%Y-%m-%d %H:%M:%S')}")
            print_process_status(process_status)
            
            # 모든 샘플이 완료되었는지 확인
            if all_finished(process_status, total_samples):
                logger.info("All samples have been processed!")
                break
            
            logger.info(f"Waiting {intv} seconds before next check...")
            time.sleep(intv)  # 지정된 간격으로 대기
        except KeyboardInterrupt:
            logger.info("Interrupted by user. Exiting...")
            break
        except Exception as e:
            logger.error(f"Error in main loop: {e}", exc_info=True)
            logger.info("Continuing after error...")
            time.sleep(intv)  # 에러 후에도 대기

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

    parser = argparse.ArgumentParser(
        description="MD Pipeline Batch Runner for Artificial Samples",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Using TSV file (like run_batch_dev.py)
  python run_batch_md.py -l cordlife -s sample_sheet.tsv -r /home/ken/ken-nipt -w md_test -m 1 -i 30
  
  # Scanning batch_output directory
  python run_batch_md.py -l cordlife -b batch_output -r /home/ken/ken-nipt -w md_test -m 1 -i 30
  
  # Generate TSV file from specific directory (e.g., batch_output)
  python run_batch_md.py --generate-tsv sample_sheet.tsv /home/ken/ken-nipt/src/md_test/batch_output
  
  # Generate TSV file from analysis/{work_dir}/ directory (scan existing BAM files)
  python run_batch_md.py --generate-tsv sample_sheet.tsv -r /home/ken/ken-nipt -w md_test
        """
    )

    # Input source: either TSV file or batch_output directory (not required if --generate-tsv is used)
    input_group = parser.add_mutually_exclusive_group(required=False)
    input_group.add_argument('-s', '--sample_list', help='Path to sample list TSV file (format: SAMPLE_NAME, WORK_DIR, BAM_PATH, GENDER)')
    input_group.add_argument('-b', '--batch_output', help='Path to batch_output directory containing artificial samples')
    
    parser.add_argument('-l', '--labcode', help='Lab code (required for daemon mode, optional for --generate-tsv)')
    parser.add_argument('-r', '--root_dir', help='Root directory (e.g., /home/ken/ken-nipt). Required for daemon mode, optional for --generate-tsv (defaults to current directory)')
    parser.add_argument('-w', '--work_dir', help='Work directory (e.g., md_test). Required for daemon mode, optional for --generate-tsv (used only if scan_dir not provided)')
    parser.add_argument('-m', '--max_samples', type=int, default=2, help='Maximum concurrent samples (default: 2)')
    parser.add_argument('-i', '--poll_intv', type=int, default=30, help='Polling interval in seconds (default: 30)')
    
    parser.add_argument('--generate-tsv', nargs='+', metavar=('OUTPUT_TSV', '[SCAN_DIR]'), 
                       help='Generate TSV file. Scans subdirectories in SCAN_DIR (or analysis/{work_dir}/ if not provided). Checks for BAM files and shows summary.')

    parser.add_argument('--ignore-zscore', action='store_true',
                       help='Ignore z-score threshold requirement (only check overlap with target)')
    parser.add_argument('--ignore-min-length', action='store_true',
                       help='Ignore minimum length requirement (only check z-score threshold)')
    parser.add_argument('--skip-npz', action='store_true',
                       help='Skip NPZ file creation for WC/WCX (use existing NPZ files, faster execution)')
    parser.add_argument('-f', '--force', action='store_true',
                       help='Force re-run even if output files already exist')

    parser.add_argument("--log_stdout", default="Y", choices=["Y", "N"], help="Pipeline logging to Standard output (default: Y)")

    args = parser.parse_args()

    # Setup logging
    if args.log_stdout == "N":
        log_file = Path(args.batch_output) / "batch_md.log"
        logging.basicConfig(filename=str(log_file), level=logging.INFO, format="%(asctime)s [%(levelname) 7s] (%(lineno)4d) | %(message)s")
    else:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname) 7s] (%(lineno)4d) | %(message)s")

    logging.addLevelName(logging.ERROR, "\033[91m%7s\033[0m" % logging.getLevelName(logging.ERROR))
    logging.addLevelName(logging.WARNING, "\033[93m%7s\033[0m" % logging.getLevelName(logging.WARNING))
    logger = logging.getLogger(__name__)

    try:
        # Generate TSV file mode
        if args.generate_tsv:
            if len(args.generate_tsv) == 1:
                # Only output_tsv provided
                output_tsv = args.generate_tsv[0]
                scan_dir = None
                # Use work_dir if provided, otherwise default to analysis/
                if not args.work_dir:
                    args.work_dir = "md_test"  # Default work_dir
                    logger.warning(f"Work directory not specified, using default: {args.work_dir}")
                if not args.root_dir:
                    args.root_dir = os.getcwd()
                    logger.warning(f"Root directory not specified, using current directory: {args.root_dir}")
                logger.info(f"Generating TSV file from analysis/{args.work_dir}/ directory...")
            elif len(args.generate_tsv) == 2:
                # Both output_tsv and scan_dir provided
                output_tsv, scan_dir = args.generate_tsv
                logger.info(f"Generating TSV file from directory: {scan_dir}")
                
                # If output_tsv is a relative path, save it in scan_dir
                output_tsv_path = Path(output_tsv)
                scan_path = Path(scan_dir)
                if not output_tsv_path.is_absolute():
                    # Relative path: save in scan_dir
                    output_tsv = str(scan_path / output_tsv)
                    logger.info(f"Output TSV will be saved to: {output_tsv}")
                
                # Auto-detect work_dir and root_dir from scan_dir if not provided
                
                if not args.work_dir:
                    # Use last directory name as work_dir
                    args.work_dir = scan_path.name
                    logger.info(f"Auto-detected WORK_DIR from scan directory: {args.work_dir}")
                
                if not args.root_dir:
                    # Try to find the root directory (usually /home/ken/ken-nipt)
                    # If scan_dir is like /home/ken/ken-nipt/analysis/md_validation/1p36
                    # Then root_dir should be /home/ken/ken-nipt
                    # If scan_dir is like /home/ken/ken-nipt/analysis/md_validation/1p36
                    # Try to go up until we find a directory that looks like root
                    current = scan_path.parent
                    # Look for common root patterns
                    while current != current.parent:  # Stop at root
                        # Check if this looks like a root directory (contains 'analysis' subdirectory)
                        if (current / "analysis").exists():
                            args.root_dir = str(current)
                            logger.info(f"Auto-detected root_dir from scan directory: {args.root_dir}")
                            break
                        current = current.parent
                    else:
                        # If not found, use parent of scan_dir (original behavior)
                        args.root_dir = str(scan_path.parent)
                        logger.info(f"Auto-detected root_dir from scan directory (fallback): {args.root_dir}")
                
                # If root_dir was auto-detected and contains 'analysis', 
                # we need to adjust work_dir to include the path from analysis/ to scan_dir
                root_path = Path(args.root_dir)
                if (root_path / "analysis").exists():
                    # Calculate relative path from root_dir/analysis to scan_dir
                    try:
                        analysis_rel = scan_path.relative_to(root_path / "analysis")
                        # If the relative path has multiple parts, use all parts as work_dir
                        # Otherwise, just use the last part
                        if len(analysis_rel.parts) > 1:
                            args.work_dir = str(analysis_rel)
                        else:
                            args.work_dir = analysis_rel.name
                        logger.info(f"Adjusted WORK_DIR to match analysis structure: {args.work_dir}")
                    except ValueError:
                        # scan_dir is not under root_dir/analysis, keep original work_dir
                        pass
            else:
                logger.error("--generate-tsv expects 1 or 2 arguments: OUTPUT_TSV [SCAN_DIR]")
                sys.exit(1)
            
            success = generate_sample_sheet(scan_dir, output_tsv, args.root_dir, args.work_dir)
            sys.exit(0 if success else 1)
        
        # Daemon mode: check required arguments
        if not args.labcode:
            parser.error("--labcode (-l) is required for daemon mode")
        if not args.root_dir:
            parser.error("--root_dir (-r) is required for daemon mode")
        if not args.work_dir:
            parser.error("--work_dir (-w) is required for daemon mode")
        if not args.sample_list and not args.batch_output:
            parser.error("Either --sample_list (-s) or --batch_output (-b) must be provided for daemon mode")
        
        logger.info("Run MD batch daemon...")

        # Build command template
        # Note: --fetal_gender GENDER will be replaced or removed based on availability
        command = f"bash {args.root_dir}/src/run_md_pipeline.sh -s SAMPLE_NAME -l LABCODE -root {args.root_dir} -work WORK_DIR --fetal_gender GENDER --detached"
        
        # Add ignore options if specified
        if args.ignore_zscore:
            command += " --ignore-zscore"
        if args.ignore_min_length:
            command += " --ignore-min-length"
        if args.skip_npz:
            command += " --skip-npz"
        if args.force:
            command += " -f"

        logger.info(f"Starting MD batch daemon with:")
        logger.info(f"  Lab code: {args.labcode}")
        if args.sample_list:
            logger.info(f"  Sample list (TSV): {args.sample_list}")
        if args.batch_output:
            logger.info(f"  Batch output: {args.batch_output}")
        logger.info(f"  Root dir: {args.root_dir}")
        logger.info(f"  Work dir: {args.work_dir}")
        logger.info(f"  Max samples: {args.max_samples}")
        logger.info(f"  Poll interval: {args.poll_intv}s")
        logger.info(f"  Ignore z-score: {args.ignore_zscore}")
        logger.info(f"  Ignore min-length: {args.ignore_min_length}")
        logger.info(f"  Skip NPZ: {args.skip_npz}")
        logger.info(f"  Force execution: {args.force}")

        run_daemon(command, args.labcode, args.sample_list, args.batch_output, args.root_dir, args.work_dir, args.max_samples, args.poll_intv, force=args.force)

    except Exception as e:
        logger.exception("Error occurred during pipeline execution.")

