import pandas as pd
import sys
import json
import logging
import os
from pathlib import Path

# Setup logger - use root logger to ensure logging works across modules
# This ensures that logs from this module are visible when imported
logger = logging.getLogger(__name__)
# If logger has no handlers or level is NOTSET, use root logger
if not logger.handlers or logger.level == logging.NOTSET:
    logger = logging.getLogger()
    # Ensure root logger has at least INFO level
    if logger.level == logging.NOTSET:
        logger.setLevel(logging.INFO)
    # Add a stream handler if none exists
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(handler)

# Function to check overlap between two regions
def check_overlap(start1, end1, start2, end2):
    return max(start1, start2) <= min(end1, end2)

def check_full_inclusion(start1, end1, start2, end2):
    return start1 <= start2 and end1 >= end2

def process_wc(input_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if "# Test results" in line:
                data_start = i + 1
                break

    df = pd.read_csv(input_file, sep=r"\s+", skiprows=data_start)
    if len(df) > 0:
        df[['chr', 'start_end']] = df['location'].str.split(':', expand=True)
        df[['start', 'end']] = df['start_end'].str.split('-', expand=True)
        df['start'] = pd.to_numeric(df['start'])
        df['end'] = pd.to_numeric(df['end'])
        df['length'] = df['end'] - df['start']
        df = df[['chr', 'start', 'end', 'effect', 'z-score', 'length']]
        df.columns = ['chr', 'start', 'end', 'effect', 'zscore', 'length']

    return df

def process_wcx(input_file):
    df = pd.read_csv(input_file, sep=r"\s+")
    df['chr'] = df['chr'].astype('str')
    df['length'] = df['end'] - df['start'] + 1
    df['effect'] = df['ratio']
    df = df[['chr', 'start', 'end', 'effect', 'zscore', 'length']]

    return df

def process_wcff(input_file):
    df = pd.read_csv(input_file, sep=r"\s+")
    df['chr'] = df['chr'].astype('str')
    df['length'] = df['end'] - df['start'] + 1
    df['effect'] = df['ratio']
    df = df[['chr', 'start', 'end', 'effect', 'zscore', 'length']]

    return df

def get_min_length(length, md_cfg):
    """Calculate minimum length requirement based on target length and config"""
    if length <= md_cfg["region1"]["less_than"]:
        return md_cfg["region1"]["min_length"]
    elif length <= md_cfg["region2"]["less_than"]:
        return md_cfg["region2"]["min_length"]
    elif length <= md_cfg["region3"]["less_than"]:
        return md_cfg["region3"]["min_length"]
    elif length > md_cfg["region4"]["greater_than"]:
        return md_cfg["region4"]["min_length"]
    else:
        return md_cfg["region3"]["min_length"]


def compare_with_target(result_df, target_file, md_cfg, threshold, ignore_min_length=False, ignore_zscore=False):
    """
    Compare detection results with target database.
    
    Args:
        result_df: DataFrame with detection results
        target_file: Path to target database BED file
        md_cfg: MD configuration dictionary
        threshold: Z-score threshold
        ignore_min_length: If True, skip minimum length check (only check z-score)
        ignore_zscore: If True, skip z-score threshold check (only check overlap with target)
    
    Returns:
        DataFrame with matched results
    """
    print(f"[{__name__}] [compare_with_target] Starting comparison: {len(result_df)} result rows, ignore_min_length={ignore_min_length}, ignore_zscore={ignore_zscore}")
    logger.info(f"[compare_with_target] Starting comparison: {len(result_df)} result rows, ignore_min_length={ignore_min_length}, ignore_zscore={ignore_zscore}")
    
    # Load the DB.bed file
    db_df = pd.read_csv(target_file, sep="\t", header=None, names=['chr', 'start', 'end', 'Disease-name', 'gain_loss', 'overlap_inside', 'minimum_length', 'cri_genes', 'cri_start', 'cri_end'])
    logger.info(f"[compare_with_target] Loaded {len(db_df)} target regions from {target_file}")

    md_output_df = pd.DataFrame()

    for idx, result_row in result_df.iterrows():
        chr1 = result_row["chr"]
        start1 = int(result_row["start"])
        end1 = int(result_row["end"])
        z = float(result_row["zscore"])
        detected_length = end1 - start1 + 1

        logger.debug(f"[compare_with_target] Processing result row {idx}: chr={chr1}, start={start1}, end={end1}, z={z}, length={detected_length}")

        match_found = False

        for db_idx, db_row in db_df.iterrows():
            chr2 = str(db_row["chr"])
            start2 = int(db_row["start"])
            end2 = int(db_row["end"])
            min_len_field = db_row["minimum_length"]
            mode = db_row["overlap_inside"]
            gain_loss = db_row["gain_loss"]
            target_length = end2 - start2 + 1

            overlap = False
            if mode == "overlap":
                detected = check_overlap(start1, end1, start2, end2)
            elif mode == "inside":
                detected = check_full_inclusion(start1, end1, start2, end2)
            else:
                continue  # skip invalid rows

            if chr1 == chr2 and detected:
                print(f"[{__name__}] [OVERLAP FOUND] chr={chr1}, target={db_row['Disease-name']}, detected_region={start1}-{end1}, target_region={start2}-{end2}, z={z}, threshold={threshold}, gain_loss={gain_loss}")
                logger.info(f"[OVERLAP FOUND] chr={chr1}, target={db_row['Disease-name']}, detected_region={start1}-{end1}, target_region={start2}-{end2}, z={z}, threshold={threshold}, gain_loss={gain_loss}")
                
                # Check zscore match with gain/loss type
                zscore_match = False
                if ignore_zscore:
                    # Skip z-score check, any overlap is considered a match
                    zscore_match = True
                    print(f"[{__name__}] [Z-SCORE CHECK SKIPPED] Target={db_row['Disease-name']}, z={z}, ignore_zscore=True")
                    logger.info(f"[Z-SCORE CHECK SKIPPED] Target={db_row['Disease-name']}, z={z}, ignore_zscore=True")
                else:
                    # Normal z-score check
                    if gain_loss == 'loss' and z < -float(threshold):
                        zscore_match = True
                    elif gain_loss == 'gain' and z > float(threshold):
                        zscore_match = True
                    elif gain_loss == 'both' and (z < -float(threshold) or z > float(threshold)):
                        zscore_match = True
                
                if not zscore_match:
                    # Threshold condition check details
                    if gain_loss == 'loss':
                        required_condition = f"z < -{threshold} (i.e., z < {-float(threshold)})"
                        actual_condition = f"{z} < {-float(threshold)} = {z < -float(threshold)}"
                    elif gain_loss == 'gain':
                        required_condition = f"z > {threshold}"
                        actual_condition = f"{z} > {float(threshold)} = {z > float(threshold)}"
                    else:  # both
                        required_condition = f"z < -{threshold} or z > {threshold}"
                        actual_condition = f"{z} < {-float(threshold)} or {z} > {float(threshold)} = {z < -float(threshold) or z > float(threshold)}"
                    
                    print(f"[{__name__}] [Z-SCORE MISMATCH] Target={db_row['Disease-name']}, z={z}, threshold={threshold}, gain_loss={gain_loss}, required={required_condition}, actual={actual_condition}")
                    logger.info(f"[Z-SCORE MISMATCH] Target={db_row['Disease-name']}, z={z}, threshold={threshold}, gain_loss={gain_loss}, required={required_condition}, actual={actual_condition}")
                    continue

                logger.info(f"[compare_with_target] Z-score match: z={z}, threshold={threshold}, gain_loss={gain_loss}")

                # Check cri_genes condition
                cri_genes = db_row.get("cri_genes", "-")
                if cri_genes != "-":
                    try:
                        cri_start = int(db_row["cri_start"])
                        cri_end = int(db_row["cri_end"])
                    except ValueError:
                        logger.debug(f"[compare_with_target] Invalid cri_genes range, skipping")
                        continue  # skip if invalid cri range

                    # Make sure the entire region is within critical range
                    detected = check_overlap(start1, end1, cri_start, cri_end)
                    if detected == True:
                        logger.info(f"[MATCHED] Critical genes match: {db_row['Disease-name']}, detected_length={detected_length}")
                        result_row["Disease-name"] = db_row["Disease-name"]
                        md_output_df = pd.concat([md_output_df, result_row.to_frame().T], ignore_index=True)
                        match_found = True
                        break  # match found
                    else:
                        logger.debug(f"[compare_with_target] Critical genes mismatch")
                        continue  # doesn't fit cri_gene bounds

                else:
                    # Determine minimum length requirement
                    if ignore_min_length or ignore_zscore:
                        # Skip length check if either ignore_min_length or ignore_zscore is True
                        # When ignore_zscore=True, we want to include any overlap regardless of length
                        logger.info(f"[MATCHED] Ignoring min_length: {db_row['Disease-name']}, detected_length={detected_length}, z={z}, ignore_min_length={ignore_min_length}, ignore_zscore={ignore_zscore}")
                        print(f"[{__name__}] [MATCHED] Ignoring min_length: {db_row['Disease-name']}, detected_length={detected_length}, z={z}, ignore_min_length={ignore_min_length}, ignore_zscore={ignore_zscore}")
                        result_row["Disease-name"] = db_row["Disease-name"]
                        md_output_df = pd.concat([md_output_df, result_row.to_frame().T], ignore_index=True)
                        match_found = True
                        break  # One good match is enough
                    else:
                        # Apply minimum length check
                        if min_len_field != "-":
                            try:
                                min_len = int(min_len_field)
                            except ValueError:
                                logger.debug(f"[compare_with_target] Invalid min_len_field: {min_len_field}")
                                continue
                        else:
                            min_len = get_min_length(target_length, md_cfg)
                            logger.debug(f"[compare_with_target] Calculated min_len={min_len} for target_length={target_length}")

                        if detected_length >= min_len:
                            logger.info(f"[MATCHED] detected_length={detected_length}, min_len={min_len}, disease={db_row['Disease-name']}")
                            result_row["Disease-name"] = db_row["Disease-name"]
                            md_output_df = pd.concat([md_output_df, result_row.to_frame().T], ignore_index=True)
                            match_found = True
                            break  # One good match is enough
                        else:
                            logger.info(f"[SHORT LENGTH] detected_length={detected_length}, min_len={min_len}, disease={db_row['Disease-name']}")

    logger.info(f"[compare_with_target] Found {len(md_output_df)} matches out of {len(result_df)} result rows")
    return md_output_df


def process_microdeletion_result(sample_name, tool_name, input_file, tag, target_file, config_file, md_select, output_dir, threshold, ignore_min_length=False, ignore_zscore=False):
    """
    Process microdeletion detection results.
    
    Args:
        sample_name: Sample name
        tool_name: Tool name (WC, WCX, WCFF)
        input_file: Input file path
        tag: Tag (orig, fetus, mom)
        target_file: Target BED file path
        config_file: Config file path
        md_select: MD target selection
        output_dir: Output directory
        threshold: Z-score threshold
        ignore_min_length: If True, skip minimum length check (only check z-score)
        ignore_zscore: If True, skip z-score threshold check (any overlap with target will be included)
    """
    print(f"[{__name__}] [process_microdeletion_result] Starting: sample={sample_name}, tool={tool_name}, tag={tag}, ignore_min_length={ignore_min_length}, ignore_zscore={ignore_zscore}")
    logger.info(f"[process_microdeletion_result] Starting: sample={sample_name}, tool={tool_name}, tag={tag}, ignore_min_length={ignore_min_length}, ignore_zscore={ignore_zscore}")
    print(f"[{__name__}] [process_microdeletion_result] Input file: {input_file}")
    logger.info(f"[process_microdeletion_result] Input file: {input_file}")
    print(f"[{__name__}] [process_microdeletion_result] Target file: {target_file}")
    logger.info(f"[process_microdeletion_result] Target file: {target_file}")
    print(f"[{__name__}] [process_microdeletion_result] Config file: {config_file}")
    logger.info(f"[process_microdeletion_result] Config file: {config_file}")
    
    if tool_name == "WC":
        result_df = process_wc(input_file)
    elif tool_name == "WCX":
        result_df = process_wcx(input_file)
    elif tool_name == "WCFF":
        result_df = process_wcff(input_file)
    else:
        logger.info(f"[ERROR] Invalid tool name: {tool_name}. Use 'WC' or 'WCX' or 'WCFF'.")
        return False

    logger.info(f"[process_microdeletion_result] Loaded {len(result_df)} result rows from {input_file}")
    
    if len(result_df) > 0:
        logger.info(f"[process_microdeletion_result] First result row: {result_df.iloc[0].to_dict()}")

    md_file_dict = {
        "MD_Target_8": "md8",
        "MD_Target_108": "md108",
        "MD_Target_320": "md320",
        "MD_Target_87": "md87"
    }

    if md_select not in md_file_dict:
        logger.info(f"[ERROR] Invalid md_select: {md_select}.")
        return False

    output_filename_md = f"{sample_name}_{tool_name}_{tag}_{md_file_dict[md_select]}.tsv"

    try:
        with open(config_file) as f:
            cfg_json = json.load(f)
    except Exception as e:
        logger.info(f"[ERROR] Failed to load config file: {e}")
        return False

    if md_select not in cfg_json:
        logger.info(f"[ERROR] {md_select} not found in config file.")
        return False

    md_cfg = cfg_json[md_select]

    if result_df.empty:
        logger.info(f"[INFO] No result found for {tool_name} in {md_select} : {sample_name}.")
        logger.info(f"[INFO] This means the input file exists but contains no detection results.")
        return False

    logger.info(f"[process_microdeletion_result] Calling compare_with_target with {len(result_df)} result rows")
    logger.info(f"[process_microdeletion_result] Target file: {target_file}, Threshold: {threshold}, IgnoreMinLength: {ignore_min_length}, IgnoreZscore: {ignore_zscore}")
    
    output_df = compare_with_target(result_df, target_file, md_cfg, threshold, ignore_min_length=ignore_min_length, ignore_zscore=ignore_zscore)

    if output_df.empty:
        logger.info(f"[INFO] No matching microdeletion targets found for {tool_name} in {md_select}:{sample_name}.")
        logger.info(f"[INFO] This means detection results exist but did not match any target regions in {target_file}")
        return False

    output_path = os.path.join(output_dir, tag, output_filename_md)
    output_df.to_csv(output_path, sep='\t', index=False)
    logger.info(f"[SUCCESS] {output_path} was generated successfully!")

    return True


def run_microdeletion_test_decision_pipeline(
    sample_name, labcode, config, analysis_dir, output_dir, bed_dir,
    types=None, md_targets=None, ignore_min_length=False, ignore_zscore=False
):
    """
    Run microdeletion decision for specified methods, types, and MD targets.
    This is a flexible version for testing/artificial samples.
    
    Args:
        sample_name: Sample ID
        labcode: Lab code
        config: Configuration dictionary
        analysis_dir: Analysis directory path
        output_dir: Output directory path
        bed_dir: BED files directory path
        types: List of BAM types to process (default: ["orig", "fetus"])
        md_targets: List of MD targets to process (default: ["MD_Target_8"])
        ignore_min_length: If True, skip minimum length check (only check z-score)
        ignore_zscore: If True, skip z-score threshold check (only check overlap with target)
    
    Returns:
        bool: True if successful
    """
    
    # Print to stdout immediately (for debugging)
    print(f"[{__name__}] [run_microdeletion_test_decision_pipeline] Starting: sample={sample_name}, ignore_min_length={ignore_min_length}, ignore_zscore={ignore_zscore}")
    logger.info(f"[run_microdeletion_test_decision_pipeline] Starting: sample={sample_name}, ignore_min_length={ignore_min_length}, ignore_zscore={ignore_zscore}")
    print(f"[{__name__}] [run_microdeletion_test_decision_pipeline] Analysis dir: {analysis_dir}")
    logger.info(f"[run_microdeletion_test_decision_pipeline] Analysis dir: {analysis_dir}")
    print(f"[{__name__}] [run_microdeletion_test_decision_pipeline] Output dir: {output_dir}")
    logger.info(f"[run_microdeletion_test_decision_pipeline] Output dir: {output_dir}")
    print(f"[{__name__}] [run_microdeletion_test_decision_pipeline] Bed dir: {bed_dir}")
    logger.info(f"[run_microdeletion_test_decision_pipeline] Bed dir: {bed_dir}")
    
    # Default values for testing
    if types is None:
        types = ["orig", "fetus"]
    
    if md_targets is None:
        md_targets = ["MD_Target_8"]
    
    logger.info(f"[run_microdeletion_test_decision_pipeline] Types: {types}, MD targets: {md_targets}")
    
    methods = ["WC", "WCX"]

    wc_input_files = {
        "orig": f"{analysis_dir}/{sample_name}/Output_WC/orig/{sample_name}.wc.orig.report.txt",
        "fetus": f"{analysis_dir}/{sample_name}/Output_WC/fetus/{sample_name}.wc.fetus.report.txt",
        "mom": f"{analysis_dir}/{sample_name}/Output_WC/mom/{sample_name}.wc.mom.report.txt",
    }

    wcx_input_files = {
        "orig": f"{analysis_dir}/{sample_name}/Output_WCX/orig/{sample_name}.wcx.orig_aberrations.bed",
        "fetus": f"{analysis_dir}/{sample_name}/Output_WCX/fetus/{sample_name}.wcx.fetus_aberrations.bed",
        "mom": f"{analysis_dir}/{sample_name}/Output_WCX/mom/{sample_name}.wcx.mom_aberrations.bed",
    }

    config_file = f"/Work/NIPT/config/{labcode}/pipeline_config.json"
    
    logger.info(f"[run_microdeletion_test_decision_pipeline] Config file: {config_file}")
    
    # Load config file to get MD target and threshold information
    try:
        with open(config_file) as f:
            config_data = json.load(f)
        logger.info(f"[run_microdeletion_test_decision_pipeline] Config loaded successfully")
    except Exception as e:
        logger.error(f"[ERROR] Failed to load config file {config_file}: {e}")
        return False

    logger.info(f"[run_microdeletion_test_decision_pipeline] Processing {len(md_targets)} MD targets, {len(types)} types, {len(methods)} methods")
    
    for md_key in md_targets:
        if md_key not in config_data:
            logger.warning(f"[SKIP] MD target {md_key} not found in config")
            continue
        
        bed_file = f"{bed_dir}/common/{config_data[md_key]['bed']}"
        logger.info(f"[run_microdeletion_test_decision_pipeline] Processing MD target: {md_key}, BED file: {bed_file}")

        for method in methods:
            for tag in types:
                logger.info(f"[run_microdeletion_test_decision_pipeline] Processing: method={method}, tag={tag}")
                
                # Determine input file and output dir
                if method == "WC":
                    input_file = wc_input_files.get(tag)
                    output_subdir = "Output_WC"
                elif method == "WCX":
                    input_file = wcx_input_files.get(tag)
                    output_subdir = "Output_WCX"
                else:
                    continue

                if input_file is None:
                    logger.warning(f"[SKIP] Unknown type: {tag}")
                    continue

                output_path = f"{analysis_dir}/{sample_name}/{output_subdir}"
                threshold = config_data.get(method, {}).get(f"{tag}_call_threshold")

                logger.info(f"[run_microdeletion_test_decision_pipeline] Checking input file: {input_file}")
                if not Path(input_file).exists():
                    logger.warning(f"[SKIP] Input file missing: {input_file}")
                    continue

                logger.info(f"[MD] {method}-{tag} on {md_key} (Threshold: {threshold}, IgnoreMinLength: {ignore_min_length}, IgnoreZscore: {ignore_zscore})")
                logger.info(f"[MD] Processing input file: {input_file}")
                logger.info(f"[MD] Target BED file: {bed_file}")
                logger.info(f"[MD] Output directory: {output_path}")
                
                # Check if target BED file exists
                if not Path(bed_file).exists():
                    logger.error(f"[ERROR] Target BED file not found: {bed_file}")
                    continue
                
                logger.info(f"[run_microdeletion_test_decision_pipeline] Calling process_microdeletion_result...")
                success = process_microdeletion_result(
                    sample_name=sample_name,
                    tool_name=method,
                    input_file=input_file,
                    tag=tag,
                    target_file=bed_file,
                    config_file=config_file,
                    md_select=md_key,
                    output_dir=output_path,
                    threshold=threshold,
                    ignore_min_length=ignore_min_length,
                    ignore_zscore=ignore_zscore
                )

                if success:
                    logger.info(f"[PASS] MD call: {method}-{tag} {md_key}")
                else:
                    logger.info(f"[NOT MATCHED] MD call: {method}-{tag} {md_key}")
    
    logger.info(f"[run_microdeletion_test_decision_pipeline] Completed processing all targets")
    return True

