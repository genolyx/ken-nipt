#!/usr/bin/env python3
"""
Microdeletion-Only Analysis Pipeline

Performs MD analysis on existing BAM files without FASTQ processing.
Supports both proper_paired.bam and sorted.bam inputs.

Usage:
    python3 md_pipeline.py --sample_id SAMPLE001 --work_dir 2507 --labcode cordlife

Author: Ken
Version: 1.0
"""

import argparse
import datetime
import json
import logging
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam


def get_lab_bed_paths(labcode):
    """Get lab-specific BED file paths with fallback to common"""
    lab_bed_dir = Path("/Work/NIPT/data/bed") / labcode
    common_bed_dir = Path("/Work/NIPT/data/bed/common")

    def get_bed_path(filename, lab_specific=False):
        lab_path = lab_bed_dir / filename
        common_path = common_bed_dir / filename

        if lab_specific:
            # Try lab-specific first, then fallback to common
            if lab_path.exists():
                return str(lab_path)
            elif common_path.exists():
                return str(common_path)
            else:
                # Neither exists, return common path as fallback
                return str(common_path)
        else:
            # Not lab-specific, use common
            if common_path.exists():
                return str(common_path)
            else:
                return str(common_path)  # fallback

    return {
        # Main filter BED files
        "of": get_bed_path("Uniform_2017_allY.bed", lab_specific=True),
        "nf08": get_bed_path(
            "hg19_mappability_0.8_clean_all_36mer.bed", lab_specific=True
        ),
        "nf09": get_bed_path(
            "hg19_mappability_0.9_clean_all_36mer.bed", lab_specific=True
        ),
        # FF calculation BED files
        "y_regions_09": get_bed_path("chrY_target_0.9.bed"),
        "y_regions_noPARs": get_bed_path("chrY_noPARs.bed"),
        "y_regions_target": get_bed_path("chrY_target.bed"),
        "autosome_control": get_bed_path("autosome_control.bed"),
    }

sys.path.append("/Work/NIPT/bin")
sys.path.append("/Work/NIPT/bin/scripts/modules")
try:
    from process_md_result import run_microdeletion_decision_pipeline, run_microdeletion_test_decision_pipeline
except ImportError as e:
    logging.warning(f"Could not import run_microdeletion functions: {e}")
    run_microdeletion_decision_pipeline = None
    run_microdeletion_test_decision_pipeline = None

# Try to import dev version with ignore_min_length support
try:
    from process_md_result_dev import run_microdeletion_test_decision_pipeline as run_microdeletion_test_decision_pipeline_dev
    DEV_MODULE_AVAILABLE = True
except ImportError:
    DEV_MODULE_AVAILABLE = False
    run_microdeletion_test_decision_pipeline_dev = None

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

logger = logging.getLogger(__name__)


def log_and_print(message):
    """Log and print message"""
    logger.info(message)
    print(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}")


def find_bam_file(analysis_dir, work_dir, sample_id):
    """
    Find BAM file for the sample.
    Priority: proper_paired.bam > sorted.bam
    
    Returns: (bam_path, bam_type)
        bam_type: 'proper_paired' or 'sorted' or None
    """
    base_dir = os.path.join(analysis_dir, work_dir, sample_id)
    
    # Try proper_paired.bam first
    proper_paired_bam = os.path.join(base_dir, f"{sample_id}.proper_paired.bam")
    if os.path.exists(proper_paired_bam):
        log_and_print(f"Found proper_paired.bam: {proper_paired_bam}")
        return proper_paired_bam, 'proper_paired'
    
    # Try sorted.bam
    sorted_bam = os.path.join(base_dir, f"{sample_id}.sorted.bam")
    if os.path.exists(sorted_bam):
        log_and_print(f"Found sorted.bam: {sorted_bam}")
        return sorted_bam, 'sorted'
    
    log_and_print(f"ERROR: No BAM file found for {sample_id} in {base_dir}")
    return None, None


def create_proper_paired_bam(sorted_bam_path, output_path):
    """
    Create proper_paired.bam from sorted.bam
    Extracts only properly paired reads (flag 2)
    """
    log_and_print(f"Creating proper_paired.bam from {sorted_bam_path}")
    
    try:
        # samtools view -b -f 2 sorted.bam > proper_paired.bam
        cmd = [
            "samtools", "view",
            "-b", "-f", "2",
            "-o", output_path,
            sorted_bam_path
        ]
        
        log_and_print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Index the BAM
        log_and_print(f"Indexing {output_path}")
        subprocess.run(["samtools", "index", output_path], check=True)
        
        log_and_print(f"Successfully created proper_paired.bam: {output_path}")
        return True
        
    except subprocess.CalledProcessError as e:
        log_and_print(f"ERROR creating proper_paired.bam: {e}")
        log_and_print(f"STDERR: {e.stderr}")
        return False


def read_gender_from_json(sample_analysis_dir, sample_id):
    """Read gender from JSON metadata file (for artificial samples)"""
    # Try to find JSON file in sample directory
    json_pattern = os.path.join(sample_analysis_dir, "*.json")
    import glob
    json_files = glob.glob(json_pattern)
    
    if json_files:
        json_file = json_files[0]  # Use first JSON file found
        try:
            log_and_print(f"Reading gender from JSON: {json_file}")
            with open(json_file, 'r') as f:
                metadata = json.load(f)
            
            gender = metadata.get("gender")
            if gender:
                gender = str(gender).upper()
                gender = gender[0]  # Take first character (M/F)
                log_and_print(f"Gender from JSON: {gender}")
                return gender
            else:
                log_and_print("WARNING: gender not found in JSON")
                return None
                
        except Exception as e:
            log_and_print(f"ERROR reading JSON file: {e}")
            return None
    
    return None


def read_gender_from_file(sample_analysis_dir, sample_id):
    """Read gender from existing gender.txt file (gd_2 value)"""
    gender_file = os.path.join(sample_analysis_dir, "Output_FF", f"{sample_id}.gender.txt")
    
    if not os.path.exists(gender_file):
        log_and_print(f"Gender file not found: {gender_file}")
        return None
    
    try:
        log_and_print(f"Reading gender from: {gender_file}")
        with open(gender_file, 'r') as f:
            lines = f.readlines()
        
        # Parse the file to find gd_2
        for line in lines:
            if line.startswith('gd_2'):
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    gender_str = parts[2]  # XX or XY
                    if gender_str == "XY":
                        gender = "M"
                    else:
                        gender = "F"
                    log_and_print(f"Gender from gd_2: {gender_str} -> {gender}")
                    return gender
        
        log_and_print("WARNING: gd_2 not found in gender.txt")
        return None
        
    except Exception as e:
        log_and_print(f"ERROR reading gender file: {e}")
        return None


def detect_gender_from_bam(bam_file):
    """Fallback: Simple gender detection based on Y chromosome coverage"""
    try:
        log_and_print("Detecting fetal gender from BAM...")
        
        # Count Y chromosome reads
        y_cmd = ["samtools", "view", "-c", "-F", "260", bam_file, "chrY"]
        y_result = subprocess.run(y_cmd, capture_output=True, text=True, check=True)
        y_reads = int(y_result.stdout.strip())
        
        # Count autosome reads (chr1)
        a_cmd = ["samtools", "view", "-c", "-F", "260", bam_file, "chr1"]
        a_result = subprocess.run(a_cmd, capture_output=True, text=True, check=True)
        a_reads = int(a_result.stdout.strip())
        
        if a_reads == 0:
            log_and_print("WARNING: No autosome reads found, defaulting to female")
            return "F", 0.0
        
        # Calculate Y/A ratio
        y_to_a_ratio = y_reads / a_reads
        
        # Threshold for gender determination (similar to gd_1 calculation)
        threshold = 0.0001  # Adjust based on your data
        
        if y_to_a_ratio > threshold:
            gender = "M"
        else:
            gender = "F"
        
        log_and_print(f"Gender detected from BAM: {'Male' if gender == 'M' else 'Female'} (Y/A ratio: {y_to_a_ratio:.6f})")
        return gender
        
    except Exception as e:
        log_and_print(f"ERROR in gender detection: {e}, defaulting to female")
        return "F"


def process_bam_filters(sample_id, sample_dir, bed_file, filter_type="of"):
    """Create orig and fetus BAM files from proper_paired.bam
    
    Args:
        sample_id: Sample ID
        sample_dir: Sample analysis directory
        bed_file: BED file for filtering regions
        filter_type: Filter type prefix for output files (e.g., 'nf08', 'of')
        
    Returns:
        tuple: (orig_bam_path, fetus_bam_path) or (None, None) on failure
    """
    log_and_print(f"=== Processing BAM filters (orig, fetus) - Filter: {filter_type} ===")
    
    sam_tools = os.environ.get("SAMTools", "samtools")
    ref_hg = "/Work/NIPT/data/refs/common/hg19/ucsc.hg19.fasta"
    fetus_awk = "/Work/NIPT/bin/scripts/fetus.awk"
    
    # File paths with filter_type
    proper_paired_bam = os.path.join(sample_dir, f"{sample_id}.proper_paired.bam")
    orig_bam = os.path.join(sample_dir, f"{sample_id}.{filter_type}_orig.bam")
    fetus_bam = os.path.join(sample_dir, f"{sample_id}.{filter_type}_fetus.bam")
    
    # Check if files already exist
    if os.path.exists(orig_bam) and os.path.exists(fetus_bam):
        log_and_print(f"✓ BAM filters already exist, skipping creation:")
        log_and_print(f"  - {orig_bam}")
        log_and_print(f"  - {fetus_bam}")
        return orig_bam, fetus_bam
    
    # Check inputs
    if not os.path.exists(proper_paired_bam):
        log_and_print(f"ERROR: proper_paired.bam not found: {proper_paired_bam}")
        return None, None
    
    if not os.path.exists(bed_file):
        log_and_print(f"ERROR: BED file not found: {bed_file}")
        return None, None
    
    if not os.path.exists(fetus_awk):
        log_and_print(f"ERROR: fetus.awk not found: {fetus_awk}")
        return None, None
    
    try:
        # Step 1: Create orig BAM (filter by BED regions)
        if not os.path.exists(orig_bam):
            log_and_print(f"Creating orig BAM: {orig_bam}")
            cmd = [sam_tools, "view", "-b", "-L", bed_file, proper_paired_bam]
            with open(orig_bam, 'wb') as f:
                subprocess.run(cmd, stdout=f, check=True)
            log_and_print(f"✓ orig BAM created: {orig_bam}")
        
        # Index orig BAM
        orig_bai = f"{orig_bam}.bai"
        if not os.path.exists(orig_bai):
            log_and_print(f"Indexing orig BAM...")
            subprocess.run([sam_tools, "index", orig_bam], check=True)
            log_and_print(f"✓ orig BAM indexed")
        
        # Step 2: Create fetus BAM (size filter using fetus.awk)
        if not os.path.exists(fetus_bam):
            log_and_print(f"Creating fetus BAM: {fetus_bam}")
            
            # samtools view orig.bam | awk -f fetus.awk | samtools view -bt ref.fai -o fetus.bam -
            view_proc = subprocess.Popen(
                [sam_tools, "view", orig_bam],
                stdout=subprocess.PIPE
            )
            awk_proc = subprocess.Popen(
                ["awk", "-f", fetus_awk],
                stdin=view_proc.stdout,
                stdout=subprocess.PIPE
            )
            view_proc.stdout.close()
            
            with open(fetus_bam, 'wb') as f:
                subprocess.run(
                    [sam_tools, "view", "-bt", f"{ref_hg}.fai", "-o", fetus_bam, "-"],
                    stdin=awk_proc.stdout,
                    check=True
                )
            awk_proc.stdout.close()
            
            log_and_print(f"✓ fetus BAM created: {fetus_bam}")
        
        # Index fetus BAM
        fetus_bai = f"{fetus_bam}.bai"
        if not os.path.exists(fetus_bai):
            log_and_print(f"Indexing fetus BAM...")
            subprocess.run([sam_tools, "index", fetus_bam], check=True)
            log_and_print(f"✓ fetus BAM indexed")
        
        log_and_print("=== BAM filter processing completed ===")
        return orig_bam, fetus_bam
        
    except subprocess.CalledProcessError as e:
        log_and_print(f"ERROR in BAM filtering: {e}")
        return None, None


def run_wisecondor(sample_id, bam_file, reference_dir, sample_dir, plots_dir, bam_type="orig", force=False, skip_npz=False):
    """Run Wisecondor analysis
    
    Args:
        sample_dir: Sample analysis directory (e.g., /Work/NIPT/analysis/2508/GNCI25080071)
        bam_type: BAM type (orig, fetus, mom)
        force: If True, force re-run even if output files exist
        skip_npz: If True, skip NPZ file creation (use existing NPZ files)
    """
    log_and_print(f"=== Starting Wisecondor Analysis ({bam_type}) ===")
    
    # Create Output_WC/bam_type structure
    wc_base_dir = os.path.join(sample_dir, "Output_WC")
    wc_type_dir = os.path.join(wc_base_dir, bam_type)
    os.makedirs(wc_type_dir, exist_ok=True)
    
    # Extract BAM filename to determine NPZ name (e.g., proper_paired, of_fetus, etc.)
    bam_basename = os.path.basename(bam_file)
    bam_name_parts = bam_basename.replace(f"{sample_id}.", "").replace(".bam", "")  # e.g., "proper_paired" or "of_fetus"
    
    # NPZ file in base dir - name based on input BAM
    npz_file = os.path.join(wc_base_dir, f"{sample_id}.wc.{bam_name_parts}.npz")
    # Output files in bam_type subdirectory
    out_npz = os.path.join(wc_type_dir, f"{sample_id}.wc.{bam_type}.out.npz")
    report_file = os.path.join(wc_type_dir, f"{sample_id}.wc.{bam_type}.report.txt")
    plot_base = os.path.join(wc_type_dir, f"{sample_id}.wc.{bam_type}")  # wisecondor adds _z.png
    
    # Check if output NPZ already exists
    if not force and os.path.exists(out_npz):
        log_and_print(f"[SKIP] WC {bam_type} NPZ file already exists: {out_npz}")
        log_and_print(f"[SKIP] Skipping Wisecondor analysis for {bam_type}")
        return True
    
    # Reference file - different for orig and fetus
    if bam_type == "fetus":
        # For fetus, use fetus-specific reference
        wc_reference = os.path.join(reference_dir, "WC", "fetus_200k_of.npz")
    else:
        # For orig (proper_paired), use orig reference
        wc_reference = os.path.join(reference_dir, "WC", "orig_200k_proper_paired.npz")
    
    cyto_file = "/Work/NIPT/data/bed/common/cytoBand.txt"
    
    log_and_print(f"Using NPZ file: {npz_file}")
    log_and_print(f"Using reference: {wc_reference}")
    
    if not os.path.exists(wc_reference):
        log_and_print(f"WARNING: WC reference not found: {wc_reference}")
        return False
    
    try:
        # Step 1: Convert BAM to NPZ (only if base NPZ doesn't exist or force=True)
        if skip_npz:
            log_and_print(f"[SKIP-NPZ] Skipping NPZ file creation, using existing: {npz_file}")
            if not os.path.exists(npz_file):
                log_and_print(f"[ERROR] NPZ file not found: {npz_file}")
                log_and_print(f"[ERROR] Cannot skip NPZ creation when file does not exist")
                return False
        elif not force and os.path.exists(npz_file):
            log_and_print(f"[SKIP] WC base NPZ file already exists: {npz_file}")
        else:
            log_and_print("WC Step 1: Converting BAM to NPZ")
            binsize = 200000  # Default binsize
            
            convert_cmd = [
                os.environ.get("PYTHON2", "python2.7"),
                os.environ.get("WC", "/opt/wisecondor/wisecondor.py"),
                "convert",
                bam_file,
                npz_file,
                "-binsize", str(binsize)
            ]
            
            log_and_print(f"Running WC convert: {' '.join(convert_cmd)}")
            subprocess.run(convert_cmd, capture_output=True, text=True, check=True)
            log_and_print(f"WC NPZ created: {npz_file}")
        
        # Step 2: Run prediction with reference
        log_and_print("WC Step 2: Running prediction")
        test_cmd = [
            os.environ.get("PYTHON2", "python2.7"),
            os.environ.get("WC", "/opt/wisecondor/wisecondor.py"),
            "test",
            npz_file,
            out_npz,
            wc_reference
        ]
        
        log_and_print(f"Running WC test: {' '.join(test_cmd)}")
        subprocess.run(test_cmd, capture_output=True, text=True, check=True)
        
        # Step 3: Generate report
        log_and_print("WC Step 3: Generating report")
        report_cmd = [
            os.environ.get("PYTHON2", "python2.7"),
            os.environ.get("WC", "/opt/wisecondor/wisecondor.py"),
            "report",
            npz_file,
            out_npz
        ]
        
        log_and_print(f"Running WC report: {' '.join(report_cmd)}")
        result = subprocess.run(report_cmd, capture_output=True, text=True, check=True)
        with open(report_file, 'w') as f:
            f.write(result.stdout)
        log_and_print(f"Report saved: {report_file}")
        
        # Step 4: Generate plot (skip if skip_npz is True)
        if skip_npz:
            log_and_print("[SKIP-NPZ] Skipping WC plot generation")
        else:
            log_and_print("WC Step 4: Generating plot")
            
            # Build plot command - include cytofile only if it exists
            plot_cmd = [
                os.environ.get("PYTHON2", "python2.7"),
                os.environ.get("WC", "/opt/wisecondor/wisecondor.py"),
                "plot"
            ]
            
            # Add cytofile if it exists
            if os.path.exists(cyto_file):
                plot_cmd.extend(["-cytofile", cyto_file])
            else:
                log_and_print(f"WARNING: cytoband file not found: {cyto_file}, generating plot without it")
            
            plot_cmd.extend([
                "-filetype", "png",
                out_npz,
                plot_base
            ])
            
            log_and_print(f"Running WC plot: {' '.join(plot_cmd)}")
            subprocess.run(plot_cmd, capture_output=True, text=True, check=True)
            log_and_print(f"Plot saved: {plot_base}_z.png")
        
        log_and_print(f"Wisecondor completed: {out_npz}")
        log_and_print(f"Report saved: {report_file}")
        return True
        
    except subprocess.CalledProcessError as e:
        log_and_print(f"ERROR in Wisecondor: {e}")
        log_and_print(f"STDERR: {e.stderr}")
        return False


def run_wisecondorx(sample_id, bam_file, reference_dir, sample_dir, plots_dir, gender="F", bam_type="orig", force=False, ignore_zscore=False, config_file=None, skip_npz=False):
    """Run WisecondorX analysis with gender-specific reference
    
    Args:
        sample_dir: Sample analysis directory (e.g., /Work/NIPT/analysis/2508/GNCI25080071)
        gender: Fetal gender (M or F) for reference selection
        bam_type: BAM type (orig, fetus, mom)
        force: If True, force re-run even if output files exist
        ignore_zscore: If True, skip z-score threshold in aberrations.bed generation
        config_file: Path to pipeline config JSON file (to read WCX thresholds)
        skip_npz: If True, skip NPZ file creation (use existing NPZ files)
    """
    log_and_print(f"=== Starting WisecondorX Analysis ({bam_type}) ===")
    
    # Create Output_WCX/bam_type structure
    wcx_base_dir = os.path.join(sample_dir, "Output_WCX")
    wcx_type_dir = os.path.join(wcx_base_dir, bam_type)
    os.makedirs(wcx_type_dir, exist_ok=True)
    
    # Extract BAM filename to determine NPZ name (e.g., proper_paired, of_fetus, etc.)
    bam_basename = os.path.basename(bam_file)
    bam_name_parts = bam_basename.replace(f"{sample_id}.", "").replace(".bam", "")  # e.g., "proper_paired" or "of_fetus"
    
    # NPZ file in base dir - name based on input BAM
    npz_file = os.path.join(wcx_base_dir, f"{sample_id}.wcx.{bam_name_parts}.npz")
    
    # Check if output NPZ already exists (skip this check if skip_npz is True)
    if not skip_npz and not force and os.path.exists(npz_file):
        # Also check if aberrations.bed exists
        out_prefix = os.path.join(wcx_type_dir, f"{sample_id}.wcx.{bam_type}")
        aberrations_bed = f"{out_prefix}_aberrations.bed"
        if os.path.exists(aberrations_bed):
            log_and_print(f"[SKIP] WCX {bam_type} NPZ file already exists: {npz_file}")
            log_and_print(f"[SKIP] WCX {bam_type} output file already exists: {aberrations_bed}")
            log_and_print(f"[SKIP] Skipping WisecondorX analysis for {bam_type}")
            return True
    
    try:
        # Step 1: Convert BAM to NPZ
        if skip_npz:
            log_and_print(f"[SKIP-NPZ] Skipping NPZ file creation, using existing: {npz_file}")
            if not os.path.exists(npz_file):
                log_and_print(f"[ERROR] NPZ file not found: {npz_file}")
                log_and_print(f"[ERROR] Cannot skip NPZ creation when file does not exist")
                return False
        else:
            log_and_print("WCX Step 1: Converting BAM to NPZ")
            binsize = 200000  # Default binsize
            
            convert_cmd = [
                os.environ.get("WCX", "wisecondorx"),
                "convert",
                bam_file,
                npz_file,
                "--binsize", str(binsize)
            ]
            
            log_and_print(f"Running WCX convert: {' '.join(convert_cmd)}")
            subprocess.run(convert_cmd, capture_output=True, text=True, check=True)
            log_and_print(f"WCX NPZ created: {npz_file}")
        
        # Step 2: Predict with gender-specific reference
        # Different reference for orig and fetus
        if bam_type == "fetus":
            # For fetus, use fetus-specific reference with gender
            wcx_reference = os.path.join(reference_dir, "WCX", f"fetus_{gender}_200k_of.npz")
        else:
            # For orig (proper_paired), use orig reference with gender
            wcx_reference = os.path.join(reference_dir, "WCX", f"orig_{gender}_200k_proper_paired.npz")
        
        log_and_print(f"Using WCX reference: {wcx_reference}")
        
        if not os.path.exists(wcx_reference):
            log_and_print(f"WARNING: WCX reference not found: {wcx_reference}")
            return False
        
        log_and_print("WCX Step 2: Running prediction")
        # Output files in bam_type subdirectory
        out_prefix = os.path.join(wcx_type_dir, f"{sample_id}.wcx.{bam_type}")
        aberrations_bed = f"{out_prefix}_aberrations.bed"
        # WCX creates plots automatically in {out_prefix}.plots/
        
        # Read threshold from config if available
        wcx_threshold = None
        if not ignore_zscore and config_file and os.path.exists(config_file):
            try:
                import json
                with open(config_file, 'r') as f:
                    config_data = json.load(f)
                    if "WCX" in config_data:
                        wcx_config = config_data["WCX"]
                        if bam_type == "fetus":
                            wcx_threshold = wcx_config.get("fetus_threshold")
                        elif bam_type == "orig":
                            wcx_threshold = wcx_config.get("orig_threshold")
                        elif bam_type == "mom":
                            wcx_threshold = wcx_config.get("mom_threshold")
                        log_and_print(f"WCX threshold from config: {wcx_threshold} (ignore_zscore={ignore_zscore})")
            except Exception as e:
                log_and_print(f"Warning: Could not read WCX threshold from config: {e}")
        
        predict_cmd = [
            os.environ.get("WCX", "wisecondorx"),
            "predict",
            npz_file,
            wcx_reference,
            out_prefix,
            "--alpha", "0.01",
            "--bed",
            "--regions", "/Work/NIPT/data/empty_regions.txt"
        ]
        
        # Add --plot option (skip if skip_npz is True)
        if not skip_npz:
            predict_cmd.insert(-2, "--plot")  # Insert before --regions
        
        # Add --zscore option
        if ignore_zscore:
            # When ignoring zscore (for testing), remove --zscore option entirely
            # This should include all bins in aberrations.bed (both positive and negative z-scores)
            log_and_print(f"WCX ignoring config threshold, removing --zscore option to include all bins (ignore_zscore=True)")
        elif wcx_threshold is not None:
            # Use config threshold value
            predict_cmd.extend(["--zscore", str(wcx_threshold)])
            log_and_print(f"WCX using threshold: {wcx_threshold} (from config)")
        else:
            log_and_print(f"WCX using default threshold (no threshold in config)")
        
        log_and_print(f"Running WCX predict: {' '.join(predict_cmd)}")
        subprocess.run(predict_cmd, capture_output=True, text=True, check=True)
        
        log_and_print(f"WisecondorX completed: {aberrations_bed}")
        return True
        
    except subprocess.CalledProcessError as e:
        log_and_print(f"ERROR in WisecondorX: {e}")
        log_and_print(f"STDERR: {e.stderr}")
        return False


def count_reads_in_regions(bam_file, bed_file):
    """Count reads in BED file regions"""
    try:
        bamfile = pysam.AlignmentFile(bam_file, "rb")
        regions = []
        counts = []
        total_bases = 0

        # Read BED file
        with open(bed_file, "r") as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                        regions.append((chrom, start, end))
                        total_bases += end - start

        # Count reads in each region
        for chrom, start, end in regions:
            count = bamfile.count(chrom, start, end)
            counts.append(count)

        bamfile.close()
        return counts, regions, total_bases

    except Exception as e:
        log_and_print(f"Error counting reads in regions: {e}")
        return [], [], 0


def calculate_normalized_coverage(counts, regions):
    """Calculate normalized coverage for regions"""
    try:
        normalized = []
        for i, (chrom, start, end) in enumerate(regions):
            region_length = end - start
            if region_length > 0:
                coverage = counts[i] / region_length
                normalized.append(coverage)
        return normalized
    except Exception as e:
        log_and_print(f"Error calculating normalized coverage: {e}")
        return []


def calculate_yff(sample_name, bam_file, ff_config, lab_bed_paths):
    """Calculate Y-based fetal fraction"""
    log_and_print(f"Calculating Y-based fetal fraction for {sample_name}")

    # Get BED file paths
    bed_dir = Path("/Work/NIPT/data/bed/common")
    y_bed = bed_dir / "chrY_target_0.9.bed"
    a_bed = bed_dir / "autosome_control.bed"

    if not os.path.exists(y_bed):
        log_and_print(f"Y chromosome BED file not found: {y_bed}")
        return {
            "sample_name": sample_name,
            "YFF1": 0,
            "gd_1_value": 0,
            "gd_1_gender": "UNKNOWN",
            "status": "FAILED",
        }

    if not os.path.exists(a_bed):
        log_and_print(f"Autosome control BED file not found: {a_bed}")
        return {
            "sample_name": sample_name,
            "YFF1": 0,
            "gd_1_value": 0,
            "gd_1_gender": "UNKNOWN",
            "status": "FAILED",
        }

    try:
        # Count reads in Y chromosome regions
        log_and_print("Counting reads in Y chromosome regions...")
        y_counts, y_regions, y_total_bases = count_reads_in_regions(bam_file, str(y_bed))

        # Count reads in autosomal regions
        log_and_print("Counting reads in autosomal control regions...")
        a_counts, a_regions, a_total_bases = count_reads_in_regions(bam_file, str(a_bed))

        if not y_counts or not a_counts:
            log_and_print("Failed to count reads in regions")
            return {
                "sample_name": sample_name,
                "YFF1": 0,
                "gd_1_value": 0,
                "gd_1_gender": "UNKNOWN",
                "status": "FAILED",
            }

        # Calculate normalized coverage
        y_normalized = calculate_normalized_coverage(y_counts, y_regions)
        a_normalized = calculate_normalized_coverage(a_counts, a_regions)

        if not y_normalized or not a_normalized:
            log_and_print("Failed to calculate normalized coverage")
            return {
                "sample_name": sample_name,
                "YFF1": 0,
                "gd_1_value": 0,
                "gd_1_gender": "UNKNOWN",
                "status": "FAILED",
            }

        # Calculate median coverage
        y_median_coverage = np.median(y_normalized)
        a_median_coverage = np.median(a_normalized)

        log_and_print(f"Median Y coverage: {y_median_coverage:.6f}")
        log_and_print(f"Median Autosome coverage: {a_median_coverage:.6f}")
        log_and_print(f"Y regions count: {len(y_regions)}")
        log_and_print(f"Autosome regions count: {len(a_regions)}")
        log_and_print(f"Total Y reads: {sum(y_counts)}")
        log_and_print(f"Total Autosome reads: {sum(a_counts)}")

        if a_median_coverage <= 0:
            log_and_print("Autosomal median coverage is zero or negative. Cannot calculate YFF.")
            return {
                "sample_name": sample_name,
                "YFF1": 0,
                "gd_1_value": 0,
                "gd_1_gender": "UNKNOWN",
                "status": "FAILED",
            }

        # Calculate the ratio of Y to autosome
        y_to_a_ratio = y_median_coverage / a_median_coverage
        log_and_print(f"Y/A ratio: {y_to_a_ratio:.6f}")

        # Determine gender based on Y coverage
        gender_threshold = ff_config.get("gd_1_threshold", 0.01)

        log_and_print(f"Gender threshold: {gender_threshold}")
        if y_to_a_ratio < gender_threshold:
            # Female fetus - Y chromosome not present, YFF not applicable
            gd_1_gender = "XX"
            gd_1_value = y_to_a_ratio
            fetal_fraction = 0.0  # YFF is 0 for females
            status = "OK"
            log_and_print(f"Detected as Female (Y/A ratio {y_to_a_ratio:.6f} < threshold {gender_threshold:.6f})")
        else:
            # Male fetus - calculate YFF
            gd_1_gender = "XY"
            gd_1_value = y_to_a_ratio
            # Calculate fetal fraction assuming a male fetus
            # FF = 2 * (Y/A ratio) since only the male fetus contributes Y chromosome
            fetal_fraction = 2 * y_to_a_ratio * 100  # Convert to percentage
            status = "OK"
            log_and_print(f"Detected as Male (Y/A ratio {y_to_a_ratio:.6f} >= threshold {gender_threshold:.6f})")
            log_and_print(f"Calculated YFF: {fetal_fraction:.2f}%")

        return {
            "sample_name": sample_name,
            "YFF1": fetal_fraction,
            "gd_1_value": gd_1_value,
            "gd_1_gender": gd_1_gender,
            "status": "OK",
        }

    except Exception as e:
        log_and_print(f"Error in calculate_yff: {e}")
        return {
            "sample_name": sample_name,
            "YFF1": 0,
            "gd_1_value": 0,
            "gd_1_gender": "UNKNOWN",
            "status": f"ERROR: {str(e)}",
        }


def calculate_seqff(sample_name, bam_path, force=False):
    """Calculate seqFF fetal fraction
    
    Args:
        force: If True, force re-run even if output file exists
    """
    # Get sample directory from bam_path
    sample_dir = os.path.dirname(bam_path)
    output_dir = os.path.join(sample_dir, "Output_FF")
    os.makedirs(output_dir, exist_ok=True)
    seqff_txt = os.path.join(output_dir, f"{sample_name}.seqff.txt")
    
    # Check if seqff.txt already exists
    if not force and os.path.exists(seqff_txt):
        log_and_print(f"[SKIP] seqFF output file already exists: {seqff_txt}")
        # Try to parse existing result
        try:
            df = pd.read_csv(seqff_txt, index_col=0)
            seqff_value = float(df.loc["SeqFF", "x"]) * 100
            log_and_print(f"[SKIP] Using existing seqFF value: {seqff_value:.2f}%")
            return {"seqff_value": round(seqff_value, 2), "seqff_file": seqff_txt}
        except Exception as e:
            log_and_print(f"[SKIP] Could not parse existing seqFF result: {e}, will re-run")
            # Continue to re-run if parsing fails

    log_and_print("Running official R-based seqFF script")
    rscript = os.environ.get("Rscript", "Rscript")
    seqff_r_script = "/opt/nipt/bin/scripts/seqFF_R/seqff.r"

    original_cwd = os.getcwd()
    seqff_dir = "/opt/nipt/bin/scripts/seqFF_R"

    cmd = f"{rscript} --vanilla {seqff_r_script} -f {bam_path} -o {seqff_txt}"
    try:
        os.chdir(seqff_dir)
        subprocess.run(cmd, shell=True, check=True)
        log_and_print(f"Official seqFF output written to: {seqff_txt}")
    except subprocess.CalledProcessError:
        log_and_print("R-based seqFF failed to run.")
        os.chdir(original_cwd)
        return {"seqff_value": 0, "seqff_file": None}
    finally:
        os.chdir(original_cwd)

    # Parse seqFF result
    try:
        df = pd.read_csv(seqff_txt, index_col=0)
        seqff_value = float(df.loc["SeqFF", "x"]) * 100
        return {"seqff_value": round(seqff_value, 2), "seqff_file": seqff_txt}
    except Exception as e:
        log_and_print(f"Could not parse seqFF result: {e}")
        return {"seqff_value": 0, "seqff_file": seqff_txt}


def calculate_fetal_fraction(sample_id, bam_path, ff_config, lab_bed_paths, force=False):
    """Calculate fetal fraction using YFF and seqFF
    
    Args:
        force: If True, force re-run even if output files exist
    """
    log_and_print(f"\n=== Calculating Fetal Fraction for {sample_id} ===")
    
    ff_results = {}
    
    # YFF calculation
    yff_result = calculate_yff(sample_id, bam_path, ff_config, lab_bed_paths)
    if yff_result["status"] == "OK":
        yff_value = round(yff_result.get("YFF1", 0), 2)
        detected_gender = "MALE" if yff_result.get("gd_1_gender") == "XY" else "FEMALE"
        ff_results["yff"] = {
            "yff_value": yff_value,
            "gender": detected_gender,
            "status": "OK"
        }
        log_and_print(f"YFF: {yff_value:.2f}% (detected gender: {detected_gender})")
        
        # Warn if detected gender doesn't match expected gender from sample_id
        if "_M" in sample_id and detected_gender == "FEMALE":
            log_and_print(f"WARNING: Sample ID suggests Male (_M) but YFF detected Female")
        elif "_F" in sample_id and detected_gender == "MALE":
            log_and_print(f"WARNING: Sample ID suggests Female (_F) but YFF detected Male")
    else:
        ff_results["yff"] = {
            "yff_value": 0,
            "gender": "UNKNOWN",
            "status": yff_result.get("status", "FAILED")
        }
        log_and_print(f"YFF calculation failed: {yff_result.get('status', 'FAILED')}")
    
    # seqFF calculation
    seqff_result = calculate_seqff(sample_id, bam_path, force=force)
    ff_results["seqff"] = {
        "seqff_value": seqff_result.get("seqff_value", 0),
        "status": "OK" if seqff_result.get("seqff_value", 0) > 0 else "FAILED"
    }
    log_and_print(f"seqFF: {ff_results['seqff']['seqff_value']:.2f}%")
    
    # Determine final FF
    min_threshold = ff_config.get("ff_min_threshold", 2.0)
    max_threshold = ff_config.get("ff_max_threshold", 40.0)
    
    final_ff = None
    final_method = None
    
    # Priority: YFF for males, seqFF for females
    if ff_results["yff"]["status"] == "OK" and ff_results["yff"]["gender"] == "MALE":
        if ff_results["yff"]["yff_value"] >= min_threshold:
            final_ff = ff_results["yff"]["yff_value"]
            final_method = "YFF"
    
    if final_ff is None and ff_results["seqff"]["status"] == "OK":
        if min_threshold <= ff_results["seqff"]["seqff_value"] <= max_threshold:
            final_ff = ff_results["seqff"]["seqff_value"]
            final_method = "seqFF"
    
    if final_ff is None:
        # Fallback to available method
        if ff_results["yff"]["status"] == "OK":
            final_ff = ff_results["yff"]["yff_value"]
            final_method = "YFF"
        elif ff_results["seqff"]["status"] == "OK":
            final_ff = ff_results["seqff"]["seqff_value"]
            final_method = "seqFF"
    
    log_and_print(f"Final FF: {final_ff:.2f}% (method: {final_method})")
    
    return {
        "yff": ff_results["yff"]["yff_value"],
        "seqff": ff_results["seqff"]["seqff_value"],
        "fragment_ff": None,  # Not implemented yet
        "final_ff": final_ff,
        "final_method": final_method,
        "yff_status": ff_results["yff"].get("status"),
        "yff_gender": ff_results["yff"].get("gender"),
        "seqff_status": ff_results["seqff"].get("status")
    }


def update_json_with_ff(sample_analysis_dir, sample_id, ff_results):
    """Update JSON metadata file with calculated FF results"""
    import glob
    json_pattern = os.path.join(sample_analysis_dir, "*.json")
    json_files = glob.glob(json_pattern)
    
    if not json_files:
        log_and_print(f"No JSON file found to update for {sample_id}")
        return False
    
    json_file = json_files[0]
    try:
        with open(json_file, 'r') as f:
            metadata = json.load(f)
        
        # Update calculated_ff section
        metadata["calculated_ff"] = {
            "yff": ff_results.get("yff"),
            "seqff": ff_results.get("seqff"),
            "fragment_ff": ff_results.get("fragment_ff"),
            "final_ff": ff_results.get("final_ff"),
            "final_method": ff_results.get("final_method")
        }
        
        # Add YFF status and detected gender for debugging
        if "yff_status" in ff_results:
            metadata["calculated_ff"]["yff_status"] = ff_results.get("yff_status")
        if "yff_gender" in ff_results:
            metadata["calculated_ff"]["yff_gender"] = ff_results.get("yff_gender")
        if "seqff_status" in ff_results:
            metadata["calculated_ff"]["seqff_status"] = ff_results.get("seqff_status")
        
        # Add comparison with target FF if available
        if "target_parameters" in metadata:
            target_ff = metadata["target_parameters"].get("ff_target_percent")
            if target_ff and ff_results.get("final_ff"):
                diff = ff_results["final_ff"] - target_ff
                metadata["ff_comparison"] = {
                    "target_ff": target_ff,
                    "calculated_ff": ff_results["final_ff"],
                    "difference": round(diff, 2),
                    "method": ff_results.get("final_method")
                }
        
        # Write updated JSON
        with open(json_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        log_and_print(f"Updated JSON file with FF results: {json_file}")
        return True
        
    except Exception as e:
        log_and_print(f"Error updating JSON file: {e}")
        return False


def run_md_detection(sample_id, labcode, analysis_dir, output_dir, bed_dir, types=None, md_targets=None, ignore_min_length=False, ignore_zscore=False):
    """Run MD detection pipeline
    
    Args:
        sample_id: Sample ID
        labcode: Lab code
        analysis_dir: Analysis directory
        output_dir: Output directory
        bed_dir: BED directory
        types: List of BAM types to process (default: ["orig", "fetus"] for test mode)
        md_targets: List of MD targets to process (default: ["MD_Target_8"] for test mode)
        ignore_min_length: If True, skip minimum length check (only check z-score)
        ignore_zscore: If True, skip z-score threshold check (only check overlap with target)
    """
    log_and_print("=== Starting MD Detection Pipeline ===")
    
    # md_pipeline.py always uses dev module (required)
    if not DEV_MODULE_AVAILABLE or run_microdeletion_test_decision_pipeline_dev is None:
        log_and_print("ERROR: Dev module (process_md_result_dev) is required but not available!")
        log_and_print("ERROR: Please ensure process_md_result_dev.py exists in /Work/NIPT/bin/scripts/modules/")
        return False
    
    log_and_print(f"Using dev module (ignore_min_length={ignore_min_length}, ignore_zscore={ignore_zscore})")
    pipeline_func = run_microdeletion_test_decision_pipeline_dev
    
    try:
        # Prepare config (minimal version for MD-only)
        config = {
            "labcode": labcode,
            "analysis_dir": analysis_dir,
            "output_dir": output_dir
        }
        
        # Always use dev module with ignore_min_length and ignore_zscore support
        success = pipeline_func(
            sample_id,
            labcode,
            config,
            analysis_dir,
            output_dir,
            bed_dir,
            types=types,
            md_targets=md_targets,
            ignore_min_length=ignore_min_length,
            ignore_zscore=ignore_zscore
        )
        
        if success:
            log_and_print("MD detection completed successfully")
        else:
            log_and_print("MD detection completed with warnings")
        
        return True
        
    except Exception as e:
        log_and_print(f"ERROR in MD detection: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(description='Microdeletion-Only Analysis Pipeline')
    parser.add_argument('--sample_id', required=True, help='Sample ID')
    parser.add_argument('--work_dir', required=True, help='Work directory (e.g., 2507)')
    parser.add_argument('--labcode', required=True, help='Lab code (e.g., cordlife)')
    parser.add_argument('--analysis_dir', default='/Work/NIPT/analysis', help='Analysis directory')
    parser.add_argument('--output_dir', default='/Work/NIPT/output', help='Output directory')
    parser.add_argument('--data_dir', default='/Work/NIPT/data', help='Data directory')
    parser.add_argument('--fetal_gender', choices=['M', 'F', 'Male', 'Female', 'male', 'female'], 
                        help='Fetal gender (for artificial samples). M/Male or F/Female')
    parser.add_argument('--types', type=str, default='orig,fetus',
                        help='Comma-separated BAM types to process (default: orig,fetus)')
    parser.add_argument('--md_targets', type=str, default='MD_Target_8',
                        help='Comma-separated MD targets to process (default: MD_Target_8)')
    parser.add_argument('--filter_type', type=str, default='nf08',
                        help='BAM filter type prefix (default: nf08)')
    parser.add_argument('--ignore-min-length', action='store_true',
                        help='Ignore minimum length requirement (only check z-score threshold)')
    parser.add_argument('--ignore-zscore', action='store_true',
                        help='Ignore z-score threshold requirement (only check overlap with target)')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Force re-run even if output files already exist')
    parser.add_argument('-ro', '--result-only', action='store_true',
                        help='Skip WC/WCX analysis and only run MD detection (use existing NPZ files)')
    parser.add_argument('--skip-npz', action='store_true',
                        help='Skip NPZ file creation for WC/WCX (use existing NPZ files, faster execution)')
    
    args = parser.parse_args()
    
    sample_id = args.sample_id
    work_dir = args.work_dir
    labcode = args.labcode
    analysis_dir = args.analysis_dir
    output_dir = args.output_dir
    data_dir = args.data_dir
    fetal_gender_arg = args.fetal_gender
    
    # Parse types and md_targets from comma-separated strings
    types = [t.strip() for t in args.types.split(',')]
    md_targets = [mt.strip() for mt in args.md_targets.split(',')]
    filter_type = args.filter_type
    ignore_min_length = args.ignore_min_length
    ignore_zscore = args.ignore_zscore
    force_execution = args.force
    result_only = args.result_only
    skip_npz = args.skip_npz
    
    # Setup directories
    sample_analysis_dir = os.path.join(analysis_dir, work_dir, sample_id)
    sample_output_dir = os.path.join(output_dir, work_dir, sample_id)
    # Don't create plots_dir unless actually needed (plots are created in Output_WC/WCX subdirs)
    # plots_dir = os.path.join(sample_analysis_dir, "plots")
    
    os.makedirs(sample_analysis_dir, exist_ok=True)
    os.makedirs(sample_output_dir, exist_ok=True)
    # os.makedirs(plots_dir, exist_ok=True)  # Not needed - plots are created in Output_WC/WCX subdirs
    
    # Log file
    log_file = os.path.join(sample_analysis_dir, f"{sample_id}_md_analysis.log")
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)
    
    log_and_print("="*80)
    log_and_print(f"MD Pipeline Started for Sample: {sample_id}")
    log_and_print(f"Work Directory: {work_dir}")
    log_and_print(f"Lab Code: {labcode}")
    log_and_print("="*80)
    
    # Step 1: Find BAM file
    log_and_print("\n=== Step 1: Finding BAM File ===")
    bam_path, bam_type = find_bam_file(analysis_dir, work_dir, sample_id)
    
    if bam_path is None:
        log_and_print("ERROR: No BAM file found. Exiting.")
        return 1
    
    # Step 2: Create proper_paired.bam if needed
    if bam_type == 'sorted':
        log_and_print("\n=== Step 2: Creating proper_paired.bam ===")
        proper_paired_path = os.path.join(sample_analysis_dir, f"{sample_id}.proper_paired.bam")
        if not create_proper_paired_bam(bam_path, proper_paired_path):
            log_and_print("ERROR: Failed to create proper_paired.bam. Exiting.")
            return 1
        bam_path = proper_paired_path
    else:
        log_and_print("\n=== Step 2: Using existing proper_paired.bam ===")
    
    # Reference and BED directories
    reference_dir = os.path.join(data_dir, "refs", labcode)
    bed_dir = os.path.join(data_dir, "bed")
    
    # Get lab-specific BED paths
    lab_bed_paths = get_lab_bed_paths(labcode)
    
    # Step 2.5: Determine Fetal Gender
    log_and_print("\n=== Step 2.5: Determining Fetal Gender ===")
    gender = None
    
    # Priority 1: --fetal_gender argument (for artificial samples)
    if fetal_gender_arg:
        gender = fetal_gender_arg[0].upper()  # M or F
        log_and_print(f"Using fetal gender from argument: {'Male' if gender == 'M' else 'Female'}")
    
    # Priority 2: Read from JSON metadata file (for artificial samples)
    if gender is None:
        gender = read_gender_from_json(sample_analysis_dir, sample_id)
        if gender:
            log_and_print(f"Using fetal gender from JSON metadata: {'Male' if gender == 'M' else 'Female'}")
    
    # Priority 3: Read from Output_FF/gender.txt (gd_2)
    if gender is None:
        gender = read_gender_from_file(sample_analysis_dir, sample_id)
        if gender:
            log_and_print(f"Using fetal gender from gender.txt: {'Male' if gender == 'M' else 'Female'}")
    
    # Priority 4: Detect from BAM file (fallback)
    if gender is None:
        log_and_print("Gender not found in existing files, detecting from BAM...")
        gender = detect_gender_from_bam(bam_path)
        log_and_print(f"Using detected fetal gender: {'Male' if gender == 'M' else 'Female'}")
    
    # Final fallback
    if gender is None:
        log_and_print("WARNING: Could not determine gender, defaulting to Female")
        gender = "F"
    
    # If result_only mode, skip WC/WCX analysis and go directly to MD detection
    if result_only:
        log_and_print("\n=== Result-Only Mode: Skipping WC/WCX Analysis ===")
        log_and_print("Using existing NPZ files and output files for MD detection")
        log_and_print("Skipping Step 2.8 (BAM filters), Step 3 (Wisecondor), Step 4 (WisecondorX)")
        
        # Step 5: Run MD Detection directly
        log_and_print("\n=== Step 5: Running MD Detection ===")
        log_and_print(f"Processing types: {types}")
        log_and_print(f"Processing MD targets: {md_targets}")
        
        # Reference and BED directories (needed for MD detection)
        reference_dir = os.path.join(data_dir, "refs", labcode)
        bed_dir = os.path.join(data_dir, "bed")
        
        # Pass work_dir level path (not sample-specific) as expected by the MD detection function
        md_analysis_dir = os.path.join(analysis_dir, work_dir)
        md_output_dir = os.path.join(output_dir, work_dir)
        run_md_detection(sample_id, labcode, md_analysis_dir, md_output_dir, bed_dir, types=types, md_targets=md_targets, ignore_min_length=ignore_min_length, ignore_zscore=ignore_zscore)
        
        log_and_print("\n" + "="*80)
        log_and_print("MD Pipeline Completed Successfully (Result-Only Mode)!")
        log_and_print(f"Results saved in: {sample_analysis_dir}")
        log_and_print(f"WC output: {sample_analysis_dir}/Output_WC/")
        log_and_print(f"WCX output: {sample_analysis_dir}/Output_WCX/")
        log_and_print("="*80)
        
        return 0
    
    # Step 2.8: Process BAM filters (orig, fetus)
    log_and_print("\n=== Step 2.8: Processing BAM Filters ===")
    # Get BED file path from lab_bed_paths
    if filter_type in lab_bed_paths:
        bed_file = lab_bed_paths[filter_type]
    else:
        # Fallback to old method if filter_type not in lab_bed_paths
        bed_file = os.path.join(data_dir, f"bed/{filter_type}.bed")
        log_and_print(f"WARNING: filter_type '{filter_type}' not found in lab_bed_paths, using fallback: {bed_file}")
    
    log_and_print(f"Using BED file: {bed_file}")
    
    # Check if BED file exists
    if not os.path.exists(bed_file):
        log_and_print(f"ERROR: BED file not found: {bed_file}")
        return 1
    
    orig_bam, fetus_bam = process_bam_filters(sample_id, sample_analysis_dir, bed_file, filter_type)
    
    # Step 3: Run Wisecondor for orig and fetus
    log_and_print("\n=== Step 3: Running Wisecondor ===")
    
    # 3.1: WiseCONDOR on orig
    log_and_print("Running WiseCONDOR on orig BAM...")
    run_wisecondor(sample_id, bam_path, reference_dir, sample_analysis_dir, None, bam_type="orig", force=force_execution, skip_npz=skip_npz)
    
    # 3.2: WiseCONDOR on fetus (if available)
    if fetus_bam and os.path.exists(fetus_bam):
        log_and_print("Running WiseCONDOR on fetus BAM...")
        run_wisecondor(sample_id, fetus_bam, reference_dir, sample_analysis_dir, None, bam_type="fetus", force=force_execution, skip_npz=skip_npz)
    
    # Step 4: Run WisecondorX (with gender-specific reference)
    log_and_print("\n=== Step 4: Running WisecondorX ===")
    
    # Get config file path for WCX thresholds
    config_file = os.path.join("/Work/NIPT/config", labcode, "pipeline_config.json")
    
    # 4.1: WiseCONDORX on orig
    log_and_print("Running WiseCONDORX on orig BAM...")
    run_wisecondorx(sample_id, bam_path, reference_dir, sample_analysis_dir, None, gender, bam_type="orig", force=force_execution, ignore_zscore=ignore_zscore, config_file=config_file, skip_npz=skip_npz)
    
    # 4.2: WiseCONDORX on fetus (if available)
    if fetus_bam and os.path.exists(fetus_bam):
        log_and_print("Running WiseCONDORX on fetus BAM...")
        run_wisecondorx(sample_id, fetus_bam, reference_dir, sample_analysis_dir, None, gender, bam_type="fetus", force=force_execution, ignore_zscore=ignore_zscore, config_file=config_file, skip_npz=skip_npz)
    
    # Step 5: Run MD Detection
    log_and_print("\n=== Step 5: Running MD Detection ===")
    log_and_print(f"Processing types: {types}")
    log_and_print(f"Processing MD targets: {md_targets}")
    # Pass work_dir level path (not sample-specific) as expected by the MD detection function
    md_analysis_dir = os.path.join(analysis_dir, work_dir)
    md_output_dir = os.path.join(output_dir, work_dir)
    run_md_detection(sample_id, labcode, md_analysis_dir, md_output_dir, bed_dir, types=types, md_targets=md_targets, ignore_min_length=ignore_min_length, ignore_zscore=ignore_zscore)
    
    # Step 6: Calculate Fetal Fraction and update JSON
    log_and_print("\n=== Step 6: Calculating Fetal Fraction ===")
    # Get FF config (default values if config file doesn't exist)
    ff_config = {
        "ff_min_threshold": 2.0,
        "ff_max_threshold": 40.0,
        "gd_1_threshold": 0.01
    }
    
    # Try to load config from config file if available
    config_file = os.path.join("/Work/NIPT/config", labcode, "pipeline_config.json")
    if os.path.exists(config_file):
        try:
            with open(config_file, 'r') as f:
                config_data = json.load(f)
                if "FF_Gender_Config" in config_data:
                    ff_config.update(config_data["FF_Gender_Config"])
        except Exception as e:
            log_and_print(f"Warning: Could not load FF config from {config_file}: {e}")
    
    # Calculate FF
    ff_results = calculate_fetal_fraction(sample_id, bam_path, ff_config, lab_bed_paths, force=force_execution)
    
    # Update JSON file with FF results
    update_json_with_ff(sample_analysis_dir, sample_id, ff_results)
    
    # Mark completion
    marker_file = os.path.join(sample_analysis_dir, f"{sample_id}.md_pipeline_completed.marker")
    with open(marker_file, 'w') as f:
        f.write(f"MD Pipeline completed at {datetime.datetime.now()}\n")
    
    log_and_print("\n" + "="*80)
    log_and_print("MD Pipeline Completed Successfully!")
    log_and_print(f"Results saved in: {sample_analysis_dir}")
    log_and_print(f"WC output: {sample_analysis_dir}/Output_WC/")
    log_and_print(f"WCX output: {sample_analysis_dir}/Output_WCX/")
    if ff_results.get("final_ff"):
        log_and_print(f"Final FF: {ff_results['final_ff']:.2f}% (method: {ff_results['final_method']})")
    log_and_print("="*80)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

