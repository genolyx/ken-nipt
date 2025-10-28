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

import pysam

sys.path.append("/Work/NIPT/bin")
sys.path.append("/Work/NIPT/bin/scripts/modules")
try:
    from process_md_result import run_microdeletion_decision_pipeline, run_microdeletion_test_decision_pipeline
except ImportError as e:
    logging.warning(f"Could not import run_microdeletion functions: {e}")
    run_microdeletion_decision_pipeline = None
    run_microdeletion_test_decision_pipeline = None

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
    fetus_awk = "/Work/NIPT/src/fetus.awk"
    
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


def run_wisecondor(sample_id, bam_file, reference_dir, sample_dir, plots_dir, bam_type="orig"):
    """Run Wisecondor analysis
    
    Args:
        sample_dir: Sample analysis directory (e.g., /Work/NIPT/analysis/2508/GNCI25080071)
        bam_type: BAM type (orig, fetus, mom)
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
        # Step 1: Convert BAM to NPZ
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
        
        # Step 4: Generate plot
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
        
        log_and_print(f"Wisecondor completed: {out_npz}")
        log_and_print(f"Report saved: {report_file}")
        log_and_print(f"Plot saved: {plot_base}_z.png")
        return True
        
    except subprocess.CalledProcessError as e:
        log_and_print(f"ERROR in Wisecondor: {e}")
        log_and_print(f"STDERR: {e.stderr}")
        return False


def run_wisecondorx(sample_id, bam_file, reference_dir, sample_dir, plots_dir, gender="F", bam_type="orig"):
    """Run WisecondorX analysis with gender-specific reference
    
    Args:
        sample_dir: Sample analysis directory (e.g., /Work/NIPT/analysis/2508/GNCI25080071)
        gender: Fetal gender (M or F) for reference selection
        bam_type: BAM type (orig, fetus, mom)
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
    
    try:
        # Step 1: Convert BAM to NPZ
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
        
        predict_cmd = [
            os.environ.get("WCX", "wisecondorx"),
            "predict",
            npz_file,
            wcx_reference,
            out_prefix,
            "--plot",
            "--bed",
            "--regions", "/Work/NIPT/data/empty_regions.txt"
        ]
        
        log_and_print(f"Running WCX predict: {' '.join(predict_cmd)}")
        subprocess.run(predict_cmd, capture_output=True, text=True, check=True)
        
        log_and_print(f"WisecondorX completed: {aberrations_bed}")
        return True
        
    except subprocess.CalledProcessError as e:
        log_and_print(f"ERROR in WisecondorX: {e}")
        log_and_print(f"STDERR: {e.stderr}")
        return False


def run_md_detection(sample_id, labcode, analysis_dir, output_dir, bed_dir, types=None, md_targets=None):
    """Run MD detection pipeline
    
    Args:
        sample_id: Sample ID
        labcode: Lab code
        analysis_dir: Analysis directory
        output_dir: Output directory
        bed_dir: BED directory
        types: List of BAM types to process (default: ["orig", "fetus"] for test mode)
        md_targets: List of MD targets to process (default: ["MD_Target_8"] for test mode)
    """
    log_and_print("=== Starting MD Detection Pipeline ===")
    
    if run_microdeletion_test_decision_pipeline is None:
        log_and_print("ERROR: MD detection module not available")
        return False
    
    try:
        # Prepare config (minimal version for MD-only)
        config = {
            "labcode": labcode,
            "analysis_dir": analysis_dir,
            "output_dir": output_dir
        }
        
        # Use test version with configurable types and targets
        success = run_microdeletion_test_decision_pipeline(
            sample_id,
            labcode,
            config,
            analysis_dir,
            output_dir,
            bed_dir,
            types=types,
            md_targets=md_targets
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
    
    # Setup directories
    sample_analysis_dir = os.path.join(analysis_dir, work_dir, sample_id)
    sample_output_dir = os.path.join(output_dir, work_dir, sample_id)
    plots_dir = os.path.join(sample_analysis_dir, "plots")
    
    os.makedirs(sample_analysis_dir, exist_ok=True)
    os.makedirs(sample_output_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)
    
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
    
    # Step 2.5: Determine Fetal Gender
    log_and_print("\n=== Step 2.5: Determining Fetal Gender ===")
    gender = None
    
    # Priority 1: --fetal_gender argument (for artificial samples)
    if fetal_gender_arg:
        gender = fetal_gender_arg[0].upper()  # M or F
        log_and_print(f"Using fetal gender from argument: {'Male' if gender == 'M' else 'Female'}")
    
    # Priority 2: Read from Output_FF/gender.txt (gd_2)
    if gender is None:
        gender = read_gender_from_file(sample_analysis_dir, sample_id)
        if gender:
            log_and_print(f"Using fetal gender from gender.txt: {'Male' if gender == 'M' else 'Female'}")
    
    # Priority 3: Detect from BAM file (fallback)
    if gender is None:
        log_and_print("Gender not found in existing files, detecting from BAM...")
        gender = detect_gender_from_bam(bam_path)
        log_and_print(f"Using detected fetal gender: {'Male' if gender == 'M' else 'Female'}")
    
    # Final fallback
    if gender is None:
        log_and_print("WARNING: Could not determine gender, defaulting to Female")
        gender = "F"
    
    # Step 2.8: Process BAM filters (orig, fetus)
    log_and_print("\n=== Step 2.8: Processing BAM Filters ===")
    bed_file = os.path.join(data_dir, f"bed/{filter_type}.bed")
    log_and_print(f"Using BED file: {bed_file}")
    orig_bam, fetus_bam = process_bam_filters(sample_id, sample_analysis_dir, bed_file, filter_type)
    
    # Step 3: Run Wisecondor for orig and fetus
    log_and_print("\n=== Step 3: Running Wisecondor ===")
    
    # 3.1: WiseCONDOR on orig
    log_and_print("Running WiseCONDOR on orig BAM...")
    run_wisecondor(sample_id, bam_path, reference_dir, sample_analysis_dir, plots_dir, bam_type="orig")
    
    # 3.2: WiseCONDOR on fetus (if available)
    if fetus_bam and os.path.exists(fetus_bam):
        log_and_print("Running WiseCONDOR on fetus BAM...")
        run_wisecondor(sample_id, fetus_bam, reference_dir, sample_analysis_dir, plots_dir, bam_type="fetus")
    
    # Step 4: Run WisecondorX (with gender-specific reference)
    log_and_print("\n=== Step 4: Running WisecondorX ===")
    
    # 4.1: WiseCONDORX on orig
    log_and_print("Running WiseCONDORX on orig BAM...")
    run_wisecondorx(sample_id, bam_path, reference_dir, sample_analysis_dir, plots_dir, gender, bam_type="orig")
    
    # 4.2: WiseCONDORX on fetus (if available)
    if fetus_bam and os.path.exists(fetus_bam):
        log_and_print("Running WiseCONDORX on fetus BAM...")
        run_wisecondorx(sample_id, fetus_bam, reference_dir, sample_analysis_dir, plots_dir, gender, bam_type="fetus")
    
    # Step 5: Run MD Detection
    log_and_print("\n=== Step 5: Running MD Detection ===")
    log_and_print(f"Processing types: {types}")
    log_and_print(f"Processing MD targets: {md_targets}")
    # Pass work_dir level path (not sample-specific) as expected by the MD detection function
    md_analysis_dir = os.path.join(analysis_dir, work_dir)
    md_output_dir = os.path.join(output_dir, work_dir)
    run_md_detection(sample_id, labcode, md_analysis_dir, md_output_dir, bed_dir, types=types, md_targets=md_targets)
    
    # Mark completion
    marker_file = os.path.join(sample_analysis_dir, f"{sample_id}.md_pipeline_completed.marker")
    with open(marker_file, 'w') as f:
        f.write(f"MD Pipeline completed at {datetime.datetime.now()}\n")
    
    log_and_print("\n" + "="*80)
    log_and_print("MD Pipeline Completed Successfully!")
    log_and_print(f"Results saved in: {sample_analysis_dir}")
    log_and_print(f"Plots saved in: {plots_dir}")
    log_and_print("="*80)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

