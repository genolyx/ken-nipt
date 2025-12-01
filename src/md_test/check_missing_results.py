#!/usr/bin/env python3
"""
Check for samples that need reprocessing

Scans sample directories and checks for missing results:
- FF: Female needs seqFF only, Male needs both YFF and seqFF
- WC_orig, WC_fetus, WCX_orig, WCX_fetus

If WC, WCX results for both orig and fetus exist, proper_paired.bam only is OK.
Otherwise, of_orig.bam and of_fetus.bam are needed.
"""

import os
import sys
import json
import pandas as pd
from pathlib import Path
import argparse


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
            print(f"Warning: Could not read gender from {json_file}: {e}", file=sys.stderr)
    return 'UNKNOWN'


def check_ff_results(sample_dir, sample_name, gender):
    """Check FF results based on gender
    
    For Male samples, checks JSON file for YFF results in calculated_ff section.
    For Female samples, checks seqFF file.
    
    Returns:
        str: 'OK' if requirements met, 'MISSING' if missing required files
    """
    ff_dir = sample_dir / "Output_FF"
    if not ff_dir.exists():
        return "MISSING"
    
    # Check seqFF
    seqff_file = ff_dir / f"{sample_name}.seqff.txt"
    seqff_exists = False
    if seqff_file.exists():
        try:
            df = pd.read_csv(seqff_file, index_col=0)
            if "SeqFF" in df.index:
                seqff_exists = True
        except Exception as e:
            # If reading fails, try to check if file has content
            try:
                with open(seqff_file, 'r') as f:
                    content = f.read()
                    if "SeqFF" in content:
                        # File exists and contains SeqFF, but parsing failed
                        # Still consider it as existing (might be format issue)
                        seqff_exists = True
            except:
                pass
    
    # Check YFF - for Male samples, check JSON file first
    yff_exists = False
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
                yff_exists = True
        except:
            pass
    
    # Fallback: also check fetal_fraction.txt file if JSON doesn't have YFF
    if not yff_exists:
        ff_txt = ff_dir / f"{sample_name}.fetal_fraction.txt"
        if ff_txt.exists():
            try:
                with open(ff_txt, 'r') as f:
                    content = f.read()
                    if "YFF" in content or "yff" in content:
                        try:
                            ff_data = json.loads(content)
                            if ff_data.get("yff", 0) > 0 or ff_data.get("YFF", 0) > 0:
                                yff_exists = True
                        except:
                            if "YFF" in content and "0.00" not in content.split("YFF")[1][:10]:
                                yff_exists = True
            except:
                pass
    
    # Check requirements based on gender
    if gender == 'FEMALE':
        # Female: seqFF only required
        return "OK" if seqff_exists else "MISSING"
    elif gender == 'MALE':
        # Male: both YFF (from JSON) and seqFF required
        return "OK" if (yff_exists and seqff_exists) else "MISSING"
    else:
        # Unknown: at least seqFF should exist
        return "OK" if seqff_exists else "MISSING"


def check_wc_orig(sample_dir, sample_name):
    """Check WC orig results"""
    wc_orig_npz = sample_dir / "Output_WC" / "orig" / f"{sample_name}.wc.orig.out.npz"
    wc_orig_report = sample_dir / "Output_WC" / "orig" / f"{sample_name}.wc.orig.report.txt"
    if wc_orig_npz.exists() and wc_orig_report.exists():
        return "OK"
    return "MISSING"


def check_wc_fetus(sample_dir, sample_name):
    """Check WC fetus results"""
    wc_fetus_npz = sample_dir / "Output_WC" / "fetus" / f"{sample_name}.wc.fetus.out.npz"
    wc_fetus_report = sample_dir / "Output_WC" / "fetus" / f"{sample_name}.wc.fetus.report.txt"
    if wc_fetus_npz.exists() and wc_fetus_report.exists():
        return "OK"
    return "MISSING"


def check_wcx_orig(sample_dir, sample_name):
    """Check WCX orig results"""
    wcx_orig_npz = sample_dir / "Output_WCX" / f"{sample_name}.wcx.proper_paired.npz"
    wcx_orig_bed = sample_dir / "Output_WCX" / "orig" / f"{sample_name}.wcx.orig_aberrations.bed"
    wcx_orig_plots = sample_dir / "Output_WCX" / "orig" / f"{sample_name}.wcx.orig.plots"
    if wcx_orig_npz.exists() and wcx_orig_bed.exists() and wcx_orig_plots.exists():
        return "OK"
    return "MISSING"


def check_wcx_fetus(sample_dir, sample_name):
    """Check WCX fetus results"""
    wcx_fetus_npz = sample_dir / "Output_WCX" / f"{sample_name}.wcx.of_fetus.npz"
    wcx_fetus_bed = sample_dir / "Output_WCX" / "fetus" / f"{sample_name}.wcx.fetus_aberrations.bed"
    wcx_fetus_plots = sample_dir / "Output_WCX" / "fetus" / f"{sample_name}.wcx.fetus.plots"
    if wcx_fetus_npz.exists() and wcx_fetus_bed.exists() and wcx_fetus_plots.exists():
        return "OK"
    return "MISSING"


def check_sample(sample_dir, sample_name):
    """Check all results for a sample
    
    Returns:
        dict: Status for each component, or None if sample should be skipped
    """
    # Get gender
    gender = get_gender_from_json(sample_dir, sample_name)
    
    # Check all components
    ff_status = check_ff_results(sample_dir, sample_name, gender)
    wc_orig_status = check_wc_orig(sample_dir, sample_name)
    wc_fetus_status = check_wc_fetus(sample_dir, sample_name)
    wcx_orig_status = check_wcx_orig(sample_dir, sample_name)
    wcx_fetus_status = check_wcx_fetus(sample_dir, sample_name)
    
    # If all WC/WCX results exist, proper_paired.bam only is OK
    # Otherwise, check if proper_paired.bam exists (pipeline will create filtered BAMs)
    all_wc_wcx_complete = (wc_orig_status == "OK" and wc_fetus_status == "OK" and 
                          wcx_orig_status == "OK" and wcx_fetus_status == "OK")
    
    proper_paired_bam = sample_dir / f"{sample_name}.proper_paired.bam"
    if not all_wc_wcx_complete and not proper_paired_bam.exists():
        # Need BAM but doesn't exist - skip this sample (can't process)
        return None
    
    # Check if any component is missing
    if (ff_status == "MISSING" or wc_orig_status == "MISSING" or 
        wc_fetus_status == "MISSING" or wcx_orig_status == "MISSING" or 
        wcx_fetus_status == "MISSING"):
        return {
            'sample_name': sample_name,
            'gender': gender,
            'FF': ff_status,
            'WC_orig': wc_orig_status,
            'WC_fetus': wc_fetus_status,
            'WCX_orig': wcx_orig_status,
            'WCX_fetus': wcx_fetus_status
        }
    
    # All complete
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Check for samples that need reprocessing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Check current directory
  python check_missing_results.py
  
  # Check specific directory
  python check_missing_results.py -d /path/to/analysis/md_validation/WBS
  
  # Output to file
  python check_missing_results.py -o missing_samples.tsv
        """
    )
    
    parser.add_argument('-d', '--dir', default='.', 
                       help='Directory to scan (default: current directory)')
    parser.add_argument('-o', '--output', 
                       help='Output TSV file (default: stdout)')
    parser.add_argument('--root-dir',
                       help='Root directory (if samples are in analysis/{work_dir} structure)')
    parser.add_argument('--work-dir',
                       help='Work directory (if samples are in analysis/{work_dir} structure)')
    
    args = parser.parse_args()
    
    # Determine base directory
    if args.root_dir and args.work_dir:
        base_dir = Path(args.root_dir) / "analysis" / args.work_dir
    else:
        base_dir = Path(args.dir).resolve()
    
    if not base_dir.exists():
        print(f"Error: Directory not found: {base_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Scanning directory: {base_dir}", file=sys.stderr)
    
    # Find all sample directories
    missing_samples = []
    sample_dirs = [d for d in base_dir.iterdir() if d.is_dir()]
    
    for sample_dir in sorted(sample_dirs):
        sample_name = sample_dir.name
        
        # Skip hidden directories
        if sample_name.startswith('.'):
            continue
        
        result = check_sample(sample_dir, sample_name)
        if result:
            missing_samples.append(result)
    
    # Output results
    if missing_samples:
        df = pd.DataFrame(missing_samples)
        # Reorder columns
        df = df[['sample_name', 'FF', 'WC_orig', 'WC_fetus', 'WCX_orig', 'WCX_fetus']]
        
        if args.output:
            df.to_csv(args.output, sep='\t', index=False)
            print(f"\nFound {len(missing_samples)} samples needing reprocessing", file=sys.stderr)
            print(f"Results saved to: {args.output}", file=sys.stderr)
        else:
            print(df.to_string(index=False))
            print(f"\nFound {len(missing_samples)} samples needing reprocessing", file=sys.stderr)
    else:
        print("All samples have complete results!", file=sys.stderr)
    
    return 0 if missing_samples else 0


if __name__ == "__main__":
    sys.exit(main())

