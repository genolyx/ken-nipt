#!/usr/bin/env python3
"""
Calculate performance metrics (Sensitivity, Specificity, PPV, NPV) by z-score threshold

Reads zscore_extraction.tsv and analyzes all detected regions from WC/WCX files
to calculate TP, FP, FN, TN for each threshold.
"""

import argparse
import json
from pathlib import Path
import pandas as pd
import numpy as np
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def parse_deletion_from_json(json_file: Path) -> dict:
    """Parse deletion information from JSON file"""
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        deletion = data.get('deletion', {})
        if not deletion:
            return None
        
        return {
            'chromosome': str(deletion.get('chromosome', '')).replace('chr', ''),
            'start': deletion.get('start'),
            'end': deletion.get('end'),
        }
    except Exception as e:
        logger.warning(f"Failed to parse JSON {json_file}: {e}")
        return None


def check_overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
    """Check if two regions overlap"""
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    return overlap_start < overlap_end


def parse_wc_report_all_regions(report_file: Path) -> list:
    """Parse WC report.txt and extract all detected regions
    
    Returns list of dicts: [{'chr': str, 'start': int, 'end': int, 'zscore': float}, ...]
    """
    regions = []
    
    if not report_file.exists():
        return regions
    
    try:
        with open(report_file, 'r') as f:
            lines = f.readlines()
        
        in_test_section = False
        
        for line in lines:
            if '# test results:' in line.lower():
                in_test_section = True
                continue
            
            if in_test_section:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if 'z-score' in line.lower() and 'effect' in line.lower():
                    continue
                
                fields = line.split('\t')
                if len(fields) < 4:
                    fields = line.split()
                    if len(fields) < 4:
                        continue
                
                try:
                    zscore = float(fields[0])
                    location = fields[3]
                    
                    if ':' in location:
                        chr_part, coord_part = location.split(':', 1)
                        if '-' in coord_part:
                            start_str, end_str = coord_part.split('-', 1)
                            chr_name = str(chr_part).replace('chr', '')
                            start = int(start_str)
                            end = int(end_str)
                            
                            regions.append({
                                'chr': chr_name,
                                'start': start,
                                'end': end,
                                'zscore': zscore
                            })
                except (ValueError, IndexError):
                    continue
    except Exception as e:
        logger.warning(f"Failed to parse WC report {report_file}: {e}")
    
    return regions


def parse_wcx_bed_all_regions(bed_file: Path) -> list:
    """Parse WCX aberrations.bed and extract all detected regions
    
    Returns list of dicts: [{'chr': str, 'start': int, 'end': int, 'zscore': float}, ...]
    """
    regions = []
    
    if not bed_file.exists():
        return regions
    
    try:
        with open(bed_file, 'r') as f:
            lines = f.readlines()
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 5:
                continue
            
            try:
                chr_name = str(fields[0]).replace('chr', '')
                start = int(fields[1])
                end = int(fields[2])
                
                # Find z-score (usually column 4 or 5)
                zscore = None
                for i in range(3, min(len(fields), 6)):
                    try:
                        val = float(fields[i])
                        if abs(val) >= 2.0:  # Likely z-score
                            zscore = val
                            break
                    except ValueError:
                        continue
                
                if zscore is not None:
                    regions.append({
                        'chr': chr_name,
                        'start': start,
                        'end': end,
                        'zscore': zscore
                    })
            except (ValueError, IndexError):
                continue
    except Exception as e:
        logger.warning(f"Failed to parse WCX bed {bed_file}: {e}")
    
    return regions


def calculate_metrics_for_sample(sample_dir: Path, sample_name: str, 
                                  expected_deletion: dict, method: str, 
                                  output_type: str, threshold: float, 
                                  min_length: int = None) -> dict:
    """Calculate TP, FP, FN, TN for a single sample at given threshold
    
    Logic:
    - TP: Target 영역이 검출됨
    - FN: Target 영역이 검출 안됨
    - FP: 다른 영역이 검출됨 (Target 영역이 아닌 다른 영역, min_length 이상인 경우만)
    - TN: 다른 영역이 검출 안됨
    
    Args:
        min_length: Minimum length (bp) for a region to be considered as FP.
                    If None, all regions are considered.
    
    Returns: {'tp': int, 'fp': int, 'fn': int, 'tn': int}
    """
    result = {'tp': 0, 'fp': 0, 'fn': 0, 'tn': 0}
    
    # Get all detected regions
    if method == 'WCX':
        bed_file = sample_dir / f"Output_WCX" / output_type / f"{sample_name}.wcx.{output_type}_aberrations.bed"
        detected_regions = parse_wcx_bed_all_regions(bed_file)
    else:  # WC
        report_file = sample_dir / f"Output_WC" / output_type / f"{sample_name}.wc.{output_type}.report.txt"
        detected_regions = parse_wc_report_all_regions(report_file)
    
    # Filter by threshold (for deletions, z-score <= threshold)
    filtered_regions = [r for r in detected_regions if r['zscore'] <= threshold]
    
    expected_chr = str(expected_deletion['chromosome']).replace('chr', '')
    expected_start = expected_deletion['start']
    expected_end = expected_deletion['end']
    
    # Check if target region is detected
    target_detected = False
    other_regions_detected = False
    
    for region in filtered_regions:
        region_length = region['end'] - region['start']
        
        # Check if this region overlaps with target deletion
        if region['chr'] == expected_chr:
            if check_overlap(expected_start, expected_end, region['start'], region['end']):
                target_detected = True
                continue
        
        # This is a different region (not target)
        # Only count as FP if min_length is None or region length >= min_length
        if min_length is None or region_length >= min_length:
            other_regions_detected = True
    
    # Assign TP/FN/FP/TN
    if target_detected:
        result['tp'] = 1
    else:
        result['fn'] = 1
    
    if other_regions_detected:
        result['fp'] = 1
    else:
        result['tn'] = 1
    
    return result


def calculate_metrics_at_threshold_or(df: pd.DataFrame, analysis_dir: Path,
                                     methods: list, threshold: float,
                                     min_length: int = None, group_by_ff_length: bool = False) -> dict:
    """Calculate performance metrics at a given threshold using OR logic across multiple methods
    
    Args:
        methods: List of (method, output_type) tuples, e.g., [('WC', 'orig'), ('WCX', 'orig')]
        OR logic: If any method detects the target, it's TP. If any method detects other regions, it's FP.
    
    Returns:
        dict or list of dicts (if group_by_ff_length=True)
    """
    
    if group_by_ff_length:
        # Group by FF and deletion_length_mb
        results = []
        
        # Get all samples (including those without z-score)
        all_samples = df.copy()
        
        # Group by FF and deletion_length_mb
        if 'ff' not in all_samples.columns or 'deletion_length_mb' not in all_samples.columns:
            logger.warning("FF or deletion_length_mb columns not found. Cannot group by FF/length.")
            group_by_ff_length = False
        else:
            # Convert to numeric
            all_samples['ff'] = pd.to_numeric(all_samples['ff'], errors='coerce')
            all_samples['deletion_length_mb'] = pd.to_numeric(all_samples['deletion_length_mb'], errors='coerce')
            
            # Remove rows with missing FF or deletion_length_mb
            all_samples = all_samples[all_samples[['ff', 'deletion_length_mb']].notna().all(axis=1)]
            
            grouped = all_samples.groupby(['ff', 'deletion_length_mb'])
            
            for (ff, del_length_mb), group_df in grouped:
                tp_total = 0
                fp_total = 0
                fn_total = 0
                tn_total = 0
                
                # Process each sample in this group
                for _, row in group_df.iterrows():
                    sample_name = row['sample_name']
                    sample_dir = analysis_dir / sample_name
                    
                    expected_deletion = {
                        'chromosome': row['expected_deletion_chr'],
                        'start': row['expected_deletion_start'],
                        'end': row['expected_deletion_end']
                    }
                    
                    # Check all methods with OR logic
                    target_detected_any = False
                    other_regions_detected_any = False
                    has_any_zscore = False
                    
                    for method, output_type in methods:
                        method_name = f'{method}_{output_type}'
                        zscore_col = method_name + '_zscore'
                        
                        # Check if sample has z-score for this method
                        has_zscore = pd.notna(row.get(zscore_col))
                        if has_zscore:
                            has_any_zscore = True
                        
                        # Calculate metrics for this method
                        metrics = calculate_metrics_for_sample(
                            sample_dir, sample_name, expected_deletion,
                            method, output_type, threshold, min_length
                        )
                        
                        # OR logic: if any method detects target, it's TP
                        if metrics['tp'] == 1:
                            target_detected_any = True
                        
                        # OR logic: if any method detects other regions, it's FP
                        if metrics['fp'] == 1:
                            other_regions_detected_any = True
                    
                    # If no method has z-score, check files directly
                    if not has_any_zscore:
                        for method, output_type in methods:
                            if method == 'WCX':
                                bed_file = sample_dir / f"Output_WCX" / output_type / f"{sample_name}.wcx.{output_type}_aberrations.bed"
                                detected_regions = parse_wcx_bed_all_regions(bed_file)
                            else:
                                report_file = sample_dir / f"Output_WC" / output_type / f"{sample_name}.wc.{output_type}.report.txt"
                                detected_regions = parse_wc_report_all_regions(report_file)
                            
                            filtered_regions = [r for r in detected_regions if r['zscore'] <= threshold]
                            if min_length is not None:
                                filtered_regions = [r for r in filtered_regions if (r['end'] - r['start']) >= min_length]
                            
                            if filtered_regions:
                                other_regions_detected_any = True
                    
                    # Assign TP/FN/FP/TN based on OR logic
                    if target_detected_any:
                        tp_total += 1
                    else:
                        fn_total += 1
                    
                    if other_regions_detected_any:
                        fp_total += 1
                    else:
                        tn_total += 1
                
                # Calculate metrics for this group
                total_samples = len(group_df)
                sensitivity = tp_total / (tp_total + fn_total) if (tp_total + fn_total) > 0 else 0.0
                specificity = tn_total / (tn_total + fp_total) if (tn_total + fp_total) > 0 else 0.0
                ppv = tp_total / (tp_total + fp_total) if (tp_total + fp_total) > 0 else 0.0
                npv = tn_total / (tn_total + fn_total) if (tn_total + fn_total) > 0 else 0.0
                
                method_names = '_OR_'.join([f'{m}_{o}' for m, o in methods])
                results.append({
                    'threshold': threshold,
                    'method': method_names,
                    'ff': ff,
                    'deletion_length_mb': del_length_mb,
                    'tp': tp_total,
                    'fp': fp_total,
                    'fn': fn_total,
                    'tn': tn_total,
                    'sensitivity': sensitivity,
                    'specificity': specificity,
                    'ppv': ppv,
                    'npv': npv,
                    'total_samples': total_samples
                })
            
            return results
    
    # Original logic: calculate for all samples together
    # Get all samples (including those without z-score)
    all_samples = df.copy()
    
    tp_total = 0
    fp_total = 0
    fn_total = 0
    tn_total = 0
    
    # Process each sample
    for _, row in all_samples.iterrows():
        sample_name = row['sample_name']
        sample_dir = analysis_dir / sample_name
        
        expected_deletion = {
            'chromosome': row['expected_deletion_chr'],
            'start': row['expected_deletion_start'],
            'end': row['expected_deletion_end']
        }
        
        # Check all methods with OR logic
        target_detected_any = False
        other_regions_detected_any = False
        has_any_zscore = False
        
        for method, output_type in methods:
            method_name = f'{method}_{output_type}'
            zscore_col = method_name + '_zscore'
            
            # Check if sample has z-score for this method
            has_zscore = pd.notna(row.get(zscore_col))
            if has_zscore:
                has_any_zscore = True
            
            # Calculate metrics for this method
            metrics = calculate_metrics_for_sample(
                sample_dir, sample_name, expected_deletion,
                method, output_type, threshold, min_length
            )
            
            # OR logic: if any method detects target, it's TP
            if metrics['tp'] == 1:
                target_detected_any = True
            
            # OR logic: if any method detects other regions, it's FP
            if metrics['fp'] == 1:
                other_regions_detected_any = True
        
        # If no method has z-score, check files directly
        if not has_any_zscore:
            for method, output_type in methods:
                if method == 'WCX':
                    bed_file = sample_dir / f"Output_WCX" / output_type / f"{sample_name}.wcx.{output_type}_aberrations.bed"
                    detected_regions = parse_wcx_bed_all_regions(bed_file)
                else:
                    report_file = sample_dir / f"Output_WC" / output_type / f"{sample_name}.wc.{output_type}.report.txt"
                    detected_regions = parse_wc_report_all_regions(report_file)
                
                filtered_regions = [r for r in detected_regions if r['zscore'] <= threshold]
                if min_length is not None:
                    filtered_regions = [r for r in filtered_regions if (r['end'] - r['start']) >= min_length]
                
                if filtered_regions:
                    other_regions_detected_any = True
        
        # Assign TP/FN/FP/TN based on OR logic
        if target_detected_any:
            tp_total += 1
        else:
            fn_total += 1
        
        if other_regions_detected_any:
            fp_total += 1
        else:
            tn_total += 1
    
    # Calculate performance metrics
    total_samples = len(all_samples)
    sensitivity = tp_total / (tp_total + fn_total) if (tp_total + fn_total) > 0 else 0.0
    specificity = tn_total / (tn_total + fp_total) if (tn_total + fp_total) > 0 else 0.0
    ppv = tp_total / (tp_total + fp_total) if (tp_total + fp_total) > 0 else 0.0
    npv = tn_total / (tn_total + fn_total) if (tn_total + fn_total) > 0 else 0.0
    
    method_names = '_OR_'.join([f'{m}_{o}' for m, o in methods])
    return {
        'threshold': threshold,
        'method': method_names,
        'tp': tp_total,
        'fp': fp_total,
        'fn': fn_total,
        'tn': tn_total,
        'sensitivity': sensitivity,
        'specificity': specificity,
        'ppv': ppv,
        'npv': npv,
        'total_samples': total_samples
    }


def calculate_metrics_at_threshold(df: pd.DataFrame, analysis_dir: Path,
                                   method: str, output_type: str, threshold: float,
                                   min_length: int = None, group_by_ff_length: bool = False) -> dict:
    """Calculate performance metrics at a given threshold
    
    Args:
        group_by_ff_length: If True, calculate metrics grouped by FF and deletion_length_mb
    """
    
    method_name = f'{method}_{output_type}'
    zscore_col = method_name + '_zscore'
    
    if group_by_ff_length:
        # Group by FF and deletion_length_mb
        results = []
        
        # Get all samples (including those without z-score)
        all_samples = df.copy()
        
        # Group by FF and deletion_length_mb
        if 'ff' not in all_samples.columns or 'deletion_length_mb' not in all_samples.columns:
            logger.warning("FF or deletion_length_mb columns not found. Cannot group by FF/length.")
            group_by_ff_length = False
        else:
            # Convert to numeric
            all_samples['ff'] = pd.to_numeric(all_samples['ff'], errors='coerce')
            all_samples['deletion_length_mb'] = pd.to_numeric(all_samples['deletion_length_mb'], errors='coerce')
            
            # Remove rows with missing FF or deletion_length_mb
            all_samples = all_samples[all_samples[['ff', 'deletion_length_mb']].notna().all(axis=1)]
            
            grouped = all_samples.groupby(['ff', 'deletion_length_mb'])
            
            for (ff, del_length_mb), group_df in grouped:
                tp_total = 0
                fp_total = 0
                fn_total = 0
                tn_total = 0
                
                # Process each sample in this group
                for _, row in group_df.iterrows():
                    sample_name = row['sample_name']
                    sample_dir = analysis_dir / sample_name
                    
                    expected_deletion = {
                        'chromosome': row['expected_deletion_chr'],
                        'start': row['expected_deletion_start'],
                        'end': row['expected_deletion_end']
                    }
                    
                    # Check if z-score exists
                    has_zscore = pd.notna(row.get(zscore_col))
                    
                    if has_zscore:
                        # Sample has z-score, check detection
                        metrics = calculate_metrics_for_sample(
                            sample_dir, sample_name, expected_deletion,
                            method, output_type, threshold, min_length
                        )
                        tp_total += metrics['tp']
                        fp_total += metrics['fp']
                        fn_total += metrics['fn']
                        tn_total += metrics['tn']
                    else:
                        # Sample has no z-score = not detected = FN
                        # Also check if there are any FP regions (need to check report/bed files)
                        if method == 'WCX':
                            bed_file = sample_dir / f"Output_WCX" / output_type / f"{sample_name}.wcx.{output_type}_aberrations.bed"
                            detected_regions = parse_wcx_bed_all_regions(bed_file)
                        else:
                            report_file = sample_dir / f"Output_WC" / output_type / f"{sample_name}.wc.{output_type}.report.txt"
                            detected_regions = parse_wc_report_all_regions(report_file)
                        
                        # Filter by threshold and min_length
                        filtered_regions = [r for r in detected_regions if r['zscore'] <= threshold]
                        if min_length is not None:
                            filtered_regions = [r for r in filtered_regions if (r['end'] - r['start']) >= min_length]
                        
                        # No target detection (FN)
                        fn_total += 1
                        
                        # Check for FP (other regions detected)
                        if filtered_regions:
                            fp_total += 1
                        else:
                            tn_total += 1
                
                # Calculate metrics for this group
                total_samples = len(group_df)
                sensitivity = tp_total / (tp_total + fn_total) if (tp_total + fn_total) > 0 else 0.0
                specificity = tn_total / (tn_total + fp_total) if (tn_total + fp_total) > 0 else 0.0
                ppv = tp_total / (tp_total + fp_total) if (tp_total + fp_total) > 0 else 0.0
                npv = tn_total / (tn_total + fn_total) if (tn_total + fn_total) > 0 else 0.0
                
                results.append({
                    'threshold': threshold,
                    'method': method_name,
                    'ff': ff,
                    'deletion_length_mb': del_length_mb,
                    'tp': tp_total,
                    'fp': fp_total,
                    'fn': fn_total,
                    'tn': tn_total,
                    'sensitivity': sensitivity,
                    'specificity': specificity,
                    'ppv': ppv,
                    'npv': npv,
                    'total_samples': total_samples
                })
            
            return results
    
    # Original logic: calculate for all samples together
    # Get samples with valid data
    valid_samples = df[df[zscore_col].notna()].copy()
    
    # Also get all samples (including those without z-score) for FN counting
    all_samples = df.copy()
    
    if len(valid_samples) == 0 and len(all_samples) == 0:
        return {
            'threshold': threshold,
            'method': method_name,
            'tp': 0, 'fp': 0, 'fn': 0, 'tn': 0,
            'sensitivity': 0.0, 'specificity': 0.0, 'ppv': 0.0, 'npv': 0.0,
            'total_samples': 0
        }
    
    tp_total = 0
    fp_total = 0
    fn_total = 0
    tn_total = 0
    
    # Process samples with z-score
    for _, row in valid_samples.iterrows():
        sample_name = row['sample_name']
        sample_dir = analysis_dir / sample_name
        
        expected_deletion = {
            'chromosome': row['expected_deletion_chr'],
            'start': row['expected_deletion_start'],
            'end': row['expected_deletion_end']
        }
        
        metrics = calculate_metrics_for_sample(
            sample_dir, sample_name, expected_deletion,
            method, output_type, threshold, min_length
        )
        
        tp_total += metrics['tp']
        fp_total += metrics['fp']
        fn_total += metrics['fn']
        tn_total += metrics['tn']
    
    # Process samples without z-score (they are FN)
    samples_without_zscore = all_samples[all_samples[zscore_col].isna()]
    for _, row in samples_without_zscore.iterrows():
        sample_name = row['sample_name']
        sample_dir = analysis_dir / sample_name
        
        # Check if there are any FP regions
        if method == 'WCX':
            bed_file = sample_dir / f"Output_WCX" / output_type / f"{sample_name}.wcx.{output_type}_aberrations.bed"
            detected_regions = parse_wcx_bed_all_regions(bed_file)
        else:
            report_file = sample_dir / f"Output_WC" / output_type / f"{sample_name}.wc.{output_type}.report.txt"
            detected_regions = parse_wc_report_all_regions(report_file)
        
        # Filter by threshold and min_length
        filtered_regions = [r for r in detected_regions if r['zscore'] <= threshold]
        if min_length is not None:
            filtered_regions = [r for r in filtered_regions if (r['end'] - r['start']) >= min_length]
        
        # No target detection (FN)
        fn_total += 1
        
        # Check for FP (other regions detected)
        if filtered_regions:
            fp_total += 1
        else:
            tn_total += 1
    
    # Calculate performance metrics
    total_samples = len(all_samples)
    sensitivity = tp_total / (tp_total + fn_total) if (tp_total + fn_total) > 0 else 0.0
    specificity = tn_total / (tn_total + fp_total) if (tn_total + fp_total) > 0 else 0.0
    ppv = tp_total / (tp_total + fp_total) if (tp_total + fp_total) > 0 else 0.0
    npv = tn_total / (tn_total + fn_total) if (tn_total + fn_total) > 0 else 0.0
    
    return {
        'threshold': threshold,
        'method': method_name,
        'tp': tp_total,
        'fp': fp_total,
        'fn': fn_total,
        'tn': tn_total,
        'sensitivity': sensitivity,
        'specificity': specificity,
        'ppv': ppv,
        'npv': npv,
        'total_samples': total_samples
    }


def main():
    parser = argparse.ArgumentParser(
        description="Calculate performance metrics by z-score threshold"
    )
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help='Input zscore_extraction.tsv file'
    )
    parser.add_argument(
        '-d', '--sample-dir',
        type=str,
        required=True,
        help='Directory containing sample directories (e.g., /data/md_validation/1p36 or /home/ken/ken-nipt/analysis/md_validation/1p36)'
    )
    parser.add_argument(
        '-o', '--output-dir',
        type=str,
        required=True,
        help='Output directory for results'
    )
    parser.add_argument(
        '-t', '--thresholds',
        type=str,
        default=None,
        help='Comma-separated list of thresholds (if not provided, uses thresholds from analyze_sensitivity_by_zscore output)'
    )
    parser.add_argument(
        '--method',
        type=str,
        choices=['WC_orig', 'WC_fetus', 'WCX_orig', 'WCX_fetus', 'all'],
        default='all',
        help='Method to analyze (default: all)'
    )
    parser.add_argument(
        '--min-length',
        type=int,
        default=None,
        help='Minimum length (bp) for a detected region to be considered as FP. Regions shorter than this will be ignored for FP counting. (default: None, all regions counted)'
    )
    parser.add_argument(
        '--group-by-ff-length',
        action='store_true',
        help='Group results by FF and deletion_length_mb. Output will include separate metrics for each FF/length combination.'
    )
    parser.add_argument(
        '--ppv',
        type=float,
        default=None,
        help='Minimum PPV threshold (0-100). Only calculate and output results where PPV >= this value. This significantly reduces computation when you only care about high PPV thresholds. (e.g., --ppv 90 for PPV >= 90%%)'
    )
    
    args = parser.parse_args()
    
    input_file = Path(args.input)
    analysis_dir = Path(args.sample_dir)
    output_dir = Path(args.output_dir)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        return
    
    if not analysis_dir.exists():
        logger.error(f"Sample directory not found: {analysis_dir}")
        return
    
    if not analysis_dir.is_dir():
        logger.error(f"Sample path is not a directory: {analysis_dir}")
        return
    
    # Read z-score extraction results
    logger.info(f"Reading z-score data from {input_file}")
    df = pd.read_csv(input_file, sep='\t')
    
    logger.info(f"Found {len(df)} samples")
    
    # Determine thresholds
    if args.thresholds:
        thresholds = [float(t.strip()) for t in args.thresholds.split(',')]
    else:
        # Use a range of thresholds based on z-score distribution
        methods_to_check = ['WCX_orig', 'WCX_fetus', 'WC_orig', 'WC_fetus']
        all_thresholds = set()
        
        # If PPV filter is specified, use two-stage sampling for efficiency
        if args.ppv is not None:
            min_ppv = args.ppv / 100.0
            logger.info(f"Using two-stage sampling to find thresholds with PPV >= {min_ppv:.0%}")
            
            # Find global min/max z-score across all methods
            global_min_z = None
            global_max_z = None
            for method in methods_to_check:
                zscore_col = method + '_zscore'
                valid_scores = df[df[zscore_col].notna()][zscore_col].values
                if len(valid_scores) > 0:
                    method_min = float(np.min(valid_scores))
                    method_max = float(np.max(valid_scores))
                    if global_min_z is None or method_min < global_min_z:
                        global_min_z = method_min
                    if global_max_z is None or method_max > global_max_z:
                        global_max_z = method_max
            
            if global_min_z is not None and global_max_z is not None:
                # Stage 1: Coarse sampling (step 5.0) to find PPV >= min_ppv range
                coarse_step = 5.0
                coarse_thresholds = np.arange(global_min_z, global_max_z + coarse_step, coarse_step)
                
                # Find the range where PPV >= min_ppv for each method
                # Use the union of all ranges to be safe
                ppv_ranges = []
                
                for method in methods_to_check:
                    zscore_col = method + '_zscore'
                    if zscore_col not in df.columns or not df[zscore_col].notna().any():
                        continue
                    
                    # Parse method name
                    if method.startswith('WCX_'):
                        method_type = 'WCX'
                        output_type = method.split('_')[1]
                    else:
                        method_type = 'WC'
                        output_type = method.split('_')[1]
                    
                    method_range_start = None
                    method_range_end = None
                    
                    # Quick check on coarse thresholds to find PPV range
                    for coarse_thresh in coarse_thresholds:
                        try:
                            temp_metrics = calculate_metrics_at_threshold(
                                df, analysis_dir, method_type, output_type, coarse_thresh,
                                args.min_length, group_by_ff_length=False
                            )
                            if temp_metrics['ppv'] >= min_ppv:
                                if method_range_start is None:
                                    method_range_start = coarse_thresh
                                method_range_end = coarse_thresh
                        except Exception as e:
                            logger.debug(f"Error calculating metrics for {method} at {coarse_thresh}: {e}")
                            continue
                    
                    if method_range_start is not None:
                        ppv_ranges.append((method_range_start, method_range_end))
                        logger.info(f"{method}: Found PPV >= {min_ppv:.0%} range: {method_range_start:.1f} to {method_range_end:.1f}")
                
                if ppv_ranges:
                    # Use the union of all ranges (min start, max end) with safety margin
                    overall_start = min(r[0] for r in ppv_ranges)
                    overall_end = max(r[1] for r in ppv_ranges)
                    
                    # Expand range for safety (add 3.0 on each side to catch edge cases)
                    safety_margin = 3.0
                    fine_start = max(global_min_z, overall_start - safety_margin)
                    fine_end = min(global_max_z, overall_end + safety_margin)
                    
                    # Stage 2: Fine sampling (step 0.5) within the PPV range
                    fine_step = 0.5
                    fine_thresholds = np.arange(fine_start, fine_end + fine_step, fine_step)
                    all_thresholds.update(fine_thresholds)
                    
                    logger.info(f"Overall PPV >= {min_ppv:.0%} range: {overall_start:.1f} to {overall_end:.1f}")
                    logger.info(f"Calculating {len(fine_thresholds)} thresholds in fine-grained range [{fine_start:.1f}, {fine_end:.1f}]")
                    logger.info(f"Reduced from ~{len(coarse_thresholds) * 10} thresholds to {len(fine_thresholds)} thresholds")
                else:
                    # No threshold meets PPV requirement
                    logger.warning(f"No threshold found with PPV >= {min_ppv:.0%} in coarse scan. Using normal sampling.")
                    for method in methods_to_check:
                        zscore_col = method + '_zscore'
                        valid_scores = df[df[zscore_col].notna()][zscore_col].values
                        if len(valid_scores) > 0:
                            min_z = float(np.min(valid_scores))
                            max_z = float(np.max(valid_scores))
                            thresholds_range = np.arange(min_z, max_z + 1, 1.0)
                            all_thresholds.update(thresholds_range)
            else:
                # Fallback to normal sampling
                logger.warning("Could not determine z-score range. Using normal sampling.")
                for method in methods_to_check:
                    zscore_col = method + '_zscore'
                    valid_scores = df[df[zscore_col].notna()][zscore_col].values
                    if len(valid_scores) > 0:
                        min_z = float(np.min(valid_scores))
                        max_z = float(np.max(valid_scores))
                        thresholds_range = np.arange(min_z, max_z + 1, 1.0)
                        all_thresholds.update(thresholds_range)
        else:
            # Normal sampling without PPV filter
            for method in methods_to_check:
                zscore_col = method + '_zscore'
                valid_scores = df[df[zscore_col].notna()][zscore_col].values
                if len(valid_scores) > 0:
                    min_z = float(np.min(valid_scores))
                    max_z = float(np.max(valid_scores))
                    thresholds_range = np.arange(min_z, max_z + 1, 1.0)
                    all_thresholds.update(thresholds_range)
        
        thresholds = sorted(all_thresholds)
        logger.info(f"Generated {len(thresholds)} thresholds from z-score distribution")
    
    # Determine methods to analyze
    if args.method == 'all':
        methods_to_analyze = [
            ('WC', 'orig'), ('WC', 'fetus'),
            ('WCX', 'orig'), ('WCX', 'fetus')
        ]
        # Also add OR combinations
        methods_to_analyze_or = [
            ('ALL_OR', [('WC', 'orig'), ('WC', 'fetus'), ('WCX', 'orig'), ('WCX', 'fetus')]),
            ('ORIG_OR', [('WC', 'orig'), ('WCX', 'orig')]),
            ('FETUS_OR', [('WC', 'fetus'), ('WCX', 'fetus')])
        ]
    else:
        method_map = {
            'WC_orig': ('WC', 'orig'),
            'WC_fetus': ('WC', 'fetus'),
            'WCX_orig': ('WCX', 'orig'),
            'WCX_fetus': ('WCX', 'fetus')
        }
        methods_to_analyze = [method_map[args.method]]
        methods_to_analyze_or = []
    
    # Calculate metrics for each threshold and method
    all_results = []
    
    print("\n" + "="*80)
    if args.group_by_ff_length:
        print("PERFORMANCE METRICS BY Z-SCORE THRESHOLD (GROUPED BY FF AND DELETION LENGTH)")
    else:
        print("PERFORMANCE METRICS BY Z-SCORE THRESHOLD")
    print("="*80)
    
    for method, output_type in methods_to_analyze:
        method_name = f'{method}_{output_type}'
        
        if args.group_by_ff_length:
            # Grouped output
            print(f"\n{method_name}:")
            print("=" * 80)
            
            # Get unique FF and deletion_length combinations
            df['ff'] = pd.to_numeric(df['ff'], errors='coerce')
            df['deletion_length_mb'] = pd.to_numeric(df['deletion_length_mb'], errors='coerce')
            valid_df = df[df[['ff', 'deletion_length_mb']].notna().all(axis=1)]
            
            if len(valid_df) == 0:
                print("  No valid FF/deletion_length data found")
                continue
            
            ff_values = sorted(valid_df['ff'].unique())
            length_values = sorted(valid_df['deletion_length_mb'].unique())
            
            # Calculate metrics for each threshold, grouped by FF/length
            for threshold in thresholds:
                grouped_results = calculate_metrics_at_threshold(
                    df, analysis_dir, method, output_type, threshold, 
                    args.min_length, group_by_ff_length=True
                )
                
                if isinstance(grouped_results, list):
                    # Filter by PPV if specified
                    if args.ppv is not None:
                        min_ppv = args.ppv / 100.0
                        filtered_results = [r for r in grouped_results if r['ppv'] >= min_ppv]
                        if len(filtered_results) == 0:
                            continue  # Skip this threshold if no group meets PPV requirement
                        grouped_results = filtered_results
                    
                    all_results.extend(grouped_results)
                    
                    # Print summary for this threshold
                    print(f"\nThreshold: {threshold:.2f}")
                    print("-" * 80)
                    print(f"{'FF':<8} {'DelLen(Mb)':<12} {'TP':<6} {'FP':<6} {'FN':<6} {'TN':<6} {'Sens':<10} {'Spec':<10} {'PPV':<10} {'NPV':<10} {'N':<6}")
                    print("-" * 80)
                    
                    for result in grouped_results:
                        print(f"{result['ff']:>6.0f}  {result['deletion_length_mb']:>10.0f}  "
                              f"{result['tp']:>4}  {result['fp']:>4}  {result['fn']:>4}  {result['tn']:>4}  "
                              f"{result['sensitivity']:>8.2%}  {result['specificity']:>8.2%}  "
                              f"{result['ppv']:>8.2%}  {result['npv']:>8.2%}  {result['total_samples']:>4}")
        else:
            # Original output (all samples together)
            print(f"\n{method_name}:")
            print("-" * 80)
            print(f"{'Threshold':<12} {'TP':<8} {'FP':<8} {'FN':<8} {'TN':<8} {'Sens':<10} {'Spec':<10} {'PPV':<10} {'NPV':<10}")
            print("-" * 80)
            
            for threshold in thresholds:
                metrics = calculate_metrics_at_threshold(
                    df, analysis_dir, method, output_type, threshold, 
                    args.min_length, group_by_ff_length=False
                )
                
                # Filter by PPV if specified
                if args.ppv is not None:
                    min_ppv = args.ppv / 100.0
                    if metrics['ppv'] < min_ppv:
                        continue  # Skip this threshold if PPV doesn't meet requirement
                
                all_results.append(metrics)
                
                # Print every 5th threshold or all if less than 20
                idx = thresholds.index(threshold)
                if len(thresholds) <= 20 or idx % max(1, len(thresholds) // 20) == 0 or idx == len(thresholds) - 1:
                    print(f"{threshold:>10.2f}  {metrics['tp']:>6}  {metrics['fp']:>6}  "
                          f"{metrics['fn']:>6}  {metrics['tn']:>6}  {metrics['sensitivity']:>8.2%}  "
                          f"{metrics['specificity']:>8.2%}  {metrics['ppv']:>8.2%}  {metrics['npv']:>8.2%}")
    
    # Calculate OR combinations
    if args.method == 'all' and methods_to_analyze_or:
        print("\n" + "="*80)
        print("OR COMBINATIONS")
        print("="*80)
        
        for or_name, or_methods in methods_to_analyze_or:
            print(f"\n{or_name}:")
            print("-" * 80)
            
            if args.group_by_ff_length:
                print(f"{'FF':<8} {'DelLen(Mb)':<12} {'TP':<6} {'FP':<6} {'FN':<6} {'TN':<6} {'Sens':<10} {'Spec':<10} {'PPV':<10} {'NPV':<10} {'N':<6}")
                print("-" * 80)
                
                for threshold in thresholds:
                    grouped_results = calculate_metrics_at_threshold_or(
                        df, analysis_dir, or_methods, threshold,
                        args.min_length, group_by_ff_length=True
                    )
                    
                    if isinstance(grouped_results, list):
                        # Filter by PPV if specified
                        if args.ppv is not None:
                            min_ppv = args.ppv / 100.0
                            filtered_results = [r for r in grouped_results if r['ppv'] >= min_ppv]
                            if len(filtered_results) == 0:
                                continue
                            grouped_results = filtered_results
                        
                        all_results.extend(grouped_results)
                        
                        # Print summary for this threshold
                        print(f"\nThreshold: {threshold:.2f}")
                        print("-" * 80)
                        for result in grouped_results:
                            print(f"{result['ff']:>6.0f}  {result['deletion_length_mb']:>10.0f}  "
                                  f"{result['tp']:>4}  {result['fp']:>4}  {result['fn']:>4}  {result['tn']:>4}  "
                                  f"{result['sensitivity']:>8.2%}  {result['specificity']:>8.2%}  "
                                  f"{result['ppv']:>8.2%}  {result['npv']:>8.2%}  {result['total_samples']:>4}")
            else:
                print(f"{'Threshold':<12} {'TP':<8} {'FP':<8} {'FN':<8} {'TN':<8} {'Sens':<10} {'Spec':<10} {'PPV':<10} {'NPV':<10}")
                print("-" * 80)
                
                for threshold in thresholds:
                    metrics = calculate_metrics_at_threshold_or(
                        df, analysis_dir, or_methods, threshold,
                        args.min_length, group_by_ff_length=False
                    )
                    
                    # Filter by PPV if specified
                    if args.ppv is not None:
                        min_ppv = args.ppv / 100.0
                        if metrics['ppv'] < min_ppv:
                            continue
                    
                    all_results.append(metrics)
                    
                    # Print every 5th threshold or all if less than 20
                    idx = thresholds.index(threshold)
                    if len(thresholds) <= 20 or idx % max(1, len(thresholds) // 20) == 0 or idx == len(thresholds) - 1:
                        print(f"{threshold:>10.2f}  {metrics['tp']:>6}  {metrics['fp']:>6}  "
                              f"{metrics['fn']:>6}  {metrics['tn']:>6}  {metrics['sensitivity']:>8.2%}  "
                              f"{metrics['specificity']:>8.2%}  {metrics['ppv']:>8.2%}  {metrics['npv']:>8.2%}")
    
    # Save results
    if all_results:
        results_df = pd.DataFrame(all_results)
        if args.group_by_ff_length:
            output_file = output_dir / 'performance_metrics_by_threshold_grouped.tsv'
        else:
            output_file = output_dir / 'performance_metrics_by_threshold.tsv'
        results_df.to_csv(output_file, sep='\t', index=False, float_format='%.2f')
        logger.info(f"Results saved to {output_file}")
        
        # Find optimal threshold (maximize sensitivity while keeping PPV high)
        if not args.group_by_ff_length:
            print("\n" + "="*80)
            print("OPTIMAL THRESHOLD ANALYSIS")
            print("="*80)
            
            for method_name in results_df['method'].unique():
                method_df = results_df[results_df['method'] == method_name].copy()
                method_df = method_df.sort_values('threshold')
                
                print(f"\n{method_name}:")
                print("-" * 80)
                
                # Find thresholds with high sensitivity and reasonable PPV
                for min_ppv in [0.9, 0.8, 0.7, 0.5]:
                    filtered = method_df[method_df['ppv'] >= min_ppv]
                    if len(filtered) > 0:
                        best = filtered.loc[filtered['sensitivity'].idxmax()]
                        print(f"  PPV >= {min_ppv:.0%}: threshold = {best['threshold']:.2f}, "
                              f"Sensitivity = {best['sensitivity']:.2%}, PPV = {best['ppv']:.2%}")
        else:
            # For grouped results, show summary by FF and length
            print("\n" + "="*80)
            print("SUMMARY BY FF AND DELETION LENGTH")
            print("="*80)
            
            # Group by FF and deletion_length_mb and show average metrics
            if 'ff' in results_df.columns and 'deletion_length_mb' in results_df.columns:
                summary_cols = ['method', 'ff', 'deletion_length_mb', 'sensitivity', 'specificity', 'ppv', 'npv', 'total_samples']
                available_cols = [col for col in summary_cols if col in results_df.columns]
                summary_df = results_df[available_cols].groupby(['method', 'ff', 'deletion_length_mb']).agg({
                    'sensitivity': 'mean',
                    'specificity': 'mean',
                    'ppv': 'mean',
                    'npv': 'mean',
                    'total_samples': 'first'
                }).reset_index()
                
                for method_name in summary_df['method'].unique():
                    method_summary = summary_df[summary_df['method'] == method_name]
                    print(f"\n{method_name} (Average across thresholds):")
                    print("-" * 80)
                    print(f"{'FF':<8} {'DelLen(Mb)':<12} {'Sens':<10} {'Spec':<10} {'PPV':<10} {'NPV':<10} {'N':<6}")
                    print("-" * 80)
                    for _, row in method_summary.iterrows():
                        print(f"{row['ff']:>6.0f}  {row['deletion_length_mb']:>10.0f}  "
                              f"{row['sensitivity']:>8.2%}  {row['specificity']:>8.2%}  "
                              f"{row['ppv']:>8.2%}  {row['npv']:>8.2%}  {row['total_samples']:>4.0f}")
    
    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)
    logger.info(f"Results saved to {output_dir}")


if __name__ == '__main__':
    main()

