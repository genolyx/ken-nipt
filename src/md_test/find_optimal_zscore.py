#!/usr/bin/env python3
"""
Find Optimal Z-Score Threshold for MD Detection

This script:
1. Evaluates different z-score thresholds (e.g., 3.0, 3.5, 4.0, ..., 7.0)
2. For each threshold, calculates TP, FP, FN
3. Computes Sensitivity, PPV, F1 score
4. Identifies optimal threshold based on metrics
"""

import json
import csv
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def parse_deletion_from_json(json_file: Path) -> Optional[Dict]:
    """Parse deletion information from JSON file"""
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        deletion = data.get('deletion', {})
        if not deletion:
            return None
        
        return {
            'chromosome': deletion.get('chromosome'),
            'start': deletion.get('start'),
            'end': deletion.get('end'),
            'size_bp': deletion.get('size_bp')
        }
    except Exception as e:
        logger.warning(f"Failed to parse JSON {json_file}: {e}")
        return None


def check_overlap(deletion: Dict, detected_start: int, detected_end: int) -> bool:
    """Check if detected region overlaps with expected deletion"""
    # Calculate overlap
    overlap_start = max(deletion['start'], detected_start)
    overlap_end = min(deletion['end'], detected_end)
    
    # Check if there's meaningful overlap (at least 50% of deletion or detected region)
    if overlap_start < overlap_end:
        overlap_size = overlap_end - overlap_start
        deletion_size = deletion['end'] - deletion['start']
        detected_size = detected_end - detected_start
        
        # At least 50% overlap with either region
        overlap_ratio_deletion = overlap_size / deletion_size if deletion_size > 0 else 0
        overlap_ratio_detected = overlap_size / detected_size if detected_size > 0 else 0
        
        return overlap_ratio_deletion >= 0.5 or overlap_ratio_detected >= 0.5
    
    return False


def parse_wcx_bed_file(bed_file: Path, zscore_threshold: float) -> List[Dict]:
    """
    Parse WCX aberrations.bed file and extract regions with z-score >= threshold
    
    Returns list of detected regions: [{'chr': str, 'start': int, 'end': int, 'zscore': float}, ...]
    """
    detected_regions = []
    
    if not bed_file.exists():
        return detected_regions
    
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                fields = line.split('\t')
                if len(fields) < 5:
                    continue
                
                chr_name = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                
                # Try to find z-score column (usually column 4 or 5)
                # Format: chr start end effect zscore ...
                zscore = None
                for i in range(3, min(len(fields), 6)):
                    try:
                        zscore_val = float(fields[i])
                        # z-score is typically negative for deletions
                        if abs(zscore_val) >= zscore_threshold:
                            zscore = abs(zscore_val)
                            break
                    except ValueError:
                        continue
                
                if zscore is not None:
                    detected_regions.append({
                        'chr': chr_name,
                        'start': start,
                        'end': end,
                        'zscore': zscore
                    })
    except Exception as e:
        logger.warning(f"Failed to parse WCX bed file {bed_file}: {e}")
    
    return detected_regions


def parse_wc_report_file(report_file: Path, zscore_threshold: float) -> List[Dict]:
    """
    Parse WC report.txt file and extract regions with z-score >= threshold
    
    Looks for "# Test results:" section at the bottom of the file.
    Format: z-score effect  mbsize  location
    Example: -51.74  -29.66  9.20    1:800000-10000000
    
    Returns list of detected regions: [{'chr': str, 'start': int, 'end': int, 'zscore': float}, ...]
    """
    detected_regions = []
    
    if not report_file.exists():
        return detected_regions
    
    try:
        with open(report_file, 'r') as f:
            lines = f.readlines()
        
        # Find "# Test results:" section
        in_test_section = False
        for i, line in enumerate(lines):
            if '# Test results:' in line or '# Test results' in line:
                in_test_section = True
                continue
            
            if in_test_section:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Skip header line if present
                if 'z-score' in line.lower() and 'effect' in line.lower():
                    continue
                
                # Parse line: z-score effect  mbsize  location
                # Example: -51.74  -29.66  9.20    1:800000-10000000
                fields = line.split()
                if len(fields) < 4:
                    continue
                
                try:
                    zscore = abs(float(fields[0]))  # z-score (absolute value)
                    location = fields[3]  # location format: chr:start-end
                    
                    if zscore < zscore_threshold:
                        continue
                    
                    # Parse location: chr:start-end
                    if ':' in location:
                        chr_part, coord_part = location.split(':', 1)
                        if '-' in coord_part:
                            start_str, end_str = coord_part.split('-', 1)
                            chr_name = chr_part
                            start = int(start_str)
                            end = int(end_str)
                            
                            detected_regions.append({
                                'chr': chr_name,
                                'start': start,
                                'end': end,
                                'zscore': zscore
                            })
                except (ValueError, IndexError) as e:
                    logger.debug(f"Failed to parse line in WC report: {line[:50]}... Error: {e}")
                    continue
    except Exception as e:
        logger.warning(f"Failed to parse WC report file {report_file}: {e}")
    
    return detected_regions


def evaluate_threshold(
    sample_dir: Path,
    sample_name: str,
    deletion: Dict,
    zscore_threshold: float,
    method: str,
    output_type: str
) -> Dict:
    """
    Evaluate detection for a single sample with given z-score threshold
    
    Args:
        sample_dir: Sample directory path
        sample_name: Sample name
        deletion: Deletion dictionary from JSON
        zscore_threshold: Z-score threshold to test
        method: 'WC' or 'WCX'
        output_type: 'orig' or 'fetus'
    
    Returns:
        Dict with 'tp', 'fp', 'fn' counts
    """
    result = {'tp': 0, 'fp': 0, 'fn': 0}
    
    # Get detected regions based on method
    if method == 'WCX':
        bed_file = sample_dir / f"Output_WCX" / output_type / f"{sample_name}.wcx.{output_type}_aberrations.bed"
        detected_regions = parse_wcx_bed_file(bed_file, zscore_threshold)
    else:  # WC
        report_file = sample_dir / f"Output_WC" / output_type / f"{sample_name}.wc.{output_type}.report.txt"
        detected_regions = parse_wc_report_file(report_file, zscore_threshold)
    
    if not detected_regions:
        # No detection: False Negative
        result['fn'] = 1
        return result
    
    # Check if any detected region overlaps with expected deletion
    has_tp = False
    has_fp = False
    
    for region in detected_regions:
        # Check chromosome match
        if region['chr'] != deletion['chromosome']:
            has_fp = True
            continue
        
        # Check overlap
        if check_overlap(deletion, region['start'], region['end']):
            has_tp = True
        else:
            has_fp = True
    
    if has_tp:
        result['tp'] = 1
    elif has_fp:
        result['fp'] = 1
    else:
        result['fn'] = 1
    
    return result


def scan_samples(root_dir: Path, work_dir: str, sample_list: Optional[List[str]] = None) -> List[Dict]:
    """Scan directories for samples"""
    samples = []
    analysis_dir = root_dir / "analysis" / work_dir
    
    if not analysis_dir.exists():
        logger.error(f"Analysis directory not found: {analysis_dir}")
        return samples
    
    # Get list of sample directories
    sample_dirs = [d for d in analysis_dir.iterdir() if d.is_dir()]
    
    if sample_list:
        sample_dirs = [d for d in sample_dirs if d.name in sample_list]
    
    for sample_dir in sample_dirs:
        sample_name = sample_dir.name
        json_file = sample_dir / f"{sample_name}.json"
        
        if not json_file.exists():
            logger.debug(f"No JSON file found for {sample_name}")
            continue
        
        deletion = parse_deletion_from_json(json_file)
        if not deletion:
            logger.debug(f"No deletion found in JSON for {sample_name}")
            continue
        
        samples.append({
            'sample_name': sample_name,
            'sample_dir': sample_dir,
            'deletion': deletion
        })
    
    logger.info(f"Found {len(samples)} samples with deletion information")
    return samples


def evaluate_all_thresholds(
    samples: List[Dict],
    thresholds: List[float],
    methods: List[str] = ['WC_orig', 'WC_fetus', 'WCX_orig', 'WCX_fetus']
) -> pd.DataFrame:
    """
    Evaluate all thresholds for all samples
    
    Returns DataFrame with columns:
    method, threshold, tp, fp, fn, sensitivity, ppv, f1_score
    """
    results = []
    
    for threshold in thresholds:
        logger.info(f"Evaluating threshold: {threshold}")
        
        for method_spec in methods:
            method, output_type = method_spec.split('_')
            
            tp_total = 0
            fp_total = 0
            fn_total = 0
            
            for sample in samples:
                result = evaluate_threshold(
                    sample['sample_dir'],
                    sample['sample_name'],
                    sample['deletion'],
                    threshold,
                    method,
                    output_type
                )
                
                tp_total += result['tp']
                fp_total += result['fp']
                fn_total += result['fn']
            
            # Calculate metrics
            sensitivity = tp_total / (tp_total + fn_total) if (tp_total + fn_total) > 0 else 0.0
            ppv = tp_total / (tp_total + fp_total) if (tp_total + fp_total) > 0 else 0.0
            f1_score = 2 * (sensitivity * ppv) / (sensitivity + ppv) if (sensitivity + ppv) > 0 else 0.0
            
            results.append({
                'method': method_spec,
                'threshold': threshold,
                'tp': tp_total,
                'fp': fp_total,
                'fn': fn_total,
                'sensitivity': sensitivity,
                'ppv': ppv,
                'f1_score': f1_score,
                'total_samples': len(samples)
            })
    
    return pd.DataFrame(results)


def find_optimal_threshold(df: pd.DataFrame, metric: str = 'f1_score') -> pd.DataFrame:
    """Find optimal threshold for each method based on specified metric"""
    optimal = []
    
    for method in df['method'].unique():
        method_df = df[df['method'] == method]
        
        if metric == 'f1_score':
            # Maximize F1 score
            best_row = method_df.loc[method_df[metric].idxmax()]
        elif metric == 'sensitivity':
            # Maximize sensitivity
            best_row = method_df.loc[method_df[metric].idxmax()]
        elif metric == 'ppv':
            # Maximize PPV
            best_row = method_df.loc[method_df[metric].idxmax()]
        else:
            # Balanced: maximize (sensitivity + ppv) / 2
            method_df['balanced'] = (method_df['sensitivity'] + method_df['ppv']) / 2
            best_row = method_df.loc[method_df['balanced'].idxmax()]
        
        optimal.append({
            'method': method,
            'optimal_threshold': best_row['threshold'],
            'sensitivity': best_row['sensitivity'],
            'ppv': best_row['ppv'],
            'f1_score': best_row['f1_score'],
            'tp': best_row['tp'],
            'fp': best_row['fp'],
            'fn': best_row['fn']
        })
    
    return pd.DataFrame(optimal)


def main():
    parser = argparse.ArgumentParser(
        description="Find optimal z-score threshold for MD detection"
    )
    parser.add_argument(
        '-r', '--root-dir',
        type=str,
        default='/home/ken/ken-nipt',
        help='Root directory path (default: /home/ken/ken-nipt)'
    )
    parser.add_argument(
        '-w', '--work-dir',
        type=str,
        required=True,
        help='Work directory (e.g., md_validation/1p36)'
    )
    parser.add_argument(
        '-o', '--output-dir',
        type=str,
        required=True,
        help='Output directory for results'
    )
    parser.add_argument(
        '--sample-sheet',
        type=str,
        help='Optional sample sheet TSV file to filter samples'
    )
    parser.add_argument(
        '--thresholds',
        type=str,
        default='3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0',
        help='Comma-separated list of z-score thresholds to test (default: 3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0)'
    )
    parser.add_argument(
        '--metric',
        type=str,
        choices=['f1_score', 'sensitivity', 'ppv', 'balanced'],
        default='f1_score',
        help='Metric to optimize (default: f1_score)'
    )
    
    args = parser.parse_args()
    
    # Parse thresholds
    thresholds = [float(t.strip()) for t in args.thresholds.split(',')]
    logger.info(f"Testing thresholds: {thresholds}")
    
    # Parse sample list if provided
    sample_list = None
    if args.sample_sheet:
        sample_list = []
        with open(args.sample_sheet, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                sample_list.append(row['SAMPLE_NAME'])
        logger.info(f"Filtering to {len(sample_list)} samples from sample sheet")
    
    # Scan samples
    root_dir = Path(args.root_dir)
    samples = scan_samples(root_dir, args.work_dir, sample_list)
    
    if len(samples) == 0:
        logger.error("No samples found!")
        return
    
    # Evaluate all thresholds
    logger.info("Evaluating all thresholds...")
    results_df = evaluate_all_thresholds(samples, thresholds)
    
    # Find optimal thresholds
    optimal_df = find_optimal_threshold(results_df, args.metric)
    
    # Save results
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results_file = output_dir / 'zscore_threshold_results.tsv'
    results_df.to_csv(results_file, sep='\t', index=False)
    logger.info(f"Saved detailed results to {results_file}")
    
    optimal_file = output_dir / 'optimal_zscore_thresholds.tsv'
    optimal_df.to_csv(optimal_file, sep='\t', index=False)
    logger.info(f"Saved optimal thresholds to {optimal_file}")
    
    # Print summary
    print("\n" + "="*80)
    print("OPTIMAL Z-SCORE THRESHOLDS")
    print("="*80)
    print(optimal_df.to_string(index=False))
    print("="*80)
    
    logger.info("Analysis complete!")


if __name__ == '__main__':
    main()

