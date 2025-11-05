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
                                  output_type: str, threshold: float) -> dict:
    """Calculate TP, FP, FN, TN for a single sample at given threshold
    
    Logic:
    - TP: Target 영역이 검출됨
    - FN: Target 영역이 검출 안됨
    - FP: 다른 영역이 검출됨 (Target 영역이 아닌 다른 영역)
    - TN: 다른 영역이 검출 안됨
    
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
        # Check if this region overlaps with target deletion
        if region['chr'] == expected_chr:
            if check_overlap(expected_start, expected_end, region['start'], region['end']):
                target_detected = True
                continue
        
        # This is a different region (not target)
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


def calculate_metrics_at_threshold(df: pd.DataFrame, root_dir: Path, work_dir: str,
                                   method: str, output_type: str, threshold: float) -> dict:
    """Calculate performance metrics at a given threshold"""
    
    method_name = f'{method}_{output_type}'
    zscore_col = method_name + '_zscore'
    
    # Get samples with valid data
    valid_samples = df[df[zscore_col].notna()].copy()
    
    if len(valid_samples) == 0:
        return {
            'threshold': threshold,
            'method': method_name,
            'tp': 0, 'fp': 0, 'fn': 0, 'tn': 0,
            'sensitivity': 0.0, 'specificity': 0.0, 'ppv': 0.0, 'npv': 0.0
        }
    
    analysis_dir = root_dir / "analysis" / work_dir
    
    tp_total = 0
    fp_total = 0
    fn_total = 0
    tn_total = 0
    
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
            method, output_type, threshold
        )
        
        tp_total += metrics['tp']
        fp_total += metrics['fp']
        fn_total += metrics['fn']
        tn_total += metrics['tn']
    
    # Calculate performance metrics
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
        'total_samples': len(valid_samples)
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
        '-r', '--root-dir',
        type=str,
        required=True,
        help='Root directory (e.g., /home/ken/ken-nipt)'
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
    
    args = parser.parse_args()
    
    input_file = Path(args.input)
    root_dir = Path(args.root_dir)
    work_dir = args.work_dir
    output_dir = Path(args.output_dir)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
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
        
        for method in methods_to_check:
            zscore_col = method + '_zscore'
            valid_scores = df[df[zscore_col].notna()][zscore_col].values
            if len(valid_scores) > 0:
                min_z = float(np.min(valid_scores))
                max_z = float(np.max(valid_scores))
                # Create thresholds from min to max with step of 1.0
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
    else:
        method_map = {
            'WC_orig': ('WC', 'orig'),
            'WC_fetus': ('WC', 'fetus'),
            'WCX_orig': ('WCX', 'orig'),
            'WCX_fetus': ('WCX', 'fetus')
        }
        methods_to_analyze = [method_map[args.method]]
    
    # Calculate metrics for each threshold and method
    all_results = []
    
    print("\n" + "="*80)
    print("PERFORMANCE METRICS BY Z-SCORE THRESHOLD")
    print("="*80)
    
    for method, output_type in methods_to_analyze:
        method_name = f'{method}_{output_type}'
        print(f"\n{method_name}:")
        print("-" * 80)
        print(f"{'Threshold':<12} {'TP':<8} {'FP':<8} {'FN':<8} {'TN':<8} {'Sens':<10} {'Spec':<10} {'PPV':<10} {'NPV':<10}")
        print("-" * 80)
        
        for threshold in thresholds:
            metrics = calculate_metrics_at_threshold(
                df, root_dir, work_dir, method, output_type, threshold
            )
            
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
        output_file = output_dir / 'performance_metrics_by_threshold.tsv'
        results_df.to_csv(output_file, sep='\t', index=False, float_format='%.2f')
        logger.info(f"Results saved to {output_file}")
        
        # Find optimal threshold (maximize sensitivity while keeping PPV high)
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
    
    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)
    logger.info(f"Results saved to {output_dir}")


if __name__ == '__main__':
    main()

