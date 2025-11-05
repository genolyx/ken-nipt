#!/usr/bin/env python3
"""
Find optimal z-score threshold for MD detection

Reads zscore_extraction.tsv and analyzes sensitivity across different thresholds.
Groups results by FF and deletion length.
"""

import csv
import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def extract_sample_info(sample_name: str) -> Dict:
    """Extract FF, deletion_length, and reads from sample name
    
    Format: mom_idx_fetus_idx_1p36deletionsyndrome_FF{ff}_15M_{del}Mb_{gender}
    Example: 10_1_1p36deletionsyndrome_FF10_15M_10Mb_F
    """
    info = {
        'ff': None,
        'deletion_length': None,  # in Mb
        'reads': None,  # in millions
        'gender': None
    }
    
    # Extract FF
    ff_match = re.search(r'FF(\d+)', sample_name)
    if ff_match:
        info['ff'] = int(ff_match.group(1))
    
    # Extract deletion length (Mb)
    del_match = re.search(r'(\d+)Mb', sample_name)
    if del_match:
        info['deletion_length'] = int(del_match.group(1))
    
    # Extract reads (M)
    reads_match = re.search(r'(\d+)M(?![a-z])', sample_name)
    if reads_match:
        info['reads'] = int(reads_match.group(1))
    
    # Extract gender (last character should be F or M)
    parts = sample_name.split('_')
    if parts and parts[-1] in ['F', 'M']:
        info['gender'] = parts[-1]
    
    return info


def calculate_sensitivity_at_threshold(df: pd.DataFrame, method: str, threshold: float) -> Dict:
    """Calculate sensitivity at a given z-score threshold
    
    Args:
        df: DataFrame with z-score data
        method: One of 'WC_orig', 'WC_fetus', 'WCX_orig', 'WCX_fetus'
        threshold: Z-score threshold (negative for deletions, e.g., -6.0)
    
    Returns:
        Dictionary with sensitivity statistics
    """
    zscore_col = method
    
    # Filter rows with valid z-score
    valid_df = df[df[zscore_col].notna()].copy()
    
    if len(valid_df) == 0:
        return {
            'total': 0,
            'detected': 0,
            'not_detected': 0,
            'sensitivity': 0.0
        }
    
    # For deletions, z-score is negative. Detection means z-score <= threshold (more negative)
    # Example: threshold = -6.0, z-score = -51.74 -> detected (because -51.74 <= -6.0)
    valid_df['detected'] = valid_df[zscore_col] <= threshold
    
    detected = valid_df['detected'].sum()
    not_detected = (~valid_df['detected']).sum()
    total = len(valid_df)
    sensitivity = detected / total if total > 0 else 0.0
    
    return {
        'total': total,
        'detected': detected,
        'not_detected': not_detected,
        'sensitivity': sensitivity
    }


def analyze_by_group(df: pd.DataFrame, method: str, threshold: float, 
                      group_by: List[str]) -> pd.DataFrame:
    """Analyze sensitivity grouped by FF, deletion_length, etc.
    
    Args:
        df: DataFrame with z-score data
        method: One of 'WC_orig', 'WC_fetus', 'WCX_orig', 'WCX_fetus'
        threshold: Z-score threshold
        group_by: List of column names to group by (e.g., ['ff', 'deletion_length'])
    
    Returns:
        DataFrame with grouped statistics
    """
    zscore_col = method
    
    # Filter rows with valid z-score and complete grouping info
    group_cols = [col for col in group_by if col in df.columns]
    valid_df = df[df[zscore_col].notna() & df[group_cols].notna().all(axis=1)].copy()
    
    if len(valid_df) == 0:
        return pd.DataFrame()
    
    # Calculate detection status
    valid_df['detected'] = valid_df[zscore_col] <= threshold
    
    # Group and calculate statistics
    grouped = valid_df.groupby(group_cols).agg({
        'detected': ['sum', 'count']
    }).reset_index()
    
    grouped.columns = group_cols + ['detected_count', 'total_count']
    grouped['not_detected_count'] = grouped['total_count'] - grouped['detected_count']
    grouped['sensitivity'] = grouped['detected_count'] / grouped['total_count']
    grouped['threshold'] = threshold
    grouped['method'] = method
    
    return grouped


def main():
    parser = argparse.ArgumentParser(
        description="Find optimal z-score threshold for MD detection"
    )
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help='Input zscore_extraction.tsv file'
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
        default='-3.0,-4.0,-5.0,-6.0,-7.0,-8.0,-9.0,-10.0',
        help='Comma-separated list of z-score thresholds to test (default: -3.0,-4.0,-5.0,-6.0,-7.0,-8.0,-9.0,-10.0)'
    )
    
    args = parser.parse_args()
    
    input_file = Path(args.input)
    output_dir = Path(args.output_dir)
    thresholds = [float(t.strip()) for t in args.thresholds.split(',')]
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        return
    
    # Read z-score extraction results
    logger.info(f"Reading z-score data from {input_file}")
    df = pd.read_csv(input_file, sep='\t')
    
    logger.info(f"Found {len(df)} samples")
    
    # Extract FF, deletion_length, reads, gender from sample names
    sample_info_list = []
    for sample_name in df['sample_name']:
        info = extract_sample_info(sample_name)
        info['sample_name'] = sample_name
        sample_info_list.append(info)
    
    sample_info_df = pd.DataFrame(sample_info_list)
    
    # Merge with original dataframe
    df = df.merge(sample_info_df, on='sample_name', how='left')
    
    # Check data availability
    methods = ['WC_orig', 'WC_fetus', 'WCX_orig', 'WCX_fetus']
    print("\n" + "="*80)
    print("DATA AVAILABILITY")
    print("="*80)
    for method in methods:
        valid_count = df[df[method].notna()].shape[0]
        print(f"{method}: {valid_count}/{len(df)} samples with z-score data")
    
    # Current threshold analysis (-6.0)
    current_threshold = -6.0
    print("\n" + "="*80)
    print(f"ANALYSIS AT CURRENT THRESHOLD: {current_threshold}")
    print("="*80)
    
    for method in methods:
        print(f"\n{method}:")
        print("-" * 80)
        
        # Overall sensitivity
        overall = calculate_sensitivity_at_threshold(df, method, current_threshold)
        print(f"  Overall: {overall['detected']}/{overall['total']} detected "
              f"(Sensitivity: {overall['sensitivity']:.2%})")
        
        # By FF
        if df['ff'].notna().any():
            print(f"\n  By FF:")
            ff_groups = analyze_by_group(df, method, current_threshold, ['ff'])
            if len(ff_groups) > 0:
                for _, row in ff_groups.iterrows():
                    print(f"    FF {row['ff']:.0f}%: {row['detected_count']:.0f}/{row['total_count']:.0f} "
                          f"(Sensitivity: {row['sensitivity']:.2%})")
        
        # By deletion length
        if df['deletion_length'].notna().any():
            print(f"\n  By Deletion Length:")
            del_groups = analyze_by_group(df, method, current_threshold, ['deletion_length'])
            if len(del_groups) > 0:
                for _, row in del_groups.iterrows():
                    print(f"    {row['deletion_length']:.0f}Mb: {row['detected_count']:.0f}/{row['total_count']:.0f} "
                          f"(Sensitivity: {row['sensitivity']:.2%})")
        
        # By FF and deletion length
        if df['ff'].notna().any() and df['deletion_length'].notna().any():
            print(f"\n  By FF and Deletion Length:")
            combined_groups = analyze_by_group(df, method, current_threshold, ['ff', 'deletion_length'])
            if len(combined_groups) > 0:
                combined_groups = combined_groups.sort_values(['ff', 'deletion_length'])
                for _, row in combined_groups.iterrows():
                    print(f"    FF {row['ff']:.0f}%, {row['deletion_length']:.0f}Mb: "
                          f"{row['detected_count']:.0f}/{row['total_count']:.0f} "
                          f"(Sensitivity: {row['sensitivity']:.2%})")
    
    # Threshold optimization
    print("\n" + "="*80)
    print("THRESHOLD OPTIMIZATION")
    print("="*80)
    
    threshold_results = []
    
    for method in methods:
        print(f"\n{method}:")
        print("-" * 80)
        print(f"{'Threshold':<12} {'Detected':<10} {'Total':<10} {'Sensitivity':<12}")
        print("-" * 80)
        
        for threshold in thresholds:
            stats = calculate_sensitivity_at_threshold(df, method, threshold)
            threshold_results.append({
                'method': method,
                'threshold': threshold,
                'detected': stats['detected'],
                'total': stats['total'],
                'sensitivity': stats['sensitivity']
            })
            print(f"{threshold:>10.1f}  {stats['detected']:>8}  {stats['total']:>8}  {stats['sensitivity']:>10.2%}")
    
    # Save threshold optimization results
    threshold_df = pd.DataFrame(threshold_results)
    threshold_file = output_dir / 'threshold_optimization.tsv'
    threshold_df.to_csv(threshold_file, sep='\t', index=False)
    logger.info(f"Threshold optimization results saved to {threshold_file}")
    
    # Detailed analysis by FF and deletion length for each threshold
    print("\n" + "="*80)
    print("DETAILED ANALYSIS BY FF AND DELETION LENGTH")
    print("="*80)
    
    detailed_results = []
    
    for method in methods:
        for threshold in thresholds:
            grouped = analyze_by_group(df, method, threshold, ['ff', 'deletion_length'])
            if len(grouped) > 0:
                detailed_results.append(grouped)
    
    if detailed_results:
        detailed_df = pd.concat(detailed_results, ignore_index=True)
        detailed_file = output_dir / 'detailed_threshold_analysis.tsv'
        detailed_df.to_csv(detailed_file, sep='\t', index=False)
        logger.info(f"Detailed analysis saved to {detailed_file}")
        
        # Print summary for current threshold and problematic cases
        print(f"\nCurrent threshold ({current_threshold}) - Problematic cases (low sensitivity):")
        current_detailed = detailed_df[
            (detailed_df['threshold'] == current_threshold) & 
            (detailed_df['sensitivity'] < 0.5)
        ].sort_values(['ff', 'deletion_length'])
        
        if len(current_detailed) > 0:
            print(f"\n{'Method':<12} {'FF':<6} {'DelLen':<8} {'Detected':<10} {'Total':<10} {'Sensitivity':<12}")
            print("-" * 80)
            for _, row in current_detailed.iterrows():
                print(f"{row['method']:<12} {row['ff']:>4.0f}%  {row['deletion_length']:>6.0f}Mb  "
                      f"{row['detected_count']:>8.0f}  {row['total_count']:>8.0f}  {row['sensitivity']:>10.2%}")
        else:
            print("  No problematic cases found at current threshold.")
    
    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)
    logger.info(f"Results saved to {output_dir}")


if __name__ == '__main__':
    main()

