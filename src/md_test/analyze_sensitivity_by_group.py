#!/usr/bin/env python3
"""
Analyze sensitivity and PPV by FF and deletion length

Reads zscore_extraction.tsv and calculates performance metrics
grouped by FF and deletion length.
"""

import argparse
from pathlib import Path
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def calculate_metrics(df: pd.DataFrame, method: str, threshold: float = None) -> pd.DataFrame:
    """Calculate sensitivity and PPV metrics
    
    Args:
        df: DataFrame with z-score data
        method: One of 'WC_orig', 'WC_fetus', 'WCX_orig', 'WCX_fetus'
        threshold: Z-score threshold (optional, if None, considers any z-score as detected)
    
    Returns:
        DataFrame with metrics grouped by FF and deletion_length_mb
    """
    zscore_col = f'{method}_zscore'
    
    # Create detected column
    if threshold is not None:
        # Detection: z-score <= threshold (more negative)
        df['detected'] = df[zscore_col].notna() & (df[zscore_col] <= threshold)
    else:
        # Detection: any z-score present
        df['detected'] = df[zscore_col].notna()
    
    # Group by FF and deletion_length_mb
    group_cols = ['ff', 'deletion_length_mb']
    
    # Filter rows with valid grouping info
    valid_df = df[df[group_cols].notna().all(axis=1)].copy()
    
    if len(valid_df) == 0:
        return pd.DataFrame()
    
    # Calculate metrics for each group
    grouped = valid_df.groupby(group_cols).agg({
        'detected': ['sum', 'count']
    }).reset_index()
    
    grouped.columns = group_cols + ['detected_count', 'total_count']
    grouped['not_detected_count'] = grouped['total_count'] - grouped['detected_count']
    
    # Sensitivity = TP / (TP + FN)
    # In this case: detected / total (since all samples have deletions)
    grouped['sensitivity'] = grouped['detected_count'] / grouped['total_count']
    
    # PPV = TP / (TP + FP)
    # Since all samples have deletions, FP = 0, so PPV = 1.0 for detected samples
    # But we can't calculate PPV without true negatives, so we'll skip it
    # grouped['ppv'] = 1.0  # Always 1.0 since all samples have deletions
    
    grouped['method'] = method
    if threshold is not None:
        grouped['threshold'] = threshold
    
    return grouped


def main():
    parser = argparse.ArgumentParser(
        description="Analyze sensitivity and PPV by FF and deletion length"
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
        '-t', '--threshold',
        type=float,
        default=None,
        help='Z-score threshold (if not provided, uses any z-score as detected)'
    )
    
    args = parser.parse_args()
    
    input_file = Path(args.input)
    output_dir = Path(args.output_dir)
    threshold = args.threshold
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        return
    
    # Read z-score extraction results
    logger.info(f"Reading z-score data from {input_file}")
    df = pd.read_csv(input_file, sep='\t')
    
    logger.info(f"Found {len(df)} samples")
    
    # Check data availability
    methods = ['WC_orig', 'WC_fetus', 'WCX_orig', 'WCX_fetus']
    
    print("\n" + "="*80)
    print("DATA AVAILABILITY")
    print("="*80)
    for method in methods:
        zscore_col = f'{method}_zscore'
        valid_count = df[df[zscore_col].notna()].shape[0]
        print(f"{method}: {valid_count}/{len(df)} samples with z-score data")
    
    # Calculate metrics for each method
    all_results = []
    
    print("\n" + "="*80)
    if threshold is not None:
        print(f"ANALYSIS AT THRESHOLD: {threshold}")
    else:
        print("ANALYSIS: Any z-score present (no threshold)")
    print("="*80)
    
    for method in methods:
        print(f"\n{method}:")
        print("-" * 80)
        
        metrics_df = calculate_metrics(df, method, threshold)
        
        if len(metrics_df) == 0:
            print("  No data available")
            continue
        
        # Sort by FF and deletion_length_mb
        metrics_df = metrics_df.sort_values(['ff', 'deletion_length_mb'])
        
        # Print results
        print(f"{'FF':<8} {'DelLen(Mb)':<12} {'Detected':<10} {'Total':<10} {'Sensitivity':<12}")
        print("-" * 80)
        
        for _, row in metrics_df.iterrows():
            print(f"{row['ff']:>6.0f}  {row['deletion_length_mb']:>10.2f}  "
                  f"{row['detected_count']:>8.0f}  {row['total_count']:>8.0f}  "
                  f"{row['sensitivity']:>10.2%}")
        
        # Summary by FF only
        print(f"\n  Summary by FF:")
        ff_summary = metrics_df.groupby('ff').agg({
            'detected_count': 'sum',
            'total_count': 'sum'
        }).reset_index()
        ff_summary['sensitivity'] = ff_summary['detected_count'] / ff_summary['total_count']
        
        print(f"{'FF':<8} {'Detected':<10} {'Total':<10} {'Sensitivity':<12}")
        print("-" * 80)
        for _, row in ff_summary.iterrows():
            print(f"{row['ff']:>6.0f}  {row['detected_count']:>8.0f}  "
                  f"{row['total_count']:>8.0f}  {row['sensitivity']:>10.2%}")
        
        # Summary by deletion length only
        print(f"\n  Summary by Deletion Length:")
        del_summary = metrics_df.groupby('deletion_length_mb').agg({
            'detected_count': 'sum',
            'total_count': 'sum'
        }).reset_index()
        del_summary['sensitivity'] = del_summary['detected_count'] / del_summary['total_count']
        
        print(f"{'DelLen(Mb)':<12} {'Detected':<10} {'Total':<10} {'Sensitivity':<12}")
        print("-" * 80)
        for _, row in del_summary.iterrows():
            print(f"{row['deletion_length_mb']:>10.2f}  {row['detected_count']:>8.0f}  "
                  f"{row['total_count']:>8.0f}  {row['sensitivity']:>10.2%}")
        
        all_results.append(metrics_df)
    
    # Save detailed results
    if all_results:
        detailed_df = pd.concat(all_results, ignore_index=True)
        output_file = output_dir / 'sensitivity_by_group.tsv'
        detailed_df.to_csv(output_file, sep='\t', index=False, float_format='%.2f')
        logger.info(f"Detailed results saved to {output_file}")
        
        # Save summary by FF
        ff_summary_all = []
        for method in methods:
            metrics_df = calculate_metrics(df, method, threshold)
            if len(metrics_df) > 0:
                ff_summary = metrics_df.groupby('ff').agg({
                    'detected_count': 'sum',
                    'total_count': 'sum'
                }).reset_index()
                ff_summary['sensitivity'] = ff_summary['detected_count'] / ff_summary['total_count']
                ff_summary['method'] = method
                if threshold is not None:
                    ff_summary['threshold'] = threshold
                ff_summary_all.append(ff_summary)
        
        if ff_summary_all:
            ff_summary_df = pd.concat(ff_summary_all, ignore_index=True)
            ff_output_file = output_dir / 'sensitivity_by_ff.tsv'
            ff_summary_df.to_csv(ff_output_file, sep='\t', index=False, float_format='%.2f')
            logger.info(f"FF summary saved to {ff_output_file}")
        
        # Save summary by deletion length
        del_summary_all = []
        for method in methods:
            metrics_df = calculate_metrics(df, method, threshold)
            if len(metrics_df) > 0:
                del_summary = metrics_df.groupby('deletion_length_mb').agg({
                    'detected_count': 'sum',
                    'total_count': 'sum'
                }).reset_index()
                del_summary['sensitivity'] = del_summary['detected_count'] / del_summary['total_count']
                del_summary['method'] = method
                if threshold is not None:
                    del_summary['threshold'] = threshold
                del_summary_all.append(del_summary)
        
        if del_summary_all:
            del_summary_df = pd.concat(del_summary_all, ignore_index=True)
            del_output_file = output_dir / 'sensitivity_by_deletion_length.tsv'
            del_summary_df.to_csv(del_output_file, sep='\t', index=False, float_format='%.2f')
            logger.info(f"Deletion length summary saved to {del_output_file}")
    
    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)
    logger.info(f"Results saved to {output_dir}")


if __name__ == '__main__':
    main()

