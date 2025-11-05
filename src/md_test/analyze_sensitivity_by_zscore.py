#!/usr/bin/env python3
"""
Analyze sensitivity by z-score threshold with adaptive binning

Reads zscore_extraction.tsv and:
1. Analyzes z-score distribution
2. Divides z-score range into bins
3. Calculates sensitivity for each bin threshold
4. Identifies regions with large sensitivity changes
5. Refines analysis in those regions
"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def calculate_sensitivity_at_threshold(df: pd.DataFrame, method: str, threshold: float) -> float:
    """Calculate sensitivity at a given z-score threshold"""
    zscore_col = f'{method}_zscore'
    
    # Filter rows with valid z-score
    valid_df = df[df[zscore_col].notna()].copy()
    
    if len(valid_df) == 0:
        return 0.0
    
    # For deletions, z-score is negative. Detection means z-score <= threshold (more negative)
    detected = (valid_df[zscore_col] <= threshold).sum()
    total = len(valid_df)
    
    sensitivity = detected / total if total > 0 else 0.0
    return sensitivity


def analyze_zscore_distribution(df: pd.DataFrame, method: str) -> dict:
    """Analyze z-score distribution for a method"""
    zscore_col = f'{method}_zscore'
    
    valid_scores = df[df[zscore_col].notna()][zscore_col].values
    
    if len(valid_scores) == 0:
        return None
    
    return {
        'min': float(np.min(valid_scores)),
        'max': float(np.max(valid_scores)),
        'mean': float(np.mean(valid_scores)),
        'median': float(np.median(valid_scores)),
        'std': float(np.std(valid_scores)),
        'q25': float(np.percentile(valid_scores, 25)),
        'q75': float(np.percentile(valid_scores, 75)),
        'count': len(valid_scores)
    }


def find_sensitivity_change_points(thresholds: list, sensitivities: list, min_change: float = 0.05) -> list:
    """Find threshold points where sensitivity changes significantly
    
    Args:
        thresholds: List of threshold values
        sensitivities: List of sensitivity values corresponding to thresholds
        min_change: Minimum change in sensitivity to consider significant
    
    Returns:
        List of threshold indices where significant changes occur
    """
    change_points = []
    
    for i in range(1, len(sensitivities)):
        change = abs(sensitivities[i] - sensitivities[i-1])
        if change >= min_change:
            change_points.append(i)
    
    return change_points


def analyze_with_adaptive_binning(df: pd.DataFrame, method: str, 
                                   initial_bins: int = 20, 
                                   refinement_factor: int = 5,
                                   min_change: float = 0.05) -> pd.DataFrame:
    """Analyze sensitivity with adaptive binning
    
    Args:
        df: DataFrame with z-score data
        method: Method name (e.g., 'WC_orig')
        initial_bins: Number of initial bins
        refinement_factor: Factor to refine bins in regions with large changes
        min_change: Minimum sensitivity change to trigger refinement
    
    Returns:
        DataFrame with threshold and sensitivity values
    """
    zscore_col = f'{method}_zscore'
    
    # Get valid z-scores
    valid_scores = df[df[zscore_col].notna()][zscore_col].values
    
    if len(valid_scores) == 0:
        return pd.DataFrame()
    
    # For deletions, z-scores are negative. We want to test thresholds from min to max
    min_zscore = float(np.min(valid_scores))
    max_zscore = float(np.max(valid_scores))
    
    # Create initial bins (from most negative to least negative)
    initial_thresholds = np.linspace(min_zscore, max_zscore, initial_bins)
    
    # Calculate sensitivity for each initial threshold
    results = []
    for threshold in initial_thresholds:
        sensitivity = calculate_sensitivity_at_threshold(df, method, threshold)
        results.append({
            'threshold': threshold,
            'sensitivity': sensitivity,
            'bin_level': 1
        })
    
    # Find regions with large sensitivity changes
    initial_df = pd.DataFrame(results)
    initial_df = initial_df.sort_values('threshold')
    
    sensitivities = initial_df['sensitivity'].values
    thresholds = initial_df['threshold'].values
    
    change_points = find_sensitivity_change_points(thresholds, sensitivities, min_change)
    
    # Refine bins in regions with large changes
    refined_results = initial_df.to_dict('records')
    
    if change_points:
        logger.info(f"Found {len(change_points)} regions with significant sensitivity changes")
        
        # Refine each region
        for idx in change_points:
            if idx == 0:
                continue
            
            # Get the region boundaries
            start_threshold = thresholds[idx - 1]
            end_threshold = thresholds[idx]
            
            # Create refined bins in this region
            refined_thresholds = np.linspace(start_threshold, end_threshold, refinement_factor + 2)
            # Remove endpoints (already in initial bins)
            refined_thresholds = refined_thresholds[1:-1]
            
            for threshold in refined_thresholds:
                sensitivity = calculate_sensitivity_at_threshold(df, method, threshold)
                refined_results.append({
                    'threshold': threshold,
                    'sensitivity': sensitivity,
                    'bin_level': 2
                })
    
    # Sort by threshold and remove duplicates
    result_df = pd.DataFrame(refined_results)
    result_df = result_df.sort_values('threshold')
    result_df = result_df.drop_duplicates(subset=['threshold'], keep='first')
    
    result_df['method'] = method
    result_df['detected_count'] = result_df['threshold'].apply(
        lambda t: (df[df[zscore_col].notna()][zscore_col] <= t).sum()
    )
    result_df['total_count'] = df[df[zscore_col].notna()].shape[0]
    
    return result_df


def main():
    parser = argparse.ArgumentParser(
        description="Analyze sensitivity by z-score threshold with adaptive binning"
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
        '--initial-bins',
        type=int,
        default=20,
        help='Number of initial bins (default: 20)'
    )
    parser.add_argument(
        '--refinement-factor',
        type=int,
        default=5,
        help='Number of refined bins in regions with large changes (default: 5)'
    )
    parser.add_argument(
        '--min-change',
        type=float,
        default=0.05,
        help='Minimum sensitivity change to trigger refinement (default: 0.05)'
    )
    
    args = parser.parse_args()
    
    input_file = Path(args.input)
    output_dir = Path(args.output_dir)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        return
    
    # Read z-score extraction results
    logger.info(f"Reading z-score data from {input_file}")
    df = pd.read_csv(input_file, sep='\t')
    
    logger.info(f"Found {len(df)} samples")
    
    methods = ['WC_orig', 'WC_fetus', 'WCX_orig', 'WCX_fetus']
    
    # Analyze z-score distribution for each method
    print("\n" + "="*80)
    print("Z-SCORE DISTRIBUTION")
    print("="*80)
    
    distributions = {}
    for method in methods:
        dist = analyze_zscore_distribution(df, method)
        if dist:
            distributions[method] = dist
            print(f"\n{method}:")
            print(f"  Count: {dist['count']}")
            print(f"  Min: {dist['min']:.2f}")
            print(f"  Max: {dist['max']:.2f}")
            print(f"  Mean: {dist['mean']:.2f}")
            print(f"  Median: {dist['median']:.2f}")
            print(f"  Q25: {dist['q25']:.2f}")
            print(f"  Q75: {dist['q75']:.2f}")
            print(f"  Std: {dist['std']:.2f}")
    
    # Analyze sensitivity by threshold with adaptive binning
    print("\n" + "="*80)
    print("SENSITIVITY ANALYSIS BY Z-SCORE THRESHOLD")
    print("="*80)
    
    all_results = []
    
    for method in methods:
        print(f"\n{method}:")
        print("-" * 80)
        
        if method not in distributions:
            print("  No data available")
            continue
        
        # Perform adaptive binning analysis
        result_df = analyze_with_adaptive_binning(
            df, method, 
            initial_bins=args.initial_bins,
            refinement_factor=args.refinement_factor,
            min_change=args.min_change
        )
        
        if len(result_df) == 0:
            print("  No results")
            continue
        
        # Print results (show key thresholds)
        print(f"{'Threshold':<12} {'Sensitivity':<12} {'Detected':<10} {'Total':<10}")
        print("-" * 80)
        
        # Show every 5th result for readability, or all if less than 20
        step = max(1, len(result_df) // 20) if len(result_df) > 20 else 1
        for idx in range(0, len(result_df), step):
            row = result_df.iloc[idx]
            print(f"{row['threshold']:>10.2f}  {row['sensitivity']:>10.2%}  "
                  f"{row['detected_count']:>8.0f}  {row['total_count']:>8.0f}")
        
        # Show last result
        if len(result_df) > 1 and (len(result_df) - 1) % step != 0:
            row = result_df.iloc[-1]
            print(f"{row['threshold']:>10.2f}  {row['sensitivity']:>10.2%}  "
                  f"{row['detected_count']:>8.0f}  {row['total_count']:>8.0f}")
        
        # Find thresholds with specific sensitivity levels
        print(f"\n  Key Thresholds:")
        for target_sens in [0.5, 0.7, 0.8, 0.9, 0.95]:
            # Find threshold closest to target sensitivity
            result_df['sens_diff'] = abs(result_df['sensitivity'] - target_sens)
            closest = result_df.loc[result_df['sens_diff'].idxmin()]
            print(f"    Sensitivity {target_sens:.0%}: threshold = {closest['threshold']:.2f}")
        
        all_results.append(result_df)
    
    # Save results
    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        output_file = output_dir / 'sensitivity_by_zscore.tsv'
        combined_df = combined_df.drop('sens_diff', axis=1, errors='ignore')
        combined_df.to_csv(output_file, sep='\t', index=False, float_format='%.2f')
        logger.info(f"Results saved to {output_file}")
        
        # Save distribution summary
        dist_df = pd.DataFrame(distributions).T
        dist_file = output_dir / 'zscore_distribution.tsv'
        dist_df.to_csv(dist_file, sep='\t', float_format='%.2f')
        logger.info(f"Distribution summary saved to {dist_file}")
    
    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)
    logger.info(f"Results saved to {output_dir}")


if __name__ == '__main__':
    main()

