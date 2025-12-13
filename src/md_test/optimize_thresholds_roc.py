#!/usr/bin/env python3
"""
Optimize Z-score Thresholds and Minimum Detection Length using ROC Analysis

This script analyzes ROC curves for each method and FF combination to find
optimal z-score threshold and minimum detection length.

Usage:
    python optimize_thresholds_roc.py \
        --input zscore_extraction.tsv \
        --sample-dir /data/md_validation/1p36 \
        --outdir roc_analysis
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from sklearn.metrics import roc_curve, auc

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)


def parse_wc_report_all_regions(report_file: Path) -> list:
    """Parse WC report.txt and extract all detected regions"""
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
                
                fields = line.split()
                if len(fields) < 4:
                    continue
                
                try:
                    zscore = float(fields[0])
                    location = fields[3]
                    if ':' in location and '-' in location:
                        chr_part, coord_part = location.split(':', 1)
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
        logger.debug(f"Failed to parse WC report {report_file}: {e}")
    
    return regions


def parse_wcx_bed_all_regions(bed_file: Path) -> list:
    """Parse WCX aberrations.bed and extract all detected regions"""
    regions = []
    if not bed_file.exists():
        return regions
    
    try:
        with open(bed_file, 'r') as f:
            lines = f.readlines()
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('track'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 5:
                fields = line.split()
                if len(fields) < 5:
                    continue
            
            try:
                chr_name = str(fields[0]).replace('chr', '')
                start = int(fields[1])
                end = int(fields[2])
                zscore = float(fields[4])
                regions.append({
                    'chr': chr_name,
                    'start': start,
                    'end': end,
                    'zscore': zscore
                })
            except (ValueError, IndexError):
                continue
    except Exception as e:
        logger.debug(f"Failed to parse WCX bed {bed_file}: {e}")
    
    return regions


def check_overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
    """Check if two regions overlap"""
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    return overlap_start < overlap_end


def calculate_roc_data(
    df: pd.DataFrame,
    sample_dir: Path,
    method: str,
    output_type: str,
    ff_value: float,
    min_length_candidates: List[int]
) -> Dict:
    """Calculate ROC curve data for a specific method and FF"""
    
    # Filter by FF
    ff_df = df[df['ff'] == ff_value].copy()
    
    if len(ff_df) == 0:
        logger.warning(f"No samples found for FF={ff_value}%")
        return None
    
    method_name = f'{method}_{output_type}'
    zscore_col = method_name + '_zscore'
    
    # Get all z-scores (for ROC x-axis)
    all_zscores = []
    
    results_by_min_length = {}
    
    for min_length in min_length_candidates:
        # Collect data for each sample
        sample_data = []
        
        for _, row in ff_df.iterrows():
            sample_name = row['sample_name']
            sample_dir_path = sample_dir / sample_name
            
            expected_deletion = {
                'chromosome': row['expected_deletion_chr'],
                'start': row['expected_deletion_start'],
                'end': row['expected_deletion_end']
            }
            
            # Get all detected regions
            if method == 'WCX':
                bed_file = sample_dir_path / f"Output_WCX" / output_type / f"{sample_name}.wcx.{output_type}_aberrations.bed"
                detected_regions = parse_wcx_bed_all_regions(bed_file)
            else:  # WC
                report_file = sample_dir_path / f"Output_WC" / output_type / f"{sample_name}.wc.{output_type}.report.txt"
                detected_regions = parse_wc_report_all_regions(report_file)
            
            if not detected_regions:
                continue
            
            expected_chr = str(expected_deletion['chromosome']).replace('chr', '')
            expected_start = expected_deletion['start']
            expected_end = expected_deletion['end']
            
            # Find target z-score (most significant in target region)
            target_zscore = None
            for region in detected_regions:
                if region['chr'] == expected_chr:
                    if check_overlap(expected_start, expected_end, region['start'], region['end']):
                        if target_zscore is None or abs(region['zscore']) > abs(target_zscore):
                            target_zscore = region['zscore']
            
            # Find FP z-scores (most significant in other regions, considering min_length)
            fp_zscore = None
            for region in detected_regions:
                region_length = region['end'] - region['start']
                
                # Skip if too short
                if region_length < min_length:
                    continue
                
                # Check if this is NOT target region
                is_target = False
                if region['chr'] == expected_chr:
                    if check_overlap(expected_start, expected_end, region['start'], region['end']):
                        is_target = True
                
                if not is_target:
                    if fp_zscore is None or abs(region['zscore']) > abs(fp_zscore):
                        fp_zscore = region['zscore']
            
            # Record for this sample
            sample_data.append({
                'sample_name': sample_name,
                'target_zscore': target_zscore,
                'fp_zscore': fp_zscore,
                'has_target': target_zscore is not None,
                'has_fp': fp_zscore is not None
            })
            
            # Collect all z-scores for threshold candidates
            if target_zscore is not None:
                all_zscores.append(target_zscore)
            if fp_zscore is not None:
                all_zscores.append(fp_zscore)
        
        results_by_min_length[min_length] = sample_data
    
    # Generate threshold candidates from collected z-scores
    if all_zscores:
        all_zscores = np.array(all_zscores)
        zcut_candidates = np.percentile(all_zscores, np.arange(0, 100, 1))
        zcut_candidates = np.unique(np.round(zcut_candidates, 1))
        zcut_candidates = np.sort(zcut_candidates)  # From least to most negative
    else:
        zcut_candidates = np.arange(-50, -2.5, 0.5)
    
    # Calculate TPR and FPR for each threshold and min_length
    roc_results = {}
    
    for min_length, sample_data in results_by_min_length.items():
        tpr_list = []
        fpr_list = []
        threshold_list = []
        sensitivity_list = []
        specificity_list = []
        ppv_list = []
        
        for zcut in zcut_candidates:
            tp = 0
            fp = 0
            fn = 0
            tn = 0
            
            for sample in sample_data:
                # True positive: target detected (z-score <= threshold)
                if sample['has_target'] and sample['target_zscore'] <= zcut:
                    tp += 1
                elif sample['has_target']:
                    fn += 1
                
                # False positive: FP detected (z-score <= threshold)
                if sample['has_fp'] and sample['fp_zscore'] <= zcut:
                    fp += 1
                else:
                    tn += 1
            
            # Calculate metrics
            sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            specificity = tn / (tn + fp) if (tn + fp) > 0 else 1.0
            ppv = tp / (tp + fp) if (tp + fp) > 0 else 0.0
            fpr = 1 - specificity
            
            tpr_list.append(sensitivity)
            fpr_list.append(fpr)
            threshold_list.append(zcut)
            sensitivity_list.append(sensitivity)
            specificity_list.append(specificity)
            ppv_list.append(ppv)
        
        # Calculate AUC
        if len(tpr_list) > 1:
            roc_auc = auc(fpr_list, tpr_list)
        else:
            roc_auc = 0.0
        
        roc_results[min_length] = {
            'fpr': fpr_list,
            'tpr': tpr_list,
            'thresholds': threshold_list,
            'sensitivity': sensitivity_list,
            'specificity': specificity_list,
            'ppv': ppv_list,
            'auc': roc_auc,
            'n_samples': len(sample_data)
        }
    
    return {
        'method': method_name,
        'ff': ff_value,
        'roc_by_min_length': roc_results
    }


def create_interactive_roc_plot(
    roc_data: Dict,
    output_file: Path
) -> None:
    """Create interactive ROC plot with threshold selection"""
    
    method = roc_data['method']
    ff = roc_data['ff']
    roc_by_min_length = roc_data['roc_by_min_length']
    
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            'ROC Curve',
            'Sensitivity vs Threshold',
            'PPV vs Threshold',
            'Sensitivity vs PPV'
        ),
        specs=[[{"secondary_y": False}, {"secondary_y": False}],
               [{"secondary_y": False}, {"secondary_y": False}]]
    )
    
    colors = px.colors.qualitative.Plotly
    
    for idx, (min_length, roc_result) in enumerate(roc_by_min_length.items()):
        color = colors[idx % len(colors)]
        min_length_mb = min_length / 1_000_000
        name = f'MinLen={min_length_mb:.1f}Mb (AUC={roc_result["auc"]:.3f})'
        
        # ROC Curve
        fig.add_trace(
            go.Scatter(
                x=roc_result['fpr'],
                y=roc_result['tpr'],
                mode='lines+markers',
                name=name,
                line=dict(color=color, width=2),
                marker=dict(size=4),
                hovertemplate=(
                    f'MinLen={min_length_mb:.1f}Mb<br>' +
                    'FPR=%{x:.3f}<br>' +
                    'TPR=%{y:.3f}<br>' +
                    'Threshold=%{text}<extra></extra>'
                ),
                text=[f'{t:.1f}' for t in roc_result['thresholds']],
                showlegend=True
            ),
            row=1, col=1
        )
        
        # Sensitivity vs Threshold
        fig.add_trace(
            go.Scatter(
                x=roc_result['thresholds'],
                y=roc_result['sensitivity'],
                mode='lines+markers',
                name=name,
                line=dict(color=color, width=2),
                marker=dict(size=4),
                hovertemplate=(
                    'Threshold=%{x:.1f}<br>' +
                    'Sensitivity=%{y:.3f}<extra></extra>'
                ),
                showlegend=False
            ),
            row=1, col=2
        )
        
        # PPV vs Threshold
        fig.add_trace(
            go.Scatter(
                x=roc_result['thresholds'],
                y=roc_result['ppv'],
                mode='lines+markers',
                name=name,
                line=dict(color=color, width=2),
                marker=dict(size=4),
                hovertemplate=(
                    'Threshold=%{x:.1f}<br>' +
                    'PPV=%{y:.3f}<extra></extra>'
                ),
                showlegend=False
            ),
            row=2, col=1
        )
        
        # Sensitivity vs PPV (trade-off)
        fig.add_trace(
            go.Scatter(
                x=roc_result['ppv'],
                y=roc_result['sensitivity'],
                mode='lines+markers',
                name=name,
                line=dict(color=color, width=2),
                marker=dict(size=4),
                hovertemplate=(
                    'PPV=%{x:.3f}<br>' +
                    'Sensitivity=%{y:.3f}<br>' +
                    'Threshold=%{text}<extra></extra>'
                ),
                text=[f'{t:.1f}' for t in roc_result['thresholds']],
                showlegend=False
            ),
            row=2, col=2
        )
    
    # Add diagonal reference line for ROC
    fig.add_trace(
        go.Scatter(
            x=[0, 1],
            y=[0, 1],
            mode='lines',
            name='Random',
            line=dict(color='gray', width=1, dash='dash'),
            showlegend=False
        ),
        row=1, col=1
    )
    
    # Update layout
    fig.update_xaxes(title_text="False Positive Rate", row=1, col=1)
    fig.update_yaxes(title_text="True Positive Rate", row=1, col=1)
    
    fig.update_xaxes(title_text="Z-score Threshold", row=1, col=2)
    fig.update_yaxes(title_text="Sensitivity", row=1, col=2)
    
    fig.update_xaxes(title_text="Z-score Threshold", row=2, col=1)
    fig.update_yaxes(title_text="PPV", row=2, col=1)
    
    fig.update_xaxes(title_text="PPV", row=2, col=2)
    fig.update_yaxes(title_text="Sensitivity", row=2, col=2)
    
    fig.update_layout(
        title_text=f'ROC Analysis: {method}, FF={ff}%',
        height=800,
        width=1200,
        hovermode='closest',
        template='plotly_white'
    )
    
    fig.write_html(output_file)
    logger.info(f"Saved interactive ROC plot: {output_file}")


def find_optimal_points(roc_data: Dict) -> pd.DataFrame:
    """Find optimal operating points for different criteria"""
    
    method = roc_data['method']
    ff = roc_data['ff']
    roc_by_min_length = roc_data['roc_by_min_length']
    
    results = []
    
    for min_length, roc_result in roc_by_min_length.items():
        thresholds = np.array(roc_result['thresholds'])
        sensitivity = np.array(roc_result['sensitivity'])
        specificity = np.array(roc_result['specificity'])
        ppv = np.array(roc_result['ppv'])
        
        # Criterion 1: Maximize Youden's J statistic (Sens + Spec - 1)
        j_stats = sensitivity + specificity - 1
        best_j_idx = np.argmax(j_stats)
        
        # Criterion 2: PPV >= 90%, maximize sensitivity
        ppv_90_mask = ppv >= 0.9
        if ppv_90_mask.any():
            best_ppv90_idx = np.argmax(sensitivity[ppv_90_mask])
            actual_idx_ppv90 = np.where(ppv_90_mask)[0][best_ppv90_idx]
        else:
            actual_idx_ppv90 = None
        
        # Criterion 3: Sensitivity >= 80%, maximize PPV
        sens_80_mask = sensitivity >= 0.8
        if sens_80_mask.any():
            best_sens80_idx = np.argmax(ppv[sens_80_mask])
            actual_idx_sens80 = np.where(sens_80_mask)[0][best_sens80_idx]
        else:
            actual_idx_sens80 = None
        
        # Record results
        results.append({
            'method': method,
            'ff': ff,
            'min_length_bp': min_length,
            'min_length_mb': min_length / 1_000_000,
            'criterion': 'Youden_J',
            'threshold': thresholds[best_j_idx],
            'sensitivity': sensitivity[best_j_idx],
            'specificity': specificity[best_j_idx],
            'ppv': ppv[best_j_idx],
            'j_statistic': j_stats[best_j_idx]
        })
        
        if actual_idx_ppv90 is not None:
            results.append({
                'method': method,
                'ff': ff,
                'min_length_bp': min_length,
                'min_length_mb': min_length / 1_000_000,
                'criterion': 'PPV>=90%',
                'threshold': thresholds[actual_idx_ppv90],
                'sensitivity': sensitivity[actual_idx_ppv90],
                'specificity': specificity[actual_idx_ppv90],
                'ppv': ppv[actual_idx_ppv90],
                'j_statistic': j_stats[actual_idx_ppv90]
            })
        
        if actual_idx_sens80 is not None:
            results.append({
                'method': method,
                'ff': ff,
                'min_length_bp': min_length,
                'min_length_mb': min_length / 1_000_000,
                'criterion': 'Sens>=80%',
                'threshold': thresholds[actual_idx_sens80],
                'sensitivity': sensitivity[actual_idx_sens80],
                'specificity': specificity[actual_idx_sens80],
                'ppv': ppv[actual_idx_sens80],
                'j_statistic': j_stats[actual_idx_sens80]
            })
    
    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(
        description="Optimize Z-score Thresholds using ROC Analysis"
    )
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Input zscore extraction TSV file'
    )
    parser.add_argument(
        '--sample-dir',
        type=str,
        required=True,
        help='Directory containing sample directories'
    )
    parser.add_argument(
        '--outdir',
        type=str,
        required=True,
        help='Output directory for ROC analysis results'
    )
    parser.add_argument(
        '--min-lengths',
        type=str,
        default='0,500000,1000000,2000000',
        help='Comma-separated minimum detection lengths in bp (default: 0,500kb,1Mb,2Mb)'
    )
    
    args = parser.parse_args()
    
    # Parse min lengths
    min_length_candidates = [int(x.strip()) for x in args.min_lengths.split(',')]
    logger.info(f"Min length candidates: {[x/1e6 for x in min_length_candidates]} Mb")
    
    # Create output directory
    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = output_dir / 'plots'
    plots_dir.mkdir(exist_ok=True)
    
    # Load data
    input_file = Path(args.input)
    logger.info(f"Loading data from {input_file}")
    df = pd.read_csv(input_file, sep='\t')
    
    # Convert columns
    df['ff'] = pd.to_numeric(df['ff'], errors='coerce')
    for col in ['WC_orig_zscore', 'WC_fetus_zscore', 'WCX_orig_zscore', 'WCX_fetus_zscore']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    sample_dir = Path(args.sample_dir)
    
    # Methods and FF values to analyze
    methods = [
        ('WC', 'orig'),
        ('WC', 'fetus'),
        ('WCX', 'orig'),
        ('WCX', 'fetus')
    ]
    ff_values = sorted(df['ff'].dropna().unique())
    
    logger.info(f"Analyzing {len(methods)} methods × {len(ff_values)} FF values")
    
    all_optimal_points = []
    
    # Analyze each method and FF combination
    for method, output_type in methods:
        for ff in ff_values:
            logger.info(f"Analyzing {method}_{output_type}, FF={ff}%...")
            
            roc_data = calculate_roc_data(
                df, sample_dir, method, output_type, ff, min_length_candidates
            )
            
            if roc_data is None:
                continue
            
            # Create interactive plot
            plot_file = plots_dir / f'roc_{method}_{output_type}_FF{int(ff)}.html'
            create_interactive_roc_plot(roc_data, plot_file)
            
            # Find optimal points
            optimal_points = find_optimal_points(roc_data)
            all_optimal_points.append(optimal_points)
    
    # Save optimal points
    if all_optimal_points:
        optimal_df = pd.concat(all_optimal_points, ignore_index=True)
        optimal_file = output_dir / 'optimal_thresholds.csv'
        optimal_df.to_csv(optimal_file, index=False, float_format='%.4f')
        logger.info(f"Saved optimal thresholds to {optimal_file}")
        
        # Create summary
        summary_file = output_dir / 'summary.txt'
        with open(summary_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("OPTIMAL THRESHOLD SUMMARY\n")
            f.write("="*80 + "\n\n")
            
            for criterion in optimal_df['criterion'].unique():
                f.write(f"\n{criterion}:\n")
                f.write("-"*80 + "\n")
                criterion_df = optimal_df[optimal_df['criterion'] == criterion]
                for _, row in criterion_df.iterrows():
                    f.write(
                        f"  {row['method']:<15} FF={row['ff']:>5.0f}% "
                        f"MinLen={row['min_length_mb']:>4.1f}Mb: "
                        f"z={row['threshold']:>6.1f} "
                        f"(Sens={row['sensitivity']:>6.2%}, PPV={row['ppv']:>6.2%})\n"
                    )
        
        logger.info(f"Saved summary to {summary_file}")
    
    logger.info("="*80)
    logger.info("ROC Analysis complete!")
    logger.info(f"Interactive plots saved to: {plots_dir}")
    logger.info(f"Results saved to: {output_dir}")
    logger.info("="*80)
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())




