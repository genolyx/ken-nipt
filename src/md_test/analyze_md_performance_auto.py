#!/usr/bin/env python3
"""
Automatic Microdeletion Detection Performance Analysis

This script automatically finds optimal z-score thresholds based on PPV requirements
and calculates sensitivity metrics for various detection modes.

Usage:
    python analyze_md_performance_auto.py \
        --input zscore_extraction.tsv \
        --sample-dir /data/md_validation/1p36 \
        --ppv 90,80,70 \
        --min-detect-length 1000000 \
        --outdir results
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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


def calculate_metrics_for_sample(
    sample_dir: Path,
    sample_name: str,
    expected_deletion: dict,
    method: str,
    output_type: str,
    zcut: float,
    min_detect_length: int
) -> dict:
    """Calculate TP, FP, FN, TN for a single sample"""
    result = {'tp': 0, 'fp': 0, 'fn': 0, 'tn': 0}
    
    # Get all detected regions
    if method == 'WCX':
        bed_file = sample_dir / f"Output_WCX" / output_type / f"{sample_name}.wcx.{output_type}_aberrations.bed"
        detected_regions = parse_wcx_bed_all_regions(bed_file)
    else:  # WC
        report_file = sample_dir / f"Output_WC" / output_type / f"{sample_name}.wc.{output_type}.report.txt"
        detected_regions = parse_wc_report_all_regions(report_file)
    
    # Filter by z-score threshold (for deletions, zscore <= threshold, more negative = detected)
    filtered_regions = [r for r in detected_regions if r['zscore'] <= zcut]
    
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
        
        # This is a different region (FP candidate)
        # Only count as FP if length >= min_detect_length
        if region_length >= min_detect_length:
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


def find_optimal_zcut_for_group(
    group_df: pd.DataFrame,
    sample_dir_root: Path,
    methods: List[Tuple[str, str]],
    target_ppv: float,
    min_detect_length: int,
    zcut_candidates: np.ndarray
) -> Optional[Dict]:
    """Find optimal z-score threshold for a specific (FF, length) group"""
    
    best_result = None
    best_sensitivity = -1
    
    for idx, zcut in enumerate(zcut_candidates):
        tp_total = 0
        fp_total = 0
        fn_total = 0
        tn_total = 0
        
        # Log progress for first and key thresholds
        if idx == 0 or idx == len(zcut_candidates) - 1 or idx % 10 == 0:
            logger.debug(f"    Testing zcut={zcut:.1f}...")
        
        # Process each sample in this group
        for _, row in group_df.iterrows():
            sample_name = row['sample_name']
            sample_dir = sample_dir_root / sample_name
            
            expected_deletion = {
                'chromosome': row['expected_deletion_chr'],
                'start': row['expected_deletion_start'],
                'end': row['expected_deletion_end']
            }
            
            # Check all methods with OR logic
            target_detected_any = False
            other_regions_detected_any = False
            
            for method, output_type in methods:
                metrics = calculate_metrics_for_sample(
                    sample_dir, sample_name, expected_deletion,
                    method, output_type, zcut, min_detect_length
                )
                
                if metrics['tp'] == 1:
                    target_detected_any = True
                if metrics['fp'] == 1:
                    other_regions_detected_any = True
            
            # OR logic results
            if target_detected_any:
                tp_total += 1
            else:
                fn_total += 1
            
            if other_regions_detected_any:
                fp_total += 1
            else:
                tn_total += 1
        
        # Calculate metrics
        sensitivity = tp_total / (tp_total + fn_total) if (tp_total + fn_total) > 0 else 0.0
        ppv = tp_total / (tp_total + fp_total) if (tp_total + fp_total) > 0 else 0.0
        
        # Log key metrics for debugging
        if idx == 0 or (tp_total > 0 and ppv >= target_ppv):
            logger.debug(
                f"      zcut={zcut:.1f}: TP={tp_total}, FP={fp_total}, "
                f"Sens={sensitivity:.2%}, PPV={ppv:.2%}"
            )
        
        # Check if PPV meets requirement
        if ppv >= target_ppv:
            if sensitivity > best_sensitivity:
                best_sensitivity = sensitivity
                best_result = {
                    'zcut': zcut,
                    'sensitivity': sensitivity,
                    'ppv': ppv,
                    'tp': tp_total,
                    'fp': fp_total,
                    'fn': fn_total,
                    'tn': tn_total,
                    'total': len(group_df)
                }
    
    return best_result


def calculate_detection_modes(
    group_df: pd.DataFrame,
    sample_dir_root: Path,
    zcut: float,
    min_detect_length: int
) -> Dict[str, Dict]:
    """Calculate metrics for all 7 detection modes"""
    
    modes = {
        'wc_orig': [('WC', 'orig')],
        'wc_fetus': [('WC', 'fetus')],
        'wcx_orig': [('WCX', 'orig')],
        'wcx_fetus': [('WCX', 'fetus')],
        'orig': [('WC', 'orig'), ('WCX', 'orig')],
        'fetus': [('WC', 'fetus'), ('WCX', 'fetus')],
        'any': [('WC', 'orig'), ('WC', 'fetus'), ('WCX', 'orig'), ('WCX', 'fetus')]
    }
    
    results = {}
    
    for mode_name, methods in modes.items():
        tp_total = 0
        fp_total = 0
        fn_total = 0
        tn_total = 0
        
        # Process each sample in this group
        for _, row in group_df.iterrows():
            sample_name = row['sample_name']
            sample_dir = sample_dir_root / sample_name
            
            expected_deletion = {
                'chromosome': row['expected_deletion_chr'],
                'start': row['expected_deletion_start'],
                'end': row['expected_deletion_end']
            }
            
            # Check all methods for this mode with OR logic
            target_detected_any = False
            other_regions_detected_any = False
            
            for method, output_type in methods:
                metrics = calculate_metrics_for_sample(
                    sample_dir, sample_name, expected_deletion,
                    method, output_type, zcut, min_detect_length
                )
                
                if metrics['tp'] == 1:
                    target_detected_any = True
                if metrics['fp'] == 1:
                    other_regions_detected_any = True
            
            # OR logic results
            if target_detected_any:
                tp_total += 1
            else:
                fn_total += 1
            
            if other_regions_detected_any:
                fp_total += 1
            else:
                tn_total += 1
        
        # Calculate metrics
        sensitivity = tp_total / (tp_total + fn_total) if (tp_total + fn_total) > 0 else 0.0
        ppv = tp_total / (tp_total + fp_total) if (tp_total + fp_total) > 0 else 0.0
        specificity = tn_total / (tn_total + fp_total) if (tn_total + fp_total) > 0 else 0.0
        
        results[mode_name] = {
            'sensitivity': sensitivity,
            'ppv': ppv,
            'specificity': specificity,
            'tp': tp_total,
            'fp': fp_total,
            'fn': fn_total,
            'tn': tn_total,
            'total': len(group_df)
        }
    
    return results


def find_optimal_thresholds(
    df: pd.DataFrame,
    sample_dir: Path,
    ppv_targets: List[float],
    min_detect_length: int
) -> pd.DataFrame:
    """Find optimal z-score thresholds for each (FF, length, PPV) combination"""
    
    # Convert to numeric
    df['ff'] = pd.to_numeric(df['ff'], errors='coerce')
    df['deletion_length_mb'] = pd.to_numeric(df['deletion_length_mb'], errors='coerce')
    
    # Remove rows with missing FF or length
    valid_df = df[df[['ff', 'deletion_length_mb']].notna().all(axis=1)]
    
    # Z-score candidates (from stringent to lenient, negative for deletions)
    # More negative = more stringent (fewer detections, higher PPV)
    # Less negative = more lenient (more detections, lower PPV)
    # Start from very stringent (-50) to lenient (-2.5)
    zcut_candidates = np.arange(-50.0, -2.5, 0.5)
    
    # Methods for 'any' mode (used for finding optimal threshold)
    methods_any = [('WC', 'orig'), ('WC', 'fetus'), ('WCX', 'orig'), ('WCX', 'fetus')]
    
    results = []
    
    # Group by FF and length
    grouped = valid_df.groupby(['ff', 'deletion_length_mb'])
    total_groups = len(grouped)
    
    logger.info(f"Processing {total_groups} (FF, length) combinations...")
    
    for idx, ((ff, length), group_df) in enumerate(grouped, 1):
        logger.info(f"[{idx}/{total_groups}] Processing FF={ff}%, Length={length}Mb ({len(group_df)} samples)...")
        
        for target_ppv in ppv_targets:
            # Find optimal z-score for this group and PPV target
            optimal = find_optimal_zcut_for_group(
                group_df,
                sample_dir,
                methods_any,
                target_ppv / 100.0,
                min_detect_length,
                zcut_candidates
            )
            
            if optimal is None:
                logger.warning(
                    f"  No z-score found for FF={ff}%, Length={length}Mb, "
                    f"PPV>={target_ppv}%"
                )
                results.append({
                    'ff': ff,
                    'length': int(length),
                    'ppv_target': target_ppv,
                    'optimal_zcut': None,
                    'sensitivity': 0.0,
                    'ppv': 0.0,
                    'tp': 0,
                    'fp': 0,
                    'fn': group_df.shape[0],
                    'tn': 0,
                    'total': group_df.shape[0],
                    'status': 'NO_VALID_ZCUT'
                })
            else:
                logger.info(
                    f"  PPV>={target_ppv}%: Optimal z={optimal['zcut']:.1f}, "
                    f"Sens={optimal['sensitivity']:.2%}, PPV={optimal['ppv']:.2%}"
                )
                results.append({
                    'ff': ff,
                    'length': int(length),
                    'ppv_target': target_ppv,
                    'optimal_zcut': optimal['zcut'],
                    'sensitivity': optimal['sensitivity'],
                    'ppv': optimal['ppv'],
                    'tp': optimal['tp'],
                    'fp': optimal['fp'],
                    'fn': optimal['fn'],
                    'tn': optimal['tn'],
                    'total': optimal['total'],
                    'status': 'OK'
                })
    
    return pd.DataFrame(results)


def calculate_sensitivity_by_mode(
    df: pd.DataFrame,
    sample_dir: Path,
    optimal_thresholds_df: pd.DataFrame,
    min_detect_length: int
) -> pd.DataFrame:
    """Calculate sensitivity for all 7 detection modes using optimal thresholds"""
    
    # Convert to numeric
    df['ff'] = pd.to_numeric(df['ff'], errors='coerce')
    df['deletion_length_mb'] = pd.to_numeric(df['deletion_length_mb'], errors='coerce')
    
    # Remove rows with missing FF or length
    valid_df = df[df[['ff', 'deletion_length_mb']].notna().all(axis=1)]
    
    results = []
    
    # Group by FF and length
    grouped = valid_df.groupby(['ff', 'deletion_length_mb'])
    
    for (ff, length), group_df in grouped:
        # Get optimal thresholds for this group
        group_thresholds = optimal_thresholds_df[
            (optimal_thresholds_df['ff'] == ff) &
            (optimal_thresholds_df['length'] == int(length))
        ]
        
        for _, threshold_row in group_thresholds.iterrows():
            ppv_target = threshold_row['ppv_target']
            zcut = threshold_row['optimal_zcut']
            
            if pd.isna(zcut):
                # No valid z-score found
                for mode in ['wc_orig', 'wc_fetus', 'wcx_orig', 'wcx_fetus', 'orig', 'fetus', 'any']:
                    results.append({
                        'ff': ff,
                        'length': int(length),
                        'ppv_target': ppv_target,
                        'zcut': None,
                        'mode': mode,
                        'sensitivity': 0.0,
                        'ppv': 0.0,
                        'specificity': 0.0,
                        'tp': 0,
                        'fp': 0,
                        'fn': len(group_df),
                        'tn': 0,
                        'total': len(group_df)
                    })
                continue
            
            # Calculate metrics for all modes
            mode_results = calculate_detection_modes(
                group_df, sample_dir, zcut, min_detect_length
            )
            
            for mode, metrics in mode_results.items():
                results.append({
                    'ff': ff,
                    'length': int(length),
                    'ppv_target': ppv_target,
                    'zcut': zcut,
                    'mode': mode,
                    'sensitivity': metrics['sensitivity'],
                    'ppv': metrics['ppv'],
                    'specificity': metrics['specificity'],
                    'tp': metrics['tp'],
                    'fp': metrics['fp'],
                    'fn': metrics['fn'],
                    'tn': metrics['tn'],
                    'total': metrics['total']
                })
    
    return pd.DataFrame(results)


def save_results(
    optimal_thresholds_df: pd.DataFrame,
    sensitivity_df: pd.DataFrame,
    output_dir: Path
) -> None:
    """Save analysis results"""
    
    # Save optimal thresholds
    optimal_file = output_dir / 'optimal_thresholds.csv'
    optimal_thresholds_df.to_csv(optimal_file, index=False, float_format='%.4f')
    logger.info(f"Saved optimal thresholds to {optimal_file}")
    
    # Save sensitivity matrix
    sensitivity_file = output_dir / 'sensitivity_matrix.csv'
    sensitivity_df.to_csv(sensitivity_file, index=False, float_format='%.4f')
    logger.info(f"Saved sensitivity matrix to {sensitivity_file}")
    
    # Save PPV-specific results
    ppv_dir = output_dir / 'by_ppv'
    ppv_dir.mkdir(parents=True, exist_ok=True)
    
    for ppv_target in sensitivity_df['ppv_target'].unique():
        ppv_data = sensitivity_df[sensitivity_df['ppv_target'] == ppv_target]
        ppv_file = ppv_dir / f'ppv{int(ppv_target)}_results.csv'
        ppv_data.to_csv(ppv_file, index=False, float_format='%.4f')
        logger.info(f"Saved PPV={ppv_target}% results to {ppv_file}")


def create_pivot_tables(
    sensitivity_df: pd.DataFrame,
    output_dir: Path
) -> None:
    """Create pivot tables for each mode, PPV, and zcut"""
    
    tables_dir = output_dir / 'tables'
    tables_dir.mkdir(parents=True, exist_ok=True)
    
    for ppv_target in sorted(sensitivity_df['ppv_target'].unique()):
        for mode in sorted(sensitivity_df['mode'].unique()):
            filtered = sensitivity_df[
                (sensitivity_df['ppv_target'] == ppv_target) &
                (sensitivity_df['mode'] == mode)
            ]
            
            if len(filtered) == 0:
                continue
            
            # Create pivot table
            pivot = filtered.pivot_table(
                index='ff',
                columns='length',
                values='sensitivity',
                aggfunc='mean'
            )
            
            # Rename columns
            pivot.columns = [f'{int(col)}Mb' for col in pivot.columns]
            
            # Save
            filename = f'table_ppv{int(ppv_target)}_mode-{mode}.csv'
            output_file = tables_dir / filename
            pivot.to_csv(output_file, float_format='%.4f')
            
            logger.info(f"Saved pivot table: {filename}")


def create_heatmaps(
    sensitivity_df: pd.DataFrame,
    output_dir: Path
) -> None:
    """Create heatmaps for each mode and PPV target"""
    
    heatmaps_dir = output_dir / 'heatmaps'
    heatmaps_dir.mkdir(parents=True, exist_ok=True)
    
    for ppv_target in sorted(sensitivity_df['ppv_target'].unique()):
        for mode in sorted(sensitivity_df['mode'].unique()):
            filtered = sensitivity_df[
                (sensitivity_df['ppv_target'] == ppv_target) &
                (sensitivity_df['mode'] == mode)
            ]
            
            if len(filtered) == 0:
                continue
            
            # Create pivot for heatmap
            pivot = filtered.pivot_table(
                index='length',
                columns='ff',
                values='sensitivity',
                aggfunc='mean'
            )
            
            # Create heatmap
            fig, ax = plt.subplots(figsize=(10, 8))
            
            sns.heatmap(
                pivot,
                annot=True,
                fmt='.2%',
                cmap='RdYlGn',
                vmin=0,
                vmax=1,
                cbar_kws={'label': 'Sensitivity'},
                ax=ax
            )
            
            ax.set_title(
                f'Sensitivity: {mode} (PPV≥{ppv_target:.0f}%)',
                fontsize=14,
                fontweight='bold'
            )
            ax.set_xlabel('Fetal Fraction (%)', fontsize=12)
            ax.set_ylabel('Deletion Length (Mb)', fontsize=12)
            
            # Save
            filename = f'heatmap_ppv{int(ppv_target)}_mode-{mode}.png'
            output_file = heatmaps_dir / filename
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Saved heatmap: {filename}")


def create_line_plots(
    sensitivity_df: pd.DataFrame,
    output_dir: Path
) -> None:
    """Create line plots for each mode and PPV target"""
    
    lineplots_dir = output_dir / 'lineplots'
    lineplots_dir.mkdir(parents=True, exist_ok=True)
    
    for ppv_target in sorted(sensitivity_df['ppv_target'].unique()):
        for mode in sorted(sensitivity_df['mode'].unique()):
            filtered = sensitivity_df[
                (sensitivity_df['ppv_target'] == ppv_target) &
                (sensitivity_df['mode'] == mode)
            ]
            
            if len(filtered) == 0:
                continue
            
            # Create line plot
            fig, ax = plt.subplots(figsize=(10, 6))
            
            ff_values = sorted(filtered['ff'].unique())
            colors = plt.cm.viridis(np.linspace(0, 1, len(ff_values)))
            
            for ff, color in zip(ff_values, colors):
                ff_data = filtered[filtered['ff'] == ff].sort_values('length')
                
                ax.plot(
                    ff_data['length'],
                    ff_data['sensitivity'],
                    marker='o',
                    linewidth=2,
                    markersize=8,
                    label=f'FF={ff:.0f}%',
                    color=color
                )
            
            ax.set_xlabel('Deletion Length (Mb)', fontsize=12)
            ax.set_ylabel('Sensitivity', fontsize=12)
            ax.set_title(
                f'Sensitivity vs. Length: {mode} (PPV≥{ppv_target:.0f}%)',
                fontsize=14,
                fontweight='bold'
            )
            ax.set_ylim(-0.05, 1.05)
            ax.grid(True, alpha=0.3)
            ax.legend(loc='best', fontsize=10)
            
            # Save
            filename = f'lineplot_ppv{int(ppv_target)}_mode-{mode}.png'
            output_file = lineplots_dir / filename
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Saved line plot: {filename}")


def create_ppv_comparison_plots(
    sensitivity_df: pd.DataFrame,
    output_dir: Path
) -> None:
    """Create comparison plots across different PPV targets"""
    
    comparison_dir = output_dir / 'ppv_comparison'
    comparison_dir.mkdir(parents=True, exist_ok=True)
    
    # For each mode, compare PPV targets
    for mode in sorted(sensitivity_df['mode'].unique()):
        mode_data = sensitivity_df[sensitivity_df['mode'] == mode]
        
        # Create comparison heatmaps for each FF
        for ff in sorted(mode_data['ff'].unique()):
            ff_data = mode_data[mode_data['ff'] == ff]
            
            if len(ff_data) == 0:
                continue
            
            # Create pivot for comparison
            pivot = ff_data.pivot_table(
                index='length',
                columns='ppv_target',
                values='sensitivity',
                aggfunc='mean'
            )
            
            # Create heatmap
            fig, ax = plt.subplots(figsize=(10, 8))
            
            sns.heatmap(
                pivot,
                annot=True,
                fmt='.2%',
                cmap='RdYlGn',
                vmin=0,
                vmax=1,
                cbar_kws={'label': 'Sensitivity'},
                ax=ax
            )
            
            ax.set_title(
                f'Sensitivity by PPV Target: {mode}, FF={ff:.0f}%',
                fontsize=14,
                fontweight='bold'
            )
            ax.set_xlabel('PPV Target (%)', fontsize=12)
            ax.set_ylabel('Deletion Length (Mb)', fontsize=12)
            
            # Save
            filename = f'ppv_comparison_{mode}_FF{int(ff)}.png'
            output_file = comparison_dir / filename
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()


def create_summary_report(
    optimal_thresholds_df: pd.DataFrame,
    sensitivity_df: pd.DataFrame,
    output_dir: Path,
    min_detect_length: int
) -> None:
    """Create summary report"""
    
    report_file = output_dir / 'summary_report.txt'
    
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("MICRODELETION DETECTION PERFORMANCE SUMMARY\n")
        f.write("(Automatic Optimal Threshold Selection)\n")
        f.write("="*80 + "\n\n")
        
        f.write("Analysis Parameters:\n")
        f.write("-" * 80 + "\n")
        f.write(f"Min detect length: {min_detect_length/1000000:.1f} Mb\n")
        f.write(f"PPV targets: {sorted(optimal_thresholds_df['ppv_target'].unique())}\n")
        f.write(f"FF values: {sorted(optimal_thresholds_df['ff'].unique())}\n")
        f.write(f"Deletion lengths (Mb): {sorted(optimal_thresholds_df['length'].unique())}\n\n")
        
        f.write("Optimal Z-score Thresholds:\n")
        f.write("-" * 80 + "\n")
        for ppv_target in sorted(optimal_thresholds_df['ppv_target'].unique()):
            f.write(f"\nPPV >= {ppv_target:.0f}%:\n")
            ppv_data = optimal_thresholds_df[
                optimal_thresholds_df['ppv_target'] == ppv_target
            ].sort_values(['ff', 'length'])
            
            for _, row in ppv_data.iterrows():
                if pd.notna(row['optimal_zcut']):
                    f.write(
                        f"  FF={row['ff']:>5.0f}% Length={row['length']:>3.0f}Mb: "
                        f"z={row['optimal_zcut']:>5.1f} "
                        f"(Sens={row['sensitivity']:>6.2%}, PPV={row['ppv']:>6.2%})\n"
                    )
                else:
                    f.write(
                        f"  FF={row['ff']:>5.0f}% Length={row['length']:>3.0f}Mb: "
                        f"NO VALID THRESHOLD\n"
                    )
        
        f.write("\n\nBest Performing Configurations (Sensitivity >= 0.95):\n")
        f.write("-" * 80 + "\n")
        high_sens = sensitivity_df[sensitivity_df['sensitivity'] >= 0.95].sort_values(
            'sensitivity', ascending=False
        )
        
        for _, row in high_sens.head(30).iterrows():
            f.write(
                f"  Mode={row['mode']:<12} FF={row['ff']:>5.0f}% "
                f"Length={row['length']:>3.0f}Mb PPV≥{row['ppv_target']:>3.0f}% "
                f"z={row['zcut']:>5.1f} Sens={row['sensitivity']:>6.2%}\n"
            )
    
    logger.info(f"Saved summary report to {report_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Automatic Microdeletion Detection Performance Analysis"
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
        '--ppv',
        type=str,
        default='90,80,70',
        help='Comma-separated PPV targets (default: 90,80,70)'
    )
    parser.add_argument(
        '--min-detect-length',
        type=int,
        default=1000000,
        help='Minimum detection length in bp (default: 1000000 = 1Mb)'
    )
    parser.add_argument(
        '--outdir',
        type=str,
        required=True,
        help='Output directory for results'
    )
    
    args = parser.parse_args()
    
    # Parse PPV targets
    ppv_targets = [float(p.strip()) for p in args.ppv.split(',')]
    logger.info(f"PPV targets: {ppv_targets}")
    logger.info(f"Min detect length: {args.min_detect_length/1000000:.1f} Mb")
    
    # Create output directory
    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    input_file = Path(args.input)
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        return 1
    
    logger.info(f"Loading data from {input_file}")
    df = pd.read_csv(input_file, sep='\t')
    
    # Convert z-score columns
    for col in ['WC_orig_zscore', 'WC_fetus_zscore', 'WCX_orig_zscore', 'WCX_fetus_zscore']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    sample_dir = Path(args.sample_dir)
    if not sample_dir.exists():
        logger.error(f"Sample directory not found: {sample_dir}")
        return 1
    
    # Step 1: Find optimal z-score thresholds
    logger.info("="*80)
    logger.info("STEP 1: Finding optimal z-score thresholds...")
    logger.info("="*80)
    optimal_thresholds_df = find_optimal_thresholds(
        df, sample_dir, ppv_targets, args.min_detect_length
    )
    
    # Step 2: Calculate sensitivity for all modes
    logger.info("="*80)
    logger.info("STEP 2: Calculating sensitivity for all detection modes...")
    logger.info("="*80)
    sensitivity_df = calculate_sensitivity_by_mode(
        df, sample_dir, optimal_thresholds_df, args.min_detect_length
    )
    
    # Step 3: Save results
    logger.info("="*80)
    logger.info("STEP 3: Saving results...")
    logger.info("="*80)
    save_results(optimal_thresholds_df, sensitivity_df, output_dir)
    
    # Step 4: Create pivot tables
    logger.info("Creating pivot tables...")
    create_pivot_tables(sensitivity_df, output_dir)
    
    # Step 5: Create visualizations
    logger.info("Creating heatmaps...")
    create_heatmaps(sensitivity_df, output_dir)
    
    logger.info("Creating line plots...")
    create_line_plots(sensitivity_df, output_dir)
    
    logger.info("Creating PPV comparison plots...")
    create_ppv_comparison_plots(sensitivity_df, output_dir)
    
    # Step 6: Create summary report
    logger.info("Creating summary report...")
    create_summary_report(
        optimal_thresholds_df, sensitivity_df, output_dir, args.min_detect_length
    )
    
    logger.info("="*80)
    logger.info("Analysis complete!")
    logger.info(f"Results saved to: {output_dir}")
    logger.info("="*80)
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())

