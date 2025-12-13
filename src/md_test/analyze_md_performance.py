#!/usr/bin/env python3
"""
Analyze Microdeletion Detection Performance

This script calculates sensitivity metrics for various detection modes
across different FF (Fetal Fraction) and deletion length combinations.

Usage:
    python analyze_md_performance.py \
        --input zscore_extraction.tsv \
        --outdir results \
        --zcut 3.0,4.0,5.0
"""

import argparse
import logging
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)


def load_data(input_file: Path) -> pd.DataFrame:
    """Load z-score extraction data"""
    logger.info(f"Loading data from {input_file}")
    df = pd.read_csv(input_file, sep='\t')
    
    # Convert FF and length to numeric
    df['ff'] = pd.to_numeric(df['ff'], errors='coerce')
    df['deletion_length_mb'] = pd.to_numeric(df['deletion_length_mb'], errors='coerce')
    
    # Convert z-scores to numeric
    for col in ['WC_orig_zscore', 'WC_fetus_zscore', 'WCX_orig_zscore', 'WCX_fetus_zscore']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    logger.info(f"Loaded {len(df)} samples")
    return df


def calculate_detection_modes(df: pd.DataFrame, zcut: float) -> pd.DataFrame:
    """Calculate 7 detection modes based on z-score threshold
    
    Args:
        df: Input dataframe with z-score columns
        zcut: Z-score threshold (absolute value)
    
    Returns:
        DataFrame with detection mode columns added
    """
    result = df.copy()
    
    # Individual detection modes (absolute value >= zcut)
    result['detect_wc_orig'] = result['WC_orig_zscore'].notna() & (
        result['WC_orig_zscore'].abs() >= zcut
    )
    result['detect_wc_fetus'] = result['WC_fetus_zscore'].notna() & (
        result['WC_fetus_zscore'].abs() >= zcut
    )
    result['detect_wcx_orig'] = result['WCX_orig_zscore'].notna() & (
        result['WCX_orig_zscore'].abs() >= zcut
    )
    result['detect_wcx_fetus'] = result['WCX_fetus_zscore'].notna() & (
        result['WCX_fetus_zscore'].abs() >= zcut
    )
    
    # Aggregated detection modes (OR logic)
    result['detect_orig'] = result['detect_wc_orig'] | result['detect_wcx_orig']
    result['detect_fetus'] = result['detect_wc_fetus'] | result['detect_wcx_fetus']
    result['detect_any'] = (
        result['detect_wc_orig'] | 
        result['detect_wc_fetus'] | 
        result['detect_wcx_orig'] | 
        result['detect_wcx_fetus']
    )
    
    return result


def calculate_sensitivity_matrix(
    df: pd.DataFrame, 
    zcut_list: List[float]
) -> pd.DataFrame:
    """Calculate sensitivity matrix for all combinations
    
    Args:
        df: Input dataframe
        zcut_list: List of z-score thresholds to test
    
    Returns:
        Long format DataFrame with sensitivity metrics
    """
    detection_modes = [
        'wc_orig', 'wc_fetus', 'wcx_orig', 'wcx_fetus',
        'orig', 'fetus', 'any'
    ]
    
    results = []
    
    for zcut in zcut_list:
        logger.info(f"Processing Zcut = {zcut}")
        
        # Calculate detection for this threshold
        df_detected = calculate_detection_modes(df, zcut)
        
        # Remove rows with missing FF or length
        valid_df = df_detected[
            df_detected[['ff', 'deletion_length_mb']].notna().all(axis=1)
        ]
        
        if len(valid_df) == 0:
            logger.warning(f"No valid data for Zcut={zcut}")
            continue
        
        # Group by FF and length
        grouped = valid_df.groupby(['ff', 'deletion_length_mb'])
        
        for (ff, length), group in grouped:
            for mode in detection_modes:
                detect_col = f'detect_{mode}'
                
                if detect_col not in group.columns:
                    continue
                
                # Calculate sensitivity
                detected = group[detect_col].sum()
                total = len(group)
                sensitivity = detected / total if total > 0 else 0.0
                
                results.append({
                    'Zcut': zcut,
                    'FF': ff,
                    'length': int(length),
                    'mode': mode,
                    'sensitivity': sensitivity,
                    'detected': detected,
                    'total': total
                })
    
    return pd.DataFrame(results)


def save_sensitivity_matrix(
    matrix_df: pd.DataFrame, 
    output_dir: Path
) -> None:
    """Save sensitivity matrix in long format"""
    output_file = output_dir / 'md_sensitivity_matrix.csv'
    matrix_df.to_csv(output_file, index=False, float_format='%.4f')
    logger.info(f"Saved sensitivity matrix to {output_file}")


def create_pivot_tables(
    matrix_df: pd.DataFrame, 
    output_dir: Path
) -> None:
    """Create pivot tables (wide format) for each mode and Zcut combination"""
    tables_dir = output_dir / 'tables'
    tables_dir.mkdir(parents=True, exist_ok=True)
    
    modes = matrix_df['mode'].unique()
    zcuts = sorted(matrix_df['Zcut'].unique())
    
    for mode in modes:
        for zcut in zcuts:
            # Filter data
            filtered = matrix_df[
                (matrix_df['mode'] == mode) & 
                (matrix_df['Zcut'] == zcut)
            ]
            
            if len(filtered) == 0:
                continue
            
            # Create pivot table
            pivot = filtered.pivot_table(
                index='FF',
                columns='length',
                values='sensitivity',
                aggfunc='mean'
            )
            
            # Rename columns to include 'Mb' suffix
            pivot.columns = [f'{int(col)}Mb' for col in pivot.columns]
            
            # Save to CSV
            filename = f'table_mode-{mode}_Z-{zcut:.1f}.csv'
            output_file = tables_dir / filename
            pivot.to_csv(output_file, float_format='%.4f')
            
            logger.info(f"Saved pivot table: {filename}")


def create_heatmaps(
    matrix_df: pd.DataFrame, 
    output_dir: Path,
    zcut_list: List[float]
) -> None:
    """Create heatmaps for each detection mode and Zcut"""
    heatmaps_dir = output_dir / 'heatmaps'
    heatmaps_dir.mkdir(parents=True, exist_ok=True)
    
    modes = sorted(matrix_df['mode'].unique())
    
    for zcut in zcut_list:
        for mode in modes:
            # Filter data
            filtered = matrix_df[
                (matrix_df['mode'] == mode) & 
                (matrix_df['Zcut'] == zcut)
            ]
            
            if len(filtered) == 0:
                logger.warning(f"No data for mode={mode}, Zcut={zcut}")
                continue
            
            # Create pivot table for heatmap
            pivot = filtered.pivot_table(
                index='length',
                columns='FF',
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
                f'Detection Sensitivity: {mode} (Zcut={zcut:.1f})',
                fontsize=14,
                fontweight='bold'
            )
            ax.set_xlabel('Fetal Fraction (%)', fontsize=12)
            ax.set_ylabel('Deletion Length (Mb)', fontsize=12)
            
            # Save figure
            filename = f'heatmap_{mode}_Z{zcut:.1f}.png'
            output_file = heatmaps_dir / filename
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Saved heatmap: {filename}")


def create_line_plots(
    matrix_df: pd.DataFrame, 
    output_dir: Path,
    zcut_list: List[float]
) -> None:
    """Create line plots for each detection mode and Zcut"""
    lineplots_dir = output_dir / 'lineplots'
    lineplots_dir.mkdir(parents=True, exist_ok=True)
    
    modes = sorted(matrix_df['mode'].unique())
    
    for zcut in zcut_list:
        for mode in modes:
            # Filter data
            filtered = matrix_df[
                (matrix_df['mode'] == mode) & 
                (matrix_df['Zcut'] == zcut)
            ]
            
            if len(filtered) == 0:
                logger.warning(f"No data for mode={mode}, Zcut={zcut}")
                continue
            
            # Create line plot
            fig, ax = plt.subplots(figsize=(10, 6))
            
            # Get unique FF values
            ff_values = sorted(filtered['FF'].unique())
            colors = plt.cm.viridis(np.linspace(0, 1, len(ff_values)))
            
            for ff, color in zip(ff_values, colors):
                ff_data = filtered[filtered['FF'] == ff].sort_values('length')
                
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
                f'Sensitivity vs. Deletion Length: {mode} (Zcut={zcut:.1f})',
                fontsize=14,
                fontweight='bold'
            )
            ax.set_ylim(-0.05, 1.05)
            ax.grid(True, alpha=0.3)
            ax.legend(loc='best', fontsize=10)
            
            # Save figure
            filename = f'lineplot_{mode}_Z{zcut:.1f}.png'
            output_file = lineplots_dir / filename
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Saved line plot: {filename}")


def create_summary_report(
    matrix_df: pd.DataFrame,
    output_dir: Path
) -> None:
    """Create a summary report with key statistics"""
    report_file = output_dir / 'summary_report.txt'
    
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("MICRODELETION DETECTION PERFORMANCE SUMMARY\n")
        f.write("="*80 + "\n\n")
        
        # Overall statistics
        f.write("Overall Statistics:\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total samples analyzed: {matrix_df['total'].max()}\n")
        f.write(f"Z-score thresholds tested: {sorted(matrix_df['Zcut'].unique())}\n")
        f.write(f"FF values: {sorted(matrix_df['FF'].unique())}\n")
        f.write(f"Deletion lengths (Mb): {sorted(matrix_df['length'].unique())}\n")
        f.write(f"Detection modes: {sorted(matrix_df['mode'].unique())}\n\n")
        
        # Best performing configurations
        f.write("Best Performing Configurations (Sensitivity >= 0.95):\n")
        f.write("-" * 80 + "\n")
        
        high_sens = matrix_df[matrix_df['sensitivity'] >= 0.95].sort_values(
            'sensitivity', ascending=False
        )
        
        for _, row in high_sens.head(20).iterrows():
            f.write(
                f"  Mode={row['mode']:<12} FF={row['FF']:>5.0f}% "
                f"Length={row['length']:>3.0f}Mb Zcut={row['Zcut']:>4.1f} "
                f"Sens={row['sensitivity']:>6.2%}\n"
            )
        
        f.write("\n")
        
        # Mode comparison at fixed conditions
        f.write("Mode Comparison (FF=10%, Length=5Mb, Zcut=3.0):\n")
        f.write("-" * 80 + "\n")
        
        comparison = matrix_df[
            (matrix_df['FF'] == 10) &
            (matrix_df['length'] == 5) &
            (matrix_df['Zcut'] == 3.0)
        ].sort_values('sensitivity', ascending=False)
        
        for _, row in comparison.iterrows():
            f.write(
                f"  {row['mode']:<12} Sensitivity: {row['sensitivity']:>6.2%} "
                f"({row['detected']}/{row['total']})\n"
            )
    
    logger.info(f"Saved summary report to {report_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze Microdeletion Detection Performance"
    )
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Input zscore extraction TSV file'
    )
    parser.add_argument(
        '--outdir',
        type=str,
        required=True,
        help='Output directory for results'
    )
    parser.add_argument(
        '--zcut',
        type=str,
        default='3.0,4.0,5.0',
        help='Comma-separated list of z-score thresholds (default: 3.0,4.0,5.0)'
    )
    
    args = parser.parse_args()
    
    # Parse z-score thresholds
    zcut_list = [float(z.strip()) for z in args.zcut.split(',')]
    logger.info(f"Z-score thresholds: {zcut_list}")
    
    # Create output directory
    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")
    
    # Load data
    input_file = Path(args.input)
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        return 1
    
    df = load_data(input_file)
    
    # Calculate sensitivity matrix
    logger.info("Calculating sensitivity matrix...")
    matrix_df = calculate_sensitivity_matrix(df, zcut_list)
    
    if len(matrix_df) == 0:
        logger.error("No sensitivity data calculated. Check input data.")
        return 1
    
    # Save sensitivity matrix (long format)
    save_sensitivity_matrix(matrix_df, output_dir)
    
    # Create pivot tables (wide format)
    logger.info("Creating pivot tables...")
    create_pivot_tables(matrix_df, output_dir)
    
    # Create visualizations
    logger.info("Creating heatmaps...")
    create_heatmaps(matrix_df, output_dir, zcut_list)
    
    logger.info("Creating line plots...")
    create_line_plots(matrix_df, output_dir, zcut_list)
    
    # Create summary report
    logger.info("Creating summary report...")
    create_summary_report(matrix_df, output_dir)
    
    logger.info("="*80)
    logger.info("Analysis complete!")
    logger.info(f"Results saved to: {output_dir}")
    logger.info("="*80)
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())




