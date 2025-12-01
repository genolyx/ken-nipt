#!/usr/bin/env python3
"""
ROC Analysis for Microdeletion Detection (One-vs-Rest)

This script performs One-vs-Rest ROC analysis where each disease is treated
as positive class and the other 7 diseases are treated as negative class.

Usage:
    python roc_analysis.py \
        --data_dir /data/md_validation/analysis_result \
        --outdir roc_results
"""

import argparse
import logging
import json
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple
import re

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

# Define diseases
DISEASES = ['1p36', '2q33', 'CDC', 'DGS', 'Jacobsen', 'PWS', 'WBS', 'WHS']

# Define disease to chromosome mapping
DISEASE_CHROMOSOMES = {
    '1p36': '1',
    '2q33': '2',
    'CDC': '5',       # Cri-du-chat syndrome
    'DGS': '22',      # DiGeorge syndrome
    'Jacobsen': '11', # Jacobsen syndrome
    'PWS': '15',      # Prader-Willi syndrome
    'WBS': '7',       # Williams-Beuren syndrome
    'WHS': '4'        # Wolf-Hirschhorn syndrome
}

# Define detection modes
INDIVIDUAL_MODES = ['wc_orig', 'wc_fetus', 'wcx_orig', 'wcx_fetus']
GROUP_MODES = ['orig', 'fetus', 'any']  # Combined modes
ALL_MODES = INDIVIDUAL_MODES + GROUP_MODES


def extract_zscore_from_aberrations_bed(bed_file: Path, target_chr: str) -> float:
    """Extract maximum absolute z-score for specific chromosome from aberrations.bed
    
    Args:
        bed_file: Path to aberrations.bed file
        target_chr: Target chromosome (e.g., '1', '22', '7')
        
    Returns:
        Maximum absolute z-score value for target chromosome, 0 if not detected
    """
    try:
        if not bed_file.exists():
            return 0.0
            
        df = pd.read_csv(bed_file, sep='\t')
        
        if len(df) == 0:
            return 0.0
        
        # Filter by target chromosome
        target_df = df[df['chr'].astype(str).isin([target_chr, f'chr{target_chr}'])]
        
        if len(target_df) == 0:
            return 0.0
        
        # Return maximum absolute z-score
        return target_df['zscore'].abs().max()
        
    except Exception as e:
        logger.debug(f"Error reading {bed_file}: {e}")
        return 0.0


def extract_zscore_for_chromosome(report_file: Path, target_chr: str) -> float:
    """Extract maximum absolute z-score for specific chromosome from Test results
    
    Args:
        report_file: Path to report txt file
        target_chr: Target chromosome (e.g., '1', '22', '7')
        
    Returns:
        Maximum absolute z-score value for target chromosome, 0 if not detected
    """
    try:
        with open(report_file, 'r') as f:
            content = f.read()
            
        # Find "Test results:" section
        test_results_match = re.search(
            r'# Test results: #\s*z-score\s+effect\s+mbsize\s+location\s*\n(.*?)(?:\n\n|\Z)',
            content,
            re.DOTALL
        )
        
        if not test_results_match:
            # No detections at all - return 0
            return 0.0
        
        results_section = test_results_match.group(1)
        
        # Extract z-scores from the target chromosome
        zscores = []
        for line in results_section.strip().split('\n'):
            if line.strip():
                parts = line.split()
                if len(parts) >= 4:  # z-score, effect, mbsize, location
                    try:
                        zscore = float(parts[0])
                        location = parts[3]  # e.g., "1:800000-10000000"
                        
                        # Extract chromosome from location
                        chrom = location.split(':')[0]
                        
                        # Check if this is the target chromosome
                        if chrom == target_chr or chrom == f'chr{target_chr}':
                            zscores.append(abs(zscore))
                    except (ValueError, IndexError):
                        continue
        
        if zscores:
            # Return the maximum absolute z-score for this chromosome
            return max(zscores)
        else:
            # No detection on target chromosome
            return 0.0
            
    except Exception as e:
        logger.warning(f"Error reading {report_file}: {e}")
        return 0.0


def extract_sample_metadata(sample_dir: Path) -> Dict:
    """Extract metadata from sample directory
    
    Args:
        sample_dir: Path to sample directory
        
    Returns:
        Dictionary with metadata (ff, length_mb, gender, etc.)
    """
    sample_name = sample_dir.name
    
    # Parse sample name: e.g., 1_1_1p36deletionsyndrome_FF10_15M_10Mb_F
    parts = sample_name.split('_')
    
    metadata = {
        'sample_id': sample_name,
        'sample_dir': str(sample_dir)
    }
    
    # Extract FF
    ff_match = re.search(r'FF(\d+)', sample_name)
    if ff_match:
        metadata['ff'] = int(ff_match.group(1))
    
    # Extract deletion length
    length_match = re.search(r'(\d+)Mb', sample_name)
    if length_match:
        metadata['deletion_length_mb'] = int(length_match.group(1))
    
    # Extract gender
    if sample_name.endswith('_F'):
        metadata['gender'] = 'F'
    elif sample_name.endswith('_M'):
        metadata['gender'] = 'M'
    
    # Try to read JSON metadata if available
    json_files = list(sample_dir.glob('*.json'))
    if json_files:
        try:
            with open(json_files[0], 'r') as f:
                json_data = json.load(f)
                if 'calculated_ff' in json_data:
                    metadata['calculated_ff'] = json_data['calculated_ff'].get('final_ff')
                if 'deletion' in json_data:
                    metadata['deletion_size_mb'] = json_data['deletion'].get('size_mb')
        except Exception as e:
            logger.debug(f"Could not read JSON from {json_files[0]}: {e}")
    
    return metadata


def collect_all_samples_for_target(
    data_dir: Path, 
    target_disease: str
) -> pd.DataFrame:
    """Collect z-scores for specific target disease from all samples (wide format)
    
    For ROC analysis, we extract z-scores for the target disease's chromosome
    from all samples (both positive and negative).
    
    Returns wide format with one row per sample and columns for each mode.
    
    Args:
        data_dir: Root directory containing disease subdirectories
        target_disease: Target disease to extract z-scores for
        
    Returns:
        DataFrame with columns: disease, sample_id, wc_orig, wc_fetus, wcx_orig, wcx_fetus, etc.
    """
    all_data = []
    target_chr = DISEASE_CHROMOSOMES[target_disease]
    
    for disease in DISEASES:
        disease_dir = data_dir / disease
        if not disease_dir.exists():
            logger.warning(f"Disease directory not found: {disease_dir}")
            continue
        
        # Get all sample directories
        sample_dirs = [d for d in disease_dir.iterdir() if d.is_dir()]
        
        for sample_dir in sample_dirs:
            # Extract metadata
            metadata = extract_sample_metadata(sample_dir)
            metadata['disease'] = disease
            metadata['target_disease'] = target_disease
            
            # Extract z-scores from each mode for target chromosome
            results_dir = sample_dir / 'results'
            if not results_dir.exists():
                continue
            
            # Collect z-scores for all modes
            zscores = {}
            for mode in INDIVIDUAL_MODES:
                mode_parts = mode.split('_')
                wc_type = mode_parts[0]  # 'wc' or 'wcx'
                correction = mode_parts[1]  # 'orig' or 'fetus'
                
                if wc_type == 'wc':
                    # WC uses report.txt files
                    report_file = results_dir / f'{wc_type}_{correction}_report.txt'
                    
                    if report_file.exists():
                        zscore = extract_zscore_for_chromosome(report_file, target_chr)
                        zscores[mode] = zscore
                    else:
                        zscores[mode] = 0.0
                        
                elif wc_type == 'wcx':
                    # WCX uses aberrations.bed files
                    bed_file = results_dir / f'{wc_type}_{correction}_aberrations.bed'
                    
                    if bed_file.exists():
                        zscore = extract_zscore_from_aberrations_bed(bed_file, target_chr)
                        zscores[mode] = zscore
                    else:
                        zscores[mode] = 0.0
            
            # Add all z-scores to the record
            record = metadata.copy()
            record.update(zscores)
            all_data.append(record)
    
    df = pd.DataFrame(all_data)
    return df


def calculate_roc_curve_manual(y_true: np.ndarray, y_score: np.ndarray) -> Tuple:
    """Calculate ROC curve manually without sklearn
    
    Args:
        y_true: True binary labels (0 or 1)
        y_score: Predicted scores (higher = more likely positive)
        
    Returns:
        Tuple of (fpr, tpr, thresholds)
    """
    # Get unique thresholds (sorted in descending order)
    thresholds = np.sort(np.unique(y_score))[::-1]
    
    # Add extreme thresholds
    thresholds = np.concatenate([[np.inf], thresholds, [-np.inf]])
    
    # Calculate TPR and FPR for each threshold
    tpr_list = []
    fpr_list = []
    
    n_positive = np.sum(y_true == 1)
    n_negative = np.sum(y_true == 0)
    
    for threshold in thresholds:
        # Predict positive if score >= threshold
        y_pred = (y_score >= threshold).astype(int)
        
        # Calculate TP, FP, TN, FN
        tp = np.sum((y_true == 1) & (y_pred == 1))
        fp = np.sum((y_true == 0) & (y_pred == 1))
        
        # Calculate TPR and FPR
        tpr = tp / n_positive if n_positive > 0 else 0
        fpr = fp / n_negative if n_negative > 0 else 0
        
        tpr_list.append(tpr)
        fpr_list.append(fpr)
    
    return np.array(fpr_list), np.array(tpr_list), thresholds


def calculate_auc(fpr: np.ndarray, tpr: np.ndarray) -> float:
    """Calculate AUC using trapezoidal rule
    
    Args:
        fpr: False positive rates
        tpr: True positive rates
        
    Returns:
        AUC value
    """
    # Sort by fpr
    sorted_indices = np.argsort(fpr)
    fpr_sorted = fpr[sorted_indices]
    tpr_sorted = tpr[sorted_indices]
    
    # Calculate AUC using trapezoidal rule
    auc_value = np.trapz(tpr_sorted, fpr_sorted)
    
    return auc_value


def calculate_roc_one_vs_rest(
    df: pd.DataFrame,
    target_disease: str,
    mode: str,
    threshold_range: np.ndarray = None
) -> Dict:
    """Calculate ROC curve for one disease vs rest
    
    The df should be in wide format with z-scores for all modes as columns.
    Supports both individual modes and group modes (OR logic).
    
    Args:
        df: DataFrame with wide format (one row per sample)
        target_disease: Target disease (positive class)
        mode: Detection mode (individual or group)
        threshold_range: Array of thresholds to test
        
    Returns:
        Dictionary with fpr, tpr, thresholds, auc
    """
    # Create binary labels: 1 if sample's disease matches target, 0 otherwise
    y_true = (df['disease'] == target_disease).astype(int).values
    
    # Get z-scores based on mode
    if mode in INDIVIDUAL_MODES:
        # Individual mode: use that mode's z-score directly
        if mode not in df.columns:
            logger.warning(f"Mode {mode} not found in data")
            return None
        y_score = df[mode].values
        
    elif mode == 'orig':
        # orig group: max of (wc_orig, wcx_orig)
        available_cols = [c for c in ['wc_orig', 'wcx_orig'] if c in df.columns]
        if not available_cols:
            logger.warning(f"No orig modes found in data")
            return None
        y_score = df[available_cols].max(axis=1).values
        
    elif mode == 'fetus':
        # fetus group: max of (wc_fetus, wcx_fetus)
        available_cols = [c for c in ['wc_fetus', 'wcx_fetus'] if c in df.columns]
        if not available_cols:
            logger.warning(f"No fetus modes found in data")
            return None
        y_score = df[available_cols].max(axis=1).values
        
    elif mode == 'any':
        # any group: max of all 4 modes
        available_cols = [c for c in INDIVIDUAL_MODES if c in df.columns]
        if not available_cols:
            logger.warning(f"No individual modes found in data")
            return None
        y_score = df[available_cols].max(axis=1).values
        
    else:
        logger.warning(f"Unknown mode: {mode}")
        return None
    
    # Check if we have both classes
    if len(np.unique(y_true)) < 2:
        logger.warning(f"Only one class present for {target_disease}, mode={mode}")
        return None
    
    # Calculate ROC curve
    fpr, tpr, thresholds = calculate_roc_curve_manual(y_true, y_score)
    roc_auc = calculate_auc(fpr, tpr)
    
    return {
        'fpr': fpr,
        'tpr': tpr,
        'thresholds': thresholds,
        'auc': roc_auc,
        'n_positive': int(y_true.sum()),
        'n_negative': int((1 - y_true).sum())
    }


def plot_roc_curve_single(
    roc_data: Dict,
    target_disease: str,
    mode: str,
    output_file: Path
) -> None:
    """Plot ROC curve for single disease vs rest
    
    Args:
        roc_data: Dictionary from calculate_roc_one_vs_rest
        target_disease: Target disease name
        mode: Detection mode
        output_file: Output file path
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    
    ax.plot(
        roc_data['fpr'],
        roc_data['tpr'],
        color='darkorange',
        lw=2,
        label=f'ROC curve (AUC = {roc_data["auc"]:.4f})'
    )
    ax.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Random')
    
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate (Sensitivity)', fontsize=12)
    ax.set_title(
        f'ROC Curve: {target_disease} vs Rest\n'
        f'Mode: {mode} (n_pos={roc_data["n_positive"]}, n_neg={roc_data["n_negative"]})',
        fontsize=14,
        fontweight='bold'
    )
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


def plot_roc_curves_all_diseases(
    all_roc_data: Dict,
    mode: str,
    output_file: Path
) -> None:
    """Plot ROC curves for all diseases on one plot
    
    Args:
        all_roc_data: Dictionary of {disease: roc_data}
        mode: Detection mode
        output_file: Output file path
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(DISEASES)))
    
    for disease, color in zip(DISEASES, colors):
        if disease not in all_roc_data or all_roc_data[disease] is None:
            continue
        
        roc_data = all_roc_data[disease]
        ax.plot(
            roc_data['fpr'],
            roc_data['tpr'],
            color=color,
            lw=2,
            label=f'{disease} (AUC={roc_data["auc"]:.3f})'
        )
    
    ax.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--', alpha=0.5)
    
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate (Sensitivity)', fontsize=12)
    ax.set_title(
        f'ROC Curves: All Diseases (One-vs-Rest)\nMode: {mode}',
        fontsize=14,
        fontweight='bold'
    )
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


def plot_roc_curves_all_modes(
    disease_mode_roc: Dict,
    target_disease: str,
    output_file: Path
) -> None:
    """Plot ROC curves for all modes for a single disease
    
    Args:
        disease_mode_roc: Dictionary of {mode: roc_data} for target disease
        target_disease: Target disease name
        output_file: Output file path
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Filter out None values
    valid_modes = [m for m in disease_mode_roc.keys() if disease_mode_roc[m] is not None]
    
    if not valid_modes:
        logger.warning(f"No valid modes for {target_disease}")
        return
    
    colors = plt.cm.Set2(np.linspace(0, 1, len(valid_modes)))
    
    for mode, color in zip(valid_modes, colors):
        roc_data = disease_mode_roc[mode]
        ax.plot(
            roc_data['fpr'],
            roc_data['tpr'],
            color=color,
            lw=2,
            label=f'{mode} (AUC={roc_data["auc"]:.4f})'
        )
    
    ax.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--', alpha=0.5)
    
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate (Sensitivity)', fontsize=12)
    ax.set_title(
        f'ROC Curves: {target_disease} (All Modes)',
        fontsize=14,
        fontweight='bold'
    )
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


def create_auc_heatmap(
    auc_matrix: pd.DataFrame,
    output_file: Path
) -> None:
    """Create heatmap of AUC values
    
    Args:
        auc_matrix: DataFrame with diseases as rows, modes as columns
        output_file: Output file path
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    sns.heatmap(
        auc_matrix,
        annot=True,
        fmt='.4f',
        cmap='RdYlGn',
        vmin=0.5,
        vmax=1.0,
        cbar_kws={'label': 'AUC'},
        ax=ax
    )
    
    ax.set_title(
        'AUC Scores: Disease Detection (One-vs-Rest)',
        fontsize=14,
        fontweight='bold'
    )
    ax.set_xlabel('Detection Mode', fontsize=12)
    ax.set_ylabel('Disease', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


def save_roc_data(
    all_roc_data: Dict,
    output_dir: Path
) -> None:
    """Save ROC data to CSV files
    
    Args:
        all_roc_data: Nested dict {disease: {mode: roc_data}}
        output_dir: Output directory
    """
    data_dir = output_dir / 'roc_data'
    data_dir.mkdir(parents=True, exist_ok=True)
    
    for disease, mode_data in all_roc_data.items():
        for mode, roc_data in mode_data.items():
            if roc_data is None:
                continue
            
            # Save ROC curve points
            df = pd.DataFrame({
                'fpr': roc_data['fpr'],
                'tpr': roc_data['tpr'],
                'threshold': roc_data['thresholds']
            })
            
            output_file = data_dir / f'roc_{disease}_{mode}.csv'
            df.to_csv(output_file, index=False, float_format='%.6f')
            logger.info(f"Saved ROC data: {output_file.name}")


def create_summary_report(
    all_roc_data: Dict,
    auc_matrix: pd.DataFrame,
    output_dir: Path
) -> None:
    """Create summary report
    
    Args:
        all_roc_data: Nested dict {disease: {mode: roc_data}}
        auc_matrix: DataFrame with AUC values
        output_dir: Output directory
    """
    report_file = output_dir / 'roc_summary.txt'
    
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("ROC ANALYSIS SUMMARY (One-vs-Rest)\n")
        f.write("="*80 + "\n\n")
        
        f.write("AUC Scores by Disease and Mode:\n")
        f.write("-" * 80 + "\n")
        f.write(auc_matrix.to_string())
        f.write("\n\n")
        
        f.write("Best AUC by Disease:\n")
        f.write("-" * 80 + "\n")
        for disease in DISEASES:
            if disease in all_roc_data:
                best_mode = auc_matrix.loc[disease].idxmax()
                best_auc = auc_matrix.loc[disease].max()
                f.write(f"{disease:<15} Best Mode: {best_mode:<15} AUC: {best_auc:.4f}\n")
        f.write("\n")
        
        f.write("Best AUC by Mode:\n")
        f.write("-" * 80 + "\n")
        for mode in auc_matrix.columns:
            best_disease = auc_matrix[mode].idxmax()
            best_auc = auc_matrix[mode].max()
            f.write(f"{mode:<15} Best Disease: {best_disease:<15} AUC: {best_auc:.4f}\n")
        f.write("\n")
        
        f.write("Overall Best:\n")
        f.write("-" * 80 + "\n")
        best_disease = auc_matrix.max(axis=1).idxmax()
        best_mode = auc_matrix.loc[best_disease].idxmax()
        best_auc = auc_matrix.loc[best_disease, best_mode]
        f.write(f"Disease: {best_disease}, Mode: {best_mode}, AUC: {best_auc:.4f}\n")
    
    logger.info(f"Saved summary report: {report_file}")


def main():
    parser = argparse.ArgumentParser(
        description="ROC Analysis for Microdeletion Detection (One-vs-Rest)"
    )
    parser.add_argument(
        '--data_dir',
        type=str,
        required=True,
        help='Root directory containing disease subdirectories'
    )
    parser.add_argument(
        '--outdir',
        type=str,
        required=True,
        help='Output directory for results'
    )
    
    args = parser.parse_args()
    
    data_dir = Path(args.data_dir)
    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("="*80)
    logger.info("Starting ROC Analysis (One-vs-Rest)")
    logger.info("="*80)
    
    # Step 1 & 2: Collect data and Calculate ROC for each target disease
    logger.info("Step 1 & 2: Collecting data and calculating ROC curves...")
    all_roc_data = {}
    auc_data = []
    all_collected_data = []
    
    for target_disease in DISEASES:
        logger.info(f"Processing target disease: {target_disease} (chr{DISEASE_CHROMOSOMES[target_disease]})")
        
        # Collect z-scores for this target disease from all samples
        df = collect_all_samples_for_target(data_dir, target_disease)
        
        if len(df) == 0:
            logger.warning(f"No data collected for {target_disease}")
            continue
        
        # Save this target's collected data
        all_collected_data.append(df)
        
        all_roc_data[target_disease] = {}
        
        for mode in ALL_MODES:
            roc_data = calculate_roc_one_vs_rest(df, target_disease, mode)
            all_roc_data[target_disease][mode] = roc_data
            
            if roc_data is not None:
                auc_data.append({
                    'disease': target_disease,
                    'mode': mode,
                    'auc': roc_data['auc'],
                    'n_positive': roc_data['n_positive'],
                    'n_negative': roc_data['n_negative']
                })
                logger.info(f"  {mode}: AUC = {roc_data['auc']:.4f}")
    
    # Save all collected data
    if all_collected_data:
        combined_df = pd.concat(all_collected_data, ignore_index=True)
        combined_df.to_csv(output_dir / 'collected_data.csv', index=False)
        logger.info(f"Saved collected data to: {output_dir / 'collected_data.csv'}")
    
    # Create AUC matrix
    auc_df = pd.DataFrame(auc_data)
    auc_matrix = auc_df.pivot(index='disease', columns='mode', values='auc')
    
    # Ensure column order (only for modes that exist)
    available_modes = [m for m in ALL_MODES if m in auc_matrix.columns]
    auc_matrix = auc_matrix[available_modes]
    
    # Save AUC matrix
    auc_matrix.to_csv(output_dir / 'auc_matrix.csv', float_format='%.4f')
    logger.info(f"Saved AUC matrix to: {output_dir / 'auc_matrix.csv'}")
    
    # Step 3: Create visualizations
    logger.info("Step 3: Creating visualizations...")
    
    # Individual ROC curves
    individual_dir = output_dir / 'individual_roc'
    individual_dir.mkdir(parents=True, exist_ok=True)
    
    for disease in DISEASES:
        for mode in available_modes:
            if mode in all_roc_data[disease] and all_roc_data[disease][mode] is not None:
                output_file = individual_dir / f'roc_{disease}_{mode}.png'
                plot_roc_curve_single(
                    all_roc_data[disease][mode],
                    disease,
                    mode,
                    output_file
                )
    
    logger.info(f"Saved individual ROC curves to: {individual_dir}")
    
    # All diseases comparison (one plot per mode)
    comparison_dir = output_dir / 'comparison_roc'
    comparison_dir.mkdir(parents=True, exist_ok=True)
    
    for mode in available_modes:
        mode_roc_data = {d: all_roc_data[d][mode] for d in DISEASES}
        output_file = comparison_dir / f'roc_all_diseases_{mode}.png'
        plot_roc_curves_all_diseases(mode_roc_data, mode, output_file)
    
    logger.info(f"Saved comparison ROC curves to: {comparison_dir}")
    
    # All modes comparison (one plot per disease)
    modes_dir = output_dir / 'modes_comparison'
    modes_dir.mkdir(parents=True, exist_ok=True)
    
    for disease in DISEASES:
        output_file = modes_dir / f'roc_{disease}_all_modes.png'
        plot_roc_curves_all_modes(all_roc_data[disease], disease, output_file)
    
    logger.info(f"Saved mode comparison plots to: {modes_dir}")
    
    # AUC heatmap
    heatmap_file = output_dir / 'auc_heatmap.png'
    create_auc_heatmap(auc_matrix, heatmap_file)
    logger.info(f"Saved AUC heatmap to: {heatmap_file}")
    
    # Step 4: Save ROC data
    logger.info("Step 4: Saving ROC data...")
    save_roc_data(all_roc_data, output_dir)
    
    # Step 5: Create summary report
    logger.info("Step 5: Creating summary report...")
    create_summary_report(all_roc_data, auc_matrix, output_dir)
    
    logger.info("="*80)
    logger.info("ROC Analysis Complete!")
    logger.info(f"Results saved to: {output_dir}")
    logger.info("="*80)
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())

