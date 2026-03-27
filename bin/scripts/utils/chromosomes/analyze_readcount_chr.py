#!/usr/bin/env python3
"""
Read Count vs Chromosome Proportion 분석
Multiplexing artifact인지 Systematic bias인지 구분
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import sys
import os

def analyze_readcount_vs_chrprop(output_dir, lab_name):
    """
    Read count와 염색체 비율의 상관관계 분석
    """
    batch_file = os.path.join(output_dir, 'sample_list_with_batch.tsv')
    cluster_file = os.path.join(output_dir, 'sample_list_with_clusters.tsv')
    
    if not os.path.exists(batch_file) or not os.path.exists(cluster_file):
        print(f"❌ Required files not found")
        return
    
    # Load data
    df_batch = pd.read_csv(batch_file, sep='\t')
    df_cluster = pd.read_csv(cluster_file, sep='\t')
    
    # Merge
    df = pd.merge(df_batch, df_cluster[['sample_id', 'cluster', 'cluster_label']], 
                  on='sample_id', how='inner')
    
    # Key metrics
    chr_metrics = ['chr13_prop', 'chr18_prop', 'chr21_prop']
    
    # Filter samples with all data
    required_cols = ['number_of_mapped_reads', 'mean_coverageData(X)'] + chr_metrics + ['cluster']
    df_analysis = df[required_cols].dropna()
    
    print(f"\n{'='*80}")
    print(f"READ COUNT vs CHROMOSOME PROPORTION ANALYSIS - {lab_name}")
    print(f"{'='*80}\n")
    print(f"Total samples: {len(df_analysis)}\n")
    
    # Convert to numeric
    df_analysis['total_reads'] = pd.to_numeric(df_analysis['number_of_mapped_reads'], errors='coerce')
    df_analysis['coverage'] = pd.to_numeric(df_analysis['mean_coverageData(X)'], errors='coerce')
    
    # Remove any NaN
    df_analysis = df_analysis.dropna()
    
    # Group info
    group1 = df_analysis[df_analysis['cluster'] == 0]
    group2 = df_analysis[df_analysis['cluster'] == 1]
    
    print(f"Group 1: {len(group1)} samples")
    print(f"Group 2: {len(group2)} samples\n")
    
    # ========================================
    # 1. Total Reads Distribution
    # ========================================
    print(f"{'='*80}")
    print(f"TOTAL READS DISTRIBUTION")
    print(f"{'='*80}\n")
    
    print(f"{'Group':<10} {'Mean Reads':>15} {'Median Reads':>15} {'Min':>15} {'Max':>15}")
    print(f"{'-'*70}")
    
    for group_id, label in [(0, 'Group 1'), (1, 'Group 2')]:
        group_data = df_analysis[df_analysis['cluster'] == group_id]
        reads = group_data['total_reads']
        print(f"{label:<10} {reads.mean():>15.0f} {reads.median():>15.0f} "
              f"{reads.min():>15.0f} {reads.max():>15.0f}")
    
    # Statistical test
    reads1 = group1['total_reads']
    reads2 = group2['total_reads']
    stat, p = stats.mannwhitneyu(reads1, reads2)
    print(f"\nMann-Whitney U test: p = {p:.2e}")
    if p < 0.001:
        print("*** Total reads significantly different between groups")
    
    # ========================================
    # 2. Correlation: Total Reads vs Chr Proportion
    # ========================================
    print(f"\n{'='*80}")
    print(f"CORRELATION: Total Reads vs Chromosome Proportions")
    print(f"{'='*80}\n")
    
    print(f"{'Chromosome':<15} {'All Samples':>20} {'Group 1':>20} {'Group 2':>20}")
    print(f"{'-'*77}")
    
    correlations = {}
    
    for chr_metric in chr_metrics:
        chr_name = chr_metric.replace('_prop', '').upper()
        
        # Overall correlation
        corr_all, p_all = stats.spearmanr(df_analysis['total_reads'], 
                                          df_analysis[chr_metric])
        
        # Group-specific correlations
        corr_g1, p_g1 = stats.spearmanr(group1['total_reads'], 
                                        group1[chr_metric])
        corr_g2, p_g2 = stats.spearmanr(group2['total_reads'], 
                                        group2[chr_metric])
        
        print(f"{chr_name:<15} r={corr_all:>6.3f} (p={p_all:.2e}) "
              f"r={corr_g1:>6.3f} (p={p_g1:.2e}) "
              f"r={corr_g2:>6.3f} (p={p_g2:.2e})")
        
        correlations[chr_name] = {
            'all': (corr_all, p_all),
            'group1': (corr_g1, p_g1),
            'group2': (corr_g2, p_g2)
        }
    
    # ========================================
    # 3. Scatter plots
    # ========================================
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    colors = {0: 'steelblue', 1: 'coral'}
    labels = {0: 'Group 1 (Low cov/dup)', 1: 'Group 2 (High cov/dup)'}
    
    for idx, chr_metric in enumerate(chr_metrics):
        chr_name = chr_metric.replace('_prop', '').upper()
        
        # Plot 1: Total reads vs chr proportion
        ax1 = axes[0, idx]
        for cluster in [0, 1]:
            cluster_data = df_analysis[df_analysis['cluster'] == cluster]
            ax1.scatter(cluster_data['total_reads'] / 1e6, 
                       cluster_data[chr_metric] * 100,
                       c=colors[cluster], 
                       label=labels[cluster],
                       alpha=0.5, 
                       s=20)
        
        ax1.set_xlabel('Total Reads (Million)', fontsize=11)
        ax1.set_ylabel(f'{chr_name} Proportion (%)', fontsize=11)
        ax1.set_title(f'{chr_name} vs Total Reads', fontsize=12, fontweight='bold')
        ax1.legend(fontsize=9)
        ax1.grid(alpha=0.3)
        
        # Add correlation text
        corr, p = correlations[chr_name]['all']
        ax1.text(0.05, 0.95, f'r={corr:.3f}\np={p:.2e}', 
                transform=ax1.transAxes, fontsize=9,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Plot 2: Coverage vs chr proportion
        ax2 = axes[1, idx]
        for cluster in [0, 1]:
            cluster_data = df_analysis[df_analysis['cluster'] == cluster]
            ax2.scatter(cluster_data['coverage'], 
                       cluster_data[chr_metric] * 100,
                       c=colors[cluster], 
                       label=labels[cluster],
                       alpha=0.5, 
                       s=20)
        
        ax2.set_xlabel('Coverage (X)', fontsize=11)
        ax2.set_ylabel(f'{chr_name} Proportion (%)', fontsize=11)
        ax2.set_title(f'{chr_name} vs Coverage', fontsize=12, fontweight='bold')
        ax2.legend(fontsize=9)
        ax2.grid(alpha=0.3)
    
    plt.suptitle(f'Multiplexing Effect: Read Count/Coverage vs Chromosome Proportions - {lab_name}',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'readcount_vs_chrprop.png'), dpi=150)
    print(f"\n✓ Saved: readcount_vs_chrprop.png")
    plt.close()
    
    # ========================================
    # 4. Interpretation
    # ========================================
    print(f"\n{'='*80}")
    print(f"INTERPRETATION")
    print(f"{'='*80}\n")
    
    # Check if chr proportions correlate with read count
    strong_corr_count = 0
    for chr_name, corr_data in correlations.items():
        corr_all, p_all = corr_data['all']
        if abs(corr_all) > 0.3 and p_all < 0.001:
            strong_corr_count += 1
            print(f"⚠️  {chr_name}: Strong correlation with total reads (r={corr_all:.3f})")
    
    if strong_corr_count == 0:
        print("✅ Chromosome proportions do NOT correlate with total reads")
        print("\nThis means:")
        print("  - The chr proportion differences are NOT due to multiplexing")
        print("  - The differences are SYSTEMATIC between the two groups")
        print("  - This is a REAL bias, not a sampling artifact")
    else:
        print(f"\n⚠️  {strong_corr_count} chromosome(s) correlate with total reads")
        print("\nThis suggests:")
        print("  - Some chr proportion variation may be due to sampling")
        print("  - But systematic differences between groups still exist")
    
    # Check read count difference
    read_ratio = reads1.median() / reads2.median()
    print(f"\nMedian read count ratio (Group1/Group2): {read_ratio:.2f}x")
    
    if abs(read_ratio - 1.0) > 0.5:
        print("  → Groups have very different sequencing depths")
        print("  → Multiplexing effect is present")
    else:
        print("  → Groups have similar sequencing depths")
        print("  → Multiplexing is NOT the main factor")
    
    print(f"\n{'='*80}\n")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python3 analyze_readcount_chr.py <output_dir> <lab_name>")
        sys.exit(1)
    
    analyze_readcount_vs_chrprop(sys.argv[1], sys.argv[2])
