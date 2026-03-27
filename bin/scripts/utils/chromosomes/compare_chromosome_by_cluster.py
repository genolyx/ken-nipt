#!/usr/bin/env python3
"""
두 클러스터 그룹 간 염색체 비율 통계 비교
실제 생물학적 신호(염색체 비율)에서 차이가 있는지 확인
"""

import pandas as pd
import numpy as np
from scipy import stats
import sys
import os

def compare_chromosome_proportions(output_dir):
    """
    두 클러스터 그룹 간 염색체 비율 비교
    """
    clustered_file = os.path.join(output_dir, 'sample_list_with_clusters.tsv')
    batch_file = os.path.join(output_dir, 'sample_list_with_batch.tsv')
    
    if not os.path.exists(clustered_file):
        print(f"❌ File not found: {clustered_file}")
        return
    
    if not os.path.exists(batch_file):
        print(f"❌ File not found: {batch_file}")
        return
    
    # Load cluster info
    df_cluster = pd.read_csv(clustered_file, sep='\t')
    # Load batch info with chromosome data
    df_batch = pd.read_csv(batch_file, sep='\t')
    
    # Merge on sample_id
    df = pd.merge(df_cluster[['sample_id', 'cluster', 'cluster_label']], 
                  df_batch, 
                  on='sample_id', 
                  how='inner')
    
    # Chromosome proportions
    chr_metrics = ['chr13_prop', 'chr18_prop', 'chr21_prop']
    
    # Filter samples with chromosome data
    df_chr = df[chr_metrics + ['cluster', 'cluster_label']].dropna()
    
    print(f"\n{'='*80}")
    print(f"CHROMOSOME PROPORTION COMPARISON BETWEEN CLUSTERS")
    print(f"{'='*80}\n")
    print(f"Total samples with chromosome data: {len(df_chr)}\n")
    
    group1 = df_chr[df_chr['cluster'] == 0]
    group2 = df_chr[df_chr['cluster'] == 1]
    
    print(f"Group 1 (Low cov/dup): {len(group1)} samples")
    print(f"Group 2 (High cov/dup): {len(group2)} samples\n")
    
    print(f"{'='*80}")
    print(f"STATISTICAL COMPARISON")
    print(f"{'='*80}\n")
    
    print(f"{'Chromosome':<15} {'Group 1':>20} {'Group 2':>20} {'Diff(%)':>10} {'p-value':>12} {'Sig':>6}")
    print(f"{'-'*85}")
    
    results = []
    
    for metric in chr_metrics:
        data1 = group1[metric].dropna()
        data2 = group2[metric].dropna()
        
        mean1 = data1.mean()
        std1 = data1.std()
        mean2 = data2.mean()
        std2 = data2.std()
        
        diff_pct = ((mean1 - mean2) / mean2) * 100 if mean2 != 0 else 0
        
        # Mann-Whitney U test
        statistic, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
        
        significant = '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else 'ns'
        
        chr_name = metric.replace('_prop', '').upper()
        
        print(f"{chr_name:<15} {mean1:>10.4f}±{std1:<7.4f} {mean2:>10.4f}±{std2:<7.4f} "
              f"{diff_pct:>9.2f} {p_value:>12.2e} {significant:>6}")
        
        results.append({
            'chromosome': chr_name,
            'group1_mean': mean1,
            'group1_std': std1,
            'group2_mean': mean2,
            'group2_std': std2,
            'diff_pct': diff_pct,
            'p_value': p_value,
            'significant': significant
        })
    
    print(f"\nSignificance: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant\n")
    
    # Also check QC metrics for comparison
    print(f"{'='*80}")
    print(f"QC METRICS COMPARISON (for reference)")
    print(f"{'='*80}\n")
    
    qc_metrics = {
        'duplication_rate(%)': 'Duplication Rate',
        'mean_coverageData(X)': 'Coverage',
        'GC_content(%)': 'GC Content',
        'mapping_rate(%)': 'Mapping Rate'
    }
    
    df_qc = df[list(qc_metrics.keys()) + ['cluster']].dropna()
    group1_qc = df_qc[df_qc['cluster'] == 0]
    group2_qc = df_qc[df_qc['cluster'] == 1]
    
    print(f"{'Metric':<20} {'Group 1':>20} {'Group 2':>20} {'Diff(%)':>10} {'p-value':>12} {'Sig':>6}")
    print(f"{'-'*85}")
    
    for metric, label in qc_metrics.items():
        data1 = group1_qc[metric].dropna()
        data2 = group2_qc[metric].dropna()
        
        if len(data1) == 0 or len(data2) == 0:
            continue
        
        mean1 = data1.mean()
        std1 = data1.std()
        mean2 = data2.mean()
        std2 = data2.std()
        
        diff_pct = ((mean1 - mean2) / mean2) * 100 if mean2 != 0 else 0
        
        statistic, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
        
        significant = '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else 'ns'
        
        if 'coverage' in metric.lower():
            print(f"{label:<20} {mean1:>10.3f}±{std1:<7.3f} {mean2:>10.3f}±{std2:<7.3f} "
                  f"{diff_pct:>9.1f} {p_value:>12.2e} {significant:>6}")
        else:
            print(f"{label:<20} {mean1:>10.2f}±{std1:<7.2f} {mean2:>10.2f}±{std2:<7.2f} "
                  f"{diff_pct:>9.1f} {p_value:>12.2e} {significant:>6}")
    
    # Summary interpretation
    print(f"\n{'='*80}")
    print(f"INTERPRETATION")
    print(f"{'='*80}\n")
    
    # Check if any chromosome proportion is significantly different
    sig_chr = [r for r in results if r['p_value'] < 0.05]
    
    if len(sig_chr) == 0:
        print("✅ **NO significant differences in chromosome proportions (chr13, 18, 21)**")
        print("\nThis means:")
        print("  - The two clusters differ ONLY in technical metrics (coverage/duplication)")
        print("  - The biological signal (chromosome ratios) is UNIFORM across clusters")
        print("  - The batch effect is purely TECHNICAL, not biological")
        print("\n💡 Implication:")
        print("  - Both groups can potentially use the SAME reference")
        print("  - Coverage/duplication differences may not affect aneuploidy detection")
        print("  - The separation in QC-only PCA is a sequencing depth artifact")
    else:
        print(f"⚠️  {len(sig_chr)} chromosome proportion(s) show significant differences:")
        for r in sig_chr:
            print(f"  - {r['chromosome']}: p={r['p_value']:.2e} (diff={r['diff_pct']:.2f}%)")
        print("\nThis suggests the batch effect may have biological implications.")
    
    print(f"\n{'='*80}\n")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python3 compare_chromosome_by_cluster.py <output_dir>")
        sys.exit(1)
    
    compare_chromosome_proportions(sys.argv[1])
