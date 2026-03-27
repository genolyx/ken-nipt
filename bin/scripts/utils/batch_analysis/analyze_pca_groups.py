#!/usr/bin/env python3
"""
PCA 그룹 분리 원인 분석
PC1/PC2 값을 기준으로 샘플을 그룹화하고 각 metric 비교
"""

import pandas as pd
import numpy as np
import sys
import os

def analyze_pca_groups(sample_list_file, output_dir):
    """
    PCA 결과에서 두 그룹을 분리하는 원인 분석
    """
    # Load enriched sample list with batch info
    enriched_file = os.path.join(output_dir, 'sample_list_with_batch.tsv')
    
    if not os.path.exists(enriched_file):
        print(f"❌ File not found: {enriched_file}")
        return
    
    df = pd.read_csv(enriched_file, sep='\t')
    
    # QC metrics
    qc_metrics = ['GC_content(%)', 'mapping_rate(%)', 'duplication_rate(%)', 
                  'mean_mapping_quality', 'mean_coverageData(X)']
    
    # Prepare data for PCA
    df_pca = df[['sample_id', 'batch_id', 'month'] + qc_metrics].dropna()
    
    print(f"\n{'='*80}")
    print(f"PCA GROUP ANALYSIS - QC Metrics")
    print(f"{'='*80}\n")
    print(f"Total samples with complete QC data: {len(df_pca)}\n")
    
    # Standardize
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    
    X = df_pca[qc_metrics].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    
    # Add PC values to dataframe
    df_pca['PC1'] = X_pca[:, 0]
    df_pca['PC2'] = X_pca[:, 1]
    
    print(f"PCA Variance Explained:")
    print(f"  PC1: {pca.explained_variance_ratio_[0]*100:.1f}%")
    print(f"  PC2: {pca.explained_variance_ratio_[1]*100:.1f}%\n")
    
    # Feature loadings
    print(f"Feature Loadings (contribution to PC axes):")
    print(f"{'Metric':<30} {'PC1':>10} {'PC2':>10}")
    print(f"{'-'*52}")
    for i, metric in enumerate(qc_metrics):
        print(f"{metric:<30} {pca.components_[0,i]:>10.3f} {pca.components_[1,i]:>10.3f}")
    print()
    
    # Group by PC1 (positive vs negative)
    print(f"\n{'='*80}")
    print(f"GROUP COMPARISON BY PC1")
    print(f"{'='*80}\n")
    
    group_pos = df_pca[df_pca['PC1'] > 0]
    group_neg = df_pca[df_pca['PC1'] < 0]
    
    print(f"Group 1 (PC1 > 0): {len(group_pos)} samples")
    print(f"Group 2 (PC1 < 0): {len(group_neg)} samples\n")
    
    print(f"{'Metric':<30} {'Group1 (PC1>0)':>20} {'Group2 (PC1<0)':>20} {'Difference':>15}")
    print(f"{'-'*87}")
    
    for metric in qc_metrics:
        mean1 = group_pos[metric].mean()
        mean2 = group_neg[metric].mean()
        std1 = group_pos[metric].std()
        std2 = group_neg[metric].std()
        diff = mean1 - mean2
        diff_pct = (diff / mean2) * 100 if mean2 != 0 else 0
        
        print(f"{metric:<30} {mean1:>10.2f} ± {std1:<6.2f} {mean2:>10.2f} ± {std2:<6.2f} {diff:>10.2f} ({diff_pct:>+6.1f}%)")
    
    # Batch distribution
    print(f"\n{'='*80}")
    print(f"BATCH DISTRIBUTION BY PC1 GROUPS")
    print(f"{'='*80}\n")
    
    print("Top batches in Group 1 (PC1 > 0):")
    batch_counts_pos = group_pos['batch_id'].value_counts().head(10)
    for batch, count in batch_counts_pos.items():
        pct = count / len(group_pos) * 100
        print(f"  {batch}: {count} samples ({pct:.1f}%)")
    
    print("\nTop batches in Group 2 (PC1 < 0):")
    batch_counts_neg = group_neg['batch_id'].value_counts().head(10)
    for batch, count in batch_counts_neg.items():
        pct = count / len(group_neg) * 100
        print(f"  {batch}: {count} samples ({pct:.1f}%)")
    
    # Month distribution
    print(f"\n{'='*80}")
    print(f"MONTH DISTRIBUTION BY PC1 GROUPS")
    print(f"{'='*80}\n")
    
    print("Month distribution in Group 1 (PC1 > 0):")
    month_counts_pos = group_pos['month'].value_counts().sort_index()
    for month, count in month_counts_pos.items():
        pct = count / len(group_pos) * 100
        print(f"  {month}: {count} samples ({pct:.1f}%)")
    
    print("\nMonth distribution in Group 2 (PC1 < 0):")
    month_counts_neg = group_neg['month'].value_counts().sort_index()
    for month, count in month_counts_neg.items():
        pct = count / len(group_neg) * 100
        print(f"  {month}: {count} samples ({pct:.1f}%)")
    
    # Group by PC2 for comparison
    print(f"\n{'='*80}")
    print(f"GROUP COMPARISON BY PC2")
    print(f"{'='*80}\n")
    
    group_pos2 = df_pca[df_pca['PC2'] > 0]
    group_neg2 = df_pca[df_pca['PC2'] < 0]
    
    print(f"Group 1 (PC2 > 0): {len(group_pos2)} samples")
    print(f"Group 2 (PC2 < 0): {len(group_neg2)} samples\n")
    
    print(f"{'Metric':<30} {'Group1 (PC2>0)':>20} {'Group2 (PC2<0)':>20} {'Difference':>15}")
    print(f"{'-'*87}")
    
    for metric in qc_metrics:
        mean1 = group_pos2[metric].mean()
        mean2 = group_neg2[metric].mean()
        std1 = group_pos2[metric].std()
        std2 = group_neg2[metric].std()
        diff = mean1 - mean2
        diff_pct = (diff / mean2) * 100 if mean2 != 0 else 0
        
        print(f"{metric:<30} {mean1:>10.2f} ± {std1:<6.2f} {mean2:>10.2f} ± {std2:<6.2f} {diff:>10.2f} ({diff_pct:>+6.1f}%)")
    
    print(f"\n{'='*80}\n")
    
    # Save analysis results
    report_file = os.path.join(output_dir, 'pca_group_analysis.txt')
    
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("PCA GROUP SEPARATION ANALYSIS\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Total samples: {len(df_pca)}\n")
        f.write(f"PC1 variance: {pca.explained_variance_ratio_[0]*100:.1f}%\n")
        f.write(f"PC2 variance: {pca.explained_variance_ratio_[1]*100:.1f}%\n\n")
        
        f.write("GROUP COMPARISON BY PC1 (Main Separation Axis)\n")
        f.write("-"*80 + "\n")
        f.write(f"Group 1 (PC1 > 0): {len(group_pos)} samples\n")
        f.write(f"Group 2 (PC1 < 0): {len(group_neg)} samples\n\n")
        
        f.write(f"{'Metric':<30} {'Group1':>15} {'Group2':>15} {'Diff(%)':>12}\n")
        f.write("-"*80 + "\n")
        
        for metric in qc_metrics:
            mean1 = group_pos[metric].mean()
            mean2 = group_neg[metric].mean()
            diff_pct = ((mean1 - mean2) / mean2) * 100 if mean2 != 0 else 0
            f.write(f"{metric:<30} {mean1:>15.2f} {mean2:>15.2f} {diff_pct:>11.1f}%\n")
    
    print(f"✓ Saved analysis report: {report_file}\n")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python3 analyze_pca_groups.py <sample_list.tsv> <output_dir>")
        sys.exit(1)
    
    analyze_pca_groups(sys.argv[1], sys.argv[2])
