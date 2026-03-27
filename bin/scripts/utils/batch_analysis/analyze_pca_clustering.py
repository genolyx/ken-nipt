#!/usr/bin/env python3
"""
PCA 기반 클러스터링 분석 - 두 그룹 구분 및 특성 비교
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy import stats
import sys
import os

def analyze_pca_clustering(sample_list_file, output_dir, lab_name):
    """
    PCA에서 보이는 두 그룹을 clustering으로 구분하고 분석
    """
    # Load enriched sample list
    enriched_file = os.path.join(output_dir, 'sample_list_with_batch.tsv')
    
    if not os.path.exists(enriched_file):
        print(f"❌ File not found: {enriched_file}")
        return
    
    df = pd.read_csv(enriched_file, sep='\t')
    
    # QC metrics for PCA
    qc_metrics = ['GC_content(%)', 'mapping_rate(%)', 'duplication_rate(%)', 
                  'mean_mapping_quality', 'mean_coverageData(X)']
    
    # Prepare data
    df_pca = df[['sample_id', 'batch_id', 'month'] + qc_metrics].dropna()
    
    print(f"\n{'='*80}")
    print(f"PCA CLUSTERING ANALYSIS - {lab_name}")
    print(f"{'='*80}\n")
    print(f"Total samples: {len(df_pca)}\n")
    
    # Standardize and perform PCA
    X = df_pca[qc_metrics].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    
    df_pca['PC1'] = X_pca[:, 0]
    df_pca['PC2'] = X_pca[:, 1]
    
    print(f"PCA Variance Explained:")
    print(f"  PC1: {pca.explained_variance_ratio_[0]*100:.1f}%")
    print(f"  PC2: {pca.explained_variance_ratio_[1]*100:.1f}%\n")
    
    # K-means clustering (k=2)
    kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
    df_pca['cluster'] = kmeans.fit_predict(X_pca)
    
    # Relabel clusters by PC1 mean (to make Group 1 = lower PC1, Group 2 = higher PC1)
    pc1_mean_0 = df_pca[df_pca['cluster'] == 0]['PC1'].mean()
    pc1_mean_1 = df_pca[df_pca['cluster'] == 1]['PC1'].mean()
    
    if pc1_mean_0 > pc1_mean_1:
        df_pca['cluster'] = df_pca['cluster'].map({0: 1, 1: 0})
    
    df_pca['cluster_label'] = df_pca['cluster'].map({0: 'Group 1 (Left)', 1: 'Group 2 (Right)'})
    
    # Count samples per cluster
    cluster_counts = df_pca['cluster_label'].value_counts()
    print(f"Cluster Sizes:")
    for label, count in cluster_counts.items():
        pct = count / len(df_pca) * 100
        print(f"  {label}: {count} samples ({pct:.1f}%)")
    print()
    
    # ======================
    # 1. PCA Plot with Clusters
    # ======================
    fig, ax = plt.subplots(figsize=(12, 8))
    
    for cluster_id, label in enumerate(['Group 1 (Left)', 'Group 2 (Right)']):
        cluster_data = df_pca[df_pca['cluster'] == cluster_id]
        color = 'steelblue' if cluster_id == 0 else 'coral'
        ax.scatter(cluster_data['PC1'], cluster_data['PC2'], 
                  c=color, label=label, alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
    
    # Plot cluster centers
    centers = kmeans.cluster_centers_
    ax.scatter(centers[:, 0], centers[:, 1], c='red', marker='X', s=300, 
              edgecolors='black', linewidth=2, label='Cluster Centers', zorder=10)
    
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=13)
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=13)
    ax.set_title(f'PCA with K-means Clustering (k=2) - {lab_name}', fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=11)
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pca_clustering_groups.png'), dpi=150)
    print(f"✓ Saved: pca_clustering_groups.png\n")
    plt.close()
    
    # ======================
    # 2. Batch Distribution by Cluster
    # ======================
    print(f"{'='*80}")
    print(f"BATCH COMPOSITION BY CLUSTER")
    print(f"{'='*80}\n")
    
    for cluster_id, label in enumerate(['Group 1 (Left)', 'Group 2 (Right)']):
        cluster_data = df_pca[df_pca['cluster'] == cluster_id]
        print(f"\n{label} - Top 10 Batches:")
        print(f"{'-'*60}")
        
        batch_counts = cluster_data['batch_id'].value_counts().head(10)
        for batch, count in batch_counts.items():
            pct = count / len(cluster_data) * 100
            print(f"  {batch:30s}: {count:4d} samples ({pct:5.1f}%)")
    
    # Month distribution
    print(f"\n{'='*80}")
    print(f"MONTH DISTRIBUTION BY CLUSTER")
    print(f"{'='*80}\n")
    
    for cluster_id, label in enumerate(['Group 1 (Left)', 'Group 2 (Right)']):
        cluster_data = df_pca[df_pca['cluster'] == cluster_id]
        print(f"\n{label}:")
        print(f"{'-'*60}")
        
        month_counts = cluster_data['month'].value_counts().sort_index()
        for month, count in month_counts.items():
            pct = count / len(cluster_data) * 100
            print(f"  {month}: {count:4d} samples ({pct:5.1f}%)")
    
    # ======================
    # 3. Statistical Comparison of QC Metrics
    # ======================
    print(f"\n{'='*80}")
    print(f"STATISTICAL COMPARISON OF QC METRICS BETWEEN CLUSTERS")
    print(f"{'='*80}\n")
    
    group1 = df_pca[df_pca['cluster'] == 0]
    group2 = df_pca[df_pca['cluster'] == 1]
    
    print(f"{'Metric':<30} {'Group 1':>15} {'Group 2':>15} {'p-value':>12} {'Significant':>12}")
    print(f"{'-'*86}")
    
    stats_results = []
    
    for metric in qc_metrics:
        data1 = group1[metric].dropna()
        data2 = group2[metric].dropna()
        
        mean1 = data1.mean()
        std1 = data1.std()
        mean2 = data2.mean()
        std2 = data2.std()
        
        # Mann-Whitney U test (non-parametric)
        statistic, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
        
        significant = '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else 'ns'
        
        print(f"{metric:<30} {mean1:>7.3f}±{std1:<5.3f} {mean2:>7.3f}±{std2:<5.3f} {p_value:>12.2e} {significant:>12}")
        
        stats_results.append({
            'metric': metric,
            'group1_mean': mean1,
            'group1_std': std1,
            'group2_mean': mean2,
            'group2_std': std2,
            'p_value': p_value,
            'significant': significant
        })
    
    print(f"\nSignificance levels: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant\n")
    
    # ======================
    # 4. Box plots for each metric
    # ======================
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes = axes.flatten()
    
    for i, metric in enumerate(qc_metrics):
        ax = axes[i]
        
        data_to_plot = [
            group1[metric].dropna().values,
            group2[metric].dropna().values
        ]
        
        bp = ax.boxplot(data_to_plot, labels=['Group 1\n(Left)', 'Group 2\n(Right)'],
                       patch_artist=True, widths=0.6)
        
        # Color boxes
        colors = ['steelblue', 'coral']
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        # Add p-value annotation
        p_val = stats_results[i]['p_value']
        sig_text = stats_results[i]['significant']
        
        y_max = max([d.max() for d in data_to_plot if len(d) > 0])
        y_min = min([d.min() for d in data_to_plot if len(d) > 0])
        y_range = y_max - y_min
        
        ax.text(1.5, y_max + y_range * 0.05, f'p={p_val:.2e}\n{sig_text}',
               ha='center', fontsize=9, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.set_ylabel(metric, fontsize=10)
        ax.set_title(f'{metric}', fontsize=11, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
    
    # Remove extra subplot
    fig.delaxes(axes[5])
    
    plt.suptitle(f'QC Metrics Comparison: Group 1 vs Group 2 - {lab_name}', 
                fontsize=14, fontweight='bold', y=1.00)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pca_groups_comparison.png'), dpi=150)
    print(f"✓ Saved: pca_groups_comparison.png\n")
    plt.close()
    
    # ======================
    # 5. Save detailed report
    # ======================
    report_file = os.path.join(output_dir, 'pca_clustering_report.txt')
    
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write(f"PCA CLUSTERING ANALYSIS - {lab_name}\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Total samples: {len(df_pca)}\n")
        f.write(f"PC1 variance: {pca.explained_variance_ratio_[0]*100:.1f}%\n")
        f.write(f"PC2 variance: {pca.explained_variance_ratio_[1]*100:.1f}%\n\n")
        
        f.write("CLUSTER SIZES\n")
        f.write("-"*80 + "\n")
        for label, count in cluster_counts.items():
            pct = count / len(df_pca) * 100
            f.write(f"{label}: {count} samples ({pct:.1f}%)\n")
        f.write("\n")
        
        f.write("STATISTICAL COMPARISON\n")
        f.write("-"*80 + "\n")
        f.write(f"{'Metric':<30} {'Group1':>12} {'Group2':>12} {'p-value':>12} {'Sig':>6}\n")
        f.write("-"*80 + "\n")
        
        for result in stats_results:
            f.write(f"{result['metric']:<30} "
                   f"{result['group1_mean']:>12.3f} "
                   f"{result['group2_mean']:>12.3f} "
                   f"{result['p_value']:>12.2e} "
                   f"{result['significant']:>6}\n")
        
        f.write("\n")
        f.write("BATCH COMPOSITION\n")
        f.write("-"*80 + "\n\n")
        
        for cluster_id, label in enumerate(['Group 1 (Left)', 'Group 2 (Right)']):
            cluster_data = df_pca[df_pca['cluster'] == cluster_id]
            f.write(f"{label} - Top 15 Batches:\n")
            
            batch_counts = cluster_data['batch_id'].value_counts().head(15)
            for batch, count in batch_counts.items():
                pct = count / len(cluster_data) * 100
                f.write(f"  {batch:30s}: {count:4d} ({pct:5.1f}%)\n")
            f.write("\n")
    
    print(f"✓ Saved: pca_clustering_report.txt\n")
    
    # Save clustered data
    output_clustered = os.path.join(output_dir, 'sample_list_with_clusters.tsv')
    df_pca.to_csv(output_clustered, sep='\t', index=False)
    print(f"✓ Saved: sample_list_with_clusters.tsv\n")
    
    print(f"{'='*80}\n")
    print("✓ Analysis complete!")
    print(f"\nGenerated files:")
    print(f"  - pca_clustering_groups.png: PCA plot with cluster labels")
    print(f"  - pca_groups_comparison.png: Box plots comparing QC metrics")
    print(f"  - pca_clustering_report.txt: Detailed statistical report")
    print(f"  - sample_list_with_clusters.tsv: Sample data with cluster assignments")
    print()

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python3 analyze_pca_clustering.py <sample_list.tsv> <output_dir> <lab_name>")
        sys.exit(1)
    
    analyze_pca_clustering(sys.argv[1], sys.argv[2], sys.argv[3])
