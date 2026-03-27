#!/usr/bin/env python3
"""
Batch Effect 분석 스크립트
시기별(month별) UCL 샘플들의 QC metrics와 분포를 비교
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import sys

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (15, 10)

def load_and_prepare_data(sample_list_file):
    """샘플 리스트 로드 및 전처리"""
    print(f"Loading: {sample_list_file}")
    df = pd.read_csv(sample_list_file, sep='\t')
    
    # Data type conversion
    numeric_cols = ['number_of_reads', 'number_of_mapped_reads', 'mapping_rate(%)', 
                    'duplication_rate(%)', 'mean_mapping_quality', 'mean_coverageData(X)',
                    'GC_content(%)', 'SeqFF', 'Fragment_FF', 'YFF_2', 'M-SeqFF']
    
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    print(f"Total samples: {len(df)}")
    print(f"Months: {sorted(df['month'].unique())}")
    print(f"Samples per month:")
    print(df['month'].value_counts().sort_index())
    
    return df

def analyze_qc_metrics(df, output_dir):
    """QC metrics 비교 (월별)"""
    print("\n" + "="*60)
    print("1. QC Metrics Analysis (by Month)")
    print("="*60)
    
    metrics = ['mapping_rate(%)', 'duplication_rate(%)', 'mean_mapping_quality', 
               'GC_content(%)', 'mean_coverageData(X)']
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    for idx, metric in enumerate(metrics):
        if metric not in df.columns:
            continue
        
        ax = axes[idx]
        
        # Box plot
        df_clean = df.dropna(subset=[metric])
        months = sorted(df_clean['month'].unique())
        data_by_month = [df_clean[df_clean['month'] == m][metric].values for m in months]
        
        bp = ax.boxplot(data_by_month, labels=months, patch_artist=True)
        
        # Color by month
        colors = plt.cm.Set3(np.linspace(0, 1, len(months)))
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
        
        ax.set_xlabel('Month', fontsize=12)
        ax.set_ylabel(metric, fontsize=12)
        ax.set_title(f'{metric} by Month', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Statistics per month
        print(f"\n{metric}:")
        for month in months:
            month_data = df_clean[df_clean['month'] == month][metric]
            print(f"  {month}: mean={month_data.mean():.2f}, std={month_data.std():.2f}, "
                  f"median={month_data.median():.2f}, n={len(month_data)}")
        
        # Kruskal-Wallis test (non-parametric)
        if len(months) > 1:
            h_stat, p_value = stats.kruskal(*data_by_month)
            print(f"  Kruskal-Wallis test: H={h_stat:.2f}, p={p_value:.4f}")
            if p_value < 0.05:
                print(f"  *** Significant difference detected (p < 0.05) ***")
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '01_qc_metrics_by_month.png'), dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/01_qc_metrics_by_month.png")
    plt.close()

def analyze_fetal_fraction(df, output_dir):
    """Fetal Fraction 비교 (월별)"""
    print("\n" + "="*60)
    print("2. Fetal Fraction Analysis (by Month)")
    print("="*60)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # SeqFF
    df_seqff = df.dropna(subset=['SeqFF'])
    months = sorted(df_seqff['month'].unique())
    data_by_month = [df_seqff[df_seqff['month'] == m]['SeqFF'].values for m in months]
    
    ax = axes[0]
    bp = ax.boxplot(data_by_month, labels=months, patch_artist=True)
    colors = plt.cm.Set3(np.linspace(0, 1, len(months)))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    
    ax.set_xlabel('Month', fontsize=12)
    ax.set_ylabel('SeqFF (%)', fontsize=12)
    ax.set_title('SeqFF by Month', fontsize=14, fontweight='bold')
    ax.axhline(y=4.0, color='red', linestyle='--', linewidth=1, label='Min threshold (4%)')
    ax.axhline(y=30.0, color='red', linestyle='--', linewidth=1, label='Max threshold (30%)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    print("\nSeqFF:")
    for month in months:
        month_data = df_seqff[df_seqff['month'] == month]['SeqFF']
        below_4 = (month_data < 4.0).sum()
        above_30 = (month_data > 30.0).sum()
        print(f"  {month}: mean={month_data.mean():.2f}, std={month_data.std():.2f}, "
              f"median={month_data.median():.2f}, <4%={below_4}, >30%={above_30}, n={len(month_data)}")
    
    if len(months) > 1:
        h_stat, p_value = stats.kruskal(*data_by_month)
        print(f"  Kruskal-Wallis test: H={h_stat:.2f}, p={p_value:.4f}")
        if p_value < 0.05:
            print(f"  *** Significant difference detected (p < 0.05) ***")
    
    # Gender distribution
    ax = axes[1]
    gender_counts = df.groupby(['month', 'fetal_gender(gd_2)']).size().unstack(fill_value=0)
    gender_counts.plot(kind='bar', ax=ax, color=['lightblue', 'lightpink'])
    ax.set_xlabel('Month', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Gender Distribution by Month', fontsize=14, fontweight='bold')
    ax.legend(title='Gender')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '02_fetal_fraction_by_month.png'), dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/02_fetal_fraction_by_month.png")
    plt.close()

def analyze_result_distribution(df, output_dir):
    """Result 분포 비교 (월별)"""
    print("\n" + "="*60)
    print("3. Result Distribution Analysis (by Month)")
    print("="*60)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Result distribution
    ax = axes[0]
    result_counts = df.groupby(['month', 'Result']).size().unstack(fill_value=0)
    result_counts.plot(kind='bar', stacked=True, ax=ax, colormap='Set3')
    ax.set_xlabel('Month', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Result Distribution by Month', fontsize=14, fontweight='bold')
    ax.legend(title='Result', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    print("\nResult Distribution:")
    for month in sorted(df['month'].unique()):
        month_df = df[df['month'] == month]
        print(f"  {month}:")
        print(f"    Total: {len(month_df)}")
        for result, count in month_df['Result'].value_counts().items():
            pct = count / len(month_df) * 100
            print(f"    {result}: {count} ({pct:.1f}%)")
    
    # MDResult distribution
    ax = axes[1]
    mdresult_counts = df.groupby(['month', 'MDResult']).size().unstack(fill_value=0)
    mdresult_counts.plot(kind='bar', stacked=True, ax=ax, colormap='Pastel1')
    ax.set_xlabel('Month', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('MDResult Distribution by Month', fontsize=14, fontweight='bold')
    ax.legend(title='MDResult', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '03_result_distribution_by_month.png'), dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/03_result_distribution_by_month.png")
    plt.close()

def analyze_correlation_matrix(df, output_dir):
    """월별 상관관계 행렬 비교"""
    print("\n" + "="*60)
    print("4. Correlation Matrix Analysis (by Month)")
    print("="*60)
    
    metrics = ['mapping_rate(%)', 'duplication_rate(%)', 'GC_content(%)', 
               'SeqFF', 'mean_coverageData(X)', 'mean_mapping_quality']
    
    months = sorted(df['month'].unique())
    n_months = len(months)
    
    fig, axes = plt.subplots(2, (n_months + 1) // 2, figsize=(6 * ((n_months + 1) // 2), 12))
    axes = axes.flatten()
    
    for idx, month in enumerate(months):
        month_df = df[df['month'] == month][metrics].dropna()
        
        if len(month_df) < 10:
            continue
        
        corr = month_df.corr()
        
        ax = axes[idx]
        sns.heatmap(corr, annot=True, fmt='.2f', cmap='coolwarm', 
                    center=0, vmin=-1, vmax=1, ax=ax,
                    cbar_kws={'label': 'Correlation'})
        ax.set_title(f'{month} (n={len(month_df)})', fontsize=14, fontweight='bold')
    
    # Hide unused subplots
    for idx in range(len(months), len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '04_correlation_by_month.png'), dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/04_correlation_by_month.png")
    plt.close()

def analyze_pca(df, output_dir):
    """PCA 분석으로 batch separation 확인"""
    print("\n" + "="*60)
    print("5. PCA Analysis (Batch Separation)")
    print("="*60)
    
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    
    metrics = ['mapping_rate(%)', 'duplication_rate(%)', 'GC_content(%)', 
               'SeqFF', 'mean_coverageData(X)', 'mean_mapping_quality']
    
    # Prepare data
    df_pca = df[['month'] + metrics].dropna()
    
    if len(df_pca) < 50:
        print("Not enough samples for PCA analysis")
        return
    
    X = df_pca[metrics].values
    months = df_pca['month'].values
    
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    
    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    unique_months = sorted(df_pca['month'].unique())
    colors = plt.cm.Set3(np.linspace(0, 1, len(unique_months)))
    
    for month, color in zip(unique_months, colors):
        mask = months == month
        ax.scatter(X_pca[mask, 0], X_pca[mask, 1], 
                  c=[color], label=month, alpha=0.6, s=50)
    
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=12)
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=12)
    ax.set_title('PCA: Batch Separation Analysis', fontsize=14, fontweight='bold')
    ax.legend(title='Month', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '05_pca_batch_separation.png'), dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/05_pca_batch_separation.png")
    
    print("\nPCA Results:")
    print(f"  PC1 explains {pca.explained_variance_ratio_[0]*100:.1f}% of variance")
    print(f"  PC2 explains {pca.explained_variance_ratio_[1]*100:.1f}% of variance")
    print(f"  Total: {sum(pca.explained_variance_ratio_)*100:.1f}%")
    
    print("\nPC1 loadings (most important features):")
    loadings = pd.DataFrame({
        'Feature': metrics,
        'PC1': pca.components_[0],
        'PC2': pca.components_[1]
    })
    loadings['PC1_abs'] = np.abs(loadings['PC1'])
    loadings = loadings.sort_values('PC1_abs', ascending=False)
    for _, row in loadings.iterrows():
        print(f"  {row['Feature']}: {row['PC1']:.3f}")
    
    plt.close()

def generate_summary_report(df, output_dir):
    """요약 리포트 생성"""
    print("\n" + "="*60)
    print("6. Generating Summary Report")
    print("="*60)
    
    report_file = os.path.join(output_dir, 'batch_effect_summary.txt')
    
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("UCL BATCH EFFECT ANALYSIS REPORT\n")
        f.write("="*80 + "\n\n")
        
        # Overall statistics
        f.write("1. OVERALL STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total samples: {len(df)}\n")
        f.write(f"Months: {', '.join(sorted(df['month'].unique()))}\n\n")
        
        # Samples per month
        f.write("Samples per month:\n")
        for month, count in df['month'].value_counts().sort_index().items():
            pct = count / len(df) * 100
            f.write(f"  {month}: {count} ({pct:.1f}%)\n")
        f.write("\n")
        
        # Gender distribution
        f.write("Gender distribution:\n")
        for month in sorted(df['month'].unique()):
            month_df = df[df['month'] == month]
            xy = (month_df['fetal_gender(gd_2)'] == 'XY').sum()
            xx = (month_df['fetal_gender(gd_2)'] == 'XX').sum()
            f.write(f"  {month}: XY={xy}, XX={xx}\n")
        f.write("\n")
        
        # QC metrics summary
        f.write("2. QC METRICS SUMMARY (by Month)\n")
        f.write("-" * 80 + "\n")
        
        metrics = ['mapping_rate(%)', 'duplication_rate(%)', 'GC_content(%)', 'SeqFF']
        
        for metric in metrics:
            if metric not in df.columns:
                continue
            
            f.write(f"\n{metric}:\n")
            for month in sorted(df['month'].unique()):
                month_data = df[df['month'] == month][metric].dropna()
                if len(month_data) > 0:
                    f.write(f"  {month}: mean={month_data.mean():.2f} ± {month_data.std():.2f}, "
                           f"median={month_data.median():.2f}, n={len(month_data)}\n")
        
        # Result distribution
        f.write("\n3. RESULT DISTRIBUTION (by Month)\n")
        f.write("-" * 80 + "\n")
        for month in sorted(df['month'].unique()):
            month_df = df[df['month'] == month]
            f.write(f"\n{month}:\n")
            for result, count in month_df['Result'].value_counts().items():
                pct = count / len(month_df) * 100
                f.write(f"  {result}: {count} ({pct:.1f}%)\n")
        
        f.write("\n" + "="*80 + "\n")
    
    print(f"✓ Saved: {report_file}")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Batch Effect Analysis for UCL samples')
    parser.add_argument('--sample-list', required=True, 
                       help='Sample list TSV file (e.g., reference_sample_list_UCL.tsv)')
    parser.add_argument('--output-dir', default='batch_effect_analysis',
                       help='Output directory (default: batch_effect_analysis)')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("="*80)
    print("UCL BATCH EFFECT ANALYSIS")
    print("="*80)
    print(f"Sample list: {args.sample_list}")
    print(f"Output dir: {args.output_dir}")
    print()
    
    # Load data
    df = load_and_prepare_data(args.sample_list)
    
    # Run analyses
    analyze_qc_metrics(df, args.output_dir)
    analyze_fetal_fraction(df, args.output_dir)
    analyze_result_distribution(df, args.output_dir)
    analyze_correlation_matrix(df, args.output_dir)
    analyze_pca(df, args.output_dir)
    generate_summary_report(df, args.output_dir)
    
    print("\n" + "="*80)
    print("✅ ANALYSIS COMPLETE!")
    print("="*80)
    print(f"\nResults saved in: {args.output_dir}/")
    print("Generated files:")
    print("  - 01_qc_metrics_by_month.png")
    print("  - 02_fetal_fraction_by_month.png")
    print("  - 03_result_distribution_by_month.png")
    print("  - 04_correlation_by_month.png")
    print("  - 05_pca_batch_separation.png")
    print("  - batch_effect_summary.txt")

if __name__ == '__main__':
    main()
