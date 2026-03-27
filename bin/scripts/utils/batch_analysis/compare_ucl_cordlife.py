#!/usr/bin/env python3
"""
UCL vs Cordlife 비교 스크립트
각각 분석된 결과를 바탕으로 직접 비교
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import sys

sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (15, 10)

def load_batch_data(ucl_file, cordlife_file):
    """UCL과 Cordlife batch 데이터 로드"""
    print("Loading UCL data...")
    df_ucl = pd.read_csv(ucl_file, sep='\t')
    print(f"  UCL: {len(df_ucl)} samples")
    
    print("Loading Cordlife data...")
    df_cordlife = pd.read_csv(cordlife_file, sep='\t')
    print(f"  Cordlife: {len(df_cordlife)} samples")
    
    # Combine
    df = pd.concat([df_ucl, df_cordlife], ignore_index=True)
    print(f"\n✓ Total: {len(df)} samples")
    
    return df

def compare_qc_metrics(df, output_dir):
    """QC metrics 비교"""
    print("\n" + "="*60)
    print("1. QC Metrics Comparison")
    print("="*60)
    
    metrics = ['mapping_rate(%)', 'duplication_rate(%)', 'mean_mapping_quality',
               'GC_content(%)', 'SeqFF', 'mean_coverageData(X)']
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    comparison_results = []
    
    for idx, metric in enumerate(metrics):
        if metric not in df.columns:
            continue
        
        ax = axes[idx]
        
        # Violin plot
        df_clean = df.dropna(subset=[metric])
        
        ucl_data = df_clean[df_clean['lab'] == 'UCL'][metric]
        cordlife_data = df_clean[df_clean['lab'] == 'Cordlife'][metric]
        
        if len(ucl_data) == 0 or len(cordlife_data) == 0:
            continue
        
        parts = ax.violinplot(
            [ucl_data.values, cordlife_data.values],
            positions=[1, 2],
            showmeans=True,
            showmedians=True,
            widths=0.7
        )
        
        # Add box plot overlay
        bp = ax.boxplot([ucl_data.values, cordlife_data.values],
                        positions=[1, 2],
                        widths=0.3,
                        patch_artist=False,
                        showfliers=False)
        
        ax.set_xticks([1, 2])
        ax.set_xticklabels(['UCL', 'Cordlife'])
        ax.set_ylabel(metric, fontsize=12)
        ax.set_title(f'{metric}: UCL vs Cordlife', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Statistics
        ucl_mean = ucl_data.mean()
        ucl_std = ucl_data.std()
        ucl_median = ucl_data.median()
        cordlife_mean = cordlife_data.mean()
        cordlife_std = cordlife_data.std()
        cordlife_median = cordlife_data.median()
        
        print(f"\n{metric}:")
        print(f"  UCL: mean={ucl_mean:.2f} ± {ucl_std:.2f}, median={ucl_median:.2f}, n={len(ucl_data)}")
        print(f"  Cordlife: mean={cordlife_mean:.2f} ± {cordlife_std:.2f}, median={cordlife_median:.2f}, n={len(cordlife_data)}")
        
        # Mann-Whitney U test
        u_stat, p_value = stats.mannwhitneyu(ucl_data, cordlife_data, alternative='two-sided')
        print(f"  Mann-Whitney U test: U={u_stat:.2f}, p={p_value:.4f}")
        
        significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
        print(f"  Significance: {significance}")
        
        # Add significance to plot
        y_max = max(ucl_data.max(), cordlife_data.max())
        y_range = y_max - min(ucl_data.min(), cordlife_data.min())
        ax.text(1.5, y_max + y_range * 0.05, significance, 
                ha='center', fontsize=16, fontweight='bold')
        
        comparison_results.append({
            'metric': metric,
            'ucl_mean': ucl_mean,
            'ucl_std': ucl_std,
            'ucl_median': ucl_median,
            'cordlife_mean': cordlife_mean,
            'cordlife_std': cordlife_std,
            'cordlife_median': cordlife_median,
            'u_stat': u_stat,
            'p_value': p_value,
            'significance': significance
        })
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ucl_vs_cordlife_comparison.png'), dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/ucl_vs_cordlife_comparison.png")
    plt.close()
    
    # Save comparison table
    comparison_df = pd.DataFrame(comparison_results)
    comparison_file = os.path.join(output_dir, 'ucl_vs_cordlife_stats.tsv')
    comparison_df.to_csv(comparison_file, sep='\t', index=False)
    print(f"✓ Saved: {comparison_file}")

def compare_result_distribution(df, output_dir):
    """Result 분포 비교"""
    print("\n" + "="*60)
    print("2. Result Distribution Comparison")
    print("="*60)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Result distribution
    ax = axes[0]
    result_counts = df.groupby(['lab', 'Result']).size().unstack(fill_value=0)
    result_pct = result_counts.div(result_counts.sum(axis=1), axis=0) * 100
    
    result_pct.plot(kind='bar', ax=ax, colormap='Set3')
    ax.set_xlabel('Lab', fontsize=12)
    ax.set_ylabel('Percentage (%)', fontsize=12)
    ax.set_title('Result Distribution', fontsize=14, fontweight='bold')
    ax.legend(title='Result', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    ax.grid(True, alpha=0.3)
    
    print("\nResult Distribution:")
    for lab in ['UCL', 'Cordlife']:
        lab_df = df[df['lab'] == lab]
        print(f"  {lab}:")
        for result, count in lab_df['Result'].value_counts().items():
            pct = count / len(lab_df) * 100
            print(f"    {result}: {count} ({pct:.1f}%)")
    
    # Gender distribution
    ax = axes[1]
    gender_counts = df.groupby(['lab', 'fetal_gender(gd_2)']).size().unstack(fill_value=0)
    gender_pct = gender_counts.div(gender_counts.sum(axis=1), axis=0) * 100
    
    gender_pct.plot(kind='bar', ax=ax, color=['lightblue', 'lightpink'])
    ax.set_xlabel('Lab', fontsize=12)
    ax.set_ylabel('Percentage (%)', fontsize=12)
    ax.set_title('Gender Distribution', fontsize=14, fontweight='bold')
    ax.legend(title='Gender')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ucl_vs_cordlife_distributions.png'), dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/ucl_vs_cordlife_distributions.png")
    plt.close()

def compare_batch_characteristics(df, output_dir):
    """Batch 특성 비교"""
    print("\n" + "="*60)
    print("3. Batch Characteristics Comparison")
    print("="*60)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Average batch size
    ax = axes[0]
    batch_sizes = df.groupby(['lab', 'batch_id']).size()
    
    ucl_batch_sizes = batch_sizes[batch_sizes.index.get_level_values(0) == 'UCL']
    cordlife_batch_sizes = batch_sizes[batch_sizes.index.get_level_values(0) == 'Cordlife']
    
    if len(ucl_batch_sizes) > 0 and len(cordlife_batch_sizes) > 0:
        bp = ax.boxplot([ucl_batch_sizes.values, cordlife_batch_sizes.values],
                        labels=['UCL', 'Cordlife'],
                        patch_artist=True)
        bp['boxes'][0].set_facecolor('lightblue')
        bp['boxes'][1].set_facecolor('lightcoral')
    
    ax.set_ylabel('Samples per Batch', fontsize=12)
    ax.set_title('Batch Size Distribution', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    print("\nBatch Size:")
    print(f"  UCL: mean={ucl_batch_sizes.mean():.1f} ± {ucl_batch_sizes.std():.1f}, "
          f"median={ucl_batch_sizes.median():.0f}, n_batches={len(ucl_batch_sizes)}")
    print(f"  Cordlife: mean={cordlife_batch_sizes.mean():.1f} ± {cordlife_batch_sizes.std():.1f}, "
          f"median={cordlife_batch_sizes.median():.0f}, n_batches={len(cordlife_batch_sizes)}")
    
    # Number of batches per month
    ax = axes[1]
    batches_per_month = df.groupby(['lab', 'month'])['batch_id'].nunique().unstack(fill_value=0)
    batches_per_month.T.plot(kind='bar', ax=ax, color=['lightblue', 'lightcoral'])
    ax.set_xlabel('Month', fontsize=12)
    ax.set_ylabel('Number of Batches', fontsize=12)
    ax.set_title('Batches per Month', fontsize=14, fontweight='bold')
    ax.legend(title='Lab')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ucl_vs_cordlife_batch_chars.png'), dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/ucl_vs_cordlife_batch_chars.png")
    plt.close()

def generate_comparison_report(df, output_dir):
    """비교 리포트 생성"""
    print("\n" + "="*60)
    print("4. Generating Comparison Report")
    print("="*60)
    
    report_file = os.path.join(output_dir, 'ucl_vs_cordlife_report.txt')
    
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("UCL vs CORDLIFE COMPARISON REPORT\n")
        f.write("="*80 + "\n\n")
        
        # Overall statistics
        f.write("1. OVERALL STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total samples: {len(df)}\n")
        f.write(f"  UCL: {(df['lab'] == 'UCL').sum()}\n")
        f.write(f"  Cordlife: {(df['lab'] == 'Cordlife').sum()}\n\n")
        
        # Batch statistics
        f.write("2. BATCH STATISTICS\n")
        f.write("-" * 80 + "\n")
        
        for lab in ['UCL', 'Cordlife']:
            lab_df = df[df['lab'] == lab]
            n_batches = lab_df['batch_id'].nunique()
            batch_sizes = lab_df.groupby('batch_id').size()
            
            f.write(f"\n{lab}:\n")
            f.write(f"  Number of batches: {n_batches}\n")
            f.write(f"  Avg batch size: {batch_sizes.mean():.1f} ± {batch_sizes.std():.1f}\n")
            f.write(f"  Median batch size: {batch_sizes.median():.0f}\n")
            f.write(f"  Min/Max batch size: {batch_sizes.min()}/{batch_sizes.max()}\n")
        
        # QC metrics comparison
        f.write("\n3. QC METRICS COMPARISON\n")
        f.write("-" * 80 + "\n")
        
        metrics = ['mapping_rate(%)', 'duplication_rate(%)', 'GC_content(%)', 'SeqFF']
        
        for metric in metrics:
            if metric not in df.columns:
                continue
            
            ucl_data = df[df['lab'] == 'UCL'][metric].dropna()
            cordlife_data = df[df['lab'] == 'Cordlife'][metric].dropna()
            
            if len(ucl_data) == 0 or len(cordlife_data) == 0:
                continue
            
            f.write(f"\n{metric}:\n")
            f.write(f"  UCL: mean={ucl_data.mean():.2f} ± {ucl_data.std():.2f}, "
                   f"median={ucl_data.median():.2f}\n")
            f.write(f"  Cordlife: mean={cordlife_data.mean():.2f} ± {cordlife_data.std():.2f}, "
                   f"median={cordlife_data.median():.2f}\n")
            
            u_stat, p_value = stats.mannwhitneyu(ucl_data, cordlife_data, alternative='two-sided')
            significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else "ns"
            f.write(f"  Mann-Whitney U test: U={u_stat:.2f}, p={p_value:.4f} ({significance})\n")
        
        # Result distribution
        f.write("\n4. RESULT DISTRIBUTION\n")
        f.write("-" * 80 + "\n")
        
        for lab in ['UCL', 'Cordlife']:
            lab_df = df[df['lab'] == lab]
            f.write(f"\n{lab}:\n")
            for result, count in lab_df['Result'].value_counts().items():
                pct = count / len(lab_df) * 100
                f.write(f"  {result}: {count} ({pct:.1f}%)\n")
        
        f.write("\n" + "="*80 + "\n")
    
    print(f"✓ Saved: {report_file}")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Compare UCL and Cordlife batch analysis results')
    parser.add_argument('--ucl', required=True,
                       help='UCL batch analysis result (sample_list_with_batch.tsv)')
    parser.add_argument('--cordlife', required=True,
                       help='Cordlife batch analysis result (sample_list_with_batch.tsv)')
    parser.add_argument('--output-dir', default='ucl_cordlife_comparison',
                       help='Output directory (default: ucl_cordlife_comparison)')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("="*80)
    print("UCL vs CORDLIFE COMPARISON")
    print("="*80)
    print(f"UCL data: {args.ucl}")
    print(f"Cordlife data: {args.cordlife}")
    print(f"Output dir: {args.output_dir}")
    print()
    
    # Load data
    df = load_batch_data(args.ucl, args.cordlife)
    
    # Run comparisons
    compare_qc_metrics(df, args.output_dir)
    compare_result_distribution(df, args.output_dir)
    compare_batch_characteristics(df, args.output_dir)
    generate_comparison_report(df, args.output_dir)
    
    print("\n" + "="*80)
    print("✅ COMPARISON COMPLETE!")
    print("="*80)
    print(f"\nResults saved in: {args.output_dir}/")
    print("Generated files:")
    print("  - ucl_vs_cordlife_comparison.png")
    print("  - ucl_vs_cordlife_stats.tsv")
    print("  - ucl_vs_cordlife_distributions.png")
    print("  - ucl_vs_cordlife_batch_chars.png")
    print("  - ucl_vs_cordlife_report.txt")

if __name__ == '__main__':
    main()
