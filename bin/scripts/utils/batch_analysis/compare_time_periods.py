#!/usr/bin/env python3
"""
시기별 샘플 품질 비교 - 2509 이전 vs 이후
"""

import pandas as pd
import numpy as np
from scipy import stats
import sys
import os

def compare_time_periods(output_dir, lab_name):
    """
    2509 이전 vs 이후 샘플 품질 비교
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
    
    print(f"\n{'='*80}")
    print(f"TIME PERIOD COMPARISON - {lab_name}")
    print(f"{'='*80}\n")
    
    # Split by month
    df['month_int'] = df['month'].astype(int)
    
    # Define periods
    period1 = df[df['month_int'] < 2509]  # 2507-2508
    period2 = df[(df['month_int'] >= 2509) & (df['month_int'] < 2510)]  # 2509
    period3 = df[df['month_int'] >= 2510]  # 2510 이후
    
    print(f"Period 1 (2507-2508): {len(period1)} samples")
    print(f"Period 2 (2509): {len(period2)} samples")
    print(f"Period 3 (2510+): {len(period3)} samples\n")
    
    # QC metrics
    qc_metrics = ['duplication_rate(%)', 'mean_coverageData(X)', 
                  'GC_content(%)', 'mapping_rate(%)']
    chr_metrics = ['chr13_prop', 'chr18_prop', 'chr21_prop']
    
    # ========================================
    # 1. QC Metrics Comparison
    # ========================================
    print(f"{'='*80}")
    print(f"QC METRICS BY PERIOD")
    print(f"{'='*80}\n")
    
    for metric in qc_metrics:
        print(f"\n{metric}:")
        print(f"{'-'*70}")
        print(f"{'Period':<15} {'Mean':>12} {'Median':>12} {'Std':>10} {'Min':>10} {'Max':>10}")
        print(f"{'-'*70}")
        
        for period_name, period_data in [
            ('2507-2508', period1),
            ('2509', period2),
            ('2510+', period3)
        ]:
            data = period_data[metric].dropna()
            if len(data) > 0:
                if 'coverage' in metric.lower():
                    print(f"{period_name:<15} {data.mean():>12.3f} {data.median():>12.3f} "
                          f"{data.std():>10.3f} {data.min():>10.3f} {data.max():>10.3f}")
                else:
                    print(f"{period_name:<15} {data.mean():>12.2f} {data.median():>12.2f} "
                          f"{data.std():>10.2f} {data.min():>10.2f} {data.max():>10.2f}")
    
    # ========================================
    # 2. Chromosome Proportion Comparison
    # ========================================
    print(f"\n{'='*80}")
    print(f"CHROMOSOME PROPORTIONS BY PERIOD")
    print(f"{'='*80}\n")
    
    for chr_metric in chr_metrics:
        chr_name = chr_metric.replace('_prop', '').upper()
        print(f"\n{chr_name}:")
        print(f"{'-'*70}")
        print(f"{'Period':<15} {'Mean':>15} {'Std':>12} {'n':>8}")
        print(f"{'-'*70}")
        
        for period_name, period_data in [
            ('2507-2508', period1),
            ('2509', period2),
            ('2510+', period3)
        ]:
            data = period_data[chr_metric].dropna()
            if len(data) > 0:
                print(f"{period_name:<15} {data.mean():>15.6f} {data.std():>12.6f} {len(data):>8}")
    
    # ========================================
    # 3. Cluster Distribution by Period
    # ========================================
    print(f"\n{'='*80}")
    print(f"CLUSTER DISTRIBUTION BY PERIOD")
    print(f"{'='*80}\n")
    
    print(f"{'Period':<15} {'Total':>8} {'Group 1':>10} {'Group 2':>10} {'G1%':>8} {'G2%':>8}")
    print(f"{'-'*65}")
    
    for period_name, period_data in [
        ('2507-2508', period1),
        ('2509', period2),
        ('2510+', period3)
    ]:
        total = len(period_data)
        g1 = len(period_data[period_data['cluster'] == 0])
        g2 = len(period_data[period_data['cluster'] == 1])
        g1_pct = (g1 / total * 100) if total > 0 else 0
        g2_pct = (g2 / total * 100) if total > 0 else 0
        
        print(f"{period_name:<15} {total:>8} {g1:>10} {g2:>10} {g1_pct:>7.1f}% {g2_pct:>7.1f}%")
    
    # ========================================
    # 4. Statistical Tests: 2509+ vs 2507-2508
    # ========================================
    print(f"\n{'='*80}")
    print(f"STATISTICAL COMPARISON: 2510+ vs 2507-2508")
    print(f"{'='*80}\n")
    
    print(f"{'Metric':<30} {'2507-2508':>15} {'2510+':>15} {'p-value':>12} {'Sig':>6}")
    print(f"{'-'*80}")
    
    for metric in qc_metrics + chr_metrics:
        data_early = period1[metric].dropna()
        data_late = period3[metric].dropna()
        
        if len(data_early) > 0 and len(data_late) > 0:
            mean_early = data_early.mean()
            mean_late = data_late.mean()
            
            stat, p = stats.mannwhitneyu(data_early, data_late, alternative='two-sided')
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
            
            if 'chr' in metric:
                print(f"{metric:<30} {mean_early:>15.6f} {mean_late:>15.6f} {p:>12.2e} {sig:>6}")
            elif 'coverage' in metric.lower():
                print(f"{metric:<30} {mean_early:>15.3f} {mean_late:>15.3f} {p:>12.2e} {sig:>6}")
            else:
                print(f"{metric:<30} {mean_early:>15.2f} {mean_late:>15.2f} {p:>12.2e} {sig:>6}")
    
    # ========================================
    # 5. 2510+ 샘플에서의 두 그룹 비교
    # ========================================
    print(f"\n{'='*80}")
    print(f"CLUSTER COMPARISON WITHIN 2510+ SAMPLES")
    print(f"{'='*80}\n")
    
    period3_g1 = period3[period3['cluster'] == 0]
    period3_g2 = period3[period3['cluster'] == 1]
    
    print(f"Group 1 (Low cov/dup): {len(period3_g1)} samples")
    print(f"Group 2 (High cov/dup): {len(period3_g2)} samples\n")
    
    print(f"{'Metric':<30} {'Group 1':>15} {'Group 2':>15} {'p-value':>12} {'Sig':>6}")
    print(f"{'-'*80}")
    
    for metric in qc_metrics + chr_metrics:
        data_g1 = period3_g1[metric].dropna()
        data_g2 = period3_g2[metric].dropna()
        
        if len(data_g1) > 0 and len(data_g2) > 0:
            mean_g1 = data_g1.mean()
            mean_g2 = data_g2.mean()
            
            stat, p = stats.mannwhitneyu(data_g1, data_g2, alternative='two-sided')
            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
            
            if 'chr' in metric:
                print(f"{metric:<30} {mean_g1:>15.6f} {mean_g2:>15.6f} {p:>12.2e} {sig:>6}")
            elif 'coverage' in metric.lower():
                print(f"{metric:<30} {mean_g1:>15.3f} {mean_g2:>15.3f} {p:>12.2e} {sig:>6}")
            else:
                print(f"{metric:<30} {mean_g1:>15.2f} {mean_g2:>15.2f} {p:>12.2e} {sig:>6}")
    
    # ========================================
    # 6. Recommendations
    # ========================================
    print(f"\n{'='*80}")
    print(f"RECOMMENDATIONS")
    print(f"{'='*80}\n")
    
    # Check duplication in 2510+
    dup_2510 = period3['duplication_rate(%)'].dropna()
    high_dup_count = len(dup_2510[dup_2510 > 15])
    
    print(f"2510+ Samples Analysis:")
    print(f"  - Total samples: {len(period3)}")
    print(f"  - Mean duplication: {dup_2510.mean():.2f}%")
    print(f"  - Samples with dup > 15%: {high_dup_count} ({high_dup_count/len(period3)*100:.1f}%)")
    
    # Check chromosome proportion variability
    chr_diff_2510 = {}
    for chr_metric in chr_metrics:
        data_g1 = period3_g1[chr_metric].dropna()
        data_g2 = period3_g2[chr_metric].dropna()
        if len(data_g1) > 0 and len(data_g2) > 0:
            diff_pct = abs(data_g1.mean() - data_g2.mean()) / data_g2.mean() * 100
            chr_diff_2510[chr_metric] = diff_pct
    
    print(f"\n2510+ Group Differences:")
    for chr_metric, diff in chr_diff_2510.items():
        chr_name = chr_metric.replace('_prop', '').upper()
        print(f"  - {chr_name}: {diff:.3f}% difference between groups")
    
    print(f"\n✅ RECOMMENDATION:")
    if dup_2510.mean() < 12 and all(d < 1.0 for d in chr_diff_2510.values()):
        print(f"  → Use 2510+ samples ({len(period3)} samples)")
        print(f"  → Low duplication, minimal chr proportion differences")
        print(f"  → Both groups can potentially be combined")
    elif dup_2510.mean() < 15:
        print(f"  → Use 2510+ samples with moderate filtering ({len(period3)} samples)")
        print(f"  → Filter out dup > 15% if needed")
        print(f"  → Consider separating groups if Z-score validation fails")
    else:
        print(f"  → Use 2510+ samples with strict filtering")
        print(f"  → Separate Group 1 and Group 2")
    
    print(f"\n{'='*80}\n")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python3 compare_time_periods.py <output_dir> <lab_name>")
        sys.exit(1)
    
    compare_time_periods(sys.argv[1], sys.argv[2])
