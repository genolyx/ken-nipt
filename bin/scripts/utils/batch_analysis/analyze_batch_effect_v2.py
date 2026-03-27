#!/usr/bin/env python3
"""
Batch Effect 분석 스크립트 v2
fastq.gz 심볼릭 링크의 날짜를 기반으로 실제 batch 구분
GNCI (Cordlife) vs GNMF (UCL) 별도 분석
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import sys
import glob
from datetime import datetime
from collections import defaultdict

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (15, 10)

def get_fastq_upload_date(analysis_dir, sample_id):
    """
    fastq 디렉토리의 생성 날짜 가져오기 (실제 업로드 날짜)
    
    Args:
        analysis_dir: 분석 디렉토리 경로 (월 추출용)
        sample_id: 샘플 ID (e.g., GNMF25070009, GNCI25120001)
    
    Returns:
        날짜 문자열 (YYYY-MM-DD) 또는 None
    """
    # 월 추출
    month = os.path.basename(analysis_dir.rstrip("/"))
    
    # fastq 디렉토리 경로 후보들
    fastq_paths = [
        f"/home/ken/ken-nipt/fastq/{month}/{sample_id}",
        f"/data/fastq_backup/{month}/{sample_id}"
    ]
    
    for fastq_dir in fastq_paths:
        if os.path.exists(fastq_dir):
            try:
                # 디렉토리의 수정 시간 가져오기
                stat_info = os.stat(fastq_dir)
                mtime = datetime.fromtimestamp(stat_info.st_mtime)
                return mtime.strftime('%Y-%m-%d')
            except Exception as e:
                continue
    
    return None

def parse_10mb_wig(wig_file):
    """
    10mb.wig.Normalization.txt 파일 파싱
    
    Returns:
        dict: {
            'chr_reads': {chr: total_reads},
            'coverage_cv': coefficient of variation,
            'valid_bins': number of valid bins
        }
    """
    try:
        df = pd.read_csv(wig_file, sep='\t')
        
        # Valid bins만 사용 (TRUE인 것들) - 대소문자 구분 없이
        valid_df = df[df['valid'].astype(str).str.upper() == 'TRUE'].copy()
        
        if len(valid_df) == 0:
            return None
        
        # Chromosome별 read count
        chr_reads = valid_df.groupby('chr')['reads'].sum().to_dict()
        
        # Coverage uniformity (CV = std/mean)
        reads = valid_df['reads'].values
        coverage_mean = np.mean(reads)
        coverage_std = np.std(reads)
        coverage_cv = coverage_std / coverage_mean if coverage_mean > 0 else np.nan
        
        return {
            'chr_reads': chr_reads,
            'coverage_cv': coverage_cv,
            'valid_bins': len(valid_df),
            'total_reads': valid_df['reads'].sum()
        }
    except Exception as e:
        print(f"  ⚠️  Error parsing {wig_file}: {e}")
        return None

def extract_chromosome_metrics(analysis_base_dir, df):
    """
    각 샘플의 chromosome-level metrics 추출
    
    Args:
        analysis_base_dir: analysis 디렉토리 경로
        df: 샘플 DataFrame
    
    Returns:
        DataFrame with additional columns for chr proportions and coverage_cv
    """
    print("\nExtracting chromosome-level metrics from 10mb.wig files...")
    
    chr_data_list = []
    missing_files = []
    
    for idx, row in df.iterrows():
        sample_id = row['sample_id']
        
        # Use sample_dir if available, otherwise construct path
        if 'sample_dir' in row and pd.notna(row['sample_dir']):
            sample_dir = os.path.join(row['sample_dir'], "Output_hmmcopy")
        else:
            month = row['month']
            sample_dir = os.path.join(analysis_base_dir, str(month), sample_id, "Output_hmmcopy")
        
        # 10mb.wig 파일 찾기 (of_orig 우선, 없으면 proper_paired, of_fetus)
        wig_patterns = [
            f"{sample_id}.of_orig.10mb.wig.Normalization.txt",
            f"{sample_id}.proper_paired.10mb.wig.Normalization.txt",
            f"{sample_id}.of_fetus.10mb.wig.Normalization.txt"
        ]
        
        wig_file = None
        for pattern in wig_patterns:
            candidate = os.path.join(sample_dir, pattern)
            if os.path.exists(candidate):
                wig_file = candidate
                break
        
        if not wig_file:
            chr_data_list.append(None)
            missing_files.append(sample_id)
            continue
        
        # Parse wig file
        metrics = parse_10mb_wig(wig_file)
        chr_data_list.append(metrics)
        
        if (idx + 1) % 50 == 0:
            print(f"  Processed {idx + 1}/{len(df)} samples...")
    
    # Extract chromosome proportions
    chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    
    for chrom in chromosomes:
        df[f'{chrom}_prop'] = [
            data['chr_reads'].get(chrom, 0) / data['total_reads'] if data and data['total_reads'] > 0 else np.nan
            for data in chr_data_list
        ]
    
    # Coverage CV
    df['coverage_cv'] = [data['coverage_cv'] if data else np.nan for data in chr_data_list]
    df['valid_bins'] = [data['valid_bins'] if data else np.nan for data in chr_data_list]
    
    print(f"\n✓ Chromosome metrics extraction complete!")
    print(f"  Total samples: {len(df)}")
    print(f"  Samples with wig data: {df['coverage_cv'].notna().sum()}")
    print(f"  Samples missing wig data: {len(missing_files)}")
    if len(missing_files) > 0 and len(missing_files) <= 20:
        print(f"  Missing samples: {', '.join(missing_files)}")
    elif len(missing_files) > 20:
        print(f"  First 20 missing samples: {', '.join(missing_files[:20])}")
    
    return df

def extract_batch_info(sample_list_file, analysis_base_dir="/home/ken/ken-nipt/analysis"):
    """
    샘플 리스트에서 batch 정보 추출
    
    Returns:
        DataFrame with additional columns: lab, upload_date, batch_id
    """
    print(f"Loading: {sample_list_file}")
    df = pd.read_csv(sample_list_file, sep='\t')
    
    print(f"\nExtracting batch information from fastq.gz timestamps...")
    
    # Lab 구분 (GNCI: Cordlife, GNMF: UCL)
    df['lab'] = df['sample_id'].apply(lambda x: 'Cordlife' if x.startswith('GNCI') else 'UCL')
    
    # Upload date 추출
    upload_dates = []
    for _, row in df.iterrows():
        month = row['month']
        sample_id = row['sample_id']
        analysis_dir = os.path.join(analysis_base_dir, str(month))
        
        upload_date = get_fastq_upload_date(analysis_dir, sample_id)
        upload_dates.append(upload_date)
        
        if len(upload_dates) % 50 == 0:
            print(f"  Processed {len(upload_dates)}/{len(df)} samples...")
    
    df['upload_date'] = upload_dates
    
    # Batch ID 생성 (month + upload_date, upload_date가 없으면 month만 사용)
    df['batch_id'] = df.apply(
        lambda x: f"{x['month']}_{x['upload_date']}" if pd.notna(x['upload_date']) else str(x['month']),
        axis=1
    )
    
    # Data type conversion
    numeric_cols = ['number_of_reads', 'number_of_mapped_reads', 'mapping_rate(%)', 
                    'duplication_rate(%)', 'mean_mapping_quality', 'mean_coverageData(X)',
                    'GC_content(%)', 'SeqFF', 'Fragment_FF', 'YFF_2', 'M-SeqFF']
    
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Summary
    print(f"\n✓ Batch extraction complete!")
    print(f"  Total samples: {len(df)}")
    print(f"  UCL samples: {(df['lab'] == 'UCL').sum()}")
    print(f"  Cordlife samples: {(df['lab'] == 'Cordlife').sum()}")
    print(f"  Unique batches: {df['batch_id'].nunique()}")
    print(f"  Samples with upload_date: {df['upload_date'].notna().sum()}")
    
    return df

def analyze_batch_structure(df, output_dir):
    """Batch 구조 분석"""
    print("\n" + "="*60)
    print("1. Batch Structure Analysis")
    print("="*60)
    
    # Batch별 샘플 수
    print("\nSamples per batch:")
    batch_counts = df['batch_id'].value_counts().sort_index()
    for batch_id, count in batch_counts.items():
        lab_dist = df[df['batch_id'] == batch_id]['lab'].value_counts()
        print(f"  {batch_id}: {count} samples ({dict(lab_dist)})")
    
    # 시각화
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    
    # 1. Batch별 샘플 수
    ax = axes[0, 0]
    batch_counts_top20 = batch_counts.head(20)
    batch_counts_top20.plot(kind='bar', ax=ax, color='skyblue')
    ax.set_xlabel('Batch ID', fontsize=12)
    ax.set_ylabel('Number of Samples', fontsize=12)
    ax.set_title('Samples per Batch (Top 20)', fontsize=14, fontweight='bold')
    ax.tick_params(axis='x', rotation=45)
    ax.grid(True, alpha=0.3)
    
    # 2. Lab별 batch 분포
    ax = axes[0, 1]
    lab_batch_counts = df.groupby(['lab', 'batch_id']).size().unstack(fill_value=0)
    lab_batch_counts.T.plot(kind='bar', stacked=True, ax=ax, color=['lightblue', 'lightcoral'])
    ax.set_xlabel('Batch ID', fontsize=12)
    ax.set_ylabel('Number of Samples', fontsize=12)
    ax.set_title('Batch Composition (Lab)', fontsize=14, fontweight='bold')
    ax.legend(title='Lab')
    ax.tick_params(axis='x', rotation=45)
    ax.grid(True, alpha=0.3)
    
    # 3. 월별 batch 수
    ax = axes[1, 0]
    month_batch_counts = df.groupby('month')['batch_id'].nunique().sort_index()
    month_batch_counts.plot(kind='bar', ax=ax, color='lightgreen')
    ax.set_xlabel('Month', fontsize=12)
    ax.set_ylabel('Number of Batches', fontsize=12)
    ax.set_title('Number of Batches per Month', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # 4. Upload date 분포 (시계열)
    ax = axes[1, 1]
    df_dated = df[df['upload_date'].notna()].copy()
    df_dated['upload_date_dt'] = pd.to_datetime(df_dated['upload_date'])
    daily_counts = df_dated.groupby(['upload_date_dt', 'lab']).size().unstack(fill_value=0)
    daily_counts.plot(ax=ax, marker='o', linewidth=2)
    ax.set_xlabel('Upload Date', fontsize=12)
    ax.set_ylabel('Number of Samples', fontsize=12)
    ax.set_title('Sample Upload Timeline', fontsize=14, fontweight='bold')
    ax.legend(title='Lab')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '01_batch_structure.png'), dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/01_batch_structure.png")
    plt.close()

def analyze_qc_by_batch(df, output_dir, top_n_batches=15):
    """Batch별 QC metrics 비교"""
    print("\n" + "="*60)
    print("2. QC Metrics by Batch")
    print("="*60)
    
    # Top N batches (샘플 수 기준)
    top_batches = df['batch_id'].value_counts().head(top_n_batches).index.tolist()
    df_top = df[df['batch_id'].isin(top_batches)]
    
    metrics = ['mapping_rate(%)', 'duplication_rate(%)', 'GC_content(%)', 'SeqFF']
    
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    axes = axes.flatten()
    
    for idx, metric in enumerate(metrics):
        if metric not in df.columns:
            continue
        
        ax = axes[idx]
        
        # Box plot
        df_clean = df_top.dropna(subset=[metric])
        batches = sorted(df_clean['batch_id'].unique())
        data_by_batch = [df_clean[df_clean['batch_id'] == b][metric].values for b in batches]
        
        bp = ax.boxplot(data_by_batch, labels=batches, patch_artist=True)
        
        # Color by batch
        colors = plt.cm.Set3(np.linspace(0, 1, len(batches)))
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
        
        ax.set_xlabel('Batch ID', fontsize=12)
        ax.set_ylabel(metric, fontsize=12)
        ax.set_title(f'{metric} by Batch (Top {top_n_batches})', fontsize=14, fontweight='bold')
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, alpha=0.3)
        
        # Statistics
        print(f"\n{metric}:")
        for batch in batches:
            batch_data = df_clean[df_clean['batch_id'] == batch][metric]
            if len(batch_data) > 0:
                print(f"  {batch}: mean={batch_data.mean():.2f} ± {batch_data.std():.2f}, "
                      f"n={len(batch_data)}")
        
        # Kruskal-Wallis test
        if len(batches) > 1:
            h_stat, p_value = stats.kruskal(*data_by_batch)
            print(f"  Kruskal-Wallis test: H={h_stat:.2f}, p={p_value:.4f}")
            if p_value < 0.05:
                print(f"  *** Significant batch effect detected (p < 0.05) ***")
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '02_qc_metrics_by_batch.png'), dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/02_qc_metrics_by_batch.png")
    plt.close()

def compare_ucl_vs_cordlife(df, output_dir):
    """UCL vs Cordlife 비교"""
    print("\n" + "="*60)
    print("3. UCL vs Cordlife Comparison")
    print("="*60)
    
    metrics = ['mapping_rate(%)', 'duplication_rate(%)', 'mean_mapping_quality',
               'GC_content(%)', 'SeqFF', 'mean_coverageData(X)']
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    for idx, metric in enumerate(metrics):
        if metric not in df.columns:
            continue
        
        ax = axes[idx]
        
        # Violin plot
        df_clean = df.dropna(subset=[metric])
        
        parts = ax.violinplot(
            [df_clean[df_clean['lab'] == 'UCL'][metric].values,
             df_clean[df_clean['lab'] == 'Cordlife'][metric].values],
            positions=[1, 2],
            showmeans=True,
            showmedians=True
        )
        
        ax.set_xticks([1, 2])
        ax.set_xticklabels(['UCL', 'Cordlife'])
        ax.set_ylabel(metric, fontsize=12)
        ax.set_title(f'{metric}: UCL vs Cordlife', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Statistics
        ucl_data = df_clean[df_clean['lab'] == 'UCL'][metric]
        cordlife_data = df_clean[df_clean['lab'] == 'Cordlife'][metric]
        
        print(f"\n{metric}:")
        print(f"  UCL: mean={ucl_data.mean():.2f} ± {ucl_data.std():.2f}, "
              f"median={ucl_data.median():.2f}, n={len(ucl_data)}")
        print(f"  Cordlife: mean={cordlife_data.mean():.2f} ± {cordlife_data.std():.2f}, "
              f"median={cordlife_data.median():.2f}, n={len(cordlife_data)}")
        
        # Mann-Whitney U test
        if len(ucl_data) > 0 and len(cordlife_data) > 0:
            u_stat, p_value = stats.mannwhitneyu(ucl_data, cordlife_data, alternative='two-sided')
            print(f"  Mann-Whitney U test: U={u_stat:.2f}, p={p_value:.4f}")
            if p_value < 0.05:
                print(f"  *** Significant difference between labs (p < 0.05) ***")
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '03_ucl_vs_cordlife.png'), dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/03_ucl_vs_cordlife.png")
    plt.close()

def analyze_batch_effect_over_time(df, output_dir):
    """시계열 batch effect 분석"""
    print("\n" + "="*60)
    print("4. Batch Effect Over Time")
    print("="*60)
    
    df_dated = df[df['upload_date'].notna()].copy()
    df_dated['upload_date_dt'] = pd.to_datetime(df_dated['upload_date'])
    df_dated = df_dated.sort_values('upload_date_dt')
    
    metrics = ['mapping_rate(%)', 'duplication_rate(%)', 'GC_content(%)', 'SeqFF']
    
    fig, axes = plt.subplots(2, 2, figsize=(18, 12))
    axes = axes.flatten()
    
    for idx, metric in enumerate(metrics):
        if metric not in df.columns:
            continue
        
        ax = axes[idx]
        
        # Daily average by lab
        daily_avg = df_dated.groupby(['upload_date_dt', 'lab'])[metric].mean().unstack(fill_value=np.nan)
        
        for lab in daily_avg.columns:
            ax.plot(daily_avg.index, daily_avg[lab], marker='o', label=lab, linewidth=2, markersize=6)
        
        ax.set_xlabel('Upload Date', fontsize=12)
        ax.set_ylabel(metric, fontsize=12)
        ax.set_title(f'{metric} Over Time', fontsize=14, fontweight='bold')
        ax.legend(title='Lab')
        ax.grid(True, alpha=0.3)
        
        # Add rolling mean (7-day window)
        for lab in daily_avg.columns:
            rolling = daily_avg[lab].rolling(window=7, min_periods=3).mean()
            ax.plot(rolling.index, rolling, linestyle='--', linewidth=2, alpha=0.6)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '04_batch_effect_over_time.png'), dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/04_batch_effect_over_time.png")
    plt.close()

def analyze_chromosome_distribution(df, output_dir):
    """Chromosome-level read distribution 분석"""
    print("\n" + "="*60)
    print("5. Chromosome-level Read Distribution Analysis")
    print("="*60)
    
    # Check if chromosome data exists
    chr_cols = [col for col in df.columns if col.endswith('_prop')]
    if not chr_cols:
        print("⚠️  No chromosome proportion data found. Skipping.")
        return
    
    df_clean = df[df['coverage_cv'].notna()].copy()
    if len(df_clean) == 0:
        print("⚠️  No valid chromosome data. Skipping.")
        return
    
    print(f"\nSamples with chromosome data: {len(df_clean)}")
    
    # 1. Key chromosome proportions by batch
    key_chrs = ['chr13_prop', 'chr18_prop', 'chr21_prop', 'chrX_prop', 'chrY_prop']
    
    fig, axes = plt.subplots(3, 2, figsize=(18, 16))
    axes = axes.flatten()
    
    for idx, chr_col in enumerate(key_chrs):
        if chr_col not in df_clean.columns:
            continue
        
        ax = axes[idx]
        
        # Violin plot by batch (top 10 batches by sample size)
        top_batches = df_clean['batch_id'].value_counts().head(10).index
        df_plot = df_clean[df_clean['batch_id'].isin(top_batches)]
        
        if len(df['lab'].unique()) > 1:
            sns.violinplot(data=df_plot, x='batch_id', y=chr_col, hue='lab', ax=ax, cut=0)
        else:
            sns.violinplot(data=df_plot, x='batch_id', y=chr_col, ax=ax, cut=0)
        
        ax.set_xlabel('Batch', fontsize=11)
        ax.set_ylabel(f'{chr_col.replace("_prop", "")} proportion', fontsize=11)
        ax.set_title(f'{chr_col.replace("_prop", "")} Read Proportion by Batch', 
                     fontsize=13, fontweight='bold')
        ax.tick_params(axis='x', rotation=90)
        ax.grid(True, alpha=0.3, axis='y')
        
        # Kruskal-Wallis test
        batch_groups = [group[chr_col].dropna().values 
                       for _, group in df_clean.groupby('batch_id')]
        batch_groups = [g for g in batch_groups if len(g) > 0]
        
        if len(batch_groups) >= 2:
            try:
                h_stat, p_value = stats.kruskal(*batch_groups)
                ax.text(0.02, 0.98, f'Kruskal-Wallis p={p_value:.4f}',
                       transform=ax.transAxes, fontsize=10,
                       verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            except ValueError as e:
                # All values are identical
                ax.text(0.02, 0.98, 'All values identical',
                       transform=ax.transAxes, fontsize=10,
                       verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))
    
    # Coverage CV plot
    ax = axes[5]
    top_batches = df_clean['batch_id'].value_counts().head(10).index
    df_plot = df_clean[df_clean['batch_id'].isin(top_batches)]
    
    if len(df['lab'].unique()) > 1:
        sns.violinplot(data=df_plot, x='batch_id', y='coverage_cv', hue='lab', ax=ax, cut=0)
    else:
        sns.violinplot(data=df_plot, x='batch_id', y='coverage_cv', ax=ax, cut=0)
    
    ax.set_xlabel('Batch', fontsize=11)
    ax.set_ylabel('Coverage CV', fontsize=11)
    ax.set_title('Coverage Uniformity (CV) by Batch', fontsize=13, fontweight='bold')
    ax.tick_params(axis='x', rotation=90)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '05_chromosome_distribution.png'), dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/05_chromosome_distribution.png")
    plt.close()
    
    # 2. Coverage CV statistics
    print("\nCoverage Uniformity (CV) Statistics:")
    for batch in df_clean['batch_id'].value_counts().head(10).index:
        batch_data = df_clean[df_clean['batch_id'] == batch]['coverage_cv']
        print(f"  {batch}: mean={batch_data.mean():.4f} ± {batch_data.std():.4f}, "
              f"median={batch_data.median():.4f}, n={len(batch_data)}")
    
    # 3. Lab comparison if multiple labs
    if len(df_clean['lab'].unique()) > 1:
        print("\nCoverage CV by Lab:")
        for lab in sorted(df_clean['lab'].unique()):
            lab_data = df_clean[df_clean['lab'] == lab]['coverage_cv']
            print(f"  {lab}: mean={lab_data.mean():.4f} ± {lab_data.std():.4f}, "
                  f"median={lab_data.median():.4f}, n={len(lab_data)}")
        
        # Mann-Whitney U test
        ucl_cv = df_clean[df_clean['lab'] == 'UCL']['coverage_cv'].dropna()
        cordlife_cv = df_clean[df_clean['lab'] == 'Cordlife']['coverage_cv'].dropna()
        if len(ucl_cv) > 0 and len(cordlife_cv) > 0:
            u_stat, p_value = stats.mannwhitneyu(ucl_cv, cordlife_cv, alternative='two-sided')
            print(f"  Mann-Whitney U test: U={u_stat:.2f}, p={p_value:.4f}")
            if p_value < 0.05:
                print(f"  *** Significant difference (p < 0.05) ***")

def analyze_pca_with_batch(df, output_dir):
    """PCA 분석 (batch 색상 구분)"""
    print("\n" + "="*60)
    print("6. PCA Analysis with Batch Coloring")
    print("="*60)
    
    try:
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
    except ImportError:
        print("⚠ sklearn not installed. Skipping PCA analysis.")
        print("  To enable PCA: pip install scikit-learn")
        return None
    
    metrics = ['mapping_rate(%)', 'duplication_rate(%)', 'GC_content(%)', 
               'SeqFF', 'mean_coverageData(X)', 'mean_mapping_quality']
    
    # Add chromosome metrics if available
    chr_metrics = ['coverage_cv', 'chr13_prop', 'chr18_prop', 'chr21_prop', 'chrX_prop', 'chrY_prop']
    available_chr_metrics = [m for m in chr_metrics if m in df.columns]
    
    if available_chr_metrics:
        print(f"Including chromosome metrics: {', '.join(available_chr_metrics)}")
        metrics = metrics + available_chr_metrics
    
    # Prepare data
    df_pca = df[['batch_id', 'lab'] + metrics].dropna()
    
    if len(df_pca) < 50:
        print("Not enough samples for PCA analysis")
        return None
    
    X = df_pca[metrics].values
    batches = df_pca['batch_id'].values
    labs = df_pca['lab'].values
    
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    
    # Plot by batch (top 20) and lab
    unique_labs = df_pca['lab'].unique()
    
    if len(unique_labs) > 1:
        # Both batch and lab plots
        fig, axes = plt.subplots(1, 2, figsize=(20, 8))
        
        # Plot 1: By Batch (top 20)
        ax = axes[0]
        top_batches = df_pca['batch_id'].value_counts().head(20).index
        
        for batch in top_batches:
            mask = batches == batch
            ax.scatter(X_pca[mask, 0], X_pca[mask, 1], 
                      label=batch, alpha=0.6, s=50)
        
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=12)
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=12)
        ax.set_title('PCA: Batch Effect (Top 20 Batches)', fontsize=14, fontweight='bold')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        ax.grid(True, alpha=0.3)
        
        # Plot 2: By Lab
        ax = axes[1]
        colors = {'UCL': 'blue', 'Cordlife': 'red'}
        for lab in unique_labs:
            mask = labs == lab
            ax.scatter(X_pca[mask, 0], X_pca[mask, 1], 
                      c=colors.get(lab, 'gray'), label=lab, alpha=0.6, s=50)
        
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=12)
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=12)
        ax.set_title('PCA: Lab Effect', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
    else:
        # Single lab: only batch plot (larger)
        fig, ax = plt.subplots(1, 1, figsize=(14, 10))
        
        top_batches = df_pca['batch_id'].value_counts().head(20).index
        colors = plt.cm.Set3(np.linspace(0, 1, len(top_batches)))
        
        for batch, color in zip(top_batches, colors):
            mask = batches == batch
            ax.scatter(X_pca[mask, 0], X_pca[mask, 1], 
                      c=[color], label=batch, alpha=0.7, s=80, edgecolors='black', linewidth=0.5)
        
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=14)
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=14)
        ax.set_title(f'PCA: Batch Effect - {unique_labs[0]}', fontsize=16, fontweight='bold')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, ncol=2)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '06_pca_batch_and_lab.png'), dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/06_pca_batch_and_lab.png")
    
    print("\nPCA Results:")
    print(f"  PC1 explains {pca.explained_variance_ratio_[0]*100:.1f}% of variance")
    print(f"  PC2 explains {pca.explained_variance_ratio_[1]*100:.1f}% of variance")
    print(f"  Total explained: {sum(pca.explained_variance_ratio_)*100:.1f}%")
    
    # Feature importance
    print("\nFeature contributions to PC1:")
    feature_importance = pd.DataFrame({
        'Feature': metrics,
        'PC1': np.abs(pca.components_[0]),
        'PC2': np.abs(pca.components_[1])
    }).sort_values('PC1', ascending=False)
    
    for _, row in feature_importance.iterrows():
        print(f"  {row['Feature']}: {row['PC1']:.3f}")
    
    plt.close()
    
    # Return PCA results for report
    return {
        'explained_variance': pca.explained_variance_ratio_,
        'feature_importance': feature_importance,
        'n_samples': len(df_pca)
    }

def analyze_pca_chromosome_only(df, output_dir):
    """PCA 분석 - Chromosome proportions만 사용 (chr13, chr18, chr21)"""
    print("\n" + "="*60)
    print("7. PCA Analysis - Chromosome Proportions Only")
    print("="*60)
    
    try:
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
    except ImportError:
        print("⚠ sklearn not installed. Skipping PCA analysis.")
        return None
    
    # Chr13, 18, 21만 사용 (성별 영향 배제)
    chr_metrics = ['chr13_prop', 'chr18_prop', 'chr21_prop']
    
    # Check if data exists
    available_metrics = [m for m in chr_metrics if m in df.columns and df[m].notna().sum() > 0]
    
    if len(available_metrics) < 3:
        print(f"⚠️  Not enough chromosome data. Available: {available_metrics}")
        return None
    
    print(f"Using chromosome proportions: {', '.join(available_metrics)}")
    
    # Prepare data
    df_pca = df[['batch_id', 'lab'] + available_metrics].dropna()
    
    if len(df_pca) < 50:
        print(f"Not enough samples for PCA analysis (n={len(df_pca)})")
        return None
    
    X = df_pca[available_metrics].values
    batches = df_pca['batch_id'].values
    labs = df_pca['lab'].values
    
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    
    # Plot
    unique_labs = df_pca['lab'].unique()
    
    if len(unique_labs) > 1:
        # Both batch and lab plots
        fig, axes = plt.subplots(1, 2, figsize=(20, 8))
        
        # Plot 1: By Batch (top 20)
        ax = axes[0]
        top_batches = df_pca['batch_id'].value_counts().head(20).index
        colors = plt.cm.Set3(np.linspace(0, 1, len(top_batches)))
        
        for batch, color in zip(top_batches, colors):
            mask = batches == batch
            ax.scatter(X_pca[mask, 0], X_pca[mask, 1], 
                      c=[color], label=batch, alpha=0.7, s=80, edgecolors='black', linewidth=0.5)
        
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=12)
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=12)
        ax.set_title('PCA: Batch Effect (Chr13/18/21 proportions)', fontsize=14, fontweight='bold')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=2)
        ax.grid(True, alpha=0.3)
        
        # Plot 2: By Lab
        ax = axes[1]
        colors_lab = {'UCL': 'blue', 'Cordlife': 'red'}
        for lab in unique_labs:
            mask = labs == lab
            ax.scatter(X_pca[mask, 0], X_pca[mask, 1], 
                      c=colors_lab.get(lab, 'gray'), label=lab, alpha=0.6, s=50)
        
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=12)
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=12)
        ax.set_title('PCA: Lab Effect (Chr13/18/21 proportions)', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
    else:
        # Single lab: only batch plot
        fig, ax = plt.subplots(1, 1, figsize=(14, 10))
        
        top_batches = df_pca['batch_id'].value_counts().head(20).index
        colors = plt.cm.Set3(np.linspace(0, 1, len(top_batches)))
        
        for batch, color in zip(top_batches, colors):
            mask = batches == batch
            ax.scatter(X_pca[mask, 0], X_pca[mask, 1], 
                      c=[color], label=batch, alpha=0.7, s=80, edgecolors='black', linewidth=0.5)
        
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=14)
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=14)
        ax.set_title(f'PCA: Batch Effect - {unique_labs[0]} (Chr13/18/21 proportions)', 
                     fontsize=16, fontweight='bold')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, ncol=2)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '07_pca_chromosome_only.png'), dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/07_pca_chromosome_only.png")
    
    print("\nPCA Results (Chromosome proportions only):")
    print(f"  PC1 explains {pca.explained_variance_ratio_[0]*100:.1f}% of variance")
    print(f"  PC2 explains {pca.explained_variance_ratio_[1]*100:.1f}% of variance")
    print(f"  Total explained: {sum(pca.explained_variance_ratio_)*100:.1f}%")
    
    # Feature importance
    print("\nFeature contributions to PC1:")
    feature_importance = pd.DataFrame({
        'Feature': available_metrics,
        'PC1': np.abs(pca.components_[0]),
        'PC2': np.abs(pca.components_[1])
    }).sort_values('PC1', ascending=False)
    
    for _, row in feature_importance.iterrows():
        print(f"  {row['Feature']}: PC1={row['PC1']:.3f}, PC2={row['PC2']:.3f}")
    
    plt.close()
    
    # Return PCA results
    return {
        'explained_variance': pca.explained_variance_ratio_,
        'feature_importance': feature_importance,
        'n_samples': len(df_pca)
    }

def analyze_pca_qc_only(df, output_dir):
    """PCA 분석 - QC metrics만 사용 (5개)"""
    print("\n" + "="*60)
    print("8. PCA Analysis - QC Metrics Only")
    print("="*60)
    
    try:
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
    except ImportError:
        print("⚠ sklearn not installed. Skipping PCA analysis.")
        return None
    
    # QC metrics만 사용
    qc_metrics = ['GC_content(%)', 'mapping_rate(%)', 'duplication_rate(%)', 
                  'mean_mapping_quality', 'mean_coverageData(X)']
    
    # Check if data exists
    available_metrics = [m for m in qc_metrics if m in df.columns and df[m].notna().sum() > 0]
    
    if len(available_metrics) < 5:
        print(f"⚠️  Not enough QC data. Available: {available_metrics}")
        return None
    
    print(f"Using QC metrics: {', '.join(available_metrics)}")
    
    # Prepare data
    df_pca = df[['batch_id', 'lab'] + available_metrics].dropna()
    
    if len(df_pca) < 50:
        print(f"Not enough samples for PCA analysis (n={len(df_pca)})")
        return None
    
    X = df_pca[available_metrics].values
    batches = df_pca['batch_id'].values
    labs = df_pca['lab'].values
    
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    
    # Plot
    unique_labs = df_pca['lab'].unique()
    
    if len(unique_labs) > 1:
        # Both batch and lab plots
        fig, axes = plt.subplots(1, 2, figsize=(20, 8))
        
        # Plot 1: By Batch (top 20)
        ax = axes[0]
        top_batches = df_pca['batch_id'].value_counts().head(20).index
        colors = plt.cm.Set3(np.linspace(0, 1, len(top_batches)))
        
        for batch, color in zip(top_batches, colors):
            mask = batches == batch
            ax.scatter(X_pca[mask, 0], X_pca[mask, 1], 
                      c=[color], label=batch, alpha=0.7, s=80, edgecolors='black', linewidth=0.5)
        
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=12)
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=12)
        ax.set_title('PCA: Batch Effect (QC metrics only)', fontsize=14, fontweight='bold')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=2)
        ax.grid(True, alpha=0.3)
        
        # Plot 2: By Lab
        ax = axes[1]
        colors_lab = {'UCL': 'blue', 'Cordlife': 'red'}
        for lab in unique_labs:
            mask = labs == lab
            ax.scatter(X_pca[mask, 0], X_pca[mask, 1], 
                      c=colors_lab.get(lab, 'gray'), label=lab, alpha=0.6, s=50)
        
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=12)
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=12)
        ax.set_title('PCA: Lab Effect (QC metrics only)', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
    else:
        # Single lab: only batch plot
        fig, ax = plt.subplots(1, 1, figsize=(14, 10))
        
        top_batches = df_pca['batch_id'].value_counts().head(20).index
        colors = plt.cm.Set3(np.linspace(0, 1, len(top_batches)))
        
        for batch, color in zip(top_batches, colors):
            mask = batches == batch
            ax.scatter(X_pca[mask, 0], X_pca[mask, 1], 
                      c=[color], label=batch, alpha=0.7, s=80, edgecolors='black', linewidth=0.5)
        
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=14)
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=14)
        ax.set_title(f'PCA: Batch Effect - {unique_labs[0]} (QC metrics only)', 
                     fontsize=16, fontweight='bold')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, ncol=2)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '08_pca_qc_only.png'), dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: {output_dir}/08_pca_qc_only.png")
    
    print("\nPCA Results (QC metrics only):")
    print(f"  PC1 explains {pca.explained_variance_ratio_[0]*100:.1f}% of variance")
    print(f"  PC2 explains {pca.explained_variance_ratio_[1]*100:.1f}% of variance")
    print(f"  Total explained: {sum(pca.explained_variance_ratio_)*100:.1f}%")
    
    # Feature importance
    print("\nFeature contributions to PC1:")
    feature_importance = pd.DataFrame({
        'Feature': available_metrics,
        'PC1': np.abs(pca.components_[0]),
        'PC2': np.abs(pca.components_[1])
    }).sort_values('PC1', ascending=False)
    
    for _, row in feature_importance.iterrows():
        print(f"  {row['Feature']}: PC1={row['PC1']:.3f}, PC2={row['PC2']:.3f}")
    
    plt.close()
    
    # Return PCA results
    return {
        'explained_variance': pca.explained_variance_ratio_,
        'feature_importance': feature_importance,
        'n_samples': len(df_pca)
    }

def generate_batch_summary_report(df, output_dir, lab_filter='all', pca_results=None, pca_chr_results=None, pca_qc_results=None):
    """Batch 분석 요약 리포트"""
    print("\n" + "="*60)
    print("7. Generating Batch Summary Report")
    print("="*60)
    
    report_file = os.path.join(output_dir, 'batch_effect_summary.txt')
    
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write(f"BATCH EFFECT ANALYSIS REPORT - {lab_filter.upper()}\n")
        f.write("="*80 + "\n\n")
        
        # Overall statistics
        f.write("1. OVERALL STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total samples: {len(df)}\n")
        f.write(f"  UCL: {(df['lab'] == 'UCL').sum()}\n")
        f.write(f"  Cordlife: {(df['lab'] == 'Cordlife').sum()}\n")
        f.write(f"Unique batches: {df['batch_id'].nunique()}\n")
        f.write(f"Months covered: {', '.join(map(str, sorted(df['month'].unique())))}\n\n")
        
        # Data completeness
        f.write("2. DATA COMPLETENESS\n")
        f.write("-" * 80 + "\n")
        
        qc_metrics = {
            'mapping_rate(%)': 'Mapping Rate',
            'duplication_rate(%)': 'Duplication Rate',
            'GC_content(%)': 'GC Content',
            'SeqFF': 'Fetal Fraction',
            'coverage_cv': 'Coverage Uniformity (10mb.wig)',
            'chr13_prop': 'Chromosome proportions (10mb.wig)'
        }
        
        for col, desc in qc_metrics.items():
            if col in df.columns:
                available = df[col].notna().sum()
                pct = available / len(df) * 100
                f.write(f"{desc:40s}: {available:4d}/{len(df):4d} ({pct:5.1f}%)\n")
        f.write("\n")
        
        # Batch sizes
        f.write("3. BATCH SIZES (Top 20)\n")
        f.write("-" * 80 + "\n")
        batch_counts = df['batch_id'].value_counts().head(20)
        for batch_id, count in batch_counts.items():
            lab_dist = df[df['batch_id'] == batch_id]['lab'].value_counts()
            f.write(f"{batch_id}: {count} samples {dict(lab_dist)}\n")
        f.write("\n")
        
        # QC metrics summary (only if single lab)
        unique_labs = df['lab'].unique()
        
        if len(unique_labs) == 1:
            f.write("4. QC METRICS SUMMARY (by Month)\n")
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
            
            # Chromosome metrics if available
            if 'coverage_cv' in df.columns and df['coverage_cv'].notna().sum() > 0:
                f.write(f"\ncoverage_cv:\n")
                for month in sorted(df['month'].unique()):
                    month_data = df[df['month'] == month]['coverage_cv'].dropna()
                    if len(month_data) > 0:
                        f.write(f"  {month}: mean={month_data.mean():.4f} ± {month_data.std():.4f}, "
                               f"median={month_data.median():.4f}, n={len(month_data)}\n")
                
                # Batch effect test for coverage_cv
                batch_groups_cv = [group['coverage_cv'].dropna().values 
                                   for _, group in df.groupby('batch_id')]
                batch_groups_cv = [g for g in batch_groups_cv if len(g) > 0]
                if len(batch_groups_cv) >= 2:
                    try:
                        h_stat, p_value = stats.kruskal(*batch_groups_cv)
                        f.write(f"  Kruskal-Wallis test: H={h_stat:.2f}, p={p_value:.4f}\n")
                        if p_value < 0.05:
                            f.write(f"  *** SIGNIFICANT BATCH EFFECT (p < 0.05) ***\n")
                    except ValueError:
                        f.write(f"  Kruskal-Wallis test: All values identical\n")
            
            # Key chromosome proportions
            chr_metrics = ['chr13_prop', 'chr18_prop', 'chr21_prop', 'chrX_prop', 'chrY_prop']
            available_chrs = [c for c in chr_metrics if c in df.columns and df[c].notna().sum() > 0]
            
            if available_chrs:
                f.write(f"\nChromosome Read Proportions (by Month, mean ± std):\n")
                for month in sorted(df['month'].unique()):
                    f.write(f"\n{month}:\n")
                    for chr_col in available_chrs:
                        month_data = df[df['month'] == month][chr_col].dropna()
                        if len(month_data) > 0:
                            f.write(f"  {chr_col.replace('_prop', '')}: {month_data.mean():.6f} ± {month_data.std():.6f} (n={len(month_data)})\n")
                
                # Batch effect tests for key chromosomes
                f.write(f"\nBatch Effect Tests (Chromosome Proportions):\n")
                for chr_col in available_chrs:
                    batch_groups_chr = [group[chr_col].dropna().values 
                                       for _, group in df.groupby('batch_id')]
                    batch_groups_chr = [g for g in batch_groups_chr if len(g) > 0]
                    if len(batch_groups_chr) >= 2:
                        try:
                            h_stat, p_value = stats.kruskal(*batch_groups_chr)
                            sig_marker = " ***" if p_value < 0.05 else ""
                            f.write(f"  {chr_col.replace('_prop', '')}: H={h_stat:.2f}, p={p_value:.4f}{sig_marker}\n")
                        except ValueError:
                            f.write(f"  {chr_col.replace('_prop', '')}: All values identical\n")
            
            # Result distribution
            f.write("\n5. RESULT DISTRIBUTION (by Month)\n")
            f.write("-" * 80 + "\n")
            for month in sorted(df['month'].unique()):
                month_df = df[df['month'] == month]
                f.write(f"\n{month}:\n")
                for result, count in month_df['Result'].value_counts().items():
                    pct = count / len(month_df) * 100
                    f.write(f"  {result}: {count} ({pct:.1f}%)\n")
        
        else:
            # Multiple labs: include comparison
            f.write("4. QC METRICS COMPARISON (UCL vs Cordlife)\n")
            f.write("-" * 80 + "\n")
            
            metrics = ['mapping_rate(%)', 'duplication_rate(%)', 'GC_content(%)', 'SeqFF']
            
            for metric in metrics:
                if metric not in df.columns:
                    continue
                
                ucl_data = df[df['lab'] == 'UCL'][metric].dropna()
                cordlife_data = df[df['lab'] == 'Cordlife'][metric].dropna()
                
                f.write(f"\n{metric}:\n")
                if len(ucl_data) > 0:
                    f.write(f"  UCL: mean={ucl_data.mean():.2f} ± {ucl_data.std():.2f}, "
                           f"median={ucl_data.median():.2f}, n={len(ucl_data)}\n")
                if len(cordlife_data) > 0:
                    f.write(f"  Cordlife: mean={cordlife_data.mean():.2f} ± {cordlife_data.std():.2f}, "
                           f"median={cordlife_data.median():.2f}, n={len(cordlife_data)}\n")
                
                if len(ucl_data) > 0 and len(cordlife_data) > 0:
                    u_stat, p_value = stats.mannwhitneyu(ucl_data, cordlife_data, alternative='two-sided')
                    f.write(f"  Mann-Whitney U test: U={u_stat:.2f}, p={p_value:.4f}\n")
                    if p_value < 0.05:
                        f.write(f"  *** SIGNIFICANT DIFFERENCE (p < 0.05) ***\n")
        
        # PCA results if available
        if pca_results:
            section_num = 5 if len(unique_labs) == 1 else 4
            f.write(f"\n{section_num+1}. PCA ANALYSIS RESULTS\n")
            f.write("-" * 80 + "\n")
            f.write(f"Samples analyzed: {pca_results['n_samples']}\n")
            f.write(f"PC1 explains: {pca_results['explained_variance'][0]*100:.1f}% of variance\n")
            f.write(f"PC2 explains: {pca_results['explained_variance'][1]*100:.1f}% of variance\n")
            f.write(f"Total explained: {sum(pca_results['explained_variance'])*100:.1f}%\n\n")
            
            f.write("Feature Contributions to PC1 (sorted by importance):\n")
            for _, row in pca_results['feature_importance'].iterrows():
                f.write(f"  {row['Feature']:30s}: {row['PC1']:.4f} (PC2: {row['PC2']:.4f})\n")
            f.write("\n")
        
        # Chromosome-only PCA results if available
        if pca_chr_results:
            section_num = 6 if len(unique_labs) == 1 else 5
            f.write(f"\n{section_num+1}. PCA ANALYSIS - CHROMOSOME PROPORTIONS ONLY\n")
            f.write("-" * 80 + "\n")
            f.write(f"Chromosomes used: Chr13, Chr18, Chr21\n")
            f.write(f"Samples analyzed: {pca_chr_results['n_samples']}\n")
            f.write(f"PC1 explains: {pca_chr_results['explained_variance'][0]*100:.1f}% of variance\n")
            f.write(f"PC2 explains: {pca_chr_results['explained_variance'][1]*100:.1f}% of variance\n")
            f.write(f"Total explained: {sum(pca_chr_results['explained_variance'])*100:.1f}%\n\n")
            
            f.write("Feature Contributions (sorted by PC1 importance):\n")
            for _, row in pca_chr_results['feature_importance'].iterrows():
                f.write(f"  {row['Feature']:30s}: PC1={row['PC1']:.4f}, PC2={row['PC2']:.4f}\n")
            f.write("\n")
        
        # QC-only PCA results if available
        if pca_qc_results:
            section_num = 7 if pca_chr_results else 6 if len(unique_labs) == 1 else 5
            if pca_chr_results and len(unique_labs) == 1:
                section_num = 7
            elif pca_chr_results:
                section_num = 6
            elif len(unique_labs) == 1:
                section_num = 6
            else:
                section_num = 5
            
            f.write(f"\n{section_num+1}. PCA ANALYSIS - QC METRICS ONLY\n")
            f.write("-" * 80 + "\n")
            f.write(f"Metrics used: GC content, Mapping rate, Duplication rate, Mapping quality, Coverage\n")
            f.write(f"Samples analyzed: {pca_qc_results['n_samples']}\n")
            f.write(f"PC1 explains: {pca_qc_results['explained_variance'][0]*100:.1f}% of variance\n")
            f.write(f"PC2 explains: {pca_qc_results['explained_variance'][1]*100:.1f}% of variance\n")
            f.write(f"Total explained: {sum(pca_qc_results['explained_variance'])*100:.1f}%\n\n")
            
            f.write("Feature Contributions (sorted by PC1 importance):\n")
            for _, row in pca_qc_results['feature_importance'].iterrows():
                f.write(f"  {row['Feature']:30s}: PC1={row['PC1']:.4f}, PC2={row['PC2']:.4f}\n")
            f.write("\n")
        
        f.write("\n" + "="*80 + "\n")
    
    print(f"✓ Saved: {report_file}")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Batch Effect Analysis')
    parser.add_argument('--sample-list', required=True,
                       help='Sample list TSV file')
    parser.add_argument('--analysis-dir', default='/home/ken/ken-nipt/analysis',
                       help='Analysis base directory (default: /home/ken/ken-nipt/analysis)')
    parser.add_argument('--output-dir', default='batch_effect_analysis_v2',
                       help='Output directory (default: batch_effect_analysis_v2)')
    parser.add_argument('--lab', choices=['UCL', 'Cordlife', 'all'], default='all',
                       help='Lab to analyze (UCL, Cordlife, or all for comparison)')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("="*80)
    print(f"BATCH EFFECT ANALYSIS v2 - {args.lab.upper()}")
    print("="*80)
    print(f"Sample list: {args.sample_list}")
    print(f"Analysis dir: {args.analysis_dir}")
    print(f"Output dir: {args.output_dir}")
    print(f"Target lab: {args.lab}")
    print()
    
    # Extract batch info
    df = extract_batch_info(args.sample_list, args.analysis_dir)
    
    # Filter by lab if specified
    if args.lab != 'all':
        original_count = len(df)
        df = df[df['lab'] == args.lab].copy()
        print(f"\n✓ Filtered to {args.lab}: {original_count} -> {len(df)} samples")
        
        if len(df) == 0:
            print(f"\n❌ ERROR: No {args.lab} samples found!")
            print(f"   Available labs: {', '.join(df['lab'].unique())}")
            sys.exit(1)
    
    # Extract chromosome-level metrics
    df = extract_chromosome_metrics(args.analysis_dir, df)
    
    # Save enriched data
    enriched_file = os.path.join(args.output_dir, 'sample_list_with_batch.tsv')
    df.to_csv(enriched_file, sep='\t', index=False)
    print(f"\n✓ Saved enriched sample list: {enriched_file}")
    
    # Run analyses
    analyze_batch_structure(df, args.output_dir)
    analyze_qc_by_batch(df, args.output_dir)
    
    # Only compare labs if analyzing all
    if args.lab == 'all' and len(df['lab'].unique()) > 1:
        compare_ucl_vs_cordlife(df, args.output_dir)
    
    analyze_batch_effect_over_time(df, args.output_dir)
    analyze_chromosome_distribution(df, args.output_dir)
    pca_results = analyze_pca_with_batch(df, args.output_dir)
    pca_chr_results = analyze_pca_chromosome_only(df, args.output_dir)
    pca_qc_results = analyze_pca_qc_only(df, args.output_dir)
    generate_batch_summary_report(df, args.output_dir, args.lab, pca_results, pca_chr_results, pca_qc_results)
    
    print("\n" + "="*80)
    print("✅ ANALYSIS COMPLETE!")
    print("="*80)
    print(f"\nResults saved in: {args.output_dir}/")
    print("Generated files:")
    print("  - sample_list_with_batch.tsv (enriched data)")
    print("  - 01_batch_structure.png")
    print("  - 02_qc_metrics_by_batch.png")
    if args.lab == 'all' and len(df['lab'].unique()) > 1:
        print("  - 03_ucl_vs_cordlife.png")
    print("  - 04_batch_effect_over_time.png")
    print("  - 05_chromosome_distribution.png")
    print("  - 06_pca_batch_and_lab.png (all metrics)")
    print("  - 07_pca_chromosome_only.png (chr13/18/21 only)")
    print("  - 08_pca_qc_only.png (QC metrics only)")
    print("  - batch_effect_summary.txt")

if __name__ == '__main__':
    main()
