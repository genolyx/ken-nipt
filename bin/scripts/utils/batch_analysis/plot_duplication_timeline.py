#!/usr/bin/env python3
"""
전체 샘플의 QC Metrics Timeline 시각화
Batch별로 그룹화하고 시간순으로 정렬
Duplication Rate와 Coverage 모두 시각화
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys
import os

def plot_metric_timeline(output_dir, lab_name, metric_name, metric_column, ylabel):
    """
    샘플들의 특정 metric을 시간순으로 시각화
    """
    # Load clustered data
    clustered_file = os.path.join(output_dir, 'sample_list_with_clusters.tsv')
    
    if not os.path.exists(clustered_file):
        print(f"❌ File not found: {clustered_file}")
        return
    
    df = pd.read_csv(clustered_file, sep='\t')
    
    # Filter samples with the metric
    df = df[df[metric_column].notna()].copy()
    
    print(f"\n{'='*80}")
    print(f"{metric_name.upper()} TIMELINE - {lab_name}")
    print(f"{'='*80}\n")
    print(f"Total samples: {len(df)}\n")
    
    # Sort by month and batch_id (chronologically)
    df['sort_key'] = df['month'].astype(str) + '_' + df['batch_id'].astype(str)
    df = df.sort_values(['month', 'batch_id']).reset_index(drop=True)
    
    # Add sample index for x-axis
    df['sample_idx'] = range(len(df))
    
    # Get batch boundaries and midpoints for x-axis labels
    batch_info = []
    current_batch = None
    batch_start_idx = 0
    
    for idx, row in df.iterrows():
        batch = row['batch_id']
        if batch != current_batch:
            if current_batch is not None:
                # Save previous batch info
                batch_end_idx = idx - 1
                batch_mid = (batch_start_idx + batch_end_idx) / 2
                batch_size = batch_end_idx - batch_start_idx + 1
                batch_info.append({
                    'batch': current_batch,
                    'start': batch_start_idx,
                    'end': batch_end_idx,
                    'mid': batch_mid,
                    'size': batch_size
                })
            current_batch = batch
            batch_start_idx = idx
    
    # Add last batch
    if current_batch is not None:
        batch_end_idx = len(df) - 1
        batch_mid = (batch_start_idx + batch_end_idx) / 2
        batch_size = batch_end_idx - batch_start_idx + 1
        batch_info.append({
            'batch': current_batch,
            'start': batch_start_idx,
            'end': batch_end_idx,
            'mid': batch_mid,
            'size': batch_size
        })
    
    print(f"Number of batches: {len(batch_info)}\n")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(20, 8))
    
    # Plot duplication rate by cluster
    colors = {0: 'steelblue', 1: 'coral'}
    labels = {0: 'Group 1 (Low cov/dup)', 1: 'Group 2 (High cov/dup)'}
    
    for cluster in [0, 1]:
        cluster_data = df[df['cluster'] == cluster]
        ax.scatter(cluster_data['sample_idx'], 
                  cluster_data[metric_column],
                  c=colors[cluster], 
                  label=labels[cluster],
                  alpha=0.6, 
                  s=20,
                  edgecolors='black',
                  linewidth=0.3)
    
    # Draw vertical lines at batch boundaries
    for i, batch in enumerate(batch_info):
        if i > 0:  # Don't draw at the very beginning
            ax.axvline(x=batch['start'] - 0.5, color='gray', linestyle='--', 
                      linewidth=0.8, alpha=0.5, zorder=0)
    
    # Add batch labels on x-axis
    batch_positions = [b['mid'] for b in batch_info]
    batch_labels = [b['batch'] for b in batch_info]
    batch_sizes = [b['size'] for b in batch_info]
    
    # Set x-axis with batch labels
    ax.set_xticks(batch_positions)
    ax.set_xticklabels(batch_labels, rotation=90, ha='center', fontsize=8)
    
    # Add secondary x-axis showing sample counts per batch
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(batch_positions)
    ax2.set_xticklabels([f'n={s}' for s in batch_sizes], rotation=90, ha='center', fontsize=7, color='gray')
    ax2.tick_params(axis='x', length=0)
    
    # Labels and title
    ax.set_xlabel('Batch ID (chronologically sorted)', fontsize=13, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=13, fontweight='bold')
    ax.set_title(f'{metric_name} Timeline by Batch - {lab_name}\n(All samples sorted chronologically)', 
                fontsize=14, fontweight='bold')
    
    # Grid
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.7)
    ax.set_axisbelow(True)
    
    # Legend
    ax.legend(loc='upper left', fontsize=11, framealpha=0.9)
    
    # Add statistics text box
    metric_data = df[metric_column]
    if 'coverage' in metric_column.lower():
        stats_text = (
            f"Total samples: {len(df)}\n"
            f"Batches: {len(batch_info)}\n"
            f"Mean: {metric_data.mean():.3f}\n"
            f"Median: {metric_data.median():.3f}\n"
            f"Range: {metric_data.min():.3f} - {metric_data.max():.3f}"
        )
    else:
        stats_text = (
            f"Total samples: {len(df)}\n"
            f"Batches: {len(batch_info)}\n"
            f"Mean: {metric_data.mean():.2f}%\n"
            f"Median: {metric_data.median():.2f}%\n"
            f"Range: {metric_data.min():.2f}% - {metric_data.max():.2f}%"
        )
    ax.text(0.98, 0.97, stats_text, transform=ax.transAxes,
           fontsize=10, verticalalignment='top', horizontalalignment='right',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    output_filename = f'{metric_name.lower().replace(" ", "_")}_timeline.png'
    plt.savefig(os.path.join(output_dir, output_filename), dpi=150, bbox_inches='tight')
    print(f"✓ Saved: {output_filename}\n")
    plt.close()
    
    # ========================================
    # Create a more detailed plot focusing on batch trends
    # ========================================
    fig, ax = plt.subplots(figsize=(18, 7))
    
    # Plot with connecting lines within each batch
    for i, batch in enumerate(batch_info):
        batch_data = df[(df['sample_idx'] >= batch['start']) & (df['sample_idx'] <= batch['end'])]
        
        # Use different colors for different clusters
        for cluster in [0, 1]:
            cluster_batch_data = batch_data[batch_data['cluster'] == cluster]
            if len(cluster_batch_data) > 0:
                ax.scatter(cluster_batch_data['sample_idx'], 
                          cluster_batch_data[metric_column],
                          c=colors[cluster],
                          alpha=0.6,
                          s=25,
                          edgecolors='black',
                          linewidth=0.3,
                          zorder=3)
        
        # Add batch mean line
        batch_mean = batch_data[metric_column].mean()
        ax.hlines(batch_mean, batch['start'] - 0.5, batch['end'] + 0.5,
                 colors='red', linestyles='-', linewidth=2, alpha=0.7, zorder=2)
    
    # Batch boundaries
    for i, batch in enumerate(batch_info):
        if i > 0:
            ax.axvline(x=batch['start'] - 0.5, color='gray', linestyle='--', 
                      linewidth=1, alpha=0.6, zorder=1)
    
    # X-axis labels
    ax.set_xticks(batch_positions)
    ax.set_xticklabels(batch_labels, rotation=90, ha='center', fontsize=8)
    
    # Labels
    ax.set_xlabel('Batch ID (chronologically sorted)', fontsize=12, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
    ax.set_title(f'{metric_name} by Batch with Batch Means (Red Lines) - {lab_name}', 
                fontsize=14, fontweight='bold')
    
    # Grid
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.7)
    ax.set_axisbelow(True)
    
    # Legend
    legend_elements = [
        mpatches.Patch(color=colors[0], label='Group 1 (Low cov/dup)', alpha=0.6),
        mpatches.Patch(color=colors[1], label='Group 2 (High cov/dup)', alpha=0.6),
        plt.Line2D([0], [0], color='red', linewidth=2, label='Batch Mean', alpha=0.7)
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10, framealpha=0.9)
    
    plt.tight_layout()
    output_filename_means = f'{metric_name.lower().replace(" ", "_")}_timeline_with_means.png'
    plt.savefig(os.path.join(output_dir, output_filename_means), dpi=150, bbox_inches='tight')
    print(f"✓ Saved: {output_filename_means}\n")
    plt.close()
    
    # ========================================
    # Print batch statistics
    # ========================================
    print(f"{'='*80}")
    print(f"BATCH STATISTICS")
    print(f"{'='*80}\n")
    print(f"{'Batch':<30} {'N':>5} {'Mean':>8} {'Median':>8} {'Std':>8} {'Min':>8} {'Max':>8}")
    print(f"{'-'*80}")
    
    for batch in batch_info:
        batch_data = df[(df['sample_idx'] >= batch['start']) & (df['sample_idx'] <= batch['end'])]
        metric_batch_data = batch_data[metric_column]
        
        if 'coverage' in metric_column.lower():
            print(f"{batch['batch']:<30} {batch['size']:>5} "
                  f"{metric_batch_data.mean():>8.3f} {metric_batch_data.median():>8.3f} "
                  f"{metric_batch_data.std():>8.3f} {metric_batch_data.min():>8.3f} {metric_batch_data.max():>8.3f}")
        else:
            print(f"{batch['batch']:<30} {batch['size']:>5} "
                  f"{metric_batch_data.mean():>8.2f} {metric_batch_data.median():>8.2f} "
                  f"{metric_batch_data.std():>8.2f} {metric_batch_data.min():>8.2f} {metric_batch_data.max():>8.2f}")
    
    print(f"\n{'='*80}\n")

def plot_all_timelines(output_dir, lab_name):
    """
    Plot timelines for both duplication rate and coverage
    """
    print(f"\n{'='*80}")
    print(f"TIMELINE ANALYSIS - {lab_name}")
    print(f"{'='*80}\n")
    
    # Plot duplication rate
    plot_metric_timeline(output_dir, lab_name, 
                        'Duplication Rate', 'duplication_rate(%)', 'Duplication Rate (%)')
    
    # Plot coverage
    plot_metric_timeline(output_dir, lab_name,
                        'Coverage', 'mean_coverageData(X)', 'Coverage (X)')
    
    print("✓ All timeline analyses complete!")
    print(f"\nGenerated files:")
    print(f"  - duplication_rate_timeline.png")
    print(f"  - duplication_rate_timeline_with_means.png")
    print(f"  - coverage_timeline.png")
    print(f"  - coverage_timeline_with_means.png")
    print()

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python3 plot_duplication_timeline.py <output_dir> <lab_name>")
        sys.exit(1)
    
    plot_all_timelines(sys.argv[1], sys.argv[2])
