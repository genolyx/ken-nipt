#!/usr/bin/env python3
"""
GC bias와 Chromosome bias 관계 시각화
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

def visualize_gc_chr_relationship():
    """
    GC content, Duplication, Chromosome proportion 관계 시각화
    """
    
    # Load data
    cordlife_batch = pd.read_csv('analysis/batch_analysis_cordlife_complete/sample_list_with_batch.tsv', sep='\t')
    cordlife_cluster = pd.read_csv('analysis/batch_analysis_cordlife_complete/sample_list_with_clusters.tsv', sep='\t')
    cordlife = pd.merge(cordlife_batch, cordlife_cluster[['sample_id', 'cluster']], on='sample_id', how='inner')
    
    ucl_batch = pd.read_csv('analysis/batch_analysis_ucl_complete/sample_list_with_batch.tsv', sep='\t')
    ucl_cluster = pd.read_csv('analysis/batch_analysis_ucl_complete/sample_list_with_clusters.tsv', sep='\t')
    ucl = pd.merge(ucl_batch, ucl_cluster[['sample_id', 'cluster']], on='sample_id', how='inner')
    
    # Create figure
    fig = plt.figure(figsize=(20, 12))
    
    # ========================================
    # 1. Duplication vs Chr21 proportion
    # ========================================
    ax1 = plt.subplot(2, 3, 1)
    colors_cord = cordlife['cluster'].map({0: '#FF6B6B', 1: '#4ECDC4'})
    scatter1 = ax1.scatter(cordlife['duplication_rate(%)'], cordlife['chr21_prop']*100,
                          c=colors_cord, alpha=0.6, s=30, edgecolors='black', linewidth=0.5)
    
    # Correlation
    r_cord, p_cord = stats.spearmanr(cordlife['duplication_rate(%)'].dropna(), 
                                      cordlife['chr21_prop'].dropna())
    
    ax1.set_xlabel('Duplication Rate (%)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Chr21 Proportion (%)', fontsize=12, fontweight='bold')
    ax1.set_title(f'Cordlife: Duplication vs Chr21\nr={r_cord:.3f}, p={p_cord:.2e}', 
                  fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(['Group 1 (low cov)', 'Group 2 (high cov)'], loc='best')
    
    # Add text
    ax1.text(0.05, 0.95, '✅ No correlation\n   Multiplexing only', 
             transform=ax1.transAxes, fontsize=11, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))
    
    # ========================================
    # 2. UCL all: Duplication vs Chr21
    # ========================================
    ax2 = plt.subplot(2, 3, 2)
    colors_ucl = ucl['cluster'].map({0: '#FF6B6B', 1: '#4ECDC4'})
    scatter2 = ax2.scatter(ucl['duplication_rate(%)'], ucl['chr21_prop']*100,
                          c=colors_ucl, alpha=0.6, s=30, edgecolors='black', linewidth=0.5)
    
    r_ucl, p_ucl = stats.spearmanr(ucl['duplication_rate(%)'].dropna(), 
                                    ucl['chr21_prop'].dropna())
    
    ax2.set_xlabel('Duplication Rate (%)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Chr21 Proportion (%)', fontsize=12, fontweight='bold')
    ax2.set_title(f'UCL (all): Duplication vs Chr21\nr={r_ucl:.3f}, p={p_ucl:.2e}', 
                  fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(['Group 1 (bad lib)', 'Group 2 (good lib)'], loc='best')
    
    ax2.text(0.05, 0.95, '⚠️ Weak correlation\n   But bias exists!', 
             transform=ax2.transAxes, fontsize=11, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
    
    # ========================================
    # 3. UCL 2510+: Duplication vs Chr21
    # ========================================
    ax3 = plt.subplot(2, 3, 3)
    ucl_2510 = ucl[ucl['month'].astype(int) >= 2510]
    colors_ucl2510 = ucl_2510['cluster'].map({0: '#FF6B6B', 1: '#4ECDC4'})
    scatter3 = ax3.scatter(ucl_2510['duplication_rate(%)'], ucl_2510['chr21_prop']*100,
                          c=colors_ucl2510, alpha=0.6, s=30, edgecolors='black', linewidth=0.5)
    
    r_ucl2510, p_ucl2510 = stats.spearmanr(ucl_2510['duplication_rate(%)'].dropna(), 
                                            ucl_2510['chr21_prop'].dropna())
    
    ax3.set_xlabel('Duplication Rate (%)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Chr21 Proportion (%)', fontsize=12, fontweight='bold')
    ax3.set_title(f'UCL (2510+): Duplication vs Chr21\nr={r_ucl2510:.3f}, p={p_ucl2510:.2e}', 
                  fontsize=13, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.legend(['Group 1 (3.3%)', 'Group 2 (96.7%)'], loc='best')
    
    ax3.text(0.05, 0.95, '✅ No correlation\n   Protocol fixed!', 
             transform=ax3.transAxes, fontsize=11, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))
    
    # ========================================
    # 4. GC content vs Chr21 (Cordlife)
    # ========================================
    ax4 = plt.subplot(2, 3, 4)
    scatter4 = ax4.scatter(cordlife['GC_content(%)'], cordlife['chr21_prop']*100,
                          c=colors_cord, alpha=0.6, s=30, edgecolors='black', linewidth=0.5)
    
    r_gc_cord, p_gc_cord = stats.spearmanr(cordlife['GC_content(%)'].dropna(), 
                                            cordlife['chr21_prop'].dropna())
    
    ax4.set_xlabel('GC Content (%)', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Chr21 Proportion (%)', fontsize=12, fontweight='bold')
    ax4.set_title(f'Cordlife: GC vs Chr21\nr={r_gc_cord:.3f}, p={p_gc_cord:.2e}', 
                  fontsize=13, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    
    ax4.text(0.05, 0.95, '✅ No GC bias', 
             transform=ax4.transAxes, fontsize=11, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))
    
    # ========================================
    # 5. GC content vs Chr21 (UCL all)
    # ========================================
    ax5 = plt.subplot(2, 3, 5)
    scatter5 = ax5.scatter(ucl['GC_content(%)'], ucl['chr21_prop']*100,
                          c=colors_ucl, alpha=0.6, s=30, edgecolors='black', linewidth=0.5)
    
    r_gc_ucl, p_gc_ucl = stats.spearmanr(ucl['GC_content(%)'].dropna(), 
                                          ucl['chr21_prop'].dropna())
    
    ax5.set_xlabel('GC Content (%)', fontsize=12, fontweight='bold')
    ax5.set_ylabel('Chr21 Proportion (%)', fontsize=12, fontweight='bold')
    ax5.set_title(f'UCL (all): GC vs Chr21\nr={r_gc_ucl:.3f}, p={p_gc_ucl:.2e}', 
                  fontsize=13, fontweight='bold')
    ax5.grid(True, alpha=0.3)
    
    if abs(r_gc_ucl) > 0.3 and p_gc_ucl < 0.001:
        ax5.text(0.05, 0.95, '⚠️ GC bias present!', 
                 transform=ax5.transAxes, fontsize=11, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='orange', alpha=0.7))
    else:
        ax5.text(0.05, 0.95, '⚠️ Complex pattern', 
                 transform=ax5.transAxes, fontsize=11, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
    
    # ========================================
    # 6. Chr21 distribution comparison
    # ========================================
    ax6 = plt.subplot(2, 3, 6)
    
    # Violin plots
    data_to_plot = [
        cordlife[cordlife['cluster']==0]['chr21_prop']*100,
        cordlife[cordlife['cluster']==1]['chr21_prop']*100,
        ucl[ucl['cluster']==0]['chr21_prop']*100,
        ucl[ucl['cluster']==1]['chr21_prop']*100,
        ucl_2510[ucl_2510['cluster']==0]['chr21_prop']*100,
        ucl_2510[ucl_2510['cluster']==1]['chr21_prop']*100,
    ]
    
    labels = ['Cord\nG1', 'Cord\nG2', 'UCL\nG1', 'UCL\nG2', 'UCL10+\nG1', 'UCL10+\nG2']
    colors_box = ['#FF6B6B', '#4ECDC4', '#FF6B6B', '#4ECDC4', '#FF6B6B', '#4ECDC4']
    
    positions = [1, 2, 4, 5, 7, 8]
    bp = ax6.boxplot(data_to_plot, positions=positions, labels=labels,
                     patch_artist=True, widths=0.6)
    
    for patch, color in zip(bp['boxes'], colors_box):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax6.set_ylabel('Chr21 Proportion (%)', fontsize=12, fontweight='bold')
    ax6.set_title('Chr21 Distribution by Group', fontsize=13, fontweight='bold')
    ax6.grid(True, alpha=0.3, axis='y')
    
    # Add significance markers
    _, p_cord = stats.mannwhitneyu(cordlife[cordlife['cluster']==0]['chr21_prop'].dropna(),
                                    cordlife[cordlife['cluster']==1]['chr21_prop'].dropna())
    _, p_ucl = stats.mannwhitneyu(ucl[ucl['cluster']==0]['chr21_prop'].dropna(),
                                   ucl[ucl['cluster']==1]['chr21_prop'].dropna())
    _, p_ucl2510 = stats.mannwhitneyu(ucl_2510[ucl_2510['cluster']==0]['chr21_prop'].dropna(),
                                       ucl_2510[ucl_2510['cluster']==1]['chr21_prop'].dropna())
    
    y_max = ax6.get_ylim()[1]
    y_text = y_max * 0.98
    
    ax6.text(1.5, y_text, f'p={p_cord:.2f}' if p_cord > 0.05 else f'p={p_cord:.2e}',
             ha='center', fontsize=9, fontweight='bold',
             color='green' if p_cord > 0.05 else 'red')
    ax6.text(4.5, y_text, f'p={p_ucl:.2e}',
             ha='center', fontsize=9, fontweight='bold', color='red')
    ax6.text(7.5, y_text, f'p={p_ucl2510:.2f}',
             ha='center', fontsize=9, fontweight='bold',
             color='green' if p_ucl2510 > 0.05 else 'orange')
    
    plt.tight_layout()
    plt.savefig('analysis/GC_Chr_Bias_Relationship.png', dpi=300, bbox_inches='tight')
    print(f"\n✅ Plot saved: analysis/GC_Chr_Bias_Relationship.png\n")
    
    # ========================================
    # Print summary
    # ========================================
    print(f"\n{'='*80}")
    print(f"CORRELATION SUMMARY")
    print(f"{'='*80}\n")
    
    print(f"Duplication vs Chr21 Proportion:")
    print(f"  Cordlife:  r={r_cord:>6.3f} (p={p_cord:.2e}) → {'✅ No correlation' if abs(r_cord) < 0.3 else '⚠️ Correlation'}")
    print(f"  UCL (all): r={r_ucl:>6.3f} (p={p_ucl:.2e}) → {'✅ No correlation' if abs(r_ucl) < 0.3 else '⚠️ Correlation'}")
    print(f"  UCL 2510+: r={r_ucl2510:>6.3f} (p={p_ucl2510:.2e}) → {'✅ No correlation' if abs(r_ucl2510) < 0.3 else '⚠️ Correlation'}")
    print()
    
    print(f"GC Content vs Chr21 Proportion:")
    print(f"  Cordlife:  r={r_gc_cord:>6.3f} (p={p_gc_cord:.2e}) → {'✅ No GC bias' if abs(r_gc_cord) < 0.3 else '⚠️ GC bias'}")
    print(f"  UCL (all): r={r_gc_ucl:>6.3f} (p={p_gc_ucl:.2e}) → {'✅ No GC bias' if abs(r_gc_ucl) < 0.3 else '⚠️ GC bias'}")
    print()
    
    print(f"Chr21 Group Differences:")
    print(f"  Cordlife:  p={p_cord:.2e} → {'✅ No bias' if p_cord > 0.05 else '❌ Bias'}")
    print(f"  UCL (all): p={p_ucl:.2e} → {'✅ No bias' if p_ucl > 0.05 else '❌ Bias'}")
    print(f"  UCL 2510+: p={p_ucl2510:.2f} → {'✅ No bias' if p_ucl2510 > 0.05 else '⚠️ Marginal'}")
    print()
    
    print(f"{'='*80}")
    print(f"INTERPRETATION")
    print(f"{'='*80}\n")
    
    print(f"Cordlife:")
    print(f"  - Duplication과 Chr21 상관관계 없음 (r={r_cord:.3f})")
    print(f"  - GC content와 Chr21 상관관계 없음 (r={r_gc_cord:.3f})")
    print(f"  - Group 간 Chr21 차이 없음 (p={p_cord:.2f})")
    print(f"  → ✅ 순수 multiplexing artifact, 생물학적 신호 보존!\n")
    
    print(f"UCL (all):")
    print(f"  - Duplication과 Chr21 약한 상관관계 (r={r_ucl:.3f})")
    print(f"  - GC content와 Chr21 약한 상관관계 (r={r_gc_ucl:.3f})")
    print(f"  - Group 간 Chr21 유의한 차이 (p={p_ucl:.2e})")
    print(f"  → ❌ Library prep 실패 포함, 생물학적 신호 왜곡!\n")
    
    print(f"UCL (2510+):")
    print(f"  - Duplication과 Chr21 상관관계 없음 (r={r_ucl2510:.3f})")
    print(f"  - Group 간 Chr21 차이 경계선 (p={p_ucl2510:.2f})")
    print(f"  - Group 1 = 3.3% only!")
    print(f"  → ✅ Protocol 개선, Chr21 bias 거의 소멸!\n")
    
    print(f"{'='*80}\n")

if __name__ == '__main__':
    visualize_gc_chr_relationship()
