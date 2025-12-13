#!/usr/bin/env python3
"""
Test Aberration-based ROC Logic

Tests the new logic that uses detected aberration length + z-score for ROC/Performance calculation.
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Load aberration data
data_file = Path("/data/md_validation/roc_results/aberration_data.csv")
print("="*80)
print("Loading aberration data...")
print("="*80)
df = pd.read_csv(data_file)

print(f"Total records: {len(df)}")
print(f"Unique samples: {df['sample_id'].nunique()}")
print(f"Columns: {df.columns.tolist()}")
print()

# Show sample data
print("Sample data:")
print(df.head(10))
print()

# Test parameters
target_disease = '1p36'
mode = 'wc_orig'
ff_value = 10.0
deletion_length_mb = 3.0
zscore_threshold = 3.0
min_detect_length = 1.0

print("="*80)
print(f"Test Case: One-vs-Rest ROC Calculation")
print("="*80)
print(f"Target Disease: {target_disease}")
print(f"Mode: {mode}")
print(f"FF: {ff_value}%")
print(f"Deletion Length: {deletion_length_mb} Mb")
print(f"Z-score Threshold: {zscore_threshold}")
print(f"Min Detect Length: {min_detect_length} Mb")
print()

# Step 1: Create One-vs-Rest dataset
print("-"*80)
print("Step 1: Create One-vs-Rest Dataset")
print("-"*80)

# Add target_disease column for One-vs-Rest
df['target_disease'] = target_disease

# Filter by mode, FF, and deletion_length
df_filtered = df[
    (df['mode'] == mode) &
    (df['ff'] == ff_value) &
    (df['deletion_length_mb'] == deletion_length_mb)
].copy()

print(f"After basic filtering: {len(df_filtered)} records")

# Separate positive and negative samples
df_positive = df_filtered[df_filtered['disease'] == target_disease].copy()
df_negative = df_filtered[df_filtered['disease'] != target_disease].copy()

print(f"Positive samples ({target_disease}): {len(df_positive)}")
print(f"Negative samples (other diseases): {len(df_negative)}")
print()

# Step 2: Apply Min Detect Length filter
print("-"*80)
print("Step 2: Apply Min Detect Length Filter")
print("-"*80)

df_positive_filtered = df_positive[df_positive['detected_mb'] >= min_detect_length].copy()
df_negative_filtered = df_negative[df_negative['detected_mb'] >= min_detect_length].copy()

print(f"Positive samples (after min detect length filter): {len(df_positive_filtered)}")
print(f"  Removed: {len(df_positive) - len(df_positive_filtered)}")
print(f"Negative samples (after min detect length filter): {len(df_negative_filtered)}")
print(f"  Removed: {len(df_negative) - len(df_negative_filtered)}")
print()

# Show examples of removed positive samples
if len(df_positive) > len(df_positive_filtered):
    removed = df_positive[df_positive['detected_mb'] < min_detect_length]
    print(f"Examples of removed positive samples (detected_mb < {min_detect_length}):")
    print(removed[['sample_id', 'zscore', 'detected_mb']].head())
    print()

# Step 3: Calculate TP, FP, TN, FN at given threshold
print("-"*80)
print("Step 3: Calculate Confusion Matrix at Threshold")
print("-"*80)

# True Positives: Positive samples with |zscore| >= threshold AND detected_mb >= min_length
tp = np.sum(df_positive_filtered['zscore'].abs() >= zscore_threshold)

# False Negatives: Positive samples with |zscore| < threshold OR detected_mb < min_length
# Option 1: Only count samples that passed min_length filter but failed threshold
fn_passed_length = np.sum(df_positive_filtered['zscore'].abs() < zscore_threshold)
# Option 2: Also count samples that failed min_length filter
fn_failed_length = len(df_positive) - len(df_positive_filtered)
fn_total = fn_passed_length + fn_failed_length

# False Positives: Negative samples with |zscore| >= threshold AND detected_mb >= min_length
fp = np.sum(df_negative_filtered['zscore'].abs() >= zscore_threshold)

# True Negatives: Negative samples with |zscore| < threshold OR detected_mb < min_length
tn_passed_length = np.sum(df_negative_filtered['zscore'].abs() < zscore_threshold)
tn_failed_length = len(df_negative) - len(df_negative_filtered)
tn_total = tn_passed_length + tn_failed_length

print(f"TP (detected correctly): {tp}")
print(f"FN (missed detection):")
print(f"  - Failed threshold (but passed length): {fn_passed_length}")
print(f"  - Failed min length: {fn_failed_length}")
print(f"  - Total FN: {fn_total}")
print()
print(f"FP (false alarm): {fp}")
print(f"TN (correctly rejected):")
print(f"  - Failed threshold (but passed length): {tn_passed_length}")
print(f"  - Failed min length: {tn_failed_length}")
print(f"  - Total TN: {tn_total}")
print()

# Step 4: Calculate metrics
print("-"*80)
print("Step 4: Calculate Performance Metrics")
print("-"*80)

sensitivity = tp / (tp + fn_total) if (tp + fn_total) > 0 else 0
specificity = tn_total / (tn_total + fp) if (tn_total + fp) > 0 else 0
ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
npv = tn_total / (tn_total + fn_total) if (tn_total + fn_total) > 0 else 0

print(f"Sensitivity (TPR): {sensitivity:.4f} ({sensitivity*100:.2f}%)")
print(f"Specificity (TNR): {specificity:.4f} ({specificity*100:.2f}%)")
print(f"PPV (Precision): {ppv:.4f} ({ppv*100:.2f}%)")
print(f"NPV: {npv:.4f} ({npv*100:.2f}%)")
print()

# Step 5: Show distribution of z-scores and detected lengths
print("-"*80)
print("Step 5: Distribution Analysis")
print("-"*80)

print("Positive samples - Z-score distribution:")
print(f"  Mean: {df_positive_filtered['zscore'].abs().mean():.2f}")
print(f"  Median: {df_positive_filtered['zscore'].abs().median():.2f}")
print(f"  Min: {df_positive_filtered['zscore'].abs().min():.2f}")
print(f"  Max: {df_positive_filtered['zscore'].abs().max():.2f}")
print()

print("Positive samples - Detected length distribution:")
print(f"  Mean: {df_positive_filtered['detected_mb'].mean():.2f} Mb")
print(f"  Median: {df_positive_filtered['detected_mb'].median():.2f} Mb")
print(f"  Min: {df_positive_filtered['detected_mb'].min():.2f} Mb")
print(f"  Max: {df_positive_filtered['detected_mb'].max():.2f} Mb")
print()

if len(df_negative_filtered) > 0:
    print("Negative samples - Z-score distribution:")
    print(f"  Mean: {df_negative_filtered['zscore'].abs().mean():.2f}")
    print(f"  Median: {df_negative_filtered['zscore'].abs().median():.2f}")
    print(f"  Min: {df_negative_filtered['zscore'].abs().min():.2f}")
    print(f"  Max: {df_negative_filtered['zscore'].abs().max():.2f}")
    print()
    
    print("Negative samples - Detected length distribution:")
    print(f"  Mean: {df_negative_filtered['detected_mb'].mean():.2f} Mb")
    print(f"  Median: {df_negative_filtered['detected_mb'].median():.2f} Mb")
    print(f"  Min: {df_negative_filtered['detected_mb'].min():.2f} Mb")
    print(f"  Max: {df_negative_filtered['detected_mb'].max():.2f} Mb")
    print()

# Step 6: Test with different min_detect_length values
print("="*80)
print("Step 6: Test Different Min Detect Length Values")
print("="*80)

for min_len in [0.5, 1.0, 2.0, 3.0]:
    df_pos_len = df_positive[df_positive['detected_mb'] >= min_len].copy()
    df_neg_len = df_negative[df_negative['detected_mb'] >= min_len].copy()
    
    tp_len = np.sum(df_pos_len['zscore'].abs() >= zscore_threshold)
    fn_len = len(df_positive) - tp_len
    fp_len = np.sum(df_neg_len['zscore'].abs() >= zscore_threshold)
    tn_len = len(df_negative) - fp_len
    
    sens_len = tp_len / (tp_len + fn_len) if (tp_len + fn_len) > 0 else 0
    spec_len = tn_len / (tn_len + fp_len) if (tn_len + fp_len) > 0 else 0
    
    print(f"Min Detect Length = {min_len} Mb:")
    print(f"  Positive samples retained: {len(df_pos_len)}/{len(df_positive)} ({len(df_pos_len)/len(df_positive)*100:.1f}%)")
    print(f"  Sensitivity: {sens_len:.4f} ({sens_len*100:.2f}%)")
    print(f"  Specificity: {spec_len:.4f} ({spec_len*100:.2f}%)")
    print()

print("="*80)
print("Logic Verification Complete!")
print("="*80)


