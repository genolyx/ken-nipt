import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import pandas as pd

def compare_matrix_npz(file1, file2, key='median'):
    a = np.load(file1)
    b = np.load(file2)
    mat_a = a[key]
    mat_b = b[key]
    diff = mat_a - mat_b
    abs_diff = np.abs(diff)

    print(f"[{key}] Matrix Comparison:")
    print(f" - {os.path.basename(file1)} vs {os.path.basename(file2)}")
    print(f" → Shape: {mat_a.shape}")
    print(f" → Max diff: {np.max(abs_diff):.6f}")
    print(f" → Mean diff: {np.mean(abs_diff):.6f}")
    print(f" → Non-zero diff count: {(abs_diff > 1e-6).sum()}")

    plt.figure(figsize=(8, 6))
    sns.heatmap(diff, cmap='bwr', center=0, linewidths=0.3, linecolor='gray')
    plt.title(f"Difference in '{key}' matrix")
    plt.xlabel("Columns")
    plt.ylabel("Rows")
    plt.tight_layout()
    out_name = f"diff_{key}_mat_{os.path.basename(file1).replace('.npz','')}__{os.path.basename(file2).replace('.npz','')}.png"
    plt.savefig(out_name)
    print(f"[INFO] Saved heatmap: {out_name}")
    plt.close()

def compare_vector_npz(file1, file2, key='median', chr_ids_key='chr_ids', bin_ids_key='bin_ids'):
    a = np.load(file1)
    b = np.load(file2)
    vec_a = a[key]
    vec_b = b[key]
    chr_ids = a[chr_ids_key]
    bin_ids = a[bin_ids_key]

    diff = vec_a - vec_b
    abs_diff = np.abs(diff)

    print(f"[{key}] Vector Comparison:")
    print(f" - {os.path.basename(file1)} vs {os.path.basename(file2)}")
    print(f" → Length: {len(vec_a)}")
    print(f" → Max diff: {np.max(abs_diff):.6f}")
    print(f" → Mean diff: {np.mean(abs_diff):.6f}")
    print(f" → Non-zero diff count: {(abs_diff > 1e-6).sum()}")

    df = pd.DataFrame({'chr': chr_ids, 'bin': bin_ids, 'diff': diff})
    pivot = df.pivot(index='chr', columns='bin', values='diff')

    plt.figure(figsize=(20, 6))
    sns.heatmap(pivot, cmap='bwr', center=0, linewidths=0.3, linecolor='gray')
    plt.title(f"Difference in '{key}' values (bin-level)")
    plt.xlabel("Genomic Bin")
    plt.ylabel("Chromosome")
    plt.tight_layout()
    out_name = f"diff_{key}_bin_{os.path.basename(file1).replace('.npz','')}__{os.path.basename(file2).replace('.npz','')}.png"
    plt.savefig(out_name)
    print(f"[INFO] Saved heatmap: {out_name}")
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare .npz reference files (matrix or bin vector).")
    parser.add_argument("file1", help="Path to first .npz file")
    parser.add_argument("file2", help="Path to second .npz file")
    parser.add_argument("--key", default="median", help="Key to compare (median or mad)")
    args = parser.parse_args()

    # Detect type: chrom-level matrix vs bin-level vector
    with np.load(args.file1) as f:
        keys = set(f.keys())
        is_matrix = 'chroms' in keys
        is_vector = 'chr_ids' in keys and 'bin_ids' in keys

    if is_matrix:
        compare_matrix_npz(args.file1, args.file2, key=args.key)
    elif is_vector:
        compare_vector_npz(args.file1, args.file2, key=args.key)
    else:
        raise ValueError("Could not determine file type (matrix or bin vector).")
