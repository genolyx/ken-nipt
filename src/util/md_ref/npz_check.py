#!/usr/bin/env python3
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import zscore

def load_reference_matrices(input_dir):
    npz_files = [f for f in os.listdir(input_dir) if f.endswith(".npz")]
    data_matrices = []
    sample_names = []
    
    for npz_file in npz_files:
        filepath = os.path.join(input_dir, npz_file)
        try:
            with np.load(filepath) as reference:
                # Check for common keys
                if "reference_matrix" in reference:
                    data_matrix = reference["reference_matrix"]
                elif "bins" in reference:
                    data_matrix = reference["bins"]
                else:
                    print(f"Skipping {npz_file}: No valid reference matrix found.")
                    continue
                
                # Ensure the matrix is 2D (reshape if 1D)
                if len(data_matrix.shape) == 1:
                    data_matrix = data_matrix.reshape(1, -1)
                
                data_matrices.append(data_matrix)
                sample_names.append(npz_file)
        except Exception as e:
            print(f"Error processing {npz_file}: {e}")
    
    return data_matrices, sample_names

def main():
    parser = argparse.ArgumentParser(
        description="Check NPZ files for a valid reference set, perform PCA, and detect outliers."
    )
    parser.add_argument("input_dir", help="Input directory containing .npz files")
    parser.add_argument("output_png", help="Output PNG file for the PCA scatter plot")
    parser.add_argument("output_outliers", help="Output text file to list outlier sample names")
    args = parser.parse_args()
    
    data_matrices, sample_names = load_reference_matrices(args.input_dir)
    
    if not data_matrices:
        print("No valid reference matrices found. Exiting.")
        return
    
    # Ensure consistency: all matrices must have the same number of columns
    num_cols = data_matrices[0].shape[1]
    consistent_data_matrices = []
    consistent_sample_names = []
    
    for mat, name in zip(data_matrices, sample_names):
        if mat.shape[1] == num_cols:
            consistent_data_matrices.append(mat)
            consistent_sample_names.append(name)
        else:
            print(f"Inconsistent shape in file {name}. Expected {num_cols} columns, got {mat.shape[1]}. Skipping file.")
    
    if not consistent_data_matrices:
        print("No consistent reference matrices found. Exiting.")
        return
    
    # Combine matrices vertically
    data_matrix = np.vstack(consistent_data_matrices)
    print(f"Loaded {len(consistent_sample_names)} samples with shape: {data_matrix.shape}")
    
    # Perform PCA to reduce dimensions to 2 components
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(data_matrix)
    pc1, pc2 = principal_components[:, 0], principal_components[:, 1]
    
    # Save PCA scatter plot
    plt.figure(figsize=(8, 6))
    plt.scatter(pc1, pc2, alpha=0.7)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("PCA of Reference Samples")
    plt.savefig(args.output_png)
    plt.close()
    print(f"PCA plot saved as {args.output_png}")
    
    # Detect outliers using z-score (threshold > 3)
    z_scores = np.abs(zscore(principal_components, axis=0))
    outlier_indices = np.where((z_scores[:, 0] > 3) | (z_scores[:, 1] > 3))[0]
    
    # Save outlier sample names to the specified file
    with open(args.output_outliers, "w") as f:
        for idx in outlier_indices:
            f.write(consistent_sample_names[idx] + "\n")
    print(f"Outliers saved in {args.output_outliers}")
    
    if outlier_indices.size > 0:
        print(f"Outlier samples: {[consistent_sample_names[i] for i in outlier_indices]}")
    else:
        print("No outliers detected.")

if __name__ == "__main__":
    main()
