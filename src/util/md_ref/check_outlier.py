#!/usr/bin/env python3
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import zscore

def load_samples(input_dir):
    """
    Loads all npz files from input_dir that contain a 'sample' key.
    Each file is expected to have a 'sample' object (a dictionary) with keys '1' to '24'.
    Returns a dictionary mapping sample filename to its sample dictionary.
    """
    sample_data_dict = {}
    npz_files = [f for f in os.listdir(input_dir) if f.endswith(".npz")]
    
    for filename in npz_files:
        filepath = os.path.join(input_dir, filename)
        try:
            data = np.load(filepath, allow_pickle=True)
            if 'sample' in data:
                # Convert the object array (holding the dictionary) to a dictionary.
                sample_dict = data['sample'].item()
                sample_data_dict[filename] = sample_dict
            else:
                print(f"Skipping {filename}: 'sample' key not found.")
        except Exception as e:
            print(f"Error processing {filename}: {e}")
    return sample_data_dict

def perform_pca_and_outlier_detection(data_by_key, sample_names, key, output_png, output_outliers):
    """
    Given a list of arrays (data_by_key) and corresponding sample names for a particular key,
    stacks the arrays, performs PCA, saves the scatter plot, detects outliers based on a z-score threshold,
    writes outlier sample names to output_outliers, and returns the list of outlier sample names.
    """
    # Stack arrays into a matrix: shape = (num_samples, length_of_array)
    data_matrix = np.vstack(data_by_key)
    
    # Perform PCA (reduce to 2 components)
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(data_matrix)
    pc1 = principal_components[:, 0]
    pc2 = principal_components[:, 1]
    
    # Plot the PCA results
    plt.figure(figsize=(8, 6))
    plt.scatter(pc1, pc2, alpha=0.7)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title(f"PCA for key {key}")
    plt.savefig(output_png)
    plt.close()
    print(f"PCA plot saved: {output_png}")
    
    # Detect outliers using z-score with a threshold of 3
    z_scores = np.abs(zscore(principal_components, axis=0))
    outlier_indices = np.where((z_scores[:, 0] > 3) | (z_scores[:, 1] > 3))[0]
    outlier_sample_names = [sample_names[i] for i in outlier_indices]
    
    # Write the outlier sample names to file
    with open(output_outliers, "w") as f:
        for name in outlier_sample_names:
            f.write(name + "\n")
    print(f"Outliers saved: {output_outliers}")
    
    return outlier_sample_names

def main():
    parser = argparse.ArgumentParser(
        description="Perform PCA outlier detection per key in npz files and aggregate outliers."
    )
    parser.add_argument("input_dir", help="Input directory containing npz files")
    parser.add_argument("output_dir", help="Output directory for png and outlier files")
    parser.add_argument("--combined_outliers", default="combined_outliers.txt",
                        help="Filename for combined outliers across all keys (default: combined_outliers.txt)")
    args = parser.parse_args()
    
    # Create output directory if it does not exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load samples from npz files
    sample_data_dict = load_samples(args.input_dir)
    if not sample_data_dict:
        print("No valid npz files with 'sample' key found. Exiting.")
        return
    
    # Assume all samples have the same set of keys. Sort keys numerically.
    first_sample = next(iter(sample_data_dict.values()))
    keys = sorted(first_sample.keys(), key=lambda x: int(x))
    
    combined_outlier_set = set()
    
    # Process each key separately
    for key in keys:
        data_by_key = []
        sample_names = []
        for sample_name, sample_dict in sample_data_dict.items():
            if key in sample_dict:
                data_by_key.append(sample_dict[key])
                sample_names.append(sample_name)
            else:
                print(f"Sample {sample_name} does not have key {key}. Skipping it for this key.")
        
        if not data_by_key:
            print(f"No data found for key {key}. Skipping.")
            continue
        
        # Define output filenames for the current key
        output_png = os.path.join(args.output_dir, f"pca_key_{key}.png")
        output_outliers = os.path.join(args.output_dir, f"outlier_key_{key}.txt")
        
        # Perform PCA and outlier detection for the current key
        outlier_samples = perform_pca_and_outlier_detection(
            data_by_key, sample_names, key, output_png, output_outliers
        )
        combined_outlier_set.update(outlier_samples)
    
    # Write combined outlier sample names (present in at least one key) to file
    combined_outliers_file = os.path.join(args.output_dir, args.combined_outliers)
    with open(combined_outliers_file, "w") as f:
        for name in sorted(combined_outlier_set):
            f.write(name + "\n")
    print(f"Combined outliers saved: {combined_outliers_file}")

if __name__ == "__main__":
    main()

