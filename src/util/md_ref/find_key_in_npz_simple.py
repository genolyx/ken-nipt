import numpy as np
import argparse
import os

def check_npz_keys(npz_file):
    """Prints the keys available in an NPZ file."""
    try:
        data = np.load(npz_file, allow_pickle=True)
        print(f"Keys in {npz_file}: {list(data.keys())}")
    except Exception as e:
        print(f"Error loading NPZ file: {e}")

def inspect_npz_sample(npz_file):
    """Extracts and prints the sample key to check its real structure."""
    try:
        data = np.load(npz_file, allow_pickle=True)
        print(f"Keys in {npz_file}: {list(data.keys())}")

        if "sample" in data:
            sample_data = data["sample"]
            print(f"'sample' raw type: {type(sample_data)}")

            # If 'sample' is a dictionary
            if isinstance(sample_data, dict):
                print(f"'sample' contains data for {len(sample_data)} chromosomes.")

                # Check for chromosome X and Y representation
                chrom_keys = sorted(sample_data.keys(), key=lambda x: int(x))
                print(f"Chromosomes present: {chrom_keys}")

                if "23" in sample_data:
                    print(f"Chromosome X detected as key '23'.")
                if "24" in sample_data:
                    print(f"Chromosome Y detected as key '24'.")

                # Print first few values for Chromosome X and Y (if present)
                for chrom in ["23", "24"]:
                    if chrom in sample_data:
                        print(f"\nSample values for Chromosome {chrom} (first 10 bins):")
                        print(sample_data[chrom][:10])

            else:
                print(f"'sample' is not a dictionary. It contains: {sample_data}")

        else:
            print("No 'sample' key found in this file.")

    except Exception as e:
        print(f"Error inspecting NPZ file: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Inspect the structure of an NPZ file.")
    parser.add_argument("npz_file", type=str, help="Path to the NPZ file")
    args = parser.parse_args()

    if not os.path.exists(args.npz_file):
        print(f"Error: File '{args.npz_file}' not found.")
    else:
        print("\n--- Checking NPZ Keys ---")
        check_npz_keys(args.npz_file)

        print("\n--- Inspecting Sample Key ---")
        inspect_npz_sample(args.npz_file)

