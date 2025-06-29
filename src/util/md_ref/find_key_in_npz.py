import numpy as np
import argparse
import os

def check_npz_keys(npz_file):
    """Prints the keys available in an npz file."""
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
            print(f"'sample' shape: {sample_data.shape}")

            # If the sample is a scalar (0-dimensional), print directly
            if sample_data.shape == ():
                print(f"'sample' contains a scalar value: {sample_data}")
            elif isinstance(sample_data, np.ndarray):
                if sample_data.size > 0:
                    print(f"'sample' first element:\n{sample_data.flat[0]}")  # Use .flat for safety
                else:
                    print("Warning: 'sample' array is empty.")
            elif isinstance(sample_data, list):
                print(f"'sample' length: {len(sample_data)}")
                if len(sample_data) > 0:
                    print(f"'sample' first 5 elements:\n{sample_data[:5]}")
            elif isinstance(sample_data, dict):
                print(f"'sample' keys:\n{list(sample_data.keys())}")
            else:
                print(f"'sample' value:\n{sample_data}")

        else:
            print("No 'sample' key found in this file.")

    except Exception as e:
        print(f"Error inspecting NPZ file: {e}")

def check_npz_structure(npz_file):
    """Prints the keys and structure of an npz file."""
    try:
        data = np.load(npz_file, allow_pickle=True)
        print(f"Keys in {npz_file}: {list(data.keys())}")

        if "sample" in data:
            sample_data = data["sample"]
            print(f"'sample' data type: {type(sample_data)}")
            print(f"'sample' shape: {sample_data.shape}")

            if sample_data.shape == ():
                print(f"'sample' contains a scalar value: {sample_data}")
            elif isinstance(sample_data, np.ndarray):
                if sample_data.size > 0:
                    print(f"'sample' first element:\n{sample_data.flat[0]}")
                else:
                    print("Warning: 'sample' array is empty.")
            elif isinstance(sample_data, list):
                print(f"'sample' length: {len(sample_data)}")
                if len(sample_data) > 0:
                    print(f"'sample' first 5 elements:\n{sample_data[:5]}")
            elif isinstance(sample_data, dict):
                print(f"'sample' keys:\n{list(sample_data.keys())}")
            else:
                print(f"'sample' contents: {sample_data}")

        else:
            print("No 'sample' key found in this file.")

    except Exception as e:
        print(f"Error checking NPZ structure: {e}")

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

        print("\n--- Checking NPZ Structure ---")
        check_npz_structure(args.npz_file)

