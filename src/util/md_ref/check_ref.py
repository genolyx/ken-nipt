import numpy as np
import sys

def check_nan_npz(npz_file):
    """Checks for NaN, Inf values, and non-numeric data in an NPZ file."""
    data = np.load(npz_file, allow_pickle=True)
    
    print(f"\nChecking file: {npz_file}\n{'-'*40}")
    
    for key in data.files:
        arr = data[key]
        if isinstance(arr, np.ndarray) and np.issubdtype(arr.dtype, np.number):
            nan_count = np.isnan(arr).sum()
            inf_count = np.isinf(arr).sum()
            print(f"Key: {key} | NaNs: {nan_count} | Infs: {inf_count} | Dtype: {arr.dtype}")
        else:
            print(f"Key: {key} | Non-numeric data detected: {arr.dtype}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python check_nan_npz.py <npz_file>")
    else:
        check_nan_npz(sys.argv[1])

