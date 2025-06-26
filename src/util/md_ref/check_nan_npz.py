import numpy as np
import sys

def check_nan_inf_in_npz(npz_file):
    """Checks for NaN, Inf values, and non-numeric data in an NPZ file."""
    data = np.load(npz_file, allow_pickle=True)
    
    print(f"\n🔍 Checking file: {npz_file}\n{'-'*50}")

    for key in data.files:
        arr = data[key]
        
        if np.issubdtype(arr.dtype, np.number):
            nan_count = np.isnan(arr).sum()
            inf_count = np.isinf(arr).sum()
            print(f"✅ Key: {key} | NaNs: {nan_count} | Infs: {inf_count} | Dtype: {arr.dtype}")
        else:
            print(f"⚠️ Key: {key} | Non-numeric data detected: {arr.dtype}")

    print("\n🔍 Checking segment_w if available...\n" + "-"*50)
    
    # Check if segment_w is stored in the NPZ file
    if "segment_w" in data.files:
        segment_w = data["segment_w"]
        print(f"✅ segment_w Found | NaNs: {np.isnan(segment_w).sum()} | Infs: {np.isinf(segment_w).sum()} | Dtype: {segment_w.dtype}")
    else:
        print("⚠️ segment_w is NOT stored in the NPZ file. It might be computed dynamically during runtime.")

    # Specifically check null_ratios.F and null_ratios.M
    print("\n🔍 Checking null_ratios for Female (F) and Male (M)...\n" + "-"*50)
    for gender in ["F", "M"]:
        key = f"null_ratios.{gender}"
        if key in data.files:
            arr = data[key]
            print(f"✅ {key} | NaNs: {np.isnan(arr).sum()} | Infs: {np.isinf(arr).sum()} | Dtype: {arr.dtype}")
        else:
            print(f"⚠️ {key} not found in NPZ.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python check_nan_npz.py <npz_file>")
    else:
        check_nan_inf_in_npz(sys.argv[1])

