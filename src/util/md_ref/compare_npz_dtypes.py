import numpy as np
import sys

def compare_npz_dtypes(file1, file2):
    """Compares data types of corresponding keys in two NPZ files."""
    old_data = np.load(file1, allow_pickle=True)
    new_data = np.load(file2, allow_pickle=True)

    print(f"\nComparing NPZ files:\n - {file1}\n - {file2}\n{'-'*50}")

    for key in old_data.files:
        old_dtype = old_data[key].dtype if key in old_data.files else "Missing"
        new_dtype = new_data[key].dtype if key in new_data.files else "Missing"
        
        if old_dtype != new_dtype:
            print(f"Key: {key} | Old dtype: {old_dtype} | New dtype: {new_dtype} ⚠️ Mismatch")
        else:
            print(f"Key: {key} | Dtypes match: {old_dtype}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare_npz_dtypes.py <old_npz> <new_npz>")
    else:
        compare_npz_dtypes(sys.argv[1], sys.argv[2])

