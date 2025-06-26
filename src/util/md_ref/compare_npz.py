import numpy as np
import sys

def load_npz(file_path):
    """Load NPZ file and return its keys and contents."""
    data = np.load(file_path, allow_pickle=True)
    return {key: data[key] for key in data.files}

def show_npz(file1):
    data1 = load_npz(file1)

    keys1 = set(data1.keys())
    print(f"keys1 : {keys1}")

def compare_npz(file1, file2):
    """Compare two NPZ files and print differences in structure and content."""
    data1 = load_npz(file1)
    data2 = load_npz(file2)

    keys1 = set(data1.keys())
    keys2 = set(data2.keys())
    print(f"keys1 : {keys1}")
    print(f"keys2 : {keys2}")

    print(f"\nComparing {file1} and {file2}...\n")

    # Check for missing or extra keys
    if keys1 != keys2:
        print(f"Keys only in {file1}: {keys1 - keys2}")
        print(f"Keys only in {file2}: {keys2 - keys1}")

    # Compare contents
    for key in keys1 & keys2:
        print(f"\nChecking key: {key}")

        value1, value2 = data1[key], data2[key]

        # Compare types
        if value1.dtype != value2.dtype:
            print(f"  - Different data types: {value1.dtype} vs {value2.dtype}")

        # Compare shapes
        if value1.shape != value2.shape:
            print(f"  - Different shapes: {value1.shape} vs {value2.shape}")

        # If both are numeric arrays, compare values
        if np.issubdtype(value1.dtype, np.number) and np.issubdtype(value2.dtype, np.number):
            if not np.allclose(value1, value2, equal_nan=True):
                print("  - Numerical differences detected.")

        # If objects, check their types
        if value1.dtype == 'O' or value2.dtype == 'O':
            print("  - Object arrays detected, comparing structures.")
            obj1, obj2 = value1.item(), value2.item()
            if type(obj1) != type(obj2):
                print(f"    - Different object types: {type(obj1)} vs {type(obj2)}")
            elif isinstance(obj1, dict):
                print(f"    - Dictionary keys comparison: {set(obj1.keys()) ^ set(obj2.keys())}")

if __name__ == "__main__":
    show_npz(sys.argv[1])
    '''
    if len(sys.argv) != 3:
        print("Usage: python compare_npz.py <old_reference.npz> <new_reference.npz>")
    else:
        compare_npz(sys.argv[1], sys.argv[2])
    '''
