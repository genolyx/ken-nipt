import numpy as np
import sys

def check_binsize(npz_file):
    """Checks the binsize used in a WisecondorX NPZ file."""
    data = np.load(npz_file, allow_pickle=True)

    print(f"\n🔍 Checking binsize for: {npz_file}\n" + "-"*50)

    # 1️⃣ Directly check if 'binsize' key exists
    if "binsize" in data.files:
        print(f"✅ Binsize found: {data['binsize'].item()} bp")
        return
    
    print("⚠️ 'binsize' key not found. Attempting to infer binsize...")

    # 2️⃣ Check bin-related keys
    genome_size = 3_100_000_000  # Approximate human genome size in base pairs
    bin_keys = ["bins_per_chr", "masked_bins_per_chr"]
    
    for key in bin_keys:
        if key in data.files:
            total_bins = np.sum(data[key])
            if total_bins > 0:
                inferred_binsize = genome_size // total_bins
                print(f"🔍 Estimated binsize based on {key}: {inferred_binsize} bp")
                return
    
    print("❌ Unable to determine binsize from the NPZ file.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python check_binsize.py <npz_file>")
    else:
        check_binsize(sys.argv[1])

