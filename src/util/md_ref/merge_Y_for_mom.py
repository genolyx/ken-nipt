import os
import numpy as np

# Define directories
mom_dir = "/Work/NIPT/refs/cordlife_new/npz/WCX/mom"
orig_dir = "/Work/NIPT/refs/cordlife_new/npz/WCX/orig"
output_dirs = {
    "F": "/Work/NIPT/refs/cordlife_new/npz/WCX/mom/F_Y",
    "M": "/Work/NIPT/refs/cordlife_new/npz/WCX/mom/M_Y"
}

# Ensure output directories exist
for out_dir in output_dirs.values():
    os.makedirs(out_dir, exist_ok=True)

def get_common_prefix(filename):
    """Extracts the common sample ID before '.proper_paired'"""
    return filename.split(".proper_paired")[0]

def merge_y_chromosome(mom_file, orig_file, output_file):
    """Replace Y chromosome data in mom_file with data from orig_file and save."""
    mom_data = np.load(mom_file, allow_pickle=True)
    orig_data = np.load(orig_file, allow_pickle=True)

    # Check if Y chromosome exists in original data
    if "Y" in orig_data:
        y_chr_data = orig_data["Y"]
    else:
        print(f"Warning: No Y chromosome data found in {orig_file}. Skipping.")
        return

    # Create merged dataset
    merged_data = {key: mom_data[key] for key in mom_data.files}  # Copy mom data
    merged_data["Y"] = y_chr_data  # Inject Y-chromosome data

    # Save new NPZ file
    np.savez(output_file, **merged_data)
    print(f"Saved merged file: {output_file}")

# Process Female and Male directories
for gender in ["F", "M"]:
    mom_path = os.path.join(mom_dir, gender)
    orig_path = os.path.join(orig_dir, gender)
    output_path = output_dirs[gender]

    if not os.path.exists(mom_path) or not os.path.exists(orig_path):
        print(f"Skipping {gender}, missing directories.")
        continue

    # Create a mapping of orig files based on their prefix
    orig_files = {get_common_prefix(f): os.path.join(orig_path, f) for f in os.listdir(orig_path) if f.endswith(".npz")}

    # Iterate over mom NPZ files
    for file in os.listdir(mom_path):
        if file.endswith(".npz"):
            mom_file = os.path.join(mom_path, file)
            common_prefix = get_common_prefix(file)

            if common_prefix in orig_files:
                orig_file = orig_files[common_prefix]
                output_file = os.path.join(output_path, file)

                merge_y_chromosome(mom_file, orig_file, output_file)
            else:
                print(f"Warning: No matching original file for {file}. Skipping.")

