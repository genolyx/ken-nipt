#!/usr/bin/env python3
import os
import csv
import argparse

try:
    import pandas as pd
except ImportError:
    pd = None

def load_age_mapping(sample_info_path, default_age=30):
    """
    Reads the Excel sheet at sample_info_path, expecting columns
    'Lab ID' and 'Age'. Returns a dict mapping Lab ID (as str) → Age.
    """
    if not sample_info_path:
        return {}
    if pd is None:
        raise RuntimeError("pandas is required to read sample_info Excel files.")
    df = pd.read_excel(sample_info_path, engine='openpyxl')
    if 'Lab ID' not in df.columns or 'Age' not in df.columns:
        raise ValueError("sample_info must have 'Lab ID' and 'Age' columns")
    # Ensure keys are strings
    return { str(lab_id): age for lab_id, age in zip(df['Lab ID'], df['Age']) }

def generate_tsv(root_dir, output_tsv, lab, sample_info_path, default_age=30):
    work_dir = os.path.basename(os.path.normpath(root_dir))
    age_mapping = load_age_mapping(sample_info_path, default_age)

    rows = []
    for subdir in sorted(os.listdir(root_dir)):
        subdir_path = os.path.join(root_dir, subdir)
        if not os.path.isdir(subdir_path):
            continue

        fq1 = fq2 = None
        for fname in os.listdir(subdir_path):
            if 'R1' in fname:
                fq1 = fname
            elif 'R2' in fname:
                fq2 = fname

        if not (fq1 and fq2):
            print(f"[WARNING] Skipping {subdir}: R1 or R2 file not found.", file=sys.stderr)
            continue

        # look up age, default if missing
        age = age_mapping.get(subdir, default_age)

        rows.append([subdir, work_dir, fq1, fq2, age, lab])

    # write out TSV
    with open(output_tsv, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['SAMPLE_NAME', 'WORK_DIR', 'FQ1', 'FQ2', 'AGE', 'LAB'])
        writer.writerows(rows)

    print(f"[INFO] TSV file saved as: {output_tsv}")

def main():
    parser = argparse.ArgumentParser(
        description="Generate a TSV sample sheet from subdirectories of FASTQ data."
    )
    parser.add_argument(
        'root_dir',
        help='Root directory (e.g. /home/ken/…/fastq/250612_01) containing sample subdirs'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output TSV filename'
    )
    parser.add_argument(
        '--lab', '-l',
        default='cordlife',
        help='Lab name to put in the LAB column (default: cordlife)'
    )
    parser.add_argument(
        '--sample_info', '-s',
        help='Path to Excel file (.xlsx) with columns "Lab ID" and "Age" for real ages'
    )
    parser.add_argument(
        '--default_age', '-d',
        type=int,
        default=30,
        help='Default age to use if not found in sample_info (default: 30)'
    )

    args = parser.parse_args()
    generate_tsv(
        root_dir=args.root_dir,
        output_tsv=args.output,
        lab=args.lab,
        sample_info_path=args.sample_info,
        default_age=args.default_age
    )

if __name__ == '__main__':
    import sys
    main()
