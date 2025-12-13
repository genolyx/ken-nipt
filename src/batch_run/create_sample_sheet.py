#!/usr/bin/env python3
"""
Create sample_sheet.tsv from old_cordlife_samples.tsv and old_sample_match_result.tsv

This script merges:
- old_cordlife_samples.tsv: SAMPLE_NAME, AGE
- old_sample_match_result.tsv: Order_id (SAMPLE_NAME), Path (fastq directory)

Output: sample_sheet.tsv with columns: SAMPLE_NAME, FQ1, FQ2, LAB, AGE, WORK_DIR
"""

import pandas as pd
import os
import sys
import glob
from pathlib import Path

def find_fastq_files(path_dir):
    """Find R1 and R2 fastq files in the given directory"""
    path_dir = Path(path_dir)
    
    # Find R1 and R2 files (excluding hidden files starting with ._)
    r1_files = sorted([f for f in path_dir.glob("*_R1_*.fastq.gz") if not f.name.startswith("._")])
    r2_files = sorted([f for f in path_dir.glob("*_R2_*.fastq.gz") if not f.name.startswith("._")])
    
    # Also check for .fq.gz extension
    if not r1_files:
        r1_files = sorted([f for f in path_dir.glob("*_R1_*.fq.gz") if not f.name.startswith("._")])
    if not r2_files:
        r2_files = sorted([f for f in path_dir.glob("*_R2_*.fq.gz") if not f.name.startswith("._")])
    
    if not r1_files or not r2_files:
        return None, None
    
    # Return the first matching pair
    return str(r1_files[0]), str(r2_files[0])

def extract_work_dir(path):
    """Extract WORK_DIR from path (e.g., /data/fastq_backup/2411/OPC241100001 -> 2411)"""
    parts = Path(path).parts
    # Find the part that looks like a work directory (4-digit number)
    for part in reversed(parts):
        if part.isdigit() and len(part) == 4:
            return part
    return None

def main():
    # File paths
    old_samples_file = "/data/fastq_backup/old_cordlife_samples.tsv"
    match_result_file = "/data/fastq_backup/old_sample_match_result.tsv"
    output_file = "sample_sheet.tsv"
    
    # Read files
    print(f"Reading {old_samples_file}...")
    samples_df = pd.read_csv(old_samples_file, sep='\t', header=None, names=['SAMPLE_NAME', 'AGE'], dtype={"SAMPLE_NAME": str})
    print(f"  Found {len(samples_df)} samples")
    
    print(f"Reading {match_result_file}...")
    match_df = pd.read_csv(match_result_file, sep='\t', dtype={"Order_id": str})
    print(f"  Found {len(match_df)} matches")
    
    # Merge dataframes
    merged_df = samples_df.merge(
        match_df[['Order_id', 'Path']],
        left_on='SAMPLE_NAME',
        right_on='Order_id',
        how='inner'
    )
    
    print(f"Merged {len(merged_df)} samples")
    
    # Find fastq files and create output
    results = []
    missing_fastq = []
    
    for idx, row in merged_df.iterrows():
        sample_name = row['SAMPLE_NAME']
        age = row['AGE']
        path = row['Path']
        
        # Find fastq files
        fq1, fq2 = find_fastq_files(path)
        
        if fq1 is None or fq2 is None:
            missing_fastq.append(sample_name)
            print(f"WARNING: Could not find fastq files for {sample_name} in {path}")
            continue
        
        # Extract work directory
        work_dir = extract_work_dir(path)
        if work_dir is None:
            print(f"WARNING: Could not extract work_dir from {path}")
            work_dir = "unknown"
        
        results.append({
            'SAMPLE_NAME': sample_name,
            'FQ1': fq1,
            'FQ2': fq2,
            'LAB': 'cordlife',
            'AGE': age,
            'WORK_DIR': work_dir
        })
    
    # Create output dataframe
    output_df = pd.DataFrame(results)
    
    # Save to file
    output_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nCreated {output_file} with {len(output_df)} samples")
    
    if missing_fastq:
        print(f"\nWARNING: {len(missing_fastq)} samples missing fastq files:")
        for sample in missing_fastq[:10]:  # Show first 10
            print(f"  - {sample}")
        if len(missing_fastq) > 10:
            print(f"  ... and {len(missing_fastq) - 10} more")
    
    print(f"\nSample sheet columns: {', '.join(output_df.columns)}")
    print(f"\nFirst few rows:")
    print(output_df.head())

if __name__ == "__main__":
    main()

