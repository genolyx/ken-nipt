#!/usr/bin/env python3
"""
Filter Incomplete Samples from Sample Sheet

Reads sample sheet and filters out completed samples based on marker files.
Outputs a new sample sheet with only incomplete samples.

Usage:
    python3 filter_incomplete_samples.py \
        -i sample_sheet.tsv \
        -o sample_sheet_incomplete.tsv \
        -r /home/ken/ken-nipt \
        -w md_validation/WBS
"""

import argparse
import sys
from pathlib import Path

def check_sample_completed(root_dir: Path, work_dir: str, sample_id: str) -> bool:
    """Check if sample has completed marker file"""
    sample_dir = root_dir / "analysis" / work_dir / sample_id
    
    # Check for completion marker
    marker_file = sample_dir / f"{sample_id}.md_pipeline_completed.marker"
    
    return marker_file.exists()


def main():
    parser = argparse.ArgumentParser(
        description="Filter incomplete samples from sample sheet"
    )
    parser.add_argument('-i', '--input', type=Path, required=True,
                       help='Input sample sheet TSV')
    parser.add_argument('-o', '--output', type=Path, required=True,
                       help='Output sample sheet TSV (incomplete only)')
    parser.add_argument('-r', '--root_dir', type=Path, required=True,
                       help='Root directory (e.g., /home/ken/ken-nipt)')
    parser.add_argument('-w', '--work_dir', type=str, required=True,
                       help='Work directory (e.g., md_validation/WBS)')
    
    args = parser.parse_args()
    
    if not args.input.exists():
        print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
        return 1
    
    print(f"Reading sample sheet: {args.input}")
    
    # Read input file
    with open(args.input, 'r') as f:
        lines = f.readlines()
    
    if len(lines) == 0:
        print("ERROR: Empty input file", file=sys.stderr)
        return 1
    
    # Separate header and data
    header = lines[0]
    data_lines = lines[1:]
    
    print(f"Total samples in input: {len(data_lines)}")
    
    # Filter incomplete samples
    incomplete_lines = []
    completed_count = 0
    
    for line in data_lines:
        if not line.strip():
            continue
        
        # Parse line to get sample_id
        parts = line.strip().split('\t')
        if len(parts) < 2:
            continue
        
        work_dir_from_line = parts[0]
        sample_id = parts[1]
        
        # Check if completed
        if check_sample_completed(args.root_dir, args.work_dir, sample_id):
            completed_count += 1
            print(f"  [Completed] {sample_id}")
        else:
            incomplete_lines.append(line)
            print(f"  [Incomplete] {sample_id}")
    
    print(f"\nSummary:")
    print(f"  Total: {len(data_lines)}")
    print(f"  Completed: {completed_count}")
    print(f"  Incomplete: {len(incomplete_lines)}")
    
    # Write output
    with open(args.output, 'w') as f:
        f.write(header)
        for line in incomplete_lines:
            f.write(line)
    
    print(f"\nWrote {len(incomplete_lines)} incomplete samples to: {args.output}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())



