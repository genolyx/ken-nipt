#!/usr/bin/env python3
"""
Create BED files for different deletion sizes for all diseases in a BED file.

Usage:
    python3 create_deletion_beds.py --base_bed ~/ken-nipt/data/bed/common/TargetDB_md8.bed --sizes 0.5,1,3,5,7,10 --output_dir bed/temp

Output:
    bed/temp/temp_1p36_0_5Mb.bed
    bed/temp/temp_1p36_1Mb.bed
    bed/temp/temp_2q33_0_5Mb.bed
    ...
"""

import argparse
from pathlib import Path
from typing import List, Tuple


# Disease name mapping to short codes
DISEASE_CODES = {
    '1p36 deletion syndrome': '1p36',
    '2q33.1 deletion syndrome': '2q33',
    'Wolf-Hirschhorn syndrome': 'WHS',
    'Cri Du Chat syndrome': 'CDC',
    'Williams-Beuren syndrome': 'WBS',
    'Jacobsen syndrome': 'Jacobsen',
    'Prader-willi/Angelman syndrome': 'PWS',
    'DiGeorge syndrome': 'DGS'
}


def parse_all_beds(bed_file: Path) -> List[Tuple[str, int, int, str, list]]:
    """Parse BED file and return all entries.
    
    Returns:
        List of (chrom, start, end, disease_name, extra_fields)
    """
    # If path is relative and doesn't exist, try looking in bed/ directory
    if not bed_file.is_absolute() and not bed_file.exists():
        script_dir = Path(__file__).parent
        bed_dir_path = script_dir / "bed" / bed_file
        if bed_dir_path.exists():
            bed_file = bed_dir_path
    
    entries = []
    with bed_file.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 4:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                disease = parts[3]
                extra_fields = parts[4:] if len(parts) > 4 else []
                entries.append((chrom, start, end, disease, extra_fields))
    
    if not entries:
        raise ValueError(f"No valid BED entries found in {bed_file}")
    
    return entries


def create_bed_file(output_path: Path, chrom: str, start: int, end: int, disease: str, extra_fields: list = None):
    """Create a BED file with the specified coordinates."""
    with output_path.open('w') as f:
        fields = [chrom, str(start), str(end), disease]
        if extra_fields:
            fields.extend(extra_fields)
        f.write('\t'.join(fields) + '\n')


def format_size_for_filename(size_mb: float) -> str:
    """Format size for filename, replacing '.' with '_' for decimals.
    
    Examples:
        0.5 -> "0_5Mb"
        1.0 -> "1Mb"
        3.5 -> "3_5Mb"
    """
    if size_mb == int(size_mb):
        # Whole number
        return f"{int(size_mb)}Mb"
    else:
        # Has decimal, replace '.' with '_'
        return f"{size_mb}Mb".replace('.', '_')


def main():
    parser = argparse.ArgumentParser(
        description="Create BED files for different deletion sizes for all diseases"
    )
    parser.add_argument(
        '--base_bed',
        required=True,
        type=Path,
        help='Base BED file containing all diseases (e.g., TargetDB_md8.bed)'
    )
    parser.add_argument(
        '--sizes',
        required=True,
        help='Comma-separated list of sizes in Mb (e.g., 0.5,1,3,5,7,10)'
    )
    parser.add_argument(
        '--output_prefix',
        default='temp',
        help='Output prefix for BED files (default: temp)'
    )
    parser.add_argument(
        '--output_dir',
        type=Path,
        default=None,
        help='Output directory (default: bed/ in script directory)'
    )
    parser.add_argument(
        '--mode',
        choices=['start', 'center'],
        default='start',
        help='Alignment mode: start (keep start position) or center (center alignment)'
    )
    
    args = parser.parse_args()
    
    # Parse all BED entries
    entries = parse_all_beds(args.base_bed)
    
    print(f"Found {len(entries)} disease(s) in {args.base_bed}")
    print(f"Mode: {args.mode}")
    print()
    
    # Parse sizes
    sizes_mb = [float(s.strip()) for s in args.sizes.split(',')]
    sizes_mb = sorted(sizes_mb)
    
    print(f"Sizes to generate: {', '.join([f'{s}Mb' for s in sizes_mb])}")
    print()
    
    # Determine output directory
    script_dir = Path(__file__).parent
    if args.output_dir:
        output_dir = args.output_dir
        if not output_dir.is_absolute():
            output_dir = script_dir / output_dir
    else:
        output_dir = script_dir / "bed"
    
    # Create output directory if needed
    output_dir.mkdir(parents=True, exist_ok=True)
    
    created_files = []
    
    # Process each disease
    for chrom, base_start, base_end, disease, extra_fields in entries:
        base_size = base_end - base_start
        
        # Get disease code
        disease_code = DISEASE_CODES.get(disease, disease.replace(' ', '_'))
        
        print(f"Processing: {disease} ({disease_code})")
        print(f"  Base region: {chrom}:{base_start}-{base_end} (size: {base_size/1_000_000:.2f} Mb)")
        
        # Create BED files for each size
        for size_mb in sizes_mb:
            size_bp = int(size_mb * 1_000_000)
            
            if args.mode == 'start':
                # Keep start position, adjust end
                new_start = base_start
                new_end = base_start + size_bp
            else:
                # Center alignment
                center = (base_start + base_end) // 2
                new_start = center - size_bp // 2
                new_end = center + size_bp // 2
                # Ensure start >= 1
                if new_start < 1:
                    new_start = 1
                    new_end = size_bp
            
            # Create filename with proper formatting
            size_str = format_size_for_filename(size_mb)
            output_file = output_dir / f"{args.output_prefix}_{disease_code}_{size_str}.bed"
            
            # Use original extra fields if available, otherwise default
            fields_to_use = extra_fields if extra_fields else ['loss', 'overlap', '-', '-', '-', '-']
            create_bed_file(output_file, chrom, new_start, new_end, disease, fields_to_use)
            
            actual_size = (new_end - new_start) / 1_000_000
            print(f"  Created: {output_file.name} -> {chrom}:{new_start}-{new_end} ({actual_size:.2f} Mb)")
            created_files.append(output_file)
        
        print()
    
    print(f"Total: Created {len(created_files)} BED files in {output_dir}")
    
    return 0


if __name__ == '__main__':
    exit(main())

