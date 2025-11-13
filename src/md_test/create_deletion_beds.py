#!/usr/bin/env python3
"""
Create BED files for different deletion sizes based on 1p36 region.

Usage:
    python3 create_deletion_beds.py --base_bed test_1p36_only.bed --sizes 1,3,5,7,10,15,20 --output_prefix test_1p36

Output:
    test_1p36_1Mb.bed
    test_1p36_3Mb.bed
    ...
"""

import argparse
from pathlib import Path


def parse_bed(bed_file: Path):
    """Parse BED file and return chromosome, start, end, disease name."""
    # If path is relative and doesn't exist, try looking in bed/ directory
    if not bed_file.is_absolute() and not bed_file.exists():
        script_dir = Path(__file__).parent
        bed_dir_path = script_dir / "bed" / bed_file
        if bed_dir_path.exists():
            bed_file = bed_dir_path
    
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
                disease = parts[3] if len(parts) > 3 else '1p36 deletion syndrome'
                return chrom, start, end, disease
    raise ValueError(f"No valid BED entry found in {bed_file}")


def create_bed_file(output_path: Path, chrom: str, start: int, end: int, disease: str, extra_fields: list = None):
    """Create a BED file with the specified coordinates."""
    with output_path.open('w') as f:
        fields = [chrom, str(start), str(end), disease]
        if extra_fields:
            fields.extend(extra_fields)
        f.write('\t'.join(fields) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description="Create BED files for different deletion sizes"
    )
    parser.add_argument(
        '--base_bed',
        required=True,
        type=Path,
        help='Base BED file (e.g., test_1p36_only.bed)'
    )
    parser.add_argument(
        '--sizes',
        required=True,
        help='Comma-separated list of sizes in Mb (e.g., 1,3,5,7,10,15,20)'
    )
    parser.add_argument(
        '--output_prefix',
        default='test_1p36',
        help='Output prefix for BED files (default: test_1p36)'
    )
    parser.add_argument(
        '--mode',
        choices=['start', 'center'],
        default='start',
        help='Alignment mode: start (keep start position) or center (center alignment)'
    )
    
    args = parser.parse_args()
    
    # Parse base BED
    chrom, base_start, base_end, disease = parse_bed(args.base_bed)
    base_size = base_end - base_start
    
    print(f"Base BED: {chrom}:{base_start}-{base_end} (size: {base_size/1_000_000:.2f} Mb)")
    print(f"Disease: {disease}")
    print(f"Mode: {args.mode}")
    print()
    
    # Parse sizes
    sizes_mb = [float(s.strip()) for s in args.sizes.split(',')]
    sizes_mb = sorted(sizes_mb)
    
    # Determine output directory: if output_prefix contains path, use bed/ directory
    script_dir = Path(__file__).parent
    output_prefix_path = Path(args.output_prefix)
    
    # If output_prefix has a parent directory, create it under bed/
    if output_prefix_path.parent != Path('.'):
        # e.g., "2q33/temp_WHS" -> "bed/2q33/temp_WHS"
        output_dir = script_dir / "bed" / output_prefix_path.parent
        output_prefix_name = output_prefix_path.name
    else:
        # e.g., "temp_WHS" -> "bed/temp_WHS"
        output_dir = script_dir / "bed"
        output_prefix_name = output_prefix_path.name
    
    # Create output directory if needed
    output_dir.mkdir(parents=True, exist_ok=True)
    
    created_files = []
    
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
        
        # Create BED file in bed/ directory
        output_file = output_dir / f"{output_prefix_name}_{int(size_mb)}Mb.bed"
        create_bed_file(output_file, chrom, new_start, new_end, disease, 
                       ['loss', 'overlap', '-', '-', '-', '-'])
        
        actual_size = (new_end - new_start) / 1_000_000
        print(f"Created: {output_file.name}")
        print(f"  Region: {chrom}:{new_start}-{new_end} (size: {actual_size:.2f} Mb)")
        created_files.append(output_file)
    
    print()
    print(f"Created {len(created_files)} BED files:")
    for f in created_files:
        print(f"  - {f}")
    
    return 0


if __name__ == '__main__':
    exit(main())

