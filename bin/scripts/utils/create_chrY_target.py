#!/usr/bin/env python3
"""
create_chrY_target.py - Create Y chromosome target regions BED file for NIPT analysis
This script generates a BED file containing Y-chromosome specific regions,
excluding pseudoautosomal regions (PARs) and other regions with mapping ambiguity.

Usage:
    python3 create_chrY_target.py --genome GENOME --output OUTPUT_FILE [--window_size SIZE]

Arguments:
    --genome      Genome reference (hg19 or hg38)
    --output      Output BED file
    --window_size Size of windows in bp (default: 10000)
"""

import argparse
import os
import sys
import logging
import random

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('create_chrY_target')

# Define pseudoautosomal regions (PARs) and other problematic regions on Y chromosome
# These regions should be excluded from the Y-specific target regions
EXCLUDE_REGIONS = {
    'hg19': [
        # PAR1
        ('chrY', 10000, 2649520, 'PAR1'),
        # PAR2
        ('chrY', 59034049, 59373566, 'PAR2'),
        # Centromere
        ('chrY', 10000000, 13800000, 'centromere'),
        # Heterochromatin
        ('chrY', 13800000, 20000000, 'heterochromatin'),
        # Other problematic regions (regions with known mapping ambiguity)
        ('chrY', 20000000, 20200000, 'ambiguous_region_1'),
        ('chrY', 40000000, 40100000, 'ambiguous_region_2'),
        ('chrY', 56000000, 56100000, 'ambiguous_region_3'),
    ],
    'hg38': [
        # PAR1
        ('chrY', 10000, 2781479, 'PAR1'),
        # PAR2
        ('chrY', 56887902, 57217415, 'PAR2'),
        # Centromere
        ('chrY', 10300000, 13400000, 'centromere'),
        # Heterochromatin
        ('chrY', 13400000, 19000000, 'heterochromatin'),
        # Other problematic regions
        ('chrY', 19000000, 19200000, 'ambiguous_region_1'),
        ('chrY', 38000000, 38100000, 'ambiguous_region_2'),
        ('chrY', 54000000, 54100000, 'ambiguous_region_3'),
    ]
}

# Y chromosome sizes
Y_SIZES = {
    'hg19': 59373566,
    'hg38': 57227415
}

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Create Y chromosome target BED file')
    parser.add_argument('--genome', required=True, choices=['hg19', 'hg38'], 
                        help='Genome reference (hg19 or hg38)')
    parser.add_argument('--output', required=True, help='Output BED file')
    parser.add_argument('--window_size', type=int, default=10000, 
                        help='Size of windows in bp (default: 10000)')
    parser.add_argument('--num_regions', type=int, default=100,
                        help='Number of target regions to generate (default: 100)')
    return parser.parse_args()

def is_excluded_region(start, end, exclude_regions):
    """Check if a region overlaps with excluded regions"""
    for _, ex_start, ex_end, _ in exclude_regions:
        if start < ex_end and end > ex_start:
            return True
    return False

def generate_y_specific_regions(genome, window_size, num_regions, exclude_regions):
    """Generate Y-specific regions, avoiding excluded regions"""
    y_regions = []
    
    # Get Y chromosome size
    y_size = Y_SIZES[genome]
    
    # Generate candidate regions
    candidate_starts = list(range(1, y_size - window_size, window_size))
    
    # Shuffle the start positions for random selection
    random.shuffle(candidate_starts)
    
    # Select regions, avoiding excluded regions
    for start in candidate_starts:
        if len(y_regions) >= num_regions:
            break
            
        end = start + window_size
        
        # Check if the region overlaps with excluded regions
        if is_excluded_region(start, end, exclude_regions):
            continue
            
        # Add the region
        region_name = f"chrY_{start}_{end}"
        y_regions.append(('chrY', start, end, region_name))
    
    # If we couldn't get enough regions with the current window size, try with a smaller one
    if len(y_regions) < num_regions and window_size > 1000:
        logger.warning(f"Could only find {len(y_regions)} regions with window size {window_size}")
        logger.warning(f"Trying with a smaller window size...")
        
        # Try with half the window size
        half_window = window_size // 2
        remaining_regions = num_regions - len(y_regions)
        
        additional_regions = generate_y_specific_regions(
            genome, half_window, remaining_regions, exclude_regions
        )
        
        y_regions.extend(additional_regions)
    
    return y_regions

def write_bed_file(regions, output_file):
    """Write regions to a BED file"""
    with open(output_file, 'w') as f:
        for chrom, start, end, name in regions:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\n")

def main():
    """Main function"""
    args = parse_args()
    
    logger.info(f"Creating Y chromosome target BED file for {args.genome}")
    logger.info(f"Window size: {args.window_size} bp")
    logger.info(f"Target number of regions: {args.num_regions}")
    
    # Generate Y-specific regions
    y_regions = generate_y_specific_regions(
        args.genome, 
        args.window_size, 
        args.num_regions, 
        EXCLUDE_REGIONS[args.genome]
    )
    
    # Write to BED file
    write_bed_file(y_regions, args.output)
    
    logger.info(f"Created {len(y_regions)} Y-specific regions")
    logger.info(f"Output written to {args.output}")
    
    # Calculate total base pairs
    total_bp = sum(end - start for _, start, end, _ in y_regions)
    logger.info(f"Total base pairs: {total_bp:,}")

if __name__ == "__main__":
    main()
