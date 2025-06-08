#!/usr/bin/env python3
"""
create_autosome_control.py - Create autosomal control regions BED file for NIPT analysis
This script generates a BED file containing autosomal control regions based on known
housekeeping genes and/or regions with consistent coverage across samples.

Usage:
    python3 create_autosome_control.py --genome GENOME --output OUTPUT_FILE [--window_size SIZE]

Arguments:
    --genome      Genome reference (hg19 or hg38)
    --output      Output BED file
    --window_size Size of windows in bp (default: 100000)
"""

import argparse
import os
import sys
import subprocess
import random
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('create_autosome_control')

# Define housekeeping genes (a small subset of commonly used ones)
# These are typically expressed at consistent levels across tissues and samples
HOUSEKEEPING_GENES = {
    'hg19': [
        # Gene name, chromosome, start, end
        ('ACTB', 'chr7', 5566782, 5570340),      # Beta-actin
        ('GAPDH', 'chr12', 6643585, 6647537),    # Glyceraldehyde-3-phosphate dehydrogenase
        ('B2M', 'chr15', 45003675, 45010357),    # Beta-2-microglobulin
        ('PPIA', 'chr7', 44836997, 44839439),    # Cyclophilin A
        ('RPL13A', 'chr19', 49993619, 49997295), # Ribosomal protein L13a
        ('LDHA', 'chr11', 18406614, 18410417),   # Lactate dehydrogenase A
        ('HPRT1', 'chrX', 133594175, 133634699), # Hypoxanthine phosphoribosyltransferase 1 (exclude for NIPT)
        ('TBP', 'chr6', 170554297, 170572870),   # TATA-box binding protein
        ('HMBS', 'chr11', 118962444, 118965359), # Hydroxymethylbilane synthase
        ('GUSB', 'chr7', 65425605, 65432600),    # Beta-glucuronidase
    ],
    'hg38': [
        # Gene name, chromosome, start, end - positions updated for hg38
        ('ACTB', 'chr7', 5527151, 5530709),
        ('GAPDH', 'chr12', 6534517, 6538371),
        ('B2M', 'chr15', 44711489, 44718171),
        ('PPIA', 'chr7', 44796585, 44798990),
        ('RPL13A', 'chr19', 49582411, 49586137),
        ('LDHA', 'chr11', 18394388, 18408379),
        ('TBP', 'chr6', 170554141, 170572907),
        ('HMBS', 'chr11', 118959400, 118962392),
        ('GUSB', 'chr7', 65960684, 65975385),
    ]
}

# Define autosome sizes (approximate)
AUTOSOME_SIZES = {
    'hg19': {
        'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276,
        'chr5': 180915260, 'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022,
        'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
        'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753,
        'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520,
        'chr21': 48129895, 'chr22': 51304566
    },
    'hg38': {
        'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
        'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
        'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
        'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
        'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
        'chr21': 46709983, 'chr22': 50818468
    }
}

# Regions to avoid (centromeres, telomeres, etc.)
PROBLEMATIC_REGIONS = {
    'hg19': [
        # Chromosome, start, end
        ('chr1', 121500000, 125000000),  # centromere
        ('chr2', 90500000, 96800000),    # centromere
        ('chr3', 87900000, 93900000),    # centromere
        ('chr4', 48200000, 52700000),    # centromere
        ('chr5', 46100000, 50700000),    # centromere
        ('chr6', 58700000, 63300000),    # centromere
        ('chr7', 58000000, 61700000),    # centromere
        ('chr8', 43100000, 48100000),    # centromere
        ('chr9', 47300000, 65900000),    # centromere
        ('chr10', 38000000, 42300000),   # centromere
        ('chr11', 51600000, 55700000),   # centromere
        ('chr12', 33300000, 38200000),   # centromere
        ('chr13', 16300000, 19500000),   # centromere
        ('chr14', 16100000, 19100000),   # centromere
        ('chr15', 15800000, 20700000),   # centromere
        ('chr16', 35300000, 38400000),   # centromere
        ('chr17', 22200000, 25800000),   # centromere
        ('chr18', 15400000, 19000000),   # centromere
        ('chr19', 24400000, 28600000),   # centromere
        ('chr20', 25600000, 29400000),   # centromere
        ('chr21', 10900000, 13200000),   # centromere
        ('chr22', 13700000, 17900000),   # centromere
    ],
    'hg38': [
        # Updated for hg38
        ('chr1', 121700000, 125100000),
        ('chr2', 91800000, 96800000),
        ('chr3', 87800000, 93800000),
        ('chr4', 48500000, 52700000),
        ('chr5', 46400000, 50700000),
        ('chr6', 58500000, 63300000),
        ('chr7', 58100000, 62100000),
        ('chr8', 43200000, 48600000),
        ('chr9', 42200000, 45500000),
        ('chr10', 38000000, 42300000),
        ('chr11', 51000000, 55800000),
        ('chr12', 33300000, 38200000),
        ('chr13', 16500000, 18900000),
        ('chr14', 16100000, 19100000),
        ('chr15', 17000000, 20500000),
        ('chr16', 35800000, 38600000),
        ('chr17', 22700000, 27400000),
        ('chr18', 15400000, 20000000),
        ('chr19', 24500000, 28600000),
        ('chr20', 25800000, 30400000),
        ('chr21', 10900000, 13200000),
        ('chr22', 13700000, 17900000),
    ]
}

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Create autosomal control BED file')
    parser.add_argument('--genome', required=True, choices=['hg19', 'hg38'], 
                        help='Genome reference (hg19 or hg38)')
    parser.add_argument('--output', required=True, help='Output BED file')
    parser.add_argument('--window_size', type=int, default=100000, 
                        help='Size of windows in bp (default: 100000)')
    parser.add_argument('--num_regions', type=int, default=200,
                        help='Number of control regions to generate (default: 200)')
    return parser.parse_args()

def is_problematic_region(chrom, start, end, problematic_regions):
    """Check if a region overlaps with problematic regions"""
    for p_chrom, p_start, p_end in problematic_regions:
        if chrom == p_chrom and start < p_end and end > p_start:
            return True
    return False

def generate_control_regions(genome, window_size, num_regions, problematic_regions):
    """Generate control regions based on known housekeeping genes and random autosomal regions"""
    control_regions = []
    
    # Add regions around housekeeping genes (excluding those on X and Y)
    for gene, chrom, start, end in HOUSEKEEPING_GENES[genome]:
        if chrom.startswith(('chrX', 'chrY')):
            continue  # Skip genes on sex chromosomes
        
        # Check if the gene region overlaps with problematic regions
        if is_problematic_region(chrom, start, end, problematic_regions):
            continue
            
        # Add the gene region
        control_regions.append((chrom, start, end, f"{gene}_gene"))
        
        # Add flanking regions (upstream and downstream)
        upstream_start = max(0, start - window_size)
        if not is_problematic_region(chrom, upstream_start, start, problematic_regions):
            control_regions.append((chrom, upstream_start, start, f"{gene}_upstream"))
            
        downstream_end = end + window_size
        if not is_problematic_region(chrom, end, downstream_end, problematic_regions):
            control_regions.append((chrom, end, downstream_end, f"{gene}_downstream"))
    
    # Add random autosomal regions to reach the desired number
    autosomes = [f"chr{i}" for i in range(1, 23)]
    
    while len(control_regions) < num_regions:
        # Randomly select an autosome
        chrom = random.choice(autosomes)
        
        # Get the chromosome size
        chrom_size = AUTOSOME_SIZES[genome][chrom]
        
        # Generate a random position, avoiding the ends of chromosomes
        max_start = chrom_size - window_size - 1000000  # Avoid the last 1Mb
        start = random.randint(1000000, max_start)  # Avoid the first 1Mb
        end = start + window_size
        
        # Check if the region overlaps with problematic regions
        if is_problematic_region(chrom, start, end, problematic_regions):
            continue
            
        # Check if the region overlaps with existing regions
        overlaps = False
        for existing_chrom, existing_start, existing_end, _ in control_regions:
            if chrom == existing_chrom and start < existing_end and end > existing_start:
                overlaps = True
                break
                
        if not overlaps:
            control_regions.append((chrom, start, end, f"random_{chrom}_{start}"))
    
    return control_regions

def write_bed_file(regions, output_file):
    """Write regions to a BED file"""
    with open(output_file, 'w') as f:
        for chrom, start, end, name in regions:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\n")

def main():
    """Main function"""
    args = parse_args()
    
    logger.info(f"Creating autosomal control BED file for {args.genome}")
    logger.info(f"Window size: {args.window_size} bp")
    logger.info(f"Target number of regions: {args.num_regions}")
    
    # Generate control regions
    control_regions = generate_control_regions(
        args.genome, 
        args.window_size, 
        args.num_regions, 
        PROBLEMATIC_REGIONS[args.genome]
    )
    
    # Write to BED file
    write_bed_file(control_regions, args.output)
    
    logger.info(f"Created {len(control_regions)} control regions")
    logger.info(f"Output written to {args.output}")

if __name__ == "__main__":
    main()
