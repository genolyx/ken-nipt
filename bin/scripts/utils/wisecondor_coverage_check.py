#!/usr/bin/env python3
"""
Wisecondor Coverage Check Tool for NIPT

Analyzes BAM read counts in bins to determine if Wisecondor/WisecondorX analysis is possible.
Designed for shallow depth WGS (0.2-0.3x average depth) used in NIPT.

Instead of per-base coverage, this tool:
1. Divides regions into bins (default 200kb, matching Wisecondor)
2. Counts reads in each bin
3. Determines if bins have sufficient reads for statistical analysis
4. Outputs coverable regions with adequate read counts

Usage:
    python3 wisecondor_coverage_check.py --bam proper_paired.bam --bed target_regions.bed --output coverable_regions.bed

Author: Ken
Version: 1.0
"""

import argparse
import sys
import os
from pathlib import Path
from collections import defaultdict

try:
    import pysam
    import numpy as np
    import pandas as pd
except ImportError as e:
    print(f"ERROR: Required library not found: {e}", file=sys.stderr)
    print("Please install required packages: pip install pysam numpy pandas", file=sys.stderr)
    sys.exit(1)


class WisecondorCoverageAnalyzer:
    """Analyzes BAM read counts in bins for Wisecondor/WisecondorX compatibility"""
    
    def __init__(self, bam_path, bin_size=200000, min_reads_per_bin=10, 
                 min_bin_fraction=0.8, min_total_reads=1000):
        """
        Initialize Wisecondor coverage analyzer
        
        Args:
            bam_path: Path to BAM file
            bin_size: Bin size for Wisecondor analysis (default: 200000)
            min_reads_per_bin: Minimum reads per bin (default: 10)
            min_bin_fraction: Minimum fraction of bins with sufficient reads (default: 0.8)
            min_total_reads: Minimum total reads in region (default: 1000)
        """
        self.bam_path = bam_path
        self.bin_size = bin_size
        self.min_reads_per_bin = min_reads_per_bin
        self.min_bin_fraction = min_bin_fraction
        self.min_total_reads = min_total_reads
        
        # Open BAM file
        try:
            self.bamfile = pysam.AlignmentFile(bam_path, "rb")
            if not self.bamfile.has_index():
                print(f"WARNING: BAM file is not indexed. Creating index...", file=sys.stderr)
                pysam.index(bam_path)
                self.bamfile = pysam.AlignmentFile(bam_path, "rb")
        except Exception as e:
            print(f"ERROR: Cannot open BAM file {bam_path}: {e}", file=sys.stderr)
            sys.exit(1)
        
        # Get total mapped reads for normalization info
        self.total_mapped_reads = self._get_total_mapped_reads()
    
    def _get_total_mapped_reads(self):
        """Get total number of mapped reads in BAM"""
        try:
            stats = pysam.idxstats(self.bam_path).split('\n')
            total = 0
            for line in stats:
                if line.strip():
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        total += int(parts[2])  # mapped reads
            return total
        except:
            return None
    
    def read_bed_file(self, bed_path):
        """
        Read BED file and return list of regions
        
        Returns:
            List of tuples (chrom, start, end, name)
        """
        regions = []
        try:
            with open(bed_path, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#') or line.startswith('track'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) < 3:
                        print(f"WARNING: Skipping invalid BED line {line_num}: {line}", file=sys.stderr)
                        continue
                    
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    name = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
                    
                    regions.append((chrom, start, end, name))
            
            print(f"Loaded {len(regions)} regions from {bed_path}")
            return regions
            
        except Exception as e:
            print(f"ERROR: Cannot read BED file {bed_path}: {e}", file=sys.stderr)
            sys.exit(1)
    
    def count_reads_in_bins(self, chrom, start, end):
        """
        Count reads in each bin for a region
        
        Returns:
            numpy array of read counts per bin
        """
        region_size = end - start
        num_bins = (region_size + self.bin_size - 1) // self.bin_size  # Ceiling division
        
        bin_counts = np.zeros(num_bins, dtype=np.int32)
        
        try:
            # Count reads for each bin
            for i in range(num_bins):
                bin_start = start + (i * self.bin_size)
                bin_end = min(start + ((i + 1) * self.bin_size), end)
                
                # Count reads in this bin
                count = self.bamfile.count(chrom, bin_start, bin_end)
                bin_counts[i] = count
            
            return bin_counts
            
        except Exception as e:
            print(f"WARNING: Error counting reads for {chrom}:{start}-{end}: {e}", file=sys.stderr)
            return np.zeros(num_bins, dtype=np.int32)
    
    def analyze_region_bins(self, chrom, start, end, name):
        """
        Analyze bin-level read counts for a region
        
        Returns:
            dict with bin analysis results
        """
        region_size = end - start
        
        # Get bin counts
        bin_counts = self.count_reads_in_bins(chrom, start, end)
        
        if len(bin_counts) == 0:
            return {
                'chrom': chrom,
                'start': start,
                'end': end,
                'name': name,
                'size': region_size,
                'num_bins': 0,
                'total_reads': 0,
                'mean_reads_per_bin': 0.0,
                'median_reads_per_bin': 0.0,
                'min_reads': 0,
                'max_reads': 0,
                'bins_with_reads': 0,
                'bins_sufficient': 0,
                'bin_fraction': 0.0,
                'is_analyzable': False,
                'reason': 'No bins',
                'bin_counts': [],
                'coverable_bins': []
            }
        
        # Calculate statistics
        num_bins = len(bin_counts)
        total_reads = np.sum(bin_counts)
        mean_reads = np.mean(bin_counts)
        median_reads = np.median(bin_counts)
        min_reads = np.min(bin_counts)
        max_reads = np.max(bin_counts)
        
        # Count bins with sufficient reads
        bins_with_reads = np.sum(bin_counts > 0)
        bins_sufficient = np.sum(bin_counts >= self.min_reads_per_bin)
        bin_fraction = bins_sufficient / num_bins if num_bins > 0 else 0.0
        
        # Find contiguous groups of sufficient bins
        coverable_bins = self._find_coverable_bin_groups(bin_counts, start)
        
        # Determine if region is analyzable
        is_analyzable = (
            total_reads >= self.min_total_reads and
            bin_fraction >= self.min_bin_fraction and
            num_bins >= 1
        )
        
        # Determine reason if not analyzable
        reason = "OK"
        if not is_analyzable:
            reasons = []
            if total_reads < self.min_total_reads:
                reasons.append(f"Low total reads ({total_reads} < {self.min_total_reads})")
            if bin_fraction < self.min_bin_fraction:
                reasons.append(f"Low bin fraction ({bin_fraction:.2%} < {self.min_bin_fraction:.2%})")
            if num_bins < 1:
                reasons.append(f"Too few bins ({num_bins})")
            reason = "; ".join(reasons)
        
        return {
            'chrom': chrom,
            'start': start,
            'end': end,
            'name': name,
            'size': region_size,
            'num_bins': num_bins,
            'total_reads': total_reads,
            'mean_reads_per_bin': mean_reads,
            'median_reads_per_bin': median_reads,
            'min_reads': min_reads,
            'max_reads': max_reads,
            'bins_with_reads': bins_with_reads,
            'bins_sufficient': bins_sufficient,
            'bin_fraction': bin_fraction,
            'is_analyzable': is_analyzable,
            'reason': reason,
            'bin_counts': bin_counts.tolist(),
            'coverable_bins': coverable_bins
        }
    
    def _find_coverable_bin_groups(self, bin_counts, region_start):
        """
        Find contiguous groups of bins with sufficient reads
        
        Args:
            bin_counts: numpy array of read counts
            region_start: genomic start position of region
        
        Returns:
            List of (start, end, num_bins) tuples for coverable groups
        """
        coverable_groups = []
        
        # Find stretches of sufficient bins
        sufficient = bin_counts >= self.min_reads_per_bin
        
        in_group = False
        group_start_idx = 0
        
        for i, is_sufficient in enumerate(sufficient):
            if is_sufficient and not in_group:
                # Start of a coverable group
                group_start_idx = i
                in_group = True
            elif not is_sufficient and in_group:
                # End of a coverable group
                group_end_idx = i
                num_group_bins = group_end_idx - group_start_idx
                
                # Only include groups with at least 1 bin
                if num_group_bins >= 1:
                    group_start = region_start + (group_start_idx * self.bin_size)
                    group_end = region_start + (group_end_idx * self.bin_size)
                    coverable_groups.append((group_start, group_end, num_group_bins))
                
                in_group = False
        
        # Handle case where group extends to the end
        if in_group:
            group_end_idx = len(sufficient)
            num_group_bins = group_end_idx - group_start_idx
            if num_group_bins >= 1:
                group_start = region_start + (group_start_idx * self.bin_size)
                group_end = region_start + (group_end_idx * self.bin_size)
                coverable_groups.append((group_start, group_end, num_group_bins))
        
        return coverable_groups
    
    def analyze_all_regions(self, regions):
        """
        Analyze all regions
        
        Returns:
            List of analysis results
        """
        results = []
        
        print(f"\nAnalyzing {len(regions)} regions...")
        print(f"Bin size: {self.bin_size:,} bp")
        print(f"Minimum reads per bin: {self.min_reads_per_bin}")
        print(f"Minimum bin fraction: {self.min_bin_fraction:.2%}")
        
        for i, (chrom, start, end, name) in enumerate(regions, 1):
            if i % 10 == 0 or i == len(regions):
                print(f"  Progress: {i}/{len(regions)} regions", end='\r')
            
            result = self.analyze_region_bins(chrom, start, end, name)
            results.append(result)
        
        print(f"  Progress: {len(regions)}/{len(regions)} regions - Complete!")
        return results
    
    def write_bed_output(self, results, output_path, include_partial=True):
        """
        Write coverable regions to BED file
        
        Args:
            results: List of analysis results
            output_path: Output BED file path
            include_partial: If True, include partially covered regions
        """
        try:
            with open(output_path, 'w') as f:
                # Write header
                f.write("# Wisecondor Coverage Analysis - Coverable Regions\n")
                f.write(f"# Bin size: {self.bin_size:,} bp\n")
                f.write(f"# Min reads per bin: {self.min_reads_per_bin}\n")
                f.write(f"# Min bin fraction: {self.min_bin_fraction:.2%}\n")
                f.write(f"# Min total reads: {self.min_total_reads}\n")
                f.write("#chrom\tstart\tend\tname\tscore\tstrand\ttotal_reads\tmean_reads\tbin_fraction\n")
                
                regions_written = 0
                
                for result in results:
                    if result['is_analyzable']:
                        # Write fully analyzable region
                        f.write(f"{result['chrom']}\t{result['start']}\t{result['end']}\t"
                               f"{result['name']}\t1000\t.\t"
                               f"{result['total_reads']}\t{result['mean_reads_per_bin']:.2f}\t"
                               f"{result['bin_fraction']:.4f}\n")
                        regions_written += 1
                    
                    elif include_partial and result['coverable_bins']:
                        # Write partially covered regions (coverable bin groups)
                        for j, (cov_start, cov_end, num_bins) in enumerate(result['coverable_bins'], 1):
                            # Count total reads in this group
                            group_reads = self.bamfile.count(result['chrom'], cov_start, cov_end)
                            mean_reads = group_reads / num_bins if num_bins > 0 else 0
                            
                            name = f"{result['name']}_partial{j}"
                            f.write(f"{result['chrom']}\t{cov_start}\t{cov_end}\t"
                                   f"{name}\t500\t.\t"
                                   f"{group_reads}\t{mean_reads:.2f}\t1.0000\n")
                            regions_written += 1
            
            print(f"\nWrote {regions_written} coverable regions to {output_path}")
            
        except Exception as e:
            print(f"ERROR: Cannot write output BED file: {e}", file=sys.stderr)
            sys.exit(1)
    
    def generate_report(self, results, report_path):
        """
        Generate detailed coverage report
        """
        try:
            # Summary statistics
            total_regions = len(results)
            analyzable_regions = sum(1 for r in results if r['is_analyzable'])
            partial_regions = sum(1 for r in results if not r['is_analyzable'] and r['coverable_bins'])
            non_analyzable = total_regions - analyzable_regions - partial_regions
            
            with open(report_path, 'w') as f:
                f.write("=" * 80 + "\n")
                f.write("Wisecondor Coverage Analysis Report (NIPT Shallow WGS)\n")
                f.write("=" * 80 + "\n\n")
                
                f.write(f"BAM File: {self.bam_path}\n")
                if self.total_mapped_reads:
                    f.write(f"Total Mapped Reads: {self.total_mapped_reads:,}\n")
                f.write(f"\nAnalysis Parameters:\n")
                f.write(f"  - Bin size: {self.bin_size:,} bp\n")
                f.write(f"  - Minimum reads per bin: {self.min_reads_per_bin}\n")
                f.write(f"  - Minimum bin fraction: {self.min_bin_fraction:.2%}\n")
                f.write(f"  - Minimum total reads: {self.min_total_reads}\n\n")
                
                f.write(f"Summary:\n")
                f.write(f"  - Total regions: {total_regions}\n")
                f.write(f"  - Fully analyzable: {analyzable_regions} ({analyzable_regions/total_regions*100:.1f}%)\n")
                f.write(f"  - Partially covered: {partial_regions} ({partial_regions/total_regions*100:.1f}%)\n")
                f.write(f"  - Not analyzable: {non_analyzable} ({non_analyzable/total_regions*100:.1f}%)\n\n")
                
                # Read count statistics
                total_reads_all = [r['total_reads'] for r in results]
                mean_reads_all = [r['mean_reads_per_bin'] for r in results]
                
                f.write(f"Read Count Statistics:\n")
                f.write(f"  - Total reads (avg per region): {np.mean(total_reads_all):.0f}\n")
                f.write(f"  - Total reads (median per region): {np.median(total_reads_all):.0f}\n")
                f.write(f"  - Mean reads per bin (avg): {np.mean(mean_reads_all):.2f}\n")
                f.write(f"  - Mean reads per bin (median): {np.median(mean_reads_all):.2f}\n\n")
                
                # Estimate average depth if we have total mapped reads
                if self.total_mapped_reads:
                    # Assume 150bp reads, hg19 genome size ~3Gb
                    estimated_depth = (self.total_mapped_reads * 150) / 3e9
                    f.write(f"Estimated Average Depth: {estimated_depth:.2f}x\n")
                    f.write(f"  (Based on {self.total_mapped_reads:,} mapped reads, 150bp reads, 3Gb genome)\n\n")
                
                # Per-region details
                f.write("=" * 80 + "\n")
                f.write("Per-Region Analysis\n")
                f.write("=" * 80 + "\n\n")
                
                for result in results:
                    status = "✓ PASS" if result['is_analyzable'] else "✗ FAIL"
                    f.write(f"{status} {result['name']} ({result['chrom']}:{result['start']}-{result['end']})\n")
                    f.write(f"  Size: {result['size']:,} bp\n")
                    f.write(f"  Number of bins: {result['num_bins']}\n")
                    f.write(f"  Total reads: {result['total_reads']:,}\n")
                    f.write(f"  Mean reads per bin: {result['mean_reads_per_bin']:.2f}\n")
                    f.write(f"  Median reads per bin: {result['median_reads_per_bin']:.2f}\n")
                    f.write(f"  Read range per bin: {result['min_reads']}-{result['max_reads']}\n")
                    f.write(f"  Bins with sufficient reads: {result['bins_sufficient']}/{result['num_bins']} "
                           f"({result['bin_fraction']:.2%})\n")
                    
                    if not result['is_analyzable']:
                        f.write(f"  Reason: {result['reason']}\n")
                        
                        if result['coverable_bins']:
                            f.write(f"  Coverable bin groups: {len(result['coverable_bins'])}\n")
                            for i, (start, end, num_bins) in enumerate(result['coverable_bins'], 1):
                                f.write(f"    {i}. {result['chrom']}:{start}-{end} ({num_bins} bins, {end-start:,} bp)\n")
                    
                    # Show bin-level details for small regions
                    if result['num_bins'] <= 20:
                        f.write(f"  Bin-level reads: {result['bin_counts']}\n")
                    
                    f.write("\n")
            
            print(f"Detailed report written to {report_path}")
            
        except Exception as e:
            print(f"ERROR: Cannot generate report: {e}", file=sys.stderr)
    
    def close(self):
        """Close BAM file"""
        if hasattr(self, 'bamfile'):
            self.bamfile.close()


def main():
    parser = argparse.ArgumentParser(
        description='Wisecondor Coverage Check Tool - For NIPT shallow depth WGS (0.2-0.3x)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (default bin size 200kb)
  python3 wisecondor_coverage_check.py --bam sample.proper_paired.bam --bed MD_Target_8.bed --output coverable.bed
  
  # With custom parameters
  python3 wisecondor_coverage_check.py --bam sample.bam --bed targets.bed --output output.bed \\
      --bin-size 200000 --min-reads-per-bin 10 --min-bin-fraction 0.8
  
  # Generate detailed report
  python3 wisecondor_coverage_check.py --bam sample.bam --bed targets.bed --output output.bed --report

Note: This tool is designed for shallow depth WGS (0.2-0.3x average depth) used in NIPT.
      It analyzes read counts per bin, not per-base coverage.
        """
    )
    
    # Required arguments
    parser.add_argument('--bam', required=True, help='Input BAM file (proper_paired.bam recommended)')
    parser.add_argument('--bed', required=True, help='Input BED file with target regions')
    parser.add_argument('--output', '-o', required=True, help='Output BED file for coverable regions')
    
    # Optional arguments
    parser.add_argument('--bin-size', type=int, default=200000,
                       help='Bin size for Wisecondor analysis in bp (default: 200000)')
    parser.add_argument('--min-reads-per-bin', type=int, default=10,
                       help='Minimum reads per bin (default: 10)')
    parser.add_argument('--min-bin-fraction', type=float, default=0.8,
                       help='Minimum fraction of bins with sufficient reads (default: 0.8)')
    parser.add_argument('--min-total-reads', type=int, default=1000,
                       help='Minimum total reads in region (default: 1000)')
    parser.add_argument('--report', action='store_true',
                       help='Generate detailed coverage report')
    parser.add_argument('--no-partial', action='store_true',
                       help='Do not include partially covered regions in output')
    
    args = parser.parse_args()
    
    # Validate arguments
    if not os.path.exists(args.bam):
        print(f"ERROR: BAM file not found: {args.bam}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.bed):
        print(f"ERROR: BED file not found: {args.bed}", file=sys.stderr)
        sys.exit(1)
    
    if args.min_bin_fraction < 0 or args.min_bin_fraction > 1:
        print(f"ERROR: Bin fraction must be between 0 and 1", file=sys.stderr)
        sys.exit(1)
    
    # Print configuration
    print("=" * 80)
    print("Wisecondor Coverage Check Tool (NIPT Shallow WGS)")
    print("=" * 80)
    print(f"BAM file: {args.bam}")
    print(f"BED file: {args.bed}")
    print(f"Output file: {args.output}")
    print(f"")
    print(f"Parameters:")
    print(f"  - Bin size: {args.bin_size:,} bp")
    print(f"  - Minimum reads per bin: {args.min_reads_per_bin}")
    print(f"  - Minimum bin fraction: {args.min_bin_fraction:.2%}")
    print(f"  - Minimum total reads: {args.min_total_reads:,}")
    print(f"  - Include partial regions: {not args.no_partial}")
    print("=" * 80)
    
    # Initialize analyzer
    analyzer = WisecondorCoverageAnalyzer(
        args.bam,
        bin_size=args.bin_size,
        min_reads_per_bin=args.min_reads_per_bin,
        min_bin_fraction=args.min_bin_fraction,
        min_total_reads=args.min_total_reads
    )
    
    try:
        # Read BED file
        regions = analyzer.read_bed_file(args.bed)
        
        # Analyze coverage
        results = analyzer.analyze_all_regions(regions)
        
        # Write output BED file
        analyzer.write_bed_output(results, args.output, include_partial=not args.no_partial)
        
        # Generate report if requested
        if args.report:
            report_path = args.output.replace('.bed', '_report.txt')
            analyzer.generate_report(results, report_path)
        
        # Print summary
        print("\n" + "=" * 80)
        print("Analysis Complete!")
        print("=" * 80)
        
        total = len(results)
        analyzable = sum(1 for r in results if r['is_analyzable'])
        partial = sum(1 for r in results if not r['is_analyzable'] and r['coverable_bins'])
        
        print(f"Total regions: {total}")
        print(f"  ✓ Fully analyzable: {analyzable} ({analyzable/total*100:.1f}%)")
        print(f"  ⚠ Partially covered: {partial} ({partial/total*100:.1f}%)")
        print(f"  ✗ Not analyzable: {total-analyzable-partial} ({(total-analyzable-partial)/total*100:.1f}%)")
        print(f"\nWisecondor/WisecondorX Analysis:")
        if analyzable == total:
            print("  ✓ All regions are suitable for analysis")
        elif analyzable > 0:
            print(f"  ⚠ {analyzable}/{total} regions are suitable for analysis")
            print(f"  ⚠ Consider using coverable regions from {args.output}")
        else:
            print("  ✗ No regions have sufficient read counts")
            print("  ✗ Check BAM quality and sequencing depth")
        
        # Show estimated depth if available
        if analyzer.total_mapped_reads:
            estimated_depth = (analyzer.total_mapped_reads * 150) / 3e9
            print(f"\nEstimated average depth: {estimated_depth:.2f}x")
            if estimated_depth < 0.1:
                print("  ⚠ Warning: Very low depth, may affect analysis quality")
        
    finally:
        analyzer.close()


if __name__ == '__main__':
    main()

