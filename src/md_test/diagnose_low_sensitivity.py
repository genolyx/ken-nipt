#!/usr/bin/env python3
"""
Diagnose Low Sensitivity Issues

This script checks:
1. Z-score distribution in detected vs non-detected samples
2. Actual deletion region coverage in BAM files
3. Validation ratios from JSON metadata
"""

import json
import argparse
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def parse_deletion_from_json(json_file: Path) -> Optional[Dict]:
    """Parse deletion information from JSON file"""
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        return data.get('deletion')
    except Exception as e:
        logger.warning(f"Failed to parse JSON {json_file}: {e}")
        return None


def parse_wcx_zscore(bed_file: Path, return_raw: bool = False) -> List[float]:
    """Extract all z-scores from WCX aberrations.bed file
    
    WCX aberrations.bed format: chr start end effect zscore ratio length [disease-name]
    z-score is in column 5 (index 4)
    Note: For deletions, z-score is typically negative. 
    
    Args:
        bed_file: Path to aberrations.bed file
        return_raw: If True, return raw z-scores (including negative). If False, return absolute values.
    
    Returns:
        List of z-scores (absolute values by default, or raw if return_raw=True)
    """
    zscores = []
    if not bed_file.exists():
        return zscores
    
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                fields = line.split('\t')
                if len(fields) < 5:
                    continue
                
                # z-score is in column 5 (index 4)
                # Format: chr start end effect zscore ratio length [disease-name]
                try:
                    if len(fields) >= 5:
                        zscore_val = float(fields[4])  # 5th column is z-score
                        # For deletions, z-score is negative
                        if return_raw:
                            zscores.append(zscore_val)
                        else:
                            zscores.append(abs(zscore_val))  # Use absolute value for magnitude
                except (ValueError, IndexError):
                    # If column 4 is not a number, try to find z-score column
                    # z-score for deletions is typically negative and large magnitude
                    # Ratio is typically between 0 and 1
                    for i in range(3, min(len(fields), 7)):
                        try:
                            val = float(fields[i])
                            # Check if it's likely z-score (large magnitude) vs ratio (0-1 range)
                            if abs(val) >= 2.0:  # Likely z-score (deletions are negative, large magnitude)
                                if return_raw:
                                    zscores.append(val)
                                else:
                                    zscores.append(abs(val))
                                break
                        except ValueError:
                            continue
    except Exception as e:
        logger.warning(f"Failed to parse WCX bed file {bed_file}: {e}")
    
    return zscores


def parse_wc_zscore(report_file: Path) -> List[float]:
    """Extract z-scores from WC report.txt file"""
    zscores = []
    if not report_file.exists():
        return zscores
    
    try:
        with open(report_file, 'r') as f:
            lines = f.readlines()
        
        in_test_section = False
        for line in lines:
            if '# Test results:' in line or '# Test results' in line:
                in_test_section = True
                continue
            
            if in_test_section:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                if 'z-score' in line.lower() and 'effect' in line.lower():
                    continue
                
                fields = line.split()
                if len(fields) >= 4:
                    try:
                        zscore = abs(float(fields[0]))
                        zscores.append(zscore)
                    except ValueError:
                        continue
    except Exception as e:
        logger.warning(f"Failed to parse WC report file {report_file}: {e}")
    
    return zscores


def check_overlap(deletion: Dict, detected_start: int, detected_end: int) -> bool:
    """Check if detected region overlaps with expected deletion"""
    overlap_start = max(deletion['start'], detected_start)
    overlap_end = min(deletion['end'], detected_end)
    
    if overlap_start < overlap_end:
        overlap_size = overlap_end - overlap_start
        deletion_size = deletion['end'] - deletion['start']
        detected_size = detected_end - detected_start
        
        overlap_ratio_deletion = overlap_size / deletion_size if deletion_size > 0 else 0
        overlap_ratio_detected = overlap_size / detected_size if detected_size > 0 else 0
        
        return overlap_ratio_deletion >= 0.5 or overlap_ratio_detected >= 0.5
    
    return False


def get_detection_status(sample_dir: Path, sample_name: str, deletion: Dict, method: str, output_type: str) -> Dict:
    """
    Check detection status and extract z-scores
    Returns: {'detected': bool, 'zscore': float or None, 'max_zscore': float or None}
    """
    result = {'detected': False, 'zscore': None, 'max_zscore': None, 'overlaps': False}
    
    if method == 'WCX':
        bed_file = sample_dir / f"Output_WCX" / output_type / f"{sample_name}.wcx.{output_type}_aberrations.bed"
        zscores = parse_wcx_zscore(bed_file, return_raw=False)
        zscores_raw = parse_wcx_zscore(bed_file, return_raw=True)  # Get raw z-scores (with sign)
        
        if zscores:
            result['max_zscore'] = max(zscores)
            # Store raw z-score info for debugging
            if zscores_raw:
                result['min_raw_zscore'] = min(zscores_raw)
                result['max_raw_zscore'] = max(zscores_raw)
            
            # Check if any detected region overlaps with expected deletion
            if bed_file.exists():
                try:
                    with open(bed_file, 'r') as f:
                        for line in f:
                            line = line.strip()
                            if not line or line.startswith('#'):
                                continue
                            
                            fields = line.split('\t')
                            if len(fields) < 3:
                                continue
                            
                            chr_name = fields[0]
                            start = int(fields[1])
                            end = int(fields[2])
                            
                            if chr_name == str(deletion['chromosome']):
                                if check_overlap(deletion, start, end):
                                    result['detected'] = True
                                    result['overlaps'] = True
                                    # Get z-score for this overlapping region
                                    # z-score is in column 5 (index 4)
                                    # Format: chr start end effect zscore ratio length [disease-name]
                                    if len(fields) >= 5:
                                        try:
                                            zscore_val = abs(float(fields[4]))  # 5th column is z-score
                                            if result['zscore'] is None or zscore_val > result['zscore']:
                                                result['zscore'] = zscore_val
                                        except (ValueError, IndexError):
                                            # Fallback: try to find z-score (value > 2 or < -2)
                                            for i in range(3, min(len(fields), 7)):
                                                try:
                                                    zscore_val = abs(float(fields[i]))
                                                    if zscore_val >= 2.0:
                                                        if result['zscore'] is None or zscore_val > result['zscore']:
                                                            result['zscore'] = zscore_val
                                                        break
                                                except ValueError:
                                                    continue
                                    break
                except Exception as e:
                    logger.debug(f"Error checking overlap: {e}")
    else:  # WC
        report_file = sample_dir / f"Output_WC" / output_type / f"{sample_name}.wc.{output_type}.report.txt"
        zscores = parse_wc_zscore(report_file)
        
        if zscores:
            result['max_zscore'] = max(zscores)
            
            # Check if any detected region overlaps
            if report_file.exists():
                try:
                    with open(report_file, 'r') as f:
                        lines = f.readlines()
                    
                    in_test_section = False
                    for line in lines:
                        if '# Test results:' in line or '# Test results' in line:
                            in_test_section = True
                            continue
                        
                        if in_test_section:
                            line = line.strip()
                            if not line or line.startswith('#'):
                                continue
                            
                            if 'z-score' in line.lower() and 'effect' in line.lower():
                                continue
                            
                            fields = line.split()
                            if len(fields) >= 4:
                                try:
                                    zscore = abs(float(fields[0]))
                                    location = fields[3]
                                    
                                    if ':' in location:
                                        chr_part, coord_part = location.split(':', 1)
                                        if '-' in coord_part:
                                            start_str, end_str = coord_part.split('-', 1)
                                            chr_name = chr_part
                                            start = int(start_str)
                                            end = int(end_str)
                                            
                                            if chr_name == str(deletion['chromosome']):
                                                if check_overlap(deletion, start, end):
                                                    result['detected'] = True
                                                    result['overlaps'] = True
                                                    if result['zscore'] is None or zscore > result['zscore']:
                                                        result['zscore'] = zscore
                                except (ValueError, IndexError):
                                    continue
                except Exception as e:
                    logger.debug(f"Error checking overlap: {e}")
    
    return result


def scan_samples(root_dir: Path, work_dir: str, sample_list: Optional[List[str]] = None) -> List[Dict]:
    """Scan directories for samples"""
    samples = []
    analysis_dir = root_dir / "analysis" / work_dir
    
    if not analysis_dir.exists():
        logger.error(f"Analysis directory not found: {analysis_dir}")
        return samples
    
    sample_dirs = [d for d in analysis_dir.iterdir() if d.is_dir()]
    
    if sample_list:
        sample_dirs = [d for d in sample_dirs if d.name in sample_list]
    
    for sample_dir in sample_dirs:
        sample_name = sample_dir.name
        json_file = sample_dir / f"{sample_name}.json"
        
        if not json_file.exists():
            continue
        
        deletion = parse_deletion_from_json(json_file)
        if not deletion:
            continue
        
        # Extract metadata from JSON
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            validation = data.get('deletion_validation', {})
            ff = data.get('target_parameters', {}).get('ff_target_percent')
            if ff is None:
                ff = data.get('target_parameters', {}).get('ff_target')
            
            samples.append({
                'sample_name': sample_name,
                'sample_dir': sample_dir,
                'deletion': deletion,
                'ff': ff,
                'deletion_length_mb': deletion.get('size_mb') or (deletion.get('size_bp', 0) / 1_000_000),
                'upstream_ratio': validation.get('upstream_ratio'),
                'deletion_ratio': validation.get('deletion_ratio'),
                'downstream_ratio': validation.get('downstream_ratio')
            })
        except Exception as e:
            logger.warning(f"Failed to parse JSON metadata for {sample_name}: {e}")
    
    logger.info(f"Found {len(samples)} samples")
    return samples


def main():
    parser = argparse.ArgumentParser(
        description="Diagnose low sensitivity issues in MD detection"
    )
    parser.add_argument(
        '-r', '--root-dir',
        type=str,
        default='/home/ken/ken-nipt',
        help='Root directory path'
    )
    parser.add_argument(
        '-w', '--work-dir',
        type=str,
        required=True,
        help='Work directory (e.g., md_validation/1p36)'
    )
    parser.add_argument(
        '-o', '--output-dir',
        type=str,
        required=True,
        help='Output directory for results'
    )
    parser.add_argument(
        '--sample-sheet',
        type=str,
        help='Optional sample sheet TSV file to filter samples'
    )
    parser.add_argument(
        '--ff-filter',
        type=int,
        help='Filter by FF value (e.g., 15)'
    )
    parser.add_argument(
        '--analyze-existing',
        type=str,
        help='Analyze z-score distribution from existing diagnosis_results.tsv file'
    )
    
    args = parser.parse_args()
    
    # Analyze existing results file if requested
    if args.analyze_existing:
        results_file = Path(args.analyze_existing)
        if not results_file.exists():
            logger.error(f"Results file not found: {results_file}")
            return
        
        output_dir = Path(args.output_dir) if args.output_dir else results_file.parent
        analyze_zscore_distribution(results_file, output_dir)
        return
    
    # Debug: Show file format if requested
    if args.debug_file:
        debug_path = Path(args.debug_file)
        if debug_path.exists():
            print(f"\n{'='*80}")
            print(f"Debug: Showing first 5 lines of: {debug_path}")
            print(f"{'='*80}")
            with open(debug_path, 'r') as f:
                for i, line in enumerate(f):
                    if i >= 5:
                        break
                    print(f"Line {i+1}: {line.rstrip()}")
                    if i == 0:
                        # Show column breakdown
                        fields = line.rstrip().split('\t')
                        print(f"  Columns: {len(fields)}")
                        for j, field in enumerate(fields[:8]):  # Show first 8 columns
                            print(f"    [{j}]: '{field}'")
            print(f"{'='*80}\n")
            return
    
    # Parse sample list if provided
    sample_list = None
    if args.sample_sheet:
        import csv
        sample_list = []
        with open(args.sample_sheet, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                sample_list.append(row['SAMPLE_NAME'])
        logger.info(f"Filtering to {len(sample_list)} samples from sample sheet")
    
    # Scan samples
    root_dir = Path(args.root_dir)
    samples = scan_samples(root_dir, args.work_dir, sample_list)
    
    if len(samples) == 0:
        logger.error("No samples found!")
        return
    
    # Filter by FF and deletion length if specified
    if args.ff_filter:
        samples = [s for s in samples if s['ff'] == args.ff_filter]
    
    if args.del_length_filter:
        samples = [s for s in samples if abs(s['deletion_length_mb'] - args.del_length_filter) < 0.1]
    
    logger.info(f"Analyzing {len(samples)} samples")
    
    # Analyze each method
    methods = ['WC_orig', 'WC_fetus', 'WCX_orig', 'WCX_fetus']
    results = []
    
    for sample in samples:
        for method_spec in methods:
            method, output_type = method_spec.split('_')
            
            detection = get_detection_status(
                sample['sample_dir'],
                sample['sample_name'],
                sample['deletion'],
                method,
                output_type
            )
            
            results.append({
                'sample_name': sample['sample_name'],
                'method': method_spec,
                'ff': sample['ff'],
                'deletion_length_mb': sample['deletion_length_mb'],
                'detected': detection['detected'],
                'overlaps': detection['overlaps'],
                'zscore': detection['zscore'],
                'max_zscore': detection['max_zscore'],
                'upstream_ratio': sample['upstream_ratio'],
                'deletion_ratio': sample['deletion_ratio'],
                'downstream_ratio': sample['downstream_ratio']
            })
    
    df = pd.DataFrame(results)
    
    # Save results
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results_file = output_dir / 'diagnosis_results.tsv'
    df.to_csv(results_file, sep='\t', index=False)
    logger.info(f"Saved results to {results_file}")
    
    # Print summary statistics
    print("\n" + "="*80)
    print("DIAGNOSIS SUMMARY")
    print("="*80)
    
    # Group by FF and deletion length
    for method in methods:
        method_df = df[df['method'] == method]
        
        print(f"\n{method}:")
        print("-" * 80)
        
        if args.ff_filter or args.del_length_filter:
            # Show detailed breakdown
            grouped = method_df.groupby(['ff', 'deletion_length_mb'])
            for (ff, del_len), group in grouped:
                total = len(group)
                detected = group['detected'].sum()
                sensitivity = detected / total if total > 0 else 0
                
                # Z-score statistics
                zscores_detected = group[group['detected']]['zscore'].dropna()
                zscores_all = group['max_zscore'].dropna()
                
                print(f"  FF={ff}%, Del={del_len}Mb: "
                      f"Sensitivity={sensitivity:.2%} ({detected}/{total})")
                
                if len(zscores_detected) > 0:
                    print(f"    Detected z-scores: min={zscores_detected.min():.2f}, "
                          f"mean={zscores_detected.mean():.2f}, max={zscores_detected.max():.2f}")
                
                if len(zscores_all) > 0:
                    print(f"    All max z-scores (abs): min={zscores_all.min():.2f}, "
                          f"mean={zscores_all.mean():.2f}, max={zscores_all.max():.2f}")
                    # Note: For deletions, z-scores are typically negative
                    # The absolute values are shown above for threshold comparison
                
                # Validation ratios
                del_ratios = group['deletion_ratio'].dropna()
                if len(del_ratios) > 0:
                    print(f"    Deletion ratios: min={del_ratios.min():.4f}, "
                          f"mean={del_ratios.mean():.4f}, max={del_ratios.max():.4f}")
        else:
            # Show overall statistics
            total = len(method_df)
            detected = method_df['detected'].sum()
            sensitivity = detected / total if total > 0 else 0
            
            zscores_detected = method_df[method_df['detected']]['zscore'].dropna()
            zscores_all = method_df['max_zscore'].dropna()
            
            print(f"  Overall Sensitivity: {sensitivity:.2%} ({detected}/{total})")
            
            if len(zscores_detected) > 0:
                print(f"  Detected z-scores: min={zscores_detected.min():.2f}, "
                      f"mean={zscores_detected.mean():.2f}, max={zscores_detected.max():.2f}")
            
            if len(zscores_all) > 0:
                print(f"  All max z-scores (abs): min={zscores_all.min():.2f}, "
                      f"mean={zscores_all.mean():.2f}, max={zscores_all.max():.2f}")
                # Note: For deletions, z-scores are typically negative
                # The absolute values are shown above for threshold comparison
                # If you see values like 17.31 with many decimals, it might be reading ratio column instead
                print(f"  Note: If z-scores look like ratios (0-1 range with many decimals),")
                print(f"        the file format might be different. Check actual file format.")
    
def analyze_zscore_distribution(results_file: Path, output_dir: Path):
    """Analyze z-score distribution from saved diagnosis results"""
    import pandas as pd
    import numpy as np
    
    # Read results
    df = pd.read_csv(results_file, sep='\t')
    
    print("\n" + "="*80)
    print("Z-SCORE DISTRIBUTION ANALYSIS")
    print("="*80)
    
    methods = ['WC_orig', 'WC_fetus', 'WCX_orig', 'WCX_fetus']
    
    for method in methods:
        method_df = df[df['method'] == method]
        
        print(f"\n{method}:")
        print("-" * 80)
        
        # Overall statistics
        zscores_all = method_df['max_zscore'].dropna()
        zscores_detected = method_df[method_df['detected']]['zscore'].dropna()
        zscores_not_detected = method_df[~method_df['detected']]['max_zscore'].dropna()
        
        print(f"  Total samples: {len(method_df)}")
        print(f"  Samples with z-score data: {len(zscores_all)}")
        print(f"  Detected samples: {len(zscores_detected)}")
        print(f"  Not detected samples: {len(zscores_not_detected)}")
        
        if len(zscores_all) > 0:
            print(f"\n  All z-scores (abs):")
            print(f"    Count: {len(zscores_all)}")
            print(f"    Min: {zscores_all.min():.2f}")
            print(f"    25th percentile: {zscores_all.quantile(0.25):.2f}")
            print(f"    Median: {zscores_all.median():.2f}")
            print(f"    75th percentile: {zscores_all.quantile(0.75):.2f}")
            print(f"    95th percentile: {zscores_all.quantile(0.95):.2f}")
            print(f"    Max: {zscores_all.max():.2f}")
            print(f"    Mean: {zscores_all.mean():.2f}")
            print(f"    Std: {zscores_all.std():.2f}")
        
        if len(zscores_detected) > 0:
            print(f"\n  Detected z-scores (abs):")
            print(f"    Count: {len(zscores_detected)}")
            print(f"    Min: {zscores_detected.min():.2f}")
            print(f"    Median: {zscores_detected.median():.2f}")
            print(f"    Max: {zscores_detected.max():.2f}")
            print(f"    Mean: {zscores_detected.mean():.2f}")
        
        if len(zscores_not_detected) > 0:
            print(f"\n  Not detected z-scores (abs):")
            print(f"    Count: {len(zscores_not_detected)}")
            print(f"    Min: {zscores_not_detected.min():.2f}")
            print(f"    Median: {zscores_not_detected.median():.2f}")
            print(f"    Max: {zscores_not_detected.max():.2f}")
            print(f"    Mean: {zscores_not_detected.mean():.2f}")
        
        # Threshold analysis
        thresholds = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0]
        print(f"\n  Threshold analysis:")
        for threshold in thresholds:
            above_threshold = (zscores_all >= threshold).sum()
            pct = (above_threshold / len(zscores_all) * 100) if len(zscores_all) > 0 else 0
            detected_above = (zscores_detected >= threshold).sum() if len(zscores_detected) > 0 else 0
            not_detected_above = (zscores_not_detected >= threshold).sum() if len(zscores_not_detected) > 0 else 0
            print(f"    z-score >= {threshold:.1f}: {above_threshold}/{len(zscores_all)} ({pct:.1f}%) "
                  f"[Detected: {detected_above}, Not detected: {not_detected_above}]")
        
        # Group by FF and deletion length
        print(f"\n  Z-score by FF and deletion length:")
        grouped = method_df.groupby(['ff', 'deletion_length_mb'])['max_zscore']
        for (ff, del_len), group in grouped:
            zscores_group = group.dropna()
            if len(zscores_group) > 0:
                detected_count = method_df[(method_df['ff'] == ff) & 
                                          (method_df['deletion_length_mb'] == del_len) & 
                                          method_df['detected']].shape[0]
                total_count = len(zscores_group)
                print(f"    FF={ff}%, Del={del_len}Mb: "
                      f"mean={zscores_group.mean():.2f}, "
                      f"median={zscores_group.median():.2f}, "
                      f"min={zscores_group.min():.2f}, "
                      f"max={zscores_group.max():.2f} "
                      f"[Detected: {detected_count}/{total_count}]")
    
    # Save detailed distribution to file
    dist_file = output_dir / 'zscore_distribution.tsv'
    dist_data = []
    
    for method in methods:
        method_df = df[df['method'] == method]
        zscores_all = method_df['max_zscore'].dropna()
        
        if len(zscores_all) > 0:
            for threshold in [3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0]:
                above_count = (zscores_all >= threshold).sum()
                dist_data.append({
                    'method': method,
                    'threshold': threshold,
                    'above_threshold': above_count,
                    'total': len(zscores_all),
                    'percentage': above_count / len(zscores_all) * 100 if len(zscores_all) > 0 else 0
                })
    
    if dist_data:
        dist_df = pd.DataFrame(dist_data)
        dist_df.to_csv(dist_file, sep='\t', index=False)
        print(f"\n  Saved threshold distribution to: {dist_file}")
    
    print("\n" + "="*80)


if __name__ == '__main__':
    main()

