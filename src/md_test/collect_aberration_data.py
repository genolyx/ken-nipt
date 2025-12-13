#!/usr/bin/env python3
"""
Collect Aberration Data from WC and WCX Results

Extracts z-scores and detected aberration lengths from WC report.txt and WCX aberrations.bed files.
This data is used for ROC and Performance analysis with Min Detect Length filtering.

Usage:
    python collect_aberration_data.py \
        --data-dir /data/md_validation/analysis_result \
        --output /data/md_validation/roc_results/aberration_data.csv
"""

import argparse
import logging
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional
import re

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

# Disease to chromosome mapping
DISEASE_CHROMOSOMES = {
    '1p36': 'chr1',
    '2q33': 'chr2',
    'CDC': 'chr5',  # Cri Du Chat syndrome: chr5p deletion
    'DGS': 'chr22',  # DiGeorge syndrome: chr22q11.2 deletion
    'Jacobsen': 'chr11',
    'PWS': 'chr15',
    'WBS': 'chr7',
    'WHS': 'chr4'
}


def parse_sample_name(sample_name: str) -> Dict:
    """Parse sample name to extract metadata"""
    # Example: 4_3_1p36deletionsyndrome_FF10_15M_3Mb_F
    parts = sample_name.split('_')
    
    metadata = {
        'sample_id': sample_name,
        'ff': None,
        'deletion_length_mb': None,
        'gender': None,
        'disease': None
    }
    
    # Extract disease name
    for part in parts:
        if 'deletionsyndrome' in part.lower() or 'syndrome' in part.lower():
            disease_part = part.lower()
            if '1p36' in disease_part:
                metadata['disease'] = '1p36'
            elif '2q33' in disease_part:
                metadata['disease'] = '2q33'
            elif 'digeorge' in disease_part or 'dgs' in disease_part:
                metadata['disease'] = 'DGS'
            elif 'cri' in disease_part or 'cdc' in disease_part or 'criduchat' in disease_part:
                metadata['disease'] = 'CDC'
            elif 'jacobsen' in disease_part:
                metadata['disease'] = 'Jacobsen'
            elif 'prader' in disease_part or 'pws' in disease_part:
                metadata['disease'] = 'PWS'
            elif 'williams' in disease_part or 'wbs' in disease_part:
                metadata['disease'] = 'WBS'
            elif 'wolf' in disease_part or 'whs' in disease_part:
                metadata['disease'] = 'WHS'
    
    # Extract FF
    for i, part in enumerate(parts):
        if part.upper().startswith('FF'):
            ff_str = part[2:]
            try:
                metadata['ff'] = float(ff_str)
            except:
                pass
    
    # Extract deletion length (Mb)
    for i, part in enumerate(parts):
        if 'Mb' in part:
            mb_str = part.replace('Mb', '').replace('_', '.')
            try:
                mb_val = float(mb_str)
                # "0Mb" in sample name actually means 0.5Mb deletion
                if mb_val == 0.0:
                    metadata['deletion_length_mb'] = 0.5
                else:
                    metadata['deletion_length_mb'] = mb_val
            except:
                pass
    
    # Extract gender (last part)
    if len(parts) > 0:
        last_part = parts[-1]
        if last_part in ['F', 'M']:
            metadata['gender'] = last_part
    
    return metadata


def parse_wc_report(report_file: Path, target_chr: str) -> Optional[Dict]:
    """Parse WC report.txt to extract z-score and detected length for target chromosome"""
    
    if not report_file.exists():
        return None
    
    try:
        with open(report_file, 'r') as f:
            content = f.read()
        
        # Find "Test results" section
        test_results_match = re.search(r'# Test results:.*?\n(.*?)(?:\n\n|$)', content, re.DOTALL)
        if not test_results_match:
            return None
        
        test_results = test_results_match.group(1)
        
        # Parse lines (skip header)
        lines = test_results.strip().split('\n')
        
        # Find header line
        header_idx = -1
        for i, line in enumerate(lines):
            if 'z-score' in line and 'effect' in line and 'mbsize' in line and 'location' in line:
                header_idx = i
                break
        
        if header_idx == -1:
            return None
        
        # Parse data lines
        best_aberration = None
        max_abs_zscore = 0
        
        for line in lines[header_idx + 1:]:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split()
            if len(parts) < 4:
                continue
            
            try:
                zscore = float(parts[0])
                effect = float(parts[1])
                mbsize = float(parts[2])
                location = parts[3]
                
                # Extract chromosome from location (e.g., "1:800000-8200000")
                chr_match = re.match(r'(\d+|X|Y):', location)
                if not chr_match:
                    continue
                
                chr_num = chr_match.group(1)
                chr_name = f'chr{chr_num}'
                
                # Check if this is the target chromosome
                if chr_name == target_chr:
                    abs_zscore = abs(zscore)
                    if abs_zscore > max_abs_zscore:
                        max_abs_zscore = abs_zscore
                        best_aberration = {
                            'zscore': zscore,
                            'detected_mb': mbsize,
                            'location': location
                        }
            except (ValueError, IndexError):
                continue
        
        return best_aberration
    
    except Exception as e:
        logger.debug(f"Error parsing WC report {report_file}: {e}")
        return None


def parse_wcx_bed(bed_file: Path, target_chr: str) -> Optional[Dict]:
    """Parse WCX aberrations.bed to extract z-score and detected length for target chromosome"""
    
    if not bed_file.exists():
        return None
    
    try:
        # Read BED file
        df = pd.read_csv(bed_file, sep='\t')
        
        if len(df) == 0:
            return None
        
        # Normalize chromosome names (e.g., "chr1" or "1" -> "chr1")
        df['chr_normalized'] = df['chr'].apply(lambda x: f'chr{x}' if not str(x).startswith('chr') else str(x))
        
        # Filter by target chromosome
        target_df = df[df['chr_normalized'] == target_chr].copy()
        
        if len(target_df) == 0:
            return None
        
        # Calculate detected length and find max absolute z-score
        target_df['detected_mb'] = (target_df['end'] - target_df['start']) / 1_000_000
        target_df['abs_zscore'] = target_df['zscore'].abs()
        
        # Get aberration with max absolute z-score
        max_idx = target_df['abs_zscore'].idxmax()
        best_row = target_df.loc[max_idx]
        
        return {
            'zscore': best_row['zscore'],
            'detected_mb': best_row['detected_mb'],
            'location': f"{best_row['chr']}:{int(best_row['start'])}-{int(best_row['end'])}"
        }
    
    except Exception as e:
        logger.debug(f"Error parsing WCX BED {bed_file}: {e}")
        return None


def process_sample(sample_dir: Path, disease: str) -> List[Dict]:
    """Process a single sample directory and extract aberration data for ALL target chromosomes"""
    
    sample_name = sample_dir.name
    metadata = parse_sample_name(sample_name)
    
    if not metadata['disease']:
        metadata['disease'] = disease
    
    results_dir = sample_dir / 'results'
    if not results_dir.exists():
        return []
    
    records = []
    
    # Process for ALL disease chromosomes (One-vs-Rest)
    for target_disease, target_chr in DISEASE_CHROMOSOMES.items():
        # Process WC orig
        wc_orig_report = results_dir / 'wc_orig_report.txt'
        wc_orig_data = parse_wc_report(wc_orig_report, target_chr)
        # ALWAYS append data, even if no aberration detected (use 0 values)
        records.append({
            'sample_id': sample_name,
            'disease': metadata['disease'],
            'target_disease': target_disease,
            'ff': metadata['ff'],
            'deletion_length_mb': metadata['deletion_length_mb'],
            'gender': metadata['gender'],
            'mode': 'wc_orig',
            'target_chr': target_chr,
            'zscore': wc_orig_data['zscore'] if wc_orig_data else 0.0,
            'detected_mb': wc_orig_data['detected_mb'] if wc_orig_data else 0.0,
            'location': wc_orig_data['location'] if wc_orig_data else ''
        })
        
        # Process WC fetus
        wc_fetus_report = results_dir / 'wc_fetus_report.txt'
        wc_fetus_data = parse_wc_report(wc_fetus_report, target_chr)
        # ALWAYS append data, even if no aberration detected (use 0 values)
        records.append({
            'sample_id': sample_name,
            'disease': metadata['disease'],
            'target_disease': target_disease,
            'ff': metadata['ff'],
            'deletion_length_mb': metadata['deletion_length_mb'],
            'gender': metadata['gender'],
            'mode': 'wc_fetus',
            'target_chr': target_chr,
            'zscore': wc_fetus_data['zscore'] if wc_fetus_data else 0.0,
            'detected_mb': wc_fetus_data['detected_mb'] if wc_fetus_data else 0.0,
            'location': wc_fetus_data['location'] if wc_fetus_data else ''
        })
        
        # Process WCX orig
        wcx_orig_bed = results_dir / 'wcx_orig_aberrations.bed'
        wcx_orig_data = parse_wcx_bed(wcx_orig_bed, target_chr)
        # ALWAYS append data, even if no aberration detected (use 0 values)
        records.append({
            'sample_id': sample_name,
            'disease': metadata['disease'],
            'target_disease': target_disease,
            'ff': metadata['ff'],
            'deletion_length_mb': metadata['deletion_length_mb'],
            'gender': metadata['gender'],
            'mode': 'wcx_orig',
            'target_chr': target_chr,
            'zscore': wcx_orig_data['zscore'] if wcx_orig_data else 0.0,
            'detected_mb': wcx_orig_data['detected_mb'] if wcx_orig_data else 0.0,
            'location': wcx_orig_data['location'] if wcx_orig_data else ''
        })
        
        # Process WCX fetus
        wcx_fetus_bed = results_dir / 'wcx_fetus_aberrations.bed'
        wcx_fetus_data = parse_wcx_bed(wcx_fetus_bed, target_chr)
        # ALWAYS append data, even if no aberration detected (use 0 values)
        records.append({
            'sample_id': sample_name,
            'disease': metadata['disease'],
            'target_disease': target_disease,
            'ff': metadata['ff'],
            'deletion_length_mb': metadata['deletion_length_mb'],
            'gender': metadata['gender'],
            'mode': 'wcx_fetus',
            'target_chr': target_chr,
            'zscore': wcx_fetus_data['zscore'] if wcx_fetus_data else 0.0,
            'detected_mb': wcx_fetus_data['detected_mb'] if wcx_fetus_data else 0.0,
            'location': wcx_fetus_data['location'] if wcx_fetus_data else ''
        })
    
    return records


def main():
    parser = argparse.ArgumentParser(
        description="Collect aberration data from WC and WCX results"
    )
    parser.add_argument(
        '--data-dir',
        type=str,
        default='/data/md_validation/analysis_result',
        help='Base directory containing disease directories'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='/data/md_validation/roc_results/aberration_data.csv',
        help='Output CSV file'
    )
    parser.add_argument(
        '--diseases',
        type=str,
        nargs='+',
        default=['1p36', '2q33', 'CDC', 'DGS', 'Jacobsen', 'PWS', 'WBS', 'WHS'],
        help='List of diseases to process'
    )
    
    args = parser.parse_args()
    
    data_dir = Path(args.data_dir)
    output_file = Path(args.output)
    
    if not data_dir.exists():
        logger.error(f"Data directory not found: {data_dir}")
        return 1
    
    # Create output directory
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    logger.info("="*80)
    logger.info("Collecting aberration data from WC and WCX results")
    logger.info("="*80)
    logger.info(f"Data directory: {data_dir}")
    logger.info(f"Output file: {output_file}")
    logger.info(f"Diseases: {', '.join(args.diseases)}")
    logger.info("")
    
    all_records = []
    
    for disease in args.diseases:
        disease_dir = data_dir / disease
        
        if not disease_dir.exists():
            logger.warning(f"Disease directory not found: {disease_dir}")
            continue
        
        logger.info(f"Processing {disease}...")
        
        # Find all sample directories
        sample_dirs = [d for d in disease_dir.iterdir() if d.is_dir()]
        
        logger.info(f"  Found {len(sample_dirs)} sample directories")
        
        processed = 0
        for idx, sample_dir in enumerate(sample_dirs, 1):
            if idx % 500 == 0 or idx == 1 or idx == len(sample_dirs):
                logger.info(f"    Processing {idx}/{len(sample_dirs)}...")
            
            records = process_sample(sample_dir, disease)
            all_records.extend(records)
            
            if records:
                processed += 1
        
        logger.info(f"  Successfully processed {processed} samples")
        logger.info("")
    
    # Create DataFrame
    if not all_records:
        logger.error("No aberration data collected!")
        return 1
    
    df = pd.DataFrame(all_records)
    
    logger.info("="*80)
    logger.info("Summary")
    logger.info("="*80)
    logger.info(f"Total records: {len(df)}")
    logger.info(f"Unique samples: {df['sample_id'].nunique()}")
    logger.info("")
    logger.info("Records by disease:")
    for disease in df['disease'].unique():
        count = len(df[df['disease'] == disease])
        logger.info(f"  {disease}: {count}")
    logger.info("")
    logger.info("Records by mode:")
    for mode in df['mode'].unique():
        count = len(df[df['mode'] == mode])
        logger.info(f"  {mode}: {count}")
    logger.info("")
    
    # Save to CSV
    df.to_csv(output_file, index=False)
    logger.info(f"Saved to: {output_file}")
    logger.info("")
    
    # Show sample data
    logger.info("Sample data (first 5 rows):")
    print(df.head().to_string())
    logger.info("")
    
    logger.info("="*80)
    logger.info("Complete!")
    logger.info("="*80)
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())

