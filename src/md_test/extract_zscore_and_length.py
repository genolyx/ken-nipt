#!/usr/bin/env python3
"""
Extract Z-scores and Deletion Length from Artificial Sample Analysis Results

Scans artificial sample analysis directories and extracts:
- Expected deletion information from JSON
- Detected z-scores from WC/WCX report and BED files
- Detected region coordinates

Usage:
    python3 extract_zscore_and_length.py \
        -i /data/md_validation/1p36 \
        -o zscore_extraction.tsv
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import sys

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def parse_json_metadata(json_file: Path, sample_dir: Path) -> Optional[Dict]:
    """Parse sample metadata JSON file"""
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        # Extract sample name from directory name
        sample_name = sample_dir.name
        
        # Extract FF
        ff = data.get('target_parameters', {}).get('ff_target_percent', 0)
        
        # Extract deletion info
        deletion = data.get('deletion', {})
        if not deletion:
            logger.warning(f"No deletion info in {json_file}")
            return None
        
        deletion_chr = str(deletion.get('chromosome', '')).replace('chr', '')
        deletion_start = deletion.get('start', 0)
        deletion_end = deletion.get('end', 0)
        deletion_size_bp = deletion.get('size_bp', 0)
        deletion_length_mb = round(deletion_size_bp / 1_000_000) if deletion_size_bp > 0 else 0
        
        return {
            'sample_name': sample_name,
            'ff': ff,
            'deletion_length_mb': deletion_length_mb,
            'expected_deletion_chr': deletion_chr,
            'expected_deletion_start': deletion_start,
            'expected_deletion_end': deletion_end,
        }
    except Exception as e:
        logger.warning(f"Failed to parse JSON {json_file}: {e}")
        return None


def check_overlap(expected_start: int, expected_end: int, detected_start: int, detected_end: int) -> bool:
    """Check if two regions overlap"""
    return not (detected_end <= expected_start or detected_start >= expected_end)


def parse_wc_report(report_file: Path, expected_chr: str, expected_start: int, expected_end: int) -> Tuple[Optional[float], Optional[int], Optional[int]]:
    """Parse WC report.txt and extract z-score for region overlapping expected deletion
    
    Format:
    # Test results: #
    
    z-score effect  mbsize  location
    -51.74  -29.66  9.20    1:800000-10000000
    
    Returns: (zscore, detected_start, detected_end) or (None, None, None)
    """
    if not report_file.exists():
        return None, None, None
    
    try:
        with open(report_file, 'r') as f:
            lines = f.readlines()
        
        in_test_section = False
        best_zscore = None
        best_start = None
        best_end = None
        
        for line in lines:
            line_lower = line.lower()
            
            # Check for test results section
            if '# test results:' in line_lower or '#test results:' in line_lower:
                in_test_section = True
                continue
            
            if not in_test_section:
                continue
            
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Skip header line
            if 'z-score' in line_lower and ('effect' in line_lower or 'mbsize' in line_lower):
                continue
            
            # Parse data line: z-score  effect  mbsize  location
            fields = line.split()
            if len(fields) < 4:
                continue
            
            try:
                zscore = float(fields[0])
                # effect = float(fields[1])
                # mbsize = float(fields[2])
                location = fields[3]  # Format: chr:start-end
                
                # Parse location
                if ':' not in location or '-' not in location:
                    continue
                
                chr_part, range_part = location.split(':')
                chr_name = chr_part.replace('chr', '')
                start_str, end_str = range_part.split('-')
                start = int(start_str)
                end = int(end_str)
                
                # Check if chromosome matches and regions overlap
                if chr_name == expected_chr and check_overlap(expected_start, expected_end, start, end):
                    # Keep the most significant (most negative) z-score
                    if best_zscore is None or zscore < best_zscore:
                        best_zscore = zscore
                        best_start = start
                        best_end = end
                        
            except (ValueError, IndexError) as e:
                logger.debug(f"Failed to parse line: {line} - {e}")
                continue
        
        return best_zscore, best_start, best_end
        
    except Exception as e:
        logger.warning(f"Failed to parse WC report {report_file}: {e}")
        return None, None, None


def parse_wcx_bed(bed_file: Path, expected_chr: str, expected_start: int, expected_end: int) -> Tuple[Optional[float], Optional[int], Optional[int]]:
    """Parse WCX aberrations.bed and extract z-score for region overlapping expected deletion
    
    Format:
    chr     start   end     ratio   zscore  type
    1       800001  10000000        -0.5131 -49.417224811792835     loss
    
    Returns: (zscore, detected_start, detected_end) or (None, None, None)
    """
    if not bed_file.exists():
        return None, None, None
    
    try:
        best_zscore = None
        best_start = None
        best_end = None
        
        with open(bed_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('track'):
                    continue
                
                # Skip header line
                if line.startswith('chr') and 'start' in line and 'end' in line:
                    continue
                
                fields = line.split('\t')
                if len(fields) < 5:
                    # Try splitting by whitespace
                    fields = line.split()
                    if len(fields) < 5:
                        continue
                
                try:
                    chr_name = str(fields[0]).replace('chr', '')
                    start = int(fields[1])
                    end = int(fields[2])
                    # ratio = float(fields[3])
                    zscore = float(fields[4])
                    # type_str = fields[5] if len(fields) > 5 else ''
                    
                    # Check if chromosome matches and regions overlap
                    if chr_name == expected_chr and check_overlap(expected_start, expected_end, start, end):
                        # Keep the most significant (most negative) z-score
                        if best_zscore is None or zscore < best_zscore:
                            best_zscore = zscore
                            best_start = start
                            best_end = end
                            
                except (ValueError, IndexError) as e:
                    logger.debug(f"Failed to parse line: {line} - {e}")
                    continue
        
        return best_zscore, best_start, best_end
        
    except Exception as e:
        logger.warning(f"Failed to parse WCX BED {bed_file}: {e}")
        return None, None, None


def extract_sample_data(sample_dir: Path) -> Optional[Dict]:
    """Extract all data for a single sample"""
    
    # Find JSON file
    json_files = list(sample_dir.glob("*.json"))
    if not json_files:
        logger.debug(f"No JSON file found in {sample_dir}")
        return None
    
    json_file = json_files[0]
    
    # Parse JSON metadata
    metadata = parse_json_metadata(json_file, sample_dir)
    if not metadata:
        return None
    
    sample_name = metadata['sample_name']
    expected_chr = metadata['expected_deletion_chr']
    expected_start = metadata['expected_deletion_start']
    expected_end = metadata['expected_deletion_end']
    
    logger.debug(f"Processing sample: {sample_name}")
    
    # Initialize result - format basic metadata values
    result = {
        'sample_name': metadata['sample_name'],
        'ff': f"{metadata['ff']:.2f}",
        'deletion_length_mb': str(metadata['deletion_length_mb']),
        'expected_deletion_chr': metadata['expected_deletion_chr'],
        'expected_deletion_start': str(metadata['expected_deletion_start']),
        'expected_deletion_end': str(metadata['expected_deletion_end']),
    }
    
    # Extract WC orig data
    wc_orig_report = sample_dir / "Output_WC" / "orig" / f"{sample_name}.wc.orig.report.txt"
    wc_orig_zscore, wc_orig_start, wc_orig_end = parse_wc_report(
        wc_orig_report, expected_chr, expected_start, expected_end
    )
    result['WC_orig_zscore'] = f"{wc_orig_zscore:.2f}" if wc_orig_zscore is not None else ""
    if wc_orig_start is not None and wc_orig_end is not None:
        result['WC_orig_detected_length_mb'] = f"{(wc_orig_end - wc_orig_start) / 1_000_000:.2f}"
        result['WC_orig_detected_start'] = str(wc_orig_start)
        result['WC_orig_detected_end'] = str(wc_orig_end)
    else:
        result['WC_orig_detected_length_mb'] = ""
        result['WC_orig_detected_start'] = ""
        result['WC_orig_detected_end'] = ""
    
    # Extract WC fetus data
    wc_fetus_report = sample_dir / "Output_WC" / "fetus" / f"{sample_name}.wc.fetus.report.txt"
    wc_fetus_zscore, wc_fetus_start, wc_fetus_end = parse_wc_report(
        wc_fetus_report, expected_chr, expected_start, expected_end
    )
    result['WC_fetus_zscore'] = f"{wc_fetus_zscore:.2f}" if wc_fetus_zscore is not None else ""
    if wc_fetus_start is not None and wc_fetus_end is not None:
        result['WC_fetus_detected_length_mb'] = f"{(wc_fetus_end - wc_fetus_start) / 1_000_000:.2f}"
        result['WC_fetus_detected_start'] = str(wc_fetus_start)
        result['WC_fetus_detected_end'] = str(wc_fetus_end)
    else:
        result['WC_fetus_detected_length_mb'] = ""
        result['WC_fetus_detected_start'] = ""
        result['WC_fetus_detected_end'] = ""
    
    # Extract WCX orig data
    wcx_orig_bed = sample_dir / "Output_WCX" / "orig" / f"{sample_name}.wcx.orig_aberrations.bed"
    wcx_orig_zscore, wcx_orig_start, wcx_orig_end = parse_wcx_bed(
        wcx_orig_bed, expected_chr, expected_start, expected_end
    )
    result['WCX_orig_zscore'] = f"{wcx_orig_zscore:.2f}" if wcx_orig_zscore is not None else ""
    if wcx_orig_start is not None and wcx_orig_end is not None:
        result['WCX_orig_detected_length_mb'] = f"{(wcx_orig_end - wcx_orig_start) / 1_000_000:.2f}"
        result['WCX_orig_detected_start'] = str(wcx_orig_start)
        result['WCX_orig_detected_end'] = str(wcx_orig_end)
    else:
        result['WCX_orig_detected_length_mb'] = ""
        result['WCX_orig_detected_start'] = ""
        result['WCX_orig_detected_end'] = ""
    
    # Extract WCX fetus data
    wcx_fetus_bed = sample_dir / "Output_WCX" / "fetus" / f"{sample_name}.wcx.fetus_aberrations.bed"
    wcx_fetus_zscore, wcx_fetus_start, wcx_fetus_end = parse_wcx_bed(
        wcx_fetus_bed, expected_chr, expected_start, expected_end
    )
    result['WCX_fetus_zscore'] = f"{wcx_fetus_zscore:.2f}" if wcx_fetus_zscore is not None else ""
    if wcx_fetus_start is not None and wcx_fetus_end is not None:
        result['WCX_fetus_detected_length_mb'] = f"{(wcx_fetus_end - wcx_fetus_start) / 1_000_000:.2f}"
        result['WCX_fetus_detected_start'] = str(wcx_fetus_start)
        result['WCX_fetus_detected_end'] = str(wcx_fetus_end)
    else:
        result['WCX_fetus_detected_length_mb'] = ""
        result['WCX_fetus_detected_start'] = ""
        result['WCX_fetus_detected_end'] = ""
    
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Extract z-scores and deletion lengths from artificial sample analysis results"
    )
    parser.add_argument(
        '-i', '--input',
        type=Path,
        required=True,
        help='Input directory containing sample directories (e.g., /data/md_validation/1p36 or /home/ken/ken-nipt/analysis/md_validation/1p36)'
    )
    parser.add_argument(
        '-o', '--output',
        type=Path,
        required=True,
        help='Output TSV file path'
    )
    
    args = parser.parse_args()
    
    # Use input directory directly
    analysis_dir = args.input
    
    if not analysis_dir.exists():
        logger.error(f"Input directory not found: {analysis_dir}")
        return 1
    
    if not analysis_dir.is_dir():
        logger.error(f"Input path is not a directory: {analysis_dir}")
        return 1
    
    logger.info(f"Scanning input directory: {analysis_dir}")
    
    # Find all sample directories (containing JSON files)
    sample_dirs = []
    for item in analysis_dir.iterdir():
        if item.is_dir():
            json_files = list(item.glob("*.json"))
            if json_files:
                sample_dirs.append(item)
    
    logger.info(f"Found {len(sample_dirs)} sample directories")
    
    if len(sample_dirs) == 0:
        logger.error("No sample directories found!")
        return 1
    
    # Extract data from all samples
    results = []
    for sample_dir in sorted(sample_dirs):
        result = extract_sample_data(sample_dir)
        if result:
            results.append(result)
    
    logger.info(f"Successfully extracted data from {len(results)} samples")
    
    if len(results) == 0:
        logger.error("No data extracted!")
        return 1
    
    # Create output directory if needed
    output_file = args.output
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Write TSV file
    logger.info(f"Writing results to: {output_file}")
    
    # Define column order
    columns = [
        'sample_name',
        'ff',
        'deletion_length_mb',
        'expected_deletion_chr',
        'expected_deletion_start',
        'expected_deletion_end',
        'WC_orig_zscore',
        'WC_orig_detected_length_mb',
        'WC_orig_detected_start',
        'WC_orig_detected_end',
        'WC_fetus_zscore',
        'WC_fetus_detected_length_mb',
        'WC_fetus_detected_start',
        'WC_fetus_detected_end',
        'WCX_orig_zscore',
        'WCX_orig_detected_length_mb',
        'WCX_orig_detected_start',
        'WCX_orig_detected_end',
        'WCX_fetus_zscore',
        'WCX_fetus_detected_length_mb',
        'WCX_fetus_detected_start',
        'WCX_fetus_detected_end',
    ]
    
    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join(columns) + '\n')
        
        # Write data
        for result in results:
            row = []
            for col in columns:
                value = result.get(col, '')
                # Format numbers - all values are already formatted as strings in extract_sample_data
                row.append(str(value) if value != '' else "")
            f.write('\t'.join(row) + '\n')
    
    logger.info(f"✓ Successfully wrote {len(results)} rows to {output_file}")
    
    # Summary statistics
    logger.info("\nSummary Statistics:")
    logger.info(f"  Total samples: {len(results)}")
    
    wc_orig_detected = sum(1 for r in results if r.get('WC_orig_zscore'))
    wc_fetus_detected = sum(1 for r in results if r.get('WC_fetus_zscore'))
    wcx_orig_detected = sum(1 for r in results if r.get('WCX_orig_zscore'))
    wcx_fetus_detected = sum(1 for r in results if r.get('WCX_fetus_zscore'))
    
    logger.info(f"  WC_orig detected: {wc_orig_detected} ({wc_orig_detected/len(results)*100:.1f}%)")
    logger.info(f"  WC_fetus detected: {wc_fetus_detected} ({wc_fetus_detected/len(results)*100:.1f}%)")
    logger.info(f"  WCX_orig detected: {wcx_orig_detected} ({wcx_orig_detected/len(results)*100:.1f}%)")
    logger.info(f"  WCX_fetus detected: {wcx_fetus_detected} ({wcx_fetus_detected/len(results)*100:.1f}%)")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
