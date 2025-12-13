#!/usr/bin/env python3
"""
Copy Analysis Results for Individual Samples (Pattern-based)

Copies analysis results for samples matching a specific pattern
from specified source directory to /data/md_validation/analysis_result.

Source locations:
- /data/md_validation: 1p36, 2q33, CDC, DGS
- ~/ken-nipt/analysis/md_validation: Jacobsen, PWS, WBS, WHS

Usage:
    # Copy from /data/md_validation (1p36, 2q33, CDC, DGS)
    python copy_analysis_results_individual.py \
        --pattern "*_0Mb_*" \
        --source-base /data/md_validation \
        --dest /data/md_validation/analysis_result

    # Copy from ~/ken-nipt/analysis/md_validation (Jacobsen, PWS, WBS, WHS)
    python copy_analysis_results_individual.py \
        --pattern "*_0Mb_*" \
        --source-base ~/ken-nipt/analysis/md_validation \
        --dest /data/md_validation/analysis_result

    # Dry-run to see what would be copied
    python copy_analysis_results_individual.py \
        --pattern "*_0Mb_*" \
        --source-base /data/md_validation \
        --dry-run
"""

import argparse
import logging
import shutil
from pathlib import Path
from typing import List, Dict, Set
import fnmatch

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)


def find_source_directories(source_base: str = None) -> List[Path]:
    """Find all source directories containing analysis results"""
    
    source_dirs = []
    
    if source_base:
        # Use specified source base directory
        base_dir = Path(source_base)
        if base_dir.exists():
            for disease_dir in base_dir.iterdir():
                if disease_dir.is_dir() and not disease_dir.name.startswith('.'):
                    # Skip analysis_result directory itself
                    if disease_dir.name != 'analysis_result':
                        source_dirs.append(disease_dir)
        else:
            logger.error(f"Source base directory not found: {source_base}")
    else:
        # Default: search both locations
        # 1. ~/ken-nipt/analysis/md_validation/* (Jacobsen, PWS, WBS, WHS)
        home_analysis = Path.home() / "ken-nipt" / "analysis" / "md_validation"
        if home_analysis.exists():
            for disease_dir in home_analysis.iterdir():
                if disease_dir.is_dir() and not disease_dir.name.startswith('.'):
                    source_dirs.append(disease_dir)
        
        # 2. /data/md_validation/* (1p36, 2q33, CDC, DGS)
        data_validation = Path("/data/md_validation")
        if data_validation.exists():
            for disease_dir in data_validation.iterdir():
                if disease_dir.is_dir() and not disease_dir.name.startswith('.'):
                    # Skip analysis_result directory itself
                    if disease_dir.name != 'analysis_result':
                        source_dirs.append(disease_dir)
    
    return source_dirs


def find_matching_samples(source_dirs: List[Path], pattern: str) -> Dict[str, Path]:
    """Find all samples matching the pattern"""
    
    matching_samples = {}
    
    for source_dir in source_dirs:
        disease_name = source_dir.name
        
        logger.info(f"Searching in {disease_name}: {source_dir}")
        
        # Find all sample directories matching pattern
        for item in source_dir.iterdir():
            if item.is_dir():
                sample_name = item.name
                
                # Check if sample name matches pattern
                if fnmatch.fnmatch(sample_name, pattern):
                    # Check if essential files exist
                    wc_dir = item / "Output_WC"
                    wcx_dir = item / "Output_WCX"
                    
                    if wc_dir.exists() or wcx_dir.exists():
                        key = f"{disease_name}/{sample_name}"
                        matching_samples[key] = item
                        logger.debug(f"  Found matching sample: {sample_name}")
    
    return matching_samples


def get_files_to_copy(sample_dir: Path, sample_name: str) -> List[Dict]:
    """Get list of essential files to copy for a sample"""
    
    files_to_copy = []
    
    # 1. JSON metadata (optional)
    json_file = sample_dir / f"{sample_name}.json"
    if json_file.exists():
        files_to_copy.append({
            'source': json_file,
            'relative_path': f"{sample_name}.json",
            'type': 'metadata'
        })
    
    # 2. WC orig report
    wc_orig_report = sample_dir / "Output_WC" / "orig" / f"{sample_name}.wc.orig.report.txt"
    if wc_orig_report.exists():
        files_to_copy.append({
            'source': wc_orig_report,
            'relative_path': "results/wc_orig_report.txt",
            'type': 'wc_orig'
        })
    
    # 3. WC fetus report
    wc_fetus_report = sample_dir / "Output_WC" / "fetus" / f"{sample_name}.wc.fetus.report.txt"
    if wc_fetus_report.exists():
        files_to_copy.append({
            'source': wc_fetus_report,
            'relative_path': "results/wc_fetus_report.txt",
            'type': 'wc_fetus'
        })
    
    # 4. WCX orig aberrations
    wcx_orig_bed = sample_dir / "Output_WCX" / "orig" / f"{sample_name}.wcx.orig_aberrations.bed"
    if wcx_orig_bed.exists():
        files_to_copy.append({
            'source': wcx_orig_bed,
            'relative_path': "results/wcx_orig_aberrations.bed",
            'type': 'wcx_orig'
        })
    
    # 5. WCX fetus aberrations
    wcx_fetus_bed = sample_dir / "Output_WCX" / "fetus" / f"{sample_name}.wcx.fetus_aberrations.bed"
    if wcx_fetus_bed.exists():
        files_to_copy.append({
            'source': wcx_fetus_bed,
            'relative_path': "results/wcx_fetus_aberrations.bed",
            'type': 'wcx_fetus'
        })
    
    return files_to_copy


def copy_sample(
    source_sample_dir: Path,
    dest_disease_dir: Path,
    sample_name: str,
    dry_run: bool = False
) -> Dict:
    """Copy essential files for a single sample"""
    
    stats = {
        'copied_files': 0,
        'copied_bytes': 0,
        'missing_files': []
    }
    
    # Get files to copy
    files_to_copy = get_files_to_copy(source_sample_dir, sample_name)
    
    if not files_to_copy:
        logger.warning(f"  No essential files found for {sample_name}")
        return stats
    
    dest_sample_dir = dest_disease_dir / sample_name
    
    # Create destination directory
    if not dry_run:
        dest_sample_dir.mkdir(parents=True, exist_ok=True)
        results_dir = dest_sample_dir / "results"
        results_dir.mkdir(exist_ok=True)
    
    # Copy files
    for file_info in files_to_copy:
        source_file = file_info['source']
        relative_path = file_info['relative_path']
        dest_file = dest_sample_dir / relative_path
        
        if dry_run:
            logger.debug(f"    [DRY-RUN] Would copy: {file_info['type']}")
            stats['copied_files'] += 1
        else:
            try:
                # Create parent directory if needed
                dest_file.parent.mkdir(parents=True, exist_ok=True)
                
                # Copy file
                shutil.copy2(source_file, dest_file)
                
                # Update stats
                file_size = source_file.stat().st_size
                stats['copied_files'] += 1
                stats['copied_bytes'] += file_size
                
                logger.debug(f"    Copied: {file_info['type']} ({file_size/1e3:.1f} KB)")
            except Exception as e:
                logger.error(f"    Failed to copy {relative_path}: {e}")
                stats['missing_files'].append(str(relative_path))
    
    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Copy analysis results for samples matching a pattern"
    )
    parser.add_argument(
        '--pattern',
        type=str,
        required=True,
        help='Pattern to match sample names (e.g., "*_0Mb_*")'
    )
    parser.add_argument(
        '--source-base',
        type=str,
        help='Source base directory (e.g., /data/md_validation or ~/ken-nipt/analysis/md_validation). If not specified, searches both locations.'
    )
    parser.add_argument(
        '--dest',
        type=str,
        default='/data/md_validation/analysis_result',
        help='Destination directory (default: /data/md_validation/analysis_result)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be done without actually copying'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    dest_root = Path(args.dest)
    pattern = args.pattern
    
    if args.dry_run:
        logger.info("="*80)
        logger.info("DRY-RUN MODE (no files will be copied)")
        logger.info("="*80)
    
    logger.info(f"Pattern: {pattern}")
    if args.source_base:
        logger.info(f"Source base: {args.source_base}")
    else:
        logger.info(f"Source base: Both /data/md_validation and ~/ken-nipt/analysis/md_validation")
    logger.info(f"Destination: {dest_root}")
    logger.info("")
    
    # Find source directories
    logger.info("Finding source directories...")
    source_dirs = find_source_directories(args.source_base)
    
    if not source_dirs:
        logger.error("No source directories found!")
        return 1
    
    logger.info(f"Found {len(source_dirs)} disease directories:")
    for source_dir in source_dirs:
        logger.info(f"  - {source_dir}")
    logger.info("")
    
    # Find matching samples
    logger.info("Finding samples matching pattern...")
    matching_samples = find_matching_samples(source_dirs, pattern)
    
    if not matching_samples:
        logger.warning(f"No samples found matching pattern: {pattern}")
        return 0
    
    logger.info(f"Found {len(matching_samples)} matching samples")
    logger.info("")
    
    # Group samples by disease
    samples_by_disease = {}
    for key, sample_dir in matching_samples.items():
        disease_name = key.split('/')[0]
        sample_name = key.split('/')[1]
        
        if disease_name not in samples_by_disease:
            samples_by_disease[disease_name] = []
        
        samples_by_disease[disease_name].append({
            'name': sample_name,
            'source_dir': sample_dir
        })
    
    # Process each disease
    total_stats = {
        'total_samples': len(matching_samples),
        'processed_samples': 0,
        'copied_files': 0,
        'copied_bytes': 0,
        'failed_samples': 0
    }
    
    disease_stats = {}
    
    for disease_name in sorted(samples_by_disease.keys()):
        samples = samples_by_disease[disease_name]
        
        logger.info(f"Processing {disease_name} ({len(samples)} samples)...")
        
        dest_disease_dir = dest_root / disease_name
        
        # Create disease directory
        if not args.dry_run:
            dest_disease_dir.mkdir(parents=True, exist_ok=True)
        
        disease_copy_stats = {
            'samples': len(samples),
            'copied_files': 0,
            'copied_bytes': 0,
            'failed': 0
        }
        
        for idx, sample_info in enumerate(samples, 1):
            sample_name = sample_info['name']
            source_sample_dir = sample_info['source_dir']
            
            logger.info(f"  [{idx}/{len(samples)}] {sample_name}")
            
            # Copy sample
            copy_stats = copy_sample(
                source_sample_dir,
                dest_disease_dir,
                sample_name,
                args.dry_run
            )
            
            if copy_stats['copied_files'] > 0:
                total_stats['processed_samples'] += 1
                disease_copy_stats['copied_files'] += copy_stats['copied_files']
                disease_copy_stats['copied_bytes'] += copy_stats['copied_bytes']
                
                logger.info(f"    ✓ Copied {copy_stats['copied_files']} files "
                          f"({copy_stats['copied_bytes']/1e6:.2f} MB)")
            else:
                total_stats['failed_samples'] += 1
                disease_copy_stats['failed'] += 1
                logger.warning(f"    ✗ No files copied")
        
        disease_stats[disease_name] = disease_copy_stats
        total_stats['copied_files'] += disease_copy_stats['copied_files']
        total_stats['copied_bytes'] += disease_copy_stats['copied_bytes']
        
        logger.info("")
    
    # Print summary
    logger.info("="*80)
    logger.info("SUMMARY")
    logger.info("="*80)
    logger.info(f"Pattern: {pattern}")
    logger.info(f"Total samples found: {total_stats['total_samples']}")
    logger.info(f"Successfully processed: {total_stats['processed_samples']}")
    logger.info(f"Failed: {total_stats['failed_samples']}")
    logger.info("")
    
    for disease_name in sorted(disease_stats.keys()):
        stats = disease_stats[disease_name]
        logger.info(f"{disease_name}:")
        logger.info(f"  Samples: {stats['samples']}")
        logger.info(f"  Files copied: {stats['copied_files']}")
        logger.info(f"  Size: {stats['copied_bytes']/1e6:.2f} MB")
        if stats['failed'] > 0:
            logger.info(f"  Failed: {stats['failed']}")
    
    logger.info("")
    logger.info("Overall:")
    logger.info(f"  Files copied: {total_stats['copied_files']}")
    logger.info(f"  Total size: {total_stats['copied_bytes']/1e6:.2f} MB")
    logger.info("")
    logger.info("="*80)
    logger.info("Complete!")
    logger.info("="*80)
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())

