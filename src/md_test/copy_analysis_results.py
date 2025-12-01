#!/usr/bin/env python3
"""
Copy Analysis Results (Essential Files Only)

Copies only essential analysis results (JSON, report.txt, aberrations.bed)
to a compact directory structure for dashboard and analysis.

Usage:
    # Single disease
    python copy_analysis_results.py \
        --source /data/md_validation/1p36 \
        --dest /data/md_validation/analysis_result/1p36

    # Multiple diseases
    python copy_analysis_results.py \
        --source-dirs /data/md_validation/1p36,/data/md_validation/2q33 \
        --dest /data/md_validation/analysis_result \
        --cleanup
"""

import argparse
import logging
import shutil
from pathlib import Path
from typing import List, Dict
import json

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)


def get_essential_files(sample_dir: Path, sample_name: str) -> List[Dict]:
    """Get list of essential files to copy"""
    
    files_to_copy = []
    
    # 1. JSON metadata
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
            'relative_path': f"results/wc_orig_report.txt",
            'type': 'wc_orig'
        })
    
    # 3. WC fetus report
    wc_fetus_report = sample_dir / "Output_WC" / "fetus" / f"{sample_name}.wc.fetus.report.txt"
    if wc_fetus_report.exists():
        files_to_copy.append({
            'source': wc_fetus_report,
            'relative_path': f"results/wc_fetus_report.txt",
            'type': 'wc_fetus'
        })
    
    # 4. WCX orig aberrations
    wcx_orig_bed = sample_dir / "Output_WCX" / "orig" / f"{sample_name}.wcx.orig_aberrations.bed"
    if wcx_orig_bed.exists():
        files_to_copy.append({
            'source': wcx_orig_bed,
            'relative_path': f"results/wcx_orig_aberrations.bed",
            'type': 'wcx_orig'
        })
    
    # 5. WCX fetus aberrations
    wcx_fetus_bed = sample_dir / "Output_WCX" / "fetus" / f"{sample_name}.wcx.fetus_aberrations.bed"
    if wcx_fetus_bed.exists():
        files_to_copy.append({
            'source': wcx_fetus_bed,
            'relative_path': f"results/wcx_fetus_aberrations.bed",
            'type': 'wcx_fetus'
        })
    
    return files_to_copy


def copy_sample_results(
    source_sample_dir: Path,
    dest_sample_dir: Path,
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
    files_to_copy = get_essential_files(source_sample_dir, sample_name)
    
    if not files_to_copy:
        logger.warning(f"No essential files found for {sample_name}")
        return stats
    
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
            logger.debug(f"  [DRY-RUN] Would copy: {relative_path}")
            stats['copied_files'] += 1
        else:
            try:
                # Create parent directory if needed
                dest_file.parent.mkdir(parents=True, exist_ok=True)
                
                # Copy file
                shutil.copy2(source_file, dest_file)
                
                # Update stats
                stats['copied_files'] += 1
                stats['copied_bytes'] += source_file.stat().st_size
                
                logger.debug(f"  Copied: {relative_path}")
            except Exception as e:
                logger.error(f"  Failed to copy {relative_path}: {e}")
                stats['missing_files'].append(str(relative_path))
    
    return stats


def cleanup_source_directory(
    sample_dir: Path,
    sample_name: str,
    dry_run: bool = False
) -> Dict:
    """Remove non-essential files from source directory"""
    
    stats = {
        'deleted_files': 0,
        'deleted_bytes': 0
    }
    
    # Files/directories to delete
    patterns_to_delete = [
        # BAM files
        f"{sample_name}.proper_paired.bam",
        f"{sample_name}.proper_paired.bam.bai",
        f"{sample_name}.of_orig.bam",
        f"{sample_name}.of_orig.bam.bai",
        f"{sample_name}.of_fetus.bam",
        f"{sample_name}.of_fetus.bam.bai",
        # NPZ files
        "*.npz",
        # Plot directories
        "*.plots",
        # Intermediate files
        "*.marker"
    ]
    
    import glob
    
    for pattern in patterns_to_delete:
        if '*' in pattern:
            # Glob pattern
            for file_path in sample_dir.rglob(pattern):
                if file_path.exists():
                    try:
                        size = 0
                        if file_path.is_file():
                            size = file_path.stat().st_size
                        elif file_path.is_dir():
                            size = sum(f.stat().st_size for f in file_path.rglob('*') if f.is_file())
                        
                        if dry_run:
                            logger.debug(f"  [DRY-RUN] Would delete: {file_path.name} ({size/1e6:.1f}MB)")
                        else:
                            if file_path.is_file():
                                file_path.unlink()
                            elif file_path.is_dir():
                                shutil.rmtree(file_path)
                            logger.debug(f"  Deleted: {file_path.name} ({size/1e6:.1f}MB)")
                        
                        stats['deleted_files'] += 1
                        stats['deleted_bytes'] += size
                    except Exception as e:
                        logger.error(f"  Failed to delete {file_path}: {e}")
        else:
            # Exact filename
            file_path = sample_dir / pattern
            if file_path.exists():
                try:
                    size = file_path.stat().st_size
                    
                    if dry_run:
                        logger.debug(f"  [DRY-RUN] Would delete: {pattern} ({size/1e6:.1f}MB)")
                    else:
                        file_path.unlink()
                        logger.debug(f"  Deleted: {pattern} ({size/1e6:.1f}MB)")
                    
                    stats['deleted_files'] += 1
                    stats['deleted_bytes'] += size
                except Exception as e:
                    logger.error(f"  Failed to delete {pattern}: {e}")
    
    return stats


def process_disease_directory(
    source_dir: Path,
    dest_dir: Path,
    disease_name: str = None,
    cleanup: bool = False,
    dry_run: bool = False
) -> Dict:
    """Process a single disease directory"""
    
    if disease_name is None:
        disease_name = source_dir.name
    
    logger.info(f"Processing disease: {disease_name}")
    logger.info(f"  Source: {source_dir}")
    logger.info(f"  Destination: {dest_dir}")
    
    # Find all sample directories (containing JSON files)
    sample_dirs = []
    for item in source_dir.iterdir():
        if item.is_dir():
            json_files = list(item.glob("*.json"))
            if json_files:
                sample_dirs.append(item)
    
    logger.info(f"  Found {len(sample_dirs)} sample directories")
    
    if len(sample_dirs) == 0:
        logger.warning(f"  No sample directories found in {source_dir}")
        return None
    
    # Create destination directory
    if not dry_run:
        dest_dir.mkdir(parents=True, exist_ok=True)
    
    # Process each sample
    total_stats = {
        'total_samples': len(sample_dirs),
        'processed_samples': 0,
        'copied_files': 0,
        'copied_bytes': 0,
        'deleted_files': 0,
        'deleted_bytes': 0
    }
    
    for idx, sample_dir in enumerate(sample_dirs, 1):
        sample_name = sample_dir.name
        
        if idx % 100 == 0 or idx == 1 or idx == len(sample_dirs):
            logger.info(f"  Processing sample {idx}/{len(sample_dirs)}: {sample_name}")
        
        dest_sample_dir = dest_dir / sample_name
        
        # Copy essential files
        copy_stats = copy_sample_results(
            sample_dir, dest_sample_dir, sample_name, dry_run
        )
        
        total_stats['copied_files'] += copy_stats['copied_files']
        total_stats['copied_bytes'] += copy_stats['copied_bytes']
        
        if copy_stats['copied_files'] > 0:
            total_stats['processed_samples'] += 1
        
        # Cleanup if requested
        if cleanup and copy_stats['copied_files'] > 0:
            cleanup_stats = cleanup_source_directory(
                sample_dir, sample_name, dry_run
            )
            total_stats['deleted_files'] += cleanup_stats['deleted_files']
            total_stats['deleted_bytes'] += cleanup_stats['deleted_bytes']
    
    return total_stats


def main():
    parser = argparse.ArgumentParser(
        description="Copy essential analysis results to compact directory"
    )
    parser.add_argument(
        '--source',
        type=str,
        help='Source directory (single disease, e.g., /data/md_validation/1p36)'
    )
    parser.add_argument(
        '--source-dirs',
        type=str,
        help='Comma-separated source directories (multiple diseases)'
    )
    parser.add_argument(
        '--dest',
        type=str,
        required=True,
        help='Destination directory (e.g., /data/md_validation/analysis_result or /data/md_validation/analysis_result/1p36)'
    )
    parser.add_argument(
        '--zscore-tsv',
        type=str,
        help='Z-score extraction TSV file to copy (optional, e.g., ~/ken-nipt/analysis/md_validation/zscore/zscore_extraction_1p36.tsv)'
    )
    parser.add_argument(
        '--zscore-dir',
        type=str,
        help='Directory containing z-score TSV files (optional, e.g., ~/ken-nipt/analysis/md_validation/zscore)'
    )
    parser.add_argument(
        '--cleanup',
        action='store_true',
        help='Remove non-essential files from source after successful copy'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be done without actually copying/deleting'
    )
    
    args = parser.parse_args()
    
    if not args.source and not args.source_dirs:
        logger.error("Either --source or --source-dirs must be specified")
        return 1
    
    dest_root = Path(args.dest)
    
    # Determine source directories
    source_dirs = []
    if args.source:
        source_dirs.append(Path(args.source))
    
    if args.source_dirs:
        for source_path in args.source_dirs.split(','):
            source_dirs.append(Path(source_path.strip()))
    
    # Validate source directories
    for source_dir in source_dirs:
        if not source_dir.exists():
            logger.error(f"Source directory not found: {source_dir}")
            return 1
    
    if args.dry_run:
        logger.info("="*80)
        logger.info("DRY-RUN MODE (no files will be copied or deleted)")
        logger.info("="*80)
    
    if args.cleanup and not args.dry_run:
        logger.warning("="*80)
        logger.warning("CLEANUP MODE: Non-essential files will be DELETED from source!")
        logger.warning("="*80)
        response = input("Are you sure? Type 'yes' to continue: ")
        if response.lower() != 'yes':
            logger.info("Aborted by user")
            return 0
    
    # Copy z-score TSV files if provided
    zscore_stats = {'copied_files': 0, 'copied_bytes': 0}
    
    if args.zscore_dir:
        zscore_src_dir = Path(args.zscore_dir)
        if zscore_src_dir.exists():
            zscore_dest_dir = dest_root / 'zscore_data'
            
            if not args.dry_run:
                zscore_dest_dir.mkdir(parents=True, exist_ok=True)
            
            logger.info("Copying z-score TSV files...")
            logger.info(f"  Source: {zscore_src_dir}")
            logger.info(f"  Destination: {zscore_dest_dir}")
            
            # Copy all TSV files
            tsv_files = list(zscore_src_dir.glob("*.tsv"))
            for tsv_file in tsv_files:
                dest_tsv = zscore_dest_dir / tsv_file.name
                
                if args.dry_run:
                    logger.info(f"  [DRY-RUN] Would copy: {tsv_file.name}")
                    zscore_stats['copied_files'] += 1
                else:
                    try:
                        shutil.copy2(tsv_file, dest_tsv)
                        file_size = tsv_file.stat().st_size
                        zscore_stats['copied_files'] += 1
                        zscore_stats['copied_bytes'] += file_size
                        logger.info(f"  Copied: {tsv_file.name} ({file_size/1e6:.1f} MB)")
                    except Exception as e:
                        logger.error(f"  Failed to copy {tsv_file.name}: {e}")
        else:
            logger.warning(f"Z-score directory not found: {zscore_src_dir}")
    
    elif args.zscore_tsv:
        zscore_src_file = Path(args.zscore_tsv)
        if zscore_src_file.exists():
            zscore_dest_dir = dest_root / 'zscore_data'
            
            if not args.dry_run:
                zscore_dest_dir.mkdir(parents=True, exist_ok=True)
            
            dest_tsv = zscore_dest_dir / zscore_src_file.name
            
            if args.dry_run:
                logger.info(f"[DRY-RUN] Would copy z-score TSV: {zscore_src_file.name}")
                zscore_stats['copied_files'] += 1
            else:
                try:
                    shutil.copy2(zscore_src_file, dest_tsv)
                    file_size = zscore_src_file.stat().st_size
                    zscore_stats['copied_files'] += 1
                    zscore_stats['copied_bytes'] += file_size
                    logger.info(f"Copied z-score TSV: {zscore_src_file.name} ({file_size/1e6:.1f} MB)")
                except Exception as e:
                    logger.error(f"Failed to copy z-score TSV: {e}")
        else:
            logger.warning(f"Z-score TSV file not found: {zscore_src_file}")
    
    # Process each disease directory
    all_stats = []
    
    for source_dir in source_dirs:
        disease_name = source_dir.name
        
        # Determine destination
        if len(source_dirs) == 1 and dest_root.name == disease_name:
            # Destination already includes disease name
            dest_dir = dest_root
        else:
            # Destination is parent, add disease name
            dest_dir = dest_root / disease_name
        
        stats = process_disease_directory(
            source_dir, dest_dir, disease_name, args.cleanup, args.dry_run
        )
        
        if stats:
            stats['disease'] = disease_name
            all_stats.append(stats)
    
    # Print summary
    logger.info("")
    logger.info("="*80)
    logger.info("SUMMARY")
    logger.info("="*80)
    
    for stats in all_stats:
        logger.info(f"\n{stats['disease']}:")
        logger.info(f"  Samples processed: {stats['processed_samples']}/{stats['total_samples']}")
        logger.info(f"  Files copied: {stats['copied_files']}")
        logger.info(f"  Total copied: {stats['copied_bytes']/1e6:.1f} MB")
        
        if args.cleanup:
            logger.info(f"  Files deleted: {stats['deleted_files']}")
            logger.info(f"  Space freed: {stats['deleted_bytes']/1e9:.2f} GB")
    
    # Overall summary
    total_copied_mb = sum(s['copied_bytes'] for s in all_stats) / 1e6
    total_freed_gb = sum(s['deleted_bytes'] for s in all_stats) / 1e9
    
    logger.info("")
    logger.info("Overall:")
    logger.info(f"  Sample results copied: {total_copied_mb:.1f} MB")
    if zscore_stats['copied_files'] > 0:
        logger.info(f"  Z-score TSV files copied: {zscore_stats['copied_files']} files ({zscore_stats['copied_bytes']/1e6:.1f} MB)")
        total_copied_mb += zscore_stats['copied_bytes'] / 1e6
    logger.info(f"  Total copied: {total_copied_mb:.1f} MB")
    if args.cleanup:
        logger.info(f"  Total space freed: {total_freed_gb:.2f} GB")
    
    logger.info("")
    logger.info("="*80)
    logger.info("Complete!")
    logger.info("="*80)
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())

