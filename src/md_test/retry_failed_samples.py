#!/usr/bin/env python3
"""
Retry failed samples from run_parallel_v3.py output.

This script:
1. Finds failed samples (no BAM file or error in log)
2. Extracts task information from original run
3. Re-runs only failed samples

Usage:
    python3 retry_failed_samples.py \
        --output_dir batch_output \
        --script make_artificial.sh \
        --max_workers 16
"""

import argparse
import subprocess
import sys
import logging
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional
import json

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)7s] | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)


def find_failed_samples(output_dir: Path, mom_tsv: Path, female_tsv: Path, male_tsv: Path) -> List[Dict]:
    """Find samples that failed by checking for missing BAM files or error logs."""
    failed_samples = []
    
    logger.info(f"Scanning {output_dir} for failed samples...")
    
    sample_dirs = [d for d in output_dir.iterdir() if d.is_dir()]
    total = len(sample_dirs)
    
    for idx, sample_dir in enumerate(sample_dirs, 1):
        if idx % 100 == 0:
            logger.info(f"Checked {idx}/{total} samples...")
        
        sample_id = sample_dir.name
        
        # Check if BAM file exists
        bam_file = sample_dir / f"{sample_id}.proper_paired.bam"
        json_file = sample_dir / f"{sample_id}.json"
        log_file = sample_dir / "generation.log"
        
        # Sample failed if:
        # 1. BAM file doesn't exist, OR
        # 2. JSON file doesn't exist (not completed), OR
        # 3. Log file contains error
        failed = False
        
        if not bam_file.exists():
            failed = True
        elif not json_file.exists():
            failed = True
        elif log_file.exists():
            # Check log for errors
            try:
                with log_file.open() as f:
                    content = f.read()
                    if "ERROR" in content.upper() or "FAILED" in content.upper() or "error" in content.lower():
                        # Check return code in log
                        if "Return code" in content or "exit code" in content.lower():
                            failed = True
            except Exception:
                pass
        
        if failed:
            # Try to extract task info from existing files
            task_info = extract_task_info(sample_dir, sample_id, mom_tsv, female_tsv, male_tsv)
            if task_info:
                failed_samples.append(task_info)
    
    logger.info(f"Found {len(failed_samples)} failed samples")
    return failed_samples


def load_ff_from_tsv(tsv_path: Path) -> Dict[str, float]:
    """Load FF values from TSV file.
    
    Expected format:
    Work_Dir | Sample_ID | Gender | FF_Method | FF_Value | ... | BAM_Path (last column)
    Column 2 = Sample_ID (index 1)
    Column 5 = FF_Value (index 4)
    """
    ff_map = {}
    if not tsv_path.exists():
        return ff_map
    
    with tsv_path.open() as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue
            sample_id = parts[1]  # Column 2: Sample_ID
            try:
                ff_value = float(parts[4])  # Column 5: FF_Value
                ff_map[sample_id] = ff_value
            except (ValueError, IndexError):
                continue
    
    return ff_map


def is_valid_ff_combination(mom_ff: float, fetus_ff: float, target_ff: float) -> bool:
    """Check if target FF is achievable with given mom and fetus FF values.
    
    Beta = (fT - fA) / (fB - fA)
    Valid if: 0.0 <= beta <= 1.0
    """
    fA = mom_ff / 100.0
    fB = fetus_ff / 100.0
    fT = target_ff / 100.0
    
    if fB == fA:
        return False  # Cannot mix if FF values are the same
    
    beta = (fT - fA) / (fB - fA)
    return 0.0 <= beta <= 1.0


def extract_task_info(sample_dir: Path, sample_id: str, mom_tsv: Path, female_tsv: Path, male_tsv: Path) -> Optional[Dict]:
    """Extract task information from sample directory."""
    # Parse sample_id: {mom_idx}_{fetus_idx}_{disease}_FF{ff}_{reads}M_{del}Mb_{gender}
    try:
        parts = sample_id.split('_')
        if len(parts) < 7:
            return None
        
        mom_idx = int(parts[0])
        fetus_idx = int(parts[1])
        
        # Find disease name (may contain underscores)
        # Format: {mom_idx}_{fetus_idx}_{disease}_FF{ff}_{reads}M_{del}Mb_{gender}
        # Find FF, reads, del, gender
        ff_idx = None
        reads_idx = None
        del_idx = None
        gender_idx = None
        
        for i, part in enumerate(parts):
            if part.startswith('FF'):
                ff_idx = i
                ff_target = int(part[2:])
            elif part.endswith('M') and part[:-1].isdigit():
                reads_idx = i
                reads_m = int(part[:-1])
                reads = reads_m * 1_000_000
            elif part.endswith('Mb') and part[:-2].isdigit():
                del_idx = i
                del_size_mb = int(part[:-2])
            elif part in ['M', 'F']:
                gender_idx = i
                gender = part
        
        if None in [ff_idx, reads_idx, del_idx, gender_idx]:
            return None
        
        # Extract disease name (between fetus_idx and FF)
        disease_name = '_'.join(parts[2:ff_idx])
        
        # Try to find BAM files from log or JSON
        log_file = sample_dir / "generation.log"
        json_file = sample_dir / f"{sample_id}.json"
        
        mom_bam = None
        fetus_bam = None
        bed_file = None
        
        # Try to read from log file
        if log_file.exists():
            try:
                with log_file.open() as f:
                    content = f.read()
                    # Extract paths from command line
                    import re
                    mom_match = re.search(r'--mom_bam\s+(\S+)', content)
                    fetus_match = re.search(r'--fetus_bam\s+(\S+)', content)
                    bed_match = re.search(r'--md_bed\s+(\S+)', content)
                    
                    if mom_match:
                        mom_bam = Path(mom_match.group(1))
                    if fetus_match:
                        fetus_bam = Path(fetus_match.group(1))
                    if bed_match:
                        bed_file = Path(bed_match.group(1))
            except Exception:
                pass
        
        # Try to read from JSON file
        if json_file.exists():
            try:
                with json_file.open() as f:
                    data = json.load(f)
                    if 'source_samples' in data:
                        if 'mom' in data['source_samples'] and 'bam_path' in data['source_samples']['mom']:
                            mom_bam = Path(data['source_samples']['mom']['bam_path'])
                        if 'fetus' in data['source_samples'] and 'bam_path' in data['source_samples']['fetus']:
                            fetus_bam = Path(data['source_samples']['fetus']['bam_path'])
                    if 'deletion' in data and 'bed_file' in data['deletion']:
                        bed_file = Path(data['deletion']['bed_file'])
            except Exception:
                pass
        
        if mom_bam and fetus_bam and bed_file and mom_bam.exists() and fetus_bam.exists() and bed_file.exists():
            # Get FF values from TSV files
            mom_id = mom_bam.parent.name
            fetus_id = fetus_bam.parent.name
            
            # Load FF maps from TSV files
            mom_ff_map = load_ff_from_tsv(mom_tsv)
            female_ff_map = load_ff_from_tsv(female_tsv)
            male_ff_map = load_ff_from_tsv(male_tsv)
            
            # Get FF values
            mom_ff = mom_ff_map.get(mom_id, 0.0)
            
            # Try both female and male TSV files for fetus
            fetus_ff = female_ff_map.get(fetus_id, 0.0)
            if fetus_ff == 0.0:
                fetus_ff = male_ff_map.get(fetus_id, 0.0)
            
            # Fallback: try to read from existing ff_map.tsv if TSV lookup fails
            if mom_ff == 0.0 or fetus_ff == 0.0:
                ff_map_file = sample_dir / "ff_map.tsv"
                if ff_map_file.exists():
                    try:
                        with ff_map_file.open() as f:
                            for line in f:
                                parts = line.strip().split('\t')
                                if len(parts) >= 2:
                                    sample_id_from_map = parts[0]
                                    ff_value = float(parts[1])
                                    if mom_id in sample_id_from_map:
                                        mom_ff = ff_value
                                    elif fetus_id in sample_id_from_map:
                                        fetus_ff = ff_value
                    except Exception:
                        pass
            
            # Note: FF combination validation is now handled by make_artificial.sh
            # It will clamp beta to [0.0, 1.0] and proceed with closest achievable FF
            
            return {
                'sample_id': sample_id,
                'mom_bam': mom_bam,
                'fetus_bam': fetus_bam,
                'bed_file': bed_file,
                'ff_target': ff_target,
                'reads': reads,
                'mom_idx': mom_idx,
                'fetus_idx': fetus_idx,
                'gender': gender,
                'disease': disease_name,
                'del_size_mb': del_size_mb,
                'mom_ff': mom_ff,
                'fetus_ff': fetus_ff,
                'sample_dir': sample_dir
            }
    
    except Exception as e:
        logger.warning(f"Failed to extract task info for {sample_id}: {e}")
        return None
    
    return None


def retry_failed_samples(failed_samples: List[Dict], script: Path, output_dir: Path, max_workers: int):
    """Retry failed samples using similar logic to run_parallel_v3.py."""
    import sys
    import time
    import subprocess
    from datetime import datetime
    
    logger.info(f"\nRetrying {len(failed_samples)} failed samples...")
    
    # Simple retry logic without importing TaskManager
    # Process samples sequentially with max_workers parallel limit
    import queue
    import threading
    
    task_queue = queue.Queue()
    for idx, sample_info in enumerate(failed_samples, 1):
        task_queue.put((idx, sample_info))
    
    completed = 0
    failed = 0
    lock = threading.Lock()
    
    def worker():
        nonlocal completed, failed
        while True:
            try:
                idx, sample_info = task_queue.get_nowait()
            except queue.Empty:
                break
            
            sample_id = sample_info['sample_id']
            sample_dir = sample_info['sample_dir']
            mom_bam = sample_info['mom_bam']
            fetus_bam = sample_info['fetus_bam']
            bed_file = sample_info['bed_file']
            ff_target = sample_info['ff_target']
            reads = sample_info['reads']
            mom_idx = sample_info['mom_idx']
            fetus_idx = sample_info['fetus_idx']
            gender = sample_info['gender']
            
            # Get sample IDs
            mom_id = mom_bam.parent.name
            fetus_id = fetus_bam.parent.name
            
            # Create FF map file
            ff_map_file = sample_dir / "ff_map.tsv"
            with ff_map_file.open('w') as f:
                f.write(f"{mom_id}\t{sample_info['mom_ff']}\n")
                f.write(f"{fetus_id}\t{sample_info['fetus_ff']}\n")
            
            # Build command
            cmd = [
                str(script.resolve()),
                "--mom_bam", str(mom_bam),
                "--fetus_bam", str(fetus_bam),
                "--ff_map", str(ff_map_file),
                "--md_bed", str(bed_file),
                "--ff_target", str(ff_target),
                "--reads", str(reads),
                "--mom_idx", str(mom_idx),
                "--fetus_idx", str(fetus_idx),
                "--gender", gender,
                "--sample_id", sample_id,
                "--outdir", str(sample_dir)
            ]
            
            # Run command
            log_file = sample_dir / "retry.log"
            script_dir = script.parent
            
            try:
                with log_file.open('w') as lf:
                    lf.write(f"Retry command: {' '.join(cmd)}\n")
                    lf.write(f"Start: {datetime.now()}\n\n")
                
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    cwd=str(script_dir)
                )
                
                stdout, stderr = process.communicate()
                returncode = process.returncode
                
                with log_file.open('a') as lf:
                    lf.write("\n=== STDOUT ===\n")
                    lf.write(stdout)
                    lf.write("\n=== STDERR ===\n")
                    lf.write(stderr)
                    lf.write(f"\nReturn code: {returncode}\n")
                
                with lock:
                    if returncode == 0:
                        completed += 1
                        logger.info(f"[{idx:04d}] ✓ Retry succeeded: {sample_id}")
                    else:
                        failed += 1
                        logger.error(f"[{idx:04d}] ✗ Retry failed: {sample_id} (rc={returncode})")
            
            except Exception as e:
                with lock:
                    failed += 1
                    logger.error(f"[{idx:04d}] ✗ Retry error: {sample_id} - {e}")
            
            finally:
                task_queue.task_done()
    
    # Start workers
    threads = []
    for _ in range(max_workers):
        t = threading.Thread(target=worker)
        t.start()
        threads.append(t)
    
    # Wait for all threads
    for t in threads:
        t.join()
    
    logger.info(f"\n{'='*60}")
    logger.info("Retry Complete!")
    logger.info(f"  Succeeded: {completed}")
    logger.info(f"  Failed: {failed}")
    logger.info(f"  Total: {len(failed_samples)}")
    logger.info(f"{'='*60}\n")
    
    return failed == 0


def main():
    parser = argparse.ArgumentParser(
        description='Retry failed samples from run_parallel_v3.py output'
    )
    
    parser.add_argument('--output_dir', type=Path, required=True,
                       help='Output directory containing sample directories')
    parser.add_argument('--mom_tsv', type=Path, required=True,
                       help='Mom samples TSV file')
    parser.add_argument('--female_tsv', type=Path, required=True,
                       help='Female fetus samples TSV file')
    parser.add_argument('--male_tsv', type=Path, required=True,
                       help='Male fetus samples TSV file')
    parser.add_argument('--script', type=Path, default=None,
                       help='make_artificial.sh path (default: same directory as this script)')
    parser.add_argument('--max_workers', type=int, default=16,
                       help='Maximum concurrent workers')
    
    args = parser.parse_args()
    
    # Set default script path
    if args.script is None:
        script_dir = Path(__file__).parent
        args.script = script_dir / 'make_artificial.sh'
    
    args.script = args.script.resolve()
    
    if not args.script.exists():
        logger.error(f"Script not found: {args.script}")
        return 1
    
    if not args.output_dir.exists():
        logger.error(f"Output directory not found: {args.output_dir}")
        return 1
    
    # Validate TSV files
    if not args.mom_tsv.exists():
        logger.error(f"Mom TSV file not found: {args.mom_tsv}")
        return 1
    if not args.female_tsv.exists():
        logger.error(f"Female TSV file not found: {args.female_tsv}")
        return 1
    if not args.male_tsv.exists():
        logger.error(f"Male TSV file not found: {args.male_tsv}")
        return 1
    
    # Find failed samples
    failed_samples = find_failed_samples(args.output_dir, args.mom_tsv, args.female_tsv, args.male_tsv)
    
    if not failed_samples:
        logger.info("No failed samples found!")
        return 0
    
    logger.info(f"\nFailed samples ({len(failed_samples)}):")
    for sample in failed_samples[:10]:  # Show first 10
        logger.info(f"  - {sample['sample_id']}")
    if len(failed_samples) > 10:
        logger.info(f"  ... and {len(failed_samples) - 10} more")
    
    # Retry failed samples
    success = retry_failed_samples(failed_samples, args.script, args.output_dir, args.max_workers)
    
    return 0 if success else 1


if __name__ == '__main__':
    sys.exit(main())

