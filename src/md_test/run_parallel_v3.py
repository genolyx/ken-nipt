#!/usr/bin/env python3
"""
Parallel Artificial Sample Generator v3
Combines features from v2 (artificial sample generation) and run_batch_dev.py (slot management)

Key Features:
- Max concurrent sample limit with slot management
- One completes → next starts immediately
- Status tracking: Init → Running → Completed/Failed
- Detailed logging and summary generation
"""

import argparse
import subprocess
import sys
import time
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple
import re

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)7s] | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)


def load_bams_from_tsv(tsv_path: Path) -> List[Path]:
    """Load BAM paths from TSV file.
    
    Expected format:
    Work_Dir | Sample_ID | Gender | FF_Method | FF_Value | ... | BAM_Path (last column)
    """
    bams = []
    if not tsv_path.exists():
        return bams
    
    with tsv_path.open() as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue
            bam_path = Path(parts[-1])
            if bam_path.exists() and str(bam_path).endswith('.bam'):
                bams.append(bam_path)
    
    return sorted(bams)


def load_bams_from_multiple_tsv(tsv_paths: List[Path]) -> List[Path]:
    """Load BAM paths from multiple TSV files."""
    all_bams = []
    seen = set()
    for tsv_path in tsv_paths:
        bams = load_bams_from_tsv(tsv_path)
        for bam in bams:
            bam_str = str(bam)
            if bam_str not in seen:
                all_bams.append(bam)
                seen.add(bam_str)
    return sorted(all_bams)


def load_ff_from_multiple_tsv(tsv_paths: List[Path]) -> Dict[str, float]:
    """Load FF values from multiple TSV files."""
    ff_map = {}
    for tsv_path in tsv_paths:
        new_map = load_ff_from_tsv(tsv_path)
        ff_map.update(new_map)
    return ff_map


def load_gender_from_multiple_tsv(tsv_paths: List[Path]) -> Dict[str, str]:
    """Load gender from multiple TSV files."""
    gender_map = {}
    for tsv_path in tsv_paths:
        new_map = load_gender_from_tsv(tsv_path)
        gender_map.update(new_map)
    return gender_map


def load_ff_from_tsv(tsv_path: Path) -> Dict[str, float]:
    """Load FF values from TSV file."""
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


def load_gender_from_tsv(tsv_path: Path) -> Dict[str, str]:
    """Load gender from TSV file (Column 3: Gender -> XX/XY -> F/M)."""
    gender_map = {}
    if not tsv_path.exists():
        return gender_map
    
    with tsv_path.open() as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            sample_id = parts[1]  # Column 2: Sample_ID
            gender = parts[2].strip().upper()  # Column 3: Gender
            # Convert XX -> F, XY -> M
            if gender == "XX":
                gender_map[sample_id] = "F"
            elif gender == "XY":
                gender_map[sample_id] = "M"
            else:
                gender_map[sample_id] = "F"  # Default to F
    
    return gender_map


def extract_sample_id(bam_path: Path) -> str:
    """Extract sample ID from BAM path."""
    return bam_path.parent.name


def generate_sample_id(mom_idx: int, fetus_idx: int, disease_name: str, 
                       ff_target: int, reads: int, del_size_mb: int, gender: str) -> str:
    """Generate sample ID in format: {mom_idx}_{fetus_idx}_{disease}_FF{ff}_{reads}M_{del}Mb_{gender}"""
    reads_m = reads // 1_000_000
    return f"{mom_idx}_{fetus_idx}_{disease_name}_FF{ff_target}_{reads_m}M_{del_size_mb}Mb_{gender}"


def get_bed_deletion_info(bed_file: Path) -> Tuple[str, int]:
    """Extract disease name and deletion size from BED file.
    Returns: (disease_name, del_size_mb)
    """
    with bed_file.open() as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                disease_name = parts[3].strip()
                # Clean disease name
                disease_name = ''.join(c for c in disease_name if c.isalnum() or c == '_').lower()
                if not disease_name:
                    disease_name = bed_file.stem
                
                # Calculate deletion size
                if len(parts) >= 3:
                    try:
                        start = int(parts[1])
                        end = int(parts[2])
                        del_size_mb = (end - start) // 1_000_000
                        return disease_name, del_size_mb
                    except ValueError:
                        pass
    
    # Fallback
    return bed_file.stem, 0


def parse_bed_diseases(bed_file: Path) -> List[str]:
    """Parse disease names from BED file (column 4)."""
    diseases = []
    with bed_file.open() as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                disease = parts[3].strip()
                if disease and disease not in diseases:
                    diseases.append(disease)
    return diseases


class TaskManager:
    """Manage artificial sample generation tasks with slot control."""
    
    def __init__(self, max_workers: int, output_dir: Path, script: Path):
        self.max_workers = max_workers
        self.output_dir = output_dir
        self.script = script
        # Note: md_bed is now per-task, not global
        
        # Task queue and status
        self.tasks: List[Dict] = []
        self.task_status: Dict[int, str] = {}  # idx -> status (Init/Running/Completed/Failed)
        self.running_processes: Dict[int, subprocess.Popen] = {}  # idx -> process
        self.task_metadata: Dict[int, Dict] = {}  # idx -> result metadata
        
        # Statistics
        self.completed = 0
        self.failed = 0
        self.start_times: Dict[int, float] = {}
    
    def add_tasks(self, tasks: List[Dict]):
        """Add tasks to the queue."""
        self.tasks = tasks
        for task in tasks:
            idx = task['idx']
            self.task_status[idx] = "Init"
    
    def get_next_task(self) -> Optional[Dict]:
        """Get next task with Init status."""
        for task in self.tasks:
            if self.task_status[task['idx']] == "Init":
                return task
        return None
    
    def get_running_count(self) -> int:
        """Count currently running tasks."""
        return sum(1 for status in self.task_status.values() if status == "Running")
    
    def can_start_new(self) -> bool:
        """Check if we can start a new task."""
        return self.get_running_count() < self.max_workers
    
    def all_finished(self) -> bool:
        """Check if all tasks are completed or failed."""
        return (self.completed + self.failed) == len(self.tasks)
    
    def start_task(self, task: Dict) -> bool:
        """Start a single task."""
        idx = task['idx']
        mom_bam = task['mom_bam']
        fetus_bam = task['fetus_bam']
        ff_target = task['ff_target']
        reads = task['reads']  # Changed from pairs to reads
        mom_idx = task['mom_idx']
        fetus_idx = task['fetus_idx']
        gender = task['gender']
        sample_id = task['sample_id']
        bed_file = task['bed_file']  # BED file is now per-task
        
        mom_id = extract_sample_id(mom_bam)
        fetus_id = extract_sample_id(fetus_bam)
        
        # Create output directory: each sample gets its own directory
        sample_dir = self.output_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        
        # Log file in sample directory
        log_file = sample_dir / "generation.log"
        
        # Create FF map file
        ff_map_file = sample_dir / "ff_map.tsv"
        with ff_map_file.open('w') as f:
            f.write(f"{mom_id}\t{task['mom_ff']}\n")
            f.write(f"{fetus_id}\t{task['fetus_ff']}\n")
        
        # Build command with new parameters
        cmd = [
            str(self.script.resolve()),  # Use absolute path
            "--mom_bam", str(mom_bam),
            "--fetus_bam", str(fetus_bam),
            "--ff_map", str(ff_map_file),
            "--md_bed", str(bed_file),  # Use per-task BED file
            "--ff_target", str(ff_target),
            "--reads", str(reads),  # Changed from --pairs
            "--mom_idx", str(mom_idx),
            "--fetus_idx", str(fetus_idx),
            "--gender", gender,
            "--sample_id", sample_id,
            "--outdir", str(sample_dir)
        ]
        
        try:
            # Start process - run from script directory to ensure relative paths work
            script_dir = self.script.parent
            
            with open(log_file, 'w') as lf:
                lf.write(f"Command: {' '.join(cmd)}\n")
                lf.write(f"Working directory: {script_dir}\n")
                lf.write(f"Start: {datetime.now()}\n\n")
            
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                cwd=str(script_dir)  # Change to script directory
            )
            
            self.running_processes[idx] = process
            self.task_status[idx] = "Running"
            self.start_times[idx] = time.time()
            
            # Store task info for later processing
            self.task_metadata[idx] = {
                'task': task,
                'log_file': log_file,
                'ff_map_file': ff_map_file,
                'sample_dir': sample_dir,
                'sample_id': sample_id,
                'mom_id': mom_id,
                'fetus_id': fetus_id
            }
            
            logger.info(f"[{idx:04d}] Started: {sample_id}")
            return True
            
        except Exception as e:
            logger.error(f"[{idx:04d}] Failed to start: {e}")
            self.task_status[idx] = "Failed"
            self.failed += 1
            if ff_map_file.exists():
                ff_map_file.unlink()
            return False
    
    def check_running_tasks(self):
        """Check status of running tasks and process completed ones."""
        finished_indices = []
        
        for idx, process in list(self.running_processes.items()):
            returncode = process.poll()
            
            if returncode is not None:  # Process finished
                finished_indices.append(idx)
                
                # Get output
                stdout, stderr = process.communicate()
                elapsed = time.time() - self.start_times[idx]
                
                metadata = self.task_metadata[idx]
                task = metadata['task']
                log_file = metadata['log_file']
                ff_map_file = metadata['ff_map_file']
                sample_dir = metadata['sample_dir']
                sample_id = metadata['sample_id']
                
                # Write log
                with log_file.open('a') as lf:
                    lf.write("\n=== STDOUT ===\n")
                    lf.write(stdout)
                    lf.write("\n=== STDERR ===\n")
                    lf.write(stderr)
                    lf.write(f"\nElapsed: {elapsed:.1f}s\n")
                
                if returncode == 0:
                    # Success - parse results and rename output
                    self.task_status[idx] = "Completed"
                    self.completed += 1
                    
                    # Parse alpha, beta from stdout
                    alpha = beta = None
                    for line in stdout.splitlines():
                        if "alpha=" in line and "beta=" in line:
                            match = re.search(r'alpha=([0-9.]+).*beta=([0-9.]+)', line)
                            if match:
                                alpha = float(match.group(1))
                                beta = float(match.group(2))
                                break
                    
                    # Parse deletion check results
                    upstream_mom = upstream_output = upstream_ratio = "N/A"
                    deletion_mom = deletion_output = deletion_ratio = "N/A"
                    downstream_mom = downstream_output = downstream_ratio = "N/A"
                    
                    for line in stdout.splitlines():
                        if "Upstream" in line and "Mom=" in line:
                            match = re.search(r'Mom=(\d+).*Output=(\d+).*Ratio=([0-9.]+|N/A)', line)
                            if match:
                                upstream_mom = match.group(1)
                                upstream_output = match.group(2)
                                upstream_ratio = match.group(3)
                        
                        if "Deletion" in line and "Mom=" in line:
                            match = re.search(r'Mom=(\d+).*Output=(\d+).*Ratio=([0-9.]+|N/A)', line)
                            if match:
                                deletion_mom = match.group(1)
                                deletion_output = match.group(2)
                                deletion_ratio = match.group(3)
                        
                        if "Downstream" in line and "Mom=" in line:
                            match = re.search(r'Mom=(\d+).*Output=(\d+).*Ratio=([0-9.]+|N/A)', line)
                            if match:
                                downstream_mom = match.group(1)
                                downstream_output = match.group(2)
                                downstream_ratio = match.group(3)
                    
                    # Check for output BAM (should be named {sample_id}.proper_paired.bam)
                    final_bam = sample_dir / f"{sample_id}.proper_paired.bam"
                    final_bai = sample_dir / f"{sample_id}.proper_paired.bam.bai"
                    
                    # If output.bam exists (old format), rename it
                    temp_bam = sample_dir / "output.bam"
                    if temp_bam.exists() and not final_bam.exists():
                        temp_bam.rename(final_bam)
                        temp_bai = sample_dir / "output.bam.bai"
                        if temp_bai.exists():
                            temp_bai.rename(final_bai)
                    
                    # Store result metadata
                    self.task_metadata[idx]['result'] = {
                        "success": True,
                        "idx": idx,
                        "sample_id": sample_id,
                        "bam_name": f"{sample_id}.proper_paired.bam",
                        "ff_target": task['ff_target'],
                        "reads": task['reads'],
                        "disease": task['disease'],
                        "mom_idx": task['mom_idx'],
                        "fetus_idx": task['fetus_idx'],
                        "mom_id": metadata['mom_id'],
                        "fetus_id": metadata['fetus_id'],
                        "mom_ff": task['mom_ff'],
                        "fetus_ff": task['fetus_ff'],
                        "gender": task['gender'],
                        "alpha": alpha,
                        "beta": beta,
                        "upstream_mom": upstream_mom,
                        "upstream_output": upstream_output,
                        "upstream_ratio": upstream_ratio,
                        "deletion_mom": deletion_mom,
                        "deletion_output": deletion_output,
                        "deletion_ratio": deletion_ratio,
                        "downstream_mom": downstream_mom,
                        "downstream_output": downstream_output,
                        "downstream_ratio": downstream_ratio,
                        "bam_path": str(final_bam),
                        "elapsed_sec": elapsed
                    }
                    
                    logger.info(f"[{idx:04d}] ✓ Completed: {sample_id} ({elapsed:.1f}s)")
                    
                else:
                    # Failed
                    self.task_status[idx] = "Failed"
                    self.failed += 1
                    
                    self.task_metadata[idx]['result'] = {
                        "success": False,
                        "idx": idx,
                        "error": f"Return code {returncode}",
                        "elapsed_sec": elapsed
                    }
                    
                    logger.error(f"[{idx:04d}] ✗ Failed: {sample_id} (rc={returncode})")
                
                # Cleanup temporary ff_map file
                if ff_map_file.exists():
                    ff_map_file.unlink()
        
        # Remove finished processes from running dict
        for idx in finished_indices:
            del self.running_processes[idx]
    
    def print_status(self):
        """Print current status."""
        running = self.get_running_count()
        init = sum(1 for s in self.task_status.values() if s == "Init")
        total = len(self.tasks)
        
        logger.info("=" * 60)
        logger.info(f"Status: Running={running}/{self.max_workers}, "
                   f"Completed={self.completed}, Failed={self.failed}, "
                   f"Pending={init}, Total={total}")
        logger.info("=" * 60)
    
    def get_results(self) -> List[Dict]:
        """Get all successful results."""
        results = []
        for idx, metadata in self.task_metadata.items():
            if 'result' in metadata and metadata['result']['success']:
                results.append(metadata['result'])
        return results


def save_summary(metadata_list: List[Dict], summary_dir: Path):
    """Save comprehensive summary of all generated samples."""
    summary_dir.mkdir(parents=True, exist_ok=True)
    
    # All samples TSV
    tsv_path = summary_dir / "all_samples.tsv"
    with tsv_path.open("w") as f:
        f.write("idx\tsample_id\tbam_name\tff_target\treads\tmom_idx\tfetus_idx\tgender\tdisease\tmom_id\tfetus_id\tmom_ff\tfetus_ff\talpha\tbeta\t"
                "upstream_mom\tupstream_out\tupstream_ratio\t"
                "deletion_mom\tdeletion_out\tdeletion_ratio\t"
                "downstream_mom\tdownstream_out\tdownstream_ratio\t"
                "bam_path\telapsed_sec\n")
        for m in sorted(metadata_list, key=lambda x: x["idx"]):
            f.write(f"{m['idx']}\t{m['sample_id']}\t{m['bam_name']}\t{m['ff_target']}\t{m['reads']}\t"
                   f"{m['mom_idx']}\t{m['fetus_idx']}\t{m['gender']}\t"
                   f"{m['disease']}\t{m['mom_id']}\t{m['fetus_id']}\t"
                   f"{m.get('mom_ff', 'N/A')}\t{m.get('fetus_ff', 'N/A')}\t"
                   f"{m.get('alpha', 'N/A')}\t{m.get('beta', 'N/A')}\t"
                   f"{m.get('upstream_mom', 'N/A')}\t{m.get('upstream_output', 'N/A')}\t{m.get('upstream_ratio', 'N/A')}\t"
                   f"{m.get('deletion_mom', 'N/A')}\t{m.get('deletion_output', 'N/A')}\t{m.get('deletion_ratio', 'N/A')}\t"
                   f"{m.get('downstream_mom', 'N/A')}\t{m.get('downstream_output', 'N/A')}\t{m.get('downstream_ratio', 'N/A')}\t"
                   f"{m['bam_path']}\t{m['elapsed_sec']:.1f}\n")
    
    # Sample mixing log
    mix_log = summary_dir / "sample_mix.log"
    with mix_log.open("w") as f:
        f.write("# Artificial Sample Mixing Log\n")
        f.write("# Format: Index | BAM File | Mom (FF%) | Fetus (FF%) | Mix Ratio (α:β) | Deletion Check\n")
        f.write("#" + "="*100 + "\n\n")
        
        for m in sorted(metadata_list, key=lambda x: x["idx"]):
            alpha = m.get('alpha', 'N/A')
            beta = m.get('beta', 'N/A')
            mom_ff = m.get('mom_ff', 'N/A')
            fetus_ff = m.get('fetus_ff', 'N/A')
            
            up_ratio = m.get('upstream_ratio', 'N/A')
            del_ratio = m.get('deletion_ratio', 'N/A')
            down_ratio = m.get('downstream_ratio', 'N/A')
            
            f.write(f"{m['idx']:4d}.\t{m['sample_id']}\n")
            f.write(f"\tMom: {m['mom_id']} (idx={m['mom_idx']}, FF={mom_ff}%)\n")
            f.write(f"\tFetus: {m['fetus_id']} (idx={m['fetus_idx']}, FF={fetus_ff}%, Gender={m['gender']})\n")
            f.write(f"\tMixing: α={alpha} (Mom) + β={beta} (Fetus) → Target FF={m['ff_target']}%\n")
            f.write(f"\tRegion Check:\n")
            f.write(f"\t  Upstream   Ratio: {up_ratio}\n")
            f.write(f"\t  Deletion   Ratio: {del_ratio} ← Should be lower!\n")
            f.write(f"\t  Downstream Ratio: {down_ratio}\n")
            f.write("\n")
    
    logger.info(f"\n✓ Summary saved to:")
    logger.info(f"  - {tsv_path}")
    logger.info(f"  - {mix_log}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate artificial MD samples with slot management (v3)'
    )
    
    # Input
    parser.add_argument('--mom_tsv', type=Path, required=True, 
                       help='Mom samples TSV file')
    parser.add_argument('--female_tsv', type=Path, required=True, 
                       help='Female fetus samples TSV file')
    parser.add_argument('--male_tsv', type=Path, required=True, 
                       help='Male fetus samples TSV file')
    parser.add_argument('--md_bed', type=str, required=True, 
                       help='Microdeletion BED file(s), comma-separated for multiple files')
    parser.add_argument('--script', type=Path, default=None,
                       help='make_artificial.sh path (default: same directory as this script)')
    
    # Output
    parser.add_argument('--output', type=Path, default=Path('test_output'),
                       help='Output directory')
    
    # Parameters
    parser.add_argument('--ff_targets', type=str, default='5,10,15',
                       help='Comma-separated FF targets (e.g., 5,10,15)')
    parser.add_argument('--coverages', type=str, default='10M',
                       help='Comma-separated coverages (e.g., 5M,10M,15M)')
    parser.add_argument('--max_workers', type=int, default=4,
                       help='Maximum concurrent workers (slot limit)')
    parser.add_argument('--poll_interval', type=int, default=5,
                       help='Status check interval in seconds (default: 5)')
    
    # Sample selection
    parser.add_argument('--n_moms', type=int, default=None,
                       help='Limit number of mom samples (default: all)')
    parser.add_argument('--n_fetuses', type=int, default=None,
                       help='Limit number of fetus samples per gender (default: all)')
    parser.add_argument('--limit_samples', type=int, default=None,
                       help='Total sample limit (for quick testing)')
    
    args = parser.parse_args()
    
    # Set default script path if not provided
    if args.script is None:
        # Use make_artificial.sh in the same directory as this script
        script_dir = Path(__file__).parent
        args.script = script_dir / 'make_artificial.sh'
    
    # Make script path absolute
    args.script = args.script.resolve()
    
    # Parse multiple BED files
    bed_files_str = args.md_bed.split(',')
    bed_files = [Path(b.strip()) for b in bed_files_str]
    for bed_file in bed_files:
        if not bed_file.exists():
            logger.error(f"BED file not found: {bed_file}")
            return 1
    
    # Validate Mom TSV file
    if not args.mom_tsv.exists():
        logger.error(f"Mom TSV file not found: {args.mom_tsv}")
        return 1
    
    # Validate
    if not args.script.exists():
        logger.error(f"Script not found: {args.script}")
        return 1
    
    # Parse BED file info for each BED file
    bed_info = {}
    for bed_file in bed_files:
        disease_name, del_size_mb = get_bed_deletion_info(bed_file)
        if not disease_name:
            logger.error(f"Could not parse disease info from {bed_file}")
            return 1
        bed_info[bed_file] = (disease_name, del_size_mb)
        logger.info(f"BED {bed_file.name}: Disease={disease_name}, Deletion size={del_size_mb}Mb")
    
    # Load BAMs
    logger.info("\nLoading samples...")
    mom_bams = load_bams_from_tsv(args.mom_tsv)
    female_bams = load_bams_from_tsv(args.female_tsv)
    male_bams = load_bams_from_tsv(args.male_tsv)
    
    # Load FF values and gender
    mom_ff_map = load_ff_from_tsv(args.mom_tsv)
    female_ff_map = load_ff_from_tsv(args.female_tsv)
    male_ff_map = load_ff_from_tsv(args.male_tsv)
    
    female_gender_map = load_gender_from_tsv(args.female_tsv)
    male_gender_map = load_gender_from_tsv(args.male_tsv)
    
    
    # Apply limits
    if args.n_moms:
        mom_bams = mom_bams[:args.n_moms]
    if args.n_fetuses:
        female_bams = female_bams[:args.n_fetuses]
        male_bams = male_bams[:args.n_fetuses]
    
    logger.info(f"  Moms: {len(mom_bams)}")
    logger.info(f"  Female fetuses: {len(female_bams)}")
    logger.info(f"  Male fetuses: {len(male_bams)}")
    
    if not mom_bams or (not female_bams and not male_bams):
        logger.error("No moms or fetuses found.")
        return 1
    
    # Parse parameters
    ff_targets = [int(x) for x in args.ff_targets.split(',')]
    coverages_str = args.coverages.split(',')
    reads_list = []
    for c in coverages_str:
        c = c.strip().upper()
        if c.endswith('M'):
            reads_list.append(int(c[:-1]) * 1_000_000)
        else:
            reads_list.append(int(c))
    
    logger.info(f"\nFF targets: {ff_targets}")
    logger.info(f"Reads targets: {[f'{r//1_000_000}M' for r in reads_list]}")
    logger.info(f"Max workers: {args.max_workers}")
    logger.info(f"Poll interval: {args.poll_interval}s")
    
    # Generate task list
    tasks = []
    global_idx = 1
    
    # Track indices for mom and fetus
    mom_idx_map = {bam: idx + 1 for idx, bam in enumerate(mom_bams)}
    female_idx_map = {bam: idx + 1 for idx, bam in enumerate(female_bams)}
    male_idx_map = {bam: idx + 1 for idx, bam in enumerate(male_bams)}
    
    # Process each BED file
    for bed_file in bed_files:
        disease_name, del_size_mb = bed_info[bed_file]
        
        for ff_target in ff_targets:
            for reads in reads_list:
                # Female fetus combinations
                for mom_bam in mom_bams:
                    for fetus_bam in female_bams:
                        mom_id = extract_sample_id(mom_bam)
                        fetus_id = extract_sample_id(fetus_bam)
                        mom_ff = mom_ff_map.get(mom_id, 0.0)
                        fetus_ff = female_ff_map.get(fetus_id, 0.0)
                        gender = female_gender_map.get(fetus_id, "F")
                        
                        mom_idx = mom_idx_map[mom_bam]
                        fetus_idx = female_idx_map[fetus_bam]
                        
                        sample_id = generate_sample_id(
                            mom_idx, fetus_idx, disease_name, 
                            ff_target, reads, del_size_mb, gender
                        )
                        
                        tasks.append({
                            'mom_bam': mom_bam,
                            'fetus_bam': fetus_bam,
                            'ff_target': ff_target,
                            'reads': reads,
                            'disease': disease_name,
                            'disease_label': disease_name,
                            'mom_idx': mom_idx,
                            'fetus_idx': fetus_idx,
                            'gender': gender,
                            'sample_id': sample_id,
                            'idx': global_idx,
                            'mom_ff': mom_ff,
                            'fetus_ff': fetus_ff,
                            'bed_file': bed_file
                        })
                        global_idx += 1
                
                # Male fetus combinations
                for mom_bam in mom_bams:
                    for fetus_bam in male_bams:
                        mom_id = extract_sample_id(mom_bam)
                        fetus_id = extract_sample_id(fetus_bam)
                        mom_ff = mom_ff_map.get(mom_id, 0.0)
                        fetus_ff = male_ff_map.get(fetus_id, 0.0)
                        gender = male_gender_map.get(fetus_id, "M")
                        
                        mom_idx = mom_idx_map[mom_bam]
                        fetus_idx = male_idx_map[fetus_bam]
                        
                        sample_id = generate_sample_id(
                            mom_idx, fetus_idx, disease_name,
                            ff_target, reads, del_size_mb, gender
                        )
                        
                        tasks.append({
                            'mom_bam': mom_bam,
                            'fetus_bam': fetus_bam,
                            'ff_target': ff_target,
                            'reads': reads,
                            'disease': disease_name,
                            'disease_label': disease_name,
                            'mom_idx': mom_idx,
                            'fetus_idx': fetus_idx,
                            'gender': gender,
                            'sample_id': sample_id,
                            'idx': global_idx,
                            'mom_ff': mom_ff,
                            'fetus_ff': fetus_ff,
                            'bed_file': bed_file
                        })
                        global_idx += 1
    
    # Apply sample limit
    if args.limit_samples and args.limit_samples < len(tasks):
        logger.info(f"\n⚠ Limiting to {args.limit_samples} samples (of {len(tasks)} total)")
        tasks = tasks[:args.limit_samples]
    
    logger.info(f"\nTotal samples to generate: {len(tasks)}")
    
    # Create task manager
    manager = TaskManager(
        max_workers=args.max_workers,
        output_dir=args.output,
        script=args.script
    )
    manager.add_tasks(tasks)
    
    # Main loop
    logger.info(f"\n{'='*60}")
    logger.info(f"Starting generation with slot management...")
    logger.info(f"{'='*60}\n")
    
    start_time = time.time()
    
    while not manager.all_finished():
        # Start new tasks if slots available
        while manager.can_start_new():
            next_task = manager.get_next_task()
            if next_task is None:
                break
            manager.start_task(next_task)
        
        # Check running tasks
        manager.check_running_tasks()
        
        # Print status
        manager.print_status()
        
        # Wait before next check
        time.sleep(args.poll_interval)
    
    # Final check
    manager.check_running_tasks()
    
    total_time = time.time() - start_time
    
    # Print final results
    logger.info(f"\n{'='*60}")
    logger.info("Generation Complete!")
    logger.info(f"  Total time: {total_time:.1f}s")
    logger.info(f"  Succeeded: {manager.completed}")
    logger.info(f"  Failed: {manager.failed}")
    logger.info(f"  Total: {len(tasks)}")
    logger.info(f"{'='*60}\n")
    
    # Save summary
    results = manager.get_results()
    if results:
        summary_dir = args.output.parent / "summary"
        save_summary(results, summary_dir)
    
    return 0 if manager.failed == 0 else 1


if __name__ == '__main__':
    sys.exit(main())

