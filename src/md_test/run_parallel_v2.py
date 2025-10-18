#!/usr/bin/env python3
"""
Parallel Artificial Sample Generator v2
Generates artificial microdeletion samples in parallel using make_artificial.sh
"""

import argparse
import subprocess
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import time
import re


def load_bams_from_tsv(tsv_path: Path):
    """Load BAM paths directly from TSV file.
    
    Expected format:
    Work_Dir | Sample_ID | Gender | FF_Method | FF_Value | ... | BAM_Path (last column)
    """
    bams = []
    if not tsv_path.exists():
        return bams
    
    with tsv_path.open() as f:
        for line in f:
            if line.startswith('#'):  # Skip comments
                continue
            parts = line.strip().split("\t")
            if len(parts) < 11:  # Need at least 11 columns
                continue
            bam_path = Path(parts[-1])  # Last column is BAM path
            if bam_path.exists() and str(bam_path).endswith('.bam'):
                bams.append(bam_path)
    
    return sorted(bams)


def load_ff_from_tsv(tsv_path: Path):
    """Load FF values from TSV file.
    
    Expected format:
    Work_Dir | Sample_ID | Gender | FF_Method | FF_Value | ... | BAM_Path
    Column 2 = Sample_ID (index 1)
    Column 5 = FF_Value (index 4)
    """
    ff_map = {}
    if not tsv_path.exists():
        return ff_map
    
    with tsv_path.open() as f:
        for line in f:
            if line.startswith('#'):  # Skip comments
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


def extract_sample_id(bam_path: Path) -> str:
    """Extract sample ID from BAM path."""
    # Example: /path/to/GNCI25080163/GNCI25080163.proper_paired.bam -> GNCI25080163
    return bam_path.parent.name


def out_bam_name(ff_target: int, pairs: int, disease_label: str, idx: int):
    """
    Generate simplified output BAM filename.
    Format: {idx}_{disease}_FF{ff}_{coverage}.bam
    """
    pairs_m = pairs // 1_000_000
    return f"{idx:04d}_{disease_label}_FF{ff_target:02d}_{pairs_m}M.bam"


def run_one(
    script: Path,
    mom_bam: Path,
    fetus_bam: Path,
    ff_target: int,
    pairs: int,
    md_bed: Path,
    out_dir: Path,
    disease_label: str,
    idx: int,
    mom_ff: float,
    fetus_ff: float
):
    """Run one artificial sample generation."""
    start_time = time.time()
    
    # Create output directory structure
    ff_dir = out_dir / f"FF{ff_target:02d}_{pairs//1_000_000}M"
    ff_dir.mkdir(parents=True, exist_ok=True)
    
    # Log file
    log_dir = out_dir.parent / "test_logs" / disease_label / f"FF{ff_target:02d}_{pairs//1_000_000}M"
    log_dir.mkdir(parents=True, exist_ok=True)
    
    final_name = out_bam_name(ff_target, pairs, disease_label, idx)
    log_file = log_dir / f"{final_name.replace('.bam', '.log')}"
    
    mom_id = extract_sample_id(mom_bam)
    fetus_id = extract_sample_id(fetus_bam)
    
    # Create temporary ff_map file for this run
    ff_map_file = ff_dir / f".ff_map_{idx}.tsv"
    with ff_map_file.open('w') as f:
        f.write(f"{mom_id}\t{mom_ff}\n")
        f.write(f"{fetus_id}\t{fetus_ff}\n")
    
    # Run make_artificial.sh with named arguments
    cmd = [
        str(script),
        "--mom_bam", str(mom_bam),
        "--fetus_bam", str(fetus_bam),
        "--ff_map", str(ff_map_file),
        "--md_bed", str(md_bed),
        "--ff_target", str(ff_target),
        "--pairs", str(pairs),
        "--outdir", str(ff_dir)
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Save log
        with log_file.open('w') as f:
            f.write(f"Command: {' '.join(cmd)}\n")
            f.write(f"Start: {datetime.now()}\n")
            f.write("\n=== STDOUT ===\n")
            f.write(result.stdout)
            f.write("\n=== STDERR ===\n")
            f.write(result.stderr)
        
        # Parse alpha, beta from log
        alpha = beta = None
        for line in result.stdout.splitlines():
            if "alpha=" in line and "beta=" in line:
                match = re.search(r'alpha=([0-9.]+).*beta=([0-9.]+)', line)
                if match:
                    alpha = float(match.group(1))
                    beta = float(match.group(2))
                    break
        
        # Parse deletion check results (upstream, deletion, downstream)
        upstream_mom = upstream_output = upstream_ratio = "N/A"
        deletion_mom = deletion_output = deletion_ratio = "N/A"
        downstream_mom = downstream_output = downstream_ratio = "N/A"
        
        for line in result.stdout.splitlines():
            # Upstream
            if "Upstream" in line and "Mom=" in line:
                match = re.search(r'Mom=(\d+).*Output=(\d+).*Ratio=([0-9.]+|N/A)', line)
                if match:
                    upstream_mom = match.group(1)
                    upstream_output = match.group(2)
                    upstream_ratio = match.group(3)
            
            # Deletion
            if "Deletion" in line and "Mom=" in line:
                match = re.search(r'Mom=(\d+).*Output=(\d+).*Ratio=([0-9.]+|N/A)', line)
                if match:
                    deletion_mom = match.group(1)
                    deletion_output = match.group(2)
                    deletion_ratio = match.group(3)
            
            # Downstream
            if "Downstream" in line and "Mom=" in line:
                match = re.search(r'Mom=(\d+).*Output=(\d+).*Ratio=([0-9.]+|N/A)', line)
                if match:
                    downstream_mom = match.group(1)
                    downstream_output = match.group(2)
                    downstream_ratio = match.group(3)
        
        # Rename output.bam to final name
        temp_bam = ff_dir / "output.bam"
        temp_bai = ff_dir / "output.bam.bai"
        final_bam = ff_dir / final_name
        final_bai = ff_dir / f"{final_name}.bai"
        
        if temp_bam.exists():
            temp_bam.rename(final_bam)
            if temp_bai.exists():
                temp_bai.rename(final_bai)
        
        # Cleanup temporary ff_map file
        if ff_map_file.exists():
            ff_map_file.unlink()
        
        elapsed = time.time() - start_time
        
        return {
            "success": True,
            "idx": idx,
            "bam_name": final_name,
            "ff_target": ff_target,
            "pairs": pairs,
            "disease": disease_label,
            "mom_id": mom_id,
            "fetus_id": fetus_id,
            "mom_ff": mom_ff,
            "fetus_ff": fetus_ff,
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
        
    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        with log_file.open('w') as f:
            f.write(f"Command: {' '.join(cmd)}\n")
            f.write(f"ERROR: {e}\n")
            f.write("\n=== STDOUT ===\n")
            f.write(e.stdout)
            f.write("\n=== STDERR ===\n")
            f.write(e.stderr)
        
        # Cleanup temporary ff_map file
        if ff_map_file.exists():
            ff_map_file.unlink()
        
        return {
            "success": False,
            "idx": idx,
            "error": str(e),
            "elapsed_sec": elapsed
        }


def save_summary(metadata_list, summary_dir: Path):
    """Save comprehensive summary of all generated samples."""
    summary_dir.mkdir(parents=True, exist_ok=True)
    
    # All samples TSV with sample mapping
    tsv_path = summary_dir / "all_samples.tsv"
    with tsv_path.open("w") as f:
        f.write("idx\tbam_name\tff_target\tpairs\tdisease\tmom_id\tfetus_id\tmom_ff\tfetus_ff\talpha\tbeta\t"
                "upstream_mom\tupstream_out\tupstream_ratio\t"
                "deletion_mom\tdeletion_out\tdeletion_ratio\t"
                "downstream_mom\tdownstream_out\tdownstream_ratio\t"
                "bam_path\telapsed_sec\n")
        for m in sorted(metadata_list, key=lambda x: x["idx"]):
            f.write(f"{m['idx']}\t{m['bam_name']}\t{m['ff_target']}\t{m['pairs']}\t"
                   f"{m['disease']}\t{m['mom_id']}\t{m['fetus_id']}\t"
                   f"{m.get('mom_ff', 'N/A')}\t{m.get('fetus_ff', 'N/A')}\t"
                   f"{m.get('alpha', 'N/A')}\t{m.get('beta', 'N/A')}\t"
                   f"{m.get('upstream_mom', 'N/A')}\t{m.get('upstream_output', 'N/A')}\t{m.get('upstream_ratio', 'N/A')}\t"
                   f"{m.get('deletion_mom', 'N/A')}\t{m.get('deletion_output', 'N/A')}\t{m.get('deletion_ratio', 'N/A')}\t"
                   f"{m.get('downstream_mom', 'N/A')}\t{m.get('downstream_output', 'N/A')}\t{m.get('downstream_ratio', 'N/A')}\t"
                   f"{m['bam_path']}\t{m['elapsed_sec']:.1f}\n")
    
    # Sample mixing log (human-readable)
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
            
            f.write(f"{m['idx']:4d}.\t{m['bam_name']}\n")
            f.write(f"\tMom: {m['mom_id']} (FF={mom_ff}%)\n")
            f.write(f"\tFetus: {m['fetus_id']} (FF={fetus_ff}%)\n")
            f.write(f"\tMixing: α={alpha} (Mom) + β={beta} (Fetus) → Target FF={m['ff_target']}%\n")
            f.write(f"\tRegion Check:\n")
            f.write(f"\t  Upstream   Ratio: {up_ratio}\n")
            f.write(f"\t  Deletion   Ratio: {del_ratio} ← Should be lower!\n")
            f.write(f"\t  Downstream Ratio: {down_ratio}\n")
            f.write("\n")
    
    print(f"\n✓ Summary saved to:")
    print(f"  - {tsv_path}")
    print(f"  - {mix_log}")


def parse_bed_diseases(bed_file: Path):
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


def main():
    parser = argparse.ArgumentParser(description='Generate artificial MD samples in parallel')
    
    # Input
    parser.add_argument('--mom_tsv', type=Path, required=True, help='Mom samples TSV file')
    parser.add_argument('--female_tsv', type=Path, required=True, help='Female fetus samples TSV file')
    parser.add_argument('--male_tsv', type=Path, required=True, help='Male fetus samples TSV file')
    parser.add_argument('--md_bed', type=Path, required=True, help='Microdeletion BED file')
    parser.add_argument('--script', type=Path, default=Path('make_artificial.sh'), 
                       help='make_artificial.sh path')
    
    # Output
    parser.add_argument('--output', type=Path, default=Path('test_output'),
                       help='Output directory')
    
    # Parameters
    parser.add_argument('--ff_targets', type=str, default='5,10,15',
                       help='Comma-separated FF targets (e.g., 5,10,15)')
    parser.add_argument('--coverages', type=str, default='10M',
                       help='Comma-separated coverages (e.g., 5M,10M,15M)')
    parser.add_argument('--workers', type=int, default=4,
                       help='Number of parallel workers')
    
    # Sample selection
    parser.add_argument('--n_moms', type=int, default=None,
                       help='Limit number of mom samples (default: all)')
    parser.add_argument('--n_fetuses', type=int, default=None,
                       help='Limit number of fetus samples per gender (default: all)')
    parser.add_argument('--limit_samples', type=int, default=None,
                       help='Total sample limit (for quick testing)')
    
    args = parser.parse_args()
    
    # Validate
    if not args.script.exists():
        print(f"ERROR: Script not found: {args.script}")
        return 1
    
    if not args.md_bed.exists():
        print(f"ERROR: BED file not found: {args.md_bed}")
        return 1
    
    # Parse diseases
    diseases = parse_bed_diseases(args.md_bed)
    if not diseases:
        print(f"ERROR: No diseases found in {args.md_bed}")
        return 1
    
    print(f"Found {len(diseases)} disease(s): {', '.join(diseases)}")
    
    # Load BAMs
    print("\nLoading samples...")
    mom_bams = load_bams_from_tsv(args.mom_tsv)
    female_bams = load_bams_from_tsv(args.female_tsv)
    male_bams = load_bams_from_tsv(args.male_tsv)
    
    # Load FF values
    mom_ff_map = load_ff_from_tsv(args.mom_tsv)
    female_ff_map = load_ff_from_tsv(args.female_tsv)
    male_ff_map = load_ff_from_tsv(args.male_tsv)
    
    # Apply limits
    if args.n_moms:
        mom_bams = mom_bams[:args.n_moms]
    if args.n_fetuses:
        female_bams = female_bams[:args.n_fetuses]
        male_bams = male_bams[:args.n_fetuses]
    
    print(f"  Moms: {len(mom_bams)}")
    print(f"  Female fetuses: {len(female_bams)}")
    print(f"  Male fetuses: {len(male_bams)}")
    
    if not mom_bams or (not female_bams and not male_bams):
        print("ERROR: no moms or fetuses found.")
        return 1
    
    # Parse parameters
    ff_targets = [int(x) for x in args.ff_targets.split(',')]
    coverages_str = args.coverages.split(',')
    coverages = []
    for c in coverages_str:
        c = c.strip().upper()
        if c.endswith('M'):
            coverages.append(int(c[:-1]) * 1_000_000)
        else:
            coverages.append(int(c))
    
    print(f"\nFF targets: {ff_targets}")
    print(f"Coverages: {[f'{c//1_000_000}M' for c in coverages]}")
    print(f"Workers: {args.workers}")
    
    # Generate task list
    tasks = []
    global_idx = 1
    
    for disease in diseases:
        disease_label = disease.replace(" ", "_").replace("/", "-")
        
        for ff_target in ff_targets:
            for pairs in coverages:
                # Female fetus combinations
                for mom_bam in mom_bams:
                    for fetus_bam in female_bams:
                        mom_id = extract_sample_id(mom_bam)
                        fetus_id = extract_sample_id(fetus_bam)
                        mom_ff = mom_ff_map.get(mom_id, 0.0)
                        fetus_ff = female_ff_map.get(fetus_id, 0.0)
                        
                        tasks.append({
                            'mom_bam': mom_bam,
                            'fetus_bam': fetus_bam,
                            'ff_target': ff_target,
                            'pairs': pairs,
                            'disease': disease,
                            'disease_label': disease_label,
                            'idx': global_idx,
                            'mom_ff': mom_ff,
                            'fetus_ff': fetus_ff
                        })
                        global_idx += 1
                
                # Male fetus combinations
                for mom_bam in mom_bams:
                    for fetus_bam in male_bams:
                        mom_id = extract_sample_id(mom_bam)
                        fetus_id = extract_sample_id(fetus_bam)
                        mom_ff = mom_ff_map.get(mom_id, 0.0)
                        fetus_ff = male_ff_map.get(fetus_id, 0.0)
                        
                        tasks.append({
                            'mom_bam': mom_bam,
                            'fetus_bam': fetus_bam,
                            'ff_target': ff_target,
                            'pairs': pairs,
                            'disease': disease,
                            'disease_label': disease_label,
                            'idx': global_idx,
                            'mom_ff': mom_ff,
                            'fetus_ff': fetus_ff
                        })
                        global_idx += 1
    
    # Apply sample limit
    if args.limit_samples and args.limit_samples < len(tasks):
        print(f"\n⚠ Limiting to {args.limit_samples} samples (of {len(tasks)} total)")
        tasks = tasks[:args.limit_samples]
    
    print(f"\nTotal samples to generate: {len(tasks)}")
    
    # Create output directory
    out_dir = args.output
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Run in parallel
    print(f"\n{'='*60}")
    print(f"Starting parallel generation with {args.workers} workers...")
    print(f"{'='*60}\n")
    
    results = []
    succeeded = 0
    failed = 0
    
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {}
        for task in tasks:
            future = executor.submit(
                run_one,
                args.script,
                task['mom_bam'],
                task['fetus_bam'],
                task['ff_target'],
                task['pairs'],
                args.md_bed,
                out_dir / task['disease_label'],
                task['disease_label'],
                task['idx'],
                task['mom_ff'],
                task['fetus_ff']
            )
            futures[future] = task
        
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            
            if result['success']:
                succeeded += 1
                print(f"[{succeeded}/{len(tasks)}] ✓ {result['bam_name']} ({result['elapsed_sec']:.1f}s)")
            else:
                failed += 1
                print(f"[FAIL] {result.get('error', 'Unknown error')}")
    
    # Save summary
    print(f"\n{'='*60}")
    print(f"Generation complete!")
    print(f"  Succeeded: {succeeded}")
    print(f"  Failed: {failed}")
    print(f"{'='*60}\n")
    
    successful_results = [r for r in results if r['success']]
    if successful_results:
        summary_dir = out_dir.parent / "summary"
        save_summary(successful_results, summary_dir)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

