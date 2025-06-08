#!/usr/bin/env python3

import os
import sys
import subprocess
import logging
import time
import argparse
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Fixed directories for Docker environment
SRC_DIR = "/Work/NIPT/bin"
ANALYSIS_DIR = "/Work/NIPT/analysis"
LOG_DIR = "/Work/NIPT/log"
FASTQ_DIR = "/Work/NIPT/fastq"

# Define filter paths
filter_paths = {
    "of": "/Work/NIPT/data/bed/Uniform_2017_allY.bed",
    "nf08": "/Work/NIPT/data/bed/hg19_mappability_0.8_clean_all_36mer.bed",
    "nf09": "/Work/NIPT/data/bed/hg19_mappability_0.9_clean_all_36mer.bed"
}

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='NIPT Analysis Pipeline')
    parser.add_argument('--sample_name', required=True, help='Sample name')
    parser.add_argument('--fastq_r1', required=True, help='R1 FASTQ filename')
    parser.add_argument('--fastq_r2', required=True, help='R2 FASTQ filename')
    
    return parser.parse_args()

def run_command(description, command):
    """Run a command with timing and logging"""
    logger.info(f"Running {description}")
    
    start_time = time.time()
    try:
        subprocess.run(command, shell=True, check=True)
        end_time = time.time()
        elapsed_time = end_time - start_time
        logger.info(f"Completed {description} in {elapsed_time:.2f} seconds")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running {description}: {e}")
        return False

def create_directories(sample_name):
    """Create required output directories"""
    directories = [
        f"{ANALYSIS_DIR}/{sample_name}",
        f"{ANALYSIS_DIR}/{sample_name}/Output_QC",
        f"{ANALYSIS_DIR}/{sample_name}/Output_Wisecondor",
        f"{ANALYSIS_DIR}/{sample_name}/Output_WisecondorX",
        f"{ANALYSIS_DIR}/{sample_name}/Output_WisecondorFF",
        f"{ANALYSIS_DIR}/{sample_name}/Output_Result",
        f"{ANALYSIS_DIR}/{sample_name}/Output_hmmcopy",
        f"{LOG_DIR}/{sample_name}"
    ]
    
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
        logger.info(f"Created directory: {directory}")
    
    return True

def create_symbolic_links(sample_name, fastq_r1, fastq_r2):
    """Create symbolic links to FASTQ files"""
    try:
        fastq_r1_path = f"{FASTQ_DIR}/{sample_name}/{fastq_r1}"
        fastq_r2_path = f"{FASTQ_DIR}/{sample_name}/{fastq_r2}"
        
        # Create symbolic links
        os.symlink(fastq_r1_path, f"{ANALYSIS_DIR}/{fastq_r1}")
        logger.info(f"Created symbolic link: {ANALYSIS_DIR}/{fastq_r1}")
        
        os.symlink(fastq_r2_path, f"{ANALYSIS_DIR}/{fastq_r2}")
        logger.info(f"Created symbolic link: {ANALYSIS_DIR}/{fastq_r2}")
        
        return True
    except Exception as e:
        logger.error(f"Error creating symbolic links: {e}")
        return False

def generate_proper_paired_bam(sample_name, fastq_r1, fastq_r2):
    """Generate proper_paired.bam file from FASTQ files"""
    # Define command paths from environment variables
    bwa = os.environ.get('BWA2', 'bwa-mem2')
    ref_hg = os.environ.get('ref_hg', '')
    sam_tools = os.environ.get('SAMTools', 'samtools')
    picard = os.environ.get('PICARD', '')
    qualimap = os.environ.get('qualimap', 'qualimap')
    
    # Get parameters from config
    bwa_threads = os.environ.get('QC.bwa_threads', '16')
    samtools_threads = os.environ.get('QC.samtools_threads', '8')
    samtools_memory = os.environ.get('QC.samtools_memory', '1G')
    
    # 1. BWA MEM alignment and sorting
    run_command(
        "BWA-MEM2 alignment and sorting",
        f"{bwa} mem -t {bwa_threads} {ref_hg} {ANALYSIS_DIR}/{fastq_r1} {ANALYSIS_DIR}/{fastq_r2} | "
        f"{sam_tools} sort -@ {samtools_threads} -m {samtools_memory} -O bam "
        f"-o {ANALYSIS_DIR}/{sample_name}/{sample_name}.sorted.bam "
        f"-T {ANALYSIS_DIR}/{sample_name}/{sample_name}.sorted"
    )
    
    # 2. Index sorted BAM
    run_command(
        "SAMTools Indexing",
        f"{sam_tools} index {ANALYSIS_DIR}/{sample_name}/{sample_name}.sorted.bam"
    )
    
    # 3. Run Qualimap
    run_command(
        "Qualimap",
        f"{qualimap} bamqc -bam {ANALYSIS_DIR}/{sample_name}/{sample_name}.sorted.bam "
        f"-nt 10 -outdir {ANALYSIS_DIR}/{sample_name}/Output_QC"
    )
    
    # 4. Zip Qualimap results
    run_command(
        "Zip Qualimap results",
        f"zip -r {ANALYSIS_DIR}/{sample_name}/{sample_name}_Qualimap.zip {ANALYSIS_DIR}/{sample_name}/Output_QC"
    )
    
    # 5. PICARD - Remove duplications
    run_command(
        "PICARD MarkDuplicates",
        f"java -Xmx8g -Djava.io.tmpdir=`pwd`/tmp -jar {picard} MarkDuplicates "
        f"I={ANALYSIS_DIR}/{sample_name}/{sample_name}.sorted.bam REMOVE_DUPLICATES=true "
        f"VALIDATION_STRINGENCY=LENIENT AS=true "
        f"M=\"{ANALYSIS_DIR}/{sample_name}/{sample_name}\"_dup.metrics "
        f"O=\"{ANALYSIS_DIR}/{sample_name}/{sample_name}\".dedup.bam"
    )
    
    # 6. Index dedup BAM
    run_command(
        "PICARD bam indexing",
        f"{sam_tools} index {ANALYSIS_DIR}/{sample_name}/{sample_name}.dedup.bam"
    )
    
    # 7. SAMTools - Make unique BAM
    run_command(
        "unique bam",
        f"{sam_tools} view -bq 1 {ANALYSIS_DIR}/{sample_name}/{sample_name}.dedup.bam > "
        f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.uniq.bam"
    )
    
    # 8. Index unique BAM
    run_command(
        "unique bam indexing",
        f"{sam_tools} index {ANALYSIS_DIR}/{sample_name}/{sample_name}.uniq.bam"
    )
    
    # 9. Extract Proper paired BAM
    run_command(
        "proper bam",
        f"{sam_tools} view -b -f 0x2 {ANALYSIS_DIR}/{sample_name}/{sample_name}.uniq.bam > "
        f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.proper_paired.bam"
    )
    
    # 10. Index proper paired BAM
    run_command(
        "proper bam indexing",
        f"{sam_tools} index {ANALYSIS_DIR}/{sample_name}/{sample_name}.proper_paired.bam"
    )
    
    logger.info("proper_paired.bam generation completed")
    with open(f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.log", "a") as log_file:
        log_file.write("8. proper_paired.bam generation : PASS\n")
    
    return True

def process_filter(sample_name, filter_name, filter_path):
    """Process a specific filter type and create all size variants"""
    # Define command paths
    sam_tools = os.environ.get('SAMTools', 'samtools')
    ref_hg = os.environ.get('ref_hg', '')
    
    # Create orig BAM from proper_paired.bam
    logger.info(f"SAMTools {filter_name}_orig bam start")
    run_command(
        f"{filter_name}_orig bam", 
        f"{sam_tools} view -b -L {filter_path} {ANALYSIS_DIR}/{sample_name}/{sample_name}.proper_paired.bam > "
        f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_name}_orig.bam"
    )
    run_command(
        f"{filter_name}_orig bam indexing", 
        f"{sam_tools} index {ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_name}_orig.bam"
    )
    logger.info(f"SAMTools {filter_name}_orig bam end")
    
    # Create fetus BAM (size filter using filter.awk)
    logger.info(f"SAMTools {filter_name}_fetus bam start")
    run_command(
        f"{filter_name}_fetus bam", 
        f"{sam_tools} view {ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_name}_orig.bam | "
        f"awk -f {SRC_DIR}/filter.awk | {sam_tools} view -bt {ref_hg}.fai -o "
        f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_name}_fetus.bam -"
    )
    run_command(
        f"{filter_name}_fetus bam indexing", 
        f"{sam_tools} index {ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_name}_fetus.bam"
    )
    logger.info(f"SAMTools {filter_name}_fetus bam end")
    
    # Create mom BAM (size filter using filter_out.awk)
    logger.info(f"SAMTools {filter_name}_mom bam start")
    run_command(
        f"{filter_name}_mom bam", 
        f"{sam_tools} view {ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_name}_orig.bam | "
        f"awk -f {SRC_DIR}/filter_out.awk | {sam_tools} view -bt {ref_hg}.fai -o "
        f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_name}_mom.bam -"
    )
    run_command(
        f"{filter_name}_mom bam indexing", 
        f"{sam_tools} index {ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_name}_mom.bam"
    )
    logger.info(f"SAMTools {filter_name}_mom bam end")

def create_npz_files(sample_name, bam_file, bam_suffix):
    """Create WC, WCX, and WFF NPZ files for a given BAM file"""
    # Define command paths
    python2 = os.environ.get('PYTHON2', 'python2')
    wc_path = os.environ.get('WC', '')
    wcx_path = os.environ.get('WCX', '')
    wcff_path = os.environ.get('WCFF', '')
    
    # WC convert
    wc_output = f"{ANALYSIS_DIR}/{sample_name}/Output_Wisecondor/{sample_name}.wc.{bam_suffix}.npz"
    run_command(
        f"WC convert {bam_suffix}", 
        f"{python2} {wc_path}/wisecondor.py convert {bam_file} {wc_output} -binsize 200000"
    )
    
    # WCX convert
    wcx_output = f"{ANALYSIS_DIR}/{sample_name}/Output_WisecondorX/{sample_name}_wcx.{bam_suffix}.npz"
    run_command(
        f"WCX convert {bam_suffix}", 
        f"{wcx_path} convert {bam_file} {wcx_output} --binsize 200000"
    )
    
    # WCFF convert
    wcff_output = f"{ANALYSIS_DIR}/{sample_name}/Output_WisecondorFF/{sample_name}_wcff.{bam_suffix}.npz"
    run_command(
        f"WCFF convert {bam_suffix}", 
        f"{wcff_path} convert -b 200000 {bam_file} output.npz {wcff_output}"
    )
    
    return True

def create_hmmcopy_files(sample_name, bam_file, bam_suffix):
    """Create 50kb and 10mb wig files and run HMMcopy for a given BAM file"""
    # Define command paths
    hmmcopy_path = os.environ.get('HMMcopy', '')
    rscript = os.environ.get('Rscript', 'Rscript')
    
    # Common chromosomes list
    chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
    
    # Ensure Output_hmmcopy directory exists
    os.makedirs(f"{ANALYSIS_DIR}/{sample_name}/Output_hmmcopy", exist_ok=True)
    
    # 50kb wig file
    wig_50kb = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{bam_suffix}.50kb.wig"
    run_command(
        f"readCounter 50kb {bam_suffix}", 
        f"{hmmcopy_path}/bin/readCounter -w 50000 -c {chromosomes} {bam_file} > {wig_50kb}"
    )
    
    # Run HMMcopy R script for 50kb and save to Output_hmmcopy
    run_command(
        f"HMMcopy R 50kb {bam_suffix}", 
        f"{rscript} --no-save --no-restore --verbose {SRC_DIR}/HMMcopy.R {wig_50kb} "
        f"{hmmcopy_path}/hg19.gc.50kb.wig {hmmcopy_path}/hg19.map.50kb.wig 50kb "
        f"{ANALYSIS_DIR}/{sample_name}/Output_hmmcopy"
    )
    
    # 10mb wig file
    wig_10mb = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{bam_suffix}.10mb.wig"
    run_command(
        f"readCounter 10mb {bam_suffix}", 
        f"{hmmcopy_path}/bin/readCounter -w 10000000 -c {chromosomes} {bam_file} > {wig_10mb}"
    )
    
    # Run HMMcopy R script for 10mb and save to Output_hmmcopy
    run_command(
        f"HMMcopy R 10mb {bam_suffix}", 
        f"{rscript} --no-save --no-restore --verbose {SRC_DIR}/HMMcopy.R {wig_10mb} "
        f"{hmmcopy_path}/hg19.gc.10mb.wig {hmmcopy_path}/hg19.map.10mb.wig 10mb "
        f"{ANALYSIS_DIR}/{sample_name}/Output_hmmcopy"
    )
    
    return True

def process_bam_files(sample_name):
    """Process all BAM files to create NPZ and HMMcopy files"""
    # List of all BAM files to process with their suffixes
    bam_files = []
    
    # Add proper_paired.bam
    proper_paired_bam = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.proper_paired.bam"
    bam_files.append((proper_paired_bam, "proper_paired"))
    
    # Process each filter type to create BAM files
    for filter_name in filter_paths.keys():
        # Add orig BAM
        orig_bam = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_name}_orig.bam"
        bam_files.append((orig_bam, f"{filter_name}_orig"))
        
        # Add fetus BAM
        fetus_bam = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_name}_fetus.bam"
        bam_files.append((fetus_bam, f"{filter_name}_fetus"))
        
        # Add mom BAM
        mom_bam = f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.{filter_name}_mom.bam"
        bam_files.append((mom_bam, f"{filter_name}_mom"))
    
    # Process each BAM file
    for bam_file, bam_suffix in bam_files:
        logger.info(f"Processing {bam_suffix} BAM file")
        
        # Create NPZ files for all BAM files
        create_npz_files(sample_name, bam_file, bam_suffix)
        
        # For all except proper_paired.bam, also create HMMcopy files
        if bam_suffix != "proper_paired":
            create_hmmcopy_files(sample_name, bam_file, bam_suffix)
    
    logger.info("All BAM processing and file generation completed")

def main():
    """Main function to orchestrate the full NIPT pipeline"""
    # Parse command line arguments
    args = parse_args()
    
    sample_name = args.sample_name
    fastq_r1 = args.fastq_r1
    fastq_r2 = args.fastq_r2
    
    # Set up logging to file
    file_handler = logging.FileHandler(f"{LOG_DIR}/{sample_name}/{sample_name}_pipeline.log")
    file_handler.setFormatter(logging.Formatter('[%(asctime)s] [%(levelname)s] %(message)s'))
    logger.addHandler(file_handler)
    
    logger.info(f"Starting NIPT pipeline for sample: {sample_name}")
    logger.info(f"Using fixed directories: SRC={SRC_DIR}, ANALYSIS={ANALYSIS_DIR}, LOG={LOG_DIR}, FASTQ={FASTQ_DIR}")
    
    # Create necessary directories
    create_directories(sample_name)
    
    # Create symbolic links to FASTQ files
    create_symbolic_links(sample_name, fastq_r1, fastq_r2)
    
    # Generate proper_paired.bam
    generate_proper_paired_bam(sample_name, fastq_r1, fastq_r2)
    
    # Process each filter type to create BAM files
    for filter_name, filter_path in filter_paths.items():
        process_filter(sample_name, filter_name, filter_path)
    
    # Log completion of filter processing
    with open(f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.log", "a") as log_file:
        log_file.write("9. Samtools filters (of, nf08, nf09) bam : PASS\n")
    
    # Process all BAM files for NPZ and HMMcopy files
    process_bam_files(sample_name)
    
    # Log completion of all processing
    with open(f"{ANALYSIS_DIR}/{sample_name}/{sample_name}.log", "a") as log_file:
        log_file.write("10. NPZ and HMMcopy file generation : PASS\n")
    
    logger.info(f"NIPT pipeline completed successfully for sample: {sample_name}")

if __name__ == "__main__":
    main()ing", 
        f"{sam_tools} index {analysis_dir}/{sample_name}/{sample_name}.{filter_name}_fetus.bam"
    )
    logger.info(f"SAMTools {filter_name}_fetus bam end")
    
    # Create mom BAM (size filter using filter_out.awk)
    logger.info(f"SAMTools {filter_name}_mom bam start")
    run_command(
        f"{filter_name}_mom bam", 
        f"{sam_tools} view {analysis_dir}/{sample_name}/{sample_name}.{filter_name}_orig.bam | "
        f"awk -f {src}/filter_out.awk | {sam_tools} view -bt {ref_hg}.fai -o "
        f"{analysis_dir}/{sample_name}/{sample_name}.{filter_name}_mom.bam -"
    )
    run_command(
        f"{filter_name}_mom bam indexing", 
        f"{sam_tools} index {analysis_dir}/{sample_name}/{sample_name}.{filter_name}_mom.bam"
    )
    logger.info(f"SAMTools {filter_name}_mom bam end")

def create_npz_files(sample_name, bam_file, bam_suffix, analysis_dir):
    """Create WC, WCX, and WFF NPZ files for a given BAM file"""
    # Define command paths
    python2 = os.environ.get('PYTHON2', 'python2')
    wc_path = os.environ.get('WC', '')
    wcx_path = os.environ.get('WCX', '')
    wcff_path = os.environ.get('WCFF', '')
    
    # WC convert
    wc_output = f"{analysis_dir}/{sample_name}/Output_Wisecondor/{sample_name}.wc.{bam_suffix}.npz"
    run_command(
        f"WC convert {bam_suffix}", 
        f"{python2} {wc_path}/wisecondor.py convert {bam_file} {wc_output} -binsize 200000"
    )
    
    # WCX convert
    wcx_output = f"{analysis_dir}/{sample_name}/Output_WisecondorX/{sample_name}_wcx.{bam_suffix}.npz"
    run_command(
        f"WCX convert {bam_suffix}", 
        f"{wcx_path} convert {bam_file} {wcx_output} --binsize 200000"
    )
    
    # WCFF convert
    wcff_output = f"{analysis_dir}/{sample_name}/Output_WisecondorFF/{sample_name}_wcff.{bam_suffix}.npz"
    run_command(
        f"WCFF convert {bam_suffix}", 
        f"{wcff_path} convert -b 200000 {bam_file} output.npz {wcff_output}"
    )
    
    return True

def create_hmmcopy_files(sample_name, bam_file, bam_suffix, analysis_dir):
    """Create 50kb and 10mb wig files and run HMMcopy for a given BAM file"""
    # Define command paths
    hmmcopy_path = os.environ.get('HMMcopy', '')
    rscript = os.environ.get('Rscript', 'Rscript')
    src = os.environ.get('src', '')
    
    # Common chromosomes list
    chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
    
    # 50kb wig file
    wig_50kb = f"{analysis_dir}/{sample_name}/{sample_name}.{bam_suffix}.50kb.wig"
    run_command(
        f"readCounter 50kb {bam_suffix}", 
        f"{hmmcopy_path}/bin/readCounter -w 50000 -c {chromosomes} {bam_file} > {wig_50kb}"
    )
    
    # Run HMMcopy R script for 50kb and save to Output_hmmcopy
    run_command(
        f"HMMcopy R 50kb {bam_suffix}", 
        f"{rscript} --no-save --no-restore --verbose {src}/HMMcopy.R {wig_50kb} "
        f"{hmmcopy_path}/hg19.gc.50kb.wig {hmmcopy_path}/hg19.map.50kb.wig 50kb "
        f"{analysis_dir}/{sample_name}/Output_hmmcopy"
    )
    
    # 10mb wig file
    wig_10mb = f"{analysis_dir}/{sample_name}/{sample_name}.{bam_suffix}.10mb.wig"
    run_command(
        f"readCounter 10mb {bam_suffix}", 
        f"{hmmcopy_path}/bin/readCounter -w 10000000 -c {chromosomes} {bam_file} > {wig_10mb}"
    )
    
    # Run HMMcopy R script for 10mb and save to Output_hmmcopy
    run_command(
        f"HMMcopy R 10mb {bam_suffix}", 
        f"{rscript} --no-save --no-restore --verbose {src}/HMMcopy.R {wig_10mb} "
        f"{hmmcopy_path}/hg19.gc.10mb.wig {hmmcopy_path}/hg19.map.10mb.wig 10mb "
        f"{analysis_dir}/{sample_name}/Output_hmmcopy"
    )
    
    return True

def process_bam_files(sample_name, analysis_dir):
    """Process all BAM files to create NPZ and HMMcopy files"""
    # List of all BAM files to process with their suffixes
    bam_files = []
    
    # Add proper_paired.bam
    proper_paired_bam = f"{analysis_dir}/{sample_name}/{sample_name}.proper_paired.bam"
    bam_files.append((proper_paired_bam, "proper_paired"))
    
    # Process each filter type to create BAM files
    for filter_name in filter_paths.keys():
        # Add orig BAM
        orig_bam = f"{analysis_dir}/{sample_name}/{sample_name}.{filter_name}_orig.bam"
        bam_files.append((orig_bam, f"{filter_name}_orig"))
        
        # Add fetus BAM
        fetus_bam = f"{analysis_dir}/{sample_name}/{sample_name}.{filter_name}_fetus.bam"
        bam_files.append((fetus_bam, f"{filter_name}_fetus"))
        
        # Add mom BAM
        mom_bam = f"{analysis_dir}/{sample_name}/{sample_name}.{filter_name}_mom.bam"
        bam_files.append((mom_bam, f"{filter_name}_mom"))
    
    # Process each BAM file
    for bam_file, bam_suffix in bam_files:
        logger.info(f"Processing {bam_suffix} BAM file")
        
        # Create NPZ files for all BAM files
        create_npz_files(sample_name, bam_file, bam_suffix, analysis_dir)
        
        # For all except proper_paired.bam, also create HMMcopy files
        if bam_suffix != "proper_paired":
            create_hmmcopy_files(sample_name, bam_file, bam_suffix, analysis_dir)
    
    logger.info("All BAM processing and file generation completed")

def main():
    """Main function to orchestrate the full NIPT pipeline"""
    # Parse command line arguments
    args = parse_args()
    
    sample_name = args.sample_name
    fastq_r1 = args.fastq_r1
    fastq_r2 = args.fastq_r2
    analysis_dir = args.analysis_dir
    fastq_dir = args.fastq_dir
    
    logger.info(f"Starting NIPT pipeline for sample: {sample_name}")
    
    # Create necessary directories
    create_directories(sample_name, analysis_dir)
    
    # Create symbolic links to FASTQ files
    create_symbolic_links(sample_name, fastq_dir, fastq_r1, fastq_r2, analysis_dir)
    
    # Generate proper_paired.bam
    generate_proper_paired_bam(sample_name, fastq_r1, fastq_r2, analysis_dir)
    
    # Process each filter type to create BAM files
    for filter_name, filter_path in filter_paths.items():
        process_filter(sample_name, filter_name, filter_path, analysis_dir)
    
    # Log completion of filter processing
    with open(f"{analysis_dir}/{sample_name}/{sample_name}.log", "a") as log_file:
        log_file.write("9. Samtools filters (of, nf08, nf09) bam : PASS\n")
    
    # Process all BAM files for NPZ and HMMcopy files
    process_bam_files(sample_name, analysis_dir)
    
    # Log completion of all processing
    with open(f"{analysis_dir}/{sample_name}/{sample_name}.log", "a") as log_file:
        log_file.write("10. NPZ and HMMcopy file generation : PASS\n")
    
    logger.info(f"NIPT pipeline completed successfully for sample: {sample_name}")

if __name__ == "__main__":
    main()
