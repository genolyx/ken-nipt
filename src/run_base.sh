#!/bin/bash

if [ $# -ne 5 ]; then
    echo "==============================================================================="
    echo "Usage: $0 <Root_dir> <Working_dir> <Sample_name> <Fastq_R1> <Fastq_R2>"
    echo "bash run_base.sh /home/ken 2504 SAMPLE SAMPLE_1.fq.gz SAMPLE_2.fq.gz"
    echo "==============================================================================="
    exit 1
fi

root_dir=$1
work_dir=$2
sample_name=$3
fastq_r1=$4
fastq_r2=$5

echo "Running: $sample_name"

docker run -v "$root_dir/fastq/$work_dir:/Work/NIPT/fastq" -v "$root_dir/analysis/$work_dir:/Work/NIPT/analysis" -v "$root_dir/log/$work_dir:/Work/NIPT/log" -v "$root_dir/output/$work_dir:/Work/NIPT/output" -v "$root_dir/config:/Work/NIPT/config" -v "$root_dir/refs:/Work/NIPT/refs" --name $sample_name -d nipt_v1.0 perl /Work/NIPT/bin/GC_NIPT_base_v2.5.pl $sample_name $fastq_r1 $fastq_r2 
