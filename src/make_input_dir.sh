#!/bin/bash

if [ $# -ne 2 ]; then
    echo "==============================================================================="
    echo "Usage: $0 <Root_dir> <Working_dir> 
    echo "bash make_input_dir.sh /data2/gc_analysis/connect 2505 
    echo "==============================================================================="
    exit 1
fi

root_dir=$1
work_dir=$2

while read -r sample_name; do
  echo "Processing $sample_name ..."

  mkdir -p "$root_dir/fastq/$work_dir/$sample_name"

  find $root_dir/fastq/$work_dir -maxdepth 1 -type f -name ${sample_name}*_R1_001.fastq.gz -exec mv {} $root_dir/fastq/$work_dir/$sample_name/ \;
  find $root_dir/fastq/$work_dir -maxdepth 1 -type f -name ${sample_name}*_R2_001.fastq.gz -exec mv {} $root_dir/fastq/$work_dir/$sample_name/ \;

done < $root_dir/fastq/$work_dir/sample_list.txt
