#!/bin/bash

if [ $# -ne 3 ]; then
    echo "==============================================================================="
    echo "Usage: $0 <Root_dir> <Working_dir> <base_image> <new_container_name>"
    echo "bash make_container_with_hostdir2.sh /data2/gc_analysis/dev_connect nipt_v1.7 nipt_v1.7_REF_dev"
    echo "==============================================================================="
    exit 1
fi

root_dir=$1
image=$2
container=$3

docker run -v $root_dir/input_data:/Work/NIPT/fastq -v $root_dir/analysis:/Work/NIPT/analysis -v $root_dir/log:/Work/NIPT/log -v $root_dir/output_data:/Work/NIPT/output -v $root_dir/config:/Work/NIPT/config -v $root_dir/refs:/Work/NIPT/refs --name $container -dt $image /bin/bash
