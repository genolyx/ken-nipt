#!/bin/bash

if [ $# -ne 4 ]; then
    echo "==============================================================================="
    echo "Usage: $0 <Root_dir> <Working_dir> <base_image> <new_container_name>"
    echo "bash make_container.sh /data2/ken/NIPT 240802 nipt_v0.95 nipt_v0.95_dev"
    echo "==============================================================================="
    exit 1
fi

root_dir=$1
work_dir=$2
image=$3
container=$4

docker run -v $root_dir/fastq/$work_dir:/Work/NIPT/fastq -v $root_dir/analysis/$work_dir:/Work/NIPT/analysis -v $root_dir/log/$work_dir:/Work/NIPT/log -v $root_dir/output/$work_dir:/Work/NIPT/output -v $root_dir/config:/Work/NIPT/config -v $root_dir/refs:/Work/NIPT/refs --name $container -dt $image /bin/bash
