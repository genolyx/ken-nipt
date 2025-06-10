# ---------------------------------------#
# Utilities
# ---------------------------------------#

1. cleanup_bam.sh 
    : After the analysis completion, remove various bams to save disk
    : Iterating subdirectories under the target_dir, remove bams listed in the script.

    ./cleanup_bam.sh /home/ken/ken-nipt/analysis/2503

2. get_threads.sh 
    : Finding the optimal parameter for bwa-mem2 and samtools

3. go_container.sh
    : Going into Docker container
    : docker exec -it $container /bin/bas

4. make_archives.sh
    : Iterating subdirectories, it compresses files excluding "bam" and "bai"

5. make_container_with_hostdir.sh
    : Making docker container with "root_dir", "work_dir"

6. make_container_with_hostdir_no_workdir.sh
    : Making docker container only with "root_dir"

7. make_input_dir.sh
    : Finding R1, R2 fastq files, Making directories and move fastq files to sampel directory
    : It needs the sample_list.txt
    <sample_list.txt>
        2504280008
        2504300001
        2504300002
        2504300003
        2504300004
        2504300006

8. make_list_base.py
    : make the sample sheet for nipt_v1.0 docker run
    <sample_sheet.tsv>
        SAMPLE_NAME	FQ1	FQ2	LAB
        2504300012	2504300012_S42_R1_001.fastq.gz	2504300012_S42_R2_001.fastq.gz	cordlife
        2504300013	2504300013_S33_R1_001.fastq.gz	2504300013_S33_R2_001.fastq.gz	cordlife
        2504300014	2504300014_S28_R1_001.fastq.gz	2504300014_S28_R2_001.fastq.gz	cordlife
        2504300015	2504300015_S23_R1_001.fastq.gz	2504300015_S23_R2_001.fastq.gz	cordlife

9. make_list.py
    : same as make_list_base.py but it adds "age" column with default value.
    
10. move_fastq.sh
    : Reading sample_name, fq1, fq2 from sample_list.csv, it create sample_name directory and move fastq files into it. This is very similar to "7. make_input_dir.sh"

11. resource_logger.sh
    : Monitor resources
    : bash resource_logger.sh 10 rsc_log &

        timestamp	cpu_usage_percent	mem_used_mb	mem_total_mb	memory_usage_percent
        2025-06-08 16:08:19	0	1538	201473	0.8
        2025-06-08 16:08:30	0	1537	201473	0.8

12. build_docker.sh
    : build docker image

    build_docker.sh -n nipt_docker_v1.0 -t latest

# ---------------------------------------#
# Sample Run script
# ---------------------------------------#

1. run_nipt_dev.sh
    : run development pipeline using nipt_docker_dev
    : bash run_nipt_dev.sh -s 2504300007 -1 2504300007_S8_R1_001.fastq.gz -2 2504300007_S8_R1_001.fastq.gz -l cordlife -a 30 -root /home/ken/ken-nipt -work 250430_01

2. run_nipt_dev_noage.sh
    : Sampe run without knwoing sample's age
    in docker run, it has default age as "30"
    --fastq_r2 "$FASTQ_R2" \
    --labcode "$LABCODE" \
    --age 3

    : related to "run_batch_dev_noage.py"

3. run_nipt_dev_ucl.sh

    : bash run_nipt_dev_ucl.sh -s DNP250200018 -1 DNP250200018_data_file_DNP250200018_4122400608_D19_S67_R1_001.fastq.gz -2 DNP250200018_data_file_DNP250200018_4122400608_D19_S67_R2_001.fastq.gz -l ucl -a 30 -root /home/ken/ken-nipt -work 2502

    : When we're going to use the separated analysis environment based on labcode
        HOST_FASTQ_DIR="$ROOT_DIR/ucl/fastq/$WORK_DIR"
        HOST_ANALYSIS_DIR="$ROOT_DIR/ucl/analysis/$WORK_DIR"
        HOST_OUTPUT_DIR="$ROOT_DIR/ucl/output/$WORK_DIR"
        HOST_LOG_DIR="$ROOT_DIR/ucl/log"
        HOST_DATA_DIR="$ROOT_DIR/data"
        HOST_CONFIG_DIR="$ROOT_DIR/config"

        --> config, data are shared

4. run_nipt.sh
    : run commercial pipeline using "nipt_docker_v1.0"

5. run_base.sh
    run the part1 of old pipeline 

# ---------------------------------------#
# Batch Run script
# ---------------------------------------#

1. run_batch_base.py

    : run batch run with run_base.sh (part1)
    python3 run_batch_base.py 5 ~/fastq/250605_01/sample_sheet.tsv "sh run_base.sh /home/ken 250605_01 SAMPLE_NAME FQ1 FQ2"

3. run_batch_dev_noage.py
    : Batch run without knowing sample's age
    bash /home/ken/ken-nipt/src/run_nipt_dev_noage.sh -s SAMPLE_NAME -1 FQ1 -2 FQ2 -l LABCODE -root /home/ken/ken-nipt -work WORK_DIR --detached
    
    in docker run, it has default age as "30"
    --fastq_r2 "$FASTQ_R2" \
    --labcode "$LABCODE" \
    --age 30

    : python3 run_batch_dev_noage.py -l cordlife -s ~/ken-nipt/src/test/cordlife_ref_female_fq_list.tsv -m 7 -i 10 --log_stdout Y

4. run_batch_dev.py
        : Using run_nipt_dev.sh, run batch 
        #command = "bash /home/ken/ken-nipt/src/run_nipt_dev.sh -s SAMPLE_NAME -1 FQ1 -2 FQ2 -l LABCODE -root /home/ken/ken-nipt -work WORK_DIR --detached"
        command = "bash /home/ken/ken-nipt/src/run_nipt_dev.sh -s SAMPLE_NAME -1 FQ1 -2 FQ2 -a AGE -l LABCODE -root /home/ken/ken-nipt -work WORK_DIR --detached"

    : python3 run_batch_dev.py -l cordlife -s ~/ken-nipt/src/test/250605_01_list.tsv -m 5 -i 10 --log_stdout Y
    : python3 run_batch_dev.py -l cordlife -s ~/ken-nipt/src/test/cordlife_ref_female_fq_list.tsv -m 7 -i 10 --log_stdout Y

    <cordlife_ref_female_fq_list.tsv>
        SAMPLE_NAME	WORK_DIR	FQ1	FQ2
    OPC241000004	2410	OPC241000004_data_file_2410250014_S7_R1_001.fastq.gz	OPC241000004_data_file_2410250014_S7_R2_001.fastq.gz
    OPC241000007	2410	OPC241000007_data_file_2410250007_S12_R1_001.fastq.gz	OPC241000007_data_file_2410250007_S12_R2_001.fastq.gz

5. run_batch.py

    : run commercial pipeline (using run_nipt.sh)

    : bash run_nipt.sh -s 2506030014 -1 2506030014_S1_R1_001.fastq.gz -2 2506030014_S1_R2_001.fastq.gz -l cordlife -a 30 -root /home/ken/ken-nipt -work 250605_01

# ---------------------------------------#
# Directories
# ---------------------------------------#

1. logs
    : bash script run log

2. rsc_log
    : resource_logger.sh log

3. test
    : samples list files for testing
