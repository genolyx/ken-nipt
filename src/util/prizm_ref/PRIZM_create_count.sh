#!/bin/bash

for bam in `ls *.uniq.good.bam`
do
    file_name="${bam%.of_orig.bam}"
    echo $bam

    samtools idxstats $bam | head -n 25 | tail -n 24 | cut -f 1,3 > $file_name.count.txt
    /Work/BIO/bin/HMMcopy/bin/readCounter -w 10000000 -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY $bam > $file_name.10mb.wig
    Rscript /Work/NICE/NICE_V4/bin/HMMcopy.R $file_name.10mb.wig /Work/BIO/bin/HMMcopy/hg19.gc.10mb.wig  /Work/BIO/bin/HMMcopy/hg19.map.10mb.wig
    sed -e 's/NA/0/g' $file_name.10mb.wig.Normalization.txt > $file_name.10mb.txt
    rm $file_name.10mb.wig.Normalization.txt
done
