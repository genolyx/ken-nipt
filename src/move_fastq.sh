#!/bin/bash

while IFS=$'\t' read -r sample_name fq1 fq2; do
  echo $sample_name

  mkdir -p $sample_name
  mv $fq1 $fq2 $sample_name
done < sample_list.csv
