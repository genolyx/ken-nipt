#!/bin/bash

image=$1

#docker run -it --rm -e TZ=Asia/Seoul --entrypoint bash $image
docker run -it --rm \
  -e TZ=Asia/Seoul \
  -v /home/ken/ken-nipt/data:/Work/NIPT/data \
  -v /home/ken/ken-nipt/analysis/2509:/Work/NIPT/analysis \
  --entrypoint bash $image
