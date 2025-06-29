#!/bin/bash

image=$1

docker run -it --rm -e TZ=Asia/Seoul --entrypoint bash $image
