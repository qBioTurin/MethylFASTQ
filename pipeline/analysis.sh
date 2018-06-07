#!/usr/bin/bash
set -x

path_r1=$(realpath $1)
path_r2=$(realpath $2)
path_ref=$(realpath $3)
path_index=$(realpath $4)

nohup ./bsmap.sh $path_r1 $path_r2 $path_ref &> /dev/null &
nohup ./bismark.sh $path_r1 $path_r2 $path_index &> /dev/null &
