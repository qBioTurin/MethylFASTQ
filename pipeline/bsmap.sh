#!/usr/bin/bash
set -x

path_bsmap=~/bioTool/bsmap-2.90/
path_samtools=~/bioTool/samtools-1.3/

r1=$(realpath $1)
r2=$(realpath $2)
folder_path=$(dirname $(realpath $r1))
ref=$(realpath $3)

data_folder=$folder_path/bsmap
mkdir $data_folder

samfile=$data_folder/alignment.sam
methfile=$data_folder/meth_calling.txt
logfile1=$data_folder/alignment.log
logfile2=$data_folder/calling.log

$path_bsmap/./bsmap -r 0 -s 16 -n 1 -a $r1 -b $r2 -d $ref -o $samfile  > $logfile1
python $path_bsmap/methratio.py -d $ref -s $path_samtools -m 1 -z -i skip -o $methfile $samfile > $logfile2
