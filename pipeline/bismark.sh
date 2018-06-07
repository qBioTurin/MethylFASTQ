#!/usr/bin/bash
set -x

path_bismark=~/bioTool/bismark_v0.19.0/
path_bowtie2=~/bioTool/bowtie2-2.2.7/
path_samtools=~/bioTool/samtools-1.3/
temp=~/Data/MethylSeq/.temp/d/

r1=$(realpath $1)
r2=$(realpath $2)
folder_path=$(dirname $(realpath $r1))
ref=$(realpath $3)

data_folder=$folder_path/bismark
mkdir $data_folder

logfile1=$data_folder/alignment.log
logfile2=$data_folder/calling.log

$path_bismark/./bismark --non_directional --bowtie2 --path_to_bowtie $path_bowtie2 \
                        --genome $ref --temp $temp -p 4 -o $data_folder -1 $r1 -2 $r2 > $logfile1

gunzip $data_folder/*.gz

$path_bismark/./bismark_methylation_extractor --comprehensive --merge_non_CpG --bedGraph --buffer_size 10G \
                        --cytosine_report --report --samtools_path $path_samtools \
                        --genome_folder $ref -o $data_folder $data_folder/*.sam > $logfile2
