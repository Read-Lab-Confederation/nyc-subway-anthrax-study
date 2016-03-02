#! /bin/bash

while read line
do
    echo "Working on ${line}"
    mkdir -p sra-fastq/${line}
    fastq-dump --split-files -O sra-fastq/${line} --gzip /data1/home/rpetit/ncbi/sra/${line}.sra
done < $1






