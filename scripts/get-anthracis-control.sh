#! /bin/bash
# Download Bacillus anthracis WGS reads corresponding to DRR014739.
#
# This project contains 300bp Illumina MiSeq reads and the NYC sample has 100bp
# reads. Use fastx_trimmer to trim the 300bp reads into 100bp reads. Then use
# seqtk to reduce the project to different coverages. Finally concatenate
# these different coverages to sample SRR1749070 (B. anthracis free).
#
# DRR014739: http://www.ebi.ac.uk/ena/data/view/DRR014739
#
# Program: fastq_trimmer - Part of FASTX Toolkit
# Version: 0.0.13.2
#
# Program: seqtk (https://github.com/lh3/seqtk)
# Version: commit 43ff625a3211b51f301cb356a34fb8d1e593d50a
#
set -x # Echo all commands
PROJECT_DIR=$(pwd)
if [ -n "$1" ]; then
    PROJECT_DIR=$1
fi
CONTROL_DIR=${PROJECT_DIR}/sra-controls/anthracis

# Download Bacillus anthracis WGS project
mkdir -p ${CONTROL_DIR}
cd ${CONTROL_DIR}
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR014/DRR014739/DRR014739_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR014/DRR014739/DRR014739_2.fastq.gz

# Trim 300bp MiSeq into 100bp reads
zcat DRR014739_1.fastq.gz | fastx_trimmer -Q 30 -f 1 -l 100 | gzip - > DRR014739_1.100bp.fastq.gz
zcat DRR014739_2.fastq.gz | fastx_trimmer -Q 30 -f 1 -l 100 | gzip - > DRR014739_2.100bp.fastq.gz

# Reduce coverage by randomly subsampling the reads
# Bacillus anthracis + PXO1 + PXO2 = 5,503,799 ~ 5,503,800
# 0.25x
${PROJECT_DIR}/bin/seqtk sample -s100 DRR014739_1.100bp.fastq.gz 6880 | gzip - > DRR014739_1.100bp.0.25x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 DRR014739_2.100bp.fastq.gz 6880 | gzip - > DRR014739_2.100bp.0.25x.fastq.gz

# 0.5x
${PROJECT_DIR}/bin/seqtk sample -s100 DRR014739_1.100bp.fastq.gz 13760 | gzip - > DRR014739_1.100bp.0.5x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 DRR014739_2.100bp.fastq.gz 13760 | gzip - > DRR014739_2.100bp.0.5x.fastq.gz

# 1x
${PROJECT_DIR}/bin/seqtk sample -s100 DRR014739_1.100bp.fastq.gz 27519 | gzip - > DRR014739_1.100bp.1x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 DRR014739_2.100bp.fastq.gz 27519 | gzip - > DRR014739_2.100bp.1x.fastq.gz

# 5x
${PROJECT_DIR}/bin/seqtk sample -s100 DRR014739_1.100bp.fastq.gz 137595 | gzip - > DRR014739_1.100bp.5x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 DRR014739_2.100bp.fastq.gz 137595 | gzip - > DRR014739_2.100bp.5x.fastq.gz

# Use SRA run SRR1749070 (100bp read lengths, does not contain B. anthracis)
# to add B. anthracis reads to.
mkdir ${CONTROL_DIR}/metagenomic
cd ${CONTROL_DIR}/metagenomic
ln -s ${PROJECT_DIR}/sra-fastq/SRR1749070/SRR1749070_1.fastq.gz ${CONTROL_DIR}/metagenomic/SRR1749070-0x_1.fastq.gz
ln -s ${PROJECT_DIR}/sra-fastq/SRR1749070/SRR1749070_2.fastq.gz ${CONTROL_DIR}/metagenomic/SRR1749070-0x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../DRR014739_1.100bp.0.25x.fastq.gz > SRR1749070-0.25x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../DRR014739_2.100bp.0.25x.fastq.gz > SRR1749070-0.25x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../DRR014739_1.100bp.0.5x.fastq.gz > SRR1749070-0.5x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../DRR014739_2.100bp.0.5x.fastq.gz > SRR1749070-0.5x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../DRR014739_1.100bp.1x.fastq.gz > SRR1749070-1x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../DRR014739_2.100bp.1x.fastq.gz > SRR1749070-1x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../DRR014739_1.100bp.5x.fastq.gz > SRR1749070-5x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../DRR014739_2.100bp.5x.fastq.gz > SRR1749070-5x_2.fastq.gz

cd ${PROJECT_DIR}
