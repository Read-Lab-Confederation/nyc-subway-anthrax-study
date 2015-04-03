#! /bin/bash
# Download Bacillus cereus VD142 WGS reads corresponding to SRR642775.
#
# This project contains 100bp Illumina HiSeq 2000 reads and the NYC sample has
# 100bp reads. Use seqtk to reduce the project to different coverages. Finally
# concatenate these different coverages to sample SRR1749070 (B. anthracis
# free).
#
# SRR642775: http://www.ebi.ac.uk/ena/data/view/SRR642775
#
# Program: seqtk (https://github.com/lh3/seqtk)
# Version: commit 43ff625a3211b51f301cb356a34fb8d1e593d50a
#
set -x # Echo all commands
PROJECT_DIR=$(pwd)
if [ -n "$1" ]; then
    PROJECT_DIR=$1
fi
CONTROL_DIR=${PROJECT_DIR}/sra-controls/cereus

# Download Bacillus cereus VD142 WGS project
mkdir -p ${CONTROL_DIR}
cd ${CONTROL_DIR}
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR642/SRR642775/SRR642775_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR642/SRR642775/SRR642775_2.fastq.gz

# Reduce coverage by randomly subsampling the reads
# Bacillus cereus VD142 = 5,923,913 ~ 5,924,000

# 0.25x
${PROJECT_DIR}/bin/seqtk sample -s100 SRR642775_1.fastq.gz 7405 | gzip - > SRR642775_1.0.25x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 SRR642775_2.fastq.gz 7405 | gzip - > SRR642775_2.0.25x.fastq.gz

# 0.5x
${PROJECT_DIR}/bin/seqtk sample -s100 SRR642775_1.fastq.gz 14810 | gzip - > SRR642775_1.0.5x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 SRR642775_2.fastq.gz 14810 | gzip - > SRR642775_2.0.5x.fastq.gz

# 1x
${PROJECT_DIR}/bin/seqtk sample -s100 SRR642775_1.fastq.gz 29620 | gzip - > SRR642775_1.1x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 SRR642775_2.fastq.gz 29620 | gzip - > SRR642775_2.1x.fastq.gz

# 5x
${PROJECT_DIR}/bin/seqtk sample -s100 SRR642775_1.fastq.gz 148100 | gzip - > SRR642775_1.5x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 SRR642775_2.fastq.gz 148100 | gzip - > SRR642775_2.5x.fastq.gz

# Use SRA run SRR1749070 (100bp read lengths, does not contain B. anthracis)
# to add B. cereus VD142 reads to.
mkdir ${CONTROL_DIR}/metagenomic
cd ${CONTROL_DIR}/metagenomic
ln -s ${PROJECT_DIR}/sra-fastq/SRR1749070/SRR1749070_1.fastq.gz ${CONTROL_DIR}/metagenomic/SRR1749070-0x_1.fastq.gz
ln -s ${PROJECT_DIR}/sra-fastq/SRR1749070/SRR1749070_2.fastq.gz ${CONTROL_DIR}/metagenomic/SRR1749070-0x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../SRR642775_1.0.25x.fastq.gz > SRR1749070-0.25x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../SRR642775_2.0.25x.fastq.gz > SRR1749070-0.25x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../SRR642775_1.0.5x.fastq.gz > SRR1749070-0.5x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../SRR642775_2.0.5x.fastq.gz > SRR1749070-0.5x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../SRR642775_1.1x.fastq.gz > SRR1749070-1x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../SRR642775_2.1x.fastq.gz > SRR1749070-1x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../SRR642775_1.5x.fastq.gz > SRR1749070-5x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../SRR642775_2.5x.fastq.gz > SRR1749070-5x_2.fastq.gz

cd ${PROJECT_DIR}
