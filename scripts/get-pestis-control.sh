#! /bin/bash
# Download Yersinia pestis WGS reads corresponding to SRR1283952.
#
# This project contains 100bp Illumina HiSeq 2000 reads and the NYC sample has
# 100bp reads. Use seqtk to reduce the project to different coverages. Finally
# concatenate these different coverages to sample SRR1749070 (Y. pestis
# free).
#
# SRR1283952: http://www.ebi.ac.uk/ena/data/view/SRR1283952
#
# Program: seqtk (https://github.com/lh3/seqtk)
# Version: commit 43ff625a3211b51f301cb356a34fb8d1e593d50a
#
set -x # Echo all commands
PROJECT_DIR=$(pwd)
if [ -n "$1" ]; then
    PROJECT_DIR=$1
fi
CONTROL_DIR=${PROJECT_DIR}/sra-controls/pestis

# Download Bacillus cereus VD142 WGS project
mkdir -p ${CONTROL_DIR}
cd ${CONTROL_DIR}
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/002/SRR1283952/SRR1283952_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/002/SRR1283952/SRR1283952_2.fastq.gz

# Reduce coverage by randomly subsampling the reads
# Yersinia pestis CO92 + plasmids ~ 4,829,855 ~ 4,830,000

# 0.25x
${PROJECT_DIR}/bin/seqtk sample -s100 SRR1283952_1.fastq.gz 6038 | gzip - > SRR1283952_1.0.25x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 SRR1283952_2.fastq.gz 6038 | gzip - > SRR1283952_2.0.25x.fastq.gz

# 0.5x
${PROJECT_DIR}/bin/seqtk sample -s100 SRR1283952_1.fastq.gz 12075 | gzip - > SRR1283952_1.0.5x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 SRR1283952_2.fastq.gz 12075 | gzip - > SRR1283952_2.0.5x.fastq.gz

# 1x
${PROJECT_DIR}/bin/seqtk sample -s100 SRR1283952_1.fastq.gz 24150 | gzip - > SRR1283952_1.1x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 SRR1283952_2.fastq.gz 24150 | gzip - > SRR1283952_2.1x.fastq.gz

# 5x
${PROJECT_DIR}/bin/seqtk sample -s100 SRR1283952_1.fastq.gz 120750 | gzip - > SRR1283952_1.5x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 SRR1283952_2.fastq.gz 120750 | gzip - > SRR1283952_2.5x.fastq.gz

# 10x
${PROJECT_DIR}/bin/seqtk sample -s100 SRR1283952_1.fastq.gz 241500 | gzip - > SRR1283952_1.10x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 SRR1283952_2.fastq.gz 241500 | gzip - > SRR1283952_2.10x.fastq.gz

# 15x
${PROJECT_DIR}/bin/seqtk sample -s100 SRR1283952_1.fastq.gz 362250 | gzip - > SRR1283952_1.15x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 SRR1283952_2.fastq.gz 362250 | gzip - > SRR1283952_2.15x.fastq.gz

# 20x
${PROJECT_DIR}/bin/seqtk sample -s100 SRR1283952_1.fastq.gz 483000 | gzip - > SRR1283952_1.20x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 SRR1283952_2.fastq.gz 483000 | gzip - > SRR1283952_2.20x.fastq.gz

# 30x
${PROJECT_DIR}/bin/seqtk sample -s100 SRR1283952_1.fastq.gz 724500 | gzip - > SRR1283952_1.30x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 SRR1283952_2.fastq.gz 724500 | gzip - > SRR1283952_2.30x.fastq.gz

# Use SRA run SRR1749070 (100bp read lengths, does not contain Y. pestis)
# to add Y. pestis reads to.
mkdir ${CONTROL_DIR}/metagenomic
cd ${CONTROL_DIR}/metagenomic
ln -s ${PROJECT_DIR}/sra-fastq/SRR1749070/SRR1749070_1.fastq.gz ${CONTROL_DIR}/metagenomic/SRR1749070-0x_1.fastq.gz
ln -s ${PROJECT_DIR}/sra-fastq/SRR1749070/SRR1749070_2.fastq.gz ${CONTROL_DIR}/metagenomic/SRR1749070-0x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../SRR1283952_1.0.25x.fastq.gz > SRR1749070-0.25x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../SRR1283952_2.0.25x.fastq.gz > SRR1749070-0.25x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../SRR1283952_1.0.5x.fastq.gz > SRR1749070-0.5x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../SRR1283952_2.0.5x.fastq.gz > SRR1749070-0.5x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../SRR1283952_1.1x.fastq.gz > SRR1749070-1x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../SRR1283952_2.1x.fastq.gz > SRR1749070-1x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../SRR1283952_1.5x.fastq.gz > SRR1749070-5x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../SRR1283952_2.5x.fastq.gz > SRR1749070-5x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../SRR1283952_1.10x.fastq.gz > SRR1749070-10x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../SRR1283952_2.10x.fastq.gz > SRR1749070-10x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../SRR1283952_1.15x.fastq.gz > SRR1749070-15x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../SRR1283952_2.15x.fastq.gz > SRR1749070-15x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../SRR1283952_1.20x.fastq.gz > SRR1749070-20x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../SRR1283952_2.20x.fastq.gz > SRR1749070-20x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../SRR1283952_1.30x.fastq.gz > SRR1749070-30x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../SRR1283952_2.30x.fastq.gz > SRR1749070-30x_2.fastq.gz

cd ${PROJECT_DIR}
