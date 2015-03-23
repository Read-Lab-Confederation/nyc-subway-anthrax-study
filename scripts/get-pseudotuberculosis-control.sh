#! /bin/bash
# Download Yersinia pseudotuberculosis WGS reads corresponding to ERR245865.
#
# This project contains 300bp Illumina MiSeq reads and the NYC sample has 100bp
# reads. Use fastx_trimmer to trim the 300bp reads into 100bp reads. Then use
# seqtk to reduce the project to different coverages. Finally concatenate
# these different coverages to sample SRR1749070 (Y. pestis free).
#
# ERR245865: http://www.ebi.ac.uk/ena/data/view/ERR245865
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
CONTROL_DIR=${PROJECT_DIR}/sra-controls/pseudotuberculosis

# Download Bacillus cereus VD142 WGS project
mkdir -p ${CONTROL_DIR}
cd ${CONTROL_DIR}
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR245/ERR245865/ERR245865_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR245/ERR245865/ERR245865_2.fastq.gz

# Trim 300bp MiSeq into 100bp reads
zcat ERR245865_1.fastq.gz | fastx_trimmer -Q 30 -f 1 -l 100 | gzip - > ERR245865_1.100bp.fastq.gz
zcat ERR245865_2.fastq.gz | fastx_trimmer -Q 30 -f 1 -l 100 | gzip - > ERR245865_2.100bp.fastq.gz

# Reduce coverage by randomly subsampling the reads
# Yersinia pseudotuberculosis + plasmid = 4,935,125 ~ 4,935,200

# 0.25x
${PROJECT_DIR}/bin/seqtk sample -s100 ERR245865_1.100bp.fastq.gz 6169 | gzip - > ERR245865_1.100bp.0.25x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 ERR245865_2.100bp.fastq.gz 6169 | gzip - > ERR245865_2.100bp.0.25x.fastq.gz

# 0.5x
${PROJECT_DIR}/bin/seqtk sample -s100 ERR245865_1.100bp.fastq.gz 12338 | gzip - > ERR245865_1.100bp.0.5x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 ERR245865_2.100bp.fastq.gz 12338 | gzip - > ERR245865_2.100bp.0.5x.fastq.gz

# 1x
${PROJECT_DIR}/bin/seqtk sample -s100 ERR245865_1.100bp.fastq.gz 24676 | gzip - > ERR245865_1.100bp.1x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 ERR245865_2.100bp.fastq.gz 24676 | gzip - > ERR245865_2.100bp.1x.fastq.gz

# 5x
${PROJECT_DIR}/bin/seqtk sample -s100 ERR245865_1.100bp.fastq.gz 123380 | gzip - > ERR245865_1.100bp.5x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 ERR245865_2.100bp.fastq.gz 123380 | gzip - > ERR245865_2.100bp.5x.fastq.gz

# 10x
${PROJECT_DIR}/bin/seqtk sample -s100 ERR245865_1.100bp.fastq.gz 246760 | gzip - > ERR245865_1.100bp.10x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 ERR245865_2.100bp.fastq.gz 246760 | gzip - > ERR245865_2.100bp.10x.fastq.gz

# 15x
${PROJECT_DIR}/bin/seqtk sample -s100 ERR245865_1.100bp.fastq.gz 370140 | gzip - > ERR245865_1.100bp.15x.fastq.gz
${PROJECT_DIR}/bin/seqtk sample -s100 ERR245865_2.100bp.fastq.gz 370140 | gzip - > ERR245865_2.100bp.15x.fastq.gz

# Use SRA run SRR1749070 (100bp read lengths, does not contain Y. pestis)
# to add Y. pseudotuberculosis reads to.
mkdir ${CONTROL_DIR}/metagenomic
cd ${CONTROL_DIR}/metagenomic
ln -s ${PROJECT_DIR}/sra-fastq/SRR1749070/SRR1749070_1.fastq.gz ${CONTROL_DIR}/metagenomic/SRR1749070-0x_1.fastq.gz
ln -s ${PROJECT_DIR}/sra-fastq/SRR1749070/SRR1749070_2.fastq.gz ${CONTROL_DIR}/metagenomic/SRR1749070-0x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../ERR245865_1.100bp.0.25x.fastq.gz > SRR1749070-0.25x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../ERR245865_2.100bp.0.25x.fastq.gz > SRR1749070-0.25x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../ERR245865_1.100bp.0.5x.fastq.gz > SRR1749070-0.5x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../ERR245865_2.100bp.0.5x.fastq.gz > SRR1749070-0.5x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../ERR245865_1.100bp.1x.fastq.gz > SRR1749070-1x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../ERR245865_2.100bp.1x.fastq.gz > SRR1749070-1x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../ERR245865_1.100bp.5x.fastq.gz > SRR1749070-5x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../ERR245865_2.100bp.5x.fastq.gz > SRR1749070-5x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../ERR245865_1.100bp.10x.fastq.gz > SRR1749070-10x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../ERR245865_2.100bp.10x.fastq.gz > SRR1749070-10x_2.fastq.gz

cat SRR1749070-0x_1.fastq.gz ../ERR245865_1.100bp.15x.fastq.gz > SRR1749070-15x_1.fastq.gz
cat SRR1749070-0x_2.fastq.gz ../ERR245865_2.100bp.15x.fastq.gz > SRR1749070-15x_2.fastq.gz

cd ${PROJECT_DIR}
