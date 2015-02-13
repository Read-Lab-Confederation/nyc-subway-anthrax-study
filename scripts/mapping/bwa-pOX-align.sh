#! /bin/bash
# Program: bwa (alignment via Burrows-Wheeler transformation)
# Version: 0.7.5a-r405
#
# Program: samtools (Tools for alignments in the SAM format)
# Version: 1.1 (using htslib 1.1)

# Map to pXO1
mkdir results/pXO1
bin/bwa mem -t 20 references/index/CP009540_pXO1 sra-pathogens/anthracis/P00134/SRR1748707_1.fastq.gz sra-pathogens/anthracis/P00134/SRR1748707_2.fastq.gz > results/pXO1/P00134_SRR1748707.sam
samtools view -bS results/pXO1/P00134_SRR1748707.sam | samtools sort - results/pXO1/P00134_SRR1748707.sorted
samtools index -b results/pXO1/P00134_SRR1748707.sorted.bam

bin/bwa mem -t 20 references/index/CP009540_pXO1 sra-pathogens/anthracis/P00134/SRR1748708_1.fastq.gz sra-pathogens/anthracis/P00134/SRR1748708_2.fastq.gz > results/pXO1/P00134_SRR1748708.sam
samtools view -bS results/pXO1/P00134_SRR1748708.sam | samtools sort - results/pXO1/P00134_SRR1748708.sorted
samtools index -b results/pXO1/P00134_SRR1748708.sorted.bam

bin/bwa mem -t 20 references/index/CP009540_pXO1 sra-pathogens/anthracis/P00497/SRR1749083_1.fastq.gz sra-pathogens/anthracis/P00497/SRR1749083_2.fastq.gz > results/pXO1/P00497.sam
samtools view -bS results/pXO1/P00497.sam | samtools sort - results/pXO1/P00497.sorted
samtools index -b results/pXO1/P00497.sorted.bam


# Map to pXO2
mkdir results/pXO1
bin/bwa mem -t 20 references/index/NC_007323_pXO2 sra-pathogens/anthracis/P00134/SRR1748707_1.fastq.gz sra-pathogens/anthracis/P00134/SRR1748707_2.fastq.gz > results/pXO2/P00134_SRR1748707.sam
samtools view -bS results/pXO2/P00134_SRR1748707.sam | samtools sort - results/pXO2/P00134_SRR1748707.sorted
samtools index -b results/pXO2/P00134_SRR1748707.sorted.bam

bin/bwa mem -t 20 references/index/NC_007323_pXO2 sra-pathogens/anthracis/P00134/SRR1748708_1.fastq.gz sra-pathogens/anthracis/P00134/SRR1748708_2.fastq.gz > results/pXO2/P00134_SRR1748708.sam
samtools view -bS results/pXO2/P00134_SRR1748708.sam | samtools sort - results/pXO2/P00134_SRR1748708.sorted
samtools index -b results/pXO2/P00134_SRR1748708.sorted.bam

bin/bwa mem -t 20 references/index/NC_007323_pXO2 sra-pathogens/anthracis/P00497/SRR1749083_1.fastq.gz sra-pathogens/anthracis/P00497/SRR1749083_2.fastq.gz > results/pXO2/P00497.sam
samtools view -bS results/pXO2/P00497.sam | samtools sort - results/pXO2/P00497.sorted
samtools index -b results/pXO2/P00497.sorted.bam
