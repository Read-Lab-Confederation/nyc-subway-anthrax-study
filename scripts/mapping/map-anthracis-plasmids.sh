#! /bin/bash
# Program: bwa (alignment via Burrows-Wheeler transformation)
# Version: 0.7.5a-r405
#
# Program: samtools (Tools for alignments in the SAM format)
# Version: 1.1 (using htslib 1.1)
#
# Program: bam2fastq - extract sequences from a BAM file
# Version: v1.1.0
#
# Program: bedtools genomecov (aka genomeCoverageBed)
# Version: v2.16.2
# Summary: Compute the coverage of a feature file among a genome.

set -x # Echo all commands

plasmids=(pXO1 pXO2)

declare -A samples
samples["SRR1748707"]="P00134"
samples["SRR1748708"]="P00134"
samples["SRR1749083"]="P00497"

declare -A reference
reference["pXO1"]="references/index/CP009540_pXO1"
reference["pXO2"]="references/index/NC_007323_pXO2"

for p in ${plasmids[@]}; do
    wd="results/${p}"
    mkdir -p ${wd}/coverage
    mkdir -p ${wd}/aligned-reads

    for s in ${!samples[@]}; do
        fq1="sra-pathogens/anthracis/${samples[$s]}/${s}_1.fastq.gz"
        fq2="sra-pathogens/anthracis/${samples[$s]}/${s}_2.fastq.gz"
        sam="${wd}/${samples[$s]}_${s}.sam"
        bam="${wd}/${samples[$s]}_${s}.bam"
        cov="${wd}/coverage/${samples[$s]}_${s}.coverage"

        # Align using BWA, sort sam to bam and index bam
        bin/bwa mem -t 20 ${reference[$p]} ${fq1} ${fq2} > ${sam}
        samtools view -bS ${sam} | samtools sort - ${wd}/${samples[$s]}_${s}
        samtools index -b ${bam}

        # Use genomeCoverageBed to get the coverage for each position and plot
        # the coverage for differening sliding windows with 0.5 overlap
        genomeCoverageBed -d -ibam ${bam} > ${cov}
        scripts/mapping/plot-coverage.R ${cov} 0.5 ${p}

        # Extract aligned reads using bam2fastq
        ofq="${wd}/aligned-reads/${samples[$s]}/${s}_#.fastq"
        bin/bam2fastq -o ${ofq} --no-unaligned ${bam}
    done
done
