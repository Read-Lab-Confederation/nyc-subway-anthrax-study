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
#
# Program: fastq_to_fasta - Part of FASTX Toolkit
# Version: 0.0.13.2
#
# Program: convert - convert SVG to PNG (ImageMagick)
# Version: 6.6.9-7 2014-03-06 Q16
#
# Program: blastn
# Version: 2.2.30+
#
set -x # Echo all commands
genomes=(anthracis cereus)

declare -A samples
samples["SRR1748707"]="P00134"
samples["SRR1748708"]="P00134"
samples["SRR1749083"]="P00497"

declare -A reference
reference["anthracis"]="references/index/CP009541-anthracis"
reference["cereus"]="references/index/NC_003909-cereus"

# BLASTN parameters
OUTFMT="6 stitle sseqid qseqid qstart qend sstart send evalue bitscore score length pident nident mismatch positive means gapopen ppos qcovs qcovhsp"
NT="/data1/home/groups/readgp/nt/nt"

for s in ${!samples[@]}; do
    for p in ${genomes[@]}; do
        wd="results/mapping/bacilli-nyc/${s}/${p}"
        mkdir -p ${wd}/coverage
        mkdir -p ${wd}/aligned-reads

        fq1="sra-pathogens/anthracis/${samples[$s]}/${s}_1.fastq.gz"
        fq2="sra-pathogens/anthracis/${samples[$s]}/${s}_2.fastq.gz"
        sam="${wd}/${s}.sam"
        bam="${wd}/${s}.bam"
        cov="${wd}/coverage/${s}.coverage.gz"

        # Align using BWA, sort sam to bam and index bam
        bin/bwa mem -t 20 ${reference[$p]} ${fq1} ${fq2} > ${sam}
        samtools view -@ 10 -bS ${sam} | samtools sort -@ 10 - ${wd}/${s}
        samtools index -b ${bam}
        rm ${sam}

        # Use genomeCoverageBed to get the coverage for each position and plot
        # the coverage for differening sliding windows with 0.5 overlap
        genomeCoverageBed -d -ibam ${bam} | gzip --best - > ${cov}
        scripts/mapping/plot-coverage.R ${cov} 0.5 ${p}

        # Extract aligned reads using bam2fastq and convert to fasta
        ofq="${wd}/aligned-reads/${samples[$s]}_${s}#.fastq"
        bin/bam2fastq -o ${ofq} --no-unaligned ${bam}
        cat ${wd}/aligned-reads/*.fastq | fastq_to_fasta -Q33 -n | gzip --best - > ${wd}/aligned-reads/${samples[$s]}_${s}.fasta.gz
        gzip --best "${wd}/aligned-reads/${samples[$s]}_${s}_1.fastq"
        gzip --best "${wd}/aligned-reads/${samples[$s]}_${s}_2.fastq"
    done
done
