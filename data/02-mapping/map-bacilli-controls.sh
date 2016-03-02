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
organisms=(anthracis cereus)
genomes=(anthracis cereus)
samples=(
    "SRR1749070-0.01x"
    "SRR1749070-0.05x"
    "SRR1749070-0.10x"
    "SRR1749070-0.25x"
    "SRR1749070-0.5x"
    "SRR1749070-1x"
    "SRR1749070-5x"
    "SRR1749070-0x"
)

declare -A reference
reference["anthracis"]="references/index/CP009541-anthracis"
reference["cereus"]="references/index/NC_003909-cereus"

# BLASTN parameters
OUTFMT="6 stitle sseqid qseqid qstart qend sstart send evalue bitscore score length pident nident mismatch positive means gapopen ppos qcovs qcovhsp"
NT="/data1/home/groups/readgp/nt/nt"

for o in ${organisms[@]}; do
    for s in ${samples[@]}; do
        for p in ${genomes[@]}; do
            fq1="sra-controls/${o}/metagenomic/${s}_1.fastq.gz"

            # anthracis has 0.01x, 0.05x, 0.10x coerages that cereus does not
            if [ -f $fq1 ] ; then
                wd="results/mapping/bacilli-controls/${o}-control/${s}/${p}"
                mkdir -p ${wd}/coverage
                mkdir -p ${wd}/aligned-genes
                mkdir -p ${wd}/aligned-reads

                fq2="sra-controls/${o}/metagenomic/${s}_2.fastq.gz"
                sam="${wd}/${s}.sam"
                bam="${wd}/${s}.bam"
                cov="${wd}/coverage/${s}.coverage.gz"

                # Align using BWA, sort sam to bam and index bam, then remove sam
                bin/bwa mem -t 20 ${reference[$p]} ${fq1} ${fq2} > ${sam}
                samtools view -@ 10 -bS ${sam} | samtools sort -@ 10 - ${wd}/${s}
                samtools index -b ${bam}
                rm ${sam}

                # Use genomeCoverageBed to get the coverage for each position and plot
                # the coverage for differening sliding windows with 0.5 overlap
                genomeCoverageBed -d -ibam ${bam} | gzip --best - > ${cov}
                scripts/mapping/plot-coverage.R ${cov} 0.5 ${p}

                # Extract aligned reads using bam2fastq and convert to fasta
                ofq="${wd}/aligned-reads/${s}#.fastq"
                bin/bam2fastq -o ${ofq} --no-unaligned ${bam}
                cat ${wd}/aligned-reads/*.fastq | fastq_to_fasta -Q33 -n | gzip --best - > ${wd}/aligned-reads/${s}.fasta.gz
                gzip --best "${wd}/aligned-reads/${s}_1.fastq"
                gzip --best "${wd}/aligned-reads/${s}_2.fastq"

            fi
        done
    done
done
