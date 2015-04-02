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
organisms=(pestis pseudotuberculosis)
plasmids=(pMT1)
samples=(
    "SRR1749070-0.25x"
    "SRR1749070-0.5x"
    "SRR1749070-1x"
    "SRR1749070-5x"
    "SRR1749070-10x"
    "SRR1749070-15x"
    "SRR1749070-20x"
    "SRR1749070-30x"
    "SRR1749070-0x"
)

declare -A reference
reference["pMT1"]="references/index/AL117211-pMT1"

declare -A gff
gff["pMT1"]="references/AL117211-pMT1.gbk.gff"

# Murine toxin gene
declare -A genes
genes["ymt"]="74612-76210"

# BLASTN parameters
OUTFMT="6 stitle sseqid qseqid qstart qend sstart send evalue bitscore score length pident nident mismatch positive means gapopen ppos qcovs qcovhsp"
NT="/data1/home/groups/readgp/nt/nt"

for o in ${organisms[@]}; do
    for s in ${samples[@]}; do
        for p in ${plasmids[@]}; do
            fq1="sra-controls/${o}/metagenomic/${s}_1.fastq.gz"

            # pestis has 20x, 30x coerages that pseudotuberculosis does not
            if [ -f $fq1 ] ; then
                fq2="sra-controls/${o}/metagenomic/${s}_2.fastq.gz"
                wd="results/mapping/${o}-control/${s}/${p}"
                mkdir -p ${wd}/coverage
                mkdir -p ${wd}/aligned-genes
                mkdir -p ${wd}/aligned-reads

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

                # Plot coverage across complete plasmid, and murine toxin gene
                scripts/mapping/plot-pmt1-murine-toxin-coverage.R ${cov} ${gff[$p]}

                # R plot to PNG did not display correctly so convert the SVG to PNG
                convert -density 300 -resize 2400x1200 ${wd}/coverage/${s}-murine-toxin.svg ${wd}/coverage/${s}-murine-toxin.png

                # Extract reads mapped to murine toxin gene and determine the top 5
                # blast hits against the NT database
                for g in ${!genes[@]}; do
                    fasta=${wd}/aligned-genes/${g}.fasta
                    blastn=${wd}/aligned-genes/${g}.blastn
                    summary=${wd}/aligned-genes/${g}.summary

                    # Extract mapped reads as FASTA
                    samtools view ${bam} "gi|5834685|emb|AL117211.1|":${genes[$g]} | awk '{OFS="\t"; print ">"$1"\n"$10}' > ${fasta}

                    # BLAST against nt
                    blastn -query ${fasta} -outfmt "${OUTFMT}" -db ${NT} -num_threads 20  -max_target_seqs 5 > ${blastn}

                    # Produce summary of hits
                    grep -c "^>" ${fasta} | awk '{print "Total Reads:",$1,"\n"}' > ${summary}
                    echo "Organism Hit Counts" >> ${summary}
                    awk '{print $1,$2}' ${blastn} | sort | uniq -c | sort -rn >> ${summary}
                done
                awk '{print $1,$2}' ${wd}/aligned-genes/*.blastn | sort | uniq -c | sort -rn > ${wd}/aligned-genes/all.summary

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
