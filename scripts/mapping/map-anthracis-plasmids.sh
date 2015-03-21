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
# Program: blastn - Part of FASTX Toolkit
# Version: 2.2.30+
#
set -x # Echo all commands

plasmids=(pXO1 pXO2)

declare -A samples
samples["SRR1748707"]="P00134"
samples["SRR1748708"]="P00134"
samples["SRR1749083"]="P00497"

declare -A reference
reference["pXO1"]="references/index/CP009540_pXO1"
reference["pXO2"]="references/index/NC_007323_pXO2"

declare -A gff
gff["pXO1"]="references/CP009540_pXO1.gbk.gff"

# Genes associated with anthrax toxin and thier positions on PXO1
declare -A genes
genes["cya"]="50445-52847"
genes["lef"]="23663-26092"
genes["pagA"]="29382-31676"
genes["pagR"]="28160-28459"

# BLASTN parameters
OUTFMT="6 stitle sseqid qseqid qstart qend sstart send evalue bitscore score length pident nident mismatch positive means gapopen ppos qcovs qcovhsp"
NT="/data1/home/groups/readgp/nt/nt"

for s in ${!samples[@]}; do
    for p in ${plasmids[@]}; do
        wd="results/anthracis/${s}/${p}"
        mkdir -p ${wd}/coverage
        mkdir -p ${wd}/aligned-genes
        mkdir -p ${wd}/aligned-reads

        fq1="sra-pathogens/anthracis/${samples[$s]}/${s}_1.fastq.gz"
        fq2="sra-pathogens/anthracis/${samples[$s]}/${s}_2.fastq.gz"
        sam="${wd}/${s}.sam"
        bam="${wd}/${s}.bam"
        cov="${wd}/coverage/${s}.coverage"

        # Align using BWA, sort sam to bam and index bam
        bin/bwa mem -t 20 ${reference[$p]} ${fq1} ${fq2} > ${sam}
        samtools view -@ 10 -bS ${sam} | samtools sort -@ 10 - ${wd}/${s}
        samtools index -b ${bam}
        rm ${sam}

        # Use genomeCoverageBed to get the coverage for each position and plot
        # the coverage for differening sliding windows with 0.5 overlap
        genomeCoverageBed -d -ibam ${bam} > ${cov}
        scripts/mapping/plot-coverage.R ${cov} 0.5 ${p}

        if [ "$p" = "pXO1" ] ; then
            # Plot coverage across complete plasmid, and lethal genes
            scripts/mapping/plot-pxo1-anthrax-toxin-coverage.R ${cov} ${gff[$p]}

            # Extract reads mapped to each lethal gene and determine the top 5
            # blast hits against the NT database
            for g in ${!genes[@]}; do
                fasta=${wd}/aligned-genes/${g}.fasta
                blastn=${wd}/aligned-genes/${g}.blastn
                summary=${wd}/aligned-genes/${g}.summary

                # Extract mapped reads as FASTA
                samtools view ${bam} "gi|753442380|gb|CP009540.1|":${genes[$g]} | awk '{OFS="\t"; print ">"$1"\n"$10}' > ${fasta}

                # BLAST against nt
                blastn -query ${fasta} -outfmt "${OUTFMT}" -db ${NT} -num_threads 20  -max_target_seqs 5 > ${blastn}

                # Produce summary of hits
                grep -c "^>" ${fasta} | awk '{print "Total Reads:",$1,"\n"}' > ${summary}
                echo "Organism Hit Counts" >> ${summary}
                awk '{print $1,$2}' ${blastn} | sort | uniq -c | sort -rn >> ${summary}
            done
        fi

        # Extract aligned reads using bam2fastq and convert to fasta
        ofq="${wd}/aligned-reads/${samples[$s]}_${s}#.fastq"
        bin/bam2fastq -o ${ofq} --no-unaligned ${bam}
        cat ${wd}/aligned-reads/*.fastq | fastq_to_fasta -Q33 -n -o ${wd}/aligned-reads/${samples[$s]}_${s}.fasta
        gzip --best "${wd}/aligned-reads/${samples[$s]}_${s}_1.fastq"
        gzip --best "${wd}/aligned-reads/${samples[$s]}_${s}_2.fastq"
    done
done
