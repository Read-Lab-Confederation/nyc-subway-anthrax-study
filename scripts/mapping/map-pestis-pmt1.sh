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

plasmids=(pMT1)

declare -A samples
samples["SRR1748611"]="P00064"
samples["SRR1748618"]="P00070"
samples["SRR1748620"]="P00073"
samples["SRR1748663"]="P00102"
samples["SRR1748670"]="P00111"
samples["SRR1748743"]="P00164"
samples["SRR1748789"]="P00193"
samples["SRR1748804"]="P00208"
samples["SRR1748810"]="P00214"
samples["SRR1748826"]="P00230"
samples["SRR1748831"]="P00235"
samples["SRR1748835"]="P00239"
samples["SRR1748838"]="P00242"
samples["SRR1748839"]="P00243"
samples["SRR1748865"]="P00275"
samples["SRR1748869"]="P00279"
samples["SRR1748876"]="P00286"
samples["SRR1748884"]="P00294"
samples["SRR1748892"]="P00302"
samples["SRR1748893"]="P00303"
samples["SRR1748896"]="P00306"
samples["SRR1748911"]="P00321"
samples["SRR1748928"]="P00338"
samples["SRR1748936"]="P00347"
samples["SRR1748967"]="P00379"
samples["SRR1748970"]="P00382"
samples["SRR1748979"]="P00391"
samples["SRR1748980"]="P00392"
samples["SRR1748994"]="P00406"
samples["SRR1749003"]="P00416"
samples["SRR1749021"]="P00435"
samples["SRR1749061"]="P00475"
samples["SRR1749062"]="P00476"
samples["SRR1749080"]="P00494"
samples["SRR1749092"]="P00506"
samples["SRR1749188"]="P00604"
samples["SRR1749213"]="P00629"
samples["SRR1749237"]="P00653"
samples["SRR1749239"]="P00655"
samples["SRR1749244"]="P00660"
samples["SRR1749249"]="P00665"
samples["SRR1749291"]="P00707"
samples["SRR1749294"]="P00710"
samples["SRR1749302"]="P00718"
samples["SRR1749309"]="P00725"
samples["SRR1749350"]="P00767"
samples["SRR1749354"]="P00771"
samples["SRR1749362"]="P00779"
samples["SRR1749366"]="P00783"
samples["SRR1749386"]="P00803"
samples["SRR1749401"]="P00818"
samples["SRR1749496"]="P00916"
samples["SRR1749571"]="P00992"
samples["SRR1749572"]="P00993"
samples["SRR1749576"]="P00997"
samples["SRR1749581"]="P01002"
samples["SRR1749596"]="P01017"
samples["SRR1749599"]="P01020"
samples["SRR1749601"]="P01022"
samples["SRR1749615"]="P01036"
samples["SRR1749627"]="P01048"
samples["SRR1749631"]="P01052"
samples["SRR1749635"]="P01056"
samples["SRR1749642"]="P01063"
samples["SRR1749651"]="P01079"
samples["SRR1749665"]="P01093"
samples["SRR1749817"]="P01245"
samples["SRR1749835"]="P01270"
samples["SRR1749873"]="P01316"
samples["SRR1749894"]="P01337"
samples["SRR1749933"]="P01376"
samples["SRR1750038"]="P01606"

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

for s in ${!samples[@]}; do
    for p in ${plasmids[@]}; do
        wd="results/yersinia/${s}/${p}"
        # mkdir -p ${wd}/coverage
        # mkdir -p ${wd}/aligned-genes
        # mkdir -p ${wd}/aligned-reads

        fq1="sra-pathogens/yersinia/${samples[$s]}/${s}_1.fastq.gz"
        fq2="sra-pathogens/yersinia/${samples[$s]}/${s}_2.fastq.gz"
        sam="${wd}/${s}.sam"
        bam="${wd}/${s}.bam"
        cov="${wd}/coverage/${s}.coverage.gz"

        # Align using BWA, sort sam to bam and index bam
        # bin/bwa mem -t 20 ${reference[$p]} ${fq1} ${fq2} > ${sam}
        # samtools view -@ 10 -bS ${sam} | samtools sort -@ 10 - ${wd}/${s}
        # samtools index -b ${bam}
        # rm ${sam}

        # Use genomeCoverageBed to get the coverage for each position and plot
        # the coverage for differening sliding windows with 0.5 overlap
        # genomeCoverageBed -d -ibam ${bam} | gzip --best - > ${cov}
        # scripts/mapping/plot-coverage.R ${cov} 0.5 ${p}

        # Plot coverage across complete plasmid, and murine toxin gene
        # scripts/mapping/plot-pmt1-murine-toxin-coverage.R ${cov} ${gff[$p]}

        # R plot to PNG did not display correctly so convert the SVG to PNG
        # convert -density 300 -resize 2400x1200 ${wd}/coverage/${s}-murine-toxin.svg ${wd}/coverage/${s}-murine-toxin.png

        # Extract reads mapped to murine toxin gene and determine the top 5
        # blast hits against the NT database
        # for g in ${!genes[@]}; do
        #     fasta=${wd}/aligned-genes/${g}.fasta
        #     blastn=${wd}/aligned-genes/${g}.blastn
        #     summary=${wd}/aligned-genes/${g}.summary

            # Extract mapped reads as FASTA
        #     samtools view ${bam} "gi|5834685|emb|AL117211.1|":${genes[$g]} | awk '{OFS="\t"; print ">"$1"\n"$10}' > ${fasta}

            # BLAST against nt
        #     blastn -query ${fasta} -outfmt "${OUTFMT}" -db ${NT} -num_threads 20  -max_target_seqs 5 > ${blastn}

            # Produce summary of hits
        #     grep -c "^>" ${fasta} | awk '{print "Total Reads:",$1,"\n"}' > ${summary}
        #     echo "Organism Hit Counts" >> ${summary}
        #     awk '{print $1,$2}' ${blastn} | sort | uniq -c | sort -rn >> ${summary}
        # done
        awk '{print $1,$2}' ${wd}/aligned-genes/*.blastn | sort | uniq -c | sort -rn > ${wd}/aligned-genes/all.summary

        # Extract aligned reads using bam2fastq and convert to fasta
        # ofq="${wd}/aligned-reads/${samples[$s]}_${s}#.fastq"
        # bin/bam2fastq -o ${ofq} --no-unaligned ${bam}
        # cat ${wd}/aligned-reads/*.fastq | fastq_to_fasta -Q33 -n | gzip --best - > ${wd}/aligned-reads/${samples[$s]}_${s}.fasta.gz
        # gzip --best "${wd}/aligned-reads/${samples[$s]}_${s}_1.fastq"
        # gzip --best "${wd}/aligned-reads/${samples[$s]}_${s}_2.fastq"

    done
done
