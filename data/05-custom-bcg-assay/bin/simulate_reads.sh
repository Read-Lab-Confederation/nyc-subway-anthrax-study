#! /bin/bash
# Use ART to Simulate reads and generate 31-mer count using Jellyfish
#
TOP_DIR=$(pwd)
SAMPLE=$1
REFERENCE=$2
OUTPUT_DIR=$3
INTERATIONS=$4
#COVERAGES=( 0.25 0.5 1 5 10 25 50 100 )
COVERAGES=( 50 )

mkdir -p ${OUTPUT_DIR}/logs
for i in "${COVERAGES[@]}"; do
    mkdir -p ${OUTPUT_DIR}/${i}x
    for j in `seq ${INTERATIONS}`; do
        RAND=`echo "$j * $i * 100 / 1" | bc`
        mkdir -p ${OUTPUT_DIR}/${i}x/${j}
        ${TOP_DIR}/bin/art_illumina -l 100 -f ${i} -na -ss HS20 -rs ${RAND} -i ${REFERENCE} -o ${OUTPUT_DIR}/${i}x/${j}/${SAMPLE} 1> ${OUTPUT_DIR}/logs/${SAMPLE}-${i}x-${j}.log 2>&1
        ${TOP_DIR}/bin/jellyfish count -C -m 31 -s 1M -o ${OUTPUT_DIR}/${i}x/${j}/${SAMPLE}.jf ${OUTPUT_DIR}/${i}x/${j}/${SAMPLE}.fq
        echo ${OUTPUT_DIR}/${i}x/${j}/${SAMPLE}.jf >> ${OUTPUT_DIR}/${i}x-jellyfish.txt

        # Don't need fastq, just the Jellyfish counts
        rm ${OUTPUT_DIR}/${i}x/${j}/${SAMPLE}.fq
    done
done
