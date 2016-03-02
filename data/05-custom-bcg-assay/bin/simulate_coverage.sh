#! /bin/bash
# Use ART to Simulate reads and generate 31-mer count using Jellyfish
#
TOP_DIR=$(pwd)
SAMPLE=$1
REFERENCE=$2
OUTPUT_DIR=$3
COVERAGE=$4
INTERATIONS=$5

mkdir -p ${OUTPUT_DIR}/logs
for i in `seq ${INTERATIONS}`; do
    mkdir -p ${OUTPUT_DIR}/
    RAND=`echo "$i * $COVERAGE * 100 / 1" | bc`
    ${TOP_DIR}/bin/art_illumina -l 100 -f ${COVERAGE} -na -ss HS20 -rs ${RAND} -i ${REFERENCE} -o ${OUTPUT_DIR}/${SAMPLE}-${i} 1> ${OUTPUT_DIR}/logs/${SAMPLE}-${i}.log 2>&1
    ${TOP_DIR}/bin/jellyfish count -C -m 31 -s 1M -o ${OUTPUT_DIR}/${SAMPLE}-${i}.jf ${OUTPUT_DIR}/${SAMPLE}-${i}.fq
    echo ${OUTPUT_DIR}/${SAMPLE}-${i}.jf >> ${OUTPUT_DIR}/${COVERAGE}x-jellyfish.txt
done
