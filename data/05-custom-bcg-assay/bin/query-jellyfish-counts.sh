#! /bin/bash
# Act as a wrapper for query-jellyfish-counts.py
#
TOP_DIR=$(pwd)
KMERS=$1
SIMULATIONS=$2
OUTPUT_DIR=$3
JELLYFISH=$4

# Assume SIMULATION_DIR is as follows SIMULATION_DIR/${coverage}/${iteration}/jellyfish-counts.jf
mkdir -p ${OUTPUT_DIR}
for i in `find ${SIMULATIONS} -mindepth 1 -maxdepth 1 ! -path "*.txt"`; do
    ITERATION=`basename ${i}`
    STATS=${OUTPUT_DIR}/${ITERATION}-kmer-stats.txt
    BY_SAMPLE=${OUTPUT_DIR}/${ITERATION}-kmer-by-sample.txt
    COUNTS=${OUTPUT_DIR}/${ITERATION}-kmer-counts.txt

    ${TOP_DIR}/bin/query-jellyfish-counts.py ${KMERS} ${i} ${JELLYFISH} ${STATS} ${BY_SAMPLE} ${COUNTS} 1>> ${OUTPUT_DIR}/stats.log 2>&1
    gzip ${BY_SAMPLE}
done
