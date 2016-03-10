#! /bin/bash
# Act as a wrapper for query-jellyfish-counts.py
#
TOP_DIR=$(pwd)
KMERS=$1
SIMULATIONS=$2
OUTPUT_DIR=$3
JELLYFISH=$4

# Assume SIMULATION_DIR is as follows SIMULATION_DIR/${coverage|/${iteration}/jellyfish-counts.jf
for c in `ls ${SIMULATIONS}`; do
    if [ $c != "logs" ]; then
        for i in `find ${SIMULATIONS}/$c -mindepth 1 -maxdepth 1 ! -path "*.txt"`; do
            mkdir -p ${OUTPUT_DIR}/${c}
            ITERATION=`basename ${i}`
            STATS=${OUTPUT_DIR}/${c}/${ITERATION}-kmer-stats.txt
            BY_SAMPLE=${OUTPUT_DIR}/${c}/${ITERATION}-kmer-by-sample.txt
            COUNTS=${OUTPUT_DIR}/${c}/${ITERATION}-kmer-counts.txt

            ${TOP_DIR}/bin/query-jellyfish-counts.py ${KMERS} ${i} ${JELLYFISH} ${STATS} ${BY_SAMPLE} ${COUNTS} 1>> ${OUTPUT_DIR}/${c}.log 2>&1
            gzip --best ${BY_SAMPLE}
        done
    fi
done
