#! /bin/bash
# Act as a wrapper for filter-kmer-by-sample.py
#
TOP_DIR=$(pwd)
STAT_DIR=$1
SAMPLES=$2
ORGANISM=$3
PREFIX=$4

rm -rf ${STAT_DIR}/${PREFIX}-kmer-stats.txt
for i in `find ${STAT_DIR} -name "*-kmer-by-sample.txt.gz"`; do
    BASENAME="${i/kmer-by-sample.txt.gz/${PREFIX}-kmer-stats.txt}"
    ${TOP_DIR}/bin/filter-kmer-by-sample.py ${i} ${SAMPLES} "${ORGANISM}" 1> ${BASENAME}
    cat ${BASENAME} >> ${STAT_DIR}/all-${PREFIX}-kmer-stats.txt
done

