#! /bin/bash
# Act as a wrapper for query-jellyfish-counts.py
#
TOP_DIR=$(pwd)
FASTA=$1
OUTPUT_DIR=$3
TAG=$4
TAG=${TAG##*/}
TAG=${TAG%-*}


JELLYFISH_DIR=$2/${TAG}
STATS=${OUTPUT_DIR}/$5-${TAG}-kmer-stats.txt
BY_SAMPLE=${OUTPUT_DIR}/$5-${TAG}-kmer-by-sample.txt
COUNTS=${OUTPUT_DIR}/$5-${TAG}-kmer-counts.txt

${TOP_DIR}/bin/query-jellyfish-counts.py ${FASTA} ${JELLYFISH_DIR} ${STATS} ${BY_SAMPLE} ${COUNTS} 1> ${OUTPUT_DIR}/${TAG}.log 2>&1
