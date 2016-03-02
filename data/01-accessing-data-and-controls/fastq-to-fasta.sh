#! /bin/bash
# Convert FASTQ to FASTA
# usage: fastq_to_fasta [-h] [-r] [-n] [-v] [-z] [-i INFILE] [-o OUTFILE]
# Part of FASTX Toolkit 0.0.13.2 by A. Gordon (gordon@cshl.edu)

fasta="fasta"
while read INFILE
do
    OUTFILE="${INFILE/fastq.gz/fasta.gz}"
    echo "Converting ${INFILE} -> ${OUTFILE}"
    zcat ${INFILE} | fastq_to_fasta -Q33 -n -z -o ${OUTFILE}
done < <(find sra-fastq/ -name "*.fastq.gz")
