---
title: Appendix
order: 6
---

## Project Data Structure

    Public
        data
            - Supplementary information from study and information about the study on SRA
        references
            - Reference genomes and plasmids used for mapping
        results
            - Mapping results of the B. anthracis and Yersinia samples against references
        scripts
            - All custom scripts used in these analyses

    In House
        bin
            - Symbolic links to programs built in src folder
        logs
            - Logged output from scripts
        sra-fastq
            - Each of the 1572 samples in FASTQ and FASTA format
        sra-pathogens
            - Symbolic links to FASTQ of samples containing B. anthracis and Yersinia
        sra-runs
            - Each of the 1572 samples in SRA format
        src
            - Programs downloaded and built for this project
