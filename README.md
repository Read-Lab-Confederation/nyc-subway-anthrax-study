# New York City Subway-Metagenome: attempts to detect *B. anthracis* and *Y. pestis*. 

Here is how to edit Markdown:  
https://guides.github.com/features/mastering-markdown/

##Introduction
Analysis of NYC Subway Metagenome data  
http://www.sciencedirect.com/science/article/pii/S2405471215000022
See follow-up apologia and discussions:  
http://microbe.net/2015/02/17/the-long-road-from-data-to-wisdom-and-from-dna-to-pathogen/

The question is - why did the software give false positive results and what exactly was found in the subway?  

##DATA
State where and when we obtained data files

Which samples were said to contain anthrax

##CONTROLS (Robert)
Mix up real *B. anthracis* and *B. cereus* genomes

##Searching for pXO1 and pXO2 using BWA (Robert)
Compare results from metagenome samples with positive controls

[pXO1-P00134 BWA plots](https://www.dropbox.com/s/5hsev0fnbe3wsmc/P00134_SRR1748707.pdf)


##Using Kraken to search for the *B. anthracis* chromosome(Matthew)
  
I basically downloaded two complete genomes each of *Bacillus anthracis*, *Bacillus cereus* and *Bacillus thurigiensis*. For each complete genome sequence, I split up the FASTA file to a multi fasta file of 100bp sequences, and fed all 6 multi-fasta files to KRAKEN.

The reuslts are in seperate folders for each whole genome sequence:

[Anthrax WGS search](https://www.dropbox.com/sh/fwfi75ft4ny1qkk/AADF16diPK-cgV-CmRHzLjWTa?dl=0)   


##*B. anthracis* specific SNPs (Sandeep)


