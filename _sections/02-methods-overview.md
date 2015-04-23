---
title: "Methods Overview"
order: 2
---
The results are organized in 4 sections:

1.  **Accessing metagenome data and controls:**  Where we obtained the data and how we constructed artificial controls by mixing recent whole genome shotgun data from pathogens and near-neighbors with NYC subway metagenome data.
2.  **Mapping metagenome data to plasmids (and chromosomes):**    Looking at the patterns of sequence coverage over the key virulence associated plasmids, pXO1 , pXO2 (and pMT of *Y. pestis*) in metagenome samples and controls.
3.  **Species identification with Kraken :** [Kraken](http://genomebiology.com/2014/15/3/R46) is a popular k-mer based software for read identification.  We showed that Kraken was a sensitive way to find *B. anthracis* when it was present in low abundance but the method also produced a small number of false positive reads on near-neighbor *B. cereus* sequence control data.
4.  **Custom [SNP (single nucleotide polymorphism)](http://en.wikipedia.org/wiki/Single-nucleotide_polymorphism) assays for *B. anthracis* core genome:** We identified 31-mer words that corresponded to SNPs in the core genome of *B. anthracis* that were not found in close relatives. This gave a rapid specific test for *B. anthracis*.  However, we still detected two potential positive SNPs in one of the NYC subway samples.
