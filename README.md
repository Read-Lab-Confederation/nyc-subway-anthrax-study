# Searching for anthrax (and plague) in the New York City Subway Metagenome*.

By Robert Petit III, Matthew Ezewudo, Sandeep Joseph and *Timothy D. Read,
Emory University School of Medicine, Atlanta, Georgia

##Introduction
In January 2015 Chris Mason and his team [published](http://www.sciencedirect.com/science/article/pii/S2405471215000022) 
an in-depth analysis of metagenomic data (environmental shotgun DNA sequence) from samples isolated from public surfaces 
in the New York City (NYC) Subway. Along with a lot of really interesting findings, the authors, rather unwisely, claimed 
to have detected DNA from the bacterial biothreat pathogens *Bacillus anthracis* (which causes anthrax) and 
*Yersinia pestis* (causes plague) in some of the samples. This predictably led to a press firestorm and skepticism from
scientistis on social media.  Chris and his team followed up with an 
[extensive re-analysis](http://microbe.net/2015/02/17/the-long-road-from-data-to-wisdom-and-from-dna-to-pathogen/) 
of the data on microbe.net, where they admitted that they overreached on the anthrax and plague claims and there was 
little to no evidence for the presence of those organisms

The question is - why did the software give false positive results and what exactly was found in the subway?  
##Methods

##Results

##Prerspective

Need - organism-sepcific detection.

Remember - "negative for anthrax" is in the context of the depth of the metagenome sequecning, which sets the sensitivity of  detection. 

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

## Acquiring NYC Subway System metagenome samples and identifying samples with evidence for *B. anthracis* and *Y. pestis* 
On February 10th 2015, we used *prefetch*, from [SRA Toolkit (v 2.4.4)](http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.4.4/), 
to download each of the 1572 runs associated with SRA study 
[SRP051511](http://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP051511). 
We then converted each run from SRA format to FASTQ format using *fastq-dump*, again from SRA Toolkit (v 2.4.4). 
Each FASTQ file was also converted to FASTA format using *fastq_to_fasta* from 
[FASTX Toolkit (v 0.0.13.2)](http://hannonlab.cshl.edu/fastx_toolkit/download.html). For each SRA run we used the 
information contained in [SRP051511_info.txt](/data/SRP051511_info.txt), acquired from SRP051511's RunInfo Table, to associate sample names used in the study with
their corresponding SRA run accessions.

We created the python script, 
[extract_pathogens.py](/scripts/extract-pathogens.py), 
to parse [DataTable5-metaphlan-metadata_v19.txt](/data/DataTable5-metaphlan-metadata_v19.txt). 
From this we were able to determine which NYC samples contained reads associated with *Bacillus antracis* and 
the *Yersinia* genus. We then used another python script, 
[map-pathogens.py](/scripts/map-pathogens.py) 
in order to associate the NYC sample names with their corresponding SRA run accession.

The following runs were identified as containing reads associated with *Bacillus antracis* and the *Yersinia* genus:

#####*Bacillus anthracis*
| NYC Sample Name | SRA Run Accession      |
|-----------------|------------------------|
| P00134          | SRR1748707, SRR1748708 |
| P00497          | SRR1749083             |

#####Yersinia genus
| NYC Sample Name | SRA Run Accession |
|-----------------|-------------------|
| P00064          | SRR1748611        |
| P00070          | SRR1748618        |
| P00073          | SRR1748620        |
| P00102          | SRR1748663        |
| P00111          | SRR1748670        |
| P00164          | SRR1748743        |
| P00193          | SRR1748789        |
| P00208          | SRR1748804        |
| P00214          | SRR1748810        |
| P00230          | SRR1748826        |
| P00235          | SRR1748831        |
| P00239          | SRR1748835        |
| P00242          | SRR1748838        |
| P00243          | SRR1748839        |
| P00275          | SRR1748865        |
| P00279          | SRR1748869        |
| P00286          | SRR1748876        |
| P00294          | SRR1748884        |
| P00302          | SRR1748892        |
| P00303          | SRR1748893        |
| P00306          | SRR1748896        |
| P00321          | SRR1748911        |
| P00338          | SRR1748928        |
| P00347          | SRR1748936        |
| P00379          | SRR1748967        |
| P00382          | SRR1748970        |
| P00391          | SRR1748979        |
| P00392          | SRR1748980        |
| P00406          | SRR1748994        |
| P00416          | SRR1749003        |
| P00435          | SRR1749021        |
| P00475          | SRR1749061        |
| P00476          | SRR1749062        |
| P00494          | SRR1749080        |
| P00506          | SRR1749092        |
| P00604          | SRR1749188        |
| P00629          | SRR1749213        |
| P00653          | SRR1749237        |
| P00655          | SRR1749239        |
| P00660          | SRR1749244        |
| P00665          | SRR1749249        |
| P00707          | SRR1749291        |
| P00710          | SRR1749294        |
| P00718          | SRR1749302        |
| P00725          | SRR1749309        |
| P00767          | SRR1749350        |
| P00771          | SRR1749354        |
| P00779          | SRR1749362        |
| P00783          | SRR1749366        |
| P00803          | SRR1749386        |
| P00818          | SRR1749401        |
| P00916          | SRR1749496        |
| P00992          | SRR1749571        |
| P00993          | SRR1749572        |
| P00997          | SRR1749576        |
| P01002          | SRR1749581        |
| P01017          | SRR1749596        |
| P01020          | SRR1749599        |
| P01022          | SRR1749601        |
| P01036          | SRR1749615        |
| P01048          | SRR1749627        |
| P01052          | SRR1749631        |
| P01056          | SRR1749635        |
| P01063          | SRR1749642        |
| P01079          | SRR1749651        |
| P01093          | SRR1749665        |
| P01245          | SRR1749817        |
| P01270          | SRR1749835        |
| P01316          | SRR1749873        |
| P01337          | SRR1749894        |
| P01376          | SRR1749933        |
| P01606          | SRR1750038        |

## Searching for *B. anthracis* pXO1 and pXO2 plasmids

[Analysis of pXO1 and pXO2](/_sections/results-anthracis.md)

##Using Kraken to search for the *B. anthracis* chromosome (Matthew)

We downloaded two complete genomes each of *Bacillus anthracis*, *Bacillus cereus* and *Bacillus thurigiensis*. For each complete genome sequence, I split up the FASTA file to a multi fasta file of 100bp sequences, and fed all 6 multi-fasta files to KRAKEN.

The reuslts are in seperate folders for each whole genome sequence:

[Anthrax WGS search](https://www.dropbox.com/sh/fwfi75ft4ny1qkk/AADF16diPK-cgV-CmRHzLjWTa?dl=0)   

 We also performed a similar approach on the actual control sequence reads for both *B. anthracis* and *B. cereus* to identify presence of *B.anthracis*.
 
 Next we ran kraken on a metagenome sample, adding varying proportions of *B.anthracis* to test for detection of the *anthracis* reads.
 
 The results from Kraken [Control Results](/results/Anthrax_control) suggests that kraken analysis could detect presence of minute amount of sequence reads from the organism.

 We used python script [parse_kraken.py](/scripts/parse_kraken.py) to extract the propotion of reads covered by *B.anthracis*, *B. cereus* group and other bacteria species respectively in the different samples.
 
###### *B.anthracis and B.cereus reads, species distribution*
 ![B. cereus](/results/Anthrax_control/cereus.png "B.cereus")
 ![B. anthracis](/results/Anthrax_control/anthracis.png "B. anthracis")
 
###### *metagenome and anthracis mixture reads, species distribution 5x,1x,0.5x and 0.25x coverage respectively*
 ![5x coverage](/results/Anthrax_control/5x_.png " anthracis 5x coverage")
 ![1x coverage](/results/Anthrax_control/1x_.png "anthracis 1x coverage")
 ![0.5x coverage](/results/Anthrax_control/0.5x_.png "anthracis 0.5x coverage")
 ![0.25x coverage](/results/Anthrax_control/0.25x_.png "anthracis 0.25x coverage")


##*B. anthracis* specific SNPs (Sandeep)

First, made a list of SNPs in B. anthracis versus B. cereus and extracted the 31-mers  <Sandeep describe how you did this>

Next, using the query:

>[txid86661\[Organism:exp\] NOT txid1392\[Organism:exp\] AND (biomol_genomic\[PROP\] AND refseq\[filter\])](http://www.ncbi.nlm.nih.gov/nuccore/?term=txid86661%5BOrganism%3Aexp%5D+NOT+txid1392%5BOrganism%3Aexp%5D+AND+(biomol_genomic%5BPROP%5D+AND+refseq%5Bfilter%5D)) 

We downloaded 6,948 *[B. cereus group](http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=86661)* 
(excluding *[Bacillus anthracis](http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=1392)*) 
sequences in FASTA format from NCBI Refseq to downselect against.
