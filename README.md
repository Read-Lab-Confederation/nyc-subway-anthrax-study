# New York City Subway-Metagenome: attempts to detect *B. anthracis* and *Y. pestis*. 

Here is how to edit Markdown:  
https://guides.github.com/features/mastering-markdown/

##Introduction
Analysis of NYC Subway Metagenome data  
http://www.sciencedirect.com/science/article/pii/S2405471215000022
See follow-up apologia and discussions:  
http://microbe.net/2015/02/17/the-long-road-from-data-to-wisdom-and-from-dna-to-pathogen/

The question is - why did the software give false positive results and what exactly was found in the subway?  

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

## Acquiring Data and Identifying Samples With Evidence For *B. anthracis* and *Y. pestis* 
On February 10th 2015, we used *prefetch*, from [SRA Toolkit (v 2.4.4)](http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.4.4/), 
to download each of the 1572 runs associated with SRA study 
[SRP051511](http://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP051511). 
We then converted each run from SRA format to FASTQ format using *fastq-dump*, again from SRA Toolkit (v 2.4.4). 
Each FASTQ file was also converted to FASTA format using *fastq_to_fasta* from 
[FASTX Toolkit (v 0.0.13.2)](http://hannonlab.cshl.edu/fastx_toolkit/download.html). For each SRA run we used the 
information contained in [SRP051511_info.txt](https://github.com/Read-Lab-Confederation/nyc-subway-metagenome/blob/master/data/SRP051511_info.txt), acquired from SRP051511's RunInfo Table, to associate sample names used in the study with
their corresponding SRA run accessions.

We created a python script 
[extract_pathogens.py](https://github.com/Read-Lab-Confederation/nyc-subway-metagenome/blob/master/scripts/extract-pathogens.py) 
which parses [DataTable5-metaphlan-metadata_v19.txt](https://github.com/Read-Lab-Confederation/nyc-subway-metagenome/blob/master/data/DataTable5-metaphlan-metadata_v19.txt) 
in order to determine which study samples contained reads associated with *Bacillus antracis* and the *Yersinia* genus. We then used another python script, [map-pathogens.py](https://github.com/Read-Lab-Confederation/nyc-subway-metagenome/blob/master/scripts/map-pathogens.py) 
in order to associate the study samples with their corresponding SRA run accession.

The following runs were identified as containing reads associated with *Bacillus antracis* and the *Yersinia* genus:

    Bacillus anthracis
      Study Sample    SRA Run Accession
      P00134          SRR1748707, SRR1748708
      P00497          SRR1749083
      
    Yersinia genus
      Study Sample    SRA Run Accession
      P00064          SRR1748611
      P00070          SRR1748618
      P00073          SRR1748620
      P00102          SRR1748663
      P00111          SRR1748670
      P00164          SRR1748743
      P00193          SRR1748789
      P00208          SRR1748804
      P00214          SRR1748810
      P00230          SRR1748826
      P00235          SRR1748831
      P00239          SRR1748835
      P00242          SRR1748838
      P00243          SRR1748839
      P00275          SRR1748865
      P00279          SRR1748869
      P00286          SRR1748876
      P00294          SRR1748884
      P00302          SRR1748892
      P00303          SRR1748893
      P00306          SRR1748896
      P00321          SRR1748911
      P00338          SRR1748928
      P00347          SRR1748936
      P00379          SRR1748967
      P00382          SRR1748970
      P00391          SRR1748979
      P00392          SRR1748980
      P00406          SRR1748994
      P00416          SRR1749003
      P00435          SRR1749021
      P00475          SRR1749061
      P00476          SRR1749062
      P00494          SRR1749080
      P00506          SRR1749092
      P00604          SRR1749188
      P00629          SRR1749213
      P00653          SRR1749237
      P00655          SRR1749239
      P00660          SRR1749244
      P00665          SRR1749249
      P00707          SRR1749291
      P00710          SRR1749294
      P00718          SRR1749302
      P00725          SRR1749309
      P00767          SRR1749350
      P00771          SRR1749354
      P00779          SRR1749362
      P00783          SRR1749366
      P00803          SRR1749386
      P00818          SRR1749401
      P00916          SRR1749496
      P00992          SRR1749571
      P00993          SRR1749572
      P00997          SRR1749576
      P01002          SRR1749581
      P01017          SRR1749596
      P01020          SRR1749599
      P01022          SRR1749601
      P01036          SRR1749615
      P01048          SRR1749627
      P01052          SRR1749631
      P01056          SRR1749635
      P01063          SRR1749642
      P01079          SRR1749651
      P01093          SRR1749665
      P01245          SRR1749817
      P01270          SRR1749835
      P01316          SRR1749873
      P01337          SRR1749894
      P01376          SRR1749933
      P01606          SRR1750038

Which samples were said to contain anthrax

##CONTROLS 
Recent B. anthracis SRA project with Illumna paired end reads.

http://trace.ncbi.nlm.nih.gov/Traces/sra/?run=DRR014739  

Recent B. cereus near-neighbor

http://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR642775

Recent PE Illumina Yersinia pestis

http://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1283952

Recent PE Illumina run of Y. pseudotuberculosis, close but less pathogenic relative

http://trace.ncbi.nlm.nih.gov/Traces/sra/?run=ERR245865


Mix up real *B. anthracis* and *Y. pestis* genomes with sample that does not contain either organism.  Also make dummy mixes of same environmental sample with *B. cereus* and *Y. pseudotuberculosis*.

## Searching for pXO1 and pXO2 using BWA
Samples P00134 and P00497 were identified as having reads from *B. anthracis*. We tested if these 
reads mapped against the PXO1 and PXO2 plasmids of *B. anthracis*. 
    [
    ADD SUMMARY OF PXO1 and PXO2
    PXO1 is responsible for encoding the components of the anthrax toxin. The genes required 
    to produce the anthrax toxin are *pagA*, *pagR*, *lef* and *cya*.
    ]

We used *BWA (v 0.7.5a-r405)* to map reads from samples P00134 and P00497 against the reference plasmid. For 
PXO1 we used reference [CP009540](http://www.ncbi.nlm.nih.gov/nuccore/CP009540.1), and for PXO2 we used reference [NC_007323](http://www.ncbi.nlm.nih.gov/nuccore/50163691). 
We also mapped each samples against PXO1-like plasmid [NC_005707](http://www.ncbi.nlm.nih.gov/nuccore/NC_005707). This process was automated using the script *[map-anthracis-plasmids.sh](https://github.com/Read-Lab-Confederation/nyc-subway-metagenome/blob/master/scripts/mapping/map-anthracis-plasmids.sh).

### PXO1 Results
For each of the plots below the top plot is against the whole PXO1 plasmid. The coverage is based on 1,000bp sliding windows with an overlap of 500bp. Each subplot against the genes *pagA*, *pagR*, *lef* and *cya* is the actual coverage in that region.

#### Sample P00134 (Run: SRR1748707)
![P00134 (Run: SRR1748707)](https://github.com/Read-Lab-Confederation/nyc-subway-metagenome/blob/master/results/pXO1/coverage/P00134_SRR1748707-antrax-toxin.png "P00134 (Run: SRR1748707)")

#### Sample P00134 (Run: SRR1748708)
![P00134 (Run: SRR1748708)](https://github.com/Read-Lab-Confederation/nyc-subway-metagenome/blob/master/results/pXO1/coverage/P00134_SRR1748708-antrax-toxin.png "P00134 (Run: SRR1748708)")

#### Sample P00497 (Run: SRR1749083)
![P00497 (Run: SRR1749083)](https://github.com/Read-Lab-Confederation/nyc-subway-metagenome/blob/master/results/pXO1/coverage/P00497_SRR1749083-antrax-toxin.png "P00497 (Run: SRR1749083)")





Compare results from metagenome samples with positive controls

[pXO1-P00134 BWA plots](https://www.dropbox.com/s/5hsev0fnbe3wsmc/P00134_SRR1748707.pdf)


##Using Kraken to search for the *B. anthracis* chromosome(Matthew)
  
I basically downloaded two complete genomes each of *Bacillus anthracis*, *Bacillus cereus* and *Bacillus thurigiensis*. For each complete genome sequence, I split up the FASTA file to a multi fasta file of 100bp sequences, and fed all 6 multi-fasta files to KRAKEN.

The reuslts are in seperate folders for each whole genome sequence:

[Anthrax WGS search](https://www.dropbox.com/sh/fwfi75ft4ny1qkk/AADF16diPK-cgV-CmRHzLjWTa?dl=0)   


##*B. anthracis* specific SNPs (Sandeep)


