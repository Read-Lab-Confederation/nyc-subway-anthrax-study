# Searching for anthrax (and plague) in the New York City Subway Metagenome.

By Robert Petit III, Matthew Ezewudo, Sandeep J. Joseph and *Timothy D. Read,
Emory University School of Medicine, Atlanta, Georgia

##Introduction
In January 2015 Chris Mason and his team [published](http://www.sciencedirect.com/science/article/pii/S2405471215000022) 
an in-depth analysis of metagenomic data (environmental shotgun DNA sequence) from samples isolated from public surfaces 
in the New York City (NYC) Subway. Along with a ton of really interesting findings, the authors claimed 
to have detected DNA from the bacterial biothreat pathogens *Bacillus anthracis* (which causes anthrax) and 
*Yersinia pestis* (causes plague) in some of the samples. This predictably led to a press firestorm and skepticism from
scientistis on social media.  Chris and his team followed up with an 
[re-analysis](http://microbe.net/2015/02/17/the-long-road-from-data-to-wisdom-and-from-dna-to-pathogen/) 
of the data on microbe.net, where they admitted that they overreached on the anthrax and plague claims and there was 
actually little to no evidence for the presence of those organisms.

We were interested in looking deeping into why the software gave false positive results and what exactly was found in the subway samples.  We stared messing around with some programs and decided to wrap up the results on this site.  The study rasied very timely questions about that hot topic of using metagenomics for detection. Mostly, we looked at *B. anthracis* but we present some results looking at *Y.pstis* and can update this later.

##Overview
The results are organized in 4 sections:  

1.  **Accessing metagenome data and controls:**  Where we obtained the data and how we constructed artificial controls by mixing recent whole genome shotgun data from pathogens and near-neighbors with NYC subway metagenome data.  
2.  [**Mapping plasmids to metagenome data:**](/_sections/results-anthracis.md)    Looked at the patterns of sequence coverage over the key virulence associated plasmids, pXO1 , pXO2 (and pMT of *Y. pstis*) in metagenome samples and controls.  
3.  **Kraken metagenome detection:** Kraken is a popular kmer based software for read identification.  We rhowed that Kraken was sensitive for *B. anthracis* detection but also produced a small n umber of false positive reads in the B. cereus genome project.  
4.  **Custom SNP assays for B. anthracis:** We identified 31-mer words that corresponded to SNPs in the core genome of *B. anthracis* that were not found in clase relative. This gave a rapid specific test for B. anthracis.  However, w still detected two potential positive SNPs in one of the NYC subway samples.

##Summary of conclusions

<finish this section when the work is completed>

##Perspective

We believe that there is no one-size-fits-all approach to bacterial species indetification in metagenome samples for several reasons.  Perhaps most importantly, there is no consistent definition for a bacterial species that can be used a cutoff for sequence dentity.  Some species (like *B. cereus* and *B. anthracis* can be 99% similar to each other).  Secondly, some species distinctions rely on the presence or absence of mobile elements (again *B. cereus* and *B. anthracis* are a great example), which are hard to model using a uniform approach like Kraken or Metaphlan. These plasmids or phage or often modular in structure with very similar 'backbones' but with key genes replaced.  Finally, the generalist species indentification programs rely on databases of sequenced genomes, whcih are not uniform in their coverage of different species, or different lineages within species.

If you were to have an organism-sepcific detection algorthim for *B. cereus* you would need to accout for the presence of the plasmids about xx times the coverage as the chromsome.  You would also expext even coverage across all the 

Even if you have developed a sophisticated algorithm you will still need to use judgement in interpreting results.  There is still much that we dont know about bacteria in the envirnemnt, especialy the conditions under which they exxchange DNA.  

Finally, the negative result is subject to a lot of nuances.  The limit of detectection will be affected by the amount of sequence generated and the "negative for anthrax" is in the context of the depth of the metagenome sequecning, which sets the sensitivity of  detection. 

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
 ![B. cereus](/results/kracken/anthracis-control/cereus.png "B.cereus")
 ![B. anthracis](/results/kracken/anthracis-control/anthracis.png "B. anthracis")
 
###### *metagenome and anthracis mixture reads, species distribution 5x,1x,0.5x and 0.25x coverage respectively*
 ![5x coverage](/results/kracken/anthracis-control/5x_.png " anthracis 5x coverage")
 ![1x coverage](/results/kracken/anthracis-control/1x_.png "anthracis 1x coverage")
 ![0.5x coverage](/results/kracken/anthracis-control/0.5x_.png "anthracis 0.5x coverage")
 ![0.25x coverage](/results/kracken/anthracis-control/0.25x_.png "anthracis 0.25x coverage")


### *B. anthracis* specific SNPs (Sandeep)

First, made a list of SNPs in B. anthracis versus B. cereus and extracted the 31-mers  <Sandeep describe how you did this>

#### Identifying B. anthracis specific SNPs from the whole genome alignment*

ProgressiveMAUVE was used to perform whole genome alignment of 10 Bacillus species genomes that were selected based on previously published whole genome phylogeny in order to capture the maximum diversity of all the Bacillus species. ProgressiveMAUVE was performed using script [progressiveMAUVE.sh](/scripts/progressiveMAUVE.sh). The B. cereus strain ATCC 10987 was used as the reference genome for SNP calling (see below) The 10 Bacillus species genomes and the progressiveMAUVE output can be accessed [here](/MAUVE).
    
#### Table 1: Genome sequences of Bacillus species used for ProgressiveMAUVE. 
| Strain name   | Species                | File name    | RefSeq/SRA  | Clade |
|---------------|------------------------|--------------|-------------|-------|
| ATCC 10987     | Bacillus cereus        | BCE.fasta    | NC_003909   | 1     |
| m1293         | Bacillus cereus        | bce1_.fasta  | SRX096996   | 1     |
| Ames ancestor | Bacillus anthracis     | GBAA_.fasta  | NC_007530   | 1     |
| Rock3-44      | Bacillus cereus        | bce22_.fasta | SRX098720   | 3     |
| AH603         | Bacillus cereus        | bce26_.fasta | SRX098634   | 3     |
| ATCC 10876    | Bacillus cereus        | bce2_.fasta  | SRX096081   | 2     |
| MM3           | Bacillus cereus        | bce6_.fasta  | SRX098606   | 1     |
| ATCC 1479     | Bacillus cereus        | BC.fasta     | NC_004721-2 | 2     |
| E33L          | Bacillus cereus        | BCZK.fasta   | NC_006274   | 1     |
| BGSC 4BD1     | Bacillus thuringiensis | bth11_.fasta | SRX098635   | 2     |

The MAUVE alignment was loaded into MAUVE Version 2.4.0 inorder to visualize the alignment and to export the SNPs using the "export SNP" option. The nucleotide positions of the reference genome (BCE.fasta, ATCC 10987) as well as the corresponding SNP position on the B. anthracis genome, along with the 10 nucleotide pattern at the position were extracted. In this nucleotide pattern, the first nuleotide will be from the reference genome and the 3rd nucleotide will be from the B. anthracis genome. Anthrax-specific SNPs (i.e SNP nucloetide positions found in only the Anthrax genome) were identified using the script [SNPPattern.awk](/scripts/SNPPattern.awk), where a number (corresponding to the position of the nucleotide variant) were assigned. Finally, those nucleotide positions that have only the number 3 assigned (position of the B. anthracis genome on the alignment) was extracted along with the corresponding nucleotide positions at the reference genome as well as at the B. anthracis genome (*[Anthrax_specific-SNPPattern.txt](/data/Anthrax_specific-SNPPattern.txt)*). Total 9538 Anthrax-specific SNPs were identified.

#### Generating all the B. cereus 31-mers to downselect against the 9538 Anthrax-specific SNPs*

Next, using the query:>[txid86661\[Organism:exp\] NOT txid1392\[Organism:exp\] AND (biomol_genomic\[PROP\] AND refseq\[filter\])](http://www.ncbi.nlm.nih.gov/nuccore/?term=txid86661%5BOrganism%3Aexp%5D+NOT+txid1392%5BOrganism%3Aexp%5D+AND+(biomol_genomic%5BPROP%5D+AND+refseq%5Bfilter%5D)), we downloaded 6,948 *[B. cereus group](http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=86661)* 
(excluding *[Bacillus anthracis](http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=1392)*) 
sequences in FASTA format from NCBI Refseq to downselect against the 9538 Anthrax-specific SNPs. 

First, we extracted all the 31-mers in and around all the the 9538 Anthrax-specific SNP positions from the B. anthracis Ames ancestor genome so that the 16th nucleotide in the 31-mer will be the SNP at that specific position. This was done using the script [extractKmers.sh](/scripts/extractKmers.sh).

We then generated all the 31-mers present in all the 6,948  *[B. cereus group](http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=86661)* 
(excluding *[Bacillus anthracis](http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=1392)*) using *[JELLYFISH](http://www.cbcb.umd.edu/software/jellyfish/)* by the script [generate-31-mers-BCereusGroup.sh](/scripts/generate-31-mers-BCereusGroup.sh). The resultant JELLYFISH database of all the B. cereus 31-mers were queried by the 9538 Anthrax-specific 31-mers using the script QueryBcereusGroup.sh(/scripts/QueryBcereusGroup.sh) and all those 31-mers with zero-counts against the Bcereus database were outputed. A total of 1793/9538 31-mers that had zero- counts were obtained and these 31-mers with the Anthrax-specific SNP at the 16th position was considered as the Anthrax-specific 31-mers and was used in all the further down stream analysis.

#### Generating all 31-mers for the B. anthracis and B. cereus control groups, and NYC Subway system B. anthracis positive samples, and querying all the Anthrax-specific 1793 31-mers against them.

Here we used the same control samples previously described for B. anthracis and B. cereus in this analysis. All the 31-mers were generated using the script [Getting_KmerCounts.sh](/scripts/Getting_KmerCounts.sh).The B. anthracis-specific 1793 31-mers were queried against the 31-mer database of all the control and real samples using [QueryAnthraxSNPs.sh](/scripts/QueryAnthraxSNPs.sh) and the 31-mers that matched the B. anthracis-specific 1793 31-mers were extracted along with their counts in each sample.

#### Table 2: Summary of the number of matched 31-mers from each of B. anthracis control samples with the 1793 Anthrax-specific 31-mers. 
|        Study       | Number of matched 31-mers with Anthrax-specific 31-mers |
|:------------------:|:--------------------------------:|
|   SRR1749070-0x_1  |                 0                |
|   SRR1749070-0x_2  |                 0                |
| SRR1749070-0.01x_1 |                 3                |
| SRR1749070-0.01x_2 |                 1                |
| SRR1749070-0.05x_1 |                20                |
| SRR1749070-0.05x_2 |                22                |
|  SRR1749070-0.1x_1 |                33                |
|  SRR1749070-0.1x_2 |                39                |
| SRR1749070-0.25x_1 |                106               |
| SRR1749070-0.25x_2 |                96                |
|  SRR1749070-0.5x_1 |                198               |
|  SRR1749070-0.5x_2 |                195               |
|  SRR1749070-1.0x_1 |                352               |
|  SRR1749070-1.0x_2 |                377               |
|  SRR1749070-5.0x_1 |               1153               |
|  SRR1749070-5.0x_2 |               1158               |

For those control samples that didn't had any Anthrax reads (0X coverage), none of the Anthrax-specific 31-mers were matched, a result that is expected. Similarily, the number of 31-mers matched to the Anthrax 31-mers were directly proportional to the coverage or the amount of B. anthracis reads present in those samples. From Table 2 we may infer that to detect the presence of Anthrax, there should be atleast 0.01X coverage of B. anthracis reads in the sample, which can be considered as the minimum threshold for Anthrax detection using 31-mers. Even at 5X coverage only 65% of the Anthrax-specific 31-mers were matched.

#### Table 3: Summary of the number of matched 31-mers from each of B. cereus control samples with the 1793 Anthrax-specific 31-mers.

|        Study       | Number of matched 31-mers with Anthrax-specific 31-mers  |
|:------------------:|:--------------------------------:|
|   SRR1749070-0x_1  |                 0                |
|   SRR1749070-0x_2  |                 0                |
| SRR1749070-0.01x_1 |                 0                |
| SRR1749070-0.01x_2 |                 0                |
| SRR1749070-0.05x_1 |                 0                |
| SRR1749070-0.05x_2 |                 0                |
|  SRR1749070-0.1x_1 |                 0                |
|  SRR1749070-0.1x_2 |                 0                |
| SRR1749070-0.25x_1 |                 0                |
| SRR1749070-0.25x_2 |                 0                |
|  SRR1749070-0.5x_1 |                 0                |
|  SRR1749070-0.5x_2 |                 0                |
|  SRR1749070-1.0x_1 |                 0                |
|  SRR1749070-1.0x_2 |                 0                |
|  SRR1749070-5.0x_1 |                 0                |
|  SRR1749070-5.0x_2 |                 0                |

As expected, none of the 31-mers from the samples in the B. cereus control group matched to the Anthrax-specific 31-mers. This indicates how reliable/sensitive are the 1793 Anthrax-specific 31-mers generated in this analysis inorder to differentiate the highly similar B. anthracis and B. cereus genomes at the nucleotide level.

#### Table 4: Summary of the number of matched 31-mers from each of the 3 NY Subway sytem Anthrax positve samples with the 1793 Anthrax-specific 31-mers.
|    Study   |Number of matched 31-mers with Anthrax-specific 31-mers |
|:----------:|:--------------------------------:|
| SRR1749083 |                 0                |
| SRR1748708 |                 2                |
| SRR1748707 |                 0                |

Only sample SRR1748708 has two 31-mers that matched Anthrax-specific 31-mers with a frequency count of 2 for each of the 31-mers. This might indicate, if some how B. anthracis is present in the sample, it would be only at amount equvalent to 0.01X coverage. 31-mers from the rest of the "Anthrax positive" samples didn't match to any of the Anthrax-specific 31-mers, indicating the absence of B. anthracis.


