# Searching for pXO1 and pXO2
  
    ADD SUMMARY OF PXO1 and PXO2
    PXO1 is responsible for encoding the components of the anthrax toxin. The genes required 
    to produce the anthrax toxin are *pagA*, *pagR*, *lef* and *cya*. 



## Identifying NYC samples with evidence of *B. anthracis*
The supplementary MetaPhlAn output (*[DataTable5-metaphlan-metadata_v19.txt](/data/DataTable5-metaphlan-metadata_v19.txt)*) 
was parsed using *[extract_pathogens.py](/scripts/extract-pathogens.py)*. From this we were able to determine which 
NYC samples contained reads associated with *Bacillus antracis*. These samples were then associated with their 
corresponding SRA run accession using *[map-pathogens.py](/scripts/map-pathogens.py)* (Table 1). For the remainder 
of these analysis, each of these samples will be referred to by their SRA run accession.

#### Table 1: NYC Subway System samples with evidence of *B. anthracis*.
| NYC Study Sample | SRA Run Accession      |
|------------------|------------------------|
| P00134           | SRR1748707, SRR1748708 |
| P00497           | SRR1749083             |

## Creating metagenomic controls for *B. anthracis* and a close relative, *B. cereus*
As a control for these analysis we randomly selected NYC sample SRR1749070 which is *B anthracis* free and added 
different amounts of known *B. anthracis* or *B. cereus* sequences to it. Using these controls we can produce
results we would expect to see if *B. anthracis* is present in the metagenomic sequences. We used *B. anthracis*
([DRR014739](http://trace.ncbi.nlm.nih.gov/Traces/sra/?run=DRR014739)) and *B. cereus* 
([SRR642775](http://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR642775)) sequencing projects as our controls. The reads from the *B. anthracis* project were trimmed from 300bp (Illumina MiSeq) down to the first 100bp using *fastq_trimmer (v0.0.13.2)*
from FASTX Toolkit. This was necessary to make the reads more similar to the NYC sample, SRR1749070, which is 100bp 
Illumina HiSeq 2000 reads. This step was not necessary for the *B. cereus* control because it too is 100bp Illumina 
HiSeq 2000 reads. 

Each control was then randomly subsampled using *seqtk (commit 43ff625a3211b51f301cb356a34fb8d1e593d50a)* to 
0.25x, 0.5x, 1x, and 5x coverage. Coverage was estimated using the total size of the *B. anthracis* genome, the 
pXO1 plasmid and the pXO2 plasmids was *B. anthracis*. For *B. cereus* we used only the total size of the *B. cereus* 
genome. These subsampled coverages were then added to the sequences of the *B. anthracis* free NYC sample SRR1749070. 
This produced5 samples for each control, one with no control sequences and 4 with the different levels of coverage. 
These samples were then subjected to the same mapping pipeline as the *B. anthracis* positive samples menthion above.
This process of downloading and preparing the controls for analysis was automated using the following scripts 
[get-anthracis-control.sh](/scripts/get-anthracis-control.sh) (*B. anthracis*) and
[get-cereus-control.sh](/scripts/get-cereus-control.sh) (*B. cereus*).


## Mapping against pXO1 and pXO2 using BWA
We used *BWA (v 0.7.5a-r405)* to map reads from each of the NYC and control samples against the reference 
pXO1 and pXO2 plasmids. For pXO1 we used reference [CP009540](http://www.ncbi.nlm.nih.gov/nuccore/CP009540.1), 
and for pXO2 we used reference [NC_007323](http://www.ncbi.nlm.nih.gov/nuccore/50163691). The SAM output was then
converted to sorted BAM and indexed suing *samtools (v 1.1)*. The per base coverage was extracted using 
*genomeCoverageBed* from *bedtools (v2.16.2)*. Coverage across the complete plasmid was then plotted for mulitple
sliding windows using the Rscript *[plot-coverage.R](/scripts/mapping/plot-coverage.R)*. Reads that mapped to 
the plasmids were extracted and saved in both FASTQ and FASTA format using *bam2fastq (v1.1.0)* and *fastq_to_fasta* 
from *FASTX Toolkit v 0.0.13.2*.

For pXO1, coverage of the genes (*cya*, *lef*, *pagA* and *pagR*) related to the anthrax toxin is most important. To
visualize this, an alternate plot was created which included subplots of coverage of these genes using the RScript
*[plot-pxo1-anthrax-toxin-coverage.R](/scripts/mapping/plot-pxo1-anthrax-toxin-coverage.R)*. The reads that mapped 
to each gene were also extracted and blasted (*blastn v2.2.30+)* against the NT database (built on Feb 9, 2015). For 
each gene, a count of the organism names of which the top five hits of each read belonged to was recorded. 

This analysis was automated using the following scripts 
*[map-anthracis-plasmids.sh](/scripts/mapping/map-anthracis-plasmids.sh)*
(NYC samples) and *[map-anthracis-controls.sh](/scripts/mapping/map-anthracis-controls.sh)* (control samples). 
Summaries of the results of mapping samples against pXO1 and pXO2 were created using the following scripts:
*[mapping-coverage-summary.py](/scripts/mapping/mapping-coverage-summary.py)*, 
*[mapping-gene-summary.py](/scripts/mapping/mapping-gene-summary.py)* and 
*[mapping-top-blast-hits-summary.py](/scripts/mapping/mapping-top-blast-hits-summary.py)*

## pXO1 Results
For each sample we calculated the summary statistics of coverage across the complete pXO1 plasmid (Table 2).
NYC sample SRR174708 has the best average coverage (2x), but each of the NYC samples have a median coverage of 
0x. These results are nearly the opposite of the *B. anthracis* control. Even at an estimated coverage of 0.25x, 
the mean coverage of pXO1 is 1.23x (median 1x). The mean coverage only improves as the estimated coverage of 
*B. anthracis* reads increases. The *B. cereus* control includes low levels of pXO1 overage.

#### Table 2: Summary statistics of per base coverage against pXO1 for NYC samples and controls.
###### *NYC Subway System*
| Study      | Min. | 1st Qu. | Median | Mean   | 3rd Qu. | Max | Coverage Plots |
|------------|------|---------|--------|--------|---------|-----|----------------|
| SRR1749083 | 0    | 0.0000  | 0.0000 | 0.4652 | 0.0000  | 44  |[pXO1](/results/anthracis/SRR1749083/pXO1/coverage/SRR1749083-coverage.pdf), [Toxin](/results/anthracis/SRR1749083/pXO1/coverage/SRR1749083-anthrax-toxin.pdf)|
| SRR1748708 | 0    | 0.0000  | 0.0000 | 2.0766 | 4.0000  | 100 |[pXO1](/results/anthracis/SRR1748708/pXO1/coverage/SRR1748708-coverage.pdf), [Toxin](/results/anthracis/SRR1748708/pXO1/coverage/SRR1748708-anthrax-toxin.pdf)|
| SRR1748707 | 0    | 0.0000  | 0.0000 | 0.2506 | 0.0000  | 44  |[pXO1](/results/anthracis/SRR1748707/pXO1/coverage/SRR1748707-coverage.pdf), [Toxin](/results/anthracis/SRR1748707/pXO1/coverage/SRR1748707-anthrax-toxin.pdf)|
###### *B. anthracis* control
| Study            | Min. | 1st Qu. | Median  | Mean    | 3rd Qu. | Max | Coverage Plots |
|------------------|------|---------|---------|---------|---------|-----|----------------|
| SRR1749070-0x    | 0    | 0.0000  | 0.0000  | 0.0103  | 0.0000  | 2   |[pXO1](/results/anthracis-control/SRR1749070-0x/pXO1/coverage/SRR1749070-0x-coverage.pdf), [Toxin](/results/anthracis-control/SRR1749070-0x/pXO1/coverage/SRR1749070-0x-anthrax-toxin.pdf)|
| SRR1749070-0.25x | 0    | 0.0000  | 1.0000  | 1.2305  | 2.0000  | 13  |[pXO1](/results/anthracis-control/SRR1749070-0.25x/pXO1/coverage/SRR1749070-0.25x-coverage.pdf), [Toxin](/results/anthracis-control/SRR1749070-0.25x/pXO1/coverage/SRR1749070-0.25x-anthrax-toxin.pdf)|
| SRR1749070-0.5x  | 0    | 1.0000  | 2.0000  | 2.4250  | 3.0000  | 20  |[pXO1](/results/anthracis-control/SRR1749070-0.5x/pXO1/coverage/SRR1749070-0.5x-coverage.pdf), [Toxin](/results/anthracis-control/SRR1749070-0.5x/pXO1/coverage/SRR1749070-0.5x-anthrax-toxin.pdf)|
| SRR1749070-1x    | 0    | 3.0000  | 4.0000  | 4.8863  | 6.0000  | 36  |[pXO1](/results/anthracis-control/SRR1749070-1x/pXO1/coverage/SRR1749070-1x-coverage.pdf), [Toxin](/results/anthracis-control/SRR1749070-1x/pXO1/coverage/SRR1749070-1x-anthrax-toxin.pdf)|
| SRR1749070-5x    | 0    | 19.0000 | 23.0000 | 24.3656 | 28.0000 | 175 |[pXO1](/results/anthracis-control/SRR1749070-5x/pXO1/coverage/SRR1749070-5x-coverage.pdf), [Toxin](/results/anthracis-control/SRR1749070-5x/pXO1/coverage/SRR1749070-5x-anthrax-toxin.pdf)|
###### *B. cereus* control
| Study            | Min. | 1st Qu. | Median | Mean   | 3rd Qu. | Max | Coverage Plots |
|------------------|------|---------|--------|--------|---------|-----|----------------|
| SRR1749070-0x    | 0    | 0.0000  | 0.0000 | 0.0103 | 0.0000  | 2   |[pXO1](/results/cereus-control/SRR1749070-0x/pXO1/coverage/SRR1749070-0x-coverage.pdf), [Toxin](/results/cereus-control/SRR1749070-0x/pXO1/coverage/SRR1749070-0x-anthrax-toxin.pdf)|
| SRR1749070-0.25x | 0    | 0.0000  | 0.0000 | 0.0368 | 0.0000  | 7   |[pXO1](/results/cereus-control/SRR1749070-0.25x/pXO1/coverage/SRR1749070-0.25x-coverage.pdf), [Toxin](/results/cereus-control/SRR1749070-0.25x/pXO1/coverage/SRR1749070-0.25x-anthrax-toxin.pdf)|
| SRR1749070-0.5x  | 0    | 0.0000  | 0.0000 | 0.0522 | 0.0000  | 9   |[pXO1](/results/cereus-control/SRR1749070-0.5x/pXO1/coverage/SRR1749070-0.5x-coverage.pdf), [Toxin](/results/cereus-control/SRR1749070-0.5x/pXO1/coverage/SRR1749070-0.5x-anthrax-toxin.pdf)|
| SRR1749070-1x    | 0    | 0.0000  | 0.0000 | 0.0912 | 0.0000  | 16  |[pXO1](/results/cereus-control/SRR1749070-1x/pXO1/coverage/SRR1749070-1x-coverage.pdf), [Toxin](/results/cereus-control/SRR1749070-1x/pXO1/coverage/SRR1749070-1x-anthrax-toxin.pdf)|
| SRR1749070-5x    | 0    | 0.0000  | 0.0000 | 0.4261 | 0.0000  | 53  |[pXO1](/results/cereus-control/SRR1749070-5x/pXO1/coverage/SRR1749070-5x-coverage.pdf), [Toxin](/results/cereus-control/SRR1749070-5x/pXO1/coverage/SRR1749070-5x-anthrax-toxin.pdf)|

To give an idea of the difference in coverages, below is coverage across the pXO1 plasmid at 1000bp
sliding windows with a 500bp overlap for three samples. SRR1748708 was selected because it has the greatest
mean coverage of the NYC samples. The selection of the *B. anthracis* 0.5x coverage control was due to 
it having a similar mean coverage (2.4x vs 2x) as SRR1748708. *B. cereus* 5x coverage control was 
selected because it has the greatest mean coverage (0.4x) of the *B. cereus* controls.
###### *NYC sample SRR1748708*
![SRR1748708](/results/anthracis/SRR1748708/pXO1/coverage/SRR1748708-coverage-1000bp.png "SRR1748708")
###### *B. anthracis* control 0.5x coverage
![B. anthracis control 0.5x coverage](/results/anthracis-control/SRR1749070-0.5x/pXO1/coverage/SRR1749070-0.5x-coverage-1000bp.png "B. anthracis control 0.5x coverage")
###### *B. cereus* control 5x coverage
![B. cereus control 5x coverage](/results/cereus-control/SRR1749070-5x/pXO1/coverage/SRR1749070-5x-coverage-1000bp.png "B. cereus control 5x coverage")

These three plots depict very different stories. Beginning with the *B. cereus* control, although there are few reads 
that mapped, there are regions of pXO1 that are similar to *B. cereus*. This does not explain why SRR1748708 has 
decent coverage accross the pOX1 plasmid. Comparing the SRR1748708 and the *B. anthracis* control, the similar mean
coverages are very different from one another. SRR1748708 has regions of high coverage and regions of low coverage.
The *B. anthracis* control, on the other hand, has a very similar coverage across the complete pXO1 plasmid. This 
suggests SRR1748708 may contain a plasmid backbone that is shared among *Bacillus* species.

Although the pXO1 plasmid has a mean coverage of 2x in SRR1748708, more importantly, are the genes associated with the 
anthrax toxin covered? Below are visualizations of SRR1748708 and the *B. anthracis* 0.5x and 5x controls that 
depict subplots of reads mapped to the *cya*, *lef*, *pagA* and *pagR* genes. For each of the plots below the 
top plot is coverage across the whole pXO1 plasmid. The coverage is based on 1,000bp sliding windows with an overlap 
of 500bp. Each subplot against the genes *pagA*, *pagR*, *lef* and *cya* is the actual coverage in that region.

###### *NYC sample SRR1748708*
![SRR1748708](/results/anthracis/SRR1748708/pXO1/coverage/SRR1748708-anthrax-toxin.png "SRR1748708")
###### *B. anthracis* control 0.5x coverage
![B. anthracis control 0.5x coverage](/results/anthracis-control/SRR1749070-0.5x/pXO1/coverage/SRR1749070-0.5x-anthrax-toxin.png "B. anthracis control 0.5x coverage")
###### *B. anthracis* control 5x coverage
![B. anthracis control 5x coverage](/results/anthracis-control/SRR1749070-5x/pXO1/coverage/SRR1749070-5x-anthrax-toxin.png "B. anthracis control 5x coverage")

The *B. anthracis* controls clearly show a number of reads mapping to each of the genes. The total number of 
reads mapped to each gene for each sample was extracted (Table 3). At 0.5x a total of 197 reads mapped to these 
genes. At 5x the coverage the total numvber of reads mapped jumped to 1,817. clearly evident. Looking at 
SRR1748708, the genes are essentially devoid of mapped reads but reads are mapped to these genes. Although only 72 
reads mapped to these genes it is necessary to investigate these reads further. 

#### Table 3: Number of reads mapped to genes (*cya*, *lef*, *pagA* and *pagR*) related to the anthrax toxin.
###### *NYC Subway System*
| Sample     | *cya* | *lef* | *pagA* | *pagR* |
|------------|-------|-------|--------|--------|
| SRR1748707 | 8     | 2     | 2      | 0      |
| SRR1748708 | 30    | 22    | 18     | 2      |
| SRR1749083 | 4     | 10    | 2      | 0      |
###### *B. anthracis* control
| Sample           | *cya* | *lef* | *pagA* | *pagR* |
|------------------|-------|-------|--------|--------|
| SRR1749070-0x    | 0     | 0     | 0      | 0      |
| SRR1749070-0.25x | 29    | 28    | 24     | 2      |
| SRR1749070-0.5x  | 57    | 76    | 55     | 9      |
| SRR1749070-1x    | 107   | 148   | 108    | 16     |
| SRR1749070-5x    | 544   | 644   | 532    | 97     |
###### *B. cereus* control
| Sample           | *cya* | *lef* | *pagA* | *pagR* |
|------------------|-------|-------|--------|--------|
| SRR1749070-0x    | 0     | 0     | 0      | 0      |
| SRR1749070-0.25x | 0     | 0     | 0      | 0      |
| SRR1749070-0.5x  | 0     | 0     | 0      | 0      |
| SRR1749070-1x    | 0     | 0     | 0      | 0      |
| SRR1749070-5x    | 0     | 0     | 0      | 0      |

The next step was to determine the likely origin of the mapped reads in the each of the genes. In order to do so, 
we extracted these reads and saved them in FASTA format. We then blasted those reads against a local NT database 
(described above). We chose to only retain the top five hits for each read. Our goal was not to determine the 
exact organism of origin, but instead to get an idea of what these reads are mapping to (Table 4). There is a 
clear difference between the results of the NYC sample SRR1748708 and the *B. anthracis* control. The reads from 
the *B. anthracis* controls most often return *B. anthracis* (~85%) as the top hit followed by *B. cereus* (~14.6%). 
The NYC samples do not have a top hit to *B. anthracis* at all. The most common is *B. thuringiensis* followed by a 
number of organisms including humans and other *Bacillus* species. It is important to state that in this case the 
total hit counts are not important. What is important is the stark difference in the organisms between the NYC
samples and the control samples.

#### Table 4: Organism names for each of the top NT blast hits of reads mapped to the anthrax toxin genes.
###### *NYC samples SRR1748707, SRR1748708, and SRR1749083*
| Organism                     |  Hit Count   |
|------------------------------|--------------|
| Bacillus thuringiensis       | 189          |
| Homo sapiens                 | 62           |
| Triticum aestivum            | 45           |
| Bacillus cereus              | 41           |
| Bacillus bombysepticus       | 22           |
| Pan troglodytes              | 17           |
| Cyprinus carpio              | 13           |
| Human DNA                    | 10           |
| Bacillus weihenstephanensis  | 7            |
| Bacillus toyonensis          | 7            |
| Staphylococcus xylosus       | 6            |
| Chlorocebus aethiops         | 4            |
| Nippostrongylus brasiliensis | 3            |
| Haemonchus placei            | 2            |
| Pongo abelii                 | 2            |
| TPA_exp: Homo                | 2            |
| Strongyloides stercoralis    | 1            |
| Diphyllobothrium latum       | 1            |
| Strongyloides papillosus     | 1            |
| Echinostoma caproni          | 1            |
###### *B. anthracis* controls (0.25x, 0.5x, 1x, and 5x)
| Organism            | Hit Count |
|---------------------|-----------|
| Bacillus anthracis  | 10486     |
| Bacillus cereus     | 1803      |
| Synthetic construct | 35        |
| Cyprinus carpio     | 2         |

## pXO2 Results
Similar to pXO1, we calculated the summary statistics of coverage across the complete pXO2 plasmid for
each sample (Table 2). We observe very similar patterns of coverage compared to pXO1. Again NYC sample SRR174708 
has the best mean coverage (4.5x), but each of the NYC samples again have a median coverage of 0x. 
The *B. anthracis* 0.5x control has a mean coverage of 1.34x (median 1x). Again the mean coverage improves 
as the estimated coverage of *B. anthracis* reads increases. Finally, as before, the *B. cereus* control 
includes low levels of pXO2 coverage.

#### Table 5: Summary statistics of per base coverage against pXO2 for NYC samples and controls.
###### *NYC Subway System*
| Study      | Min. | 1st Qu. | Median | Mean   | 3rd Qu. | Max  | Coverage Plots |
|------------|------|---------|--------|--------|---------|------|----------------|
| SRR1749083 | 0    | 0.0000  | 0.0000 | 2.1786 | 0.0000  | 3706 |[pXO2](/results/anthracis/SRR1749083/pXO2/coverage/SRR1749083-coverage.pdf)|
| SRR1748708 | 0    | 0.0000  | 0.0000 | 4.4961 | 0.0000  | 9863 |[pXO2](/results/anthracis/SRR1748708/pXO2/coverage/SRR1748708-coverage.pdf)|
| SRR1748707 | 0    | 0.0000  | 0.0000 | 0.3622 | 0.0000  | 707  |[pXO2](/results/anthracis/SRR1748707/pXO2/coverage/SRR1748707-coverage.pdf)|
###### *B. anthracis* control
| Study            | Min. | 1st Qu. | Median  | Mean    | 3rd Qu. | Max | Coverage Plots |
|------------------|------|---------|---------|---------|---------|-----|----------------|
| SRR1749070-0x    | 0    | 0.0000  | 0.0000  | 0.2483  | 0.0000  | 662 |[pXO2](/results/anthracis-control/SRR1749070-0x/pXO2/coverage/SRR1749070-0x-coverage.pdf)|
| SRR1749070-0.25x | 0    | 0.0000  | 0.0000  | 0.7916  | 1.0000  | 667 |[pXO2](/results/anthracis-control/SRR1749070-0.25x/pXO2/coverage/SRR1749070-0.25x-coverage.pdf)|
| SRR1749070-0.5x  | 0    | 0.0000  | 1.0000  | 1.3448  | 2.0000  | 667 |[pXO2](/results/anthracis-control/SRR1749070-0.5x/pXO2/coverage/SRR1749070-0.5x-coverage.pdf)|
| SRR1749070-1x    | 0    | 1.0000  | 2.0000  | 2.4348  | 3.0000  | 667 |[pXO2](/results/anthracis-control/SRR1749070-1x/pXO2/coverage/SRR1749070-1x-coverage.pdf)|
| SRR1749070-5x    | 0    | 8.0000  | 10.0000 | 11.4254 | 13.0000 | 670 |[pXO2](/results/anthracis-control/SRR1749070-5x/pXO2/coverage/SRR1749070-5x-coverage.pdf)|
###### *B. cereus* control
| Study            | Min. | 1st Qu. | Median | Mean   | 3rd Qu. | Max | Coverage Plots |
|------------------|------|---------|--------|--------|---------|-----|----------------|
| SRR1749070-0x    | 0    | 0.0000  | 0.0000 | 0.2483 | 0.0000  | 662 |[pXO2](/results/cereus-control/SRR1749070-0x/pXO2/coverage/SRR1749070-0x-coverage.pdf)|
| SRR1749070-0.25x | 0    | 0.0000  | 0.0000 | 0.2969 | 0.0000  | 667 |[pXO2](/results/cereus-control/SRR1749070-0.25x/pXO2/coverage/SRR1749070-0.25x-coverage.pdf)|
| SRR1749070-0.5x  | 0    | 0.0000  | 0.0000 | 0.3577 | 0.0000  | 676 |[pXO2](/results/cereus-control/SRR1749070-0.5x/pXO2/coverage/SRR1749070-0.5x-coverage.pdf)|
| SRR1749070-1x    | 0    | 0.0000  | 0.0000 | 0.4805 | 0.0000  | 688 |[pXO2](/results/cereus-control/SRR1749070-1x/pXO2/coverage/SRR1749070-1x-coverage.pdf)|
| SRR1749070-5x    | 0    | 0.0000  | 0.0000 | 1.3704 | 0.0000  | 782 |[pXO2](/results/cereus-control/SRR1749070-5x/pXO2/coverage/SRR1749070-5x-coverage.pdf)|

The coverage across the pXO2 plasmid at 1000bp sliding windows with a 500bp overlap for three 
samples is depicted below. SRR1748708 was selected because it has the greatest
mean coverage of the NYC samples. The selection of the *B. anthracis* 1x coverage control was due to 
it having a similar mean coverage (2.4x vs 4x) as SRR1748708. *B. cereus* 5x coverage control was 
selected because it has the greatest mean coverage (1.37x) of the *B. cereus* controls.

###### *NYC sample SRR1748708*
![SRR1748708](/results/anthracis/SRR1748708/pXO2/coverage/SRR1748708-coverage-1000bp.png "SRR1748708")
###### *B. anthracis* control 1x coverage
![B. anthracis control 1x coverage](/results/anthracis-control/SRR1749070-1x/pXO2/coverage/SRR1749070-1x-coverage-1000bp.png "B. anthracis control 1x coverage")
###### *B. cereus* control 5x coverage
![B. cereus control 5x coverage](/results/cereus-control/SRR1749070-5x/pXO2/coverage/SRR1749070-5x-coverage-1000bp.png "B. cereus control 5x coverage")

There is clearly a peak in each of the plots above. This region is likely a repeat region shared about the 
*Bacillus* species. SRR1748708 mostly only maps to this region. As with pXO1, the *B. anthracis* control 
maintains similar coverage accross the complete pXO2 plasmid. Unlike pXO1, there seems to be a number of 
regions in pXO2 that are very similar to *B. cereus* as depicted by the plot of the *B. cereus* control.
