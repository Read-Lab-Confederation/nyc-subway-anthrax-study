# Searching for pXO1 and pXO2
  
    ADD SUMMARY OF PXO1 and PXO2
    PXO1 is responsible for encoding the components of the anthrax toxin. The genes required 
    to produce the anthrax toxin are *pagA*, *pagR*, *lef* and *cya*. 



## Identifying NYC Samples samples with evidence of *B. anthracis*
Below is a table of samples from the NYC Subway System dataset which were determined to have contained
evidence for *B. anthracis*. For the remainder of these analysis, each of these samples will be 
referred to by their SRA run accession.

#### Table 1: NYC Subway System samples with evidence of *B. anthracis*.
| NYC Study Sample | SRA Run Accession      |
|------------------|------------------------|
| P00134           | SRR1748707, SRR1748708 |
| P00497           | SRR1749083             |

## Mapping against pXO1 and pXO2 using BWA
We used *BWA (v 0.7.5a-r405)* to map reads from each NYC sample and control sample against the reference 
pXO1 and pXO2 plasmids. For pXO1 we used reference [CP009540](http://www.ncbi.nlm.nih.gov/nuccore/CP009540.1), 
and for pXO2 we used reference [NC_007323](http://www.ncbi.nlm.nih.gov/nuccore/50163691). This process was 
automated using the following scripts *[map-anthracis-plasmids.sh](/scripts/mapping/map-anthracis-plasmids.sh)* 
(NYC samples) and *[map-anthracis-plasmids.sh](/scripts/mapping/map-anthracis-controls.sh)* (control samples).

## pXO1 Results
For each sample we calcualted the summary statistics of coverage accross the complete pXO1 plasmid (Table 2).
NYC sample SRR174708 had the best average coverage (2x), but each of the NYC samples had a median coverage of 
0x. 

    For each of the plots below the top plot is against the whole pXO1 plasmid. The coverage is based on 1,000bp 
    sliding windows with an overlap of 500bp. Each subplot against the genes *pagA*, *pagR*, *lef* and *cya* is 
    the actual coverage in that region.

#### Table 2: Summary statistics of per base coverage against pXO1 for NYC Subway System samples and controls.
###### *NYC Subway System*
| Study      | Min. | 1st Qu. | Median | Mean   | 3rd Qu. | Max | Coverage Plots |
|------------|------|---------|--------|--------|---------|-----|----------------|
| SRR1749083 | 0    | 0.0000  | 0.0000 | 0.4652 | 0.0000  | 44  |[pXO1](/results/anthracis/SRR1749083/pXO1/coverage/SRR1749083-coverage.pdf), [Toxins](/results/anthracis/SRR1749083/pXO1/coverage/SRR1749083-anthrax-toxin.pdf)|
| SRR1748708 | 0    | 0.0000  | 0.0000 | 2.0766 | 4.0000  | 100 |[pXO1](/results/anthracis/SRR1748708/pXO1/coverage/SRR1748708-coverage.pdf), [Toxins](/results/anthracis/SRR1748708/pXO1/coverage/SRR1748708-anthrax-toxin.pdf)|
| SRR1748707 | 0    | 0.0000  | 0.0000 | 0.2506 | 0.0000  | 44  |[pXO1](/results/anthracis/SRR1748707/pXO1/coverage/SRR1748707-coverage.pdf), [Toxins](/results/anthracis/SRR1748707/pXO1/coverage/SRR1748707-anthrax-toxin.pdf)|
###### *B. anthracis* control
| Study            | Min. | 1st Qu. | Median  | Mean    | 3rd Qu. | Max | Coverage Plots |
|------------------|------|---------|---------|---------|---------|-----|----------------|
| SRR1749070-0x    | 0    | 0.0000  | 0.0000  | 0.0103  | 0.0000  | 2   |[pXO1](/results/anthracis-control/SRR1749070-0x/pXO1/coverage/SRR1749070-0x-coverage.pdf), [Toxins](/results/anthracis-control/SRR1749070-0x/pXO1/coverage/SRR1749070-0x-anthrax-toxin.pdf)|
| SRR1749070-0.25x | 0    | 0.0000  | 1.0000  | 1.2305  | 2.0000  | 13  |[pXO1](/results/anthracis-control/SRR1749070-0.25x/pXO1/coverage/SRR1749070-0.25x-coverage.pdf), [Toxins](/results/anthracis-control/SRR1749070-0.25x/pXO1/coverage/SRR1749070-0.25x-anthrax-toxin.pdf)|
| SRR1749070-0.5x  | 0    | 1.0000  | 2.0000  | 2.4250  | 3.0000  | 20  |[pXO1](/results/anthracis-control/SRR1749070-0.5x/pXO1/coverage/SRR1749070-0.5x-coverage.pdf), [Toxins](/results/anthracis-control/SRR1749070-0.5x/pXO1/coverage/SRR1749070-0.5x-anthrax-toxin.pdf)|
| SRR1749070-1x    | 0    | 3.0000  | 4.0000  | 4.8863  | 6.0000  | 36  |[pXO1](/results/anthracis-control/SRR1749070-1x/pXO1/coverage/SRR1749070-1x-coverage.pdf), [Toxins](/results/anthracis-control/SRR1749070-1x/pXO1/coverage/SRR1749070-1x-anthrax-toxin.pdf)|
| SRR1749070-5x    | 0    | 19.0000 | 23.0000 | 24.3656 | 28.0000 | 175 |[pXO1](/results/anthracis-control/SRR1749070-5x/pXO1/coverage/SRR1749070-5x-coverage.pdf), [Toxins](/results/anthracis-control/SRR1749070-5x/pXO1/coverage/SRR1749070-5x-anthrax-toxin.pdf)|
###### *B. cereus* control
| Study            | Min. | 1st Qu. | Median | Mean   | 3rd Qu. | Max | Coverage Plots |
|------------------|------|---------|--------|--------|---------|-----|----------------|
| SRR1749070-0x    | 0    | 0.0000  | 0.0000 | 0.0103 | 0.0000  | 2   |[pXO1](/results/cereus-control/SRR1749070-0x/pXO1/coverage/SRR1749070-0x-coverage.pdf), [Toxins](/results/cereus-control/SRR1749070-0x/pXO1/coverage/SRR1749070-0x-anthrax-toxin.pdf)|
| SRR1749070-0.25x | 0    | 0.0000  | 0.0000 | 0.0368 | 0.0000  | 7   |[pXO1](/results/cereus-control/SRR1749070-0.25x/pXO1/coverage/SRR1749070-0.25x-coverage.pdf), [Toxins](/results/cereus-control/SRR1749070-0.25x/pXO1/coverage/SRR1749070-0.25x-anthrax-toxin.pdf)|
| SRR1749070-0.5x  | 0    | 0.0000  | 0.0000 | 0.0522 | 0.0000  | 9   |[pXO1](/results/cereus-control/SRR1749070-0.5x/pXO1/coverage/SRR1749070-0.5x-coverage.pdf), [Toxins](/results/cereus-control/SRR1749070-0.5x/pXO1/coverage/SRR1749070-0.5x-anthrax-toxin.pdf)|
| SRR1749070-1x    | 0    | 0.0000  | 0.0000 | 0.0912 | 0.0000  | 16  |[pXO1](/results/cereus-control/SRR1749070-1x/pXO1/coverage/SRR1749070-1x-coverage.pdf), [Toxins](/results/cereus-control/SRR1749070-1x/pXO1/coverage/SRR1749070-1x-anthrax-toxin.pdf)|
| SRR1749070-5x    | 0    | 0.0000  | 0.0000 | 0.4261 | 0.0000  | 53  |[pXO1](/results/cereus-control/SRR1749070-5x/pXO1/coverage/SRR1749070-5x-coverage.pdf), [Toxins](/results/cereus-control/SRR1749070-5x/pXO1/coverage/SRR1749070-5x-anthrax-toxin.pdf)|

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

#### Table 4: Organisms of top blast hits of reads mapped to the anthrax toxin genes.
###### *NYC Subway System*
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
###### *B. anthracis* control
| Organism            | Hit Count |
|---------------------|-----------|
| Bacillus anthracis  | 10486     |
| Bacillus cereus     | 1803      |
| Synthetic construct | 35        |
| Cyprinus carpio     | 2         |


    ADD PLOTS OF CONTROL AND COMPARE/CONTRAST
![P00134 (Run: SRR1748708)](/results/anthracis/SRR1748708/pXO1/coverage/SRR1748708-anthrax-toxin.png "SRR1748708")


