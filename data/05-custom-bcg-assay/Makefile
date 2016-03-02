TOP_DIR := $(shell pwd)
BIN_DIR := $(TOP_DIR)/bin/
DATA_DIR := $(TOP_DIR)/data/
SRC_DIR := $(TOP_DIR)/src/
NUM_CPU ?= 1
TOTAL_SIMS ?= 50

# NCBI Nucleotide queries
# BACILLUS_QUERY - All the genus Bacillus, excluding Bacillus cereus group (BCG)
# BCG_QUERY - All sequences from BCG, excluding B. anthracis
BACILLUS_QUERY = 'txid1386[Organism:exp] NOT txid86661[Organism:exp] AND "complete"[Title] NOT "ribosomal"[Title] NOT "WGS"[Keyword] AND (biomol_genomic[PROP] AND refseq[filter])'
BCG_QUERY = 'txid86661[Organism:exp] NOT txid1392[Organism:exp] AND "complete"[Title] NOT "ribosomal"[Title] NOT "WGS"[Keyword] AND (biomol_genomic[PROP] AND refseq[filter])'

# BCG_GENOMES - All BCG completed genomes
# ANTHACIS_GENOMES - All B. anthracis completed genomes
BCG_GENOMES = 'txid86661[Organism:exp] AND "complete genome"[Title] NOT "plasmid"[Title] AND refseq[filter]'
BACILLUS_GENOMES = 'txid1386[Organism:exp] NOT txid86661[Organism:exp] "complete"[Title]NOT "ribosomal"[Title] NOT "WGS"[Keyword] AND refseq[filter]'

# RRNA_QUERY - Filter out ribosomal RNA from the sequences
RRNA_QUERY := 'txid1386[Organism:exp] AND "ribosomal RNA"  NOT "whole genome" NOT "plasmid" NOT "genome" '


all: programs dl-data mash count-31mers specific-kmers simulate-reads query

programs: jellyfish art

dl-data: dl-genomes dl-sequences dl-rrna

###############################################################################
#### Download Jellyfish
JELLYFISH_DIR := $(SRC_DIR)/jellyfish-2.2.4/
$(SRC_DIR)/jellyfish-2.2.4.tar.gz:
	mkdir -p $(SRC_DIR)
	wget -O $@ https://github.com/gmarcais/Jellyfish/releases/download/v2.2.4/jellyfish-2.2.4.tar.gz

$(BIN_DIR)/jellyfish: $(SRC_DIR)/jellyfish-2.2.4.tar.gz
	tar -C $(SRC_DIR) -xzvf $^
	find $(JELLYFISH_DIR) -exec touch {} \;
	cd $(JELLYFISH_DIR) && autoreconf -f -i && ./configure && cd $(TOP_DIR)
	make -C $(JELLYFISH_DIR)
	ln -s $(JELLYFISH_DIR)/bin/jellyfish $(BIN_DIR)/jellyfish

jellyfish: $(BIN_DIR)/jellyfish


$(SRC_DIR)/art.tar.gz:
	mkdir -p $(SRC_DIR)
	wget -O $@ http://www.niehs.nih.gov/research/resources/assets/docs/artbinchocolatecherrycake031915linux64tgz.tgz

$(BIN_DIR)/art_illumina: $(SRC_DIR)/art.tar.gz
	tar -C $(SRC_DIR) -xzvf $^
	ln -s $(SRC_DIR)/art_bin_ChocolateCherryCake/art_illumina $@

art: $(BIN_DIR)/art_illumina


###############################################################################
#### Download Completed Genomes for BCG, B. anthracis, and non BCG
GENOMES_DIR := $(DATA_DIR)/completed-genomes
BCG_DIR := $(GENOMES_DIR)/bcg
BACILLUS_DIR := $(GENOMES_DIR)/bacillus
MAX_GENOMES ?= 1000

$(BCG_DIR)/completed-genomes.list:
	mkdir -p $(BCG_DIR)/fasta
	$(BIN_DIR)/download-genomes.py $(BCG_GENOMES) $(BCG_DIR)/fasta --retmax=$(MAX_GENOMES)
	find $(BCG_DIR)/fasta -empty -type f -delete
	ls $(BCG_DIR)/fasta/ | grep fasta | sed 's/.fasta//' > $(BCG_DIR)/completed-genomes.list

dl-bcg-genomes: $(BCG_DIR)/completed-genomes.list

$(BACILLUS_DIR)/completed-genomes.list:
	mkdir -p $(BACILLUS_DIR)/fasta
	$(BIN_DIR)/download-genomes.py $(BACILLUS_GENOMES) $(BACILLUS_DIR)/fasta --retmax=$(MAX_GENOMES)
	find $(BACILLUS_DIR)/fasta -empty -type f -delete
	ls $(BACILLUS_DIR)/fasta/ | grep fasta | sed 's/.fasta//' > $(BACILLUS_DIR)/completed-genomes.list

dl-bacillus-genomes: $(BACILLUS_DIR)/completed-genomes.list

dl-genomes: dl-bcg-genomes dl-bacillus-genomes


###############################################################################
#### Identify potential members of BCG using MASH
mash-sketch: $(BCG_DIR)/completed-genomes.msh

$(BCG_DIR)/completed-genomes.msh: $(BCG_DIR)/completed-genomes.list
	mash sketch -k 31 -u -s 10000 -o $(BCG_DIR)/completed-genomes $(BCG_DIR)/fasta/*.fasta

mash-self: $(BCG_DIR)/completed-genomes.mash

$(BCG_DIR)/completed-genomes.mash: $(BCG_DIR)/completed-genomes.msh
	mash dist -u $^ $^ > $@

mash-bacillus: $(BACILLUS_DIR)/completed-genomes.mash

$(BACILLUS_DIR)/completed-genomes.mash: $(BACILLUS_DIR)/completed-genomes.list $(BCG_DIR)/completed-genomes.msh
	mkdir -p $(BACILLUS_DIR)/mash
	cat $(BACILLUS_DIR)/completed-genomes.list | xargs -I {} -P $(NUM_CPU) -n 1 sh -c 'mash dist -u $(BCG_DIR)/completed-genomes.msh $(BACILLUS_DIR)/fasta/{}.fasta > $(BACILLUS_DIR)/mash/{}.mash'
	ls $(BACILLUS_DIR)/mash/* > $(BACILLUS_DIR)/completed-genomes.mash

mash-compare: $(BACILLUS_DIR)/completed-genomes.bcg

$(BACILLUS_DIR)/completed-genomes.bcg: $(BCG_DIR)/completed-genomes.mash $(BACILLUS_DIR)/completed-genomes.mash $(BCG_DIR)/completed-genomes.list
	$(BIN_DIR)/compare-mash-distances.py $(BCG_DIR)/completed-genomes.mash $(BACILLUS_DIR)/completed-genomes.mash $@
	cat $@ | xargs -I {} cp $(BACILLUS_DIR)/fasta/{}.fasta $(BCG_DIR)/fasta/{}.fasta
	cat $@ >> $(BCG_DIR)/completed-genomes.list

mash: mash-sketch mash-self mash-bacillus mash-compare


###############################################################################
#### Create a list of accessions with the organisms name
$(DATA_DIR)/bcg-accession.txt: $(BCG_DIR)/completed-genomes.list
	$(BIN_DIR)/download-genomes.py $(BCG_DIR)/completed-genomes.list NONE --dry_run > $@

bcg-accessions: $(DATA_DIR)/bcg-accession.txt


###############################################################################
#### Download non BCG sequences and non B. anthracis sequences
$(BACILLUS_DIR)/sequences.list: $(BACILLUS_DIR)/completed-genomes.bcg
	mkdir -p $(BACILLUS_DIR)/sequences
	$(BIN_DIR)/download-genomes.py $(BACILLUS_QUERY) $(BACILLUS_DIR)/sequences --skip $^
	find $(BACILLUS_DIR)/sequences -empty -type f -delete
	ls $(BACILLUS_DIR)/sequences/ | grep fasta | sed 's/.fasta//' > $(BACILLUS_DIR)/sequences.list

$(DATA_DIR)/non-bcg.fasta: $(BACILLUS_DIR)/sequences.list
	cat $(BACILLUS_DIR)/sequences/*.fasta > $(DATA_DIR)/non-bcg.fasta

dl-bacillus-seqeunces: $(DATA_DIR)/non-bcg.fasta

$(DATA_DIR)/non-ba.fasta: $(DATA_DIR)/non-bcg.fasta
	$(BIN_DIR)/entrez-batch-download.py $(BCG_QUERY) $(DATA_DIR)/bcg.fasta
	cat $(DATA_DIR)/bcg.fasta $^ > $@

dl-bcg-seqeunces: $(DATA_DIR)/non-ba.fasta

dl-sequences: dl-bacillus-seqeunces dl-bcg-seqeunces


###############################################################################
#### Download ribosomal rRNA sequnces
$(DATA_DIR)/ncbi-rrna.fasta:
	$(BIN_DIR)/entrez-batch-download.py $(RRNA_QUERY) $@

SILVA_URL = 'http://www.arb-silva.de/fileadmin/silva_databases/release_123/Exports'
$(DATA_DIR)/SILVA_123_LSUParc_tax_silva.fasta.gz:
	wget -O $@ $(SILVA_URL)/SILVA_123_LSUParc_tax_silva.fasta.gz

$(DATA_DIR)/SILVA_123_SSUParc_tax_silva.fasta.gz:
	wget -O $@ $(SILVA_URL)/SILVA_123_SSUParc_tax_silva.fasta.gz

$(DATA_DIR)/SILVA.fasta: $(DATA_DIR)/SILVA_123_LSUParc_tax_silva.fasta.gz $(DATA_DIR)/SILVA_123_SSUParc_tax_silva.fasta.gz
	zcat $(DATA_DIR)/SILVA_123_LSUParc_tax_silva.fasta.gz $(DATA_DIR)/SILVA_123_SSUParc_tax_silva.fasta.gz > $@

$(DATA_DIR)/rrna.fasta: $(DATA_DIR)/ncbi-rrna.fasta $(DATA_DIR)/SILVA.fasta
	cat $(DATA_DIR)/ncbi-rrna.fasta $(DATA_DIR)/SILVA.fasta > $(DATA_DIR)/rrna.fasta

dl-rrna: $(DATA_DIR)/rrna.fasta


###############################################################################
#### Use Jellyfish to count 31-mers in BCG
$(BCG_DIR)/completed-genomes.jf: $(BCG_DIR)/completed-genomes.list
	mkdir -p $(BCG_DIR)/jf
	cat $^ | xargs -I {} $(BIN_DIR)/jellyfish count -C -m 31 -s 10M -t $(NUM_CPU) -o $(BCG_DIR)/jf/{}.jf $(BCG_DIR)/fasta/{}.fasta
	ls $(BCG_DIR)/jf/* > $@

$(BCG_DIR)/completed-genomes.dump: $(BCG_DIR)/completed-genomes.jf
	mkdir -p $(BCG_DIR)/dump
	cat $(BCG_DIR)/completed-genomes.list | xargs -I {} $(BIN_DIR)/jellyfish dump -c -o $(BCG_DIR)/dump/{}.dump $(BCG_DIR)/jf/{}.jf
	ls $(BCG_DIR)/dump/* > $@

bcg-genomes: $(BCG_DIR)/completed-genomes.dump

jf-genomes: bcg-genomes


###############################################################################
#### Use Jellyfish to count 31-mers in non BCG and non B. anthracis seqeunces
$(DATA_DIR)/non-bcg.jf: $(DATA_DIR)/non-bcg.fasta
	$(BIN_DIR)/jellyfish count -C -m 31 -s 100M -t $(NUM_CPU) -o $@ $^

non-bcg-sequences: $(DATA_DIR)/non-bcg.jf

$(DATA_DIR)/non-ba.jf: $(DATA_DIR)/non-ba.fasta
	$(BIN_DIR)/jellyfish count -C -m 31 -s 100M -t $(NUM_CPU) -o $@ $^

non-ba-sequences: $(DATA_DIR)/non-ba.jf

jf-sequences: non-bcg-sequences non-ba-sequences


###############################################################################
#### Use Jellyfish to count 31-mers in rRNA sequences
$(DATA_DIR)/rrna.jf: $(DATA_DIR)/rrna.fasta
	$(BIN_DIR)/jellyfish count -C -m 31 -s 1M -t $(NUM_CPU) -o $@ $^

$(DATA_DIR)/rrna.dump: $(DATA_DIR)/rrna.jf
	$(BIN_DIR)/jellyfish dump -c -o $@ $^

jf-rrna: $(DATA_DIR)/rrna.dump

count-31mers: jf-genomes jf-sequences jf-rrna


###############################################################################
#### Identify BCG spcific 31-mers
$(BCG_DIR)/completed-genomes.kmers: $(BCG_DIR)/completed-genomes.dump $(DATA_DIR)/rrna.dump
	$(BIN_DIR)/common-kmers.py $(BCG_DIR)/completed-genomes.dump $(DATA_DIR)/rrna.dump $@

bcg-filter-rrna: $(BCG_DIR)/completed-genomes.kmers

$(DATA_DIR)/bcg-comparison.kmers: $(BCG_DIR)/completed-genomes.kmers
	$(BIN_DIR)/jellyfish query -s $^ --load $(DATA_DIR)/non-bcg.jf | sort -n -k2,2 > $@

bcg-compare-kmers: $(DATA_DIR)/bcg-comparison.kmers

$(DATA_DIR)/bcg-specific-kmers.txt: $(DATA_DIR)/bcg-comparison.kmers
	awk '{if($$2==0){print $$0}}' $^ > $@

$(DATA_DIR)/bcg-specific-kmers.fasta: $(DATA_DIR)/bcg-specific-kmers.txt
	awk '{print ">"$$1"\n"$$1}' $^ > $@

bcg-kmers-fasta: $(DATA_DIR)/bcg-specific-kmers.fasta

bcg-specific-kmers: bcg-filter-rrna bcg-compare-kmers bcg-kmers-fasta


###############################################################################
#### Identify B. anthracis spcific 31-mers
$(DATA_DIR)/ba-rrna-filtered.kmers: $(DATA_DIR)/rrna.dump
	$(BIN_DIR)/common-kmers.py $(TOP_DIR)/ba-snp-kmers.txt $(DATA_DIR)/rrna.dump $@ --single_input

ba-filter-rrna: $(DATA_DIR)/ba-rrna-filtered.kmers

$(DATA_DIR)/ba-comparison.kmers: $(DATA_DIR)/ba-rrna-filtered.kmers
	$(BIN_DIR)/jellyfish query -s $^ --load $(DATA_DIR)/non-ba.jf | sort -n -k2,2 > $@

$(DATA_DIR)/ba-filtered-kmers.txt: $(DATA_DIR)/ba-comparison.kmers
	awk '{if($$2==0){print $$0}}' $^ > $@

$(DATA_DIR)/ba-filtered-kmers.fasta: $(DATA_DIR)/ba-filtered-kmers.txt
	awk '{print ">"$$1"\n"$$1}' $^ > $@

$(DATA_DIR)/ba-specific-kmers.txt: $(DATA_DIR)/ba-filtered-kmers.fasta $(BCG_DIR)/completed-genomes.jf $(DATA_DIR)/bcg-accession.txt
	$(BIN_DIR)/query-jellyfish-counts.py $(DATA_DIR)/ba-filtered-kmers.fasta $(BCG_DIR)/jf $(DATA_DIR)/ba-filtered-stats.txt $(DATA_DIR)/ba-filtered-kmer-by-sample.txt $(DATA_DIR)/ba-filtered-kmer-counts.txt
	gzip -f $(DATA_DIR)/ba-filtered-kmer-by-sample.txt
	$(BIN_DIR)/filter-kmer-by-sample.py $(DATA_DIR)/ba-filtered-kmer-by-sample.txt.gz $(DATA_DIR)/bcg-accession.txt 'bacillus anthracis' | awk '{if ($$11 == 1.00) {print $$1}}' > $@

$(DATA_DIR)/ba-specific-kmers.fasta: $(DATA_DIR)/ba-specific-kmers.txt
	awk '{print ">"$$1"\n"$$1}' $^ > $@

ba-kmers-fasta: $(DATA_DIR)/ba-specific-kmers.fasta

ba-specific-kmers: ba-filter-rrna ba-kmers-fasta

specific-kmers: bcg-specific-kmers ba-specific-kmers

###############################################################################
#### Query NYC jellyfish counts
$(DATA_DIR)/nyc-jellyfish/nyc-bcg-kmer-stats.txt: $(DATA_DIR)/bcg-specific-kmers.fasta
	mkdir -p $(DATA_DIR)/nyc-jellyfish
	$(BIN_DIR)/query-jellyfish-counts.py $^ $(TOP_DIR)/../00-data/sra-jellyfish $@ $(DATA_DIR)/nyc-jellyfish/nyc-bcg-kmer-by-sample.txt $(DATA_DIR)/nyc-jellyfish/nyc-bcg-kmer-counts.txt

$(DATA_DIR)/nyc-jellyfish/nyc-ba-kmer-stats.txt: $(DATA_DIR)/ba-specific-kmers.fasta
	mkdir -p $(DATA_DIR)/nyc-jellyfish
	$(BIN_DIR)/query-jellyfish-counts.py $^ $(TOP_DIR)/../00-data/sra-jellyfish $@ $(DATA_DIR)/nyc-jellyfish/nycba-kmer-by-sample.txt $(DATA_DIR)/nyc-jellyfish/nyc-ba-kmer-counts.txt

query-nyc: $(DATA_DIR)/nyc-jellyfish/nyc-ba-kmer-stats.txt $(DATA_DIR)/nyc-jellyfish/nyc-bcg-kmer-stats.txt

###############################################################################
#### Simulate Illumina reads for completed BCG genomes
$(BCG_DIR)/completed-genomes.sims: $(BCG_DIR)/completed-genomes.list
	mkdir -p $(BCG_DIR)/simulation
	cat $^ | xargs -I {} -P $(NUM_CPU) -n 1 $(BIN_DIR)/simulate_reads.sh {} $(BCG_DIR)/fasta/{}.fasta $(BCG_DIR)/simulation
	ls $(BCG_DIR)/simulation > $@

simulate-bcg: $(BCG_DIR)/completed-genomes.sims

$(BACILLUS_DIR)/completed-genomes.sims: $(BACILLUS_DIR)/completed-genomes.list
	mkdir -p $(BACILLUS_DIR)/simulation
	cat $^ | xargs -I {} -P $(NUM_CPU) -n 1 $(BIN_DIR)/simulate_reads.sh {} $(BACILLUS_DIR)/fasta/{}.fasta $(BACILLUS_DIR)/simulation
	ls $(BACILLUS_DIR)/simulation > $@

simulate-ba: $(BACILLUS_DIR)/completed-genomes.sims

simulate-reads: simulate-ba simulate-bcg


###############################################################################
#### Query BCG and Bacillus whole genome simulated read jellyfish counts
$(DATA_DIR)/bcg-sims-jellyfish/bcg.out: $(DATA_DIR)/bcg-specific-kmers.fasta
	mkdir -p $(DATA_DIR)/bcg-sims-jellyfish
	ls $(BCG_DIR)/simulation/*.txt | xargs -I {} -P ${NUM_CPU} -n 1 $(BIN_DIR)/query-jellyfish-counts.sh $^ $(BCG_DIR)/simulation $(DATA_DIR)/bcg-sims-jellyfish {} bcg
	ls $(DATA_DIR)/bcg-sims-jellyfish/bcg* > $@

$(DATA_DIR)/bcg-sims-jellyfish/ba.out: $(DATA_DIR)/ba-specific-kmers.fasta
	mkdir -p $(DATA_DIR)/bcg-sims-jellyfish
	ls $(BCG_DIR)/simulation/*.txt | xargs -I {} -P ${NUM_CPU} -n 1 $(BIN_DIR)/query-jellyfish-counts.sh $^ $(BCG_DIR)/simulation $(DATA_DIR)/bcg-sims-jellyfish {} ba
	ls $(DATA_DIR)/bcg-sims-jellyfish/ba* > $@

query-bcg-sims: $(DATA_DIR)/bcg-sims-jellyfish/bcg.out $(DATA_DIR)/bcg-sims-jellyfish/ba.out

$(DATA_DIR)/bacillus-sims-jellyfish/bcg.out: $(DATA_DIR)/bcg-specific-kmers.fasta
	mkdir -p $(DATA_DIR)/bacillus-sims-jellyfish
	ls $(BACILLUS_DIR)/simulation/*.txt | xargs -I {} -P ${NUM_CPU} -n 1 $(BIN_DIR)/query-jellyfish-counts.sh $^ $(BACILLUS_DIR)/simulation $(DATA_DIR)/bacillus-sims-jellyfish {} bcg
	ls $(DATA_DIR)/bacillus-sims-jellyfish/bcg* > $@

$(DATA_DIR)/bacillus-sims-jellyfish/ba.out: $(DATA_DIR)/ba-specific-kmers.fasta
	mkdir -p $(DATA_DIR)/bacillus-sims-jellyfish
	ls $(BACILLUS_DIR)/simulation/*.txt | xargs -I {} -P ${NUM_CPU} -n 1 $(BIN_DIR)/query-jellyfish-counts.sh $^ $(BACILLUS_DIR)/simulation $(DATA_DIR)/bacillus-sims-jellyfish {} ba
	ls $(DATA_DIR)/bacillus-sims-jellyfish/ba* > $@

query-bacillus-sims: $(DATA_DIR)/bacillus-sims-jellyfish/bcg.out $(DATA_DIR)/bacillus-sims-jellyfish/ba.out

query-sims: query-bacillus-sims query-bcg-sims

query: query-nyc query-sims


###############################################################################
#### Simulate Illumina reads for completed BCG genomes multiple times
resimulate-bcg-0.25x: $(BCG_DIR)/bcg-0.25x-simulations.txt
$(BCG_DIR)/bcg-0.25x-simulations.txt: $(BCG_DIR)/completed-genomes.list
	mkdir -p $(BCG_DIR)/bcg-0.25x-simulations
	cat $^ | xargs -I {} -P $(NUM_CPU) -n 1 $(BIN_DIR)/simulate_coverage.sh {} $(BCG_DIR)/fasta/{}.fasta $(BCG_DIR)/bcg-0.25x-simulations 0.25 $(TOTAL_SIMS)
	ls $(BCG_DIR)/bcg-0.25x-simulations > $@

resimulate-bcg-0.5x: $(BCG_DIR)/bcg-0.5x-simulations.txt
$(BCG_DIR)/bcg-0.5x-simulations.txt: $(BCG_DIR)/completed-genomes.list
	mkdir -p $(BCG_DIR)/bcg-0.5x-simulations
	cat $^ | xargs -I {} -P $(NUM_CPU) -n 1 $(BIN_DIR)/simulate_coverage.sh {} $(BCG_DIR)/fasta/{}.fasta $(BCG_DIR)/bcg-0.5x-simulations 0.5 $(TOTAL_SIMS)
	ls $(BCG_DIR)/bcg-0.5x-simulations > $@

resimulate-bcg-1x: $(BCG_DIR)/bcg-1x-simulations.txt
$(BCG_DIR)/bcg-1x-simulations.txt: $(BCG_DIR)/completed-genomes.list
	mkdir -p $(BCG_DIR)/bcg-1x-simulations
	cat $^ | xargs -I {} -P $(NUM_CPU) -n 1 $(BIN_DIR)/simulate_coverage.sh {} $(BCG_DIR)/fasta/{}.fasta $(BCG_DIR)/bcg-1x-simulations 1 $(TOTAL_SIMS)
	ls $(BCG_DIR)/bcg-1x-simulations > $@

resimulate-bcg-5x: $(BCG_DIR)/bcg-5x-simulations.txt
$(BCG_DIR)/bcg-5x-simulations.txt: $(BCG_DIR)/completed-genomes.list
	mkdir -p $(BCG_DIR)/bcg-5x-simulations
	cat $^ | xargs -I {} -P $(NUM_CPU) -n 1 $(BIN_DIR)/simulate_coverage.sh {} $(BCG_DIR)/fasta/{}.fasta $(BCG_DIR)/bcg-5x-simulations 5 $(TOTAL_SIMS)
	ls $(BCG_DIR)/bcg-5x-simulations > $@

resimulate-bcg-10x: $(BCG_DIR)/bcg-10x-simulations.txt
$(BCG_DIR)/bcg-10x-simulations.txt: $(BCG_DIR)/completed-genomes.list
	mkdir -p $(BCG_DIR)/bcg-10x-simulations
	cat $^ | xargs -I {} -P $(NUM_CPU) -n 1 $(BIN_DIR)/simulate_coverage.sh {} $(BCG_DIR)/fasta/{}.fasta $(BCG_DIR)/bcg-10x-simulations 10 $(TOTAL_SIMS)
	ls $(BCG_DIR)/bcg-10x-simulations > $@

resimulate-bcg-25x: $(BCG_DIR)/bcg-25x-simulations.txt
$(BCG_DIR)/bcg-25x-simulations.txt: $(BCG_DIR)/completed-genomes.list
	mkdir -p $(BCG_DIR)/bcg-25x-simulations
	cat $^ | xargs -I {} -P $(NUM_CPU) -n 1 $(BIN_DIR)/simulate_coverage.sh {} $(BCG_DIR)/fasta/{}.fasta $(BCG_DIR)/bcg-25x-simulations 25 $(TOTAL_SIMS)
	ls $(BCG_DIR)/bcg-25x-simulations > $@

resimulate-bcg: resimulate-bcg-0.25x resimulate-bcg-0.5x resimulate-bcg-1x
