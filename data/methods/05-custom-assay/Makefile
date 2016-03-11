TOP_DIR := $(shell pwd)
BIN_DIR := $(TOP_DIR)/bin/
DATA_DIR := $(TOP_DIR)/data/
SRC_DIR := $(TOP_DIR)/src/
NUM_CPU ?= 1

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


all: programs dl-data mash count-31mers specific-kmers

programs: jellyfish

dl-data: dl-genomes dl-sequences dl-rrna

###############################################################################
#### Download Jellyfish
JELLYFISH_DIR := $(SRC_DIR)/jellyfish-2.2.4/
$(SRC_DIR)/jellyfish-2.2.4.tar.gz:
	mkdir -p $(SRC_DIR)
	wget -O $@ https://github.com/gmarcais/Jellyfish/releases/download/v2.2.4/jellyfish-2.2.4.tar.gz

$(JELLYFISH_DIR)/bin/jellyfish: $(SRC_DIR)/jellyfish-2.2.4.tar.gz
	tar -C $(SRC_DIR) -xzvf $^
	find $(JELLYFISH_DIR) -exec touch {} \;
	cd $(JELLYFISH_DIR) && autoreconf -f -i && ./configure && cd $(TOP_DIR)
	make -C $(JELLYFISH_DIR)
	ln -s $(JELLYFISH_DIR)/bin/jellyfish $(BIN_DIR)/jellyfish

jellyfish: $(JELLYFISH_DIR)/bin/jellyfish


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
	awk '{if($2==0){print $0}}' $^ > $@

$(DATA_DIR)/bcg-specific-kmers.fasta: $(DATA_DIR)/bcg-specific-kmers.txt
	awk '{print ">"$1"\n"$2}' $^ > $@

bcg-kmers-fasta: $(DATA_DIR)/bcg-specific-kmers.fasta

bcg-specific-kmers: bcg-filter-rrna bcg-compare-kmers bcg-kmers-fasta


###############################################################################
#### Identify B. anthracis spcific 31-mers
$(DATA_DIR)/ba-rrna-filtered.kmers: $(DATA_DIR)/rrna.dump
	$(BIN_DIR)/common-kmers.py $(TOP_DIR)/ba-specific.kmers $(DATA_DIR)/rrna.dump $@ --single_input

ba-filter-rrna: $(DATA_DIR)/ba-rrna-filtered.kmers

$(DATA_DIR)/ba-comparison.kmers: $(DATA_DIR)/ba-rrna-filtered.kmers
	$(BIN_DIR)/jellyfish query -s $^ --load $(DATA_DIR)/non-ba.jf | sort -n -k2,2 > $@

ba-compare-kmers: $(DATA_DIR)/ba-comparison.kmers

$(DATA_DIR)/ba-specific-kmers.txt: $(DATA_DIR)/ba-comparison.kmers
	awk '{print $0}' $^ > $@

$(DATA_DIR)/ba-specific-kmers.fasta: $(DATA_DIR)/ba-specific-kmers.txt
	awk '{if($2==0){print ">"$1"\n"$2}}' $^ > $@

ba-kmers-fasta: $(DATA_DIR)/ba-specific-kmers.fasta

ba-specific-kmers: ba-filter-rrna ba-compare-kmers ba-kmers-fasta

specific-kmers: bcg-specific-kmers ba-specific-kmers

###############################################################################
#### Query NYC jellyfish counts
$(DATA_DIR)/nyc-jellyfish/bcg-kmer.stats: $(DATA_DIR)/bcg-specific-kmers.fasta
	mkdir -p $(DATA_DIR)/nyc-jellyfish
	$(BIN_DIR)


