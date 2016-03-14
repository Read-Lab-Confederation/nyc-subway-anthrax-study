# k-mer based typing schema for Bacillus cereus group (BCG) and Bacillus anthracis
#
# Analysis Steps
#
# Step 0: Download a few programs for analysis
# Step 1: Download completed genomes BCG and non BCG Bacilli
# Step 2: Use MASH to determine if a non BCG Bacilli should actually be in BCG
# Step 3: Download plasmids and rRNA for Bacillus
# Step 4: Count 31-mers for each of the downloaded sequences
# Step 5: Identify conserved k-mers in BCG
# Step 6: Use k-mers identified in the custom SNP assay, to identify conserved
#         k-mers in Bacillus anthracis
# Step 7: Simulate Illumina reads for both BCG and B. anthracis genomes
# Step 8: Down-select BCG and B. anthracis k-mers based on simulations
# Step 9: Test the down-selected k-mers against NYC controls
# Step 10: Run the down-selected k-mers against NYC samples

TOP_DIR := $(shell pwd)
BIN_DIR := $(TOP_DIR)/bin/
DATA_DIR := $(TOP_DIR)/data/
SRC_DIR := $(TOP_DIR)/src/
NUM_CPU ?= 1
MAX_GENOMES ?= 1000
TOTAL_SIMS ?= 30

all: programs dl-genomes mash dl-sequences count-31mers specific-kmers simulate-reads


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Step 0: Download a few programs for analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
programs: jellyfish art

jellyfish: $(BIN_DIR)/jellyfish
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

art: $(BIN_DIR)/art_illumina
$(SRC_DIR)/art.tar.gz:
	mkdir -p $(SRC_DIR)
	wget -O $@ http://www.niehs.nih.gov/research/resources/assets/docs/artbinchocolatecherrycake031915linux64tgz.tgz

$(BIN_DIR)/art_illumina: $(SRC_DIR)/art.tar.gz
	tar -C $(SRC_DIR) -xzvf $^
	touch $(SRC_DIR)/art_bin_ChocolateCherryCake/art_illumina
	ln -s $(SRC_DIR)/art_bin_ChocolateCherryCake/art_illumina $@




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Step 1: Download completed genomes BCG and non BCG Bacilli
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
GENOMES_DIR := $(DATA_DIR)/completed-genomes
BCG_DIR := $(GENOMES_DIR)/bcg
BACILLUS_DIR := $(GENOMES_DIR)/bacillus
BCG_GENOMES = 'txid86661[Organism:exp] AND "complete genome"[Title] NOT "plasmid"[Title] AND refseq[filter]'
BACILLUS_GENOMES = 'txid1386[Organism:exp] NOT txid86661[Organism:exp] "complete"[Title]NOT "ribosomal"[Title] NOT "WGS"[Keyword] AND refseq[filter]'

dl-genomes: dl-bcg-genomes dl-bacillus-genomes

dl-bcg-genomes: $(BCG_DIR)/completed-genomes.list
$(BCG_DIR)/completed-genomes.list:
	mkdir -p $(BCG_DIR)/fasta
	$(BIN_DIR)/download-genomes.py $(BCG_GENOMES) $(BCG_DIR)/fasta --retmax=$(MAX_GENOMES)
	find $(BCG_DIR)/fasta -empty -type f -delete
	ls $(BCG_DIR)/fasta/ | grep fasta | sed 's/.fasta//' > $(BCG_DIR)/completed-genomes.list

dl-bacillus-genomes: $(BACILLUS_DIR)/completed-genomes.list
$(BACILLUS_DIR)/completed-genomes.list:
	mkdir -p $(BACILLUS_DIR)/fasta
	$(BIN_DIR)/download-genomes.py $(BACILLUS_GENOMES) $(BACILLUS_DIR)/fasta --retmax=$(MAX_GENOMES)
	find $(BACILLUS_DIR)/fasta -empty -type f -delete
	ls $(BACILLUS_DIR)/fasta/ | grep fasta | sed 's/.fasta//' > $(BACILLUS_DIR)/completed-genomes.list



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Step 2: Use MASH to determine if a non BCG Bacilli should actually be in BCG
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mash: mash-sketch mash-self mash-bacillus mash-compare get-bcg-accessions

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
	touch $(BACILLUS_DIR)/completed-genomes.mash
	touch $(BCG_DIR)/completed-genomes.msh
	touch $(BCG_DIR)/completed-genomes.mash
	touch $(BACILLUS_DIR)/completed-genomes.bcg

get-bcg-accessions: $(DATA_DIR)/bcg-accession.txt
$(DATA_DIR)/bcg-accession.txt: $(BACILLUS_DIR)/completed-genomes.bcg
	$(BIN_DIR)/download-genomes.py $(BCG_DIR)/completed-genomes.list NONE --dry_run > $@




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Step 3: Download plasmids and rRNA for Bacillus
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
BACILLUS_QUERY = 'txid1386[Organism:exp] NOT txid86661[Organism:exp] AND "complete"[Title] NOT "ribosomal"[Title] NOT "WGS"[Keyword] AND (biomol_genomic[PROP] AND refseq[filter])'
BCG_QUERY = 'txid86661[Organism:exp] NOT txid1392[Organism:exp] AND "complete"[Title] NOT "ribosomal"[Title] NOT "WGS"[Keyword] AND (biomol_genomic[PROP] AND refseq[filter])'
RRNA_QUERY := 'txid1386[Organism:exp] AND "ribosomal RNA"  NOT "whole genome" NOT "plasmid" NOT "genome" '
SILVA_URL = 'http://www.arb-silva.de/fileadmin/silva_databases/release_123/Exports'

dl-sequences: dl-bacillus-seqeunces dl-bcg-seqeunces dl-rrna

dl-bacillus-seqeunces: $(DATA_DIR)/non-bcg.fasta
$(BACILLUS_DIR)/sequences.list: $(BACILLUS_DIR)/completed-genomes.bcg
	mkdir -p $(BACILLUS_DIR)/sequences
	$(BIN_DIR)/download-genomes.py $(BACILLUS_QUERY) $(BACILLUS_DIR)/sequences --skip $^
	find $(BACILLUS_DIR)/sequences -empty -type f -delete
	ls $(BACILLUS_DIR)/sequences/ | grep fasta | sed 's/.fasta//' > $(BACILLUS_DIR)/sequences.list

$(DATA_DIR)/non-bcg.fasta: $(BACILLUS_DIR)/sequences.list
	cat $(BACILLUS_DIR)/sequences/*.fasta > $(DATA_DIR)/non-bcg.fasta

dl-bcg-seqeunces: $(DATA_DIR)/non-ba.fasta
$(DATA_DIR)/non-ba.fasta: $(DATA_DIR)/non-bcg.fasta
	$(BIN_DIR)/entrez-batch-download.py $(BCG_QUERY) $(DATA_DIR)/bcg.fasta
	cat $(DATA_DIR)/bcg.fasta $^ > $@

dl-rrna: merge-all

dl-ncbi-rrna: $(DATA_DIR)/ncbi-rrna.fasta
$(DATA_DIR)/ncbi-rrna.fasta:
	$(BIN_DIR)/entrez-batch-download.py $(RRNA_QUERY) $@

dl-silva: $(DATA_DIR)/SILVA_123_LSUParc_tax_silva.fasta.gz $(DATA_DIR)/SILVA_123_SSUParc_tax_silva.fasta.gz
$(DATA_DIR)/SILVA_123_LSUParc_tax_silva.fasta.gz:
	wget -O $@ $(SILVA_URL)/SILVA_123_LSUParc_tax_silva.fasta.gz

$(DATA_DIR)/SILVA_123_SSUParc_tax_silva.fasta.gz:
	wget -O $@ $(SILVA_URL)/SILVA_123_SSUParc_tax_silva.fasta.gz

merge-silva: $(DATA_DIR)/SILVA.fasta
$(DATA_DIR)/SILVA.fasta: $(DATA_DIR)/SILVA_123_LSUParc_tax_silva.fasta.gz $(DATA_DIR)/SILVA_123_SSUParc_tax_silva.fasta.gz
	zcat $(DATA_DIR)/SILVA_123_LSUParc_tax_silva.fasta.gz $(DATA_DIR)/SILVA_123_SSUParc_tax_silva.fasta.gz > $@

merge-all: dl-ncbi-rrna merge-silva $(DATA_DIR)/rrna.fasta
$(DATA_DIR)/rrna.fasta: $(DATA_DIR)/ncbi-rrna.fasta $(DATA_DIR)/SILVA.fasta
	cat $(DATA_DIR)/ncbi-rrna.fasta $(DATA_DIR)/SILVA.fasta > $(DATA_DIR)/rrna.fasta




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Step 4: Count 31-mers for each of the downloaded sequences
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
count-31mers: jf-genomes jf-sequences jf-rrna


jf-genomes: jf-bcg-genomes
jf-bcg-genomes: $(BCG_DIR)/completed-genomes.jf
$(BCG_DIR)/completed-genomes.jf: $(BCG_DIR)/completed-genomes.list
	mkdir -p $(BCG_DIR)/jf
	cat $^ | xargs -I {} $(BIN_DIR)/jellyfish count -C -m 31 -s 10M -t $(NUM_CPU) -o $(BCG_DIR)/jf/{}.jf $(BCG_DIR)/fasta/{}.fasta
	ls $(BCG_DIR)/jf/* > $@

jf-sequences: jf-non-bcg-sequences jf-non-ba-sequences
jf-non-bcg-sequences: $(DATA_DIR)/non-bcg.jf
$(DATA_DIR)/non-bcg.jf: $(DATA_DIR)/non-bcg.fasta
	$(BIN_DIR)/jellyfish count -C -m 31 -s 100M -t $(NUM_CPU) -o $@ $^

jf-non-ba-sequences: $(DATA_DIR)/non-ba.jf
$(DATA_DIR)/non-ba.jf: $(DATA_DIR)/non-ba.fasta
	$(BIN_DIR)/jellyfish count -C -m 31 -s 100M -t $(NUM_CPU) -o $@ $^


jf-rrna: $(DATA_DIR)/rrna.jf
$(DATA_DIR)/rrna.jf: $(DATA_DIR)/rrna.fasta
	$(BIN_DIR)/jellyfish count -C -m 31 -s 1M -t $(NUM_CPU) -o $@ $^




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Step 5: Identify conserved k-mers in BCG
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
specific-kmers: bcg-specific-kmers ba-specific-kmers


bcg-specific-kmers: bcg-filter-rrna bcg-compare-kmers bcg-kmers-fasta
bcg-filter-rrna: $(BCG_DIR)/completed-genomes.kmers
$(BCG_DIR)/completed-genomes.kmers: $(BCG_DIR)/completed-genomes.jf $(DATA_DIR)/rrna.jf
	$(BIN_DIR)/common-kmers.py $(BCG_DIR)/jf $(BIN_DIR)/jellyfish $(DATA_DIR)/rrna.jf $@

bcg-compare-kmers: $(DATA_DIR)/bcg-comparison.kmers
$(DATA_DIR)/bcg-comparison.kmers: $(BCG_DIR)/completed-genomes.kmers
	$(BIN_DIR)/jellyfish query -s $^ --load $(DATA_DIR)/non-bcg.jf | sort -n -k2,2 > $@

bcg-kmers-fasta: $(DATA_DIR)/bcg-specific-kmers.fasta
$(DATA_DIR)/bcg-specific-kmers.txt: $(DATA_DIR)/bcg-comparison.kmers
	awk '{if($$2==0){print $$0}}' $^ > $@

$(DATA_DIR)/bcg-specific-kmers.fasta: $(DATA_DIR)/bcg-specific-kmers.txt
	awk '{print ">"$$1"\n"$$1}' $^ > $@




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Step 6: Use k-mers identified in the custom SNP assay, to identify conserved
#         k-mers in Bacillus anthracis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ba-specific-kmers: ba-filter-rrna ba-compare-kmers ba-kmers-fasta
ba-filter-rrna: $(DATA_DIR)/ba-rrna-filtered.kmers
$(DATA_DIR)/ba-rrna-filtered.kmers: $(DATA_DIR)/rrna.jf
	$(BIN_DIR)/common-kmers.py $(TOP_DIR)/../04-custom-ba-snp-assay/results/ba-snp-kmers.txt $(BIN_DIR)/jellyfish $(DATA_DIR)/rrna.jf $@ --single_input

ba-compare-kmers: $(DATA_DIR)/ba-comparison.kmers
$(DATA_DIR)/ba-comparison.kmers: $(DATA_DIR)/ba-rrna-filtered.kmers
	$(BIN_DIR)/jellyfish query -s $^ --load $(DATA_DIR)/non-ba.jf | sort -n -k2,2 > $@

ba-kmers-fasta: $(DATA_DIR)/ba-specific-kmers.fasta
$(DATA_DIR)/ba-filtered-kmers.txt: $(DATA_DIR)/ba-comparison.kmers
	awk '{if($$2==0){print $$0}}' $^ > $@

$(DATA_DIR)/ba-filtered-kmers.fasta: $(DATA_DIR)/ba-filtered-kmers.txt
	awk '{print ">"$$1"\n"$$1}' $^ > $@

$(DATA_DIR)/ba-specific-kmers.txt: $(DATA_DIR)/ba-filtered-kmers.fasta $(BCG_DIR)/completed-genomes.jf $(DATA_DIR)/bcg-accession.txt
	$(BIN_DIR)/query-jellyfish-counts.py $(DATA_DIR)/ba-filtered-kmers.fasta $(BCG_DIR)/jf $(BIN_DIR)/jellyfish $(DATA_DIR)/ba-filtered-stats.txt $(DATA_DIR)/ba-filtered-kmer-by-sample.txt $(DATA_DIR)/ba-filtered-kmer-counts.txt
	gzip -f $(DATA_DIR)/ba-filtered-kmer-by-sample.txt
	$(BIN_DIR)/filter-kmer-by-sample.py $(DATA_DIR)/ba-filtered-kmer-by-sample.txt.gz $(DATA_DIR)/bcg-accession.txt 'bacillus anthracis' | awk '{if ($$11 == 1.00) {print $$1}}' > $@

$(DATA_DIR)/ba-specific-kmers.fasta: $(DATA_DIR)/ba-specific-kmers.txt
	awk '{print ">"$$1"\n"$$1}' $^ > $@




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 7: Simulate Illumina reads for both BCG and B. anthracis genomes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
simulate-reads: simulate-bcg simulate-bacillus

simulate-bcg: $(BCG_DIR)/completed-genomes.sims
$(BCG_DIR)/completed-genomes.sims: $(BCG_DIR)/completed-genomes.list
	mkdir -p $(BCG_DIR)/simulation
	cat $^ | xargs -I {} -P $(NUM_CPU) -n 1 $(BIN_DIR)/simulate_reads.sh {} $(BCG_DIR)/fasta/{}.fasta $(BCG_DIR)/simulation $(TOTAL_SIMS)
	ls $(BCG_DIR)/simulation > $@

simulate-bacillus: $(BACILLUS_DIR)/completed-genomes.sims
$(BACILLUS_DIR)/completed-genomes.sims: $(BACILLUS_DIR)/completed-genomes.list
	mkdir -p $(BACILLUS_DIR)/simulation
	cat $^ | xargs -I {} -P $(NUM_CPU) -n 1 $(BIN_DIR)/simulate_reads.sh {} $(BACILLUS_DIR)/fasta/{}.fasta $(BACILLUS_DIR)/simulation $(TOTAL_SIMS)
	ls $(BACILLUS_DIR)/simulation > $@




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 8: Down-select BCG and B. anthracis k-mers based on simulations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

query-bcg: $(BCG_DIR)/ba-kmer-stats/bcg-sims.txt $(BCG_DIR)/bcg-kmer-stats/bcg-sims.txt
$(BCG_DIR)/ba-kmer-stats/bcg-sims.txt: $(DATA_DIR)/ba-specific-kmers.fasta $(BCG_DIR)/completed-genomes.sims
	mkdir -p $(BCG_DIR)/ba-kmer-stats
	grep "x$$" $(BCG_DIR)/completed-genomes.sims | xargs -I {} -P $(NUM_CPU) -n 1 $(BIN_DIR)/query-jellyfish-counts.sh $(DATA_DIR)/ba-specific-kmers.fasta $(BCG_DIR)/simulation/{} $(BCG_DIR)/ba-kmer-stats/{} $(BIN_DIR)/jellyfish
	touch $@

$(BCG_DIR)/bcg-kmer-stats/bcg-sims.txt: $(DATA_DIR)/bcg-specific-kmers.fasta $(BCG_DIR)/completed-genomes.sims
	mkdir -p $(BCG_DIR)/bcg-kmer-stats
	grep "x$$" $(BCG_DIR)/completed-genomes.sims | xargs -I {} -P $(NUM_CPU) -n 1 $(BIN_DIR)/query-jellyfish-counts.sh $(DATA_DIR)/bcg-specific-kmers.fasta $(BCG_DIR)/simulation/{} $(BCG_DIR)/bcg-kmer-stats/{} $(BIN_DIR)/jellyfish
	touch $@

ba-kmer-stats: $(BCG_DIR)/ba-kmer-stats/filter-kmer-by-sample.txt
$(BCG_DIR)/ba-kmer-stats/filter-kmer-by-sample.txt: $(BCG_DIR)/ba-kmer-stats/bcg-sims.txt
	find $(BCG_DIR)/ba-kmer-stats -type d -name "*x" | xargs -I {} -P $(NUM_CPU) -n 1 $(BIN_DIR)/filter-kmer-by-sample.sh {} $(DATA_DIR)/bcg-accession.txt "bacillus anthracis" ba
	touch $@

aggregate-ba-kmer-stats: $(BCG_DIR)/ba-kmer-stats/aggregate-ba-kmer-stats.txt
$(BCG_DIR)/ba-kmer-stats/aggregate-ba-kmer-stats.txt: $(BCG_DIR)/ba-kmer-stats/filter-kmer-by-sample.txt
	find $(BCG_DIR)/ba-kmer-stats/ -name "*all-ba-kmer-stats.txt" | xargs -I {} -P $(NUM_CPU) -n 1 $(BIN_DIR)/aggregate-kmer-stats.py {}
	touch $@


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 9: Test the down-selected k-mers against NYC controls
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 10: Run the down-selected k-mers against NYC samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
query-nyc: $(DATA_DIR)/nyc-jellyfish/nyc-ba-kmer-stats.txt $(DATA_DIR)/nyc-jellyfish/nyc-bcg-kmer-stats.txt
$(DATA_DIR)/nyc-jellyfish/nyc-bcg-kmer-stats.txt: $(DATA_DIR)/bcg-specific-kmers.fasta
	mkdir -p $(DATA_DIR)/nyc-jellyfish
	$(BIN_DIR)/query-jellyfish-counts.py $^ $(TOP_DIR)/../00-data/sra-jellyfish $@ $(DATA_DIR)/nyc-jellyfish/nyc-bcg-kmer-by-sample.txt $(DATA_DIR)/nyc-jellyfish/nyc-bcg-kmer-counts.txt
	gzip --best $(DATA_DIR)/nyc-jellyfish/nyc-bcg-kmer-by-sample.txt
	gzip --best $(DATA_DIR)/nyc-jellyfish/nyc-bcg-kmer-counts.txt

$(DATA_DIR)/nyc-jellyfish/nyc-ba-kmer-stats.txt: $(DATA_DIR)/ba-specific-kmers.fasta
	mkdir -p $(DATA_DIR)/nyc-jellyfish
	$(BIN_DIR)/query-jellyfish-counts.py $^ $(TOP_DIR)/../00-data/sra-jellyfish $@ $(DATA_DIR)/nyc-jellyfish/nyc-ba-kmer-by-sample.txt $(DATA_DIR)/nyc-jellyfish/nyc-ba-kmer-counts.txt
	gzip --best $(DATA_DIR)/nyc-jellyfish/nyc-ba-kmer-by-sample.txt
	gzip --best $(DATA_DIR)/nyc-jellyfish/nyc-ba-kmer-counts.txt


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 9: Test the down-selected k-mers against NYC controls
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###############################################################################
#### Query BCG and Bacillus whole genome simulated read jellyfish counts
$(DATA_DIR)/bcg-sims-jellyfish/bcg.out: $(DATA_DIR)/bcg-specific-kmers.fasta
	mkdir -p $(DATA_DIR)/bcg-sims-jellyfish
	ls $(BCG_DIR)/simulation/*.txt | xargs -I {} -P ${NUM_CPU} -n 1 $(BIN_DIR)/query-jellyfish-counts.sh $^ $(BCG_DIR)/simulation $(DATA_DIR)/bcg-sims-jellyfish {} bcg $(BIN_DIR)/jellyfish
	ls $(DATA_DIR)/bcg-sims-jellyfish/bcg* > $@

$(DATA_DIR)/bcg-sims-jellyfish/ba.out: $(DATA_DIR)/ba-specific-kmers.fasta
	mkdir -p $(DATA_DIR)/bcg-sims-jellyfish
	ls $(BCG_DIR)/simulation/*.txt | xargs -I {} -P ${NUM_CPU} -n 1 $(BIN_DIR)/query-jellyfish-counts.sh $^ $(BCG_DIR)/simulation $(DATA_DIR)/bcg-sims-jellyfish {} ba $(BIN_DIR)/jellyfish
	ls $(DATA_DIR)/bcg-sims-jellyfish/ba* > $@

query-bcg-sims: $(DATA_DIR)/bcg-sims-jellyfish/bcg.out $(DATA_DIR)/bcg-sims-jellyfish/ba.out

$(DATA_DIR)/bacillus-sims-jellyfish/bcg.out: $(DATA_DIR)/bcg-specific-kmers.fasta
	mkdir -p $(DATA_DIR)/bacillus-sims-jellyfish
	ls $(BACILLUS_DIR)/simulation/*.txt | xargs -I {} -P ${NUM_CPU} -n 1 $(BIN_DIR)/query-jellyfish-counts.sh $^ $(BACILLUS_DIR)/simulation $(DATA_DIR)/bacillus-sims-jellyfish {} bcg $(BIN_DIR)/jellyfish
	ls $(DATA_DIR)/bacillus-sims-jellyfish/bcg* > $@

$(DATA_DIR)/bacillus-sims-jellyfish/ba.out: $(DATA_DIR)/ba-specific-kmers.fasta
	mkdir -p $(DATA_DIR)/bacillus-sims-jellyfish
	ls $(BACILLUS_DIR)/simulation/*.txt | xargs -I {} -P ${NUM_CPU} -n 1 $(BIN_DIR)/query-jellyfish-counts.sh $^ $(BACILLUS_DIR)/simulation $(DATA_DIR)/bacillus-sims-jellyfish {} ba $(BIN_DIR)/jellyfish
	ls $(DATA_DIR)/bacillus-sims-jellyfish/ba* > $@

query-bacillus-sims: $(DATA_DIR)/bacillus-sims-jellyfish/bcg.out $(DATA_DIR)/bacillus-sims-jellyfish/ba.out

query-sims: query-bacillus-sims query-bcg-sims

query: query-nyc query-sims



