# Directories
mkdir references
mkdir results
mkdir src

# Map SRA run accessions to sample ids used in paper
awk -F"\t" '{print $6"\t"$8}' data/SRP051511_info.txt > data/runs-to-samples.txt

# Map Yersinia and B. anthracis to sample ids
./scripts/extract-pathogens.py data/DataTable5-metaphlan-metadata_v19.txt > data/pathogens-to-samples.txt

# Create Symbolic links fro the pathogens in a separate folder
./scripts/map-pathogens.py data/pathogens-to-samples.txt data/runs-to-samples.txt

# SRA to FASTQ
cat data/SRP051511.txt | xargs -I {} mkdir -p sra-fastq/{}
cat data/SRP051511.txt | xargs -I {} fastq-dump --split-files -O sra-fastq/{} -n 1 -P 15 --gzip /data1/home/rpetit/ncbi/sra/{}.sra

# Convert FASTQ to FASTA
nohup ./scripts/fastq-to-fasta.sh 1> logs/fastq-to-fasta.out 2> logs/fastq-to-fasta.err &

# Build BWA v 0.7.12
cd src/
tar -xjvf bwa-0.7.12.tar.bz2
cd bwa-0.7.12
make
ln -s /data1/home/rpetit/readgp/nyc-subway-metagenome/src/bwa-0.7.12/bwa /data1/home/rpetit/readgp/nyc-subway-metagenome/bin/bwa

# Build bam2fastq v 1.1.0
wget http://gsl.hudsonalpha.org/static/software/bam2fastq-1.1.0.tgz
tar xzvf bam2fastq-1.1.0.tgz
cd bam2fastq-1.1.0/
make
ln -s /data1/home/rpetit/readgp/nyc-subway-metagenome/src/bam2fastq-1.1.0/bam2fastq /data1/home/rpetit/readgp/nyc-subway-metagenome/bin/bam2fastq

# Make bwa index
bin/bwa index -p references/index/CP009540_pXO1 references/CP009540_pXO1.fasta
bin/bwa index -p references/index/NC_005707_pXO1-like references/NC_005707_pXO1-like.fasta
bin/bwa index -p references/index/NC_007323_pXO2 references/NC_007323_pXO2.fasta

# Covert GenBank to GFF3 (From BioPerl)
bp_genbank2gff3 -o references/ references/CP009540_pXO1.gbk

# Map to pXO1 and pXO2
nohup scripts/mapping/map-anthracis-plasmids.sh 1> logs/map-anthracis-plasmids.out 2> logs/map-anthracis-plasmids.err &

# Download Bacillus anthracis
mkdir -p sra-controls/anthracis
prefetch -a "/data1/home/rpetit/.aspera/connect/bin/ascp|/data1/home/rpetit/.aspera/connect/etc/asperaweb_id_dsa.openssh" DRR014739
fastq-dump --split-files -O sra-controls/anthracis --gzip /data1/home/rpetit/ncbi/sra/DRR014739.sra

# Download Bacillus cereus VD142
mkdir -p sra-controls/cereus
prefetch -a "/data1/home/rpetit/.aspera/connect/bin/ascp|/data1/home/rpetit/.aspera/connect/etc/asperaweb_id_dsa.openssh" SRR642775
fastq-dump --split-files -O sra-controls/cereus --gzip /data1/home/rpetit/ncbi/sra/SRR642775.sra


