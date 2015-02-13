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

# Convert FASTQ to FASTA
nohup ./scripts/fastq-to-fasta.sh 1> logs/fastq-to-fasta.out 2> logs/fastq-to-fasta.err &


# BWA
cd src/
tar -xjvf bwa-0.7.12.tar.bz2
cd bwa-0.7.12
make
ln -s /data1/home/rpetit/readgp/nyc-subway-metagenome/src/bwa-0.7.12/bwa /data1/home/rpetit/readgp/nyc-subway-metagenome/bin/bwa

# Make bwa index
bin/bwa index -p references/index/CP009540_pXO1 references/CP009540_pXO1.fasta
bin/bwa index -p references/index/NC_007323_pXO2 references/NC_007323_pXO2.fasta

# Map to pXO1 and pXO2
nohup ./scripts/mapping/bwa-pOX-align.sh 1> logs/bwa-pOX-align.out 2> logs/bwa-pOX-align.err &
