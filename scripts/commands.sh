# Directories
mkdir bin
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

# Clone fastq-stats
cd src/
git clone git@github.com:rpetit3/fastq-stats.git
cd fastq-stats
git reset --hard 8460e6bbdb09152048b957a8d9dfd1619ff7dbfd
make
ln -s /data1/home/rpetit/readgp/nyc-subway-metagenome/src/fastq-stats/fastq-stats /data1/home/rpetit/readgp/nyc-subway-metagenome/bin/fastq-stats
cd ../

# clone seqtk
git clone git@github.com:lh3/seqtk.git
cd seqtk/
git reset --hard 43ff625a3211b51f301cb356a34fb8d1e593d50a
make
ln -s /data1/home/rpetit/readgp/nyc-subway-metagenome/src/seqtk/seqtk /data1/home/rpetit/readgp/nyc-subway-metagenome/bin/seqtk
cd ../

# Build BWA v 0.7.12
tar -xjvf bwa-0.7.12.tar.bz2
cd bwa-0.7.12
make
ln -s /data1/home/rpetit/readgp/nyc-subway-metagenome/src/bwa-0.7.12/bwa /data1/home/rpetit/readgp/nyc-subway-metagenome/bin/bwa
cd ../

# Build bam2fastq v 1.1.0
wget http://gsl.hudsonalpha.org/static/software/bam2fastq-1.1.0.tgz
tar xzvf bam2fastq-1.1.0.tgz
cd bam2fastq-1.1.0/
make
ln -s /data1/home/rpetit/readgp/nyc-subway-metagenome/src/bam2fastq-1.1.0/bam2fastq /data1/home/rpetit/readgp/nyc-subway-metagenome/bin/bam2fastq
cd ../../

# Make bwa index
bin/bwa index -p references/index/CP009540-pXO1 references/CP009540-pXO1.fasta
bin/bwa index -p references/index/NC_005707-pXO1-like references/NC_005707-pXO1-like.fasta
bin/bwa index -p references/index/NC_007323-pXO2 references/NC_007323-pXO2.fasta
bin/bwa index -p references/index/AL117211-pMT1 references/AL117211-pMT1.fasta

# Covert GenBank to GFF3 (From BioPerl)
bp_genbank2gff3 -o references/ references/CP009540-pXO1.gbk
bp_genbank2gff3 -o references/ references/AL117211-pMT1.gbk

# Map to pXO1 and pXO2
nohup scripts/mapping/map-anthracis-plasmids.sh 1> logs/map-anthracis-plasmids.out 2> logs/map-anthracis-plasmids.err &

# Map Bacillus controls to pXO1 and pXO2
nohup scripts/mapping/map-anthracis-controls.sh 1> logs/map-anthracis-controls.out 2> logs/map-anthracis-controls.err &

# Map to pMT1
nohup scripts/mapping/map-pestis-pmt1.sh 1> logs/map-pestis-pmt1.out 2> logs/map-pestis-pmt1.err &
