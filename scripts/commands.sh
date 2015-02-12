# Map SRA run accessions to sample ids used in paper
awk -F"\t" '{print $6"\t"$8}' data/SRP051511_info.txt > data/runs-to-samples.txt

# Map Yersinia and B. anthracis to sample ids
./scripts/extract-pathogens.py data/DataTable5-metaphlan-metadata_v19.txt > data/pathogens-to-samples.txt

# Create Symbolic links fro the pathogens in a separate folder
./scripts/map-pathogens.py data/pathogens-to-samples.txt data/runs-to-samples.txt
