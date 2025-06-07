# set path
cd samovar
SAMOVAR=./
output_dir="samovar_10bac"

# create output directory
mkdir -p $output_dir/.database

# path to databases
DB_KRAKEN=/mnt/metagenomics/kraken/minikraken_4GB_oct2017
DB_KAIJU=/mnt/metagenomics/kaiju/kaiju_db_refseq.fmi
DB_KRAKEN2=/mnt/metagenomics/kraken2/db_standard_03_2023

# Get genomes
python -m samovar.genome_fetcher \
    --output-dir $output_dir/.genomes \
    --N 10 \
    --group "Bacteria" \
    --silent

# Prepare generation config
samovar generate \
    --genome_dir $output_dir/.genomes \
    --host_genome $SAMOVAR/data/test_genomes/host/9606.fna \
    --output_dir $output_dir

# Prepare run config
samovar preprocess \
    #--input_dir tests/data/reads \
    --output_dir $output_dir \
    --kraken2 "kraken2 $DB_KRAKEN2" \
    --kaiju "kaiju $DB_KAIJU" \
    --kraken "kraken $DB_KRAKEN" \
    --krakenuniq "krakenuniq $DB_KRAKEN"


# Run samovar
samovar exec --output_dir $output_dir
