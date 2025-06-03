# set path
cd samovar
SAMOVAR=./
output_dir="samovar_toy"

# create output directory
rm -r $output_dir/
mkdir -p $output_dir/.database

# prepare databases
## optional: prepare toy databases
### write config file
echo """
input_dir:
  - "$SAMOVAR/data/test_genomes/meta"
  - "$SAMOVAR/data/test_genomes/host"
output_dir: "$output_dir/database_prep"
mutation_rate: 0.02 # similarity at ~ species level
include_percent: 70.0 # reads from 70% of genome will be included in the metagenome
""" > $output_dir/.database/config.yaml

### build databases
samovar build_database --type kraken2 --config_path $output_dir/.database/config.yaml --db_path $output_dir/.database/kraken2_db

## or fetch databases
mkdir -p $output_dir/.database/kaiju_db
wget -O $output_dir/.database/kaiju_db/kaiju_fungi.tar.gz https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_fungi_2024-08-16.tgz
tar -xzf $output_dir/.database/kaiju_db/kaiju_fungi.tar.gz -C $output_dir/.database/kaiju_db

# Prepare generation config
samovar generate \
    --genome_dir $SAMOVAR/data/test_genomes/meta \
    --host_genome $SAMOVAR/data/test_genomes/host/9606.fna \
    --output_dir $output_dir

# Prepare run config
samovar preprocess \
    --output_dir $output_dir \
    --kraken2-test "kraken2 $output_dir/.database/kraken2_db" \
    --kaiju-test "kaiju $output_dir/.database/kaiju_db"

# Run samovar
samovar exec --output_dir $output_dir
