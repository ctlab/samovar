# set path
cd samovar
SAMOVAR=./
output_dir="samovar_out"

# create output directory
mkdir -p $output_dir/.database

# fetch databases
# correct with https://benlangmead.github.io/aws-indexes/k2
wget -O $output_dir/.database/kraken_minusb.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20250402.tar.gz
wget -O $output_dir/.database/kraken_standard_8.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20250402.tar.gz
wget -O $output_dir/.database/kraken_standard_16.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20250402.tar.gz
wget -O $output_dir/.database/kraken_pluspf_8.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20250402.tar.gz
wget -O $output_dir/.database/kraken_pluspfp_8.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_08gb_20250402.tar.gz

tar -xzf $output_dir/.database/*.tar.gz

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
