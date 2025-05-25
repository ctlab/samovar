from samovar.build_database import build_database_kraken2, add_database_kraken2, get_taxonomy_db

input_files = [
    "data/test_genomes/meta/562.fna",
    "data/test_genomes/meta/2886930.fna",
    "data/test_genomes/meta/4932.fna",
    "data/test_genomes/host/9606.fna"
]

taxids = [
    "562",
    "2886930",
    "4932",
    "9606"
]

for input_file, taxid in zip(input_files, taxids):
    add_database_kraken2(input_file, taxid, db_path="tests_outs/kraken_db")

get_taxonomy_db(db_path="tests_outs/kraken_db")

build_database_kraken2(db_path="tests_outs/kraken_db", threads=1, kmer_len=35, minimizer_len=31, minimizer_spaces=7, skip_maps=True)
