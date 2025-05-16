from samovar.build_database import build_database, add_database, get_taxonomy_db

input_files = [
    "data/test_genomes/meta/Ecoli.fna",
    "data/test_genomes/meta/Phix.fna",
    "data/test_genomes/meta/Scer.fna",
    "data/test_genomes/host/Hsap.fna"
]

taxids = [
    "562",
    "2886930",
    "4932",
    "9606"
]

for input_file, taxid in zip(input_files, taxids):
    add_database(input_file, taxid, db_path="tests_outs/kraken_db")

get_taxonomy_db(db_path="tests_outs/kraken_db")

build_database(db_path="tests_outs/kraken_db/")
