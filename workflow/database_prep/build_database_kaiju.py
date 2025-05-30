from samovar.build_database import build_database_kaiju, add_database_kaiju, get_taxonomy_db

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

# Add each genome to the database
for input_file, taxid in zip(input_files, taxids):
    add_database_kaiju(input_file, taxid, db_path="tests_outs/kaiju_db")

# Get taxonomy database
get_taxonomy_db(db_path="tests_outs/kaiju_db")

# Build the final database
build_database_kaiju(db_path="tests_outs/kaiju_db", threads=1, protein=False)
