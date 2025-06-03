import argparse
from samovar.build_database import build_database_from_config

parser = argparse.ArgumentParser(description='Build Kaiju database from config file')
parser.add_argument('--config_path', "-c", help='Path to config YAML file')
parser.add_argument('--db-path', "-o", default="tests_outs/kaiju_db", help='Path to store the database')
args = parser.parse_args()
    
build_database_from_config(args.config_path, db_type="kaiju", db_path=args.db_path)

