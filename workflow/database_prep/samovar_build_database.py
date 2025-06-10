#!/usr/bin/env python3

import argparse
from samovar.build_database import build_database_from_config

def main():
    parser = argparse.ArgumentParser(description='Build Kaiju, Kraken2, Kraken, or KrakenUniq database from config file')
    parser.add_argument('--type', choices=['kaiju', 'kraken2', 'kraken', 'krakenunique'], required=True,
                      help='Type of database to build (kaiju, kraken2, kraken, or krakenunique)')
    parser.add_argument('--config_path', required=True,
                      help='Path to config YAML file')
    parser.add_argument('--db_path', required=True,
                      help='Path to store the database')
    
    args = parser.parse_args()
    # If type is 'kraken', use krakenunique processing
    db_type = 'krakenunique' if args.type == 'kraken' else args.type
    build_database_from_config(args.config_path, db_type=db_type, db_path=args.db_path)

if __name__ == '__main__':
    main() 