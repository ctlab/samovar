#!/usr/bin/env python3

import argparse
from samovar.build_database import build_database_from_config

def main():
    parser = argparse.ArgumentParser(description='Build Kaiju or Kraken2 database from config file')
    parser.add_argument('--type', choices=['kaiju', 'kraken2'], required=True,
                      help='Type of database to build (kaiju or kraken2)')
    parser.add_argument('--config_path', required=True,
                      help='Path to config YAML file')
    parser.add_argument('--db_path', required=True,
                      help='Path to store the database')
    
    args = parser.parse_args()
    build_database_from_config(args.config_path, db_type=args.type, db_path=args.db_path)

if __name__ == '__main__':
    main() 