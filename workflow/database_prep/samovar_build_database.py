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
    parser.set_defaults(example_omit=None)
    omit = parser.add_mutually_exclusive_group()
    omit.add_argument(
        '--example-omit',
        dest='example_omit',
        action='store_true',
        help='Toy/test only: omit Escherichia from Kraken2 and Phage Phi X from Kaiju',
    )
    omit.add_argument(
        '--no-example-omit',
        dest='example_omit',
        action='store_false',
        help='Do not apply toy taxon gaps even if input_dir is data/test_genomes',
    )
    parser.add_argument('--index', default=None, help='Name to record in the install config databases block')
    parser.add_argument('--flags', default="", help='Extra CLI flags stored with --index')

    args = parser.parse_args()
    # If type is 'kraken', use krakenunique processing
    db_type = 'krakenunique' if args.type == 'kraken' else args.type
    build_database_from_config(
        args.config_path,
        db_type=db_type,
        db_path=args.db_path,
        example_omit=args.example_omit,
    )
    if args.index:
        from samovar.genome_index import register_database
        from samovar.paths import absolute_path

        register_database(db_type, args.index, absolute_path(args.db_path), args.flags or "")
        print(f"Indexed database {db_type}/{args.index} -> {args.db_path}")

if __name__ == '__main__':
    main() 