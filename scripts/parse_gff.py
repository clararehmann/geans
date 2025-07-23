import gffutils
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Process GFF files to TSV format with relevant features")
    parser.add_argument("input_file", help="Path to the input GFF file")
    parser.add_argument("output_file", help="Path to save GFF database to")
    return parser.parse_args()

def main():
    args = parse_args()
    # Create a database from the GFF file
    db = gffutils.create_db(args.input_file, args.output_file, force=True, keep_order=True, merge_strategy="create_unique")

if __name__ == "__main__":
    main()