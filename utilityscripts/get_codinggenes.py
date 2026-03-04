from Bio import SeqIO
import numpy as np
import argparse

"""
Utility script to extract coding gene IDs from a FASTA file containing protein sequences.
"""


def parse_args():
    parser = argparse.ArgumentParser(description="Get coding gene IDs from FASTA")
    parser.add_argument("--input_file", help="Path to the input FASTA file (protein coding)")
    parser.add_argument("--output_file", help="Path to save coding gene IDs to")
    return parser.parse_args()

def main():
    args = parse_args()
    coding_genes = []
    for record in SeqIO.parse(args.input_file, "fasta"):
        recdesc = record.description.split('|')
        coding_genes.append([i for i in recdesc if 'gene=' in i][0].split('=')[1].strip())
    coding_genes = np.unique(coding_genes)

    with open(args.output_file, "w") as out_f:
        for gene in coding_genes:
            out_f.write(f"{gene}\n")

if __name__ == "__main__":
    main()