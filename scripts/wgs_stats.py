"""
Calculate popgen statistics for a given gene in a VCF file
"""

import numpy as np
import argparse
from Bio import Entrez, SeqIO, codonalign
import gffutils
import allel

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Calculate popgen statistics for a gene in a VCF file")
    parser.add_argument('--vcf', type=str, help='Path to the VCF file')
    parser.add_argument('--annotation', help='Path to genome annotation DB (created by parse_gff.py)', type=str, default=None)
    parser.add_argument('--fasta', type=str, help='Path to the reference genome FASTA file')
    parser.add_argument('--gene', type=str, help='Gene name to analyze')
    parser.add_argument('--chromosome', type=str, help='Chromosome of the gene')
    parser.add_argument('--start', type=int, help='Start position of the gene')
    parser.add_argument('--end', type=int, help='End position of the gene')
    parser.add_argument('--output', type=str, default='stats_output.txt', help='Output file for statistics')
    return parser.parse_args()


def fetch_gene_coordinates(gene, annotation):
    """Fetch gene coordinates from a GFF database."""
    db = gffutils.FeatureDB(annotation, keep_order=True)
    feature = db[gene]
    return feature.chrom, feature.start, feature.end

def fetch_gene_transcripts(gene, annotation, transcriptID='mRNA'):
    """Fetch transcripts for a given gene from a GFF database."""
    db = gffutils.FeatureDB(annotation, keep_order=True)
    transcripts = list(db.children(gene, featuretype=transcriptID))
    return transcripts

def fetch_coding_sequence(transcript, annotation, fasta, cdsID='CDS'):
    """Fetch the coding sequence for a gene from a GFF database and reference FASTA, maintaining genome-wide coordinate references."""
    db = gffutils.FeatureDB(annotation, keep_order=True)
    cds = list(db.children(transcript.id, featuretype=cdsID))
    cdsid = [c.id for c in cds]
    cdsstart = [c.start for c in cds]
    cdsend = [c.end for c in cds]
    cdsseq = [c.sequence(fasta, use_strand=False) for c in cds]
    cdsstrand = transcript.strand
    coding_seq = {
        'id': transcript.id,
        'cds_id': cdsid,
        'cds_start': cdsstart,
        'cds_end': cdsend,
        'cds_seq': cdsseq,
        'strand': cdsstrand
    }
    return coding_seq

def calc_pi(callset):
    """Calculate nucleotide diversity (pi) for a given genotype data."""
    positions = callset['variants/POS']
    # Convert genotype data to allele counts
    ac = allel.GenotypeArray(callset['calldata/GT']).count_alleles()
    # Calculate pi
    pi = allel.sequence_diversity(positions, ac)
    return pi

def calculate_statistics(vcf_file, chromosome, gene_start, gene_end):
    """Calculate popgen statistics for a gene in a VCF file."""
    # Load the VCF file
    callset = allel.read_vcf(vcf_file, region=f'{chromosome}:{gene_start}-{gene_end}')
    # Calculate basic statistics
    stats = {
        'num_samples': callset['calldata/GT'].shape[1],
        'pi': calc_pi(callset)
    }
    return stats

def main():
    """Popgen statistics for a gene..."""
    args = parse_args()
    print(args)
    if not args.start or not args.end:
        # If gene coordinates are not provided, find them
        args.chromosome, args.start, args.end = fetch_gene_coordinates(args.gene, args.annotation)
    print(args.chromosome.strip('AgamP4_'), args.start, args.end)
    stats = calculate_statistics(
        vcf_file=args.vcf,
        chromosome=args.chromosome.strip('AgamP4_'),
        gene_start=args.start,
        gene_end=args.end
    )
    print(stats)

if __name__ == "__main__":
    main()