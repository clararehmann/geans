import argparse
import sys
from .gene_stats import Gene

_VALID_STATS = ['pi', 'theta', 'tajD', 'Fst', 'pairwiseFst']
_DEFAULT_STATS = ['pi', 'theta', 'tajD', 'Fst']


def parse_args():
    parser = argparse.ArgumentParser(
        prog='geans',
        description='Calculate per-gene population genetics statistics from a VCF file.',
    )
    parser.add_argument('--id', required=True, metavar='GENE_ID',
                        help='Gene identifier (required; must match annotation database)')
    parser.add_argument('--vcf', required=True, metavar='VCF',
                        help='Path to VCF file (required)')
    parser.add_argument('--annotation', required=True, metavar='DB',
                        help='Path to gffutils annotation database (required)')
    parser.add_argument('--fasta', required=True, metavar='FASTA',
                        help='Path to reference genome FASTA (required)')
    parser.add_argument('--metadata', default=None, metavar='TSV',
                        help='Path to sample metadata TSV (must contain sample_id column)')
    parser.add_argument('--filter', default=None, metavar='VCF',
                        dest='filter_vcf',
                        help='Path to a filter VCF; only sites present in this file are kept')
    parser.add_argument('--protein-fasta', default=None, metavar='FASTA',
                        dest='protein_fasta',
                        help='Path to protein FASTA for amino acid sequence validation')
    parser.add_argument('--longest', default=False, action='store_true',
                        help='Only analyze the longest transcript for the gene (default: False)')
    parser.add_argument('--stats', nargs='+', default=_DEFAULT_STATS,
                        choices=_VALID_STATS, metavar='STAT',
                        help=(
                            f'Statistics to calculate (default: {" ".join(_DEFAULT_STATS)}). '
                            f'Choices: {", ".join(_VALID_STATS)}'
                        ))
    parser.add_argument('--output', default=None, metavar='FILE',
                        help='Write statistics to this tab-separated file')
    return parser.parse_args()


def main():
    args = parse_args()
    stat_flags = {s: True for s in args.stats}

    gene = Gene(
        args.id,
        vcf=args.vcf,
        annotation=args.annotation,
        fasta=args.fasta,
        metadata=args.metadata,
        statistics=stat_flags,
    )

    gene.fetch_gene_coordinates()
    gene.fetch_gene_transcripts()

    if args.longest:
        gene.keep_longest_transcript()

    gene.fetch_variation(filter=args.filter_vcf, protein_fasta=args.protein_fasta)
    gene.calculate_statistics()

    print(gene)
    print()
    for tx in gene:
        print(tx)
        print()

    if args.output:
        first = True
        for tx in gene:
            gene.save_to_df(tx.id, output_file=args.output, append=not first)
            first = False
        print(f'Statistics written to {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
