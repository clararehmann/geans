"""
Calculate popgen statistics for a given gene in a VCF file
"""

from matplotlib.pyplot import fill
import numpy as np
import argparse
from Bio import Entrez, SeqIO, codonalign
from Bio.Seq import Seq
import gffutils
import allel
from allel.util import ignore_invalid

codontab = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G'     # Glicina
}

class Gene:
    def __init__(self, name, chromosome, start, end):
        self.name = name
        self.chromosome = chromosome
        self.start = start
        self.end = end

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

def count_sites(cds_seq, aa_seq):
    """ 
    calculate number of silent and non-silent sites in a coding sequence 
    returns array with shape (n_codons, 2) where:
    - first column is the number of silent sites per codon (divided by 3)
    - second column is the number of non-silent sites per codon (divided by 3)
    stop codons are excluded from the count
    """
    codons = np.split(cds_seq, len(aa_seq))
    codon_ss_sites = np.empty((len(codons), 2), dtype=float)

    for i, codon in enumerate(codons):
        silent_sites = 0
        non_silent_sites = 0
        originalAA = codontab[''.join(codon)]
        assert originalAA == aa_seq[i], f"Original AA {originalAA} does not match expected {aa_seq[i]} at position {i}."
        states = ["A", "C", "G", "T"]
        for codonpos in range(3):
            for state in states:
                new_codon = codon.copy()
                new_codon[codonpos] = state
                if ''.join(new_codon) == ''.join(codon):
                    continue
                newAA = codontab[''.join(new_codon)]
                if newAA == originalAA:
                    silent_sites += 1
                elif newAA == '*':
                    # stop codon, we don't count these as silent or non-silent
                    continue
                else:
                    non_silent_sites += 1
        codon_ss_sites[i][0] = silent_sites/3
        codon_ss_sites[i][1] = non_silent_sites/3
    return codons, codon_ss_sites

def mean_ss_sites(cds_seq_array, aa_seq):
    """

    """
    ss_sites = np.empty((len(cds_seq_array[0,:]), 2), dtype=float)
    print(ss_sites.shape)
    for i in range(len(cds_seq_array[0,:])):
        for j in [0, 1]:
            #print(cds_seq_array[:, i, j])
            aa_seq_arr = codonalign.codonseq.CodonSeq("".join(cds_seq_array[:, i, j])).translate()
            c_a, s_a, = count_sites(cds_seq_array[:, i, j], aa_seq_arr)
            #iprint(s_a)
            #print(np.mean(s_a, axis=0))
            ss_sites[i, :] = np.sum(s_a, axis=0)
    #print(ss_sites)
    mean_silent_sites = np.mean(ss_sites[:, 0])
    mean_non_silent_sites = np.mean(ss_sites[:, 1])
    return mean_silent_sites, mean_non_silent_sites

def calc_pi(callset, start=None, stop=None, coding_seq=None, translate=False):
    """Calculate nucleotide diversity (pi) for a given genotype data."""
    positions = callset['variants/POS']

    if coding_seq:
        # If coding sequence is provided, filter positions to those within the coding sequence
        positionmask = np.repeat(False, len(positions))
        for start, end in zip(coding_seq['cds_start'], coding_seq['cds_end']):
            positionmask[(positions >= start) & (positions <= end)] = True
        positions = positions[positionmask]
        
        if not translate:
            # Convert genotype data to allele counts
            ac = allel.GenotypeArray(callset['calldata/GT'][positionmask]).count_alleles()
            # Calculate pi by hand to deal with # of possible sites glitch

            # mean pairwise difference
            an = np.sum(ac, axis=1)
            n_pairs = an * (an - 1) / 2
            n_same = np.sum(ac * (ac - 1) / 2, axis=1)
            n_diff = n_pairs - n_same
            mpd = n_diff / n_pairs
            mpd_sum = np.sum(mpd)
            n_bases = np.sum([x+1 - y for x, y in zip(coding_seq['cds_end'], coding_seq['cds_start'])])
            pi = mpd_sum / n_bases if n_bases > 0 else 0

    if translate:
        # Calculate pi for translated sequences

        # deal with data...
        alt = callset['variants/ALT'][positionmask] # alternate alleles for variant positions
        gt = allel.GenotypeArray(callset['calldata/GT'][positionmask]) # genotype calls for variant positions
        print(coding_seq['cds_id'], coding_seq['strand'])
        #cds_index = np.arange(len(cds_positions)) # indexer for CDS positions
        if coding_seq['strand'] == '-':
            cds_seq = [str(Seq(x).reverse_complement()) for x in coding_seq['cds_seq']] # reverse complement the coding sequence if on the reverse strand
            cds_seq = np.fromiter("".join(cds_seq), dtype='U1') # numpy array of coding sequence
            cds_positions = np.concatenate([np.flip(np.arange(s, e+1)) for s, e in zip(coding_seq['cds_start'], coding_seq['cds_end'])])
            complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        else:
            cds_seq = np.fromiter("".join(coding_seq['cds_seq']), dtype='U1') # numpy array of coding sequence
            cds_positions = np.concatenate([np.arange(s, e+1) for s, e in zip(coding_seq['cds_start'], coding_seq['cds_end'])])
        assert(np.all(cds_seq[0:3] == ['A','T','G'])) # ensure start codon

        cds_seq_array = np.tile(cds_seq, (2, gt.shape[1], 1)).T
        variant_cds_positions, variant_cds_index, variant_position_index = np.intersect1d(cds_positions, positions, return_indices=True) # find intersection of CDS and variant positions
        for vcp, vci, vpi in zip(variant_cds_positions, variant_cds_index, variant_position_index):
            for i in np.unique(gt[vpi]):
                if i > 0:
                    if coding_seq['strand'] == '-':
                        altNT = complement_dict[alt[vpi][i-1]] 
                    else:
                        altNT = alt[vpi][i-1]
                    assert(altNT != cds_seq[vci]) # ensure alternate allele is not the same as the reference
                    cds_seq_array[vci][gt[vpi]==i] = altNT

        wt_AA_seq = np.fromiter(codonalign.codonseq.CodonSeq("".join(cds_seq)).translate(), dtype='U1')  # wild-type coding sequence
        aa_array = np.apply_along_axis(lambda x: np.fromiter(codonalign.codonseq.CodonSeq("".join(x)).translate(), dtype='U1'), 0, cds_seq_array) # array of AA sequences (translated)
        pdct = SeqIO.index('VectorBase-68_AgambiaePEST_AnnotatedProteins.fasta', 'fasta')
        assert(np.all(np.fromiter(pdct[coding_seq['id'].replace('R','P')].seq, dtype='U1') == wt_AA_seq[:-1]))
        nmean_ss, mean_ns = mean_ss_sites(cds_seq_array, wt_AA_seq)

        """
        with open('aa_array.txt', 'w') as f:
            for i in range(aa_array.shape[0]):
                f.write(' '.join(aa_array[i, :, :].flatten()) + '\n')
        with open('cds_seq_array.txt', 'w') as f:
            for i in range(cds_seq_array.shape[0]):
                f.write(' '.join(cds_seq_array[i, :, :].flatten()) + '\n')
        with open('genotypes.txt', 'w') as f:
            for i in range(gt.shape[0]):
                f.write(' '.join(map(str, gt[i, :])) + '\n')
        np.savetxt('wt_AA_seq.txt', wt_AA_seq, fmt='%s')
        np.savetxt('wt_cds_seq.txt', cds_seq, fmt='%s')
        np.savetxt('positions.txt', positions)
        np.savetxt('alt.txt', alt, fmt='%s')
        """
        # convert aa_array to 0/1/2/3... for GT array format
        for i in range(aa_array.shape[0]):
            wt_aa = wt_AA_seq[i]
            all_aa = np.unique(aa_array[i, :, :])
            all_aa = np.roll(all_aa, -1*np.where(all_aa == wt_aa)[0])
            for a in range(len(all_aa)):
                AA_id = all_aa[a]
                msk = aa_array[i, :, :] == AA_id
                aa_array[i, msk] = a  # replace amino acid with its index
        aa_array = aa_array.astype(int)  # convert to integer type
        aa_array = allel.GenotypeArray(aa_array)
        position_array = np.arange(len(aa_array))
        nonvariant_mask = aa_array.count_hom_ref(axis=1) < aa_array.shape[1]
        aa_array = aa_array[nonvariant_mask]
        position_array = position_array[nonvariant_mask]

        # manually calculate pi 
        an = np.sum(aa_array.count_alleles(), axis=1)
        ac = aa_array.count_alleles()
        n_pairs = an * (an - 1) / 2
        n_same = np.sum(ac * (ac - 1) / 2, axis=1) 
        n_diff = n_pairs - n_same
        mpd = n_diff / n_pairs
        mpd_sum = np.sum(mpd)
        n_pos = mean_ns
        pi = mpd_sum / n_pos if n_pos > 0 else 0
    else:
        # Convert genotype data to allele counts
        ac = allel.GenotypeArray(callset['calldata/GT']).count_alleles()
        # Calculate pi
        pi = allel.sequence_diversity(positions, ac, start=start, stop=stop)
    return pi

def calculate_statistics(vcf_file, gene, annotation, fasta, chromosome, gene_start, gene_end):
    """Calculate popgen statistics for a gene in a VCF file."""
    # Load the VCF file
    callset = allel.read_vcf(vcf_file, region=f'{chromosome}:{gene_start}-{gene_end}')
    transcripts = fetch_gene_transcripts(gene, annotation=annotation)
    # Calculate basic statistics
    stats = {
        'num_samples': callset['calldata/GT'].shape[1],
        'pi': calc_pi(callset, start=gene_start, stop=gene_end),
        'coding sequence pi': calc_pi(callset, coding_seq=fetch_coding_sequence(transcripts[0], annotation, fasta)),
        'amino acid pi': calc_pi(callset, coding_seq=fetch_coding_sequence(transcripts[0], annotation, fasta), translate=True)
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
        vcf_file=args.vcf.replace('CHROMOSOME', args.chromosome.split('_')[-1]),
        gene=args.gene,
        annotation=args.annotation,
        fasta=args.fasta,
        chromosome=args.chromosome.split('_')[-1],
        gene_start=args.start,
        gene_end=args.end
    )
    print(stats)


if __name__ == "__main__":
    main()