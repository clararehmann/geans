from ctypes import alignment
from wgs_stats import Gene
import numpy as np
from Bio import SeqIO, AlignIO, Align
from Bio.Align import CodonAligner
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from heapq import heapify, heappop, heappush
from math import floor
from scipy.stats import chi2_contingency
import argparse, subprocess, sys



def load_gene(gene_name, annotation, fasta, vcf, filter):
    """Load gene data including coordinates and variations."""
    gene = Gene(name=gene_name)
    gene.fetch_gene_coordinates(annotation, fasta)
    gene.fetch_variation(vcf_file=vcf, annotation=annotation, fasta=fasta, filter=filter)
    return gene

def load_records(species1, species2, gene1, gene2, out_prefix):
    """Generate SeqRecord-parseable objects for both genes."""
    species_list = []
    nt_records = []
    aa_records = []

    for gene, species in zip([gene1, gene2], [species1, species2]):
        for i in range(gene.transcripts[list(gene.transcripts.keys())[0]]['sequences']['nt_seq_array'].shape[1]):
            species_list.append(species)
            species_list.append(species)
            nt_seq = gene.transcripts[list(gene.transcripts.keys())[0]]['sequences']['nt_seq_array'][:,i]
            nt_records.append(SeqRecord(Seq(''.join(nt_seq[:,0])), id=f"{gene.name}_{species}_nt_{i}.0", description=""))
            nt_records.append(SeqRecord(Seq(''.join(nt_seq[:,1])), id=f"{gene.name}_{species}_nt_{i}.1", description=""))

        for i in range(gene.transcripts[list(gene.transcripts.keys())[0]]['sequences']['aa_seq_array'].shape[1]):
            aa_seq = gene.transcripts[list(gene.transcripts.keys())[0]]['sequences']['aa_seq_array'][:,i]
            aa_records.append(SeqRecord(Seq(''.join(aa_seq[:,0])), id=f"{gene.name}_{species}_aa_{i}.0", description=""))
            aa_records.append(SeqRecord(Seq(''.join(aa_seq[:,1])), id=f"{gene.name}_{species}_aa_{i}.1", description=""))
    SeqIO.write(nt_records, out_prefix + f"{gene1.name}-{gene2.name}-temp_nt.fasta", "fasta")
    SeqIO.write(aa_records, out_prefix + f"{gene1.name}-{gene2.name}-temp_aa.fasta", "fasta")
    return species_list, nt_records, aa_records

def align_codons(gene1, gene2, out_prefix):
    """Align codon sequences using amino acid alignment as a guide."""
    cmd = f"muscle -super5 {out_prefix}{gene1}-{gene2}-temp_aa.fasta -output {out_prefix}{gene1}-{gene2}-temp_aligned_aa.fasta"
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, text=True)
    nt_recs = SeqIO.index(out_prefix + f"{gene1}-{gene2}-temp_nt.fasta", "fasta")
    aligner = CodonAligner()
    protein_aln = Align.read(out_prefix + f"{gene1}-{gene2}-temp_aligned_aa.fasta", "fasta")
    codon_alignments = []
    for prec in protein_aln.sequences:
        ntrec = nt_recs[prec.id.replace('_aa_', '_nt_')]
        alignments = aligner.align(prec, ntrec)
        assert len(alignments) == 1
        codon_alignment = next(alignments)
        codon_alignments.append(codon_alignment)
    alignment = protein_aln.mapall(codon_alignments)
    return alignment

########################################
# this section is a bunch of stuff stolen from biopython
# ## https://github.com/biopython/biopython/blob/master/Bio/Align/analysis.py
#########################################

def _get_codon2codon_matrix(codon_table):
    """Get codon codon substitution matrix (PRIVATE).
    Elements in the matrix are number of synonymous and nonsynonymous
    substitutions required for the substitution.
    """
    bases = ("A", "T", "C", "G")
    codons = [
        codon
        for codon in list(codon_table.forward_table.keys()) + codon_table.stop_codons
        if "U" not in codon
    ]
    # set up codon_dict considering stop codons
    codon_dict = codon_table.forward_table.copy()
    for stop in codon_table.stop_codons:
        codon_dict[stop] = "stop"
    # count site
    num = len(codons)
    G = {}  # graph for substitution
    nonsyn_G = {}  # graph for nonsynonymous substitution
    graph = {}
    graph_nonsyn = {}
    for i, codon in enumerate(codons):
        graph[codon] = {}
        graph_nonsyn[codon] = {}
        for p in range(3):
            for base in bases:
                tmp_codon = codon[0:p] + base + codon[p + 1 :]
                if codon_dict[codon] != codon_dict[tmp_codon]:
                    graph_nonsyn[codon][tmp_codon] = 1
                    graph[codon][tmp_codon] = 1
                else:
                    if codon != tmp_codon:
                        graph_nonsyn[codon][tmp_codon] = 0.1
                        graph[codon][tmp_codon] = 1
    for codon1 in codons:
        nonsyn_G[codon1] = {}
        G[codon1] = {}
        for codon2 in codons:
            if codon1 == codon2:
                nonsyn_G[codon1][codon2] = 0
                G[codon1][codon2] = 0
            else:
                nonsyn_G[codon1][codon2] = _dijkstra(graph_nonsyn, codon1, codon2)
                G[codon1][codon2] = _dijkstra(graph, codon1, codon2)
    return G, nonsyn_G

def _dijkstra(graph, start, end):
    """Dijkstra's algorithm Python implementation (PRIVATE).

    Algorithm adapted from
    http://thomas.pelletier.im/2010/02/dijkstras-algorithm-python-implementation/.
    However, an obvious bug in::

        if D[child_node] >(<) D[node] + child_value:

    is fixed.
    This function will return the distance between start and end.

    Arguments:
     - graph: Dictionary of dictionary (keys are vertices).
     - start: Start vertex.
     - end: End vertex.

    Output:
       List of vertices from the beginning to the end.

    """
    D = {}  # Final distances dict
    P = {}  # Predecessor dict
    # Fill the dicts with default values
    for node in graph.keys():
        D[node] = 100  # Vertices are unreachable
        P[node] = ""  # Vertices have no predecessors
    D[start] = 0  # The start vertex needs no move
    unseen_nodes = list(graph.keys())  # All nodes are unseen
    while len(unseen_nodes) > 0:
        # Select the node with the lowest value in D (final distance)
        shortest = None
        node = ""
        for temp_node in unseen_nodes:
            if shortest is None:
                shortest = D[temp_node]
                node = temp_node
            elif D[temp_node] < shortest:
                shortest = D[temp_node]
                node = temp_node
        # Remove the selected node from unseen_nodes
        unseen_nodes.remove(node)
        # For each child (ie: connected vertex) of the current node
        for child_node, child_value in graph[node].items():
            if D[child_node] > D[node] + child_value:
                D[child_node] = D[node] + child_value
                # To go to child_node, you have to go through node
                P[child_node] = node
        if node == end:
            break
    # Set a clean path
    path = []
    # We begin from the end
    node = end
    distance = 0
    # While we are not arrived at the beginning
    while not (node == start):
        if path.count(node) == 0:
            path.insert(0, node)  # Insert the predecessor of the current node
            node = P[node]  # The current node becomes its predecessor
        else:
            break
    path.insert(0, start)  # Finally, insert the start vertex
    for i in range(len(path) - 1):
        distance += graph[path[i]][path[i + 1]]
    return distance

def _count_replacement(codons, G):
    """Count replacement needed for a given codon_set (PRIVATE)."""
    if len(codons) == 1:
        return 0, 0
    elif len(codons) == 2:
        codons = list(codons)
        return floor(G[codons[0]][codons[1]])
    else:
        subgraph = {
            codon1: {codon2: G[codon1][codon2] for codon2 in codons if codon1 != codon2}
            for codon1 in codons
        }
        return _prim(subgraph)


def _prim(G):
    """Prim's algorithm to find minimum spanning tree (PRIVATE).

    Code is adapted from
    http://programmingpraxis.com/2010/04/09/minimum-spanning-tree-prims-algorithm/
    """
    nodes = []
    edges = []
    for i in G.keys():
        nodes.append(i)
        for j in G[i]:
            if (i, j, G[i][j]) not in edges and (j, i, G[i][j]) not in edges:
                edges.append((i, j, G[i][j]))
    conn = defaultdict(list)
    for n1, n2, c in edges:
        conn[n1].append((c, n1, n2))
        conn[n2].append((c, n2, n1))
    mst = []  # minimum spanning tree
    used = set(nodes[0])
    usable_edges = conn[nodes[0]][:]
    heapify(usable_edges)
    while usable_edges:
        cost, n1, n2 = heappop(usable_edges)
        if n2 not in used:
            used.add(n2)
            mst.append((n1, n2, cost))
            for e in conn[n2]:
                if e[2] not in used:
                    heappush(usable_edges, e)
    length = 0
    for p in mst:
        length += floor(p[2])
    return length

##########################################

def calculate_mk_stats(alignment, species_list):
    # based on biopython code in Bio.Align.analysis
    codon_table = CodonTable.generic_by_id[1]
    G, nonsyn_G = _get_codon2codon_matrix(codon_table=codon_table)  
    unique_species = list(set(species_list))
    sequences = []
    for sequence in alignment.sequences:
        sequences.append(str(sequence.seq))
    syn_fix, nonsyn_fix, syn_poly, nonsyn_poly = 0, 0, 0, 0
    starts = sys.maxsize
    for ends in alignment.coordinates.transpose():
        step = min(ends - starts)
        for j in range(0, step, 3):
            codons = {key: [] for key in unique_species}
            for key, sequence, start in zip(species_list, sequences, starts):
                codon = sequence[start + j : start + j + 3]
                codons[key].append(codon)
            fixed = True
            all_codons = set()
            for value in codons.values():
                value = set(value)
                if len(value) > 1:
                    fixed = False
                all_codons.update(value)
            if len(all_codons) == 1:
                continue
            nonsyn = _count_replacement(all_codons, nonsyn_G)
            syn = _count_replacement(all_codons, G) - nonsyn
            if fixed is True:
                # fixed
                nonsyn_fix += nonsyn
                syn_fix += syn
            else:
                # not fixed
                nonsyn_poly += nonsyn
                syn_poly += syn    
        starts = ends
    alpha = 1 - ((syn_fix * nonsyn_poly) / (nonsyn_fix * syn_poly))
    tab = np.array([[syn_fix, syn_poly], [nonsyn_fix, nonsyn_poly]])
    chi2, p, dof, ex = chi2_contingency(tab)
    return syn_fix, nonsyn_fix, syn_poly, nonsyn_poly, alpha, p

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Calculate MK test statistics between two orthologous genes from two different species.")
    parser.add_argument('--gene1', type=str, help='Gene name to analyze', required=True)
    parser.add_argument('--gene2', type=str, help='Orthologous gene name to analyze', required=True)    
    parser.add_argument('--species1', type=str, help='Species 1 (e.g. gamb)', required=True)
    parser.add_argument('--species2', type=str, help='Species 2 (e.g. afun)', required=True)    
    parser.add_argument('--vcf1', type=str, help='VCF file for species 1', required=True)
    parser.add_argument('--vcf2', type=str, help='VCF file for species 2', required=True)
    parser.add_argument('--filter1', type=str, help='Site filter VCF file for species 1', required=True)
    parser.add_argument('--filter2', type=str, help='Site filter VCF file for species 2', required=True)    
    parser.add_argument('--fasta1', type=str, help='Reference FASTA file for species 1', required=True)
    parser.add_argument('--fasta2', type=str, help='Reference FASTA file for species 2', required=True)
    parser.add_argument('--annotation1', type=str, help='Annotation DB file for species 1', required=True)
    parser.add_argument('--annotation2', type=str, help='Annotation DB file for species 2', required=True)
    parser.add_argument('--outprefix', type=str, help='Output prefix for result files', required=True)
    return parser.parse_args()

def main():
    args = parse_args()
    gene1 = load_gene(args.gene1, args.annotation1, args.fasta1, args.vcf1, args.filter1)
    gene2 = load_gene(args.gene2, args.annotation2, args.fasta2, args.vcf2, args.filter2)
    species_list, nt_records, aa_records = load_records(args.species1, args.species2, gene1, gene2, out_prefix=args.outprefix)
    alignment = align_codons(args.gene1, args.gene2, args.outprefix)
    syn_fix, nonsyn_fix, syn_poly, nonsyn_poly, alpha, p = calculate_mk_stats(alignment, species_list)
    print(f"Gene1: {args.gene1}, Gene2: {args.gene2}")
    print(f"Synonymous fixed: {syn_fix}, Nonsynonymous fixed: {nonsyn_fix}")
    print(f"Synonymous polymorphic: {syn_poly}, Nonsynonymous polymorphic: {nonsyn_poly}")
    print(f"Alpha: {alpha}")
    print(f"Fisher's exact test p-value: {p}")

if __name__ == "__main__":
    main()