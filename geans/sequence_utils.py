import os
import json
import numpy as np
from Bio import SeqIO, codonalign
from Bio.Seq import Seq
import allel

_CODON_TABLE_PATH = os.path.join(os.path.dirname(__file__), 'codontable.json')


def _sort_vcf(callset, samples):
    """Filter and reorder a VCF callset to match a samples metadata DataFrame.

    Samples present in the metadata but absent from the VCF (or vice versa) are
    dropped. Returns the filtered callset, matching sample IDs, and a (n_samples, 2)
    array of geographic coordinates if the metadata contains longitude/latitude or
    x/y columns, otherwise None.
    """
    sampleIDs = np.array(samples['sample_id'].tolist())
    sample_indices = np.array([[-1, np.nan] if len(np.where(callset['samples'] == i)[0]) == 0
                    else [[-1 if len(np.where(sampleIDs == i)[0]) == 0
                    else np.where(sampleIDs == i)[0][0]][0], np.where(callset['samples'] == i)[0][0]]
                    for i in sampleIDs])
    sample_indices = sample_indices[sample_indices[:,0] != -1]
    sample_indices = sample_indices.astype(int)
    sampleIDs = sampleIDs[sample_indices[:,0]]
    callset['calldata/GT'] = allel.GenotypeArray(callset['calldata/GT'][:,sample_indices[:,1],:])
    callset['samples'] = sampleIDs
    samples = samples.set_index('sample_id').reindex(sampleIDs)
    try:
        sampleLOCs = np.array([samples['longitude'], samples['latitude']]).T
    except KeyError:
        try:
            sampleLOCs = np.array([samples['x'], samples['y']]).T
        except KeyError:
            sampleLOCs = None
    return callset, sampleIDs, sampleLOCs


def count_sites(cds_seq, aa_seq):
    """
    Calculate number of silent and non-silent sites in a coding sequence.
    Returns array with shape (n_codons, 2) where:
    - first column is the number of silent sites per codon (divided by 3)
    - second column is the number of non-silent sites per codon (divided by 3)
    Stop codons are excluded from the count.
    """
    codons = np.split(cds_seq, len(aa_seq))
    codon_ss_sites = np.empty((len(codons), 2), dtype=float)
    with open(_CODON_TABLE_PATH, 'r') as f:
        codontab = json.load(f)
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
                    continue
                else:
                    non_silent_sites += 1
        codon_ss_sites[i][0] = silent_sites / 3
        codon_ss_sites[i][1] = non_silent_sites / 3
    return codons, codon_ss_sites


def mean_ss_sites(cds_seq_array, aa_seq):
    """Calculate mean number of silent and non-silent sites across all sequences."""
    ss_sites = np.empty((len(cds_seq_array[0,:]), 2), dtype=float)
    for i in range(len(cds_seq_array[0,:])):
        for j in [0, 1]:
            aa_seq_arr = codonalign.codonseq.CodonSeq("".join(cds_seq_array[:, i, j])).translate()
            c_a, s_a = count_sites(cds_seq_array[:, i, j], aa_seq_arr)
            ss_sites[i, :] = np.sum(s_a, axis=0)
    mean_silent_sites = np.mean(ss_sites[:, 0])
    mean_non_silent_sites = np.mean(ss_sites[:, 1])
    return mean_silent_sites, mean_non_silent_sites


def get_nucleotide_sequences(callset, coding_seq):
    """Get coding sequences from a VCF file and create an array of nucleotide sequences.

    Parameters
    ----------
    callset : dict
        VCF callset from allel.read_vcf.
    coding_seq : dict
        Coding sequence dict with keys 'cds_start', 'cds_end', 'cds_seq', 'strand'.
        Produced by Gene.fetch_coding_sequence or Transcript.coding_seq.

    Returns
    -------
    cds_seq : np.ndarray
        Reference coding sequence as a 1-D character array.
    cds_seq_array : np.ndarray
        Per-sample coding sequences with shape (cds_length, n_samples, 2).
    """
    positions = callset['variants/POS']
    positionmask = np.repeat(False, len(positions))
    for start, end in zip(coding_seq['cds_start'], coding_seq['cds_end']):
        positionmask[(positions >= start) & (positions <= end)] = True
    positions = positions[positionmask]
    alt = callset['variants/ALT'][positionmask]
    gt = allel.GenotypeArray(callset['calldata/GT'][positionmask])

    if coding_seq['strand'] == '-':
        cds_seq = [str(Seq(x).reverse_complement()) for x in coding_seq['cds_seq']]
        cds_seq = np.fromiter("".join(cds_seq), dtype='U1')
        cds_positions = np.concatenate([np.flip(np.arange(s, e+1))
                                        for s, e in zip(coding_seq['cds_start'], coding_seq['cds_end'])])
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    else:
        cds_seq = np.fromiter("".join(coding_seq['cds_seq']), dtype='U1')
        cds_positions = np.concatenate([np.arange(s, e+1)
                                        for s, e in zip(coding_seq['cds_start'], coding_seq['cds_end'])])
    assert np.all(cds_seq[0:3] == ['A', 'T', 'G']), "Coding sequence does not start with ATG."
    cds_seq_array = np.tile(cds_seq, (2, gt.shape[1], 1)).T
    variant_cds_positions, variant_cds_index, variant_position_index = np.intersect1d(
        cds_positions, positions, return_indices=True)
    for vcp, vci, vpi in zip(variant_cds_positions, variant_cds_index, variant_position_index):
        for i in np.unique(gt[vpi]):
            if i > 0:
                if coding_seq['strand'] == '-':
                    altNT = complement_dict[alt[vpi][i-1]]
                else:
                    altNT = alt[vpi][i-1]
                assert altNT != cds_seq[vci], "Alternate allele matches reference."
                cds_seq_array[vci][gt[vpi] == i] = altNT
    return cds_seq, cds_seq_array


def get_aa_sequences(cds_seq_array, cds_seq, ID=None, fasta_check=None):
    """Translate array of coding sequences into protein sequences.

    Returns
    -------
    wt_AA_seq : np.ndarray
        Wild-type amino acid sequence.
    aa_array : np.ndarray
        Per-sample amino acid sequences.
    """
    wt_AA_seq = np.fromiter(codonalign.codonseq.CodonSeq("".join(cds_seq)).translate(), dtype='U1')
    aa_array = np.apply_along_axis(
        lambda x: np.fromiter(codonalign.codonseq.CodonSeq("".join(x)).translate(), dtype='U1'),
        0, cds_seq_array)
    if fasta_check:
        pdct = SeqIO.index(fasta_check, 'fasta')
        assert np.all(np.fromiter(pdct[ID.replace('R', 'P')].seq, dtype='U1') == wt_AA_seq[:-1])
    return wt_AA_seq, aa_array


def identify_ss_ns_sites(cds_seq_array, cds_seq, aa_array, wt_AA_seq):
    """Identify synonymous and non-synonymous mutations in the array of coding sequences."""
    codon_array = np.split(cds_seq, len(cds_seq) // 3)
    synonymous_array = np.zeros(cds_seq_array.shape, dtype='U10')
    nonsynonymous_array = np.zeros(cds_seq_array.shape, dtype='U10')
    for i in range(len(aa_array[0, :, 0])):
        for j in [0, 1]:
            codons = np.split(cds_seq_array[:, i, j], len(cds_seq) // 3)
            mutcodons = np.where([codons[x] != codon_array[x] for x in range(len(codon_array))])
            mutloci = mutcodons[0] * 3 + mutcodons[1]
            nscodons = np.where(aa_array[:, i, j] != wt_AA_seq)[0]
            nsloci = mutloci[np.isin(mutcodons[0], nscodons)]
            sscodons = np.array([x for x in mutcodons[0] if x not in nscodons])
            ssloci = mutloci[np.isin(mutcodons[0], sscodons)]
            if len(ssloci) > 0:
                synonymous_array[ssloci, i, j] = cds_seq_array[ssloci, i, j]
            if len(nsloci) > 0:
                nonsynonymous_array[nsloci, i, j] = cds_seq_array[nsloci, i, j]
    return synonymous_array, nonsynonymous_array


def convert_to_binaryarray(sequence_array, wt_sequence):
    """Convert nt or amino acid sequences to binary array format for genotype calls."""
    for i in range(sequence_array.shape[0]):
        wt_val = wt_sequence[i]
        all_val = np.unique(sequence_array[i, :, :])
        all_val = np.roll(all_val, -1 * np.where(all_val == wt_val)[0])
        for a in range(len(all_val)):
            val_id = all_val[a]
            msk = sequence_array[i, :, :] == val_id
            sequence_array[i, msk] = a
    sequence_array = sequence_array.astype(int)
    sequence_array = allel.GenotypeArray(sequence_array)
    nonvariant_mask = sequence_array.count_hom_ref(axis=1) < sequence_array.shape[1]
    sequence_array = sequence_array[nonvariant_mask]
    if len(sequence_array) == 0:
        return None
    else:
        return sequence_array
