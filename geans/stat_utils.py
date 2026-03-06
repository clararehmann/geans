import pandas as pd
import numpy as np
import allel
import itertools
from sklearn.cluster import HDBSCAN


def triang2df_fst(mtx1):
    """Convert triangular Fst matrix to a long-format DataFrame."""
    mtx1 = mtx1.where(np.triu(np.ones(mtx1.shape)).astype(bool))
    mtx1 = mtx1.stack().reset_index()
    mtx1.columns = ['loc1', 'loc2', 'fst']
    mtx1 = mtx1.dropna()
    return mtx1


def calc_pi(sequence_array, nsites):
    """Calculate nucleotide diversity (pi) for a given sequence array."""
    if sequence_array is None or len(sequence_array) == 0:
        return 0
    ac = sequence_array.count_alleles()
    an = np.sum(ac, axis=1)
    n_pairs = an * (an - 1) / 2
    n_same = np.sum(ac * (ac - 1) / 2, axis=1)
    n_diff = n_pairs - n_same
    mpd = n_diff / n_pairs
    pi = np.sum(mpd) / nsites if nsites > 0 else 0
    return pi


def calc_theta(sequence_array, nsites):
    """Calculate Watterson's theta (theta hat per base) for a given sequence array."""
    if sequence_array is None or len(sequence_array) == 0:
        return 0
    ac = sequence_array.count_alleles()
    S = ac.count_segregating()
    n = ac.sum(axis=1).max()
    a1 = np.sum(1 / np.arange(1, n))
    theta_hat_w = (S / a1) / nsites if nsites > 0 else 0
    return theta_hat_w


def calc_taj(sequence_array, nsites, min_sites=3):
    """Calculate Tajima's D."""
    if sequence_array is None or len(sequence_array) == 0:
        return 0
    ac = sequence_array.count_alleles()
    S = ac.count_segregating()
    if S < min_sites:
        return np.nan
    n = ac.sum(axis=1).max()
    a1 = np.sum(1 / np.arange(1, n))
    theta_hat_w_abs = S / a1
    an = np.sum(ac, axis=1)
    n_pairs = an * (an - 1) / 2
    n_same = np.sum(ac * (ac - 1) / 2, axis=1)
    n_diff = n_pairs - n_same
    mpd = n_diff / n_pairs
    theta_hat_pi_abs = np.sum(mpd)
    d = theta_hat_pi_abs - theta_hat_w_abs
    a2 = np.sum(1 / (np.arange(1, n)**2))
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n**2 + n + 3) / (9 * n * (n - 1))
    c1 = b1 - (1 / a1)
    c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1**2))
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)
    d_stdev = np.sqrt((e1 * S) + (e2 * S * (S - 1)))
    D = d / d_stdev if d_stdev > 0 else np.nan
    return D


def calc_fst_btwn(gene, seqs, loc1, loc2):
    """Calculate Hudson's Fst between two locations."""
    locs = gene.locs
    assert sum(seqs.count_hom_ref(axis=1) == seqs.shape[1]) == 0
    ind1 = np.where((locs == loc1).all(1))[0]
    ind2 = np.where((locs == loc2).all(1))[0]
    ac1 = seqs[:, ind1].count_alleles()
    ac2 = seqs[:, ind2].count_alleles()
    n, d = allel.hudson_fst(ac1, ac2, fill=np.nan)
    with np.errstate(invalid='ignore'):
        fst = np.nansum(n) / np.nansum(d)
    return fst


def calc_fst_pairwise(gene, seqs):
    """Calculate pairwise Hudson's Fst between all locations."""
    locs = np.unique(gene.locs, axis=0)
    fstmat = pd.DataFrame(index=[tuple(i) for i in locs], columns=[tuple(i) for i in locs])
    for i, j in itertools.combinations(locs, 2):
        itup = tuple(i)
        jtup = tuple(j)
        fst = calc_fst_btwn(gene, seqs, i, j)
        fstmat[itup][jtup] = fst
        fstmat[jtup][itup] = fst
    return fstmat


def calc_fst_wc(gene, seqs, minsamples=2):
    """Calculate Weir and Cockerham's Fst using HDBSCAN-derived population clusters."""
    locs = gene.locs
    cl = HDBSCAN(min_cluster_size=minsamples, metric='haversine').fit_predict(np.radians(locs))
    locs = locs[cl != -1]
    seqs = seqs[:, cl != -1]
    cl = cl[cl != -1]
    c_ = cl[0]
    clinds = []
    clind = [0]
    for i in range(1, len(cl)):
        if cl[i] == c_:
            clind.append(i)
        else:
            clinds.append(clind)
            clind = [i]
            c_ = cl[i]
    clinds.append(clind)
    a, b, c = allel.weir_cockerham_fst(seqs, clinds)
    fst = np.nansum(a) / (np.nansum(a) + np.nansum(b) + np.nansum(c))
    return fst
