"""
Integration tests for the Gene class and its full load-and-calculate workflow.

Requires fixture files in testing/:
    test.vcf.gz, filter.vcf.gz, agam.db, agam.fasta, metadata.txt
"""
import pathlib
import numpy as np
import pandas as pd
import allel
import pytest
from sklearn.cluster import HDBSCAN
from geans import Gene, Transcript, Variants
from tests.conftest import GENE_ID, TRANSCRIPT, VCF, ANNOTATION, FASTA, METADATA


# ---------------------------------------------------------------------------
# Gene initialisation
# ---------------------------------------------------------------------------

class TestGeneInit:
    def test_default_statistics_flags(self):
        g = Gene(GENE_ID, vcf=VCF, annotation=ANNOTATION, fasta=FASTA, metadata=METADATA)
        assert isinstance(g.statistics, dict)
        assert g.statistics['pi'] is True
        assert g.statistics['theta'] is True
        assert g.statistics['tajD'] is True
        assert g.statistics['Fst'] is True
        assert g.statistics['pairwiseFst'] is False

    def test_custom_statistics_flags(self):
        g = Gene(GENE_ID, vcf=VCF, annotation=ANNOTATION, fasta=FASTA, metadata=METADATA,
                 statistics={'pi': True, 'theta': False})
        assert g.statistics['pi'] is True
        assert g.statistics['theta'] is False

    def test_metadata_loaded_as_dataframe(self):
        g = Gene(GENE_ID, vcf=VCF, annotation=ANNOTATION, fasta=FASTA, metadata=METADATA)
        assert isinstance(g.metadata, pd.DataFrame)
        assert 'sample_id' in g.metadata.columns

    def test_no_metadata(self):
        g = Gene(GENE_ID, vcf=VCF, annotation=ANNOTATION, fasta=FASTA)
        assert g.metadata is None

    def test_initial_state(self):
        g = Gene(GENE_ID, vcf=VCF, annotation=ANNOTATION, fasta=FASTA)
        assert g.name == GENE_ID
        assert g.transcripts == {}
        assert g.chromosome is None
        assert g.samples is None
        assert g.locs is None


# ---------------------------------------------------------------------------
# fetch_gene_coordinates
# ---------------------------------------------------------------------------

class TestFetchGeneCoordinates:
    def test_chromosome_contains_2L(self, gene_loaded):
        assert gene_loaded.chromosome is not None
        assert '2L' in gene_loaded.chromosome

    def test_start_end_set(self, gene_loaded):
        assert gene_loaded.start is not None
        assert gene_loaded.end is not None
        assert gene_loaded.start > 0
        assert gene_loaded.end > gene_loaded.start

    def test_length_computed(self, gene_loaded):
        assert gene_loaded.length == gene_loaded.end - gene_loaded.start + 1

    def test_gc_content_valid(self, gene_loaded):
        assert gene_loaded.gc_content is not None
        assert 0.0 <= gene_loaded.gc_content <= 1.0

    def test_returns_chrom_start_end_tuple(self):
        g = Gene(GENE_ID, vcf=VCF, annotation=ANNOTATION, fasta=FASTA)
        chrom, start, end = g.fetch_gene_coordinates()
        assert '2L' in chrom
        assert start < end


# ---------------------------------------------------------------------------
# fetch_gene_transcripts + Gene accessor methods
# ---------------------------------------------------------------------------

class TestFetchTranscripts:
    def test_transcripts_populated(self, gene_loaded):
        assert len(gene_loaded.transcripts) >= 1

    def test_known_transcript_present(self, gene_loaded):
        assert TRANSCRIPT in gene_loaded.transcripts

    def test_transcript_type(self, gene_loaded):
        tx = gene_loaded.transcripts[TRANSCRIPT]
        assert isinstance(tx, Transcript)

    def test_gene_len(self, gene_loaded):
        assert len(gene_loaded) == len(gene_loaded.transcripts)

    def test_gene_iter_yields_transcripts(self, gene_loaded):
        items = list(gene_loaded)
        assert len(items) == len(gene_loaded.transcripts)
        assert all(isinstance(tx, Transcript) for tx in items)

    def test_gene_getitem(self, gene_loaded):
        tx = gene_loaded[TRANSCRIPT]
        assert isinstance(tx, Transcript)
        assert tx.id == TRANSCRIPT

    def test_gene_getitem_missing_key_raises(self, gene_loaded):
        with pytest.raises(KeyError):
            _ = gene_loaded['NONEXISTENT_TRANSCRIPT']


# ---------------------------------------------------------------------------
# Transcript attributes (derived from annotation + FASTA only)
# ---------------------------------------------------------------------------

class TestTranscriptAttributes:
    def test_cds_length_positive_and_in_frame(self, gene_loaded):
        tx = gene_loaded[TRANSCRIPT]
        assert tx.cds_length > 0
        assert tx.cds_length % 3 == 0

    def test_aa_seq_matches_cds_length(self, gene_loaded):
        tx = gene_loaded[TRANSCRIPT]
        # translate() includes stop codon, so len == cds_length // 3
        assert len(tx.AAseq) == tx.cds_length // 3

    def test_nt_seq_starts_with_atg(self, gene_loaded):
        tx = gene_loaded[TRANSCRIPT]
        assert ''.join(tx.NTseq[:3]) == 'ATG'

    def test_coding_seq_dict_keys(self, gene_loaded):
        tx = gene_loaded[TRANSCRIPT]
        for key in ('id', 'cds_id', 'cds_start', 'cds_end', 'cds_seq', 'strand'):
            assert key in tx.coding_seq

    def test_transcript_id_matches(self, gene_loaded):
        tx = gene_loaded[TRANSCRIPT]
        assert tx.id == TRANSCRIPT
        assert tx.gene == GENE_ID


# ---------------------------------------------------------------------------
# fetch_variation / Variants structure
# ---------------------------------------------------------------------------

class TestFetchVariation:
    def test_variants_loaded(self, gene_with_variation):
        tx = gene_with_variation[TRANSCRIPT]
        assert tx.variants is not None
        assert isinstance(tx.variants, Variants)

    def test_samples_set(self, gene_with_variation):
        assert gene_with_variation.samples is not None
        assert len(gene_with_variation.samples) > 0
        assert len(gene_with_variation.samples) == len(set(gene_with_variation.samples))

    def test_locs_shape(self, gene_with_variation):
        locs = gene_with_variation.locs
        assert locs is not None
        assert locs.ndim == 2
        assert locs.shape[1] == 2  # (longitude, latitude)
        assert locs.shape[0] == len(gene_with_variation.samples)

    def test_positions_sorted(self, gene_with_variation):
        pos = gene_with_variation[TRANSCRIPT].variants.positions
        assert len(pos) > 0
        assert np.all(np.diff(pos) >= 0)

    def test_gt_array_shape(self, gene_with_variation):
        tx = gene_with_variation[TRANSCRIPT]
        gt = tx.variants.gt_array
        n_samples = len(gene_with_variation.samples)
        assert gt.ndim == 3
        assert gt.shape[1] == n_samples
        assert gt.shape[2] == 2  # diploid genotype

    def test_wt_nt_seq_matches_cds_length(self, gene_with_variation):
        tx = gene_with_variation[TRANSCRIPT]
        assert len(tx.variants.wt_nt_seq) == tx.cds_length

    def test_wt_aa_seq_matches_transcript_aa(self, gene_with_variation):
        tx = gene_with_variation[TRANSCRIPT]
        assert len(tx.variants.wt_aa_seq) == len(tx.AAseq)

    def test_mean_sites_positive(self, gene_with_variation):
        tx = gene_with_variation[TRANSCRIPT]
        assert tx.variants.mean_ssites > 0
        assert tx.variants.mean_nssites > 0

    def test_cds_seq_array_shape(self, gene_with_variation):
        tx = gene_with_variation[TRANSCRIPT]
        v = tx.variants
        n_samples = len(gene_with_variation.samples)
        assert v.cds_seq_array.shape[0] == tx.cds_length
        assert v.cds_seq_array.shape[1] == n_samples
        assert v.cds_seq_array.shape[2] == 2

    def test_synarr_nonarr_shape_matches_cds(self, gene_with_variation):
        tx = gene_with_variation[TRANSCRIPT]
        v = tx.variants
        assert v.synarr.shape == v.cds_seq_array.shape
        assert v.nonarr.shape == v.cds_seq_array.shape

    def test_statistics_dataframe_initialized(self, gene_with_variation):
        tx = gene_with_variation[TRANSCRIPT]
        assert isinstance(tx.variants.statistics, pd.DataFrame)
        assert list(tx.variants.statistics.index) == ['pi', 'theta', 'tajD', 'Fst']
        assert list(tx.variants.statistics.columns) == ['gene', 'cds', 'aa', 'ss', 'ns']


# ---------------------------------------------------------------------------
# calculate_statistics / statistics DataFrame
# ---------------------------------------------------------------------------

class TestCalculateStatistics:
    def test_statistics_is_dataframe(self, gene_with_stats):
        tx = gene_with_stats[TRANSCRIPT]
        assert isinstance(tx.statistics, pd.DataFrame)

    def test_statistics_index(self, gene_with_stats):
        assert list(gene_with_stats[TRANSCRIPT].statistics.index) == ['pi', 'theta', 'tajD', 'Fst']

    def test_statistics_columns(self, gene_with_stats):
        assert list(gene_with_stats[TRANSCRIPT].statistics.columns) == ['gene', 'cds', 'aa', 'ss', 'ns']

    def test_pi_non_negative(self, gene_with_stats):
        row = gene_with_stats[TRANSCRIPT].statistics.loc['pi']
        for col, val in row.items():
            assert float(val) >= 0, f"pi_{col} is negative: {val}"

    def test_theta_non_negative(self, gene_with_stats):
        row = gene_with_stats[TRANSCRIPT].statistics.loc['theta']
        for col, val in row.items():
            assert float(val) >= 0, f"theta_{col} is negative: {val}"

    ## compare pi to scikit allel's implementation on the same data
    def test_pi_matches_skal(self, gene_with_stats):
        tx = gene_with_stats[TRANSCRIPT]
        v = tx.variants
        pi_calculated = float(tx.statistics.loc['pi', 'gene'])
        pi_skal = allel.sequence_diversity(v.positions, v.gt_array.count_alleles())
        assert np.isclose(pi_calculated, pi_skal, rtol=1e-5), f"Calculated pi {pi_calculated} does not match scikit-allel pi {pi_skal}"

    def test_tajD_matches_skal(self, gene_with_stats):
        tx = gene_with_stats[TRANSCRIPT]
        v = tx.variants
        tajD_calculated = float(tx.statistics.loc['tajD', 'gene'])
        tajD_skal = allel.tajima_d(v.gt_array.count_alleles())
        assert np.isclose(tajD_calculated, tajD_skal, rtol=1e-5), f"Calculated Tajima's D {tajD_calculated} does not match scikit-allel Tajima's D {tajD_skal}"

    def test_theta_matches_skal(self, gene_with_stats):
        tx = gene_with_stats[TRANSCRIPT]
        v = tx.variants
        theta_calculated = float(tx.statistics.loc['theta', 'gene'])
        theta_skal = allel.watterson_theta(v.positions, v.gt_array.count_alleles())
        assert np.isclose(theta_calculated, theta_skal, rtol=1e-5), f"Calculated theta {theta_calculated} does not match scikit-allel Watterson's theta {theta_skal}"

    def test_fst_matches_skal(self, gene_with_stats):
        tx = gene_with_stats[TRANSCRIPT]
        v = tx.variants
        fst_calculated = float(tx.statistics.loc['Fst', 'gene'])
        # group samples by population from metadata
        locs = gene_with_stats.locs
        cl = HDBSCAN(min_cluster_size=2, metric='haversine').fit_predict(np.radians(locs))
        locs = locs[cl != -1]
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
        a, b, c = allel.weir_cockerham_fst(v.gt_array, clinds)
        fst_skal = np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))
        assert np.isclose(fst_calculated, fst_skal, rtol=1e-5), f"Calculated Fst {fst_calculated} does not match scikit-allel Fst {fst_skal}"
# ---------------------------------------------------------------------------
# Statistics preservation across repeated calculate_statistics calls
# ---------------------------------------------------------------------------

class TestStatsPreservation:
    def test_incremental_stats_preserved(self):
        """Pi calculated first should survive a subsequent call that only calculates theta."""
        g = Gene(
            GENE_ID, vcf=VCF, annotation=ANNOTATION, fasta=FASTA, metadata=METADATA,
            statistics={'pi': True, 'theta': False, 'tajD': False, 'Fst': False, 'pairwiseFst': False},
        )
        g.fetch_gene_coordinates()
        g.fetch_gene_transcripts()
        g.fetch_variation()
        g.calculate_statistics()

        tx = g[TRANSCRIPT]
        pi_val = float(tx.statistics.loc['pi', 'gene'])
        assert np.isfinite(pi_val)
        assert pd.isna(tx.statistics.loc['theta', 'gene'])

        # Calculate only theta on the same object
        g.calculate_statistics(statistics=['theta'])

        assert float(tx.statistics.loc['pi', 'gene']) == pi_val   # pi unchanged
        assert np.isfinite(float(tx.statistics.loc['theta', 'gene']))  # theta now filled

    def test_second_call_does_not_wipe_first(self):
        """Calling calculate_statistics twice with different flags should accumulate rows."""
        g = Gene(
            GENE_ID, vcf=VCF, annotation=ANNOTATION, fasta=FASTA, metadata=METADATA,
            statistics={'pi': False, 'theta': False, 'tajD': False, 'Fst': False, 'pairwiseFst': False},
        )
        g.fetch_gene_coordinates()
        g.fetch_gene_transcripts()
        g.fetch_variation()

        g.calculate_statistics(statistics=['pi'])
        g.calculate_statistics(statistics=['tajD'])

        tx = g[TRANSCRIPT]
        assert np.isfinite(float(tx.statistics.loc['pi', 'gene']))
        assert np.isfinite(float(tx.statistics.loc['tajD', 'gene']))


# ---------------------------------------------------------------------------
# keep_longest_transcript
# ---------------------------------------------------------------------------

class TestKeepLongestTranscript:
    def _fresh_gene_with_transcripts(self):
        g = Gene(GENE_ID, vcf=VCF, annotation=ANNOTATION, fasta=FASTA)
        g.fetch_gene_coordinates()
        g.fetch_gene_transcripts()
        return g

    def test_reduces_to_one_transcript(self):
        g = self._fresh_gene_with_transcripts()
        g.keep_longest_transcript()
        assert len(g.transcripts) == 1

    def test_kept_transcript_is_transcript_type(self):
        g = self._fresh_gene_with_transcripts()
        g.keep_longest_transcript()
        tx = list(g.transcripts.values())[0]
        assert isinstance(tx, Transcript)

    def test_kept_transcript_has_max_cds_length(self):
        g = self._fresh_gene_with_transcripts()
        all_lengths = {tid: tx.cds_length for tid, tx in g.transcripts.items()}
        g.keep_longest_transcript()
        kept_id = list(g.transcripts.keys())[0]
        assert all_lengths[kept_id] == max(all_lengths.values())


# ---------------------------------------------------------------------------
# save_to_df
# ---------------------------------------------------------------------------

class TestSaveToDf:
    def test_output_file_created(self, gene_with_stats, tmp_path):
        out = str(tmp_path / "output.txt")
        gene_with_stats.save_to_df(TRANSCRIPT, output_file=out, append=False)
        assert pathlib.Path(out).exists()

    def test_output_has_header_and_data_row(self, gene_with_stats, tmp_path):
        out = str(tmp_path / "output_header.txt")
        gene_with_stats.save_to_df(TRANSCRIPT, output_file=out, append=False)
        lines = [l for l in pathlib.Path(out).read_text().splitlines() if l.strip()]
        assert len(lines) >= 2
        assert 'Gene' in lines[0]
        assert 'pi_gene' in lines[0]

    def test_header_contains_all_stat_columns(self, gene_with_stats, tmp_path):
        out = str(tmp_path / "output_cols.txt")
        gene_with_stats.save_to_df(TRANSCRIPT, output_file=out, append=False)
        header = pathlib.Path(out).read_text().splitlines()[0]
        for stat in ('pi', 'theta', 'tajD', 'Fst'):
            for level in ('gene', 'cds', 'aa', 'ss', 'ns'):
                assert f'{stat}_{level}' in header, f"Missing column {stat}_{level} in header"

    def test_append_adds_rows(self, gene_with_stats, tmp_path):
        out = str(tmp_path / "output_append.txt")
        gene_with_stats.save_to_df(TRANSCRIPT, output_file=out, append=False)
        gene_with_stats.save_to_df(TRANSCRIPT, output_file=out, append=True)
        lines = [l for l in pathlib.Path(out).read_text().splitlines() if l.strip()]
        assert len(lines) == 3  # header + 2 data rows

    def test_gene_name_in_data_row(self, gene_with_stats, tmp_path):
        out = str(tmp_path / "output_name.txt")
        gene_with_stats.save_to_df(TRANSCRIPT, output_file=out, append=False)
        content = pathlib.Path(out).read_text()
        assert GENE_ID in content


# ---------------------------------------------------------------------------
# __repr__ methods
# ---------------------------------------------------------------------------

class TestRepr:
    def test_gene_repr_contains_name(self, gene_with_stats):
        assert GENE_ID in repr(gene_with_stats)

    def test_gene_repr_contains_chromosome(self, gene_with_stats):
        assert '2L' in repr(gene_with_stats)

    def test_transcript_repr_contains_id(self, gene_with_stats):
        tx = gene_with_stats[TRANSCRIPT]
        assert TRANSCRIPT in repr(tx)

    def test_transcript_repr_contains_stats_when_calculated(self, gene_with_stats):
        tx = gene_with_stats[TRANSCRIPT]
        assert 'pi' in repr(tx)

    def test_variants_repr(self, gene_with_variation):
        tx = gene_with_variation[TRANSCRIPT]
        r = repr(tx.variants)
        assert 'Variants' in r
        assert 'positions' in r
