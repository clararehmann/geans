"""Pytest configuration and fixtures for geans integration tests.

Fixture files live in testing/ at the repository root:
  - test.vcf.gz      : VCF for AGAP002866 region (2R:28480000-28495650), 11 samples
  - filter.vcf.gz    : site-filter VCF for the same region
  - agam.db          : gffutils annotation database (VectorBase Ag. gambiae)
  - agam.fasta       : reference FASTA (requires agam.fasta.fai index)
  - metadata.txt     : sample metadata with latitude/longitude columns
"""
import pathlib
import pytest
from geans import Gene

TESTING_DIR = pathlib.Path(__file__).parent.parent / "example_data"

GENE_ID     = "AGAP002866"
TRANSCRIPT  = "AGAP002866-RA"
VCF         = str(TESTING_DIR / "test.vcf.gz")
FILTER_VCF  = str(TESTING_DIR / "filter.vcf.gz")
ANNOTATION  = str(TESTING_DIR / "agam.db")
FASTA       = str(TESTING_DIR / "agam.fasta")
METADATA    = str(TESTING_DIR / "metadata.txt")


@pytest.fixture(scope="session")
def gene_loaded():
    """Session-scoped Gene that has run fetch_gene_coordinates and fetch_gene_transcripts."""
    g = Gene(
        GENE_ID,
        vcf=VCF,
        annotation=ANNOTATION,
        fasta=FASTA,
        metadata=METADATA,
    )
    g.fetch_gene_coordinates()
    g.fetch_gene_transcripts()
    return g


@pytest.fixture(scope="session")
def gene_with_variation(gene_loaded):
    """Session-scoped Gene that has also loaded variation (no site filter)."""
    gene_loaded.fetch_variation()
    return gene_loaded


@pytest.fixture(scope="session")
def gene_with_stats(gene_with_variation):
    """Session-scoped Gene that has calculated all default statistics."""
    gene_with_variation.calculate_statistics()
    return gene_with_variation
