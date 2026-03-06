# geans: gene-based population genetics statistics

`geans` calculates per-gene population genetics statistics ( nucleotide diversity (π), Watterson's θ, Tajima's D, and Weir-Cockerham Fst) from a VCF file.
Statistics are calculated for the whole gene, the coding sequence, the amino acid sequence, synonymous sites, and nonsynonymous sites. It can be used as a Python module or run directly from the command line.

---

## Installation

**Prerequisites**: Python ≥ 3.9, a (preferably bgzipped and tabix-indexed) VCF, a gffutils annotation database (which can be generated using `utilityscripts/parse_gff.py`), and a reference genome FASTA.

```bash
git clone https://github.com/clararehmann/geans
cd geans
uv pip install -e .
```

An example VCF for *Anopheles gambiae* gene AGAP005023 with its associated filter, GFF file, and metadata are available in the `example_data` subdirectory. To run example tests, you will need to obtain the *Anopheles gambiae* reference genome FASTA:

```bash
wget -P example_data/ https://vectorbase.org/a/service/raw-files/release-68/AgambiaePEST/fasta/data/VectorBase-68_AgambiaePEST_Genome.fasta
```

and construct a gffutils database from the GFF file using the provided utility script:

```bash
python utilityscripts/parse_gff.py --input_file example_data/test.gff --output_file example_data/test.db
```

---

## Command-line usage

```bash
geans --id GENE_ID \
      --vcf variants.vcf.gz \
      --annotation genome.db \
      --fasta genome.fa
```

### Arguments

| Flag | Required | Description |
|------|----------|-------------|
| `--id` | yes | Gene identifier (must match the annotation database) |
| `--vcf` | yes | Path to bgzipped, tabix-indexed VCF |
| `--annotation` | yes | Path to gffutils annotation database |
| `--fasta` | yes | Path to reference genome FASTA |
| `--metadata` | no | Tab-separated sample metadata file (must contain a `sample_id` column) |
| `--filter` | no | Filter VCF; only sites present in this file are retained |
| `--protein-fasta` | no | Protein FASTA for amino acid sequence validation |
| `--stats` | no | Statistics to calculate (default: `pi theta tajD Fst`). Choices: `pi theta tajD Fst pairwiseFst` |
| `--longest` | no | Only analyze the longest transcript (flag; default: all transcripts) |
| `--output` | no | Write statistics to this tab-separated file |

### Examples

```bash
# Minimal — calculate pi, theta, tajD, Fst for all transcripts
geans --id AGAP005023 \
      --vcf ag1000g.vcf.gz \
      --annotation VectorBase.db \
      --fasta AgamP4.fa \
      --metadata samples.txt

# Calculate only pi and Tajima's D for the longest transcript
geans --id AGAP005023 \
      --vcf ag1000g.vcf.gz \
      --annotation VectorBase.db \
      --fasta AgamP4.fa \
      --metadata samples.txt \
      --stats pi tajD \
      --longest

# Save output to a TSV
geans --id AGAP005023 \
      --vcf ag1000g.vcf.gz \
      --annotation VectorBase.db \
      --fasta AgamP4.fa \
      --metadata samples.txt \
      --output results.txt
```

### Output

Statistics are printed to stdout as a formatted table per transcript. If `--output` is given, a tab-separated file is written with one row per transcript and columns:

```
Gene  Chromosome  Transcript  Length  GC_Content  pi_gene  pi_cds  pi_aa  pi_ss  pi_ns  theta_gene  ...
```

Pairwise Fst (if calculated) is written to a separate file named `<output>_<transcriptID>_pairwiseFst.txt` where each row and column 
represents the Fst between index locations.

---

## Python API

### Basic workflow

```python
import geans

gene = geans.Gene(
    'AGAP005023',
    vcf='ag1000g.vcf.gz',
    annotation='VectorBase.db',
    fasta='AgamP4.fa',
    metadata='samples.txt',
)

gene.fetch_gene_coordinates()
gene.fetch_gene_transcripts()
gene.fetch_variation()
gene.calculate_statistics()

print(gene)
```

```
Gene: AGAP005023
  Location: AgamP4_2L:8375887-8376198  |  Length: 312 bp  |  GC: 0.596
  Transcripts: 1 (AGAP005023-RA)
  Statistics: pi, theta, tajD, Fst [pairwiseFst disabled]
```

### Accessing transcripts and statistics

```python
# Count transcripts
len(gene)

# Iterate over transcripts
for tx in gene:
    print(tx)
    print(tx.statistics)     # pandas DataFrame

# Index by transcript ID
tx = gene['AGAP005023-RA']

# Direct DataFrame indexing
tx.statistics.loc['pi', 'cds']    # CDS-level nucleotide diversity
tx.statistics.loc['tajD']         # Full Tajima's D row across all levels
```

### Incremental statistics

Statistics are preserved across calls, so you can add a new statistic without recomputing existing ones:

```python
gene.calculate_statistics()               # computes pi, theta, tajD, Fst
gene.calculate_statistics(['pairwiseFst']) # adds pairwise Fst; pi/theta/tajD/Fst unchanged
```

### Only the longest transcript

```python
gene.fetch_gene_transcripts()
gene.keep_longest_transcript()
gene.fetch_variation()
gene.calculate_statistics()
```

### Saving to file

```python
# Write header + first transcript
gene.save_to_df('AGAP005023-RA', output_file='results.txt', append=False)

# Append additional transcripts (if they exist)
# this code won't run because AGAP005023 only has one transcript
gene.save_to_df('AGAP005023-RB', output_file='results.txt', append=True)
```

### Controlling which statistics are calculated

Pass a `statistics` dict at construction or a list to `calculate_statistics`:

```python
# At construction — only pi and theta
gene = geans.Gene('AGAP005023', ..., statistics={'pi': True, 'theta': True, 'tajD': False, 'Fst': False})

# Or selectively
gene.calculate_statistics(['pi', 'theta'])
```

---

## Statistics

| Statistic | Description |
|-----------|-------------|
| `pi` | Nucleotide diversity — mean pairwise differences per site |
| `theta` | Watterson's θ — diversity estimated from segregating sites |
| `tajD` | Tajima's D — neutrality test statistic |
| `Fst` | Weir-Cockerham Fst — genetic differentiation between sampling locations (requires metadata with coordinates) |
| `pairwiseFst` | Hudson's pairwise Fst between all pairs of sampling locations (requires metadata with coordinates) |

Each statistic is computed at five levels:

| Level | Description |
|-------|-------------|
| `gene` | All variant sites within the gene region |
| `cds` | Coding sequence sites only |
| `aa` | Amino acid (protein) level |
| `ss` | Synonymous sites |
| `ns` | Nonsynonymous sites |

---

## Input requirements

### VCF
Chromosome names in the VCF must match those in the annotation database (the last `_`-delimited field is stripped when querying, e.g. `AgamP4_2L` → `2L`).

### Annotation database
A gffutils SQLite database built from a GFF3 file. Create one with:

```bash
python scripts/parse_gff.py --input_file genome.gff3 --output_file genome.db
```

### Metadata TSV
Required for Fst calculations. Must contain a `sample_id` column matching sample IDs in the VCF. Geographic coordinates should be in `longitude`/`latitude` columns (or `x`/`y` for projected coordinates). If no coordinate columns are found, Fst is skipped silently.

### Reference FASTA
Standard indexed reference genome. Used to retrieve wild-type sequences for each transcript's coding sequence.

---

## Utility scripts

| Script | Description |
|--------|-------------|
| `scripts/parse_gff.py` | Build a gffutils annotation database from a GFF3 file |
| `scripts/get_codinggenes.py` | Extract a list of coding gene IDs from a protein FASTA |