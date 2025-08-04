# Population genetics comparison of Anopheles malaria vectors

## Basic script usage:

The `wgs_stats.py` script calculates pi (nucleotide diversity), theta hat (Watterson's estimator), and Tajima's D
for a given gene, considering:
- the entire gene
- the coding sequence
- the amino acid sequence
- silent sites in the coding sequence
- nonsynonymous sites in the coding sequence

The script requires a gene ID, a VCF file, an annotation in `.db` format (which can be created using `parse_gff.py`), and a genome-level FASTA file. 

__Example execution:__

    python wgs_stats.py --gene <gene_name> \
                        --vcf <vcf_path>\
                        --annotation <annotation.db path> \
                        --fasta <fasta path> \

__Example printed output:__

    gene-wide pi: 0.026461257859243054 cds_pi: 0.004484740740526833 aa_pi: 0.0029415135832686753 ss_pi: 0.009924957573901851 ns_pi: 0.0030197803483858632
    gene-wide theta: 0.11677863181375378 cds_theta: 0.022045760862240794 aa_theta: 0.016981474980921216 ss_theta: 0.03224173578054607 ns_theta: 0.020261158990112015
    gene-wide Tajima's D: -2.1546476407943578 cds_D: -2.215106564978898 aa_D: -2.2707703885891455 ss_D: -1.8660249387103294 ns_D: -2.3478453225818217

Tab-formatted output can be saved to a file using the `--output` flag, which will either create a 
file if one does not exist or append statistics to an existing file. The output columns are:

    Gene	Chromosome	Transcript	gene_pi	cds_pi	aa_pi	ss_pi	ns_pi	gene_theta	cds_theta	aa_theta	ss_theta	ns_theta	gene_d	cds_d	aa_d	ss_d	ns_d	

## Using the `Gene` class:

Further information about nucleotide and genetic variation is contained within the `Gene` class
and can be accessed in the `.transcripts` object after using `.fetch_variation()`.

    >>> from wgs_stats import Gene
    >>> GENE = Gene(name="AGAP007036")
    >>> GENE.fetch_gene_coordinates(annotation)
    ('AgamP4_2L', 41271509, 41272901)
    >>> GENE.fetch_variation(
    ...     vcf_file=vcf_file.replace('CHROMOSOME', GENE.chromosome.split('_')[-1]),
    ...     annotation=annotation,
    ...     fasta=fasta
    ... )
    Fetching transcripts for gene AGAP007036 on chromosome 2L
    Fetching coding sequence for transcript AGAP007036-RA


    >>> for k, v in GENE.transcripts['AGAP007036-RA'].items():
    ...     print(k, type(v))
    ... 
    id <class 'str'> # ID of transcript
    cds_id <class 'list'> # list of coding sequence IDs
    cds_start <class 'list'> # list of coding sequence start positions
    cds_end <class 'list'> # list of coding sequence stop positions
    cds_seq <class 'list'> # list of strings of coding sequences as nucleotides
    strand <class 'str'> # + or -
    sequences <class 'dict'> # dictionary of sequence variation data
    mean_nssites <class 'numpy.float64'> # mean number of nonsynonymous sites across all coding sequences
    mean_sssites <class 'numpy.float64'> # # mean number of synonymous sites across all coding sequences

Information about loaded genotype variants can be accessed in the `.transcripts['<TRANSCRIPT_ID>']['sequences']` dictionary.
By default, all coding transcripts are loaded into `.transcripts`. (TODO: add method to just load one transcript)

    >>> for k, v in GENE.transcripts['AGAP007036-RA']['sequences'].items():
    ...     print(k, type(v))
    ... 
    genotype_array <class 'allel.model.ndarray.GenotypeArray'> # scikit-allel GenotypeArray of all data
    positions <class 'numpy.ndarray'> # array of variant positions (genome-based)
    wt_nt_seq <class 'numpy.ndarray'> # wild-type coding sequence (nucleotide)
    wt_aa_seq <class 'numpy.ndarray'> # wild-type amino acid sequence
    nt_seq_array <class 'numpy.ndarray'> # array of genotypes as coding sequence nucleotide variation 
                                         # (second axis corresponds to sample locations in genotype_array)
    aa_seq_array <class 'numpy.ndarray'> # array of genotypes as amino acid variation
                                         # (second axis corresponds to sample locations in genotype_array)
    synonymous_array <class 'numpy.ndarray'> # array of variant genotypes at synonymous sites
    nonsynonymous_array <class 'numpy.ndarray'> # array of variant genotypes at nonsynonymous sites

Statistics can be directly calculated on loaded genotypes using `.calculate_statistics('<TRANSCRIPT_ID>'):

    ```
    >>> GENE.calculate_statistics('AGAP007036-RA')
    gene-wide pi: 0.0023784062791499727 cds_pi: 2.639262834664575e-06 aa_pi: 3.571741754680041e-06 ss_pi: 0 ns_pi: 3.571741754680041e-06
    gene-wide theta: 0.11677863181375378 cds_theta: 9.073708765637434e-05 aa_theta: 0.00012279544137237239 ss_theta: 0 ns_theta: 0.00012279544137237239
    gene-wide Tajima's D: -0.9842000160717197 cds_D: nan aa_D: nan ss_D: 0 ns_D: nan
    ```