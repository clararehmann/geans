# Population genetics comparison of Anopheles malaria vectors

Basic script usage:

The `wgs_stats.py` script calculates pi (nucleotide diversity), theta hat (Watterson's estimator), and Tajima's D
for a given gene, considering:
- the entire gene
- the coding sequence
- the amino acid sequence
- silent sites in the coding sequence
- nonsynonymous sites in the coding sequence

The script requires a gene ID, a VCF file, an annotation in `.db` format (which can be created using `parse_gff.py`), and a genome-level FASTA file. 

    Example execution:
    ```
    python wgs_stats.py --gene <gene_name> \
                        --vcf <vcf_path>\
                        --annotation <annotation.db path> \
                        --fasta <fasta path> \
    ```

    Example printed output:
        ```
        gene-wide pi: 0.026461257859243054 cds_pi: 0.004484740740526833 aa_pi: 0.0029415135832686753 ss_pi: 0.009924957573901851 ns_pi: 0.0030197803483858632
        gene-wide theta: 0.11677863181375378 cds_theta: 0.022045760862240794 aa_theta: 0.016981474980921216 ss_theta: 0.03224173578054607 ns_theta: 0.020261158990112015
        gene-wide Tajima's D: -2.1546476407943578 cds_D: -2.215106564978898 aa_D: -2.2707703885891455 ss_D: -1.8660249387103294 ns_D: -2.3478453225818217
        ```