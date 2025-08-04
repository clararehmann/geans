#!/bin/bash

vcf="/sietch_colab/data_share/Ag1000G/Ag3.0/vcf/phased_vcf/gamb/gamb.CHROMOSOME.phased.n1470.derived.vcf.gz"
annotation="data/VectorBase-68_AgambiaePest.db"
fasta="data/VectorBase-68_AgambiaePEST_Genome.fasta"

for gene in AGAP007033 AGAP006184 AGAP006914 AGAP007035 AGAP007036 AGAP010489 AGAP004370 AGAP006914 AGAP007033 AGAP003879 AGAP001587 AGAP029989 AGAP006184; do
    python scripts/wgs_stats.py \
        --vcf $vcf \
        --gene $gene \
        --annotation $annotation \
        --fasta $fasta \
        --proteinfasta "data/VectorBase-68_AgambiaePEST_AnnotatedProteins.fasta" \
        --output "ag_gene_stats.tsv"
done