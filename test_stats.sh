#!/bin/bash

vcf="/sietch_colab/data_share/Ag1000G/Ag3.0/vcf/phased_vcf/gamb_colu_arab/ag3_gamb1470_colu507_arab368.CHROMOSOME.phased.recode.vcf.gz"
annotation="data/VectorBase-68_AgambiaePest.db"
fasta="data/VectorBase-68_AgambiaePEST_Genome.fasta"
samples=$(echo $(cat "/sietch_colab/data_share/Ag1000G/Ag3.0/vcf/sample_logs/gamb_inds.n1470.keep.txt"))

for gene in AGAP007033; do # AGAP006184 AGAP006914 AGAP007035 AGAP007036 AGAP010489 AGAP004370 AGAP006914 AGAP007033 AGAP003879 AGAP001587 AGAP029989 AGAP006184; do
    python scripts/wgs_stats.py \
        --vcf $vcf \
        --gene $gene \
        --annotation $annotation \
        --fasta $fasta \
        --proteinfasta "data/VectorBase-68_AgambiaePEST_AnnotatedProteins.fasta" \
        --samples $samples #\
        #--output "ag_gene_stats.tsv"
done