#!/bin/bash

vcf="/sietch_colab/data_share/Ag1000G/Ag3.0/vcf/unphased_vcf/gamb/ag1000g.agam_n1470.merged_allsites.vcf.gz"
filter="/sietch_colab/data_share/Ag1000G/Ag3.0/vcf/unphased_vcf/gamb/ag1000g.agam_n1470.merged_variants.sitefilt.vcf.gz"
annotation="data/gamb/VectorBase-68_AgambiaePEST.db"
fasta="data/gamb/VectorBase-68_AgambiaePEST_Genome.fasta"
samples=$(echo $(cat "/sietch_colab/data_share/Ag1000G/Ag3.0/vcf/sample_logs/gamb_inds.n1470.keep.txt"))

for gene in AGAP007033; do # AGAP006184 AGAP006914 AGAP007035 AGAP007036 AGAP010489 AGAP004370 AGAP006914 AGAP007033 AGAP003879 AGAP001587 AGAP029989 AGAP006184; do
    python scripts/wgs_stats.py \
        --vcf $vcf \
        --gene $gene \
        --annotation $annotation \
        --fasta $fasta \
        --proteinfasta "data/gamb/VectorBase-68_AgambiaePEST_AnnotatedProteins.fasta" \
        --filter $filter
        #--samples $samples #\
        #--output "ag_gene_stats.tsv"
done