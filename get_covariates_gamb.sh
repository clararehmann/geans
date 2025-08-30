#!/bin/bash

vcf="/sietch_colab/data_share/Ag1000G/Ag3.0/vcf/unphased_vcf/gamb/ag1000g.agam_n1470.merged_allsites.vcf.gz"
annotation="data/gamb/VectorBase-68_AgambiaePEST.db"
fasta="data/gamb/VectorBase-68_AgambiaePEST_Genome.fasta"
filter="/sietch_colab/data_share/Ag1000G/Ag3.0/vcf/unphased_vcf/gamb/ag1000g.agam_n1470.merged_variants.sitefilt.vcf.gz"

echo -e "Gene\tChromosome\tStart\tEnd\tLength\tGC_Content" > "out/gamb/gene_covariates.txt"

genes="data/gamb/agPEST_geneIDlist.txt"

while read gene; do
    python scripts/get_covariates.py --name $gene --annotation $annotation --fasta $fasta >> "out/gamb/gene_covariates.txt"
done < $genes

#for gene in AGAP007033; do #; do #AGAP006184 AGAP006914 AGAP007035 AGAP007036 AGAP010489 AGAP004370 AGAP006914 AGAP007033 AGAP003879 AGAP001587 AGAP029989 AGAP006184; do
#    python scripts/wgs_stats.py \
#        --vcf $vcf \
#        --gene $gene \
#        --annotation $annotation \
#        --fasta $fasta \
#        --proteinfasta "data/gamb/VectorBase-68_AgambiaePEST_AnnotatedProteins.fasta" \
#        --filter $filter
#        --output "test_out.txt"
        #--samples $samples #\
        #--output "ag_gene_stats.tsv"
#done