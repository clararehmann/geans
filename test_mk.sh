#!/bin/bash

gene1='AGAP029767'
species1='gamb'
vcf1='/sietch_colab/data_share/Ag1000G/Ag3.0/vcf/unphased_vcf/gamb/ag1000g.agam_n1470.merged_allsites.vcf.gz'
filter1='/sietch_colab/data_share/Ag1000G/Ag3.0/vcf/unphased_vcf/gamb/ag1000g.agam_n1470.merged_variants.sitefilt.vcf.gz'
fasta1='data/gamb/VectorBase-68_AgambiaePEST_Genome.fasta'
annotation1='data/gamb/VectorBase-68_AgambiaePEST.db'

gene2='AFUN2_003877'
species2='afun'
vcf2='/sietch_colab/crehmann/vo_afun_release/v1.0/merged_vcf/afun_merged.vcf.gz'
filter2='/sietch_colab/crehmann/vo_afun_release/v1.0/site_filters/dt_20200416/vcf/funestus/merged_funestus_sitefilters.vcf.gz'
fasta2='data/afun/VectorBase-68_AfunestusAfunGA1_Genome.fasta'
annotation2='data/afun/VectorBase-68_AfunestusAfunGA1.db'

outprefix='test_'

python scripts/mk_test.py --gene1 $gene1 --species1 $species1 --vcf1 $vcf1 --filter1 $filter1 --fasta1 $fasta1 --annotation1 $annotation1 --gene2 $gene2 --species2 $species2 --vcf2 $vcf2 --filter2 $filter2 --fasta2 $fasta2 --annotation2 $annotation2 --outprefix $outprefix