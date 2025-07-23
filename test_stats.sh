#!/bin/bash

vcf="/sietch_colab/data_share/Ag1000G/Ag3.0/vcf/phased_vcf/gamb/gamb.2L.phased.n1470.derived.vcf.gz"
chrom="2L"

#AGAP005948
#AgamP4_2L:23,985,166..23,988,714(+)


annotation="VectorBase-68_AgambiaePest.db"

python scripts/wgs_stats.py \
    --vcf $vcf \
    --gene AGAP005948\
    --annotation $annotation