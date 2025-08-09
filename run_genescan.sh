#!/bin/bash

# initialize input/output array
# argument order:
# vcf, site filter, annotation db, genome fasta, protein fasta, output base, metadata dir

declare -A mosqs=(
    ["gamb"]="Anopheles gambiae"
    ["colu"]="Anopheles coluzzii"
    ["arab"]="Anopheles arabiensis"
    #["afun"]="Anopheles funestus"
    #["amin"]="Anopheles minimus"
)

declare -A vcfs=(
    ["gamb"]="/projects/kernlab/shared/Ag3.0/vcf/unphased_vcf/gamb/ag1000g.agam_n1470.merged_allsites.vcf.gz"
    ["colu"]="/projects/kernlab/shared/Ag3.0/vcf/unphased_vcf/colu/ag1000g.colu_n507.merged_allsites.vcf.gz"
    ["arab"]="/projects/kernlab/shared/Ag3.0/vcf/unphased_vcf/arab/ag1000g.arab_n368.merged_allsites.vcf.gz"
)

declare -A filters=(
    ["gamb"]="/projects/kernlab/shared/Ag3.0/vcf/unphased_vcf/gamb/ag1000g.agam_n1470.merged_variants.sitefilt.vcf.gz"
    ["colu"]="/projects/kernlab/shared/Ag3.0/vcf/unphased_vcf/gamb/ag1000g.agam_n1470.merged_variants.sitefilt.vcf.gz"
    ["arab"]="data/arab/filters/merged_arab_sitefilters.vcf.gz"
)

declare -A annotations=(
    ["gamb"]="data/gamb/VectorBase-68_AgambiaePest.db"
    ["colu"]="data/gamb/VectorBase-68_AgambiaePest.db"
    ["arab"]="data/gamb/VectorBase-68_AgambiaePest.db"
)

declare -A fastas=(
    ["gamb"]="data/gamb/VectorBase-68_AgambiaePEST_Genome.fasta"
    ["colu"]="data/gamb/VectorBase-68_AgambiaePEST_Genome.fasta"
    ["arab"]="data/gamb/VectorBase-68_AgambiaePEST_Genome.fasta"
)

declare -A proteinfastas=(
    ["gamb"]="data/gamb/VectorBase-68_AgambiaePEST_AnnotatedProteins.fasta"
    ["colu"]="data/gamb/VectorBase-68_AgambiaePEST_AnnotatedProteins.fasta"
    ["arab"]="data/gamb/VectorBase-68_AgambiaePEST_AnnotatedProteins.fasta"
)

declare -A outputs=(
    ["gamb"]="out/gamb/ag_gene_stats"
    ["colu"]="out/colu/colu_gene_stats"
    ["arab"]="out/arab/arab_gene_stats"
)

declare -A metadatas=(
    ["gamb"]="data/gamb/metadata"
    ["colu"]="data/colu/metadata"
    ["arab"]="data/arab/metadata"
)

declare -A genes=(
    ["gamb"]="data/gamb/VectorBase-68_AgambiaePEST_ProteinList.txt"
    ["colu"]="data/gamb/VectorBase-68_AgambiaePEST_ProteinList.txt"
    ["arab"]="data/gamb/VectorBase-68_AgambiaePEST_ProteinList.txt"
)

for m in "${!mosqs[@]}"; do
    if [ ! -d "out/${m}" ]; then
        genes=${genes[$m]}
        # initialize output directory and headers
        mkdir -pv "out/${m}"
        echo -e "Gene\tChromosome\tTranscript\tgene_pi\tcds_pi\taa_pi\tss_pi\tns_pi\tgene_theta\tcds_theta\taa_theta\tss_theta\tns_theta\tgene_d\tcds_d\taa_d\tss_d\tns_d" > "out/${m}/${m}_gene_stats.txt"
        for S in data/${m}/metadata/*; do
            if [[ "all" == *$S* ]]; then
                continue # skip the full metadata file
            fi
            SN=$(basename $S)
            SN=${SN%.txt}  # Remove the .txt extension
            echo -e "Gene\tChromosome\tTranscript\tgene_pi\tcds_pi\taa_pi\tss_pi\tns_pi\tgene_theta\tcds_theta\taa_theta\tss_theta\tns_theta\tgene_d\tcds_d\taa_d\tss_d\tns_d" > "out/${m}/${m}_gene_stats_${SN}.txt"
        done
        for gene in $(cat $genes); do
            echo $gene
            #sbatch run_genescan.batch $gene ${mosqs[$m]}
        done
    fi
done


#mkdir -pv out/arab
#echo -e "Gene\tChromosome\tTranscript\tgene_pi\tcds_pi\taa_pi\tss_pi\tns_pi\tgene_theta\tcds_theta\taa_theta\tss_theta\tns_theta\tgene_d\tcds_d\taa_d\tss_d\tns_d" > out/arab/arab_gene_stats.txt
#for S in data/arab/metadata/*; do
#    SN=$(basename $S)
#    SN=${SN%.txt}  # Remove the .txt extension
#    echo -e "Gene\tChromosome\tTranscript\tgene_pi\tcds_pi\taa_pi\tss_pi\tns_pi\tgene_theta\tcds_theta\taa_theta\tss_theta\tns_theta\tgene_d\tcds_d\taa_d\tss_d\tns_d" > out/arab/arab_gene_stats_${SN}.txt
#done
# submit jobs
#for gene in $(cat data/gamb/VectorBase-68_AgambiaePEST_ProteinList.txt); do
#    sbatch run_genescan_arab.batch $gene
#done