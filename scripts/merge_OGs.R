## script to merge ortholog tables from different sources

library(tidyverse)
library(data.table)

# establish file paths and variables
sp1genes <- "data/gamb/GenesByTaxon_Summary_gamb.txt"
sp2genes <- "data/afun/GenesByTaxon_Summary_afun.txt"
outpath <- "data/gamb_afun_orthologs.txt"
sp1 <- "gamb"
sp2 <- "afun"

# read files and merge orthogroups1
sp1genes <- fread(sp1genes)
sp2genes <- fread(sp2genes)
orths <- inner_join(sp1genes, sp2genes, by="Ortholog Group")

# clean up so that there's only 1:1 orthologs
count1 <- as.data.frame(table(orths$`Gene ID.x`))
excld1 <- count1[count1$Freq>1,]$Var1
count2 <- as.data.frame(table(orths$`Gene ID.y`))
excld2 <- count2[count2$Freq>1,]$Var1
count3 <- as.data.frame(table(orths$`Ortholog Group`))
ogexcld <- count3[count3$Freq>1,]$Var1
orths1 <-orths %>% filter(!(`Gene ID.x` %in% excld1) & 
                          !(`Gene ID.y` %in% excld2) & 
                          !(`Ortholog Group` %in% ogexcld))

# rename columns
names(orths1)[names(orths1) == 'Gene ID.x'] <- paste0(sp1, "_GENE_ID")
names(orths1)[names(orths1) == 'Gene ID.y'] <- paste0(sp2, "_GENE_ID")

# write output
fwrite(orths1, outpath, sep="\t", quote=F, row.names=F)