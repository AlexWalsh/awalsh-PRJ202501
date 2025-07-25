# R 4.3.1 ---------------------------------------------------

setwd("C:/Users/alex_/Desktop/Work/2025/Human-microglia-work-submission/final_analysis/human/mRNA-miRNA-network")

#### corr.test  ####

library(dplyr)
library(tidyverse)
library(psych)

### load in rds file ###


miRNA.lcpm <- readRDS("C:/Users/alex_/Desktop/Work/2025/Human-microglia-work-submission/final_analysis/human/smRNA/output-files/hsa_smRNA_lcpm_resid.rds")
mRNA.lcpm <- readRDS("C:/Users/alex_/Desktop/Work/2025/Human-microglia-work-submission/final_analysis/human/mRNA/output-files/hsa_mRNA_lcpm_resid.rds")


###format data frames ---------------------------------------------

miRNAcounts <- miRNA.lcpm %>% dplyr::select(matches("CD45|gene"))

mRNAcounts <- mRNA.lcpm %>% dplyr::select(matches("CD45|gene"))

### convert mRNAs from ENSEMBL to gene names and reformat  ------------------------------------------

GRCh38 <- read_rds("GRCh38_genes.rds")

genes <- as_tibble(mRNAcounts$geneid)

colnames(genes)[1] <- "gene_id" # rename column

annotate <- left_join(genes, GRCh38, by = "gene_id") |> dplyr::select("gene_id", "gene_name")

annotate_final <- annotate %>% mutate(full_name = paste(gene_id, gene_name, sep = "|"))

mRNAcounts <- cbind(annotate_final$full_name, mRNAcounts[,-c(1)])

colnames(mRNAcounts)[1] <- "geneid"

## format so that gene ids are rownames ------------------------

geneid <- mRNAcounts$geneid
mRNAcounts <- dplyr::select(mRNAcounts, -geneid)
rownames(mRNAcounts) <- geneid

miRNAid <- miRNAcounts$geneid
miRNAcounts <- dplyr::select(miRNAcounts, -geneid)
rownames(miRNAcounts) <- miRNAid

colnames(miRNAcounts) == colnames(mRNAcounts) # check that structure is correct

#### using corr.test ####

### requires two matrices, with same row numbers ###

dim(miRNAcounts)
dim(mRNAcounts) 

# transpose both matrices #

t.miRNAcounts <- t(miRNAcounts)
t.mRNAcounts <- t(mRNAcounts)


myCorr <- corr.test(t.miRNAcounts, t.mRNAcounts, method = 'pearson', adjust = 'none', ci = FALSE) 

saveRDS(myCorr, file = "hsa_miRNA_mRNA_microglia.corr.rds")

writeLines(capture.output(sessionInfo()), "sessioninfo_corr.txt")
