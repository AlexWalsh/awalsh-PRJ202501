# R 4.4.1 #

# prepare cibersort data ---------------------------------

# control data from (control cortex tissue) Zhou et al., 2020 , https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html #


library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)

theme_set(theme_cowplot())


# set local directory --------------------------------------------

setwd("C:/Users/alex_/Desktop/Work/2025/Human-microglia-work-submission/final_analysis/human/mRNA/cibersort_prepare")

# process Zhou (2020) seurat object ------------------------------


zhou_seurat <- readRDS("Zhou_2020_control.rds")

DimPlot(zhou_seurat, group.by='cell_type', label=TRUE) + ggtitle('Zhou et al Control Cortex') + NoLegend()

cell_counts_df <- GetAssayData(zhou_seurat, assay = "RNA", layer = "counts")

set.seed(40)

downsampled_columns <- sample(ncol(cell_counts_df), 5000) # sampling 5000 cells
downsampled_df <- data.frame(cell_counts_df[, downsampled_columns])

zhou_seurat@meta.data$cell_type

zhou_seurat_cell_type_labels <- t(zhou_seurat@meta.data$cell_type[downsampled_columns])

colnames(downsampled_df) <- zhou_seurat_cell_type_labels

# annotate gene_names to ensembl id ------------------------------

GRCh38 <- read_rds("GRCh38_genes.rds")

zhou_mix <- cbind(gene_name = rownames(downsampled_df), downsampled_df)

zhou_annotate <- merge(x = zhou_mix, y = GRCh38[, c("gene_id", "gene_name")], by = "gene_name", all.x = TRUE) |> relocate(gene_name, gene_id)

zhou_annotate <- zhou_annotate[,-c(1)] # remove gene_name column

sum((is.na(zhou_annotate$gene_id))) # 3506 

zhou_annotate <- zhou_annotate[!is.na(zhou_annotate$gene_id), ]


# final formatting -------------------------------------------

zhou_mix_final <- rbind(colnames(zhou_mix), zhou_annotate)
zhou_mix_final[1,1] <- "gene_id"

rownames(zhou_mix_final) <- NULL
colnames(zhou_mix_final) <- NULL

write.table(zhou_mix_final, file = "zhou_reference.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)




# cibersort code - run as docker container ---------------------


#--refsample zhou_reference.txt \
#--fraction 0.25 \
#--rmbatchSmode TRUE \
#--mixture hsa_rawcounts_cibersort.txt \
#--sampling 1 \
#--filter TRUE \
#--remake True
