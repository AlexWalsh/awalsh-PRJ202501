# R.4.3.1
#### network analysis, miRNA, mRNA #####


library(dplyr)
library(tidyverse)
library(psych)
library(ggplot2)
library(multiMiR)

#### load in rds objects ###

corr.data <- readRDS("hsa_miRNA_mRNA_microglia.corr.rds")

pearson.corr <- corr.data$r

pearson.corr <- t(pearson.corr) #transpose to get miRNAs as header names

### read in enriched microglia miRNAs ---------------------------

lm_res <- readRDS("C:/Users/alex_/Desktop/Work/2025/Human-microglia-work-submission/final_analysis/human/smRNA/output-files/hsa_sRNA_lm_res.rds")

## subset for significantly enriched miRNAs in microglia - logFC > 1 and adj.p.val <0.05 -----------

tt_treat <- lm_res$tt_treat |> dplyr::filter(logFC > 1 & adj.P.Val < 0.05)

## multimiR -------------------------------------------

CD45Bulk_multimir_results <- get_multimir(org = 'hsa',
                                       mirna = tt_treat$genes,
                                       table = 'validated',
                                       summary = TRUE)

CD45Bulk.multiMir.data <- CD45Bulk_multimir_results@data

unique(CD45Bulk.multiMir.data$experiment)

head(pearson.corr)

### subset pearson correlation object for miRNAs enriched in microglia ------------

CD45BulkDeList_corr <- pearson.corr[, tt_treat[[1]]]

# extract p values -----------------------

pearson.corr.p.values <- corr.data$p

pearson.corr.p.values <- t(pearson.corr.p.values) #transpose to get miRNAs as header names


### reformat multimiR results ----------------------------

CD45Bulk.validated.targets <- as_tibble(CD45Bulk_multimir_results@summary)

CD45Bulk.validated.targets <- CD45Bulk.validated.targets %>% dplyr::select(-mature_mirna_acc, -target_ensembl)

CD45Bulk.validated.targets <- CD45Bulk.validated.targets %>% 
  dplyr::rename(
    miRNA = mature_mirna_id,
    Gene = target_symbol)

#### join all corr results for validated - CD45Bulk_genes_type --------------------

CD45BulkDeList_corr_t <- t(CD45BulkDeList_corr)

CD45BulkDeList_corr_t <- cbind(rownames(CD45BulkDeList_corr_t), data.frame(CD45BulkDeList_corr_t))

colnames(CD45BulkDeList_corr_t)[1] <- c("miRNA") # rename first column

# remove ENS ID from colnames for genes #

colnames(CD45BulkDeList_corr_t) <- sub("^[^.]*\\.", "", colnames(CD45BulkDeList_corr_t))

#### match validated targets with columns of CD45Bulk_miRNAs_type_corr_t

CD45Bulkoutput <- vector("double", length = nrow(CD45Bulk.validated.targets))
for (i in 1:nrow(CD45Bulk.validated.targets)) {
  if(length(paste(CD45BulkDeList_corr_t[as.character(CD45Bulk.validated.targets[i,"miRNA"]), as.character(CD45Bulk.validated.targets[i,"Gene"])])) > 0) { 
    CD45Bulkoutput[[i]] <- paste(CD45BulkDeList_corr_t[as.character(CD45Bulk.validated.targets[i,"miRNA"]), as.character(CD45Bulk.validated.targets[i,"Gene"])])
  }
  
}

CD45Bulkoutput

### append p values ###

CD45BulkDe_p_values <- pearson.corr.p.values[, tt_treat[[1]]]

CD45BulkDe_p_values_t <- t(CD45BulkDe_p_values)

CD45BulkDe_p_values_t <- cbind(rownames(CD45BulkDe_p_values_t), data.frame(CD45BulkDe_p_values_t))

colnames(CD45BulkDe_p_values_t)[1] <- c("miRNA") #rename first column

# remove ENS ID from colnames for genes #

colnames(CD45BulkDe_p_values_t) <- sub("^[^.]*\\.", "", colnames(CD45BulkDe_p_values_t))

### for loop p values to matrix ###

CD45Bulkoutput.p.values <- vector("double", length = nrow(CD45Bulk.validated.targets))
for (i in 1:nrow(CD45Bulk.validated.targets)) {
  if(length(paste(CD45BulkDe_p_values_t[as.character(CD45Bulk.validated.targets[i,"miRNA"]), as.character(CD45Bulk.validated.targets[i,"Gene"])])) > 0) { 
    CD45Bulkoutput.p.values[[i]] <- paste(CD45BulkDe_p_values_t[as.character(CD45Bulk.validated.targets[i,"miRNA"]), as.character(CD45Bulk.validated.targets[i,"Gene"])])
  }
  
}

CD45Bulkoutput.p.values

#adjust p values --------------------

CD45BulkAdj.p.val <- p.adjust(CD45Bulkoutput.p.values, method = c("fdr"))

# combine results for final output --------------

CD45Bulk.validated.targets.pearson <- cbind(CD45Bulk.validated.targets, as.numeric(CD45Bulkoutput), CD45Bulkoutput.p.values, CD45BulkAdj.p.val)

CD45Bulk.validated.targets.pearson <- CD45Bulk.validated.targets.pearson %>% 
  dplyr::rename(
    r_value = "as.numeric(CD45Bulkoutput)",
    Pvalue = CD45Bulkoutput.p.values,
    Adj.P.value = CD45BulkAdj.p.val)


CD45Bulk.validated.targets.pearson <- as_tibble(CD45Bulk.validated.targets.pearson)


###  convert to R^2 ###

CD45Bulk.validated.targets.pearson.final <- CD45Bulk.validated.targets.pearson %>% 
  mutate(R_squared = r_value ^ 2)


### filter on adj. p.value (0.05)  ###

#CD45Bulkfilter <- CD45Bulk.validated.targets.pearson %>% filter(Adj.P.value < 0.05)

write.csv(CD45Bulk.validated.targets.pearson.final, file = "Microglia-enriched-miRNA-mRNA-correlations-no_filter.csv", quote = FALSE)

CD45Bulk.validated.targets.pearson.final <- read.csv("Microglia-enriched-miRNA-mRNA-correlations-no_filter.csv")


# repeat analysis with predicted targets ----------------------------------

## multimiR -------------------------------------------

CD45Bulk_multimir_results <- get_multimir(org = 'hsa',
                                          mirna = tt_treat$genes,
                                          table = 'predicted',
                                          summary = TRUE) # run at default cutoff so only top 20% of interactions

CD45Bulk.multiMir.data <- CD45Bulk_multimir_results@data


head(pearson.corr)

### subset pearson correlation object for miRNAs enriched in microglia ------------

CD45BulkDeList_corr <- pearson.corr[, tt_treat[[1]]]

# extract p values -----------------------

pearson.corr.p.values <- corr.data$p

pearson.corr.p.values <- t(pearson.corr.p.values) #transpose to get miRNAs as header names


### reformat multimiR results ----------------------------

CD45Bulk.predicted.targets <- as_tibble(CD45Bulk_multimir_results@summary)

CD45Bulk.predicted.targets <- CD45Bulk.predicted.targets %>% dplyr::select(-mature_mirna_acc, -target_ensembl)

CD45Bulk.predicted.targets <- CD45Bulk.predicted.targets %>% 
  dplyr::rename(
    miRNA = mature_mirna_id,
    Gene = target_symbol)

#### join all corr results for predicted - CD45Bulk_genes_type --------------------

CD45BulkDeList_corr_t <- t(CD45BulkDeList_corr)

CD45BulkDeList_corr_t <- cbind(rownames(CD45BulkDeList_corr_t), data.frame(CD45BulkDeList_corr_t))

colnames(CD45BulkDeList_corr_t)[1] <- c("miRNA") # rename first column

# remove ENS ID from colnames for genes #

colnames(CD45BulkDeList_corr_t) <- sub("^[^.]*\\.", "", colnames(CD45BulkDeList_corr_t))

#### match predicted targets with columns of CD45Bulk_miRNAs_type_corr_t

CD45Bulkoutput <- vector("double", length = nrow(CD45Bulk.predicted.targets))
for (i in 1:nrow(CD45Bulk.predicted.targets)) {
  if(length(paste(CD45BulkDeList_corr_t[as.character(CD45Bulk.predicted.targets[i,"miRNA"]), as.character(CD45Bulk.predicted.targets[i,"Gene"])])) > 0) { 
    CD45Bulkoutput[[i]] <- paste(CD45BulkDeList_corr_t[as.character(CD45Bulk.predicted.targets[i,"miRNA"]), as.character(CD45Bulk.predicted.targets[i,"Gene"])])
  }
  
}

CD45Bulkoutput

### append p values ###

CD45BulkDe_p_values <- pearson.corr.p.values[, tt_treat[[1]]]

CD45BulkDe_p_values_t <- t(CD45BulkDe_p_values)

CD45BulkDe_p_values_t <- cbind(rownames(CD45BulkDe_p_values_t), data.frame(CD45BulkDe_p_values_t))

colnames(CD45BulkDe_p_values_t)[1] <- c("miRNA") #rename first column

# remove ENS ID from colnames for genes #

colnames(CD45BulkDe_p_values_t) <- sub("^[^.]*\\.", "", colnames(CD45BulkDe_p_values_t))

### for loop p values to matrix ###

CD45Bulkoutput.p.values <- vector("double", length = nrow(CD45Bulk.predicted.targets))
for (i in 1:nrow(CD45Bulk.predicted.targets)) {
  if(length(paste(CD45BulkDe_p_values_t[as.character(CD45Bulk.predicted.targets[i,"miRNA"]), as.character(CD45Bulk.predicted.targets[i,"Gene"])])) > 0) { 
    CD45Bulkoutput.p.values[[i]] <- paste(CD45BulkDe_p_values_t[as.character(CD45Bulk.predicted.targets[i,"miRNA"]), as.character(CD45Bulk.predicted.targets[i,"Gene"])])
  }
  
}

CD45Bulkoutput.p.values

#adjust p values --------------------

CD45BulkAdj.p.val <- p.adjust(CD45Bulkoutput.p.values, method = c("fdr"))

# combine results for final output --------------

CD45Bulk.predicted.targets.pearson <- cbind(CD45Bulk.predicted.targets, as.numeric(CD45Bulkoutput), CD45Bulkoutput.p.values, CD45BulkAdj.p.val)

CD45Bulk.predicted.targets.pearson <- CD45Bulk.predicted.targets.pearson %>% 
  dplyr::rename(
    r_value = "as.numeric(CD45Bulkoutput)",
    Pvalue = CD45Bulkoutput.p.values,
    Adj.P.value = CD45BulkAdj.p.val)


CD45Bulk.predicted.targets.pearson <- as_tibble(CD45Bulk.predicted.targets.pearson)


###  convert to R^2 ###

CD45Bulk.predicted.targets.pearson.final <- CD45Bulk.predicted.targets.pearson %>% 
  mutate(R_squared = r_value ^ 2)



write.csv(CD45Bulk.predicted.targets.pearson.final, file = "Microglia-enriched-miRNA-mRNA-correlations-no_filter_predicted_networks.csv", quote = FALSE)
CD45Bulk.predicted.targets.pearson.final <- read.csv("Microglia-enriched-miRNA-mRNA-correlations-no_filter_predicted_networks.csv")

# KEGG/GO analysis

CD45Bulkfilter.validated <- CD45Bulk.validated.targets.pearson.final %>% dplyr::filter(Adj.P.value < 0.05)
CD45Bulkfilter.predicted <- CD45Bulk.predicted.targets.pearson.final %>% dplyr::filter(Adj.P.value < 0.05)
CD45Bulkfilter.validated.r5 <- CD45Bulk.validated.targets.pearson.final %>% dplyr::filter(Adj.P.value < 0.05 & r_value < -0.5)
CD45Bulkfilter.predicted.r5 <- CD45Bulk.predicted.targets.pearson.final %>% dplyr::filter(Adj.P.value < 0.05 & r_value < -0.5)


validated.go <- goana(CD45Bulkfilter.validated$target_entrez, species = "Hs")
validated.go$p.adj <- p.adjust(validated.go$P.DE, method = c("fdr"))
  
predicted.go <- goana(CD45Bulkfilter.predicted$target_entrez, species = "Hs")
predicted.go$p.adj <- p.adjust(predicted.go$P.DE, method = c("fdr"))

validated.kegg <- kegga(CD45Bulkfilter.validated$target_entrez, species = "Hs")
validated.kegg$p.adj <- p.adjust(validated.kegg$P.DE, method = c("fdr"))

predicted.kegg <- kegga(CD45Bulkfilter.predicted$target_entrez, species = "Hs")
predicted.kegg$p.adj <- p.adjust(predicted.kegg$P.DE, method = c("fdr"))


validated.go.r5 <- goana(CD45Bulkfilter.validated.r5$target_entrez, species = "Hs")
validated.go.r5$p.adj <- p.adjust(validated.go.r5$P.DE, method = c("fdr"))

predicted.go.r5 <- goana(CD45Bulkfilter.predicted.r5$target_entrez, species = "Hs")
predicted.go.r5$p.adj <- p.adjust(predicted.go.r5$P.DE, method = c("fdr"))

validated.kegg.r5 <- kegga(CD45Bulkfilter.validated.r5$target_entrez, species = "Hs")
validated.kegg.r5$p.adj <- p.adjust(validated.kegg.r5$P.DE, method = c("fdr"))

predicted.kegg.r5 <- kegga(CD45Bulkfilter.predicted.r5$target_entrez, species = "Hs")
predicted.kegg.r5$p.adj <- p.adjust(predicted.kegg.r5$P.DE, method = c("fdr"))



## plots for figure

## validated GO/KEGG



topGO.validated.plot <- topGO(validated.go, n = 20, ontology = c("BP")) %>% mutate(log10p = -log10(p.adj)) %>% 
  arrange(desc(log10p)) %>% 
  mutate(Gene_ratio = DE / N)

pdf("Validated_GO_BP_tophits.pdf", width = 12, height = 10)

ggplot(data = topGO.validated.plot, aes(x = reorder(Term, Gene_ratio), y= Gene_ratio, size = DE, color = log10p)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range= c(2,8)) +
  coord_flip() + 
  xlab("GO Terms") +
  ylab("Gene Enrichment") +
  labs(size = 'Gene Count', color = '-log10(FDR)') +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 13)) +
  guides(color = guide_colorbar(order = 2),
         size = guide_legend(order=1))

dev.off()

## KEGG plot ##

topKEGG.validated.plot <- topKEGG(validated.kegg, n = 20) %>% mutate(log10p = -log10(p.adj)) %>% 
  arrange(desc(log10p)) %>% 
  mutate(Gene_ratio = DE / N)

pdf("Validated_KEGG_tophits.pdf", width = 12, height = 10)

ggplot(data = topKEGG.validated.plot, aes(x = reorder(Pathway, Gene_ratio), y= Gene_ratio, size = DE, color = log10p)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range= c(2,8)) +
  coord_flip() + 
  xlab("KEGG Terms") +
  ylab("Gene Enrichment") +
  labs(size = 'Gene Count', color = '-log10(FDR)') +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 13)) +
  guides(color = guide_colorbar(order = 2),
         size = guide_legend(order=1))

dev.off()



# SAVE GO/KEGG results

library(openxlsx)
wb <- createWorkbook()

# Add a worksheet named "Sheet1"
addWorksheet(wb, "Validated_targets_GO")
addWorksheet(wb, "Validated_targets_KEGG")
addWorksheet(wb, "Predicted_targets_GO")
addWorksheet(wb, "Predicted_targets_KEGG")
addWorksheet(wb, "Validated_targets_GO_r0.5")
addWorksheet(wb, "Validated_targets_KEGG_r0.5")
addWorksheet(wb, "Predicted_targets_GO_r0.5")
addWorksheet(wb, "Predicted_targets_KEGG_r0.5")

writeData(wb,sheet='Validated_targets_GO', x= validated.go)
writeData(wb,sheet = "Validated_targets_KEGG", x = validated.kegg)
writeData(wb,sheet='Predicted_targets_GO', x= predicted.go ) 
writeData(wb,sheet='Predicted_targets_KEGG', x= predicted.kegg)
writeData(wb,sheet='Validated_targets_GO_r0.5', x= validated.go.r5)
writeData(wb,sheet = "Validated_targets_KEGG_r0.5", x = validated.kegg.r5)
writeData(wb,sheet='Predicted_targets_GO_r0.5', x= predicted.go.r5 ) 
writeData(wb,sheet='Predicted_targets_KEGG_r0.5', x= predicted.kegg.r5) 

saveWorkbook(wb,
             file= file.path("hsa_miRNA_mRNA_GO_KEGG_output.xlsx"),overwrite = T) 



writeLines(capture.output(sessionInfo()), "sessioninfo_networks.txt")



## Analysis of age regulated miRNAs (up/down)

age_res <- readRDS("C:/Users/alex_/Desktop/Work/2025/Human-microglia-work-submission/final_analysis/human/smRNA/output-files/hsa_miRNA_age_res.Rds")


## Upregulated -----------

up.list <- age_res$microglia_age_res %>% dplyr::filter(logFC > 0 & adj.P.Val < 0.05) %>% dplyr::select(genes)
down.list <- age_res$microglia_age_res %>% dplyr::filter(logFC < 0 & adj.P.Val < 0.05) %>% dplyr::select(genes)

## multimiR -------------------------------------------

microglia.up_multimir_results <- get_multimir(org = 'hsa',
                                          mirna = up.list$genes,
                                          table = 'validated',
                                          summary = TRUE)

microglia.up.multiMir.data <- microglia.up_multimir_results@data

unique(microglia.up.multiMir.data$experiment)

head(pearson.corr)

### subset pearson correlation object for miRNAs enriched in microglia ------------

microglia.upDeList_corr <- pearson.corr[, up.list[[1]]]

# extract p values -----------------------

pearson.corr.p.values <- corr.data$p

pearson.corr.p.values <- t(pearson.corr.p.values) #transpose to get miRNAs as header names


### reformat multimiR results ----------------------------

microglia.up.validated.targets <- as_tibble(microglia.up_multimir_results@summary)

microglia.up.validated.targets <- microglia.up.validated.targets %>% dplyr::select(-mature_mirna_acc, -target_ensembl)

microglia.up.validated.targets <- microglia.up.validated.targets %>% 
  dplyr::rename(
    miRNA = mature_mirna_id,
    Gene = target_symbol)

#### join all corr results for validated - microglia.up_genes_type --------------------

microglia.upDeList_corr_t <- t(microglia.upDeList_corr)

microglia.upDeList_corr_t <- cbind(rownames(microglia.upDeList_corr_t), data.frame(microglia.upDeList_corr_t))

colnames(microglia.upDeList_corr_t)[1] <- c("miRNA") # rename first column

# remove ENS ID from colnames for genes #

colnames(microglia.upDeList_corr_t) <- sub("^[^.]*\\.", "", colnames(microglia.upDeList_corr_t))

#### match validated targets with columns of microglia.up_miRNAs_type_corr_t

microglia.upoutput <- vector("double", length = nrow(microglia.up.validated.targets))
for (i in 1:nrow(microglia.up.validated.targets)) {
  if(length(paste(microglia.upDeList_corr_t[as.character(microglia.up.validated.targets[i,"miRNA"]), as.character(microglia.up.validated.targets[i,"Gene"])])) > 0) { 
    microglia.upoutput[[i]] <- paste(microglia.upDeList_corr_t[as.character(microglia.up.validated.targets[i,"miRNA"]), as.character(microglia.up.validated.targets[i,"Gene"])])
  }
  
}

microglia.upoutput

### append p values ###

microglia.upDe_p_values <- pearson.corr.p.values[, up.list[[1]]]

microglia.upDe_p_values_t <- t(microglia.upDe_p_values)

microglia.upDe_p_values_t <- cbind(rownames(microglia.upDe_p_values_t), data.frame(microglia.upDe_p_values_t))

colnames(microglia.upDe_p_values_t)[1] <- c("miRNA") #rename first column

# remove ENS ID from colnames for genes #

colnames(microglia.upDe_p_values_t) <- sub("^[^.]*\\.", "", colnames(microglia.upDe_p_values_t))

### for loop p values to matrix ###

microglia.upoutput.p.values <- vector("double", length = nrow(microglia.up.validated.targets))
for (i in 1:nrow(microglia.up.validated.targets)) {
  if(length(paste(microglia.upDe_p_values_t[as.character(microglia.up.validated.targets[i,"miRNA"]), as.character(microglia.up.validated.targets[i,"Gene"])])) > 0) { 
    microglia.upoutput.p.values[[i]] <- paste(microglia.upDe_p_values_t[as.character(microglia.up.validated.targets[i,"miRNA"]), as.character(microglia.up.validated.targets[i,"Gene"])])
  }
  
}

microglia.upoutput.p.values

#adjust p values --------------------

microglia.upAdj.p.val <- p.adjust(microglia.upoutput.p.values, method = c("fdr"))

# combine results for final output --------------

microglia.up.validated.targets.pearson <- cbind(microglia.up.validated.targets, as.numeric(microglia.upoutput), microglia.upoutput.p.values, microglia.upAdj.p.val)

microglia.up.validated.targets.pearson <- microglia.up.validated.targets.pearson %>% 
  dplyr::rename(
    r_value = "as.numeric(microglia.upoutput)",
    Pvalue = microglia.upoutput.p.values,
    Adj.P.value = microglia.upAdj.p.val)


microglia.up.validated.targets.pearson <- as_tibble(microglia.up.validated.targets.pearson)


### convert to R^2 ###

microglia.up.validated.targets.pearson.final <- microglia.up.validated.targets.pearson %>% 
  mutate(R_squared = r_value ^ 2)


### filter on adj. p.value (0.05)  ###

#microglia.upfilter <- microglia.up.validated.targets.pearson %>% filter(Adj.P.value < 0.05)

write.csv(microglia.up.validated.targets.pearson.final, file = "Microglia-age-upregulated-miRNA-mRNA-correlations-no_filter.csv", quote = FALSE)

microglia.up.validated.targets.pearson.final <- read.csv("Microglia-age-upregulated-miRNA-mRNA-correlations-no_filter.csv")


# repeat analysis with predicted targets ----------------------------------

## multimiR -------------------------------------------

microglia.up_multimir_results <- get_multimir(org = 'hsa',
                                          mirna = up.list$genes,
                                          table = 'predicted',
                                          summary = TRUE) # run at default cutoff so only top 20% of interactions

microglia.up.multiMir.data <- microglia.up_multimir_results@data


head(pearson.corr)

### subset pearson correlation object for miRNAs enriched in microglia ------------

microglia.upDeList_corr <- pearson.corr[, up.list[[1]]]

# extract p values -----------------------

pearson.corr.p.values <- corr.data$p

pearson.corr.p.values <- t(pearson.corr.p.values) #transpose to get miRNAs as header names


### reformat multimiR results ----------------------------

microglia.up.predicted.targets <- as_tibble(microglia.up_multimir_results@summary)

microglia.up.predicted.targets <- microglia.up.predicted.targets %>% dplyr::select(-mature_mirna_acc, -target_ensembl)

microglia.up.predicted.targets <- microglia.up.predicted.targets %>% 
  dplyr::rename(
    miRNA = mature_mirna_id,
    Gene = target_symbol)

#### join all corr results for predicted - microglia.up_genes_type --------------------

microglia.upDeList_corr_t <- t(microglia.upDeList_corr)

microglia.upDeList_corr_t <- cbind(rownames(microglia.upDeList_corr_t), data.frame(microglia.upDeList_corr_t))

colnames(microglia.upDeList_corr_t)[1] <- c("miRNA") # rename first column

# remove ENS ID from colnames for genes #

colnames(microglia.upDeList_corr_t) <- sub("^[^.]*\\.", "", colnames(microglia.upDeList_corr_t))

#### match predicted targets with columns of microglia.up_miRNAs_type_corr_t

microglia.upoutput <- vector("double", length = nrow(microglia.up.predicted.targets))
for (i in 1:nrow(microglia.up.predicted.targets)) {
  if(length(paste(microglia.upDeList_corr_t[as.character(microglia.up.predicted.targets[i,"miRNA"]), as.character(microglia.up.predicted.targets[i,"Gene"])])) > 0) { 
    microglia.upoutput[[i]] <- paste(microglia.upDeList_corr_t[as.character(microglia.up.predicted.targets[i,"miRNA"]), as.character(microglia.up.predicted.targets[i,"Gene"])])
  }
  
}

microglia.upoutput

### append p values ###

microglia.upDe_p_values <- pearson.corr.p.values[, up.list[[1]]]

microglia.upDe_p_values_t <- t(microglia.upDe_p_values)

microglia.upDe_p_values_t <- cbind(rownames(microglia.upDe_p_values_t), data.frame(microglia.upDe_p_values_t))

colnames(microglia.upDe_p_values_t)[1] <- c("miRNA") #rename first column

# remove ENS ID from colnames for genes #

colnames(microglia.upDe_p_values_t) <- sub("^[^.]*\\.", "", colnames(microglia.upDe_p_values_t))

### for loop p values to matrix ###

microglia.upoutput.p.values <- vector("double", length = nrow(microglia.up.predicted.targets))
for (i in 1:nrow(microglia.up.predicted.targets)) {
  if(length(paste(microglia.upDe_p_values_t[as.character(microglia.up.predicted.targets[i,"miRNA"]), as.character(microglia.up.predicted.targets[i,"Gene"])])) > 0) { 
    microglia.upoutput.p.values[[i]] <- paste(microglia.upDe_p_values_t[as.character(microglia.up.predicted.targets[i,"miRNA"]), as.character(microglia.up.predicted.targets[i,"Gene"])])
  }
  
}

microglia.upoutput.p.values

#adjust p values --------------------

microglia.upAdj.p.val <- p.adjust(microglia.upoutput.p.values, method = c("fdr"))

# combine results for final output --------------

microglia.up.predicted.targets.pearson <- cbind(microglia.up.predicted.targets, as.numeric(microglia.upoutput), microglia.upoutput.p.values, microglia.upAdj.p.val)

microglia.up.predicted.targets.pearson <- microglia.up.predicted.targets.pearson %>% 
  dplyr::rename(
    r_value = "as.numeric(microglia.upoutput)",
    Pvalue = microglia.upoutput.p.values,
    Adj.P.value = microglia.upAdj.p.val)


microglia.up.predicted.targets.pearson <- as_tibble(microglia.up.predicted.targets.pearson)


###  convert to R^2 ###

microglia.up.predicted.targets.pearson.final <- microglia.up.predicted.targets.pearson %>% 
  mutate(R_squared = r_value ^ 2)



write.csv(microglia.up.predicted.targets.pearson.final, file = "Microglia-age-upregulated-mRNA-correlations-no_filter_predicted_networks.csv", quote = FALSE)

microglia.up.predicted.targets.pearson.final <- read.csv("Microglia-age-upregulated-mRNA-correlations-no_filter_predicted_networks.csv")


# KEGG/GO analysis

microglia.upfilter.validated <- microglia.up.validated.targets.pearson.final %>% dplyr::filter(Adj.P.value < 0.05)
microglia.upfilter.predicted <- microglia.up.predicted.targets.pearson.final %>% dplyr::filter(Adj.P.value < 0.05)
microglia.upfilter.validated.r5 <- microglia.up.validated.targets.pearson.final %>% dplyr::filter(Adj.P.value < 0.05 & r_value < -0.5)
microglia.upfilter.predicted.r5 <- microglia.up.predicted.targets.pearson.final %>% dplyr::filter(Adj.P.value < 0.05 & r_value < -0.5)

up.validated.go <- goana(microglia.upfilter.validated$target_entrez, species = "Hs")
up.validated.go$p.adj <- p.adjust(up.validated.go$P.DE, method = "fdr")

up.predicted.go <- goana(microglia.upfilter.predicted$target_entrez, species = "Hs")
up.predicted.go$p.adj <- p.adjust(up.predicted.go $P.DE, method = "fdr")

up.validated.kegg <- kegga(microglia.upfilter.validated$target_entrez, species = "Hs")
up.validated.kegg$p.adj <- p.adjust(up.validated.kegg$P.DE, method = "fdr")

up.predicted.kegg <- kegga(microglia.upfilter.predicted$target_entrez, species = "Hs")
up.predicted.kegg$p.adj <- p.adjust(up.predicted.kegg$P.DE, method = "fdr")

####

up.validated.go.r5 <- goana(microglia.upfilter.validated.r5$target_entrez, species = "Hs")
up.validated.go.r5$p.adj <- p.adjust(up.validated.go.r5$P.DE, method = "fdr")

up.predicted.go.r5 <- goana(microglia.upfilter.predicted.r5$target_entrez, species = "Hs")
up.predicted.go.r5$p.adj <- p.adjust(up.predicted.go.r5$P.DE, method = "fdr")

up.validated.kegg.r5 <- kegga(microglia.upfilter.validated.r5$target_entrez, species = "Hs")
up.validated.kegg.r5$p.adj <- p.adjust(up.validated.kegg.r5$P.DE, method = "fdr")

up.predicted.kegg.r5 <- kegga(microglia.upfilter.predicted.r5$target_entrez, species = "Hs")
up.predicted.kegg.r5$p.adj <- p.adjust(up.predicted.kegg.r5$P.DE, method = "fdr")

# SAVE GO/KEGG results

library(openxlsx)
wb <- createWorkbook()

# Add a worksheet named "Sheet1"
addWorksheet(wb, "Validated_targets_GO")
addWorksheet(wb, "Validated_targets_KEGG")
addWorksheet(wb, "Predicted_targets_GO")
addWorksheet(wb, "Predicted_targets_KEGG")
addWorksheet(wb, "Validated_targets_GO_r0.5")
addWorksheet(wb, "Validated_targets_KEGG_r0.5")
addWorksheet(wb, "Predicted_targets_GO_r0.5")
addWorksheet(wb, "Predicted_targets_KEGG_r0.5")


writeData(wb,sheet='Validated_targets_GO', x= up.validated.go)
writeData(wb,sheet = "Validated_targets_KEGG", x = up.validated.kegg)
writeData(wb,sheet='Predicted_targets_GO', x= up.predicted.go) 
writeData(wb,sheet='Predicted_targets_KEGG', x= up.predicted.kegg ) 
writeData(wb,sheet='Validated_targets_GO_r0.5', x= up.validated.go.r5)
writeData(wb,sheet = "Validated_targets_KEGG_r0.5", x = up.validated.kegg.r5)
writeData(wb,sheet='Predicted_targets_GO_r0.5', x= up.predicted.go.r5 ) 
writeData(wb,sheet='Predicted_targets_KEGG_r0.5', x= up.predicted.kegg.r5) 

saveWorkbook(wb,
             file="hsa_upregulated_miRNA_mRNA_GO_KEGG_output.xlsx",overwrite = T)



## Downregulated -----------


## multimiR -------------------------------------------

microglia.down_multimir_results <- get_multimir(org = 'hsa',
                                              mirna = down.list$genes,
                                              table = 'validated',
                                              summary = TRUE)

microglia.down.multiMir.data <- microglia.down_multimir_results@data

unique(microglia.down.multiMir.data$experiment)

head(pearson.corr)

### subset pearson correlation object for miRNAs enriched in microglia ------------

microglia.downDeList_corr <- pearson.corr[, down.list[[1]]]

# extract p values -----------------------

pearson.corr.p.values <- corr.data$p

pearson.corr.p.values <- t(pearson.corr.p.values) #transpose to get miRNAs as header names


### reformat multimiR results ----------------------------

microglia.down.validated.targets <- as_tibble(microglia.down_multimir_results@summary)

microglia.down.validated.targets <- microglia.down.validated.targets %>% dplyr::select(-mature_mirna_acc, -target_ensembl)

microglia.down.validated.targets <- microglia.down.validated.targets %>% 
  dplyr::rename(
    miRNA = mature_mirna_id,
    Gene = target_symbol)

#### join all corr results for validated - microglia.down_genes_type --------------------

microglia.downDeList_corr_t <- t(microglia.downDeList_corr)

microglia.downDeList_corr_t <- cbind(rownames(microglia.downDeList_corr_t), data.frame(microglia.downDeList_corr_t))

colnames(microglia.downDeList_corr_t)[1] <- c("miRNA") # rename first column

# remove ENS ID from colnames for genes #

colnames(microglia.downDeList_corr_t) <- sub("^[^.]*\\.", "", colnames(microglia.downDeList_corr_t))

#### match validated targets with columns of microglia.down_miRNAs_type_corr_t

microglia.downoutput <- vector("double", length = nrow(microglia.down.validated.targets))
for (i in 1:nrow(microglia.down.validated.targets)) {
  if(length(paste(microglia.downDeList_corr_t[as.character(microglia.down.validated.targets[i,"miRNA"]), as.character(microglia.down.validated.targets[i,"Gene"])])) > 0) { 
    microglia.downoutput[[i]] <- paste(microglia.downDeList_corr_t[as.character(microglia.down.validated.targets[i,"miRNA"]), as.character(microglia.down.validated.targets[i,"Gene"])])
  }
  
}

microglia.downoutput

### append p values ###

microglia.downDe_p_values <- pearson.corr.p.values[, down.list[[1]]]

microglia.downDe_p_values_t <- t(microglia.downDe_p_values)

microglia.downDe_p_values_t <- cbind(rownames(microglia.downDe_p_values_t), data.frame(microglia.downDe_p_values_t))

colnames(microglia.downDe_p_values_t)[1] <- c("miRNA") #rename first column

# remove ENS ID from colnames for genes #

colnames(microglia.downDe_p_values_t) <- sub("^[^.]*\\.", "", colnames(microglia.downDe_p_values_t))

### for loop p values to matrix ###

microglia.downoutput.p.values <- vector("double", length = nrow(microglia.down.validated.targets))
for (i in 1:nrow(microglia.down.validated.targets)) {
  if(length(paste(microglia.downDe_p_values_t[as.character(microglia.down.validated.targets[i,"miRNA"]), as.character(microglia.down.validated.targets[i,"Gene"])])) > 0) { 
    microglia.downoutput.p.values[[i]] <- paste(microglia.downDe_p_values_t[as.character(microglia.down.validated.targets[i,"miRNA"]), as.character(microglia.down.validated.targets[i,"Gene"])])
  }
  
}

microglia.downoutput.p.values

#adjust p values --------------------

microglia.downAdj.p.val <- p.adjust(microglia.downoutput.p.values, method = c("fdr"))

# combine results for final output --------------

microglia.down.validated.targets.pearson <- cbind(microglia.down.validated.targets, as.numeric(microglia.downoutput), microglia.downoutput.p.values, microglia.downAdj.p.val)

microglia.down.validated.targets.pearson <- microglia.down.validated.targets.pearson %>% 
  dplyr::rename(
    r_value = "as.numeric(microglia.downoutput)",
    Pvalue = microglia.downoutput.p.values,
    Adj.P.value = microglia.downAdj.p.val)


microglia.down.validated.targets.pearson <- as_tibble(microglia.down.validated.targets.pearson)


### convert to R^2 ###

microglia.down.validated.targets.pearson.final <- microglia.down.validated.targets.pearson %>% 
  mutate(R_squared = r_value ^ 2)


### filter on adj. p.value (0.05)  ###

#microglia.downfilter <- microglia.down.validated.targets.pearson %>% filter(Adj.P.value < 0.05)

write.csv(microglia.down.validated.targets.pearson.final, file = "Microglia-age-downregulated-miRNA-mRNA-correlations-no_filter.csv", quote = FALSE)

microglia.down.validated.targets.pearson.final <- read.csv("Microglia-age-downregulated-miRNA-mRNA-correlations-no_filter.csv")

# repeat analysis with predicted targets ----------------------------------

## multimiR -------------------------------------------

microglia.down_multimir_results <- get_multimir(org = 'hsa',
                                              mirna = down.list$genes,
                                              table = 'predicted',
                                              summary = TRUE) # run at default cutoff so only top 20% of interactions

microglia.down.multiMir.data <- microglia.down_multimir_results@data


head(pearson.corr)

### subset pearson correlation object for miRNAs enriched in microglia ------------

microglia.downDeList_corr <- pearson.corr[, down.list[[1]]]

# extract p values -----------------------

pearson.corr.p.values <- corr.data$p

pearson.corr.p.values <- t(pearson.corr.p.values) #transpose to get miRNAs as header names


### reformat multimiR results ----------------------------

microglia.down.predicted.targets <- as_tibble(microglia.down_multimir_results@summary)

microglia.down.predicted.targets <- microglia.down.predicted.targets %>% dplyr::select(-mature_mirna_acc, -target_ensembl)

microglia.down.predicted.targets <- microglia.down.predicted.targets %>% 
  dplyr::rename(
    miRNA = mature_mirna_id,
    Gene = target_symbol)

#### join all corr results for predicted - microglia.down_genes_type --------------------

microglia.downDeList_corr_t <- t(microglia.downDeList_corr)

microglia.downDeList_corr_t <- cbind(rownames(microglia.downDeList_corr_t), data.frame(microglia.downDeList_corr_t))

colnames(microglia.downDeList_corr_t)[1] <- c("miRNA") # rename first column

# remove ENS ID from colnames for genes #

colnames(microglia.downDeList_corr_t) <- sub("^[^.]*\\.", "", colnames(microglia.downDeList_corr_t))

#### match predicted targets with columns of microglia.down_miRNAs_type_corr_t

microglia.downoutput <- vector("double", length = nrow(microglia.down.predicted.targets))
for (i in 1:nrow(microglia.down.predicted.targets)) {
  if(length(paste(microglia.downDeList_corr_t[as.character(microglia.down.predicted.targets[i,"miRNA"]), as.character(microglia.down.predicted.targets[i,"Gene"])])) > 0) { 
    microglia.downoutput[[i]] <- paste(microglia.downDeList_corr_t[as.character(microglia.down.predicted.targets[i,"miRNA"]), as.character(microglia.down.predicted.targets[i,"Gene"])])
  }
  
}

microglia.downoutput

### append p values ###

microglia.downDe_p_values <- pearson.corr.p.values[, down.list[[1]]]

microglia.downDe_p_values_t <- t(microglia.downDe_p_values)

microglia.downDe_p_values_t <- cbind(rownames(microglia.downDe_p_values_t), data.frame(microglia.downDe_p_values_t))

colnames(microglia.downDe_p_values_t)[1] <- c("miRNA") #rename first column

# remove ENS ID from colnames for genes #

colnames(microglia.downDe_p_values_t) <- sub("^[^.]*\\.", "", colnames(microglia.downDe_p_values_t))

### for loop p values to matrix ###

microglia.downoutput.p.values <- vector("double", length = nrow(microglia.down.predicted.targets))
for (i in 1:nrow(microglia.down.predicted.targets)) {
  if(length(paste(microglia.downDe_p_values_t[as.character(microglia.down.predicted.targets[i,"miRNA"]), as.character(microglia.down.predicted.targets[i,"Gene"])])) > 0) { 
    microglia.downoutput.p.values[[i]] <- paste(microglia.downDe_p_values_t[as.character(microglia.down.predicted.targets[i,"miRNA"]), as.character(microglia.down.predicted.targets[i,"Gene"])])
  }
  
}

microglia.downoutput.p.values

#adjust p values --------------------

microglia.downAdj.p.val <- p.adjust(microglia.downoutput.p.values, method = c("fdr"))

# combine results for final output --------------

microglia.down.predicted.targets.pearson <- cbind(microglia.down.predicted.targets, as.numeric(microglia.downoutput), microglia.downoutput.p.values, microglia.downAdj.p.val)

microglia.down.predicted.targets.pearson <- microglia.down.predicted.targets.pearson %>% 
  dplyr::rename(
    r_value = "as.numeric(microglia.downoutput)",
    Pvalue = microglia.downoutput.p.values,
    Adj.P.value = microglia.downAdj.p.val)


microglia.down.predicted.targets.pearson <- as_tibble(microglia.down.predicted.targets.pearson)


###  convert to R^2 ###

microglia.down.predicted.targets.pearson.final <- microglia.down.predicted.targets.pearson %>% 
  mutate(R_squared = r_value ^ 2)



write.csv(microglia.down.predicted.targets.pearson.final, file = "Microglia-age-downregulated-mRNA-correlations-no_filter_predicted_networks.csv", quote = FALSE)

microglia.down.predicted.targets.pearson.final <- read.csv("Microglia-age-downregulated-mRNA-correlations-no_filter_predicted_networks.csv")

# KEGG/GO analysis

microglia.downfilter.validated <- microglia.down.validated.targets.pearson.final %>% dplyr::filter(Adj.P.value < 0.05)
microglia.downfilter.predicted <- microglia.down.predicted.targets.pearson.final %>% dplyr::filter(Adj.P.value < 0.05)
microglia.downfilter.validated.r5 <- microglia.down.validated.targets.pearson.final %>% dplyr::filter(Adj.P.value < 0.05 & r_value < -0.5)
microglia.downfilter.predicted.r5 <- microglia.down.predicted.targets.pearson.final %>% dplyr::filter(Adj.P.value < 0.05 & r_value < -0.5)

down.validated.go <- goana(microglia.downfilter.validated$target_entrez, species = "Hs")
down.validated.go$p.adj <- p.adjust(down.validated.go$P.DE, method = "fdr")

down.predicted.go <- goana(microglia.downfilter.predicted$target_entrez, species = "Hs")
down.predicted.go$p.adj <- p.adjust(down.predicted.go $P.DE, method = "fdr")

down.validated.kegg <- kegga(microglia.downfilter.validated$target_entrez, species = "Hs")
down.validated.kegg$p.adj <- p.adjust(down.validated.kegg$P.DE, method = "fdr")

down.predicted.kegg <- kegga(microglia.downfilter.predicted$target_entrez, species = "Hs")
down.predicted.kegg$p.adj <- p.adjust(down.predicted.kegg$P.DE, method = "fdr")

####

down.validated.go.r5 <- goana(microglia.downfilter.validated.r5$target_entrez, species = "Hs")
down.validated.go.r5$p.adj <- p.adjust(down.validated.go.r5$P.DE, method = "fdr")

down.predicted.go.r5 <- goana(microglia.downfilter.predicted.r5$target_entrez, species = "Hs")
down.predicted.go.r5$p.adj <- p.adjust(down.predicted.go.r5$P.DE, method = "fdr")

down.validated.kegg.r5 <- kegga(microglia.downfilter.validated.r5$target_entrez, species = "Hs")
down.validated.kegg.r5$p.adj <- p.adjust(down.validated.kegg.r5$P.DE, method = "fdr")

down.predicted.kegg.r5 <- kegga(microglia.downfilter.predicted.r5$target_entrez, species = "Hs")
down.predicted.kegg.r5$p.adj <- p.adjust(down.predicted.kegg.r5$P.DE, method = "fdr")

# SAVE GO/KEGG results

library(openxlsx)
wb <- createWorkbook()

# Add a worksheet named "Sheet1"
addWorksheet(wb, "Validated_targets_GO")
addWorksheet(wb, "Validated_targets_KEGG")
addWorksheet(wb, "Predicted_targets_GO")
addWorksheet(wb, "Predicted_targets_KEGG")
addWorksheet(wb, "Validated_targets_GO_r0.5")
addWorksheet(wb, "Validated_targets_KEGG_r0.5")
addWorksheet(wb, "Predicted_targets_GO_r0.5")
addWorksheet(wb, "Predicted_targets_KEGG_r0.5")


writeData(wb,sheet='Validated_targets_GO', x= down.validated.go)
writeData(wb,sheet = "Validated_targets_KEGG", x = down.validated.kegg)
writeData(wb,sheet='Predicted_targets_GO', x= down.predicted.go) 
writeData(wb,sheet='Predicted_targets_KEGG', x= down.predicted.kegg ) 
writeData(wb,sheet='Validated_targets_GO_r0.5', x= down.validated.go.r5)
writeData(wb,sheet = "Validated_targets_KEGG_r0.5", x = down.validated.kegg.r5)
writeData(wb,sheet='Predicted_targets_GO_r0.5', x= down.predicted.go.r5 ) 
writeData(wb,sheet='Predicted_targets_KEGG_r0.5', x= down.predicted.kegg.r5) 

saveWorkbook(wb,
             file="hsa_downregulated_miRNA_mRNA_GO_KEGG_output.xlsx",overwrite = T)


writeLines(capture.output(sessionInfo()), "sessioninfo_networks.txt")

