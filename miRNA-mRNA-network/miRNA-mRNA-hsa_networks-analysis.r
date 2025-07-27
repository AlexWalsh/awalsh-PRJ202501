# R.4.3.1
#### network analysis, miRNA, mRNA #####


library(dplyr)
library(tidyverse)
library(psych)
library(ggplot2)
library(multiMiR)

#### load in rds objects ###

corr.data <- readRDS("hsa_miRNA_mRNA_microglia.corr.rds")

pearson.corr <- t(corr.data$r) #transpose to get miRNAs as header names

### read in enriched microglia miRNAs ---------------------------

lm_res <- readRDS("hsa_sRNA_lm_res.rds")

## subset for significantly enriched miRNAs in microglia - logFC > 1 and adj.p.val <0.05 -----------

tt_treat <- lm_res$tt_treat |> filter(logFC > 1 & adj.P.Val < 0.05)

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
  rename(
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
  rename(
    r_value = "as.numeric(CD45Bulkoutput)",
    Pvalue = CD45Bulkoutput.p.values,
    Adj.P.value = CD45BulkAdj.p.val)


CD45Bulk.validated.targets.pearson <- as_tibble(CD45Bulk.validated.targets.pearson)


### remove r values of zero and then convert to R^2 ###

CD45Bulk.validated.targets.pearson.final <- CD45Bulk.validated.targets.pearson %>% 
  filter(r_value != 0) %>%
  mutate(R_squared = r_value ^ 2)


### filter on adj. p.value (0.05)  ###

#CD45Bulkfilter <- CD45Bulk.validated.targets.pearson %>% filter(Adj.P.value < 0.05)

write.csv(CD45Bulk.validated.targets.pearson.final, file = "Microglia-enriched-miRNA-mRNA-correlations-no_filter.csv", quote = FALSE)





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
  rename(
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
  rename(
    r_value = "as.numeric(CD45Bulkoutput)",
    Pvalue = CD45Bulkoutput.p.values,
    Adj.P.value = CD45BulkAdj.p.val)


CD45Bulk.predicted.targets.pearson <- as_tibble(CD45Bulk.predicted.targets.pearson)


### remove r values of zero and then convert to R^2 ###

CD45Bulk.predicted.targets.pearson.final <- CD45Bulk.predicted.targets.pearson %>% 
  filter(r_value != 0) %>%
  mutate(R_squared = r_value ^ 2)


### filter on adj. p.value (0.05)  ###

#CD45Bulkfilter <- CD45Bulk.predicted.targets.pearson %>% filter(Adj.P.value < 0.05)

write.csv(CD45Bulk.predicted.targets.pearson.final, file = "Microglia-enriched-miRNA-mRNA-correlations-no_filter_predicted_networks.csv", quote = FALSE)


writeLines(capture.output(sessionInfo()), "sessioninfo_networks.txt")

