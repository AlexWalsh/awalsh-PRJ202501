# R.4.3.1

# mmu, hsa, xla conserved signature ---------------------

rm(list=ls()); gc()

library(dplyr)
library(tidyverse)
library(corrr)
library(readxl)
library(GGally)

setwd("C:/Users/alex_/Desktop/Work/2025/Human-microglia-work-submission/final_analysis/blast")

# read in blast results ----------------------------

blast <- read.delim("blastn/allVall_noxtr.out") |> mutate(across(where(is.character), ~ gsub("_", "-", .))) #convert _ in name to -

# remove rows where species in qseqid and sseqid are the same -  ---------------------

blast_filter <- blast[ sub("-.*", "", blast[["qseqid"]]) != sub("-.*", "", blast[["sseqid"]]), ] # filters out same species comparisons

# filter for e value < 0.01 and annotate species names for query and subject species -------------

blast_filter_e.01_summary <- blast_filter |>
  dplyr::filter(evalue < 0.01) |>
  mutate(q_species = sub("-.*", "", qseqid),
         q_miRNA = sub("^[^-]+-", "", qseqid),
         s_miRNA = sub("^[^-]+-", "", sseqid),
         s_species = sub("-.*", "", sseqid)
         ) |>
  relocate(q_species)


### keep only matching family miRNA names -------------------

blast_hsa_family <- blast_filter_e.01_summary |> dplyr::filter(q_miRNA == s_miRNA) # keeps only matching family names


## conserved miRNAs - keep only miRNAs that are conserved in all 3 species  ------------------------------------

blast_hsa_family |> 
  group_by(q_miRNA) |> 
  dplyr::summarise(species_count = n_distinct(q_species)) |>
  print() #505 miRNA families

conserved_miRNAs <- blast_hsa_family |>
  group_by(q_miRNA) |>
  dplyr::summarise(species_count = n_distinct(q_species)) |>
  dplyr::filter(species_count == 3) |>
  pull(q_miRNA) 

conserved_table <- blast_hsa_family |>
 dplyr:: filter(q_miRNA %in% conserved_miRNAs)


### pairwise comparisons -------------------------------

blast_df_unique <- conserved_table |>
  dplyr::filter(
    (q_species == "hsa" & s_species == "mmu") |
      (q_species == "hsa" & s_species == "xla") |
      (q_species == "mmu" & s_species == "xla")
  ) # this filters for unique instances of each homology pair

## read in DE datasets - limma treat  datasets ------------------------

hsa_RDS <- readRDS("hsa_sRNA_lm_res.rds") 

xla_RDS <- readRDS("xeno_sRNA_lm_res.rds")

hsa_treat <- hsa_RDS$tt_treat

xla_treat <- xla_RDS$tt_treat

mmu_treat <- read.csv("mmu_Walsh2023_CD45vsBulk_treat_DE.csv")


# drop species prefix and keep only miRNAs matching to conserved_miRNAs -----------------

hsa_DE_treat_filter <- hsa_treat |> 
  mutate(miRNA = sub("^[^-]*-", "", genes)) |>
  dplyr::filter(miRNA %in% conserved_miRNAs) |>
  dplyr::select(miRNA, t, adj.P.Val) |>
  rename_with(.fn = ~ paste0("hsa_", .x), .cols = c("t", "adj.P.Val"))

xla_DE_treat_filter <- xla_treat |>
  mutate(miRNA = sub("^[^_]*_", "", genes)) |>
  mutate(miRNA = sub("^[^-]*-", "", miRNA)) |> 
  dplyr::filter(miRNA %in% conserved_miRNAs) |> 
  dplyr::select(miRNA, t, adj.P.Val) |>
  rename_with(.fn = ~ paste0("xla_", .x), .cols = c("t", "adj.P.Val"))

mmu_DE_treat_filter <- mmu_treat |>
  mutate(miRNA = sub("^[^-]*-", "", ID)) |>
  dplyr::filter(miRNA %in% conserved_miRNAs) |>
  dplyr::select(miRNA, t, adj.P.Val) |>
  rename_with(.fn = ~ paste0("mmu_", .x), .cols = c("t", "adj.P.Val"))

# join table ---------------------------

treat_result <- as_tibble(conserved_miRNAs) |> 
  rename(miRNA = value) |> 
  left_join(hsa_DE_treat_filter, by = "miRNA") |>
  left_join(mmu_DE_treat_filter, by = "miRNA") |>
  left_join(xla_DE_treat_filter, by = "miRNA")


## ggpairs and correlation matrix --------------------


ggpairs(treat_result[, c("hsa_t", "mmu_t", "xla_t")])

### filter for positive t values ----------

treat_result_final <- treat_result |> dplyr::filter(hsa_t > 0 & mmu_t > 0 & xla_t > 0 & mmu_adj.P.Val < 0.05 & hsa_adj.P.Val < 0.05 & xla_adj.P.Val < 0.05)

write_tsv(treat_result, file = "hsa_mmu_xla_TREAT_t_corr.tsv")

#---------------------------------------------------
### analysis between hsa vs mmu only ---------------


conserved_miRNAs_mmu_hsa <- blast_hsa_family |>
  dplyr::filter(q_species == "mmu" | q_species == "hsa") |>
  group_by(q_miRNA) |>
  summarise(species_count = n_distinct(q_species)) |>
  dplyr::filter(species_count == 2) |>
  pull(q_miRNA)

conserved_table_mmu_hsa <- blast_hsa_family |>
  dplyr::filter(q_miRNA %in% conserved_miRNAs_mmu_hsa ) |>
  dplyr::filter(q_species != "xla" & s_species != "xla")


#### pairwise comparisons - filter for unique instances of homology -------------

blast_df_unique_mmu_hsa <- conserved_table_mmu_hsa |>
  dplyr::filter(
    (q_species == "hsa" & s_species == "mmu")
  )


# filter DEs by new conserved list mmu vs hsa -------------------------

mmu_DE_treat2 <- mmu_treat |>
  mutate(miRNA = sub("^[^-]*-", "", ID)) |>
  dplyr::filter(miRNA %in% conserved_miRNAs_mmu_hsa) |>
  dplyr::select(miRNA, t, adj.P.Val) |>
  rename_with(.fn = ~ paste0("mmu_", .x), .cols = c("t", "adj.P.Val"))

hsa_DE_treat2 <- hsa_treat |> 
  mutate(miRNA = sub("^[^-]*-", "", genes)) |>
  dplyr::filter(miRNA %in% conserved_miRNAs_mmu_hsa) |>
  dplyr::select(miRNA, t, adj.P.Val) |>
  rename_with(.fn = ~ paste0("hsa_", .x), .cols = c("t", "adj.P.Val"))


conserved_blast_core_mmu_hsa <- as_tibble(conserved_miRNAs_mmu_hsa) |> rename(miRNA = value)

treat_result_mmu_hsa <- conserved_blast_core_mmu_hsa |>
  left_join(hsa_DE_treat2, by = "miRNA") |>
  left_join(mmu_DE_treat2, by = "miRNA") 


ggpairs(treat_result_mmu_hsa[, c("hsa_t", "mmu_t")])


treat_result_final_mmu_hsa <- treat_result_mmu_hsa |> dplyr::filter(hsa_t > 0 & mmu_t > 0 & mmu_adj.P.Val < 0.05 & hsa_adj.P.Val < 0.05)

write_tsv(treat_result_mmu_hsa, file = "hsa_mmu_TREAT_t_corr.tsv")

#-----------------------------------------------------------------------------
### Repeat analysis with NO TREAT results ###

hsa_no_treat <- hsa_RDS$tt

xla_no_treat <- xla_RDS$tt

mmu_no_treat <- read.csv("mmu_Walsh2023_CD45vsBulk_notreat_DE.csv")


# drop species prefix and keep only miRNAs matching to conserved_miRNAs -----------------

hsa_DE_no_treat_filter <- hsa_no_treat |> 
  mutate(miRNA = sub("^[^-]*-", "", genes)) |>
  dplyr::filter(miRNA %in% conserved_miRNAs) |>
  dplyr::select(miRNA, t, adj.P.Val) |>
  rename_with(.fn = ~ paste0("hsa_", .x), .cols = c("t", "adj.P.Val"))

xla_DE_no_treat_filter <- xla_no_treat |>
  mutate(miRNA = sub("^[^_]*_", "", genes)) |>
  mutate(miRNA = sub("^[^-]*-", "", miRNA)) |> 
  dplyr::filter(miRNA %in% conserved_miRNAs) |> 
  dplyr::select(miRNA, t, adj.P.Val) |>
  rename_with(.fn = ~ paste0("xla_", .x), .cols = c("t", "adj.P.Val"))

mmu_DE_no_treat_filter <- mmu_no_treat |>
  mutate(miRNA = sub("^[^-]*-", "", ID)) |>
  dplyr::filter(miRNA %in% conserved_miRNAs) |>
  dplyr::select(miRNA, t, adj.P.Val) |>
  rename_with(.fn = ~ paste0("mmu_", .x), .cols = c("t", "adj.P.Val"))

# join table ---------------------------

no_treat_result <- as_tibble(conserved_miRNAs) |> 
  rename(miRNA = value) |>
  left_join(hsa_DE_no_treat_filter, by = "miRNA") |>
  left_join(mmu_DE_no_treat_filter, by = "miRNA") |>
  left_join(xla_DE_no_treat_filter, by = "miRNA")

## ggpairs and correlation matrix --------------------


ggpairs(no_treat_result[, c("hsa_t", "mmu_t", "xla_t")])


### filter for positive t values ----------

no_treat_result_final <- no_treat_result |> dplyr::filter(hsa_t > 0 & mmu_t > 0 & xla_t > 0 & mmu_adj.P.Val < 0.05 & hsa_adj.P.Val < 0.05 & xla_adj.P.Val < 0.05)

write_tsv(no_treat_result, file = "hsa_mmu_xla_t_corr.tsv")


### analysis between hsa vs mmu only ---------------


conserved_miRNAs_mmu_hsa <- blast_hsa_family |>
  dplyr::filter(q_species == "mmu" | q_species == "hsa") |>
  group_by(q_miRNA) |>
  summarise(species_count = n_distinct(q_species)) |>
  dplyr::filter(species_count == 2) |>
  pull(q_miRNA)

conserved_table_mmu_hsa <- blast_hsa_family |>
  dplyr::filter(q_miRNA %in% conserved_miRNAs_mmu_hsa ) |>
  dplyr::filter(q_species != "xla" & s_species != "xla")


## pairwise comparison - filter for unique instances of homology -----------------

blast_df_unique_mmu_hsa <- conserved_table_mmu_hsa |>
  dplyr::filter(
    (q_species == "hsa" & s_species == "mmu")
  ) 


# filter DEs by new conserved list mmu vs hsa -------------------------

mmu_DE_no_treat2 <- mmu_no_treat |>
  mutate(miRNA = sub("^[^-]*-", "", ID)) |>
  dplyr::filter(miRNA %in% conserved_miRNAs_mmu_hsa) |>
  dplyr::select(miRNA, t, adj.P.Val) |>
  rename_with(.fn = ~ paste0("mmu_", .x), .cols = c("t", "adj.P.Val"))

hsa_DE_no_treat2 <- hsa_no_treat |> 
  mutate(miRNA = sub("^[^-]*-", "", genes)) |>
  dplyr::filter(miRNA %in% conserved_miRNAs_mmu_hsa) |>
  dplyr::select(miRNA, t, adj.P.Val) |>
  rename_with(.fn = ~ paste0("hsa_", .x), .cols = c("t", "adj.P.Val"))

no_treat_result_mmu_hsa <- as_tibble(conserved_miRNAs_mmu_hsa) |> 
  rename(miRNA = value) |>
  left_join(hsa_DE_no_treat2, by = "miRNA") |>
  left_join(mmu_DE_no_treat2, by = "miRNA") 


ggpairs(no_treat_result_mmu_hsa[, c("hsa_t", "mmu_t")])


no_treat_result_final_mmu_hsa <- no_treat_result_mmu_hsa |> dplyr::filter(hsa_t > 0 & mmu_t > 0 & mmu_adj.P.Val < 0.05 & hsa_adj.P.Val < 0.05)

write_tsv(no_treat_result_mmu_hsa, file = "hsa_mmu_t_corr.tsv")

writeLines(capture.output(sessionInfo()), "sessioninfo.txt")


