## R v4.3.1 

# sRNA - hsa --------------------------------------------------------------------

## Aim: preprocess mature miRNA counts to deduplicate those with highly similar counts (within 10%)  and retain those with higher variance, together with their  precursor mi_rnas.


##

rm(list=ls()); gc()
#renv::init()

library(tidyverse); library(readxl); library(janitor)
library(glue);      library(scales); library(here) ; library(tidysteps)

if(!require(tidyverse)){
  devtools::install_github('bansell/tidyverse'); 
  library(tidyverse)}

theme_set(theme_minimal())
fix_tidyverse_conflicts()

default_grey <- "#595959"


library(edgeR)
library(strex)


#for local work#
# Define paths #

setwd("")
Input.path <- ""
Output.path <- ""

# source fns --------------------------------------------------------------

source('input-files/as_tibble_rnm.R')
neglog <- function(x){-log10(x)}

# read data ---------------------------------------------------------------

colnames <- read_tsv(file.path(Input.path, "colnames.txt"), col_names = FALSE)

dat <- readRDS(paste(Input.path, "hsa_known_miRNA_counts.rds", sep = "/")) |> select(-c(total_read_count))
colnames(dat) <- colnames

### add up reads for read_count ###

dat <- dat |> mutate(read_count = rowSums(across(3:48))) |>
  relocate(read_count, .after = mi_rna)


# check redundant ---------------------------------------------------------


dat |> sort_prop(mi_rna) |> count(n=1)  #1793 unique
dat |> sort_prop(mi_rna) |> count(n>1) #181


#subset to duplicate/multiplicates------------------------------

dups <- dat |> sort_prop(mi_rna) |> dplyr::filter(n>1)  #181 duplicate/multiplicates.

dups |>  
  left_join(dat,by='mi_rna') |> 
  select(mi_rna,read_count) |> distinct() |> sort_prop(mi_rna) #181


# compare duplicate mature counts -----------------------------------------

dups |>
  left_join(dat,by='mi_rna') |> 
  select(mi_rna,read_count) |> distinct() |> group_by(mi_rna) |>
  mutate(rnm=paste0('dup',row_number())) |> spread(rnm,read_count) |> 
  #mutate(across(contains('dup'),log10)) |> 
  #ggplot(aes(x=dup1,y=dup2)) + 
  ggplot(aes(x=.panel_x,y=.panel_y)) +
  geom_point() + 
  scale_x_log10() + scale_y_log10() +
  ggforce::facet_matrix(vars(contains('dup'))) +
  #geom_smooth(aes(x=x,y=x),method = "lm", formula = y ~ x, se = FALSE, color = "red")
  geom_abline() #dont work

# identify miRNAs with greater than 10% variance between duplicates

var9 <- dups |>
  left_join(dat, by='mi_rna') |> 
  select(mi_rna,read_count) |> distinct() |> 
  group_by(mi_rna) |> 
  mutate(var=read_count/(sum(read_count)/n())) |> 
  dplyr::filter(var>1.1 | var<0.9)
  
var9 |> sort_prop(mi_rna) |> print() #38 unique miRNAs


#pull in precursors
dat |> dplyr::filter(mi_rna %in% var9$mi_rna) |> 
  group_by(mi_rna) |> mutate(var=read_count/(sum(read_count)/n())) |> 
  relocate(var,.after=read_count) |> 
  arrange(mi_rna,var)
  


#Results
#181 duplicated mature miRNAs (right panel)


#82 dup entries have variance > 10% between duplicates for total reads mapped.
#38 unique miRNAs with variance >10% between duplicates

#Approach
#Deduplicate low-variance genes, but retain those varying > 10%, together with precursor annotation.


# label high variance miRNAs ----------------------------------------------

var9_label <- dat |> dplyr::filter(mi_rna %in% var9$mi_rna) |>
  #append precursor to mature miRNA for unique label:
  mutate(mi_rna=paste0(mi_rna,'_',precursor)) 



# dedup low variance miRNAs -------------------------------------------------------------

### slice dups ----------------------------------------------------------------


#sliced i.e. selecting first row per group - for each low variance miRNA keep first occurence

dat_dup_slice <- dat |> dplyr::filter(!mi_rna %in% var9$mi_rna) |> 
  select(-c(precursor,read_count)) |> group_by(mi_rna) |> 
  slice_head(n=1) |> ungroup()


### sum dups ----------------------------------------------------------------

dat_dup_sum <- dat |> dplyr::filter(!mi_rna %in% var9$mi_rna) |> 
  select(-c(1,3))|> print() |> group_by(mi_rna) |> 
  reframe(across(everything() , sum)) |> print()

### median dups ----------------------------------------------------------------


dat_dup_med <- dat |> dplyr::filter(!mi_rna %in% var9$mi_rna) |> 
  select(-c(precursor,read_count))|> print() |> group_by(mi_rna) |>  
  reframe(across(everything() , median)) |> print()

### mean dups ----------------------------------------------------------------


dat_dup_mean <- dat |> dplyr::filter(!mi_rna %in% var9$mi_rna) |> 
  select(-c(precursor,read_count))|> print() |> group_by(mi_rna) |>  
  reframe(across(everything() , mean)) |> print()




# bind rows ---------------------------------------------------------------

#default to use the sampled (slice) 
dat_dedup <- bind_rows(
  dat_dup_slice,
  var9_label) |> 
  distinct() |>
  select(-c(precursor,read_count))




# save --------------------------------------------------------------------

saveRDS(dat_known_dedup, 'counts_dedup.Rds')












dat_unique_knownonly <- unique(dat_knownonly$mi_rna)

dat_unique1 <- dat |> dplyr::filter(!mi_rna %in% dat_unique_knownonly)
