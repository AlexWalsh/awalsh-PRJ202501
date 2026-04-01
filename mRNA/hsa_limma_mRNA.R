## R v4.3.1

# mRNA hsa -------------------------------------------

rm(list=ls()); gc()
#renv::init()

library(tidyverse); library(readxl); library(janitor)
library(glue);      library(scales); library(here) 

if(!require(tidysteps)){
  devtools::install_github('bansell/tidysteps'); 
  library(tidysteps)}

theme_set(theme_minimal())
fix_tidyverse_conflicts()

default_grey <- "#595959"

here()

library(edgeR)
library(limma)



## local paths - define input+output directories ##

setwd("")

Input.path <- ""
Output.path <- ""

# source fns --------------------------------------------------------------

source('input-files/as_tibble_rnm.R')
neglog <- function(x){-log10(x)}

# read in count data -------------------------------------------

dat <- read.delim(file.path(Input.path, "hsa_mRNA_genecounts.txt"), header = TRUE, sep = "", row.names=1, skip = 1, stringsAsFactors = FALSE)

# clean up data --------------------------------------------

dat <- dat |> select(-c(Chr,Start,End,Strand,Length)) # removes columns Chr to Length

colnames(dat) <- c("J01_CD45", "J02_CD45", "J03_CD45", "J04_CD45", "J05_CD45", "J06_CD45", "J07_CD45", "J08_CD45", "J09_CD45", "J10_CD45", "J11_CD45", "J12_CD45", "J13_CD45", "J14_CD45", "J15_CD45", "J16_CD45", "J17_CD45", "A01_CD45", "A02_CD45", "A03_CD45", "A04_CD45", "A05_CD45", "A06_CD45", "J01_Bulk", "J02_Bulk", "J03_Bulk", "J04_Bulk", "J05_Bulk", "J06_Bulk", "J07_Bulk", "J08_Bulk", "J09_Bulk", "J10_Bulk", "J11_Bulk", "J12_Bulk", "J13_Bulk", "J14_Bulk", "J15_Bulk", "J16_Bulk", "J17_Bulk", "A01_Bulk", "A02_Bulk", "A03_Bulk", "A04_Bulk", "A05_Bulk", "A06_Bulk")

dat <- dat[, c(1,24,2,25,3,26,4,27,5,28,6,29,7,30,8,31,9,32,10,33,11,34,12,35,13,36,14,37,15,38,16,39,17,40,18,41,19,42,20,43,21,44,22,45,23,46)] # clunky but reorders columns to be grouped by individual, J01_CD45, J01_Bulk, J02_CD45, J02_Bulk etc.  through to A06


# extract gene ids -----------------------------------------------------------------

rownames(dat) <- gsub("\\.[0-9]*$", "", rownames(dat)) # remove ensembl id version #

dat <- rownames_to_column(as.data.frame(dat), var = "gene_id") |> as_tibble()

genes <- dat[, "gene_id"]

# write raw counts file for cibersortX input -------------------------------------

write_delim(dat, file.path(Output.path, "hsa_rawcounts_cibersort.txt"), delim = "\t")

#---------------------------------------------------------------
# read in metadata -------------------------------------------

jkey <- read_tsv(file.path(Input.path, "hsa_sampleinfo.txt"))
jkey$Age_sqrt <- sqrt(jkey$Age)

# annotate gene IDs --------------------------------------

GRCh38 <- read_rds(file.path(Input.path, "GRCh38_genes.rds"))

annotate <- left_join(genes, GRCh38, by = "gene_id") |> select("gene_id", "gene_name", "description")


# remove samples based on metrics - read depth, low cibersort purity, MDS ---------------
# J05, A02, J14, J03, J16, J02, J01, CD45 --> <50% Microglia cibersort fraction #

outliers_list <- c("J05_CD45","A02_CD45", "J14_CD45", "J03_CD45", "J16_CD45", "J02_CD45", "J01_CD45")

dat <- dat |> select(-all_of(outliers_list))
jkey <- jkey |> dplyr::filter(!jkey$`sample ID` %in% outliers_list)

# create DGE ------------------------------------------

d <- DGEList(counts = dat[,-1],genes=annotate, samples = jkey)

rownames(d) <- d$genes$gene_id

# lib size --------------------------------------------

tibble(d$samples) |> mutate(id=rownames(d$samples)) |> 
  ggplot(aes(x=id,y=lib.size/1e6)) +
  geom_col() + x_angle() + ylab('Million read counts')

# MDS -------------------------------------------------

plotMDS(cpm(d)) # general clustering shows CD45/bulk grouping -- some Bulk outliers #

# filter ----------------------------------------------

keep.exprs <- filterByExpr(d, group=d$samples$Type) # by Type
d_filt <- d[keep.exprs, , keep.lib.sizes=FALSE] 
dim(d)
dim(d_filt) # 23378 geens past expression threshold 



# normalise ----------------------------------------

d_norm <- normLibSizes(d_filt)
d_norm$samples


# mds ----------------------------------------------

plotMDS(d_filt) # clustering by CD45/bulk # normalised plot looks better
plotMDS(d_filt, gene.selection = "common")


plotMDS(d_norm, col=d_norm$samples$group_num,dim=c(1,2))
plotMDS(d_norm, col=d_norm$samples$group_num,dim=c(1,3))
plotMDS(d_norm, col=d_norm$samples$group_num,dim=c(2,3))
plotMDS(d_norm, col=d_norm$samples$group_num,dim=c(3,4))
plotMDS(d_norm, col=d_norm$samples$group_num,dim=c(5,6))

# corrplot ----------------------------------------------------------------

library(corrr)
library(reshape2)



corr_matrix <- corrr::correlate(cpm(d_norm,log=T))
corr_matrix <- corrr::as_matrix(corr_matrix)

reorder_matrix <- function(corr_matrix){
  # Use correlation between variables as distance
  dd <- as.dist((1-corr_matrix)/2)
  hc <- hclust(dd)
  corr_matrix <- corr_matrix[hc$order, hc$order]
}

corr_matrix <- reorder_matrix(corr_matrix)

melted_corr_matrix <- melt(corr_matrix)

ggcorrplot <- ggplot(melted_corr_matrix, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed()

print(ggcorrplot)

# heatmap -----------------------------------------------------------------

library(tidyheatmaps)

#cpm(d_norm,log=T) |> as_tibble() |> 
#   mutate(id=rownames(d_norm)) |> relocate(id) |> 
#   gather(key,value,-1) |> 
#   tidyheatmap(
#   rows=id,
#   columns=key,values=value,
#   cluster_rows = T,
#   cluster_cols = T ,annotation_names_row = F)   

# design -------------------------------------------------------------------

Age <- jkey$Age_sqrt
Type <- as.factor(jkey$Type)
Sex <- as.factor(jkey$Sex)
Region <- gsub(" ", "_", jkey$Region) |> as.factor()
Disease <- as.factor(jkey$Disease)
ID <- as.factor(jkey$ID)
Batch <- as.factor(jkey$Batch)

design <- model.matrix(~0 + Type + Age + Region + Disease + Sex)

# wsva surrogate variable ---------------------------------------------------------------

svs <- wsva(d_norm, design, plot=T,n.sv = 10) # only use svs[,1]

# design + svs ------------------------------------------------------------

design_svs <- bind_cols(design,"sv1"=svs[,1])

# contrasts --------------------------------------------------------------

contrast.type <- makeContrasts(TypeMicroglia - TypeBulk, levels = design_svs)


# vfit -------------------------------------------------------------------

v <- voomWithQualityWeights(d_norm, design_svs, plot=T)


corr_v <- duplicateCorrelation(v, block = v$targets$ID, design = design_svs) # duplicateCorr
corr_v$consensus.correlation # 0.16 -incorporate into vfit

vfit <- lmFit(v, design_svs, block = ID, correlation = corr_v$consensus.correlation)
vfit <- contrasts.fit(vfit, contrast.type)

efit <- eBayes(vfit)

plotSA(efit)

summary(decideTests(efit)) # 8735 up 9170 down

tfit.type <- treat(vfit, lfc=1) 
summary(decideTests(tfit.type)) # 3254 up 4341 down

DGE_results <- topTable(efit, coef=1, n=Inf) |>  as_tibble()

# summary plots ------------------------------------------------------------

topTable(efit,coef=1,n=Inf) |> as_tibble() |> print() |> 
  ggplot(aes(x=logFC,y=neglog(P.Value)))  + 
  geom_point(cex=0.2,col=default_grey) +
  geom_point(data = . %>% dplyr::filter(adj.P.Val<0.01,logFC>0), cex=0.2,col="red") +
  geom_point(data = . %>% dplyr::filter(adj.P.Val<0.01,logFC<0), cex=0.2,col="blue")


topTreat(tfit.type, coef=1, n=Inf) |> as_tibble() |> arrange(logFC,adj.P.Val) |> print() |> 
  ggplot(aes(x=logFC,y=neglog(P.Value)))  + 
  geom_point(cex=0.2,col=default_grey) +
  geom_point(data = . %>% dplyr::filter(adj.P.Val<0.01,logFC>0), cex=0.5,col="red") +
  geom_point(data = . %>% dplyr::filter(adj.P.Val<0.01,logFC<0), cex=0.5,col="blue") +
  annotate("text", x = 3, y = 1, label = "enriched \nin microglia", hjust = 0, vjust = 1, color = "red") +
  annotate("text", x = -3, y = 1, label = "depleted \nin microglia", hjust = 1, vjust = 1, color = "blue") +
  xlim(-10,10) + ylab('-log10(p value)') + xlab('log2(fold-change)')

## ggsave -----------------------------------------------

ggsave(filename = file.path(Output.path, "mrna_treat_volcano.pdf"))

# save -------------------------------------------------

# log cpm ---------------------------------

saveRDS(v$E |> as_tibble_rnm('geneid'),
        file.path(Output.path, "hsa_mRNA_lcpm_resid.Rds"))

## lists -----------------------------------------------

res <- list()
res$d_norm <- d_norm
res$design_svs <- design_svs
res$lcpm_resid <- v$E |> as_tibble() |> mutate(geneID = v$genes$gene_name) |> mutate(ensembl=rownames(v$E)) |> relocate(ensembl) |> relocate(geneID)  

res$tt <- topTable(efit,coef = 1, n=Inf) |> as_tibble()

res$tt_treat <- topTreat(tfit.type,coef = 1, n=Inf) |> as_tibble()

saveRDS(res, file.path(Output.path, "hsa_mRNA_lm_resid.Rds"))

# as xlsx -----------------------------------------------

library(openxlsx)

wb <- createWorkbook()

# Add a worksheet named "Sheet1"
addWorksheet(wb, "log_counts")
addWorksheet(wb, "DGEsummary_notreat")
addWorksheet(wb, "DEGs_FC_treat")

writeData(wb,sheet='log_counts', x=res$lcpm_resid)
writeData(wb,sheet = "DGEsummary_notreat", x =res$tt)
writeData(wb,sheet='DEGs_FC_treat', x= res$tt_treat) #|> filter(adj.P.Val<0.05))


saveWorkbook(wb,file=file.path(Output.path, "hsa_mRNA_output.xlsx"),overwrite = T) 

writeLines(capture.output(sessionInfo()), "sessioninfo.txt")


## MDS plots

# MDS ------------------------

### MDS plots ####

plotMDS(d_norm)
title("MDS")

par(mar = c(5,4,4,8))

Type <- factor(jkey$Type)
#Age <- factor(jkey$Age.group)
#Age <- factor(Age, levels = c("Childhood", "Adolescence", "Adulthood"))
Sex <- factor(jkey$Sex)
Disease <- factor(jkey$Disease)
Region <- factor(jkey$Region)
ID <- factor(jkey$ID)


## Type -  ##

col.type <- c("red", "blue")[Type]
pch.type <- c(16,16)[Type]

pdf(paste(Output.path, "MDS_AgeGroup_Type_Plot.pdf", sep = "/"), width = 10, height = 10)

par(mar = c(5,4,4,8))

plotMDS(d_norm, pch=pch.type, col = col.type, cex=1.5)
legend("topright", legend = levels(as.factor(jkey$Type)), title = "Type", pch = c(15), col = c("red", "blue"), xpd = TRUE, inset = c(-0.15, 0))
title("MDS - Type")

dev.off()

## Sex ##
col.type <- c("red", "blue")[Type]
pch.sex <- c(16,17)[Sex]


pdf(paste(Output.path, "MDS_SexType_Plot.pdf", sep = "/"), width = 10, height = 10)

par(mar = c(5,4,4,8))

plotMDS(d_norm, col=col.type, pch=pch.sex , cex=1.5)
legend("topright", legend = levels(as.factor(jkey$Type)), title = "Type", pch = 15, col = c("red", "blue"), xpd = TRUE, inset = c(-0.15, 0))
legend("bottomright", legend = levels(as.factor(jkey$Sex)), title = "Sex", pch = c(16,17), xpd = TRUE, inset = c(-0.15, 0))
title("MDS - Sex/Type")

dev.off()


## Disease ##

col.type <- c("red", "blue")[Type]
pch.disease <- c(15,16,17,18)[Disease]

pdf(paste(Output.path, "MDS_Disease_Plot.pdf", sep = "/"), width = 10, height = 10)

par(mar = c(5,4,4,8))

plotMDS(d_norm, col=col.type, pch=pch.disease, cex=1.5)
legend("topright", legend = levels(as.factor(jkey$Type)), title = "Type", pch = 15, col = c("red", "blue"), xpd = TRUE, inset = c(-0.15, 0))
legend("bottomright", legend = levels(as.factor(jkey$Disease)), title = "Disease", pch = c(15,16,17,18), xpd = TRUE, inset = c(-0.15, 0))
title("MDS - Disease/Type")

dev.off()

## Region ##
col.type <- c("red", "blue")[Type]
pch.region <- c(15,16,17,18,21)[Region]


pdf(paste(Output.path, "MDS_Region_Plot.pdf", sep = "/"), width = 10, height = 10)

par(mar = c(5,4,4,8))

plotMDS(d_norm, col=col.type , pch=pch.region, cex=1.5)
legend("topright", legend = levels(as.factor(jkey$Type)), title = "Type", pch = 15, col = c("red", "blue"), xpd = TRUE, inset = c(-0.15, 0))
legend("bottomright", legend = levels(as.factor(jkey$Region)), title = "Region", pch = c(15,16,17,18,21), xpd = TRUE, inset = c(-0.15, 0))
title("MDS - Region/Type")

dev.off()

### PCA

## Type -  ##

col.type <- c("red", "blue")[Type]
pch.type <- c(16,16)[Type]

pdf(paste(Output.path, "PCA_AgeGroup_Type_Plot.pdf", sep = "/"), width = 10, height = 10)

par(mar = c(5,4,4,8))

plotMDS(d_norm, pch=pch.type, col = col.type, cex=1.5, gene.selection = "common")
legend("topright", legend = levels(as.factor(jkey$Type)), title = "Type", pch = c(15), col = c("red", "blue"), xpd = TRUE, inset = c(-0.15, 0))
title("PCA - Type")

dev.off()

## Sex ##
col.type <- c("red", "blue")[Type]
pch.sex <- c(16,17)[Sex]


pdf(paste(Output.path, "PCA_SexType_Plot.pdf", sep = "/"), width = 10, height = 10)

par(mar = c(5,4,4,8))

plotMDS(d_norm, col=col.type, pch=pch.sex , cex=1.5, gene.selection = "common")
legend("topright", legend = levels(as.factor(jkey$Type)), title = "Type", pch = 15, col = c("red", "blue"), xpd = TRUE, inset = c(-0.15, 0))
legend("bottomright", legend = levels(as.factor(jkey$Sex)), title = "Sex", pch = c(16,17), xpd = TRUE, inset = c(-0.15, 0))
title("PCA - Sex/Type")

dev.off()


## Disease ##

col.type <- c("red", "blue")[Type]
pch.disease <- c(15,16,17,18)[Disease]

pdf(paste(Output.path, "PCA_Disease_Plot.pdf", sep = "/"), width = 10, height = 10)

par(mar = c(5,4,4,8))

plotMDS(d_norm, col=col.type, pch=pch.disease, cex=1.5, gene.selection = "common")
legend("topright", legend = levels(as.factor(jkey$Type)), title = "Type", pch = 15, col = c("red", "blue"), xpd = TRUE, inset = c(-0.15, 0))
legend("bottomright", legend = levels(as.factor(jkey$Disease)), title = "Disease", pch = c(15,16,17,18), xpd = TRUE, inset = c(-0.15, 0))
title("PCA - Disease/Type")

dev.off()

## Region ##
col.type <- c("red", "blue")[Type]
pch.region <- c(15,16,17,18,21)[Region]


pdf(paste(Output.path, "PCA_Region_Plot.pdf", sep = "/"), width = 10, height = 10)

par(mar = c(5,4,4,8))

plotMDS(d_norm, col=col.type , pch=pch.region, cex=1.5, gene.selection = "common")
legend("topright", legend = levels(as.factor(jkey$Type)), title = "Type", pch = 15, col = c("red", "blue"), xpd = TRUE, inset = c(-0.15, 0))
legend("bottomright", legend = levels(as.factor(jkey$Region)), title = "Region", pch = c(15,16,17,18,21), xpd = TRUE, inset = c(-0.15, 0))
title("PCA - Region/Type")

dev.off()
