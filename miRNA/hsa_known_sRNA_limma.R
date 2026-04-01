## R v4.3.1 

# hsa sRNA --------------------------------------------------------------------

## Aim: use limma to perform DE analysis on hsa miRNAs  ##


rm(list=ls()); gc()
#renv::init()

library(tidyverse); library(readxl); library(janitor)
library(glue);      library(scales); library(here); library(tidyext)

if(!require(tidysteps)){
  devtools::install_github('bansell/tidysteps'); 
  library(tidysteps)}

theme_set(theme_minimal())
fix_tidyverse_conflicts()

default_grey <- "#595959"

library(edgeR)
library(strex)
library(limma)

# source fns --------------------------------------------------------------

setwd("")

source('input-files/as_tibble_rnm.R')
neglog <- function(x){-log10(x)}

#for local work#----------------------------------------

Input.path <- "C:/Users/alex_/Desktop/Work/2025/Human-microglia-work-submission/final_analysis/human/smRNA/input-files"
Output.path <- "C:/Users/alex_/Desktop/Work/2025/Human-microglia-work-submission/final_analysis/human/smRNA/output-files"

####

# metadata ---------------------------------

jkey <- read_tsv(file.path(Input.path, "hsa_sampleinfo.txt"))
jkey$Age_sqrt <- sqrt(jkey$Age)


# append cibersortX results -------------------------

cibersort_fractions <- read_tsv("C:/Users/alex_/Desktop/Work/2025/Human-microglia-work-submission/final_analysis/human/cibersort/CIBERSORTx_Adjusted.txt")

jkey <- bind_cols(jkey, tibble(MG_fraction = cibersort_fractions$MG))



# read in dedeuplicated counts file from deduplicate.R ------------------------

dat_dedup <- readRDS('counts_dedup.Rds')

# remove samples based on metrics - read depth, low cibersort purity, MDS #
# J06 bulk --> less than 1 million mapped reads - remove #
# J05, A02, J14, J03, J16, J02, J01, J05 CD45 --> <50% Microglia cibersort fraction #

outliers_list <- c("J06_Bulk", "J05_CD45", "A02_CD45", "J14_CD45", "J03_CD45", "J16_CD45", "J02_CD45", "J01_CD45")

dat_dedup <- dat_dedup |> dplyr::select(-all_of(outliers_list))
jkey <- jkey |> dplyr::filter(!jkey$`sample ID` %in% outliers_list)

# dgelist --------------------------------

dge <- DGEList(counts = dat_dedup %>% dplyr::select(-mi_rna),
               genes=dat_dedup$mi_rna,
               samples = jkey)
rownames(dge) <- dat_dedup$mi_rna


## limma ------------------------------

boxplot(cpm(dge,log=T))

dge$samples |> as_tibble(rownames='sample') |> 
  ggplot(aes(x=sample,y=lib.size)) + geom_col() + x_angle()

##filter ------------------------------

keep.exprs <- filterByExpr(dge, group=dge$samples$Type) # by group type - CD45/Bulk
d_filt <- dge[keep.exprs, , keep.lib.sizes=FALSE] 
dim(dge) #1881
dim(d_filt) #838 miRNAs pass expression threshold

d_norm <- normLibSizes(d_filt)

d_norm$samples


#mds -----------------------------------

plotMDS(d_norm, col=d_norm$samples$group_num,dim=c(1,2))
plotMDS(d_norm, col=d_norm$samples$group_num,dim=c(1,3))
plotMDS(d_norm, col=d_norm$samples$group_num,dim=c(2,3))
plotMDS(d_norm, col=d_norm$samples$group_num,dim=c(3,4))
plotMDS(d_norm, col=d_norm$samples$group_num,dim=c(5,6))


boxplot(cpm(dge, log=T))
boxplot(cpm(d_norm, log=T))

plotMDS(d_filt) 
plotMDS(d_filt, gene.selection = "common")


# design ------------------------------



Age <- d_norm$samples$Age_sqrt
Type <- as.factor(d_norm$samples$Type)
Sex <- as.factor(d_norm$samples$Sex)
Region <- gsub(" ", "_", d_norm$samples$Region) |> as.factor()
Disease <- as.factor(d_norm$samples$Disease)
ID <- as.factor(d_norm$samples$ID)
Batch <- as.factor(d_norm$samples$Batch)
Mg_fraction <- d_norm$samples$MG_fraction

design <- model.matrix(~0 + Type + Age + Region + Disease + Sex) 


# sva - surrogate variable ------------------------------------

svs <- wsva(d_norm, design, plot=T,n.sv = 10)


# design + svs ------------------------------------------------------------

design_svs <- bind_cols(design,"sv1"=svs[,1]) #use only svs[,1] for design

# contrasts --------------------------------------------------------------

contrast.type <- makeContrasts(TypeMicroglia - TypeBulk, levels = design_svs)

# vfit -------------------------------------------------------------------

v <- voomWithQualityWeights(d_norm, design_svs, plot=T)

corr_v <- duplicateCorrelation(v, block = v$targets$ID, design = design_svs)
corr_v$consensus.correlation # 0.023 --> can't include in design because it is not full rank, however it is minimal and will not have a large effect



vfit <- lmFit(v, design_svs)

vfit <- contrasts.fit(vfit, contrast.type)
efit <- eBayes(vfit)
plotSA(efit) 

summary(decideTests(efit)) # 170 up and 342 down

DGE_results <- topTable(efit, coef=1, n=Inf) |>  as_tibble()


# toptreat - logFC > 1 | logFC < -1 ---------------------------------------------------------------

tfit <- treat(vfit, lfc=1) 
summary(decideTests(tfit)) #97 up 245 down


ttreat <- topTreat(tfit, coef=1, n=Inf) |>  as_tibble() |>
  as_tibble() |> arrange(logFC,adj.P.Val) |> print()



# plot candidates - eg. 511-3p -----------------------------------------------------------------

mygene <- 'hsa-miR-511-3p'

geneplots <- as.data.frame(v$E) |> 
  rownames_to_column(var = "Geneid") |> 
  as_tibble()

geneplots |>
  dplyr::filter(Geneid == mygene) |> 
  gather(key,value,-Geneid) |> 
  mutate(id= str_before_first(key, '_'))|> 
  mutate(id=str_remove(key,'[0-9]')) |> 
  ggplot(aes(x=id,y=value)) + geom_point(point_size = 1) +
  ggtitle(mygene)

# plot by group -------

gene_long <- geneplots %>%
  dplyr::filter(Geneid == mygene) %>%
  pivot_longer(
    cols = -Geneid,
    names_to = "Sample",
    values_to = "value"
  ) %>%
  mutate(Type = ifelse(str_detect(Sample, "CD45"), "CD45", "Bulk"))

ggplot(gene_long, aes(x = Type, y = value, fill = Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  ggtitle(mygene) +
  theme_minimal()


## volcano ---------------------------------------------------------------

#efit
topTable(efit, coef=1, n=Inf) |>  as_tibble() |> 
  arrange(logFC,adj.P.Val) |> print() |> 
  ggplot(aes(x=logFC,y=neglog(P.Value)))  + 
  geom_point(cex=0.2) +
  geom_point(data = . %>% dplyr::filter(adj.P.Val<0.1,logFC>0), cex=0.5,col="red") +
  geom_point(data = . %>% dplyr::filter(adj.P.Val<0.1,logFC<0), cex=0.5,col="blue") +
  annotate("text", x = 3, y = 1, label = "enriched \nin microglia", hjust = 0, vjust = 1, color = "red") +
  annotate("text", x = -3, y = 1, label = "depleted \nin microglia", hjust = 1, vjust = 1, color = "blue") +
  xlim(-6,6) + 
  ylab('-log10(p value)') + xlab('log2(fold-change)')

#tfit

topTreat(tfit, coef=1, n=Inf) |>  as_tibble() |> 
  arrange(adj.P.Val) |> print() |> 
  ggplot(aes(x=logFC,y=neglog(P.Value)))  + 
  geom_point(cex=0.2) +
  geom_point(data = . %>% dplyr::filter(adj.P.Val<0.1,logFC>0), cex=0.5,col="red") +
  geom_point(data = . %>% dplyr::filter(adj.P.Val<0.1,logFC<0), cex=0.5,col="blue") +
  annotate("text", x = 3, y = 1, label = "enriched \nin microglia", hjust = 0, vjust = 1, color = "red") +
  annotate("text", x = -3, y = 1, label = "depleted \nin microglia", hjust = 1, vjust = 1, color = "blue") +
  xlim(-6,6) + 
  ylab('-log10(p value)') + xlab('log2(fold-change)')


ggsave(filename = file.path(Output.path,'hsa_srna_treat_volcano.pdf'), width=7, height=5)


# save --------------------------------------------------------------------


## save v$E voom counts - logCPM ----------------------------------------------------

saveRDS(v$E |> as_tibble_rnm('geneid'),
        file.path(Output.path, "hsa_smRNA_lcpm_resid.Rds"))


### as list -----------------------------------------------------------------

res <- list()
res$d_norm <- d_norm
res$design_svs <- design_svs
res$lcpm_resid <- v$E |> as_tibble()|> mutate(geneID=rownames(v$E)) |> relocate(geneID)

res$tt <- topTable(efit,coef = 1, n=Inf) |> as_tibble()

res$tt_treat <- topTreat(tfit,coef = 1, n=Inf) |> as_tibble()

saveRDS(res, file.path(Output.path, "hsa_sRNA_lm_res.Rds"))

### as xlsx -----------------------------------------------------------------


library(openxlsx)
wb <- createWorkbook()

# Add a worksheet named "Sheet1"
addWorksheet(wb, "log_counts")
addWorksheet(wb, "DGEsummary_notreat")
addWorksheet(wb, "DGEsummary_FC_treat")

writeData(wb,sheet='log_counts', x=res$lcpm_resid)
writeData(wb,sheet = "DGEsummary_notreat", x =res$tt)
writeData(wb,sheet='DGEsummary_FC_treat', x= res$tt_treat ) #|> filter(adj.P.Val<0.05))

saveWorkbook(wb,
             file= file.path(Output.path, "hsa_smRNA_output.xlsx"),overwrite = T) 





## plots for paper ------------------------------

#Enhanced volano -------------------------

library(EnhancedVolcano)

CD45vsBulklabels <- gsub("hsa-miR-", "", ttreat$genes)
CD45vsBulk.DE.graph <- cbind(CD45vsBulklabels, ttreat)
colnames(CD45vsBulk.DE.graph)[1] <- "label"

pdf(paste(Output.path, "CD45vsBulk.DE.volcano.pdf", sep = "/"), width = 12, height = 10)

p10 <- EnhancedVolcano(CD45vsBulk.DE.graph, lab = rownames(CD45vsBulk.DE.graph$ID), x= "logFC", y = "P.Value",
                       selectLab = c("31-3p", "511-3p", "511-5p", "4645-3p", "205-5p", "5193", "653-3p", "31-5p", "3074-3p", "223-3p", "489-3p", "500b-3p", "223-5p", "142-3p", "424-5p", "424-5p", "146a-3p", "3690", "450a-1-3p", "424-3p"),
                       FCcutoff = 1,
                       labSize = 4,
                       pCutoff = 0.05,
                       pCutoffCol = "adj.P.Val",
                       xlim = c(-8, 8),
                       vline = c(5,-5),
                       title = 'CD45+ Microglia vs Bulk miRNA expression',
                       legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                      'p-value & Log (base 2) FC'),
                       legendPosition = "top",
                       legendLabSize = 13,
                       gridlines.major = FALSE, gridlines.minor = FALSE,
                       pointSize = 3.0,
                       drawConnectors = TRUE,
                       arrowheads = FALSE,
                       typeConnectors = "open",
                       boxedLabels = TRUE,
                       directionConnectors = "both",
                       lengthConnectors = 20,
                       maxoverlapsConnectors = 15,
                       widthConnectors =  0.8,
                       colAlpha = 1)

p10 + ggplot2::coord_cartesian(xlim=c(-8, 8)) + ggplot2::scale_x_continuous(breaks=seq(-8,8, 2)) 

dev.off()

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


###############################################
##### interaction Type + Age - all samples - December 2025 - suggested by Sarrabeth ####

design.age.type.int <- model.matrix(~Type:Age + Sex + Region + Disease)

v.age.type <- voomWithQualityWeights(d_norm, design.age.type.int, plot = TRUE)

corr_v_age_type <- duplicateCorrelation(v.age.type, block = v.age.type$targets$ID, design = design.age.type.int)
corr_v_age_type$consensus.correlation # 0.04

fit.age.type <- lmFit(v.age.type, block = v.age.type$targets$ID, correlation = corr_v_age_type$consensus.correlation)

fit.age.type <- eBayes(fit.age.type)

Toptable.age.Microglia <- topTreat(fit.age.type, coef= 11, sort.by ="p", number = "inf")
Toptable.age.Bulk <- topTreat(fit.age.type, coef= 10, sort.by ="p", number = "inf")

Toptable.age.Microglia.DE <- topTreat(fit.age.type, coef= 11, sort.by ="p", number = "inf", p= 0.05)
Toptable.age.Bulk.DE <- topTreat(fit.age.type, coef= 10, sort.by ="p", number = "inf", p = 0.05)


microglia_age_res <- topTreat(fit.age.type, coef = 11, n=Inf) |> as_tibble()
bulk_age_res <- topTreat(fit.age.type, coef = 10, n=Inf) |> as_tibble()


age_res_rds <- list()

age_res_rds$microglia_age_res <- topTreat(fit.age.type, coef = 11, n=Inf) |> as_tibble()
age_res_rds$bulk_age_res <- topTreat(fit.age.type, coef = 10, n=Inf) |> as_tibble()

saveRDS(age_res_rds, file.path(Output.path, "hsa_miRNA_age_res.Rds"))


## save workbook


wb_age <- createWorkbook()

# Add a worksheet named "Sheet1"
addWorksheet(wb_age, "microglia_age_coeff")
addWorksheet(wb_age, "bulk_age_coeff")

writeData(wb_age,sheet= "microglia_age_coeff", x=age_res_rds$microglia_age_res)
writeData(wb_age,sheet = "bulk_age_coeff", x = age_res_rds$bulk_age_res)

saveWorkbook(wb_age,
             file= file.path(Output.path, "hsa_smRNA_age_analysis_output.xlsx"),overwrite = T) 

### filter microglia age.DEs by those that are enriched ####

hsa.microglia.age.up <- Toptable.age.Microglia.DE %>% dplyr::filter(logFC > 0 & adj.P.Val < 0.05)

hsa.microglia.age.down <- Toptable.age.Microglia.DE %>% dplyr::filter(logFC < 0 & adj.P.Val < 0.05)

write(hsa.microglia.age.up$genes, file = file.path(Output.path, "hsa_microglia_age_up.txt"))
write(hsa.microglia.age.down$genes, file = file.path(Output.path, "hsa_microglia_age_down.txt"))

stripchart.age.plot <- subset(v.age.type$E, select = grep("CD45", colnames(v.age.type$E)))
indices <- grep("CD45", colnames(v.age.type$E))
microglia.sample.age <- jkey$Age[c(indices)]

stripchart(stripchart.age.plot["hsa-miR-9-5p",]~microglia.sample.age, vertical=TRUE, las =2, cex.axis = 0.5, pch = 16, cex = 1.3, ylab = "Normalised log2 expression", xlab = "Age", main ="1") 
lm(stripchart.age.plot["hsa-miR-9-5p", ]~microglia.sample.age)


stripchart(v.age.type$E["hsa-miR-9-5p",]~Age, vertical=TRUE, las =2, cex.axis = 0.5, pch = 16, cex = 1.3, ylab = "Normalised log2 expression", xlab = "Age", main ="1")
lm(v.age.type$E["hsa-miR-9-5p",]~Age)


stripchart.age.ggplot <- as.data.frame(t(rbind(stripchart.age.plot, microglia.sample.age)))

lm(stripchart.age.plot["hsa-miR-511-3p", ]~microglia.sample.age)
ggplot(stripchart.age.ggplot, aes(x = stripchart.age.ggplot[,'microglia.sample.age'], stripchart.age.ggplot[,'hsa-miR-511-3p'])) + geom_point(colour = 'gray32', fill = 'gray20', stroke = 1.5, alpha = 1) + geom_point(alpha=0.4, size=3) + theme_classic(base_size = 17) + geom_abline(intercept = 2.70006, slope = 0.01978, size = 1) + ggtitle("hsa-miR-511-3p") + xlab("Age (Years)") + ylab("Log2 Expression") 


lm(stripchart.age.plot["hsa-miR-124-5p", ]~microglia.sample.age)
ggplot(stripchart.age.ggplot, aes(x = stripchart.age.ggplot[,'microglia.sample.age'], stripchart.age.ggplot[,'hsa-miR-124-5p'])) + geom_point(colour = 'gray32', fill = 'gray20', stroke = 1.5, alpha = 1) + geom_point(alpha=0.4, size=3) + theme_classic(base_size = 17) + geom_abline(intercept =  1.760992, slope = -0.005928, size = 1) + ggtitle("hsa-miR-124-5p") + xlab("Age (Years)") + ylab("Log2 Expression") 

lm(stripchart.age.plot["hsa-miR-223-5p", ]~microglia.sample.age)
ggplot(stripchart.age.ggplot, aes(x = stripchart.age.ggplot[,'microglia.sample.age'], stripchart.age.ggplot[,'hsa-miR-223-5p'])) + geom_point(colour = 'gray32', fill = 'gray20', stroke = 1.5, alpha = 1) + geom_point(alpha=0.4, size=3) + theme_classic(base_size = 17) + geom_abline(intercept =  6.811108 , slope = 0.001139, size = 1) + ggtitle("hsa-miR-223-5p") + xlab("Age (Years)") + ylab("Log2 Expression") 



## Volcano ##

EnhancedVolcano(Toptable.age.Microglia, lab = NA, x= "logFC", y = "P.Value",
                FCcutoff = 0,
                labSize = 5,
                pCutoff = 0.05,
                xlab = "Regression coefficient",
                pCutoffCol = "adj.P.Val",
                title = 'Age miRNAs Microglia',
                legendLabels=c('','Not Sig','', 'p < 0.05'),
                legendPosition = "top",
                legendLabSize = 15,
                gridlines.major = FALSE, gridlines.minor = FALSE,
                pointSize = 3.0,
                colAlpha = 1)

EnhancedVolcano(Toptable.age.Bulk, lab = NA, x= "logFC", y = "P.Value",
                FCcutoff = 0,
                labSize = 5,
                pCutoff = 0.05,
                xlab = "Regression coefficient",
                pCutoffCol = "adj.P.Val",
                title = 'Age Bulk Microglia',
                legendLabels=c('','Not Sig','', 'p < 0.05'),
                legendPosition = "top",
                legendLabSize = 15,
                gridlines.major = FALSE, gridlines.minor = FALSE,
                pointSize = 3.0,
                colAlpha = 1)


## GO AND KEGG

multimir_results_up <- get_multimir(org = 'hsa',
                                    mirna   = hsa.microglia.age.up$genes,
                                    table   = 'validated',
                                    summary = TRUE)


age.up.enriched.microglia.targets <- multimir_results_up@summary


go.microglia.up <- goana(age.up.enriched.microglia.targets$target_entrez, species = "Hs")
go.microglia.up$p.adj <- p.adjust(go.microglia.up$P.DE, method = "fdr")

kegg.microglia.up <- kegga(age.up.enriched.microglia.targets$target_entrez, species = "Hs")
kegg.microglia.up$p.adj <- p.adjust(kegg.microglia.up$P.DE, method = "fdr")



multimir_results_down <- get_multimir(org     = 'hsa',
                                      mirna   = hsa.microglia.age.down$genes,
                                      table   = 'validated',
                                      summary = TRUE)


age.down.enriched.microglia.targets <- multimir_results_down@summary

go.microglia.down <- goana(age.down.enriched.microglia.targets$target_entrez, species = "Hs")
go.microglia.down$p.adj <- p.adjust(go.microglia.down$P.DE, method = "fdr")

kegg.microglia.down <- kegga(age.down.enriched.microglia.targets$target_entrez, species = "Hs")
kegg.microglia.down$p.adj <- p.adjust(kegg.microglia.down$P.DE, method = "fdr")


