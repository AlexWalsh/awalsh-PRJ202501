# awalsh-PRJ202501
human microglia miRNA mRNA analysis


DE analysis

**mRNA analysis**  
Input file: hsa_mRNA_genecounts
1.hsa_limma_mRNA.R
Format raw gene counts

**Cibersort analysis**  
Input files: hsa_rawcounts_cibersort.txt, zhou_reference.txt
1. preprocess_cibersort.R
2. cibersort.sh

**mRNA analysis**  
Input files: hsa_mRNA_genecounts.txt, CIBERSORTx_Adjusted.txt
1. hsa_limma_mRNA.R
- filter out samples based on <50% microglia purity
2. limma-voom DE analysis

**miRNA analysis**  
Input files: hsa_known_miRNA_counts.txt, CIBERSORTx_Adjusted.txt
1. deduplicate.R
2. hsa_known_sRNA_limma.R
- filter out samples based on <50% microglia purity
-limma-voom.DE analysis

**miRNA-mRNA network**  
Input files: hsa_smRNA_lcpm_resid.rds, hsa_mRNA_lcpm_resid.rds
1. generate_corr.RDS
2. miRNA-mRNA-hsa_networks-analysis.R


   
