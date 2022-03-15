##### install required packages ####
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("data.table", quietly = TRUE))
    install.packages("BiocManager")
if (!require("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")

##### load required packages ####
library(data.table)
library(DESeq2)

##### input the phenotype table ####
pd <- fread("phenotype/TcgaTargetGTEX_phenotype.txt.gz") #download from UCSC Xena cohort: TCGA TARGET GTEx
pd <- as.data.frame.matrix(pd)

##### input the expression table and data convert####
exp <- fread("gene_expression_RANseq/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz") #download from UCSC Xena cohort: TCGA TARGET GTEx
exp <- as.data.frame.matrix(exp)
rownames(exp) <- exp[,1]
exp <- exp[,2:ncol(exp)]

##### select the normal and cancer sample id with target disease from phenotype #### 
pd_N <- pd[pd[,5] == "Normal Tissue",]
table(pd_N[,4]) #view the primary site of the samples to determine the label
pd_N <- pd_N[pd_N[,4] == "<your sample type>",] #select the normal change <your sample type> into your target tissue type
#pd_N <- pd_N[pd_N[,4] == "Thyroid",]
pd_n <- as.character(pd_N$sample)
pd_C <- pd[pd[,5] == "Primary Tumor",]
table(pd_C[,4]) #view the primary site of the samples to determine the label
pd_C <- pd_C[pd_C[,4] == "<your sample type>",] #select the normal change <your sample type> into your target tissue type
# pd_C <- pd_C[pd_C[,4] == "Thyroid Gland",]
pd_c <- as.character(pd_C$sample)

##### extract the expression matrix of normal and cancer samples ####
exp_n <- exp[,colnames(exp) %in% pd_n]
exp_c <- exp[,colnames(exp) %in% pd_c]
exp_nc <- cbind(exp_n,exp_c) #merge normal and cancer sample
#expression matrix convert
exp_nc <- (2^exp_nc - 1)
exp_nc <- round(exp_nc,0)

##### DEGs analysis by DESeq2 packages ####
condition <- factor(ifelse(substr(colnames(exp_nc),14,15) == "01","Tumor","Normal")) #generate the factor of sample info
colData <- data.frame(row.names = colnames(exp_nc),condition) #defined the info of expression matrix
dds <- DESeqDataSetFromMatrix(countData = exp_nc, colData = colData, design = ~ condition) #Generate the DESeq2 input file
dds <- dds[rowSums(counts(dds)) > 1,] #filter the low express genes
dds <- DESeq(dds) #differential analysis
res <- results(dds) #obtain the analysis result of DESeq2
#extract the DEGs
keep <- res$padj < 0.05 & (res$log2FoldChange < -2 | res$log2FoldChange > 2 ) & res$padj != 0
dds_k <- dds[keep,]
res_k <- results(dds_k)
exp_k <- assay(dds_k)

##### save the analysis result ####
write.csv(res,"Thyroid_normal_cancer.csv",quote = F)
write.csv(res_k,"Thyroid_normal_cancer_k.csv",quote = F)
