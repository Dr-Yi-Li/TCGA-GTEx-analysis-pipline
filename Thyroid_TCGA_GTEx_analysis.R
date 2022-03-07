library(data.table)
library(DESeq2)
#读入表型文件
pd <- fread("phenotype/TcgaTargetGTEX_phenotype.txt.gz")
pd <- as.data.frame.matrix(pd)
table(pd[,4])
#甲状腺正常组织样品ID
pd_N <- pd[pd[,4] == "Thyroid",]
pd_n <- as.character(pd_N$sample)
#甲状腺癌组织样品ID
pd_C <- pd[pd[,4] == "Thyroid Gland",]
pd_C <- pd_C[substr(pd_C[,1],14,15) == "01",]
pd_c <- as.character(pd_C$sample)
#读入表达矩阵文件
exp <- fread("gene_expression_RANseq/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2.gz")
exp <- as.data.frame.matrix(exp)
#表达矩阵预处理
exp[1:5,1:5]
rownames(exp) <- exp[,1]
exp[1:5,1:5]
exp <- exp[,2:ncol(exp)]
exp[1:5,1:5]
#获取正常组织和肿瘤组织表达矩阵
exp_n <- exp[,colnames(exp) %in% pd_n]
exp_c <- exp[,colnames(exp) %in% pd_c]
exp_nc <- cbind(exp_n,exp_c)
#表达矩阵数值转换
exp_nc <- (2^exp_nc - 1)
exp_nc[1:5,1:5]
exp_nc <- round(exp_nc,0)
exp_nc[1:5,1:5]
#DESeq2包进行差异分析
#生成样本信息的factor
condition <- factor(ifelse(substr(colnames(exp_nc),14,15) == "01","Tumor","Normal"))
#定义表达矩阵的信息
colData <- data.frame(row.names = colnames(exp_nc),condition)
colData[1:5]
#生成DESeq2差异分析输入文件
dds <- DESeqDataSetFromMatrix(countData = exp_nc, colData = colData, design = ~ condition)
#过滤低表达的gene
dds <- dds[rowSums(counts(dds)) > 1,]
#进行差异分析
dds <- DESeq(dds)
#获取DESeq2差异分析结果
res <- results(dds)
#获取差异基因表达矩阵
keep <- res$padj < 0.05 & (res$log2FoldChange < -2 | res$log2FoldChange > 2 ) & res$padj != 0
dds_k <- dds[keep,]
res_k <- results(dds_k)
exp_k <- assay(dds_k)
#保存结果
write.csv(res,"Thyroid_normal_cancer.csv",quote = F)
write.csv(res_k,"Thyroid_normal_cancer_k.csv",quote = F)
