#load libraries
library(edgeR)
library(DESeq2)
library(limma)
library(Biobase)
library(dplyr)

#reading RSEM Expected Count 
exprs <- read.delim('RSEM_ExpectedCount GTEX Pancreas TCGA PAAD.txt', sep ='\t')

#reading which sample belongs to which study: TCGA or GTEx
sample_study <- read.delim('sample_study.txt',sep='\t')

#Creating a DGElist 
dge <- DGEList(counts = (2 ^ exprs[,2:345]) - 1, lib.size = colSums((2 ^ exprs[,2:345]) - 1), norm.factors = rep(1,ncol(exprs[,2:345])), sample = NULL, group = NULL, genes = exprs[,1], remove.zeros = FALSE)

#Defining model matrix based on study: TCGA or GTEx
design <- model.matrix(~sample_study$study)

#Filtering out genes with low expression 
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge_counts <- dge$counts

#Voom normalisation followed by quantile normalisation 
voom_expression <- voom(dge, design, plot=TRUE, normalize.method = 'quantile')

#Fiting model
fit.voom <- lmFit(voom_expression, design)

#Calculating t-statistics
fit.voom <- eBayes(fit.voom)

#Extracting data
TopGenes <- topTable(fit.voom, coef=ncol(design), number = 20000 ,sort.by = "p")
write.table(voom_expression, file= "voom_expression_quantile.txt", sep='\t')