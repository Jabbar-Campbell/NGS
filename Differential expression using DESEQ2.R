# We will inspect a subset of the experiment from the 48h time point to compare cells 
# treated with DPN (an ER- selective agonist) to control cells. Abbreviate this dataset as pth.
install.packages(pa)

BiocManager::install("parathyroidSE")
library(parathyroidSE)
data(parathyroidGenesSE)
pth <- parathyroidGenesSE    # abbreviate
pth <- pth[, pth$time == "48h" & pth$treatment != "OHT"]    # remove 24h and OHT samples
pth
# Define a DESeqDataSet with the design ~ patient + treatment. Call the resulting object dds:
  
library(DESeq2)
colData(pth) #choose the columns you'd like to group by in design
dds <- DESeqDataSet(pth, design = ~ patient + treatment)
# We need to do one last bit of cleanup before analysis: two of the columns are technical 
# replicates of the same sample. We can collapse these replicates, combining reads from the two runs, 
# into a single sample and store the resulting DESeqDataSeq as dds_coll:
  
dds_coll <- collapseReplicates(dds, groupby = dds$sample, run = dds$run)



#  In the context of DESeq2, a size factor is a crucial normalization factor used to account for 
# differences in library size (total read counts) across samples in RNA-seq data analysis. 
# These factors ensure that the expression levels of genes are appropriately adjusted for variations in sequencing depth

# Calculate the size factors for samples in dds_coll.
library(tidyr)
size_factors <- estimateSizeFactors(dds_coll) 


# What is the size factor for the first sample, "SRS308866"?

dds_coll <- estimateSizeFactors(dds_coll)
sizeFactors(dds_coll)[1]

#Which sample has the largest size factor?
#Report the name of the sample starting with "SRS".
sizeFactors(dds_coll) %>% which.max()




# PCA plot

# Use rlog() to perform a regularized log transformation on dds_coll, 
# then make a PCA plot of the transformed data using plotPCA(). The samples should cluster by patient.
#What percent of the variance is represented by the first principal component (PC1)?
#Which patient's samples have the largest values of PC1?
rld<-rlog(dds_coll)
colData(rld)
plotPCA(rld, intgroup="patient")




# Use DESeq() to perform differential expression analysis on dds_coll and use summary() 
# on the results to determine the number of differentially expressed genes.

dds <- DESeq(dds_coll)
res <- results(dds)
head(res)

summary(res)





# How many genes are upregulated (LFC > 0) at a FDR of 0.1?
# How many genes are downregulated (LFC < 0) at a FDR of 0.1?
# How many genes have counts which are too low to calculate differential expression?
  summary(res)







  
  
# What is the ENSEMBL ID of the gene with the lowest adjusted p-value?
# Report the entire ID starting with "ENSG".
  
index<-which.min(res$padj)
rownames(res)[index] #ENSG00000044574
head(res)
  

#What is the log2 fold change in expression of this gene in DPN-treated samples versus controls?
which (rownames(res)== "ENSG00000044574")
x<-res %>% as.data.frame(.)
x[615,]
res[615]
which(dds_coll %>% names() == "ENSG00000044574") 
dds_coll[615]
 



# Create an MA plot of the genes in the dataset. (You can also answer the questions or verify your responses by manually 
# filtering the results table.)
# Which of the following statements are true about the data?
plotMA(res, ylim=c(-4,4))


Check ALL correct answers. Note that "significantly differentially expressed" refers to an adjusted p-value < 0.1, which is the default value used by plotMA.


# There are no significantly differentially expressed genes with a log fold change larger than 2.

# There are no significantly differentially expressed genes with log fold change less than -2.

# There are no significantly differentially expressed genes with mean normalized counts less than 1e+02.

# There are no significantly differentially expressed genes with mean normalized counts greater than 1e+05.




# Use plotCounts() and ggplot2 to visualize counts for the gene with the lowest adjusted p-value. 
# Your x-axis should be the treatment type, your y-axis should be the count, and the point shape should be the patient.
# Which of these is the correct plot for the gene with the lowest adjusted p-value?

library(ggplot2)
colData(dds_coll)
data <- plotCounts(dds_coll, gene=which.min(res$padj), intgroup=c("treatment","patient"), returnData=TRUE)
ggplot(data, aes(x=treatment, y=count, col=patient)) +
  geom_point(position=position_jitter(width=.1,height=0)) +
  scale_y_log10()


# Connecting by lines shows the differences which are actually being tested by *results* given that our design includes `cell + dex`


ggplot(data, aes(x=treatment, y=count, col=patient, group=patient)) +
  geom_point() + geom_line() + scale_y_log10() 


A heatmap of the top genes:
  
  ```{r}
library(pheatmap)
topgenes <- head(rownames(resSort),20)
mat <- assay(rld)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("dex","cell")])
pheatmap(mat, annotation_col=df)

Gene expression increases in all DPN-treated samples relative to control samples. Patient 1 has the highest counts under both conditions.

Gene expression increases in all DPN-treated samples relative to control samples. Patient 3 has the highest counts under both conditions.

Gene expression decreases in all DPN-treated samples relative to control samples. Patient 3 has the highest counts under both conditions.

Gene expression decreases in all DPN-treated samples relative to control samples. Patient 3 has the highest counts under both conditions.
unanswered
SaveSave your answer
Submit
You have used 0 of 2 attempts