genes<- read.table(file = "C:/Users/jabba/Downloads/SRR1039508.genes.results" , header = TRUE, stringsAsFactors = TRUE) 
isoforms<- read.table(file = "C:/Users/jabba/Downloads/SRR1039508.isoforms.results", header = TRUE, stringsAsFactors = TRUE)


# isoform table is a spread out format of genes 
# notice the FPKM in the genes table  is the sum of all transcripts for that gene in the Isoforms
# table  
library(dplyr)
isoforms %>% group_by(gene_id) %>% summarize(sum = sum(FPKM)) 



#Make a histogram of the FPKM values in genes. Make a histogram after transforming by: log10(x+1). Before moving on, we will first remove genes with FPKM values of 0.
genes2 <- genes[genes$FPKM > 0,]
genes2$gene_id <- droplevels(genes2$gene_id)
isoforms2 <- isoforms[isoforms$gene_id %in% genes2$gene_id,]
isoforms2$gene_id <- droplevels(isoforms2$gene_id)

#Run the following code to verify that the gene_id column in genes2 is equal to the levels of the gene_id column in isoforms2:
stopifnot(all(genes2$gene_id == levels(isoforms2$gene_id)))


library(tidyverse)
library(ggplot2)

genes2 %>% group_by(.,gene_id) %>% 
ggplot(data=., aes(x=effective_length,y=expected_count))+
  geom_point()+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")


genes2 %>% group_by(.,gene_id) %>% 
  ggplot(data=., aes(x=FPKM))+
  geom_histogram()+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")

genes2 %>% group_by(.,gene_id) %>% 
  ggplot(data=., aes(x=FPKM))+
  geom_bar() +
  scale_x_binned()+
  #scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")

summary(genes2)


isoforms2 %>% group_by(.,gene_id) %>% summarise(., max.Iso = max(IsoPct))

quantile(x$max.Iso,.95) 

isoforms2 %>% group_by(.,gene_id) %>% summarise(., max.Iso = max(IsoPct)) %>% dim()

isoforms2 %>% group_by(.,gene_id) %>% summarise(., max.Iso = max(IsoPct)) %>% filter(max.Iso >95.00) 

(1835/3045)*100

isoforms2 %>% group_by(.,gene_id) %>% summarise(., max.Iso = max(IsoPct)) %>% 
  ggplot(data= ., aes(x=max.Iso))+
  geom_histogram()

 









merge(isoforms2[c(2,8)] %>% group_by(.,gene_id) %>% summarise(., max.Iso = max(IsoPct))
      ,genes2[,c(1,7)],by = "gene_id") %>% 
 #group_by(.,gene_id) %>% 
  ggplot(data=., aes(x=max.Iso,y=FPKM))+
  geom_point()+
  #scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")


filter(merge(isoforms2[c(2,8)] %>% group_by(.,gene_id) %>% summarise(., max.Iso = max(IsoPct))
             ,genes2[,c(1,7)],by = "gene_id") , max.Iso >=12.5 & max.Iso <= 30.1 ) %>% summary()

filter(merge(isoforms2[c(2,8)] %>% group_by(.,gene_id) %>% summarise(., max.Iso = max(IsoPct))
             ,genes2[,c(1,7)],by = "gene_id") , max.Iso >=30.1 & max.Iso <= 47.6 ) %>% summary()

filter(merge(isoforms2[c(2,8)] %>% group_by(.,gene_id) %>% summarise(., max.Iso = max(IsoPct))
             ,genes2[,c(1,7)],by = "gene_id"), max.Iso >=47.6 & max.Iso <= 65  ) %>% summary()

filter(merge(isoforms2[c(2,8)] %>% group_by(.,gene_id) %>% summarise(., max.Iso = max(IsoPct))
             ,genes2[,c(1,7)],by = "gene_id") , max.Iso >=65 & max.Iso <=82.5 ) %>% summary()

filter(merge(isoforms2[c(2,8)] %>% group_by(.,gene_id) %>% summarise(., max.Iso = max(IsoPct))
             ,genes2[,c(1,7)],by = "gene_id") , max.Iso >=82.5 & max.Iso <= 100 ) %>% summary()

#####################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pasillaBamSubset" , "Rsamtools")
BiocManager::install("Rsamtools")
BiocManager::install("TxDb.Dmelanogaster.UCSC.dm3.ensGene")

library(pasillaBamSubset)
bam.file <- untreated3_chr4()
library(Rsamtools)
bf <- BamFile(bam.file)

library(Rsamtools)
bam.list <- BamFileList(bam.file , yieldSize = 1)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene) # this is already an assembled genome from a gtf file as a txdb object




txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
exons.by.gene <-GenomicFeatures::exonsBy(txdb, by="gene")


  
exons.by.gene[[1]]


#We can subset to only the genes on chromosome 4 using the following R command:
  
chr4.idx <- all(seqnames(exons.by.gene) == "chr4")
ebg.sub <- exons.by.gene[chr4.idx]

ebg.sub[[2]]

#this subset of exon by gene construction of our genome can now be aligned


se<-GenomicAlignments::summarizeOverlaps(ebg.sub,bam.list, 
                                     mode = "Union",
                                     ignore.strand=TRUE,    # experiment was not strand-specific
                                     singleEnd=FALSE,       # the experiment is paired-end
                                     fragments=FALSE)

table(assays(se)$counts>0)
table(se)
se[1]
colData(se)
library(GenomicAlignments)
x<-rowRanges(se)
subset(se, subset, select)
y<-x$FBgn0002521



BiocManager::install("DESeq2")
BiocManager::install("Biobase")

library(Biobase)
library(DESeq2)
download.file("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/wang_eset.RData", "wang_eset.RData")
load("wang_eset.RData")





count.matrix <- exprs(wang.eset)[,10:21]
col.data <- pData(wang.eset)[10:21,]

dds <- DESeqDataSetFromMatrix(count.matrix, col.data, design=~cell.type)

counts(dds)
assays(dds)
x<-estimateSizeFactors(dds)
colData(x)

#What tissue has the highest size factor (as calculated by DESeq's estimateSizeFactors())?

estimateSizeFactors(dds, cerebellum)
assays(dds)
DESeqSummarizedExperiment(dds)
#Copy the name exactly as it appears in your R session, only without quote marks.




#Run the varianceStabilizingTransformation() transformation (blind=FALSE) on the counts and save this to an object vsd.

vsd<-varianceStabilizingTransformation(counts(dds), blind = FALSE)
library(tidyverse)
df<-vsd %>% as.data.frame() %>% t() %>% as.data.frame() 



colnames(df)[52581]
rownames(df)

df$sample.id<-rownames(df)
 



x<-col.data[,1:2] %>% droplevels()
x$sample.id<-rownames(df)


x_1$cell.type<-x_1$cell.type %>% droplevels()


DESeqTransform(dds) %>% plotPCA(.,intgroup = "cell.type")
plotPCA(dds, intgroup="dex")

#Make a PCA plot with this object using "cell.type" as the color of the points.


#What sample clusters with cerebellum?
# Give the name of the cell-type exactly as it appears in the plot legend. Be careful - the colors can be similar.


rmeans <- rowMeans(assay(vsd)) # row mean of rlog-transformed data
idx <- c(1,2,10,7,8,9,12) # pick some samples for visualization
mat <- assay(vsd)[rmeans > 1,idx] # pull out a small matrix of rlog-transformed counts
colnames(mat) <- vsd$cell.type[idx] # name the columns of matrix by cell type

# Now we could already call pairs() on this matrix, and it would make a matrix of scatterplots. But we will add an extra bit of code that will make a fast, subsetted scatter plot on top and add Pearson correlations below (the correlation code is directly copied from ?pairs):

panel.sub <- function(x,y,...) points(cbind(x,y)[sample(length(x),1000),],...)
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)  {
	usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Now make our scatterplot matrix:

pairs(mat, asp=1, col=rgb(0,0,0,.3), lower.panel=panel.cor, upper.panel=panel.sub)





####  DE assessment I

download.file("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/wang_eset.RData", "wang_eset.RData")

load("wang_eset.RData")
library(Biobase)

count.matrix <- exprs(wang.eset)[,10:21]
col.data <- pData(wang.eset)[10:21,]

library(DESeq2)
dds <- DESeqDataSetFromMatrix(count.matrix, col.data, design=~cell.type)



dds$type <- factor(ifelse(dds$cell.type %in% c("cerebellum","mixed.brain"),
                          "brain", "other"))

dds$type <- relevel(dds$type, "other")


design(dds) <- ~ type


dds <- DESeq(dds)
res <- results(dds)
head(res)
table(res$padj < 0.1)
plotCounts(dds, gene=which.min(res$padj < 1), intgroup="type")

table(res$log2FoldChange >=2)

res.thr <- results(dds, lfcThreshold=2)
plotMA(res.thr, ylim=c(-10,10))



table(res.thr$padj < 0.1 & res.thr$log2FoldChange >=0)


#What is the ENSEMBL ID of the gene with the smallest adjusted p-value in res?
# Report the full gene ID starting with "ENSG". The rownames of res give the Ensembl ID.



 


#Now create a results table res2 but this time with a lfcThreshold of 2. You do not need to re-run DESeq(). Make an MA plot of res2 with ylim=c(-10,10).
#Call summary() on res2.
#How many genes with adjusted p-value less than 0.1 have a positive LFC?




#Earlier, we made a single counts plot for the gene with the smallest adjusted p-value:
  
plotCounts(dds, which.min(res$padj), intgroup="type") 

  
#With the following code, make normalized counts plots for the top 9 genes (smallest adjusted p-values):
  
par(mfrow=c(3,3))
for (i in 1:9) {
  plotCounts(dds, order(res$padj)[i], intgroup="type")
}


#These genes all look similar, with large counts for brain samples and small counts for the "other" group. Here we have empirically found a set of genes which seem specific (at least in our dataset) to brain. We examine their annotation below. Run the following code to return to showing one plot at a time (rather than 9):
  
par(mfrow=c(1,1))

  
#Build a vector of the names of the top 20 genes by the test statistic (largest first). The test statistic is the LFC divided by its standard error, and uniquely determines the p-value. So this is the same as asking for the genes with the smallest p-value and positive LFC (so higher expression in brain).

top <- rownames(res)[head(order(res$stat, decreasing=TRUE), 20)]



# Use org.Hs.eg.db to determine the gene symbol of the top gene in this list.

# What is the SYMBOL of the top gene?


install.packages("Rtools")
R.Version()

BiocManager::install("AnnotationDbi", force = TRUE)
 
BiocManager::install("org.Hs.eg.db")

install.packages("C:/Users/jabba/Downloads/org.Hs.eg.db_3.18.0.tar.gz", repos = NULL, type="source")

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
anno <- select(org.Hs.eg.db, keys=top,
               columns=c("SYMBOL","GENENAME"), 
               keytype="ENSEMBL")
anno[match(topgenes, anno$ENSEMBL),]
 
#####





# The code below downloads the data of Bottomly et al. 
# from the ReCount project External link:
  
download.file("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData", "bottomly_eset.RData")

Note: if the file above is not available, you can also download bottomly_eset.RData External link hosted as part of this course.

load("bottomly_eset.RData")
library(Biobase)

#Now build a DESeqDataSet from this object:
  
count.matrix <- exprs(bottomly.eset)
col.data <- pData(bottomly.eset) 
install.packages('stringr')
library(stringr)
strain <- str_replace(col.data$strain, "/", "")
col.data$strain <- strain

BiocManager::install("DESeq2")
library(DESeq2)
dds <- DESeqDataSetFromMatrix(count.matrix, col.data, design=~strain)
# The experiment.number and lane.number columns are numeric, 
# so make sure to turn them into factors:
  
  dds$experiment.number <- factor(dds$experiment.number)
dds$lane.number <- factor(dds$lane.number)

#Estimate the size factors so we can get normalized counts later:
  
  dds <- estimateSizeFactors(dds)
  
EDA
# Run the varianceStabilizingTransformation() on the dds 
# and then make a PCA plot with c("strain","experiment.number") as the intgroup to label.
varianceStabilizingTransformation(dds) %>%
plotPCA(.,intgroup = c("strain","experiment.number"))  
  
  


#The strains have names "C57BL6J" and "DBA2J".
#The experiments are numbered 4, 6 and 7. 
#We can see that both strain and experimental batch have an effect on the normalized, transformed counts.

#Because we know the experimental batches, we could just use DESeq() 
# with ~ experiment.number + strain to look for strain specific differences 
# controlling for batch. 
# But suppose we were given this data without the batch information. 
# We could use SVA to try to identify the hidden structure.

Run SVA-seq
#Run SVA-seq to find two surrogate variables using the code shown in the previous video or corresponding Rmd file External link.
#Use a design of ~ strain for the full model "mod" and ~ 1 for the reduced model "mod0".
BiocManager::install("sva")
library(sva)
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~ strain , colData(dds))

mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv=2)
plot(svseq$sv[,1], svseq$sv[,2], col=dds$experiment.number, pch=16)

#Did SVA-seq recover batch?
#Plot the surrogate variables and color by the true batches 
#(remember, normally we wouldn't have this information if we needed to run SVA):

plot(svseq$sv[,1], svseq$sv[,2], col=dds$experiment.number, pch=16)
legend("bottom", levels(dds$experiment.number), pch=16, col=1:3)
#Now add numbers to each sample:

#Which sample (the number in the plot) from experiment number 6 
#is closest to the samples from experiment number 4, 
#according to these two surrogate variables?


text(svseq$sv[,1], svseq$sv[,2], 1:ncol(dds), pos=1)
 