


# In this assessment, we will examine the isoform-level abundances which are saved as an output from Cufflinks 
# and accessible via the Bioconductor package cummeRbund External link.
# We will eventually ask, for each gene: how often is the most highly expressed isoform the same across two biological conditions?
# Create a CuffSet object:

BiocManager::install("cummeRbund")
library(cummeRbund)

myDir <- system.file("extdata", package="cummeRbund") 

gtfFile <- system.file("extdata/chr1_snippet.gtf",package="cummeRbund")

cuff <- readCufflinks(dir=myDir,gtfFile=gtfFile,genome="hg19",rebuild=TRUE)

#Extract annotation information
#Extract the annotation information with the annotation() function. 
# This gives exon-level information. We just want to know the gene which each isoform belongs to, 
# so we can remove duplicate rows belonging to the same isoform:
  
gene.features <- annotation(genes(cuff))

head(gene.features)

isoforms.per.gene <- gene.features[!duplicated(gene.features$isoform_id), c("gene_id","isoform_id")]

isoforms.per.gene <- isoforms.per.gene[order(isoforms.per.gene$isoform_id),]

head(isoforms.per.gene)

#Number of isoforms per gene
#A single call to table() gives us a vector, one value for each gene, and the number of isoforms per gene:
  
gene.tab <- table(isoforms.per.gene$gene_id)



#How many genes have only 1 isoform?

x<-gene.tab %>% as.data.frame 
table(x)
which(x$Freq > 1)


# FPKM per isoform
# The fpkm() function returns a data.frame of the FPKM estimates for each isoform and sample:
  
isoform.fpkm <- fpkm(isoforms(cuff))

head(isoform.fpkm)

table(isoform.fpkm$sample_name)

# Extract out tables for the iPS and hESC samples:
  
ips <- isoform.fpkm[isoform.fpkm$sample_name == "iPS",]
hesc <- isoform.fpkm[isoform.fpkm$sample_name == "hESC",]

# Now check that the isoform_id from our FPKM tables and our isoforms-per-gene table are identical:

# If these functions run without error, the columns are equal.

stopifnot(all(ips$isoform_id == isoforms.per.gene$isoform_id))


#Which isoform has the highest expression?
# Use sapply(), split() and which.max() to identify, for each sample, the index of the isoform with the largest FPKM. For example:
  
ips.max <- sapply(split(ips$fpkm, isoforms.per.gene$gene_id), which.max)
head(ips.max)
hesc.max <- sapply(split(hesc$fpkm, isoforms.per.gene$gene_id), which.max)
head(hesc.max)

install.packages("caret")
library(caret)
hesc.max %>% as.factor() %>% levels()

ips.max = sapply(split(ips$fpkm, isoforms.per.gene$gene_id ), which.max)    
hesc.max = sapply(split(hesc$fpkm, isoforms.per.gene$gene_id), which.max)
mean(ips.max[x1]  == hesc.max[x1]) 

x1<-ips.max %>%  names() %in% x1 %>%  which()
hesc.max %>%  names() %in% x1 %>%  which()

ips.max[ips.max %>%  names() %in% x1 %>%  which()]
filter(ips.max )

x1<-filter(x, Freq>1 )$Var1




library(dplyr)
my_indices<-c("XLOC_000001", "XLOC_000002", "XLOC_000003", "XLOC_000004", "XLOC_000005", "XLOC_000006")
merge(filter(isoform.fpkm,sample_name != "Fibroblasts" ), isoforms.per.gene, by= "isoform_id") %>%  filter(., .$gene_id %in% c(my_indices))


#Across all genes in our tables, how often is the highest expressed isoform the same one in iPS and hESC cells? 
#Give a proportion (number between 0 and 1). 
# If your R code is throwing an error, try loading the DESeq2 package before library(cummeRbund) in the above code.
XLOC_000001 XLOC_000002 XLOC_000003 XLOC_000004 XLOC_000005 XLOC_000006 

genes(cuff)

head(gene.tab)
x
isoform.fpkm []
table(isoform.fpkm$isoform_id)

dispersionPlot(genes(cuff))
sigGeneIds<-getSig(cuff,alpha=0.05,level="isoforms")
sigGeneIds<-getSig(cuff,"iPS","hESC",alpha=0.05,level="isoforms")

sigGenes<-getGenes(cuff,sigGeneIds)
            

# Subsetting to only the genes that have more than one isoform, how often is the highest expressed isoform the same one in iPS and hESC cells?
# (Hint: you already have gene.tab calculated)

x1<-filter(x, Freq >1)
x1$Var1