# Up to now, the assessments have focused on general characteristics of methylation and CpG sites across the human genome. 
#We are now ready to analyze methylation measurement data. We will be working with a small subset of TCGA data that we have 
# created for illustrative purposes. This is real DNA methylation data from cancer and normal samples. Note that we are going 
# to use some of the material we learned in course 3 (Advanced Statistics) and course 5 (Introduction to Bioconductor). If you 
# have not done so already, you will need to install the following library from the GitHub repository:


BiocManager::install("devtools")  
library(devtools)
install_github("genomicsclass/coloncancermeth")
# Now we can load the library as well as the needed data objects:
  
library(coloncancermeth)
data(coloncancermeth)
dim(meth)
dim(pd)
print(gr)
head(meth)
head(pd)
# Finding Differentially Methylated Regions in R Assessment Q1

# From dim(pd) we can see that there are a total of 26 samples.

#How many are cancer samples?
table(pd$Status)

# Which column of the meth matrix is a cancer sample and has BCR patient barcode "TCGA-A6-4107"?
  pd$bcr_patient_barcode
meth 
 

# Finding Differentially Methylated Regions in R Assessment Q2
 
# Use the methylation profiles to compute a Euclidean distance between each sample:
  
  
d <- dist(t(meth))
pd$Status

# Now use the cmdscale() function to create an MDS plot External link that graphically shows approximate distances between the samples, 
# using color to distinguish cancer and normal samples.
my_pca<-cmdscale(d) %>% as.data.frame

my_pca$status<-pd$Status

ggplot(my_pca, aes(x= V1,y= V2, color = status))+
  geom_point(aes(size=3, alpha = .5))


# Which of the following best describes the MDS plot?
  
 # x 1)The DNA methylation profiles appear similar for cancer and normal since the MDS plot shows random scatter.

 # x 2) The MDS plot shows perfect separation between cancer and normal samples in both dimensions.

 # * 3) The MDS plot shows separation between cancer and normal samples, but only in the first dimension. The second dimension seems to be associated with a large variability within the cancers.

 # x 4) We can't compute distances because methylation data is binary.

 
# For each CpG compute p-values for the cancer versus normal comparison using the limma package:

    
library(limma)
X <- model.matrix(~pd$Status) # model.matrix just models things in terms of numbers...Binary in this case
fit <- lmFit(meth,X) # create a linear model using X as the out come and , gene expression data (kinda like pytorch)!!!!
eb <- eBayes(fit) # perform stats on all the elements of fit which is a linear model based on meth data
pvals <- eb$p.value[,2]  # extract pvalues

# when many features are under consideration pvalues cutoffs will not be enough we need to calculate a qvalue(s)  
# Now use the qvalue() function in the qvalue package to obtain q-values aka False Discovery rates.
# q-value is the estimated pFDR for a list of all the features 
BiocManager::install("qvalue")
library(qvalue)
qvalue(pvals)
qvalue(pvals) %>% summary
# What proportion of CpG sites have q-values smaller than 0.05?
qvalue(pvals)$qvalue %>% summary

 qvalue(pvals)$qvalue %>% mean()
 
 plot(pvals,qvalue(pvals)$qvalue , ylab= "qvalue(FDR)", xlab = "pvalue")

# To refresh your memory on q-values, read the end of the multiple testing section External link.

 
#Finding Differentially Methylated Regions in R Assessment Q4
 


# Before high-throughput technologies were available, cancer epigenetics research focused on searching for CpG islands showings 
# higher levels of methylation in cancer (hypermethylated). Let's explore the data at hand in this regard.

# What proportion of the CpGs showing statistically significant differences (defined with q-values .23  in the previous question) are, 
# on average, higher in cancer compared to normal samples?
 
 qvals = qvalue(pvals)$qvalue     
 index = which(qvals<=0.05)
 diffs = fit$coef[index,2]
 mean(diffs > 0)

 
# Now let's determine which of the differentially methylated CpGs are in CpG islands.
 index2<- rownames(fit)[index]
 names(index)
 
# Let's redefine CpG islands as cgi using the code from a previous assessment:
  
  
library(AnnotationHub)
ah <- AnnotationHub()
cgi <- ah[["AH5086"]]


# What proportion of the differentially methylated CpGs are inside islands?
islands=gr[index]%over%cgi
hypermethylated=fit$coef[index,2]>0
prop.table( table(islands,hypermethylated) )
  




#Finding Differentially Methylated Regions in R Assessment Q6
 
# Now we will use the bumphunter package to separate the differentially methylated CpGs intro groups.

BiocManager::install("bumphunter")
library(bumphunter)
X <- model.matrix(~pd$Status)
chr <- as.character(seqnames(gr))
res <- bumphunter(meth,X,chr=chr,pos=start(gr),cutoff=0.1)


# From here we get a table of regions:
  
  
  head(res$table)


# Note that the bumphunter() function has options to assess uncertainty, which are turned on through the B argument. 
# However, these options make this function computationally intensive. 
  

dmrs <- res$table[ res$table$L>=3, ]


# Note that this table is not a GenomicRanges object, but we can turn it into one easily:
  
  
dmrs <- makeGRangesFromDataFrame(dmrs)


#For the regions in dmrs, find the distance to the closest island. (Hint: use distanceToNearest().)

x<-distanceToNearest(dmrs,  makeGRangesFromDataFrame(cgi))


# What proportion of DMRs overlap a CpG island?
  Hint: if a DMR overlaps an island, the distance to the closest island will be 0.

 x1<-x@elementMetadata@listData$distance == 0
 table(x1)
 4825/ (3078+4825) 

  
#What proportion of DMRs are within 2000 basepairs from a CpG island, but do not overlap?
  y<- x@elementMetadata@listData$distance <= 2000 & x@elementMetadata@listData$distance > 0
  table(y) 
 972/(972 +6931)
  unanswered 

