#To complete this week's assessments, we will be analyzing a much larger dataset than what we considered in the previous sections. 
# We will be using a subset of the methylation data available from The Cancer Genome Atlas project (TCGA). 
# This is not a toy example, but rather a dataset similar to those being analyzed in academic research labs. 
#  We have created a GitHub repository External link with all the data. You can download a ZIP External link file here. 
# The file is 801 MB in size, so even on a fast internet connection it will take several minutes to download. 
# We recommend you start downloading the file before you engage in the videos and assessments. 
# You will need to download and unzip this ZIP file to follow the videos and complete the exercises.
# In the cell composition assessment, we will be downloading a dataset from GEO directly from R. 
# a ratio set is not just methylation data but an object that has metrics on a comparison to baseline methylation

library(minfi)
#preprocessed data for Illumina methylation microarrays, mapped to a genomic location.
grset <- getGenomicRatioSetFromGEO("GSE32148") 
options(timeout = max(1000, getOption("timeout")))
getOption('timeout')
options(timeout = 100)
##after download  you can save the object like this:
save(grset,file="grset.rda")
##then to load later
load("grset.rda")


#CpG Island Shores Assessment
# We assume you have downloaded the TCGA data methylation dataset in the tcgaMethylationSubset External link GitHub repo. (~840 MB zip file).
# Create a variable named path specifying the path to the folder where you have downloaded this data.
path <- "C:/Users/jabba/Downloads"

# Make sure the following Bioconductor packages are downloaded if you haven't done so already:
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(tidyr)

# CpGIslandShores Assessment Q1
#Read in the sample annotation table (where path is as described above).
targets <- read.delim(file.path(path, "targets.txt"), as.is=TRUE) #file.path is like a shorter version of paste0

    
#How many samples are represented in this table?
  
targets$sample_type


 
# How many samples are from normal colon samples?
which(targets$normal_control_type == "Normal Tissue" & targets$patient.tumor_tissue_site == "COLON")




# For the following questions we will perform a differential methylation analysis using normal breast and colon samples from TCGA:

index <- which(targets$Status=="normal" & targets$Tissue%in%c("colon", "breast"))
targets <- targets[index, ]
#Now we are ready to read in the data (this will take about 2 minutes):

library(minfi)
dat <- read.metharray.exp(base=path, targets=targets, verbose=TRUE) # Creates  a RGChannel set Object

# The dat object includes the raw data. To convert this into an object that includes methylation values, 
# as well as the location of CpGs, we run the following commands (we show you the class of dat as we transform it):
## preprocess the data
dat <- preprocessIllumina(dat) #Normilization and produces a methyl set object from a RGChannel set Object
dat <- mapToGenome(dat)## assign locations to each CpG 
dat <- ratioConvert(dat, type="Illumina")## precompute methylation values from U and M values making it a Genomic Ratio Object


# Before we start we can create some quality assessment plots. First look at the distribution of each sample:
install.packages("rafalib")
library(rafalib)
mypar(1,1)
y <- getBeta(dat)##extract methylation values
shist(y)

# Note that the distributions seem similar. Nothing stands out.
# We also create an MDS plot to search for outlier samples. The first PC splits the data by tissue as 
# expected and no sample stands out as an outlier.

mds <- cmdscale(dist(t(y)))
tissue <- as.factor(pData(dat)$Tissue)
plot(mds, col=tissue)

# Now we are ready to use statistical inference to find differentially methylated regions. 
# Let's start by using the limma package to perform a site-by-site analysis.

library(limma)

      
#CpGIslandShores Assessment Q3
#Find the CpG with the largest effect size when comparing the normal breast and colon samples.

topTreat(eb) #cg25753567
index = which.max(abs( fit$coef[,2]))
seqnames(dat)[index]
start(dat)[index]


#Which chromosome is this CpG on?
seqnames(dat)[index]
  

#CpGIslandShores Assessment Q4
# when many features are under consideration pvalues cutoffs will not be enough we need to calculate a qvalue(s)  
# Now we will use the qvalue() function to determine the q-value for the CpG found in the previous question.
library(qvalue)

##create design matrix
tissue <- as.factor(pData(dat)$Tissue)

X <- model.matrix(~tissue)
   y <- getBeta(dat)##extract methylation values
       fit <- lmFit(y,X)## obtain effect sizes and pvals with limma
          eb <- eBayes(fit)

## obtain q-values
qvals <- qvalue(eb$p.value[,2])$qvalue


#What is the q-value for this CpG? Since the value is incredibly small, report the -log10(q-value) as your answer.
which(names(qvals) %in% "cg22365276")
qvals[which(names(qvals) %in% "cg22365276")] 
-log10(1.269936e-27 )


########################################I couldnt get this answer  ##############################################################################
# Find all the CpGs within 5000 basepairs of the CpG identified in the previous question. 
# Create a plot showing the methylation values across all of the samples at these CpGs. 
# Use color to distinguish breast from colon samples. 
# Additionally, create plots of the estimated effect sizes and the -log10(q-values) along these CpG sites. You should have a total of three plots.
library(rafalib)
library(AnnotationHub)
cgi<-AnnotationHub()[["AH5086"]] # this number is the id for CPG islands
index<- which(seqnames(dat)== "chr15")
# :(
# this makes 1000 permutations using a Linear model function and estimates regions where methylation deviates from base line
# from a genomic ratio set
res<-bumphunter(dat[index], X, cutoff=.1, B= 1000) 
tab<- res$tab[res$tab$fwer<=.05,] # take only significant bumps
tab<-makeGRangesFromDataFrame(tab,keep.extra.columns = TRUE)# makes a data frame from it
# :(
map<- distanceToNearest(tab,cgi) #computes distances from bumps to cg islands
d<- mcols(map)$distance #retrieve distances as a list of numbers
# :(
# :(
# distances are binned based on distance reveal who is in the island who is just out side the island (shore) 
# and who is way off (shelf)
prop.table(   table ( cut(d,c(0,1,2000,5000,Inf), include.lowest = TRUE, right= FALSE  )))
# :(
# :(
# the same data can be retrieved from the raw illumina data without bumphunting...
# :(
null<-granges(dat)
nulltab<- makeGRangesFromDataFrame(null,keep.extra.columns = TRUE)
nullmap<- distanceToNearest(nulltab,cgi)
prop.table(   table ( cut(nulld,c(0,1,2000,5000,Inf), include.lowest = TRUE, right= FALSE  )))
#or
getIslandStatus(dat) %>% table()
# :(
# :(
# :(
tab<-tab+5000
# :(
dat[329905]
which(start(ranges(granges(dat))) == 60688622) # 360649 cg17495912 
# :(
# :(
ranges(granges(cgi))
# :(
dataindex<-which(granges(dat)[360649] %over% tab )#  retrieve which of our data is covered in signigicant bumps from  bumbhunter result
cgindex<-which(granges(cgi) %over% tab)# retrieve which of our annotation is covered in  significant bumps from bumbhunter result
thecgi<- cgi[cgindex]# which ever cgi index we choose use start and stop as xlim
#use those indexs to plot
pos<- start(dat)[dataindex] #position for x
colors<- as.factor(pData(dat)$Tissue) #color code by tissue
beta<-getBeta(dat)[dataindex]# methylation values for y
xlim<- range( c(pos,start(thecgi),end(thecgi))) #limits are based on island region
#plot
matplot(pos,beta, col = as.numeric(colors), xlim = xlim  , ylim= c(0,1)  )  
# :(
######################################################################################################################




###########################BELOW IS THE PROFESSORS SOLUTION ABOVE CODE GIVES WRONG ANSWERS#################################
# Based on these plots, which best describes this region?
# :)
# :) 
library(rafalib)
mypar(3,1)
index = which.max(abs( fit$coef[,2]))
gr=granges(dat)[index]+5000
index=which(granges(dat)%over%gr)
pos= start(dat)[index]
# :) 
matplot(pos,y[index,],ylab="Methylation",col=as.numeric(tissue))
# :) 
plot(pos, fit$coef[index,2],ylab="Effect Size")
# :) 
plot(pos, -log10(qvals[index]) ,ylab="-log10 q-value")
# :) 
# :) 
# :) 
# Repeat the above exercise for the top 10 CpGs ranked by absolute value of the effect size and identified with the following code:  
o <- order(abs(fit$coef[,2]), decreasing = TRUE)[1:10]
library(rafalib)
mypar(3,1)
o <- order(abs(fit$coef[,2]), decreasing = TRUE)[1:10]
for(i in o) {
  index <- i
  gr=granges(dat)[index]+5000
  index=which(granges(dat) %over% gr)
  pos <- start(dat)[index]
  # :)  
  matplot(pos, y[index,,drop=FALSE], ylab="Methylation", col=as.numeric(tissue))
  plot(pos, fit$coef[index,2], ylab="Effect Size")
  plot(pos, -log10(qvals[index]), ylab="-log10 q-value")
}
# :) 
###########################ABOVE IS THE PROFESSORS SOLUTION MINE GAVE WRONG ANSWERS#################################



# Now we are going to explicitly search for regions using the bumphunter() function. 
# We will use permutation to assess statistical significance. Because the function is slow, we will restrict our analysis to chromosome 15.
index <- which(seqnames(dat)=="chr15")
dat2 <- dat[index,]



# If your computer has more than one core, you can use parallel computing to speed up the procedure.
install.packages("doParallel")
library(doParallel)
ncores <- detectCores()
registerDoParallel(cores = ncores)



# We can now run the bumphunter() function to find differentially methylated regions (DMRs).
# For this assessment question we will use 100 permutations, 
# although we recommend more in practice. Here we will use a cutoff of 0.1. 
# The permutations are random so make sure you set the seed to 1 to obtain exact results in the assessment question.



##create design matrix
tissue <- as.factor(pData(dat)$Tissue)
X <- model.matrix(~tissue)

##extract methylation values
set.seed(1)
res <- bumphunter(dat2,X,cutoff=0.1,B=100)
head(res$tab)
table(res$table$fwer < .05 )
 





# Previously we performed a CpG by CpG analysis and obtained qvalues.
# Create an index for the CpGs that achieve qvalues smaller than 0.05 and a large effect size larger than 0.5 (in absolute value):
##fit and qvals were defined in a previous answer
index <- which(qvals < 0.05 & abs(fit$coef[,2]) > 0.5 & seqnames(dat)=="chr15")


#Now create a table of the DMRs returned by bumphunter that had 3 or more probes and convert the table into GRanges:
tab <- res$tab[ res$tab$L >= 3,]
tab <- makeGRangesFromDataFrame(tab)


#What proportion of the CpGs indexed by index are inside regions found in tab?
# Hint: use the findOverlaps() function in the GenomicRanges package.
countOverlaps(dat[index],tab)
12/21
unanswered 


# Now download the table of CpG islands using AnnotationHub.
library(AnnotationHub)
cgi <- AnnotationHub()[["AH5086"]]


# Now we create a GRanges object from the list of DMRs we computed in the previous questions:
tab <- res$tab[res$tab$fwer <= 0.05,]
tab <- makeGRangesFromDataFrame(tab)


# What proportion of the regions represented in tab do not overlap islands, but overall CpG islands shores (within 2000 basepairs) ?
# Hint: use the distanceToNearest() function.
X<-distanceToNearest(cgi ,tab) 
d<-X@elementMetadata@listData$distance
prop.table(   table ( cut(d,c(0,2000), include.lowest = TRUE, right= FALSE  )))
getIslandStatus(X) %>% table()
#or
map <- distanceToNearest(tab,cgi)
d <- mcols(map)$distance
mean(d<=2000 & d>0)
