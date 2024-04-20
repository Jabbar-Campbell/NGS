#Cell Composition Assessment Q1
 
#To examine the importance of accounting for cellular composition in DNA methylation 
# analysis, we are going to download a GEO dataset used in an analysis of whole blood data.
# The minfi package has a function to read data directly from GEO. 
# Run the following commands. Note that this command downloads 64.7 MB of data and 
# can take several minutes, depending on your download speed.

NOTE: Occasionally downloads from GEO fail. If you receive a "cannot open URL" error running the code below, try again in a couple minutes.


library(minfi)
grset <- getGenomicRatioSetFromGEO("GSE32148")
gset <- getGEO("GSE32148", GSEMatrix =TRUE, getGPL=FALSE)


getOption("timeout")
options(timeout = 1000000)


#This creates an object of class:
  
  
  class(grset)


# which is quite convenient. Use the function pData to examine the sample 
# information table and determine which column includes the age of the individual.
grset@colData@listData$characteristics_ch1.1
# What is the name of this column?
grset@colData@listData$contact_name
  
  


library(minfi)
path <- "C:/Users/jabba/Downloads"    # same as CpGIslandShores assessment
targets <- read.delim(file.path (path,"targets.txt"),as.is=TRUE)
index <- which( targets$Tissue=="colon") #find the index  for colon data
targets <- targets[index,] #subset the targets file by this index
dat <- read.metharray.exp(base=path, targets=targets, verbose=TRUE) # read in only the base name of the subsetted targets


# The dat object again includes the raw data. 
# To convert this into an object that includes methylation values, 
# as well as the location of CpGs, we run the following lines of code:
  
  
dat <- preprocessIllumina(dat)
dat <- mapToGenome(dat)


#Now we can collapse the data as described in the video:
  
  
cdat <- cpgCollapse(dat)


#The original data includes nrow(dat) CpGs.
#How many regions are represented in the collapsed object?
library(tidyr)
cdat$blockInfo$pns %>% unique() %>% dim
  
  


# We can see the type of regions that are represented in this collapsed object:
  
head(granges(cdat$obj))


 
# OpenSea regions refer to genomic areas where DNA methylation occurs but is not specifically associated with CpG islands or 
# other well-defined regulatory elements like a  shore or a shelf
# What proportion of the regions are OpenSea regions?
granges(cdat$obj)$type %>% table()

113857/(26571  +  35725 +  47344 +113857)
 

# Now we use the blockFinder() function to find differentially methylated regions between cancer and normal:
# similar to ebayes linear modelfrom limma
status <- factor(pData(cdat$obj)$Status,
                   level=c("normal","cancer"))
X <- model.matrix(~status)
  res <- blockFinder(cdat$obj,X,cutoff=0.05) # similar to ebayes linear model it uses bumhunter under the hood
  
#Since blockFinder() calls bumphunter(), the returned results should look familiar. We can take a peek at the blocks:
  head(res$table)
# What proportion of the blocks reported in res$table are hypomethyated in cancer (lower methylation in cancer versus normal)?
table(res$table$value < .05 )
3017/(167  +3017 )

