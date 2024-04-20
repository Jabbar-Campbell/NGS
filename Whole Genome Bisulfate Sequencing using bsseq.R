# Measuring Methylation from Sequencing Q1
# Reduced Representation Bisulfite Sequencing or RRBS is an experimental technique widely used to manipulate the regions of the 
# genome we measure. An enzyme is used to cut DNA at CCGG and the general idea is to filter out small or large molecules once DNA is cut.
# We can use Bioconductor tools to predict the size of these regions.
# Load the genome package and create an object with the sequence for chr22:
  
library(tidyr) 
library("BSgenome.Hsapiens.UCSC.hg19")
chr22 <- Hsapiens[["chr22"]]


# Now use the matchPattern() function to find all the locations in which CCGG occurs on chr22.
matchPattern(pattern = "CCGG", subject = chr22) %>% length()

How many CCGG do we find on chr22?
  unanswered 


# Plot a histogram of the DNA fragment sizes after we cut at CCGG positions. How would you describe this distribution?
  
#The fragment sizes follow a normal distribution.

#The fragment size follow an exponential distribution.

# The distribution has a long right tail with most values from 0-1000, but some very large values.

# The fragment sizes follow a uniform distribution.

x<-matchPattern(pattern = "CCGG", subject = chr22)
table(x@ranges@width > 20)





# A typical size to filter are DNA regions between 40 and 220 basepairs.
#What proportion of the fragments created for chr22 are between 40 and 220 basepairs?

recognition_site <- "CCGG"  # Example recognition site for some restriction enzyme

# Define the DNA sequence for chr22 (example)
dna_sequence_chr22 <- as.character(chr22)

# Find the positions of recognition sites in the DNA sequence
positions_chr22 <- unlist(gregexpr(recognition_site, dna_sequence_chr22))

# Calculate fragment sizes for chr22
fragment_sizes_chr22 <- c(positions_chr22[1], diff(positions_chr22), nchar(dna_sequence_chr22) - positions_chr22[length(positions_chr22)] + 1)

# Calculate the proportion of fragments between 40 and 220 base pairs
fragments_between_40_and_220 <- sum(fragment_sizes_chr22 >= 40 & fragment_sizes_chr22 <= 220)
total_fragments <- length(fragment_sizes_chr22)

proportion <- fragments_between_40_and_220 / total_fragments

cat("Proportion of fragments between 40 and 220 base pairs:", proportion, "\n")


 
# If we try to sequence all of chromosome 22, we would need to sequence 51,304,566 bases. 
# If instead we only sequence fragments of sizes between 40 and 220 basepairs, how many bases would we need to sequence?
 
index<-which(fragment_sizes_chr22 >= 40 & fragment_sizes_chr22 <= 220 )


fragment_sizes_chr22[index] %>% sum()



# Whole-genome bisulfite sequencing (WGBS) is another sequencing-based protocol for measuring DNA methylation in samples. 
# In the remaining problems, we will take a look at WGBS data from a set of paired tumor and normal colon samples. 
# You can download the data from our GitHub repo External link (~40 Mb). 
# More information about these data are available on GEO GSE46644 External link. Here we use the .cov files from Bismark as 
# input to the bsseq R/Bioconductor package to create a bsseq object. 
# The data come from this paper External link.
# Step-by-step instruction on how to covert the raw FASTQ files to the files that we will be work on are available here External link.
# Let's start by the reading in the target information.

    
## path to the downloaded file collection
path <- "C:/Users/jabba/Downloads"

## read in sample metadata table
targets <- read.table(file.path(path,"targets.txt"), header = TRUE, sep = "	")
targets

  
# Now you will need the bsseq package to read in the sequencing data.

# We load the methylation calls from our alignments into R. Once the data are loaded into R, 
# we can use this package for further downstream analyses such as finding differentially methylated 
# regions between our paired tumor and normal samples. This package assumes the following data have been 
# extracted from the alignments.

# 1. genomic positions (chromosome and location) for methylation sites
# 2. M (Methylation) values, the number of reads supporting methylation covering each site
# 3. Cov (Coverage) values, the total number of reads covering each site

# For illustrative purposes, we only consider the methylation sites on chromosome 22 from the .cov file.

# Here, read in the six files:

BiocManager::install("bsseq")   
library("bsseq")

## turn metadata into DataFrame w/ sample names as rownames
targets <- DataFrame(targets, row.names = as.character(targets$Run))

## specify path to files in same order as targets table
covfiles <- file.path(path, paste0(rownames(targets), ".chr22.cov"))

## read coverage files
colonCancerWGBS <- read.bismark(files = covfiles, rmZeroCov = TRUE,
                                colData = targets)

  
#Take a look at the bsseq object and the phenotypic information about for sample:
colonCancerWGBS

# phenotypic information
pData(colonCancerWGBS)

# granges object
granges(colonCancerWGBS)

  
#Now, extract the coverage and the number of reads with evidence for methylation:
cov <- getCoverage(colonCancerWGBS,type = "Cov")
m <- getCoverage(colonCancerWGBS,type = "M")


# Check if each CpG position has coverage greater than zero in all samples
cpg_with_coverage_in_all_samples <- apply(cov, 1, function(row) all(row > 0))

# Calculate the proportion of CpGs with coverage in all samples
proportion_coverage_all_samples <- sum(cpg_with_coverage_in_all_samples) / length(cpg_with_coverage_in_all_samples)

cat("Proportion of CpGs with coverage in all samples:", proportion_coverage_all_samples, "\n")


 
# Compute the total coverage (across all samples) for each CpG site. Plot it against location.
# How would you describe coverage?
# Assuming you have already extracted coverage data
tot = rowSums(cov)
##there are some very large values
hist(tot)
loc= start(colonCancerWGBS)
##plot by pieces
for(i in 1:11){
  index=1:100000+100000*i ##very ad-hoc
  plot(loc[index],tot[index],cex=.5,ylim=c(0,300))
}
# Appears uniform across chr22
# Has some very large values (>200) as well as general varaibility
#Is 0 for most location   x
# Is above 300 for most locations   x











#Assuming cov>0, which of the following pairs gives us:
#1) an estimate of methylation status at each CpG and
# 2) a quantity proportional to the standard error of this estimate?



#m and cov
#m/cov and 1/sqrt(cov) :)
# m/cov and sqrt(cov)
# cov and m
