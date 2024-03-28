

# CpG dinucleotides have long been observed to occur with a much lower frequency in the sequence of vertebrate genomes 
# than would be expected due to random chance. For example, in the human genome, which has a 42% GC content,[4] a pair 
# of nucleotides consisting of cytosine followed by guanine would be expected to occur 
# 0.21×0.21=4.41 of the time. The frequency of CpG dinucleotides in human genomes is less than one-fifth of the expected frequency.[5]

# This underrepresentation is a consequence of the high mutation rate of methylated CpG sites: the spontaneously occurring 
# deamination of a methylated cytosine results in a thymine, and the resulting G:T mismatched bases are often improperly resolved 
# to A:T; whereas the deamination of unmethylated cytosine results in a uracil, which as a foreign base is quickly replaced by a 
# cytosine by the base excision repair mechanism. The C to T transition rate at methylated CpG sites is ~10 fold higher than at unmethylated site
# percentage of C's * percentage of G's * region length =  expected
# vcountPattern("CG") = observed gives us mutation rate???????????????

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
chr22 <- Hsapiens[["chr22"]]
s <- subseq(chr22,start=23456789,width=1000)
print( as.character(s) )
letterFrequency(s,"CG", as.prob = TRUE)
matchPattern("GC",s)






#observed divided by expected is and index for ?????????????


# In the video we briefly described an algorithm that defines CpG islands. 
# This algorithm has been used to create a list that is available in most genomic annotation data bases.

# The Bioconductor AnnotationHub package, which we covered in the Introduction to Bioconductor course, 
# permits us to explore available annotations. By typing the following commands we can see what is available:
  
BiocManager::install("AnnotationHub") 
library(AnnotationHub)
ah <- AnnotationHub()
head(ah)



# We can subset these to just the databases related to the hg19 genome:
  
  
ah <- subset(ah,ah$genome=="hg19")



#We can then use the query function to search the available annotations in this "hub".
#For example, we can search for "gene" annotations by calling:
  
  
query(ah,"genes")


# On the left you see the record IDs for AnnotationHub.
# What is the record ID used by AnnotationHub for hg19 CpG Islands?
# Hint: use query() on the object ah we created above.
  
query(ah, "CpG")

 which(ah$title == "Cpg Islands")
# If the CpG Island record ID from above is stored as the variable cgi_id, 
# we can retrieve the annotations from AnnotationHub with the following command.

cgi <- ah[[cpg_id]]
 

#What is the class of the object cgi (hint: use class function)?
cgi <- ah[["AH5086"]]
class(cgi)
 

#How many CpG islands are represented in the object cgi defined in the previous question?
cgi$name %>% length()
    

# We can extract the sequence of each CpG Island with the following command:
  
  
library(BSgenome.Hsapiens.UCSC.hg19)
cgiseq <- getSeq(Hsapiens, cgi)


# Note that we must verify the same genome builds are being used here:
  
  
genome(cgi)[1:24]  # From the annotationhub
genome(Hsapiens)[1:24] # from the BSgenome


# We will use these sequences to determine the observed to expected ratios at CpG Islands as described in the video.
#Compute the proportion of Cs for each island and report the median of these proportions.
#Hint: use the letterFrequency() function.



s<- subseq(cgiseq) 
class(s)

letterFrequency(s,"C", as.prob = TRUE) %>% median()
 




 
#Compute the proportion of Gs for each island and report the median of these proportions.

letterFrequency(s,"G", as.prob = TRUE) %>% median()


 # CpG islands Assessment Q6
 
#For each CpG island, we now have the proportion of Cs () and the proportion of Gs (). 
# Using these values, we can compute the proportion of CpGs we expect to see by chance 
# if all dinucleotides have the same probability of appearing. To a close approximation this expected proportion is simply .
# The number of CpGs that we expect to see in a genomic interval of size "L" is then probability or porportion of each letter times the length
# proportion_C * proportion_G * Interval Length. 
# Once we have this expectation we can compute the observed to expected ratio.

# Compute the expected number of CpGs in each CpG island using the formula. 
# For each island divide the observed number of CpGs by the expected number of CpGs.
# Hint: use the vcountPattern() function to get the number of CpGs in each island.
 
width(s)
 
df<-letterFrequency(s,"C", as.prob = TRUE) 
df<-as.data.frame(df) 
df$g<-letterFrequency(s,"G", as.prob = TRUE)
colnames(df)[2]<-"G"
df$widths<- width(s) 
df<-mutate(df, expected =  df$C *df$G * df$widths)
 

colnames(df)[1]<-"expected"
df$observed<- vcountPattern("CG",s, algorithm = "auto")

df<- mutate(df, ratio =  observed/expected)


summary(df)
df<-mutate(expected, ratio = observed/expected)

summary(df) 

 
 

#Repeat the previous exercise for GpCs instead of CpGs.
#What is the median observed to expected ratio?

df<-letterFrequency(s,"C", as.prob = TRUE) 
df<-as.data.frame(df) 
df$G<-letterFrequency(s,"G", as.prob = TRUE)
colnames(df)    
df$widths<- width(s) 
df<-mutate(df, expected =  df$C *df$G * df$widths)


colnames(df) 
df$observed_GC<- vcountPattern("GC",s, algorithm = "auto")

df<- mutate(df, ratio =  observed_GC/expected)


summary(df)
 
 
 
 
# Note that the CpG observed to expected ratio is below 1 and that few islands actually surpass a ratio of 1 or more. 
# However, for the rest of the genome, the observed to expected ratio is substantially smaller. 
# To look at regions that are not islands let's shift the islands we have by 20,000 nucleotides.

# To avoid problems with short sequences, we will restrict our analysis to the mapped chromosomes:

    
chr2use <- seqlevels(cgi)[1:24]
index = which( seqnames(cgi) %in% chr2use)

  
# Define the non-CpG islands by shifting the known ones by 20,000 nucleotides.

    
noncgi <- shift(cgi[index],20000)

  
# Some of these regions contain repeats or are unmapped so we remove regions that have 0 Cs or 0 Gs:

    
library(BSgenome.Hsapiens.UCSC.hg19)
noncgiseq <- getSeq(Hsapiens,noncgi)

nullres <- alphabetFrequency(noncgiseq)
keepIndex <- nullres[,"G"]>0 &  nullres[,"C"]>0 & nullres[,"N"]==0
nullres <- nullres[keepIndex,]
noncgiseq <- noncgiseq[keepIndex]

  
# Use nullres and noncgiseq defined above to compute the expected number of CpGs we should see in each of the regions. 
# Report the median observed to expected ratio for these regions.
L2 =  
cprop2 = nullres[,"C"]/L2
gprop2 = nullres[,"G"]/L2
expected2 = L2 * cprop2 * gprop2
observed2 = vcountPattern("CG",noncgiseq)
noncgioe=observed2/expected2
median(noncgioe)
