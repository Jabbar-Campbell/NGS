
#We will look for differential exon usage in the same experimental data as shown in the video, 
#but we will use a different subset of the genes (to speed up the time required for analysis).

#Build a DEXSeq dataset object
BiocManager::install("pasilla")
library("pasilla")
inDir <- system.file("extdata", package="pasilla", mustWork=TRUE)
countFiles <- list.files(inDir, pattern="fb.txt$", full.names=TRUE)    
flattenedFile <- list.files(inDir, pattern="gff$", full.names=TRUE) 

sampleTable <- data.frame(
  row.names = c("treated1", "treated2", "treated3",
                "untreated1", "untreated2", "untreated3", "untreated4"), 
  condition = c("knockdown", "knockdown", "knockdown", 
                "control", "control", "control", "control"))                 

library("DEXSeq")  
dxd <- DEXSeqDataSetFromHTSeq(countFiles, sampleData=sampleTable, 
                              design= ~ sample + exon + condition:exon, 
                              flattenedfile=flattenedFile)
Subset
# Now we will subset to 1000 genes on chr2L which do not have low counts 
# (just for demonstration, typically you would just run the analysis on the whole dataset):
  
rmean <- rowMeans(counts(dxd)) 
# We use rowRanges() to pull out chr2L:
  
dxd2L <- dxd[seqnames(rowRanges(dxd)) == "chr2L" & rmean > 10,] 
  
#Now subset to first 1000:
  
dxd2L <- dxd2L[1:1000,] 

Exon usage Q1
1 point possible (graded)
#What is the name of the gene containing the exon with the smallest adjusted p-value for differntial exon usage?
#The gene name is the part starting with "FBgn" before the ':' and the exon number, E001. So, stop before the ':'
ranges(dxd2L)
dxd2L[1,]
x<-DEXSeq(dxd2L)
  



#Make a DEXSeq plot of the DEXSeqResults object for this gene.
esf <- estimateSizeFactors( dxd2L )
esd <- estimateDispersions(esf)
tfd<- testForDEU(esd)
efc <- estimateExonFoldChanges( tfd )
dsr<-DEXSeqResults(efc)

library("tidyr")
x<-which(dsr$featureID == "E009")
c(x)
c(30,48,63)


plotMA(dsr[which(dsr$featureID == "E009"),]  , 
       foldChangeColumn = "log2fold_knockdown_control", 
       norCounts=TRUE, 
       displayTranscripts=TRUE,
       cex = 2)

plotDEXSeq( dsr, "FBgn0000256",norCounts=TRUE, displayTranscripts=TRUE )

Use the settings: norCounts=TRUE, displayTranscripts=TRUE


# You should see two exonic parts with differential exon usage (colored pink). 
# Exonic part E009 was the one with the smallest p-value in the dataset,
# consistently knocked down in control samples.


#The transcripts are shown in the white boxes below the plot 
#(the grey boxes just show all the exonic parts). 
#Make the plot as wide as possible to help see the exonic parts. 
#(In RStudio, click the zoom button)

How many transcripts does exonic part E009 appear in?