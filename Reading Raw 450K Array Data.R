# Reading Raw 450K Array Data Assessment Q1
 
#For this assessment you will need to install the following packages:
  
  
BiocManager::install(c("minfi",
                      "IlluminaHumanMethylation450kmanifest",
                      "IlluminaHumanMethylation450kanno.ilmn12.hg19"))


# You will also need to download some raw data from this link External link. https://github.com/genomicsclass/rawdata/tree/master/idats
# These are idat files produced by the Illumina 450K DNA methylation array, 
# which is one of the most widely used technologies for the measurement of DNA methylation. You should download all of the files in this directory 
# to a directory on your machine. (You can also download the entire repository by clicking "Download Zip" here External link. 
# Note that this will be ~100 Mb.)
# The first step to loading the data is to determine the basename of the idat files you have just downloaded. 
# Note that for each sample we have two files: one each for the red and green channels. If you determine the path correctly and type this:
  
path<-"https://github.com/genomicsclass/rawdata/blob/master/idats/"
  
paths<- c(paste0(path,"5775041065_R01C02_Grn.idat"),
paste0(path,"5775041065_R01C02_Red.idat"), 
paste0(path,"5775041065_R04C01_Grn.idat"),
paste0(path,"5775041065_R04C01_Red.idat"),
paste0(path,"5775041068_R01C02_Grn.idat"),
paste0(path,"5775041068_R01C02_Red.idat"),
paste0(path,"5775041068_R04C01_Grn.idat"),
paste0(path,"5775041068_R04C01_Red.idat"),
paste0(path,"5775041068_R06C01_Grn.idat"),
paste0(path,"5775041068_R06C01_Red.idat"),
paste0(path,"5775041084_R01C01_Grn.idat"),
paste0(path,"5775041084_R01C01_Red.idat"))
 
 

path <- "C:/Users/jabba/Downloads/"
dest<- c(paste0(path,"5775041065_R01C02_Grn.idat"),
          paste0(path,"5775041065_R01C02_Red.idat"), 
          paste0(path,"5775041065_R04C01_Grn.idat"),
          paste0(path,"5775041065_R04C01_Red.idat"),
          paste0(path,"5775041068_R01C02_Grn.idat"),
          paste0(path,"5775041068_R01C02_Red.idat"),
          paste0(path,"5775041068_R04C01_Grn.idat"),
          paste0(path,"5775041068_R04C01_Red.idat"),
          paste0(path,"5775041068_R06C01_Grn.idat"),
          paste0(path,"5775041068_R06C01_Red.idat"),
          paste0(path,"5775041084_R01C01_Grn.idat"),
          paste0(path,"5775041084_R01C01_Red.idat"))




for (i in seq_along(paths)) {
  download.file(paths[i], dest[i], mode = "w")
}

download.file("https://github.com/genomicsclass/rawdata/blob/master/idats/targets.csv",
              paste0(path,"targets.csv"))



#you should see a list of idat files and one csv file (targets.csv).

# Let's start by reading in the csv file which contains clinical information. This has one row for each sample and one of the columns includes the "basenames" for the files:
list.files(path, full.names = TRUE)
    
targets <- read.csv("C:/Users/jabba/Downloads/targets.csv", sep = ",")


names(targets)
targets$Basename

  
#How many cancer samples are included in this dataset? Try using the Status column.
targets$Status
 

 
#To make this script work in any working directory we can edit that column to contain the absolute paths.
targets$Basename <- file.path("C:/Users/jabba/Downloads",targets$Basename)


#Then we are ready to read in the raw data with the read.metharray() function:
list.files("C:/Users/jabba/Downloads/", full.names = TRUE)
    
 
library(minfi)
rgset <- read.metharray(targets$Basename,verbose=TRUE)
rownames(targets) <- sampleNames(rgset)
targets <- as(targets, "DataFrame")
pData(rgset) <- targets

  
# We now have the raw data, and we have access to red and green intensities:
dim(getRed(rgset))
dim(getGreen(rgset))

  
# If you are not interested in developing preprocessing algorithms then you can use the built in preprocessing 
# algorithm and go straight to object that give you access to methylation estimates:
mset <- preprocessIllumina(rgset)

  
# This performs the default preprocessing algorithm developed by Illumina. 
# However, for this to be useful we want to have the locations of each CpG and to do that we need map the CpGs to genome.
# minfi keeps this information modular so that when the genome annotation gets updated one can easily change the mapping.
mset <- mapToGenome(mset)

  
#Now we are ready to obtain the methylation values and CpG locations.
dim(getBeta(mset,type="Illumina")) ##the argument type="Illumina" gives us default procedure
head(granges(mset))

  
# If we use the Illumina approach to estimating methylation values, what is the estimated level of methylation 
# for the CpG at location 153807318 on chr4 for sample "5775041068_R04C01"?
which(mset@colData@rownames == "5775041068_R04C01")
mset@assays@data@listData$Meth[123456,4] 
x<-granges(mset)
which(x@ranges@start == 153807318) # the cpg at this location number 123456
x[123456] # this corresponds to cg09689478
y<-mset[123456] # some where in this object for cg09689478 is a methylation value
getBeta(y) # I didn't see a "beta" slot in the object and didn't know this was the command
# or 
i <- which(seqnames(granges(mset))=="chr4" & start(granges(mset))==153807318)
j <- which(rownames(pData(mset))=="5775041068_R04C01")
getBeta(mset,type="Illumina")[i,j]





#Reading Raw 450K Array Data Assessment Q3
 
#Load the bumphunter package:

    
library(bumphunter)

  
#Check the class of the mset object from above:
class(mset)

  
#The bumphunter() function needs one of the following classes:
showMethods("bumphunter")

  
# We need to convert mset to a GenomicRatioSet using the following command:
grset <- ratioConvert(mset,what="beta",type="Illumina")
 
  
#Read the help file for the bumphunter() method:
# Which of the following is the best way of finding DMRs between cancer and normal samples?
help("bumphunter,GenomicRatioSet-method")

bumphunter(grset, model.matrix(~pData(grset)$Status), cutoff=0.1)

bumphunter(mset, pData(grset)$Status, cutoff=0.1)

bumphunter(getBeta(mset), model.matrix(~pData(grset)$Status), cutoff=0.1)

bumphunter(grset)

        
 
#Run bumphunter() as determined in the previous assessment question.
#What is the "area" of the first DMR listed in the DMR table returned by bumphunter()?

x<-bumphunter(grset, model.matrix(~pData(grset)$Status), cutoff=0.1)
x$table 
 
#Reading Raw 450K Array Data Assessment Q5
#The default behvior for bumphunter() is not to smooth the data. Here we will learn how to run bumphunter() with smoothing. However to make the code run faster we will only run it on chr22. To do this, we first subset grset:

    
index <- which(seqnames(grset)=="chr22")
grset2 <- grset[index,]

  
# Now we run bumphunter() without smoothing:
library(bumphunter)
X <- model.matrix(~pData(grset2)$Status)
res <- bumphunter(grset2,X,cutoff=0.25)

  
# To use smoothing we change the smooth argument to TRUE. 
# For details on how to control how the smoothing is done, we refer you to the helpfile:
?bumphunter

  
#Now we run it with smoothing:
res2 <- bumphunter(grset2,X,cutoff=0.25,smooth=TRUE)

  
#Which best characterizes the difference between the DMRs in res and res2?
Hint: the information needed is stored in res$table and res2$table.
res$table
res2$table

#res has more DMRs and res2 has longer DMRs.

 