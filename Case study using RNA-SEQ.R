# for issues regarding updates vist...
https://stackoverflow.com/questions/41839214/installation-path-not-writable-r-unable-to-update-packages

# Scripts for visualization can be found here.....
https://github.com/genomicsclass/labs External link

# Written Tutorial can be found here...
http://genomicsclass.github.io/book/ External link


if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GenomicRanges","Biostrings","AnnotationHub"))

