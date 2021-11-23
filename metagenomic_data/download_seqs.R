##
# Script to download NEON metagenomes

rm(list = ls())
# load packages
library(neonUtilities)
library(raster)

#one way to download metadata and compile it into a stacked tables folder
# this requires you to download data from the website 
stackByTable("/Users/aho/Downloads/NEON_seq-metagenomic-microbe-benthic")

# import metadata into an R object
metadata <- loadByProduct(dpID = 'DP1.20279.001',
                          #startdate = "2016-09", enddate = "2016-09",
                          #site = "MAYF",
                          check.size = FALSE, package = 'expanded')

# change to your specific RawDataFiles name
rawFileInfo <- metadata$mms_benthicRawDataFiles

# get unique file paths to each sample
download.dir <- paste0(getwd(), "/raw_metagenomes") #started 5:22pm
already.have <- list.files(download.dir, pattern = ".fastq", recursive = T)
already.have <- gsub("_R[12].fastq|_R[12].fastq.gz", "", already.have)
already.have <- basename(already.have)
#not.have <- rawFileInfo[!(rawFileInfo$dnaSampleID %in% already.have | rawFileInfo$internalLabID %in% already.have),]
not.have <- rawFileInfo[!rawFileInfo$dnaSampleID %in% already.have,]
not.have <- not.have[!not.have$internalLabID %in% already.have,]

url_list <- unique(not.have$rawDataFilePath) 

# function to download files associated with each URL
lapply(1:length(url_list), function(i) {
  download.file(url_list[i],  destfile = paste0(download.dir, basename(url_list[i])))
})

# alternatively you can export the url list and run a for loop with wget on the command line.
