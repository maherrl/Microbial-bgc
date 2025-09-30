## script for processing Anvio gene functions and coverage data into a matrix format

library("tidyverse")
library("plyr")
library("dplyr")
library("reshape")

#############################
## functions
############################

swap_accession_for_key <- function(coverage, functions) {
  functions %>% 
    select(gene_callers_id, accession, `function`) %>% 
    left_join(coverage, by = "gene_callers_id") %>% 
    select(-gene_callers_id) %>% 
    mutate(`function` = word(`function`, sep = "!!!")) %>% 
    mutate(accession = word(accession, sep = "!!!")) %>% 
    drop_na()%>%
    group_by(accession, `function`) %>% 
    summarize_all(sum)
}


#####
## Analysis
counts <- read.csv(file = "./data/reads_per_sample.csv")
x <- list.files("./data/blitz")

lsTab <- list()
for(i in x) {
  
  # transform the read counts to RPKM from Konwar et al.
  reads <- read.table(file = paste0("./data/blitz/", i), header = TRUE)
  samp <- sub('\\.txt$', '', i)
  RPK <- reads$num_mapped_reads / (reads$length / 1000)
  samp_reads <- counts[counts$sample == samp, 2]
  scal_fac <- samp_reads / 1e6
  reads$RPKM <- RPK / scal_fac
  reads <- reads %>% select(gene_callers_id, RPKM)
  names(reads)[names(reads) == 'RPKM'] <- samp
  
  # get functional data
  samp2 <- sub('\\.[0-9]$', '', samp) # for benthic replicates
  func <- read_tsv(file = paste0("./data/functions/",samp2,"_functions.txt"), col_types = cols()) %>%
    filter(source == "COG20_CATEGORY") %>% 
    select(-e_value)
  mat <- swap_accession_for_key(reads, func)
  lsTab[[i]] <- mat
  
}

tab_all <- merge_recurse(lsTab, by = c("accession", "function"))

write.csv(tab_all, file = "./results/gene_tab_COG20_CATEGORY.csv", row.names = FALSE)

colnames(tab_all)[-2]
