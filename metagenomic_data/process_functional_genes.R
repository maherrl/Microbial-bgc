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
    select(key = gene_callers_id, accession, `function`) %>% 
    left_join(coverage, by = "key") %>% 
    select(-key) %>% 
    mutate(`function` = word(`function`, sep = "!!!")) %>% 
    mutate(accession = word(accession, sep = "!!!")) %>% 
    drop_na()%>%
    group_by(accession, `function`) %>% 
    summarize_all(sum)
}


#####
## Analysis
setwd("/gpfs/home/rmaher2/crobe/neon/surface/genes")

sample_names <- read.table(file = "../databases/names.txt")
x <- sample_names$V1

lsTab <- list()
for(i in x) {
  
  cov <- read_table(file = paste0(i,"_gene_table-GENE-COVERAGES.txt"), col_types = cols())
  func <- read_tsv(file = paste0(i,"_functions.txt"), col_types = cols()) %>%
    filter(source == "Pfam") %>% 
    select(-e_value)
  mat <- swap_accession_for_key(cov, func)
  lsTab[[i]] <- mat
  
}

tab_all <- merge_recurse(lsTab, by = c("accession", "function"))

write.csv(tab_all, file = "./Surface_gene_by_sample_Pfam_coverages.csv", row.names = FALSE)
