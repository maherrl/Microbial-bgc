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
setwd("/gpfs/home/rmaher2/crobe/neon")

x <- list.files("./profiles")

lsTab <- list()
for(i in x) {
  
  cov <- read_table(file = paste0("./gene_coverages/",i,"-gene_coverages.txt"), col_types = cols())
  cov[-1] <- sapply(cov[-1], function(x) {x/(sum(x)/1e6)}) # transcripts per million transformation
  func <- read_tsv(file = paste0("./functions/",i,"_functions.txt"), col_types = cols()) %>%
    filter(source == "Pfam") %>% 
    select(-e_value)
  mat <- swap_accession_for_key(cov, func)
  lsTab[[i]] <- mat
  
}

tab_all <- merge_recurse(lsTab, by = c("accession", "function"))

write.csv(tab_all, file = "./results/gene_tab_Pfam.csv", row.names = FALSE)

colnames(tab_all)[-2]
