library(phyloseq)
library(dplyr)
library(ggplot2)
library(vegan)
library(cowplot)
library(RColorBrewer)
library(tidyverse)
library(wrMisc)


#### BENTHIC 16S data ordinations ####

# ### 1. Prepare environmental table for ordinations ### 
# # (If using the prepared "envtab_benthic.rds" from the box, jump to item 2 below)
# 
# envtab <- read.csv("./EnvironmentalData-Edit.csv") # this table can be replaced for a newer version for the RDA, 
# # but for NMDS and PCOA it won't make a difference since we'll use only Site, Substrate type and Water Quality groups.
# 
# # We need to put the sample names as written in the seq table as row names in the environmental table
# # Extract sample names from exportfile
# exportfile_benthic <- read.csv("./exportfile_Benthic16S_noLakes_Nov2018-Dec2019_filledNAs.csv")
# exportfile_benthic$sample.names <- sapply(strsplit(basename(exportfile_benthic$rawDataFilePath), "_"), `[`, 2)
# exportfile_benthic$sample.names
# 
# # Inner join of env table and exportfile table to match sample names with dnaSampleID
# envtab_benthic <- inner_join(exportfile_benthic, envtab, by = "dnaSampleID") # 980 rows
# envtab_benthic <- subset(envtab_benthic, select = -rawDataFilePath) # remove column with raw data file path (because it is unique for each paired end sequence (R1 and R2) and we don't want 2 rows per sample)
# envtab_benthic <- distinct(envtab_benthic)
# 
# setdiff(envtab_benthic$dnaSampleID, exportfile_benthic$dnaSampleID)
# setdiff(exportfile_benthic$dnaSampleID, envtab_benthic$dnaSampleID)
# 
# rownames(envtab_benthic) <- envtab_benthic$sample.names
# saveRDS(envtab_benthic, "./envtab_benthic.rds")

### 2. Import seqtab, taxtab, and environmental table ###

seqtab_bact_arch_benthic_new <- readRDS("./seqtab_bact_arch_benthic_new.rds")
taxtab_bact_arch_benthic_new <- readRDS("./taxtab_bact_arch_benthic_new.rds")
envtab_benthic <- readRDS("./envtab_benthic.rds") # this table can be replaced for a newer version, 
# but for NMDS and PCOA it doesn't matter since we'll use only Site, Substrate type and Water Quality groups, which won't change.

# Remove ASVs that were not assigned to the Kingdom level (2567 ASVs)
seqtab_bact_arch_benthic_new <- seqtab_bact_arch_benthic_new[,-which(is.na(taxtab_bact_arch_benthic_new$Kingdom))]
taxtab_bact_arch_benthic_new <- taxtab_bact_arch_benthic_new[-which(is.na(taxtab_bact_arch_benthic_new$Kingdom)),]

# some more prep:
envtab_benthic$aquMicrobeType.x[envtab_benthic$aquMicrobeType.x == "epilithon_largeSubstrate"] <- "epilithon"
envtab_benthic$aquMicrobeType.x[is.na(envtab_benthic$aquMicrobeType.x)] <- "epilithon"
colnames(envtab_benthic)[colnames(envtab_benthic) == "aquMicrobeType.x"] <- "Type"
unique(envtab_benthic$Type)
envtab_benthic$WQGroup <- as.factor(envtab_benthic$WQGroup)

seqtab_bact_arch_benthic_new <- t(seqtab_bact_arch_benthic_new) # now taxa are rows
seqtab_bact_arch_benthic_new <- as.matrix(seqtab_bact_arch_benthic_new)
taxtab_bact_arch_benthic_new <- as.matrix(taxtab_bact_arch_benthic_new)

### 3. Make phyloseq objects ###

psb <- phyloseq(otu_table(seqtab_bact_arch_benthic_new, taxa_are_rows = TRUE),
                tax_table(taxtab_bact_arch_benthic_new),
                sample_data(envtab_benthic))
psb # 429 samples, 63813 ASVs

# Subset to only Epilithon and Epipsammon
psb_epilit <- subset_samples(psb, Type == "epilithon" | Type == "epipsammon" )
psb_epilit # 306 samples

# set colors

col <- c("ARIK" = "pink1", "BIGC" =  "slateblue1", "BLDE" = "purple3",
         "BLUE" = "turquoise2", "CARI" =  "steelblue", "COMO" = "navyblue",
         "CUPE" = "orange", "GUIL" = "tomato", "HOPB" = "palevioletred",
         "KING" = "violetred", "LECO" = "mediumpurple1",  "LEWI" = "coral2", 
         "MART" = "yellowgreen", "MAYF" = "red2", "MCDI" = "springgreen2",
         "MCRA" = "violet", "OSKR" = "skyblue", "POSE" = "grey70", 
         "PRIN" = "grey50", "REDB" = "palegreen4", "SYCA" = "tan", 
         "TECR" = "grey30", "WALK" = "blue2", "WLOU" = "brown")

cbPalette <- c("1" = "#999999", "2" = "#F0E442", "3" = "#E69F00", "5" = "#56B4E9", "6" = "#009E73")


### 4. ORDINATIONS WITH THE PHYLOSEQ OBJECT ###

## NMDS

# all benthic samples
# set.seed(1994)
# neon.ord <- ordinate(psb, "NMDS", "bray", trymax = 50)
# neon.ord$stress # 0.1839158
# 
# plot_ordination(psb, neon.ord, color = "Site", shape = "Type", 
#                 title = "NMDS of 16S rRNA gene for benthic samples (stress = 0.18)") +
#   theme_bw()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   theme(axis.title.x=element_text(size=16,color="black"))+
#   theme(axis.title.y=element_text(size=16,color="black"))+
#   theme(axis.text.y=element_text(size=14,color="black"))+
#   theme(axis.text.x=element_text(size=14,color="black")) +
#   scale_color_manual(values = col) +
#   geom_point(size = 3) +
#   guides(color = guide_legend(order = 1))
#   #scale_shape_discrete(labels=c("EPILITHON", "EPIPELON", "EPIPHYTON", "EPIPSAMMON", "EPIXYLON"))
# ggsave("ProcessedDataForAnalysis/Figures/nmds_bray_16s_site_type.jpg", dpi = 300, width = 7.5, height = 6, units = "in")
# 
# plot_ordination(psb, neon.ord, color = "WQGroup", shape = "Type", 
#                 title = "NMDS of 16S rRNA gene for benthic samples (stress = 0.18)") +
#   theme_bw()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   theme(axis.title.x=element_text(size=16,color="black"))+
#   theme(axis.title.y=element_text(size=16,color="black"))+
#   theme(axis.text.y=element_text(size=14,color="black"))+
#   theme(axis.text.x=element_text(size=14,color="black")) +
#   scale_color_manual(values = cbPalette) +
#   geom_point(size = 3) +
#   guides(color = guide_legend(order = 1))
#   #scale_shape_discrete(labels=c("EPILITHON", "EPIPELON", "EPIPHYTON", "EPIPSAMMON", "EPIXYLON"))
# ggsave("ProcessedDataForAnalysis/Figures/nmds_bray_16s_group_type.jpg", dpi = 300, width = 7.5, height = 6, units = "in")

# epilithon and epipsammon only
set.seed(1994)
neon.ord2 <- ordinate(psb_epilit, "NMDS", "bray", trymax = 50) # no convergence

plot1a=plot_ordination(psb_epilit, neon.ord2, color = "Site", shape = "Type",
                title = "NMDS of 16S rRNA gene for benthic samples (stress = 0.18)") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black")) +
  scale_color_manual(values = col) +
  scale_shape_manual(values=c(16, 3))+
  geom_point(size = 3)
ggsave("ProcessedDataForAnalysis/Figures/nmds_bray_16s_epilithon_site.jpg", dpi = 300, width = 7.5, height = 6, units = "in")

plot1b=plot_ordination(psb_epilit, neon.ord2, color = "WQGroup", shape = "Type",
                title = "NMDS of 16S rRNA gene for benthic samples (stress = 0.18)") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black")) +
  scale_color_manual(values = cbPalette) +
  scale_shape_manual(values=c(16, 3))+
  geom_point(size = 3)
ggsave("ProcessedDataForAnalysis/Figures/nmds_bray_16s_epilithon_group.jpg", dpi = 300, width = 7.5, height = 6, units = "in")


plot_grid(plot1a, plot1b, nrow=1, labels=c("A", "B"), align="hv")


# ## PCoA
# 
# # all benthic samples
# set.seed(1994)
# neon.ord.pcoa <- ordinate(psb, "PCoA", "bray", trymax = 50)
# 
# plot_ordination(psb, neon.ord.pcoa, color = "Site", shape = "Type", 
#                 title = "PCoA of 16S rRNA gene for benthic samples") +
#   scale_color_manual(values = col) +
#   geom_point(size = 3) + 
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   theme(axis.title.x=element_text(size=16,color="black"))+
#   theme(axis.title.y=element_text(size=16,color="black"))+
#   theme(axis.text.y=element_text(size=14,color="black"))+
#   theme(axis.text.x=element_text(size=14,color="black")) +
#   guides(color = guide_legend(order = 1))
#   #scale_shape_discrete(labels=c("EPILITHON", "EPIPELON", "EPIPHYTON", "EPIPSAMMON", "EPIXYLON"))
# ggsave("ProcessedDataForAnalysis/Figures/pcoa_bray_16s_site_type.jpg", dpi = 300, width = 7.5, height = 6, units = "in")
# 
# plot_ordination(psb, neon.ord.pcoa, color = "Group", shape = "Type", 
#                 title = "PCoA of 16S rRNA gene for benthic samples") +
#   scale_color_manual(values = cbPalette) +
#   geom_point(size = 3) + 
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   theme(axis.title.x=element_text(size=16,color="black"))+
#   theme(axis.title.y=element_text(size=16,color="black"))+
#   theme(axis.text.y=element_text(size=14,color="black"))+
#   theme(axis.text.x=element_text(size=14,color="black")) +
#   guides(color = guide_legend(order = 1))
#   #scale_shape_discrete(labels=c("EPILITHON", "EPIPELON", "EPIPHYTON", "EPIPSAMMON", "EPIXYLON"))
# ggsave("ProcessedDataForAnalysis/Figures/pcoa_bray_16s_group_type.jpg", dpi = 300, width = 7.5, height = 6, units = "in")
# 
# # epilithon only
# set.seed(1994)
# neon.ord.pcoa2 <- ordinate(psb_epilit, "PCoA", "bray", trymax = 50)
# 
# plot_ordination(psb_epilit, neon.ord.pcoa2, color = "Site", 
#                 title = "PCoA of 16S rRNA gene for epilithon samples") +
#   scale_color_manual(values = col) +
#   geom_point(size = 3) + 
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   theme(axis.title.x=element_text(size=16,color="black"))+
#   theme(axis.title.y=element_text(size=16,color="black"))+
#   theme(axis.text.y=element_text(size=14,color="black"))+
#   theme(axis.text.x=element_text(size=14,color="black"))
# ggsave("ProcessedDataForAnalysis/Figures/pcoa_bray_16s_epilithon_site.jpg", dpi = 300, width = 7.5, height = 6, units = "in")
# 
# plot_ordination(psb_epilit, neon.ord.pcoa2, color = "Group", 
#                 title = "PCoA of 16S rRNA gene for epilithon samples") +
#   scale_color_manual(values = cbPalette) +
#   geom_point(size = 3) + 
#   theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   theme(axis.title.x=element_text(size=16,color="black"))+
#   theme(axis.title.y=element_text(size=16,color="black"))+
#   theme(axis.text.y=element_text(size=14,color="black"))+
#   theme(axis.text.x=element_text(size=14,color="black"))
# ggsave("ProcessedDataForAnalysis/Figures/pcoa_bray_16s_epilithon_group.jpg", dpi = 300, width = 7.5, height = 6, units = "in")
# 

## PERMANOVAs using function adonis2 (vegan package)
sampledf <- data.frame(sample_data(psb))
dist <- distance(psb, method = "bray")

adonis2(dist ~ Type, sampledf) # clustering by substrate types
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist ~ Type, data = sampledf)
#            Df SumOfSqs    R2      F   Pr(>F)    
# Type       4   12.541 0.06384 7.2283  0.001 ***
# Residual 424  183.912 0.93616                  
# Total    428  196.453 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis2(dist ~ Group, sampledf) # clustering by water quality groups
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist ~ Group, data = sampledf)
#            Df SumOfSqs   R2      F    Pr(>F)    
# Group      4   17.037 0.08672 10.065  0.001 ***
# Residual 424  179.416 0.91328                  
# Total    428  196.453 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis2(dist ~ Site, sampledf) # clustering by sites
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist ~ Site, data = sampledf)
#            Df SumOfSqs   R2      F    Pr(>F)    
# Site      23   61.776 0.31446 8.0771  0.001 ***
# Residual 405  134.677 0.68554                  
# Total    428  196.453 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis2(dist ~ Group*Type, sampledf) # clustering by water quality groups and substrate types combined
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Blocks:  strata 
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist ~ Group * Type, data = sampledf, strata = sampledf$Site)
#             Df   SumOfSqs   R2       F    Pr(>F)    
# Group        4   17.005 0.08654 11.1041  0.001 ***
# Type         4   11.512 0.05858  7.5171  0.001 ***
# Group:Type   6    9.488 0.04829  4.1305  0.001 ***
# Residual   414  158.497 0.80660                   
# Total      428  196.502 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis2(dist ~ Site*Type, sampledf) # clustering by sites and substrate types combined
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist ~ Site * Type, data = sampledf)
#            Df   SumOfSqs  R2      F    Pr(>F)    
# Site       23   61.691 0.31395 8.9267  0.001 ***
# Type        4    6.813 0.03467 5.6685  0.001 ***
# Site:Type  10   10.514 0.05351 3.4992  0.001 ***
# Residual  391  117.484 0.59788                  
# Total     428  196.502 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## PERMUTEST: testing groups dispersion
bd <- betadisper(dist, sampledf$Type)
permutest(bd) # p = 0.001
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#             Df  Sum Sq  Mean Sq   F     N.Perm Pr(>F)    
# Groups      4 0.50801 0.127002 90.885    999  0.001 ***
# Residuals 424 0.59249 0.001397                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

bd2 <- betadisper(dist, sampledf$Group)
permutest(bd2)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#             Df  Sum Sq  Mean Sq   F     N.Perm Pr(>F)    
# Groups      4 0.15088 0.037720 19.724    999  0.001 ***
# Residuals 424 0.81085 0.001912                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

bd3 <- betadisper(dist, sampledf$Site)
permutest(bd3)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#            Df  Sum Sq  Mean Sq    F   N.Perm Pr(>F)    
# Groups     23 0.89844 0.039063 7.077    999  0.001 ***
# Residuals 405 2.23548 0.005520                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#### BENTHIC Metagenomics data ordination ####
setwd("/Users/brittnibertolet/Documents/BertoletLab/ResearchProjects/Microbial-BGC/")
gen.tab <- read.csv("Data/Metagenomic Process Data/gene_tab_KOfam.csv", header = TRUE)
gen.tab <- read.csv("Data/Metagenomic Process Data/gene_tab_KEGG_Class.csv", header = TRUE)

## MAKING A PHYLOSEQ OBJECT
# make taxonomy table
# for gene_tab_KOfam.csv
gen.tax <- as.matrix(gen.tab[,1:2])
rownames(gen.tax) <- gen.tab$accession
#gen.tax[gen.tax == "Module set"] <- "Metabolic capacity"
#colnames(gen.tax) <- c("function")
TAX <- tax_table(gen.tax)
#gen.tax <- matrix(gen.tab[,2])
#rownames(gen.tax) <- gen.tab$accession

# for gene_tab_KEGG_Class
gen.tax <- gen.tab %>% separate(function., c('Module_1','Module_2','Module_3'), sep = "; ")
gen.tax <- as.matrix(gen.tax[,1:4])
rownames(gen.tax) <- gen.tab$accession
colnames(gen.tax) <- c("Kingdom", "Phylum", "Class", "Order")
TAX <- tax_table(gen.tax)



# make otu (gene) table
gen.otu <- as.matrix(gen.tab[,3:ncol(gen.tab)])
rownames(gen.otu) <- gen.tab$accession
gen.otu[is.na(gen.otu)] <- 0
OTU <- otu_table(gen.otu, taxa_are_rows = TRUE)

# make metadata
gen.meta <- read.csv("Data/Metagenomic Process Data/metag_map.csv", header = TRUE, row.names = 1)
gen.meta$Group <- as.character(gen.meta$Group)
gen.meta$Site <- as.factor(gen.meta$Site)
gen.meta$Site_Type <- factor(gen.meta$Site_Type, levels = c("ARIK_EPIPSAMMON", "REDB_EPIPSAMMON", "WLOU_EPIPSAMMON", 
                                                            "BIGC_EPIPSAMMON", "HOPB_EPIPSAMMON", "COMO_EPIPSAMMON", 
                                                            "SYCA_EPIPSAMMON", "CARI_EPILITHON", "KING_EPILITHON", 
                                                            "MCDI_EPILITHON", "BIGC_EPILITHON", "BLUE_EPILITHON", 
                                                            "CUPE_EPILITHON", "HOPB_EPILITHON", "OKSR_EPILITHON", 
                                                            "COMO_EPILITHON", "SYCA_EPILITHON", "BLDE_EPILITHON", 
                                                            "GUIL_EPILITHON", "REDB_EPILITHON", "WLOU_EPILITHON", 
                                                            "MAYF_EPIXYLON", "BLDE_SS", "COMO_SS", "GUIL_SS", 
                                                            "SYCA_SS", "CARI_SS", "KING_SS", "REDB_SS", "TECR_SS", 
                                                            "WLOU_SS", "BLUE_SS", "CUPE_SS", "LECO_SS", "OKSR_SS", 
                                                            "TOMB_C0", "BLWA_C0", "HOPB_SS"))
sampledata <- sample_data(gen.meta)

# make phyloseq object
physeq = phyloseq(OTU, TAX, sampledata)
physeq

#saveRDS(physeq, file = "./data/Physeq_KOfam.rds")
#saveRDS(physeq, file = "./data/Physeq_KEGG_Module.rds")
#physeq <- Genes_phyloseq_obj_trans

## CALCULATE ALPHA DIVERSITY FROM GENE TABLE
# alpha diversity
H <- diversity(t(gen.otu))
simp <- diversity(t(gen.otu), "simpson")
invsimp <- diversity(t(gen.otu), "inv")
S <- specnumber(t(gen.otu))
J <- H/log(S)
alpha <- mergeVectors(Shannon = H, Simpson = simp, Invsimpson = invsimp, Richness = S, Evenness = J)
alpha <- as.data.frame(t(alpha))
stats <- merge(gen.meta, alpha, by = "row.names")
stats_ben <- stats[stats$Source == "benthic",] # subset to benthic
stats_ben <- stats_ben[stats_ben$Site != "OKSR",] # exclude OKSR

cbPalette <- c("#009E73", "#E69F00", "#56B4E9")
ggplot(data = stats_ben, aes(x = Site, y = Invsimpson, color = Type)) + 
  geom_boxplot() +
  geom_point() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_color_manual(values = cbPalette) +
  ylab("Inverse simpson")+ggtitle("KOfam functions")


# ORDINATIONS WITH THE PHYLOSEQ OBJECT
# subset to benthic and exclude OKSR
physeq_ben <- subset_samples(physeq, Source == "benthic") %>%
  subset_samples(Site != "OKSR")

col <- c("ARIK" = "pink1", "BIGC" =  "slateblue1", "BLDE" = "purple3",
         "BLUE" = "turquoise2", "CARI" =  "steelblue", "COMO" = "navyblue",
         "CUPE" = "orange", "GUIL" = "tomato", "HOPB" = "palevioletred",
         "KING" = "violetred", "MAYF" = "red2", "MCDI" = "springgreen2",
         "REDB" = "palegreen4", "SYCA" = "tan", "WLOU" = "brown")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# NMDS
set.seed(1994)
neon.ord <- ordinate(physeq_ben, "NMDS", "bray", trymax = 50)
neon.ord$stress
plot1c=plot_ordination(physeq_ben, neon.ord, color = "Site", shape = "Type", 
                title = "NMDS of KEGG Modules for benthic samples (stress = 0.10)") +
  scale_color_manual(values = col) +
  geom_point(size = 3) + 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black"))+
  scale_shape_manual(values=c(16, 3))
  

# Subset to only Epilithon and Epipsammon
physeq_ben_epil <- subset_samples(physeq_ben, Type == "EPILITHON" | Type == "EPIPSAMMON" )
physeq_ben_epil # 306 samples

set.seed(1993)
neon.ord <- ordinate(physeq_ben_epil, "NMDS", "bray", trymax = 50)
neon.ord$stress
plot1d=plot_ordination(physeq_ben_epil, neon.ord, color = "Group", shape="Type", 
                title = "NMDS of KEGG Modules benthic samples (stress = 0.16)") +
  scale_color_manual(values = cbPalette) +
  geom_point(size = 3) +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black"))+
  scale_shape_manual(values=c(16, 3))



plot_grid(plot1a, plot1b, plot1c, plot1d, align="hv", labels=c("A", "B", "C","D"))
ggsave("Figure1.pdf", height=10, width=15)


 