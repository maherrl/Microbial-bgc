## NEON 16S analyses
## Rebecca L. Maher
## 07/07/2023
rm(list = ls())

library(phyloseq)
library(dplyr)
library(ggplot2)
library(vegan)
library(stringr)
library(tidyverse)



### 1. Import seqtab, taxtab, and environmental table ###

seqtab_bact_arch_benthic_new <- readRDS("data/16S/seqtab_bact_arch_benthic_new.rds")
taxtab_bact_arch_benthic_new <- readRDS("data/16S/taxtab_bact_arch_benthic_new.rds")
envtab_benthic <- readRDS("data/16S/envtab_benthic.rds")
#envtab_benthic["WQGroup"][envtab_benthic["WQGroup"] == "4"] <- "3" # fix error in env data
#saveRDS(envtab_benthic, file = "data/16S/envtab_benthic.rds")
# Remove ASVs that were not assigned to the Kingdom level (2567 ASVs)
seqtab_bact_arch_benthic_new <- seqtab_bact_arch_benthic_new[,-which(is.na(taxtab_bact_arch_benthic_new$Kingdom))]
taxtab_bact_arch_benthic_new <- taxtab_bact_arch_benthic_new[-which(is.na(taxtab_bact_arch_benthic_new$Kingdom)),]
# transform into format for phyloseq
seqtab_bact_arch_benthic_new <- t(seqtab_bact_arch_benthic_new) # now taxa are rows
seqtab_bact_arch_benthic_new <- as.matrix(seqtab_bact_arch_benthic_new)
taxtab_bact_arch_benthic_new <- as.matrix(taxtab_bact_arch_benthic_new)



### 2. New phyloseq object
psb <- phyloseq(otu_table(seqtab_bact_arch_benthic_new, taxa_are_rows = TRUE),
                tax_table(taxtab_bact_arch_benthic_new),
                sample_data(envtab_benthic))
psb # 395 samples, 63813 ASVs
psb <- prune_samples(sample_sums(psb)>2000, psb)
psb


### 3. Performing rarefaction to get the average Bray Curtis Dissimilarity matrix from
###.   1,000 rarefied tables
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
} # function phyloseq object -> vegan

psb.veg <- psotu2veg(psb)

# Average Bray Curtis Dissimilarity
min_depth <- min(sample_sums(psb))
bc.dist <- avgdist(psb.veg, sample = min_depth, dmethod = "bray")
#write.table(bc.dist, file= "data/16S/Bray_Curtis_Dist.txt")
bc.dist <- read.table(file= "data/16S/Bray_Curtis_Dist.txt")



### 4. Alpha Diversity metrics averaged from 1,000 rarefied tables
# Observed richness
S <- specnumber(psb.veg)
Srare <- rarefy(psb.veg, min_depth, se = TRUE)
plot(S, Srare.t$S, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(psb.veg, step = 20, sample = min_depth, col = "blue", cex = 0.6)
Srare <- as.data.frame(t(Srare))

# Shannon
num_rarefactions <- 1000
nsamp <- nrow(psb.veg)
shannon_diversity <- matrix(nrow = nsamp, ncol = num_rarefactions)
rownames(shannon_diversity) <- rownames(psb.veg)

for (i in 1:num_rarefactions) {
  rarefied_matrix <- rrarefy(psb.veg, sample = min_depth)
  shannon_diversity[,i] <- diversity(rarefied_matrix, index = "shannon")
}
average_shannon_diversity <- rowMeans(shannon_diversity)

# Simpson (1-D)
simpson_diversity <- matrix(nrow = nsamp, ncol = num_rarefactions)
rownames(simpson_diversity) <- rownames(psb.veg)

for (i in 1:num_rarefactions) {
  rarefied_matrix <- rrarefy(psb.veg, sample = min_depth)
  simpson_diversity[,i] <- diversity(rarefied_matrix, index = "simpson")
}
average_simpson_diversity <- rowMeans(simpson_diversity)

# Invsimpson (1/D)
invsimp_diversity <- matrix(nrow = nsamp, ncol = num_rarefactions)
rownames(invsimp_diversity) <- rownames(psb.veg)

for (i in 1:num_rarefactions) {
  rarefied_matrix <- rrarefy(psb.veg, sample = min_depth)
  invsimp_diversity[,i] <- diversity(rarefied_matrix, index = "invsimpson")
}
average_invsimp_diversity <- rowMeans(invsimp_diversity)

alpha <- as.data.frame(cbind(Srare$S, average_shannon_diversity, average_simpson_diversity, average_invsimp_diversity))
colnames(alpha) <- c("obs", "shannon", "simpson","invsimp")
alpha$sample.id <- rownames(alpha)
# function
# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}
env_meta <- pssd2veg(psb)
alpha <- merge(alpha, env_meta, by.x = "sample.id", by.y = "sample.names", all = TRUE)
#write.csv(alpha, file = "data/16S/alpha_div.csv")



## 5. Aggregating samples by SITE and TIME to match environmental data
# Note this action deletes siteID and time variables from the environmental data
# If you want to use these variables in the future (ie for plotting as in #6 below)
# you have to add these variables back
# epilithon only
psb <- subset_samples(psb, aquMicrobeType.y == "epilithon")
psb_epil
nitr.epil <- subset_samples(nitr, Type =="EPILITHON")
psb <- nitr
# get meta data
env_meta <- pssd2veg(psb)
env_meta$env_grp <- paste(env_meta$siteID,env_meta$CollectDate, sep = "_")
env_meta$env_grp <- paste(env_meta$Site,env_meta$Date, sep = "_")
# move sample data from vegan to phyloseq
sample_data(psb) <- as.data.frame(env_meta)
# merging
mergedPSB <- merge_samples(psb, "env_grp", fun = sum) # counts are summed by SITE_TIME combinations
mergedPSB
#saveRDS(mergedPSB, "data/16S/seqtab_merged.rds")
mergedPSB <- readRDS(file = "data/16S/seqtab_epip_merged.rds")

# rarefying to get Bray Curtis distance matrix
min_depth <- min(sample_sums(mergedPSB))
mergedPSB.veg <- psotu2veg(mergedPSB)
merged.bc.dist <- avgdist(mergedPSB.veg, sample = min_depth, dmethod = "bray") # avg dist matrix from 1,000 rarefied tables
merged.bc.dist <- vegdist(mergedPSB.veg, dmethod = "bray")
#saveRDS(merged.bc.dist, file = "data/16S/Merged_Bray_Curtis_Dist.rds")
merged.bc.dist <- readRDS(file = "data/16S/Merged_Bray_Curtis_Dist_epip.rds")


## 6. NMDS with envfit
merged.nmds <- metaMDSiter(merged.bc.dist, k=2, trymax = 1000, maxit = 1000, autotransform = FALSe)
mergedPSB <- subset_samples(mergedPSB, Type == "EPILITHON")
merged.cca <- ordinate(mergedPSB, "CCA")
env_merged <- pssd2veg(mergedPSB)
env_merged <- env_merged[,8:ncol(env_merged)]
env_merged <- env_merged[,-ncol(env_merged)]
env_merged$site_date <- rownames(env_merged)
env_merged$site_date_type <- rownames(env_merged)
env_merged <- within(env_merged, rm(site_date))

en = envfit(.cca, env_merged, permutations = 999, na.rm = TRUE)
en

ef.adj <- en 
pvals.adj <- p.adjust (en$vectors$pvals, method = 'fdr')
ef.adj$vectors$pvals <- pvals.adj
ef.adj

plot(merged.cca)
plot(en)
# extract nmds score
data.scores = as.data.frame(scores(merged.cca))
# add 'site' column
data.scores$Site = env_merged$Site
#data.scores$Type = env_merged$Type
head(data.scores)
# add 'WQGroup' column
wq_group <- unique(env_meta[, c('Site',"Type")])
data.scores <- inner_join(data.scores, wq_group, by = join_by(Site))
data.scores$WQGroup <- as.character(data.scores$WQGroup)
data.scores$Type <- as.character(data.scores$Type)
#data.scores$Group = env_nmds$WQGroup
#data.scores$Group = as.character(data.scores$Group)
# extract segments for each environmental variable
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cont$Site <- rownames(en_coord_cont)
# only significant pvalues
A <- as.list(en$vectors)
pvals <- as.data.frame(A$pvals)
arrows<-as.data.frame(A$arrows*sqrt(A$r))
C<-cbind(arrows, pvals)
#subset
Cred<-subset(C,pvals<0.05)
Cred <- cbind(Cred, Site = rownames(Cred))
# plot

# site color assignments
col <- c("ARIK" = "pink1", "BIGC" =  "slateblue1", "BLDE" = "purple3",
         "BLUE" = "turquoise2", "CARI" =  "steelblue", "COMO" = "navyblue",
         "CUPE" = "orange", "GUIL" = "tomato", "HOPB" = "palevioletred",
         "KING" = "violetred", "LECO" = "mediumpurple1",  "LEWI" = "coral2", 
         "MART" = "yellowgreen", "MAYF" = "red2", "MCDI" = "springgreen2",
         "MCRA" = "violet", "OSKR" = "skyblue", "POSE" = "grey70", 
         "PRIN" = "grey50", "REDB" = "palegreen4", "SYCA" = "tan", 
         "TECR" = "grey30", "WALK" = "blue2", "WLOU" = "brown")
cbPalette <- c("1" = "#999999", "2" = "#F0E442", "3" = "#E69F00", "5" = "#56B4E9", "6" = "#009E73")
# plot with sites
gg = ggplot(data = data.scores, aes(x = MDS1, y = MDS2)) + 
  geom_point(data = data.scores, aes(colour = WQGroup), size = 3) + 
  scale_colour_manual(values = cbPalette) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "Water\nQuality\nGroup")

gg
# plot with sites and environmental variables
gg = ggplot(data = data.scores, aes(x = MDS1, y = MDS2)) + 
  geom_point(data = data.scores, aes(colour = Site, shape= Type), size = 3) + 
  scale_colour_manual(values = col)  + 
  geom_segment(aes(x = 0, y = 0, xend = MDS1, yend = MDS2), 
               data = Cred, linewidth =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = Cred, aes(x = MDS1, y = MDS2), colour = "grey30", 
            fontface = "bold", label = row.names(Cred)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Site")

gg
