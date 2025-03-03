rm(list = ls())

library(vegan)
library(phyloseq)
library(tidyverse)
library("wrMisc")
library("cowplot")

####
#### FUNCTIONS
# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
####



physeq <- readRDS("./data/Physeq_KOfam.rds")
physeq

physeq.veg <- psotu2veg(physeq)
sampledf.veg <- pssd2veg(physeq)
# calculate alpha diversity metrics
H <- diversity(physeq.veg, "shannon")
simp <- diversity(physeq.veg, "simpson")
invsimp <- diversity(physeq.veg, "inv")
S <- specnumber(physeq.veg)
alpha <- as.data.frame(t(mergeVectors(Shannon = H, Simpson = simp, Invsimpson = invsimp, Richness = S)))
stats <- merge(alpha, sampledf.veg, by = "row.names")
write.csv(stats, file = "results/KOfam_alpha.csv")
stats_ben <- stats[stats$Source == "benthic",] # subset to benthic
stats_ben <- stats_ben[stats_ben$Type != "EPIXYLON",]



## ANALYSES
# all data
kruskal.test(Richness ~ Site, data = stats_ben)
kruskal.test(Richness ~ Type, data = stats_ben)
kruskal.test(Richness ~ Group, data = stats_ben)

kruskal.test(Shannon ~ Site, data = stats_ben)
kruskal.test(Shannon ~ Type, data = stats_ben)
kruskal.test(Shannon ~ Group, data = stats_ben)

kruskal.test(Simpson ~ Site, data = stats_ben)
kruskal.test(Simpson ~ Type, data = stats_ben)
kruskal.test(Simpson ~ Group, data = stats_ben)
# epilithon
kruskal.test(Richness ~ Site, data = stats_ben_epil)
pairwise.wilcox.test(stats_ben_epil$Richness, stats_ben_epil$Site, p.adjust.method = 'fdr')
kruskal.test(Richness ~ Group, data = stats_ben_epil)
pairwise.wilcox.test(stats_ben_epil$Richness, stats_ben_epil$Group, p.adjust.method = 'fdr')

kruskal.test(Shannon ~ Site, data = stats_ben_epil)
kruskal.test(Shannon ~ Group, data = stats_ben_epil)

kruskal.test(Simpson ~ Site, data = stats_ben_epil)
kruskal.test(Simpson ~ Group, data = stats_ben_epil)

cbPalette <- c("#009E73", "#E69F00", "#56B4E9")
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
type <- c("#000000","#0072B2", "#CC79A7")
col <- c("ARIK" = "pink1", "BIGC" =  "slateblue1", "BLDE" = "purple3",
         "BLUE" = "turquoise2", "CARI" =  "steelblue", "COMO" = "navyblue",
         "CUPE" = "orange", "GUIL" = "tomato", "HOPB" = "palevioletred",
         "KING" = "violetred", "MAYF" = "red2", "MCDI" = "springgreen2",
         "REDB" = "palegreen4", "SYCA" = "tan", "WLOU" = "brown")

names <- c(`EPILITHON` = "L", `EPIPSAMMON` = "P")

a <- ggplot(data = stats_ben, aes(x = Site, y = Richness, color = Site)) + 
  geom_boxplot() +
  geom_jitter(width = 0.25) + 
  facet_grid(~Type, scale="free", space = "free", labeller = as_labeller(names)) +
  scale_color_manual(values = col) +
  ylab("Observed\nKOfam Functions") +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black", angle = 90, vjust = 0.5, hjust = 1)) +
  theme(strip.text.x=element_text(size=14,color="black"))+
  theme(strip.text.y=element_text(size=14,color="black"))

b <- ggplot(data = stats_ben, aes(x = Group, y = Richness, color = Group)) + 
  geom_boxplot() +
  geom_jitter(width = 0.25) + 
  facet_grid(~Type, scale="free", space = "free", labeller = as_labeller(names)) +
  scale_color_manual(values = colorBlindGrey8) +
  ylab("Observed\nKOfam Functions") +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black")) +
  theme(strip.text.x=element_text(size=14,color="black"))+
  theme(strip.text.y=element_text(size=14,color="black"))

plot_grid(a,b, ncol = 1, labels = c("A","B"))

ggplot(data = stats_ben, aes(x = Type, y = Richness, color = Type)) + 
  geom_boxplot() +
  geom_jitter(width = 0.25) + 
  scale_color_manual(values = type) +
  ylab("Observed\nKOfam Functions") +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black", angle = 35, vjust = 1, hjust = 1)) +
  theme(legend.position = "none")
