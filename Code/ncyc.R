## NCyc Analyses
rm(list = ls())

library(microbiome)
library(vegan)
library("tidyverse")
library("Hmisc")
library("corrplot")
library("lme4")
library("nortest")
library("car")
library("glmm")
library("lmerTest")
library("lsmeans")

ncyc <- read.csv(file = "results/gene_tab_NCyc.csv", header = TRUE)

gen.tax <- as.matrix(ncyc[,c(2:4,1)])
rownames(gen.tax) <- ncyc$accession
colnames(gen.tax) <- c("Kingdom", "Phylum", "Class", "Order")
TAX <- tax_table(gen.tax)

# make otu (gene) table
gen.otu <- as.matrix(ncyc[,5:ncol(ncyc)])
rownames(gen.otu) <- ncyc$accession
gen.otu[is.na(gen.otu)] <- 0
OTU <- otu_table(gen.otu, taxa_are_rows = TRUE)

gen.meta <- read.csv("./data/metag_map.csv", header = TRUE, row.names = 1)
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


## PLOTTING
## SPECIFIC KEGG MODULE ABUNDANCE PLOTS
## Nitrogen plot
nitr <- physeq %>%
  subset_samples(Source == "benthic") %>%
  subset_samples(Site != "OKSR") %>%
  subset_taxa(Kingdom != "Organic degradation and synthesis") %>%
  subset_taxa(Kingdom != "Others") %>%
  tax_glom(taxrank = "Kingdom")

psmelt(nitr) %>%
  aggregate(Abundance ~ Site_Type + Kingdom, FUN = mean) -> nitr.agr
nitr.agr$Kingdom <- as.factor(nitr.agr$Kingdom)

cbPalette <- c("1" = "#999999", "2" = "#F0E442", "3" = "#E69F00", "5" = "#56B4E9", "6" = "#009E73")
col2 <- c("Anammox" = "#555555", "Assimilatory nitrate reduction" = "#364B9A", 
         "Denitrification" = "#4A7BB7", 
         "Denitrification/Dissimilatory nitrate reduction" = "#6EA6CD", 
         "Dissimilatory nitrate reduction" = "#98CAE1",
         "Nitrification" = "#FDB366", "Nitrogen fixation"="#EAECCC")
col <- c("ARIK" = "pink1", "BIGC" =  "slateblue1", "BLDE" = "purple3",
         "BLUE" = "turquoise2", "CARI" =  "steelblue", "COMO" = "navyblue",
         "CUPE" = "orange", "GUIL" = "tomato", "HOPB" = "palevioletred",
         "KING" = "violetred", "LECO" = "mediumpurple1",  "LEWI" = "coral2", 
         "MART" = "yellowgreen", "MAYF" = "red2", "MCDI" = "springgreen2",
         "MCRA" = "violet", "OSKR" = "skyblue", "POSE" = "grey70", 
         "PRIN" = "grey50", "REDB" = "palegreen4", "SYCA" = "tan", 
         "TECR" = "grey30", "WALK" = "blue2", "WLOU" = "brown")

ggplot(nitr.agr, aes(x=Site_Type, y = Abundance, fill= Kingdom)) +
  geom_bar(stat = "identity") +
  scale_fill_manual("Nitrogen Cycling Pathways", values = col) +
  #facet_grid(~Type, scale="free", space = "free") +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Average RPKM") +
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=10,color="#555555", angle = 90, vjust = 0.5, hjust = 1))
# Makes the same graph values! yay!
nitr %>% plot_composition(average_by = "Site_Type") + scale_fill_brewer("Kingdom", palette = "Paired")

######################################################################
## ORDINATION
nitr <- physeq %>%
  subset_samples(Source == "benthic") %>%
  subset_samples(Site != "OKSR") %>%
  subset_samples(Type != "EPIXYLON") %>%
  subset_taxa(Kingdom != "Organic degradation and synthesis") %>%
  subset_taxa(Kingdom != "Others")

# permanova
sampledf <- data.frame(sample_data(nitr))
dist <- distance(nitr, method = "bray")
adonis2(dist ~ Site*Type, sampledf)
# permdisp
anova(betadisper(dist, sampledf$Type, bias.adjust = TRUE))
anova(betadisper(dist, sampledf$Group, bias.adjust = TRUE))

env <- read.csv(file = "./data/EnvironmentalData-v5-Benthic.csv", header = TRUE, row.names = 1)
rownames(env) <- str_replace(rownames(env), "-DNA[1-9]", "")
rownames(env) <-  str_replace(rownames(env), ".DNA", "")
env$ef_grp <- paste(env$Site,env$Date, sep = "_")
sample_data(nitr) <- env
nitr

physeq.cca <- ordinate(nitr, "CCA")

plot_ordination(nitr, physeq.cca,
                color = "Site", shape = "Type") +
  geom_point(size = 3) +
  scale_color_manual(values = col) +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black"))
  
 
# envfit
merged.nitr <- merge_samples(nitr, "ef_grp", fun = sum)
merged.cca <- ordinate(merged.nitr, "CCA")


env_merged <- pssd2veg(merged.nitr)
env_merged <- env_merged[,8:ncol(env_merged)]
env_merged <- env_merged[,-ncol(env_merged)]
env_merged$site_date <- rownames(env_merged)

en = envfit(merged.cca, env_merged, permutations = 999, na.rm = TRUE)
en
ordiplot(merged.cca, display = "sites")
plot(en)
# only significant pvalues
A <- as.list(en$vectors)
pvals <- as.data.frame(A$pvals)
arrows<-as.data.frame(A$arrows*sqrt(A$r))
C<-cbind(arrows, pvals)
#subset
Cred<-subset(C,pvals<0.05)
Cred <- cbind(Cred, Site = rownames(Cred))

data.scores <- as.data.frame(merged.cca$CA$u)
data.scores$Site <- rownames(data.scores)
data.scores[c('Site','Date')] <- str_split_fixed(data.scores$Site, '_', 2)

head(data.scores)

gg = plot_ordination(merged.nitr, merged.cca) + 
  geom_point(data = data.scores, aes(colour = Site),size = 3) + 
  scale_colour_manual(values = col)  + 
  geom_segment(aes(x = 0, y = 0, xend = CA1, yend = CA2), 
               data = Cred, linewidth =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = Cred, aes(x = CA1, y = CA2), colour = "grey30", 
            fontface = "bold", label = row.names(Cred)) + 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=14,color="black"))+
  theme(axis.title.y=element_text(size=14,color="black"))+
  theme(axis.text.y=element_text(size=12,color="black"))+
  theme(axis.text.x=element_text(size=12,color="black")) +
  labs(colour = "Site")

gg




## Correlations

nitr <- physeq %>%
  subset_samples(Source == "benthic") %>%
  subset_samples(Site != "OKSR") %>%
  subset_samples(Type == "EPILITHON") %>%
  psmelt()

# Now I want to duplicate the dissimilatory nitrate reduction and denitrification genes
dd.genes <- nitr[nitr$Kingdom == "Denitrification/Dissimilatory nitrate reduction",]
dd.genes$Kingdom <- "Denitrification"
de.genes <- dd.genes
dd.genes$Kingdom <- "Dissimilatory nitrate reduction"
nitr.no.dd <- nitr[nitr$Kingdom != "Denitrification/Dissimilatory nitrate reduction",]
nitr.dup <- rbind(nitr.no.dd, dd.genes, de.genes)
unique(nitr.dup$Kingdom)
# remove unwanted pathways
nitr.dup <- subset(nitr.dup, nitr.dup$Kingdom != "Organic degradation and synthesis")
nitr.dup <- subset(nitr.dup, nitr.dup$Kingdom != "Others")
nitr.dup <- subset(nitr.dup, nitr.dup$Kingdom != "Anammox")
# aggregate by sample and pathway
nitr.dup %>%
  aggregate(Abundance ~ Sample + Kingdom, FUN = mean) -> nitr.agr2

# environmental data
env <- read.csv(file = "./data/EnvironmentalData-v5-Benthic.csv", header = TRUE)
env$dnaSampleID <- str_replace(env$dnaSampleID, "-DNA[1-9]", "")
env$dnaSampleID <-  str_replace(env$dnaSampleID, ".DNA", "")
env_sub <- env[,-c(44:46)]
# remove nas from both dataframes
env_nona <- na.omit(env_sub)
env_nitr <- env_nona[env_nona$dnaSampleID %in% nitr.agr2$Sample,]
nitr.agr3 <- nitr.agr2[nitr.agr2$Sample %in% env_nitr$dnaSampleID,]

nitr.agr3 %>%
  pivot_wider(names_from = Kingdom, values_from = Abundance) -> nitr.wide

nitr.wide <- nitr.wide[,-1]
nitr.wide <- as.matrix(nitr.wide)
env_nitr <- env_nitr[,-c(1:8)]
env.mat <- as.matrix(env_nitr)
env.standard <- decostand(env.mat, method = "standardize", MARGIN = 2)

out <- cor(nitr.wide, env.standard, method = "spearman")

corrplot(out, type = "full",
         tl.col = "black", tl.srt = 45)

# to get the p values:
# doesn't seem to work for EPIPSAMMON only samples...
p.out <- t(sapply(1:nrow(t(nitr.wide)), function(x) {
  sapply(1:nrow(t(env.standard)), function(y) {
    rcorr(nitr.wide[x,],env.standard[y,], type = "spearman")[[3]][1,2]
  })
}))

## statistics
nitr.stat <- physeq %>%
  subset_samples(Source == "benthic") %>%
  subset_samples(Site != "OKSR") %>%
  subset_samples(Type != "EPIXYLON") %>%
  subset_taxa(Kingdom != "Organic degradation and synthesis") %>%
  subset_taxa(Kingdom != "Others") %>%
  psmelt()

# test for normality
hist(log(nitr.stat$Abundance))
qqnorm(nitr.stat$Abundance)
abline(0,1)
ModelData = nitr.stat$Abundance
shapiro.test(ModelData[0:5000])
ad.test(ModelData)
leveneTest(nitr.stat$Abundance ~ nitr.stat$Type)

# epilithon vs epipsammon
model1 <- lmer(rank(Abundance) ~ Type + (1|Site), data = nitr.stat)
anova(model1)
summary(model1)
model2 <- lmer(rank(Abundance) ~ Group + (1|Site), data = nitr.stat)
anova(model2)
summary(model2)
model3 <- lmer(rank(Abundance) ~ Type * Kingdom + (1|Site), data = nitr.stat)
anova(model3)
summary(model3)
lsmeans(model3, pairwise ~ Type * Kingdom)
nitr.stat.1 <- nitr.stat[nitr.stat$Kingdom == "Anammox",]
anova(lmer(rank(Abundance) ~ Type + (1|Site), data = nitr.stat.1))
nitr.stat.2 <- nitr.stat[nitr.stat$Kingdom == "Assimilatory nitrate reduction",]
anova(lmer(rank(Abundance) ~ Type + (1|Site), data = nitr.stat.2))
nitr.stat.3 <- nitr.stat[nitr.stat$Kingdom == "Denitrification",]
anova(lmer(rank(Abundance) ~ Type + (1|Site), data = nitr.stat.3))
nitr.stat.4 <- nitr.stat[nitr.stat$Kingdom == "Denitrification/Dissimilatory nitrate reduction",]
anova(lmer(rank(Abundance) ~ Type + (1|Site), data = nitr.stat.4))
nitr.stat.5 <- nitr.stat[nitr.stat$Kingdom == "Dissimilatory nitrate reduction",]
anova(lmer(rank(Abundance) ~ Type + (1|Site), data = nitr.stat.5))
nitr.stat.6 <- nitr.stat[nitr.stat$Kingdom == "Nitrification",]
anova(lmer(rank(Abundance) ~ Type + (1|Site), data = nitr.stat.6))
nitr.stat.7 <- nitr.stat[nitr.stat$Kingdom == "Nitrogen fixation",]
anova(lmer(rank(Abundance) ~ Type + (1|Site), data = nitr.stat.7))
