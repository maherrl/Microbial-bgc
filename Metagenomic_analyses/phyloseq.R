library("phyloseq")
library("data.table")
library("ggplot2")
library("tidyverse")
library("vegan")
library("microbiome")
library("wrMisc")
library(RColorBrewer)


gen.tab <- read.csv("./results/gene_tab_KOfam.csv", header = TRUE)
gen.tab <- read.csv("./results/gene_tab_KEGG_Class.csv", header = TRUE)

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
physeq <- subset_samples(physeq, Site != "OKSR")

#saveRDS(physeq, file = "./data/Physeq_KOfam.rds")
readRDS(file = "./data/Physeq_KOfam.rds")
#saveRDS(physeq, file = "./data/Physeq_KEGG_Module.rds")
#physeq <- Genes_phyloseq_obj_trans

# SUMMARIZE SEQUENCING DEPTHS
sdt = data.table(as(sample_data(physeq), "data.frame"),
                 TotalReads = sample_sums(physeq), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
pSeqDepth + facet_wrap(~Type)

ggplot(sdt, aes(x = Site, y = TotalReads)) + 
  geom_boxplot() +
  geom_point() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

# TAXA TOTAL COUNTS HISTOGRAM
tdt = data.table(tax_table(physeq),
                 TotalCounts = taxa_sums(physeq),
                 OTU = taxa_names(physeq))
ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() + 
  ggtitle("Histogram of Total Counts")

# TAXA PREVALENCE HISTOGRAM
mdt = psmelt(physeq)
prevdt = mdt[, list(Prevalence = sum(Abundance > 0), 
                    TotalCounts = sum(Abundance)),
             by = "OTU"] # getting an error here
ggplot(prevdt, aes(Prevalence)) + 
  geom_histogram() + 
  ggtitle("Histogram of Taxa Prevalence")

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
ggplot(data = stats_ben, aes(x = Site, y = Richness, color = Site)) + 
  geom_boxplot() +
  geom_point() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_color_manual(values = col) +
  ylab("Inverse simpson")
  ggtitle("KOfam functions")
  

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
plot_ordination(physeq_ben, neon.ord, color = "Site", shape = "Type", 
                title = "NMDS of KEGG Modules for benthic samples (stress = 0.10)") +
  scale_color_manual(values = col) +
  geom_point(size = 3) + 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black"))

# subset to epilithon only
physeq_ben_epil <- subset_samples(physeq_ben, Type == "EPILITHON")
set.seed(1993)
neon.ord <- ordinate(physeq_ben_epil, "NMDS", "bray", trymax = 50)
neon.ord$stress
plot_ordination(physeq_ben_epil, neon.ord, color = "Group", 
                title = "NMDS of KEGG Modules benthic epilithon samples (stress = 0.16)") +
  scale_color_manual(values = cbPalette) +
  geom_point(size = 3) +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black"))

# PCoA
set.seed(1994)
neon.ord <- ordinate(physeq_ben, "PCoA", "bray", trymax = 50)
plot_ordination(physeq_ben, neon.ord, color = "Site", shape = "Type", 
                title = "PCoA of KEGG Modules for benthic samples") +
  scale_color_manual(values = col) +
  geom_point(size = 3) + 
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black"))

# subset to epilithon only
physeq_ben_epil <- subset_samples(physeq_ben, Type == "EPILITHON")
set.seed(1993)
neon.ord <- ordinate(physeq_ben_epil, "PCoA", "bray", trymax = 50)
plot_ordination(physeq_ben_epil, neon.ord, color = "Site", 
                title = "PCoA of KEGG Modules benthic epilithon samples") +
  scale_color_manual(values = col) +
  geom_point(size = 3) +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black"))

#  KEGG Module stacked barplot
kegg <- physeq %>%
  subset_samples(Source == "benthic") %>%
  subset_samples(Site != "OKSR") %>%
  tax_glom(taxrank = "Class") # change to Order for high resolution

psmelt(kegg) %>%
  aggregate(Abundance ~ Site + Type + Class, FUN = mean) -> kegg.agr

psmelt(kegg) -> kegg.agr2
  aggregate(Abundance ~  + Type + Class, FUN = mean) -> kegg.agr2

Palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                      "navyblue", "purple3", "turquoise2","#000000")

n <- 41
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col=sample(col_vector, n)

ggplot(kegg.agr, aes(x=Site, y = Abundance, fill= Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual("KEGG Modules", values = Palette) +
  facet_grid(~Type, scale="free", space = "free") +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Average Reads Per Kilobase Million") +
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black", angle = 90, vjust = 0.5, hjust = 1))

ggplot(kegg.agr, aes(x=, y = Abundance, fill= Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual("KEGG Modules", values = Palette) +
  facet_grid(~Type, scale="free", space = "free") +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Average Reads Per Kilobase Million") +
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black", angle = 90, vjust = 0.5, hjust = 1))

# permanova
sampledf <- data.frame(sample_data(physeq_ben))
dist <- distance(physeq_ben, method = "bray")
adonis2(dist ~ Site*Type, sampledf)

## SPECIFIC KEGG MODULE ABUNDANCE PLOTS
## Nitrogen plot
nitr <- physeq %>%
  subset_samples(Source == "benthic") %>%
  subset_samples(Site != "OKSR") %>%
  #transform("compositional") %>%
  subset_taxa(accession == "M00175" | accession == "M00529" | accession == "M00530" | 
                accession == "M00615" | accession == "M00531" | accession == "M00528")
nitr
plot_bar(nitr, "Site", "Abundance", "pathway")
cbPalette <- c("M00175" = "#999999", "M00529" = "#E69F00", "M00530" = "#56B4E9", 
               "M00615" = "#009E73", "M00531" = "#F0E442", "M00531" = "#0072B2",
               "M00528" = "#D55E00")
legend_title <- "Nitrogen cycling pathway"
plot_composition(nitr, average_by = "Site") +
                scale_fill_manual(legend_title, values = cbPalette, labels=c("Nitrogen fixation, nitrogen => ammonia",
                                                         "Nitrification, ammonia => nitrite",
                                                         "Denitrification, nitrate => nitrogen",
                                                         "Dissimilatory nitrate reduction, nitrate => ammonia",
                                                         "Nitrate assimilation",
                                                         "Assimilatory nitrate reduction, nitrate => ammonia")) +
  ylab("Pathway Gene Coverage (RPKM)") +
  xlab("Site") +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black",angle=90))
  
# Sulfur plot
sulf <- physeq %>%
  subset_samples(Source == "benthic") %>%
  subset_samples(Site != "OKSR") %>%
  subset_taxa(accession == "M00176" | accession == "M00595" | accession == "M00596" | 
                accession == "M00616")
cbPalette <- c("M00176" = "#999999", "M00595" = "#E69F00", "M00596" = "#56B4E9", 
               "M00616" = "#009E73")
legend_title <- "Sulfur cycling pathway"
plot_composition(sulf, average_by = "Site") +
    scale_fill_manual(legend_title, values = cbPalette, labels=c("Assimilatory sulfate reduction, sulfate => H2S",
                                           "Thiosulfate oxidation by SOX complex, thiosulfate => sulfate",
                                           "Dissimilatory sulfate reduction, sulfate => H2S",
                                           "Sulfate-sulfur assimilation")) +
  ylab("Pathway Gene Coverage (RPKM)") +
  xlab("Site") +
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black",angle=90))
# Carbon plot
carb <- physeq %>%
  subset_samples(Source == "benthic") %>%
  subset_samples(Site != "OKSR") %>%
  #transform("compositional") %>%
  subset_taxa(accession == "M00174" | accession == "M00567" | accession == "M00001" | 
                accession == "M00009" | accession == "M00307" )
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
legend_title <- "Carbon cycling pathway"
plot_composition(carb,average_by = "Site") +
  scale_fill_manual(legend_title, labels=c("Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate",
                                           "Citrate cycle (TCA cycle, Krebs cycle)",
                                           "Methane oxidation, methanotroph, methane => formaldehyde",
                                           "Pyruvate oxidation, pyruvate => acetyl-CoA",
                                           "Methanogenesis, CO2 => methane"),
                    values = cbPalette) +
  ylab("Pathway Gene Coverage (RPKM)") +
  xlab("Site" )+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black",angle=90))
