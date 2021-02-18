#################################################################################
#                 Making Relative Abundance Plots   
#                      Phan Biogeo 2019 Project
#                        June 19th, 2020
################################################################################

library(ggplot2)
library(phyloseq)
library(tidyverse)

load(file = "ps_rare_diss_shann_.RData") # files in this: psotu2veg, pssd2veg, gut.r.3507.1k, gut.r.data, gut.jacc.dist, soil.r.26336.1k, soil_shann, 
#gut.r.phyloseq, guts.wUF.dist, guts.uwUF.dist, PV.r.data, PD.r.data, PV.jacc.dist, PD.jacc.dist, PV.r.phyloseq, PD.r.phyloseq, PV.wUF.dist, 
# PV.uwUF.dist, PD.wUF.dist, PD.uwUF.dist

############# PHYLUM LEVEL ##############

# TURN ASVs INTO PHYLUM LEVEL
gut.phylum.glom <-  tax_glom(gut.r.phyloseq, taxrank = "Phylum") 
dim(tax_table(gut.phylum.glom)) # good, this is only phyla

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
relabun.phyla.0 <- transform_sample_counts(gut.phylum.glom, function(x) x / sum(x) )
rownames(otu_table(relabun.phyla.0))

# MERGE SAMPLES so that we only have two combined abundances for PV and PD
relabun.phyla.1 <- merge_samples(relabun.phyla.0, group = "Species")
sample_data(relabun.phyla.1)

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabun.phyla.2 <- transform_sample_counts(relabun.phyla.1, function(x) x / sum(x))
sample_data(relabun.phyla.2)

# NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE 

relabun.phyla.df <-psmelt(relabun.phyla.2)
dim(relabun.phyla.df) #54 rows because 27 phyla per beetle. 
sum(relabun.phyla.df[,3]) 
colnames(relabun.phyla.df) 

relabun.phyla.df$Phylum[relabun.phyla.df$Abundance < 0.01] <- "< 1% abund."

#see names of phyla
uniphyl <- unique(relabun.phyla.df$Phylum)
uniphyl


# PLOT
quartz()
gut.phy.plot <- ggplot(data=relabun.phyla.df, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 18, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 18, face = "bold"))
gut.phy.plot + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("#4575b4", "#d73027", "#fc8d59", "#fee090", "#91bfdb", "grey"), 
                    name= "Phylum", breaks= c("D_1__Firmicutes", "D_1__Proteobacteria", "D_1__Bacteroidetes", "D_1__Actinobacteria", "D_1__Fusobacteria", "< 1% abund."), 
                    labels =c("Firmicutes", "Proteobacteria", "Bacteroidetes", "Actinobacteria", "Fusobacteria", "< 1% abund.")) + theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=2)) + theme(legend.text = element_text(colour="black", size = 16))  + theme(legend.title = element_blank())


# SEE ABUNDANCES OF EACH PHYLA
# based on discussion found here: https://github.com/joey711/phyloseq/issues/1089#issuecomment-471334036

library(dplyr)
library(forcats)

colnames(relabun.phyla.df)
relabun.phyla.df[,2] #Species (i.e. PV versus PD) was renames automatically to "sample" to distinguish between Phanaeus spp and bacterial species

top_phyla <- relabun.phyla.df %>%
  group_by(Sample, Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) %>%
  View()

########### REPEAT FOR FAMILY LEVEL ############

# TURN ASVs INTO FAM LEVEL
gut.fam.glom <-  tax_glom(gut.r.phyloseq, taxrank = "Family") #179 taxa and 199 samples
dim(tax_table(gut.fam.glom)) 
tax_table(gut.fam.glom)

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
relabun.fam.0 <- transform_sample_counts(gut.fam.glom, function(x) x / sum(x) )
rownames(otu_table(relabun.fam.0)) #rownames are 179 taxa

# MERGE SAMPLES so that we only have two combined abundances for PV and PD
relabun.fam.1 <- merge_samples(relabun.fam.0, group = "Species")
sample_data(relabun.fam.1)

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabun.fam.2 <- transform_sample_counts(relabun.fam.1, function(x) x / sum(x))
sample_data(relabun.fam.2)

# NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE: 

relabun.fam.df <-psmelt(relabun.fam.2)
dim(relabun.fam.df) #358
sum(relabun.fam.df[,3])

relabun.fam.df$Family[relabun.fam.df$Abundance < 0.01] <- "< 1% abund." #only 23 families if we do it this way

#see names of families
unifam <- unique(relabun.fam.df$Family)
unifam
length(unifam)

# PLOT!
quartz()
gut.fam.plot <- ggplot(data=relabun.fam.df, aes(x=Sample, y=Abundance, fill=Family)) + theme(axis.title.y = element_text(size = 20, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 18, face = "bold"))
gut.fam.plot + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("#a50026", "#fee090",  "#74add1", "#313695","#5e4fa2", "#f46d43", "#fdae61", "#4575b4", "#FFFF00", "#EE82EE","#8B008B", "#9ACD32", "#CD5C5C", "grey"), 
                    name= "Family", breaks= c("D_4__Enterococcaceae", "D_4__Moraxellaceae", "D_4__Porphyromonadaceae", "D_4__Enterobacteriaceae", "D_4__Planococcaceae", "D_4__Flavobacteriaceae", "D_4__Micrococcaceae", "D_4__Streptococcaceae", "D_4__Comamonadaceae", "D_4__Nocardiaceae", "D_4__Fusobacteriaceae", "D_4__Family XI", "D_4__Rickettsiales Incertae Sedis", "< 1% abund."), 
                    labels =c("Enterococcaceae", "Moraxellaceae", "Porphyromonadaceae", "Enterobacteriaceae", "Planococcaceae", "Flavobacteriaceae", "Micrococcaceae", "Streptococcaceae", "Comamonadaceae", "Nocardiaceae", "Fusobacteriaceae", "Family XI", "Rickettsiales Incertae Sedis", "< 1% abund.")) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5)) + theme(legend.text = element_text(colour="black", size = 16))  + theme(legend.title = element_blank())


# PERCENTAGES OF EACH FAMILY

library(dplyr)
library(forcats)

colnames(relabun.fam.df)
relabun.fam.df[,2] #Species (i.e. PV versus PD) was renames automatically to "sample" to distinguish between Phanaeus spp and bacterial species

top_fam <- relabun.fam.df %>%
  group_by(Sample, Family) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) %>%

# LET'S GO AHEAD AND LOOK AT ARCHAEA

# TURN ASVs INTO ARCAHEA LEVEL
gut.king.glom <-  tax_glom(gut.r.phyloseq, taxrank = "Kingdom") #2 taxa and 199 samples
dim(tax_table(gut.king.glom)) 
tax_table(gut.king.glom) #good only archaea and bacteria

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
relabun.king.0 <- transform_sample_counts(gut.king.glom, function(x) x / sum(x) )
rownames(otu_table(relabun.king.0)) #rownames are 2 taxa

# MERGE SAMPLES so that we only have two combined abundances for PV and PD
relabun.king.1 <- merge_samples(relabun.king.0, group = "Species")
sample_data(relabun.king.1)

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabun.king.2 <- transform_sample_counts(relabun.king.1, function(x) x / sum(x))
sample_data(relabun.king.2)

relabun.king.df <-psmelt(relabun.king.2)
dim(relabun.king.df) 

top_king <- relabun.king.df %>%
  group_by(Sample, Kingdom) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) %>%
  View()

save(relabun.phyla.df, uniphyl, gut.phy.plot, relabun.fam.df, unifam, gut.fam.plot, file ="barplots_June20")
