#############################################################################################
#                           INDICATOR SPECIES ANALYSIS 
#                       JUNE 18, 2020- (UPDATED JUL 20) final version
#############################################################################################
##### 

setwd("/Users/clairewinfrey/Desktop/2019_biogeo_phan_GM/_All_pops")

library(indicspecies)
library(vegan)
library(phyloseq) # necessary to load all of the objects in the below RData file

load("ps_rare_diss_shann_.RData")

# load(file = "indicator_sp_analysis_June18") # this was made in this file
# load(file="cattle.indi.sp.Sept7") #also made in this file!

# USEFUL TO DO INDICATOR SPECIES ANALYSIS TO FIND:
# 1. INDICATOR SPECIES FOR 1) P. VINDEX IN SYMP, 2) P. VINDEX IN ALLO, 3) P. DIFFORMIS IN SYMP, AND 4) P. DIFFORMIS IN ALLO
# 2. FOR BEETLES WITH DIFFERENT LEVELS OF CATTLE PRESENT IN ENVIRONMENT

######## ALL PV VERSUS ALL PD

# 1. NEED MATRIX WITH SITES IN ROWS AND SPECIES AS COLUMNS 

gut.r.data
dim(gut.r.data) # 199 x 1358 so site (i.e. beetles) are rows and OTUs are columns. 

########## MAKE 4 GROUPS- I.E. EACH PHAN SPECIES IN SYMP AND ALLO

PV_PD_symp.allo <- c(rep("PV_allo", 6), rep("PV_symp", 8), rep("PD_symp", 12), rep("PV_symp", 6), rep("PD_symp", 20), rep("PV_symp", 7), rep("PD_allo", 7),
                     rep("PD_symp", 12), rep("PV_symp", 9), rep("PD_symp", 11), rep("PV_symp", 11), rep("PD_allo", 6), rep("PD_symp", 12), rep("PV_symp", 11),
                     rep("PV_allo", 3), rep("PD_allo", 18), rep("PV_allo", 22), rep("PD_allo", 12), rep("PV_allo", 6))

set.seed(9)
PV_PD_symp.allo.rg <- multipatt(x = gut.r.data, cluster = PV_PD_symp.allo, func = "r.g", control = how(nperm = 99999))
summary(PV_PD_symp.allo.rg)
#added to table in results 

save(PV_PD_symp.allo, file = "indicator_sp_analysis_final")

##### GET RELATIVE ABUNDANCES #####
####### should be able to do this for each ASV by sample type 

sample_data(gut.r.phyloseq)


gut_merge_by_overlap <- merge_samples(gut.r.phyloseq, group = "Range.Overlap.Spec") # = need to be P.V. symp, PV allo, PD symp, PD ALLO #then do transform by sample type
# now we have samples combined by overlap - i.e. Symp. PV, Allo. PV, Symp. PD, Allo PD

sample_data(gut_merge_by_overlap) #looks good!

#Relative species abundance is calculated by dividing the number of species from one group by the total number of species from all groups.
byrange.rel.abund <- transform_sample_counts(gut_merge_by_overlap, function(x) x / sum(x) )
rownames(otu_table(byrange.rel.abund)) #YEP! Groups that we expect!
colnames(otu_table(byrange.rel.abund))
otu_table(byrange.rel.abund)

##### write this table to Excel object to look for each one of our ASVs
library(xlsx)
write.xlsx(otu_table(byrange.rel.abund), "rel.abund.indicator.analysis.xlsx")
getwd
####### JULY 20, 2020 ##### ADDING LEVEL OF CATTLE TO THIS ANALYSIS ######

## SAVE OBJECTS MADE BELOW:
save(cattle_level.all, cattle.indicspec, file = "cattle.indi.sp.Sept7")
cattle_level.all <- c(rep("high", 44), rep("med", 22), rep("high", 21), rep("low", 22), rep("high", 6), rep("low", 23), rep("high", 12), 
                      rep("med", 9), rep("low", 22), rep("med", 12), rep("high", 6))

set.seed(93)
cattle.indicspec <- multipatt(x = gut.r.data, cluster = cattle_level.all, func = "r.g", control = how(nperm = 99999))
summary(cattle.indicspec)


##### GET RELATIVE ABUNDANCES FOR CATTLE #####

# Because cattle abundance is an ordinal variable in our current gut.r.phyloseq
# object, we cannot use phyloseq's merge_samples function to combine samples
# based on cattle level. So, below, we add another matadata variable, nominal
# cattle abundance level, so that we can do so and then get the relative abundance
# of the different ASVs by level of cattle

# Make levels of cattle a nominal variable in a dataframe
catt_nominal <- data.frame(cattle_level.all) 

samps <- sample_names(gut.r.phyloseq) #extract sample names from phyloseq object

rownames(catt_nominal) <- samps #make sample names rownames
catt_nominal

# Turn into phyloseq-style sample data
catt_nominal <- sample_data(catt_nominal)
head(catt_nominal)

# Merge with original phyloseq object
catt_nom_ps <- merge_phyloseq(gut.r.phyloseq, catt_nominal)
head(sample_data(catt_nom_ps)) #now has metadata variable "cattle_level.all"

# Now, combine samples based on amount of cattle in their environment:
gut_merge_by_catt <- merge_samples(catt_nom_ps, group = "cattle_level.all")

sample_data(gut_merge_by_catt) #looks good!

#Relative species abundance is calculated by dividing the number of species from one group by the total number of species from all groups.
bycatt.rel.abund <- transform_sample_counts(gut_merge_by_catt, function(x) x / sum(x) )
rownames(otu_table(bycatt.rel.abund)) #YEP! Groups that we expect!
colnames(otu_table(bycatt.rel.abund))
otu_table(bycatt.rel.abund)

##### write this table to Excel object to look for each one of our ASVs
library(xlsx)
write.xlsx(otu_table(bycatt.rel.abund), "rel.abund.catt.indic.analysis.xlsx")
