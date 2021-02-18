#############################################################################################
#                   BOXPLOTS COMPARING DISTANCES AMONG INDIVIDUALS 
#                                 OCTOBER 31, 2020
#############################################################################################
##### 

# Description of script: This script makes boxplots to compare pairwise
# distances between six categories. The comparisons are made for quantitative
# Jaccard and weighed UniFrac dissimilarities.

## Between species:
# 1) PV Allo vs PD Allo (i.e. pairwise distances between all allopatric P. vindex and
# all allopatric P. difformis)
# 2) PV Symp vs PD Symp (i.e. pairwise distances of all sympatric
# P. vindex with all sympatric P. difformis)

## Within species:
# 3) Among all allopatric PV individuals 
# 4) Among all allopatric PD individuals 
# 5) Among all sympatric PV individuals 
# 6) Among all sympatric PD individuals 

library(vegan)
library(phyloseq)

load("ps_rare_diss_shann_.RData") # this has all the distances, phyloseq objects, etc. 

# Pull out sample metdata to use for indexing
sample_data(gut.r.phyloseq) #gut.r.phyloseq is in file loaded above.
meta_dat <- sample_data(gut.r.phyloseq)

# Make matrixes out of dissimilarities
# Jaccard
gut.jacc.dist 
gut.jacc.mat <- as.matrix(gut.jacc.dist)
diag(gut.jacc.mat) <- NA
gut.jacc.mat[lower.tri(gut.jacc.mat)] <- NA 

# Weighed UniFrac
guts.wUF.dist
gut.wUF.mat <- as.matrix(guts.wUF.dist)
diag(gut.wUF.mat) <- NA
gut.wUF.mat[lower.tri(gut.wUF.mat)] <- NA 

# Make data frame for indexing matrix
indiv <- colnames(gut.jacc.mat)
patry_spec <- meta_dat$Range.Overlap.Spec
index.df <- data.frame(indiv, patry_spec)

#Get row and column numbers of allopatric or sympatric samples
allo_PV.num <- which(index.df$patry == "Allo_vindex")
allo_PD.num <- which(index.df$patry == "Allo_difformis")
symp_PV.num <- which(index.df$patry == "Symp_vind")
symp_PD.num <- which(index.df$patry == "Symp_difformis")

######## GETTING DISTANCES "BETWEEN SPECIES" ########
### 1) PV Allo vs PD Allo 
# Jaccard
bw_allo_jacc <- gut.jacc.mat[allo_PV.num, allo_PD.num]

# Weighed UniFrac
bw_allo_wUF <- gut.wUF.mat[allo_PV.num, allo_PD.num]

### 2) PV Symp vs PD Symp 
# Jaccard
bw_symp_jacc <- gut.jacc.mat[symp_PV.num, symp_PD.num]

# Weighed UniFrac
bw_symp_wUF <- gut.wUF.mat[symp_PV.num, symp_PD.num]

######## GETTING DISTANCES "WITHIN SPECIES" ########
### 3) Among all allopatric PV individuals 
#Jaccard
wiPV_allo_jacc <- gut.jacc.mat[allo_PV.num, allo_PV.num]

#Weighted UniFrac
wiPV_allo_wUF <- gut.wUF.mat[allo_PV.num, allo_PV.num]

### 4) Among all allopatric PD individuals 
#Jaccard
wiPD_allo_jacc <- gut.jacc.mat[allo_PD.num, allo_PD.num]

#Weighted UniFrac
wiPD_allo_wUF <- gut.wUF.mat[allo_PD.num, allo_PD.num]

### 5) Among all sympatric PV individuals 
#Jaccard
wiPV_symp_jacc <- gut.jacc.mat[symp_PV.num, symp_PV.num]

#Weighted UniFrac
wiPV_symp_wUF <- gut.wUF.mat[symp_PV.num, symp_PV.num]

### 6) Among all sympatric PD individuals 
#Jaccard
wiPD_symp_jacc <- gut.jacc.mat[symp_PD.num, symp_PD.num]

#Weighted UniFrac
wiPD_symp_wUF <- gut.wUF.mat[symp_PD.num, symp_PD.num]

######## MAKE BOXPLOTS ######## 
bold_a <- expression(bold("(a)"))
bold_b <- expression(bold("(b)"))
quartz()
par(mfrow=c(2,1))
jacc_box <- boxplot(list(bw_allo_jacc, bw_symp_jacc, wiPV_allo_jacc,
                         wiPD_allo_jacc, wiPV_symp_jacc, wiPD_symp_jacc),
                    ylab = "Quantitative Jaccard Dissimilarity",
                    names = c("PV v. PD allo", "PV v. PD sym", "PV allo",
                              "PD allo", "PV sym", "PD sym"), cex.axis = 0.8,
                    cex.lab = 1, col = c("orange", "purple", "orange", "purple", "orange", "purple"),
                    ylim=c(0.0, 1.0))
mtext(text=bold_a, side=3, adj = -0.065, line = 2)
wUF_box <- boxplot(list(bw_allo_wUF, bw_symp_wUF, wiPV_allo_wUF, wiPD_allo_wUF,
                        wiPV_symp_wUF, wiPD_symp_wUF), 
                   ylab = "Weighted UniFrac Dissimilarity",
                   names = c("PV v. PD allo", "PV v. PD sym", "PV allo",
                             "PD allo", "PV sym", "PD sym"), cex.axis = 0.8,
                   cex.lab = 1, col = c("orange", "purple", "orange", "purple", "orange", "purple"),
                   ylim=c(0.0, 1.0))
mtext(text=bold_b, side=3, adj = -0.065, line = 2)

