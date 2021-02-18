#################################################################################
# PS_Rarefying_Dissim.R  
# Setting up phyloseq objects, rarefying, Shannon, and getting dissimilarities     
# Phan Biogeo 2019 Project
# Previous script was q2_to_PS2.R
################################################################################

##### NOTE: all_samp_physeq WAS RENAMED TO all_samp_physeq_raw HERE (SO DIFFERENT THAN  q2_to_phyloseq2.R)

library("phyloseq")
library("ape")
library("vegan")

#################################################################################
#                 PHYLOSEQ OBJECTS BY SAMPLE TYPE
################################################################################

library(vegan)

# CREATE SEPARATE PHYLOSEQ OBJECTS FOR GUT AND SOIL SAMPLES
# guts, controls, and blanks
only.gut.samples <-c("AN.PV.276", "AN.PV.277", "AN.PV.278", "AN.PV.279", "AN.PV.280", "AN.PV.281", "AN.PV.282", "AN.PV.283", "AN.PV.284", "BR2.PD.324", "BR2.PD.325", "BR2.PD.326", "BR2.PD.327", "BR2.PD.328", "BR2.PD.329", "BR2.PD.330", "BR2.PD.332", "BR2.PD.333", "BR2.PD.334", "BR2.PD.335", "BR2.PD.336", "B2R.PV.309", "B2R.PV.310", "B2R.PV.311", "B2R.PV.312", "B2R.PV.313", "B2R.PV.316", "B2R.PV.317", "B2R.PV.318", "B2R.PV.319", "B2R.PV.320", "B2R.PV.321", "B2R.PV.322", "BCC.PD.285", "BCC.PD.288", "BCC.PD.291", "BCC.PD.292", "BCC.PD.295", "BCC.PD.296", "BCC.PD.298", "BCC.PD.302", "BCC.PD.303", "BCC.PD.304", "BCC.PD.307", "BCC.PD.308", "BCC.PV.286", "BCC.PV.287", "BCC.PV.289", "BCC.PV.290", "BCC.PV.297", "BCC.PV.300", "BCC.PV.301", "PCRblankA", "PCRblankB", "PCRblankC", "EXT.BLANK", "CTW.PD.153", "CTW.PD.154", "CTW.PD.155", "CTW.PD.156", "CTW.PD.161", "CTW.PD.163", "CTW.PD.164", "CTW.PD.186", "CTW.PD.188", "CTW.PD.190", "CTW.PD.191", "CTW.PD.192", "CTW.PV.150", "CTW.PV.151", "CTW.PV.152", "CTW.PV.157", "CTW.PV.158", "CTW.PV.162", "CTW.PV.183", "CTW.PV.184", "CTW.PV.185", "CTW.PV.187", "CTW.PV.193", "CWE.PD.1", "CWE.PD.2", "CWE.PD.3", "CWE.PD.4", "CWE.PD.5", "CWE.PD.14", "CWE.PD.16", "CWE.PD.17", "FBR.PD.123", "FBR.PD.124", "FBR.PD.125", "FBR.PD.126", "FBR.PD.127", "FBR.PD.128", "FBR.PD.129", "FBR.PD.130", "FBR.PD.131", "FBR.PD.132", "FBR.PD.134", "FBR.PD.135", "FBR.PV.111", "FBR.PV.112", "FBR.PV.113", "FBR.PV.114", "FBR.PV.115", "FBR.PV.116", "FBR.PV.120", "FBR.PV.121", "FBR.PV.122", "FBR.PV.133", "FBR.PV.136", "FBR.PV.137", "FW.PD.142", "FW.PD.145", "FW.PD.146", "FW.PD.149", "FW.PD.165", "FW.PD.169", "FW.PD.173", "FW.PD.177", "FW.PD.178", "FW.PD.179", "FW.PD.180", "FW.PD.181", "FW.PV.138", "FW.PV.139", "FW.PV.140", "FW.PV.141", "FW.PV.144", "FW.PV.168", "FW.PV.170", "FW.PV.171", "FW.PV.174", "FW.PV.175", "FW.PV.176", "FW.PV.182", "GCR.PD.45", "GCR.PD.46", "GCR.PD.47", "GCR.PD.48", "GCR.PD.49", "GCR.PD.50", "GCR.PD.51", "GCR.PD.52", "GEW.PD.235", "GEW.PD.236", "GEW.PD.237", "GEW.PD.241", "GEW.PD.242", "GEW.PD.243", "GEW.PD.247", "GEW.PD.250", "GEW.PD.251", "GEW.PD.255", "GEW.PD.256", "GEW.PD.258", "GEW.PV.233", "GEW.PV.234", "GEW.PV.239", "GEW.PV.244", "GEW.PV.245", "GEW.PV.246", "GEW.PV.249", "GEW.PV.252", "GEW.PV.253", "GEW.PV.254", "GEW.PV.259", "GEW.PV.260", "KP.PV.74", "KP.PV.75", "KP.PV.76", "KP.PV.77", "KP.PV.78", "KP.PV.79", "MR.PD.18", "MR.PD.19", "MR.PD.20", "MR.PD.21", "MR.PD.23", "MR.PD.24", "MR.PD.25", "MR.PD.26", "MR.PD.27", "MR.PD.28", "MR.PD.29", "MR.PD.30", "MR.PD.31", "NW.PD.33", "NW.PD.34", "NW.PD.35", "NW.PD.36", "NW.PD.37", "NW.PD.38", "NW.PD.39", "NW.PD.41", "NW.PD.42", "NW.PD.43", "NW.PD.44", "PME.PV.224", "PME.PV.225", "PME.PV.226", "PME.PV.227", "PME.PV.228", "PME.PV.230", "PME.PV.231", "PMW.PV.212", "PMW.PV.213", "PMW.PV.214", "PMW.PV.215", "PMW.PV.216", "PMW.PV.217", "PMW.PV.218", "PMW.PV.219", "PMW.PV.220", "PMW.PV.221", "PMW.PV.222", "PMW.PV.223", "RCW.PV.262", "RCW.PV.263", "RCW.PV.264", "RCW.PV.265", "RCW.PV.266", "RCW.PV.267", "RCW.PV.269", "RCW.PV.270", "RCW.PV.271", "RCW.PV.272", "RCW.PV.273", "RCW.PV.274", "RCW.PV.275", "STR.PD.58", "STR.PD.59", "STR.PD.60", "STR.PD.61", "STR.PD.62", "STR.PD.63", "STR.PD.64", "STR.PD.65", "STR.PD.66", "STR.PD.67", "STR.PD.68", "STR.PD.69", "STR.PD.71", "WLR.PV.96", "WLR.PV.98", "WLR.PV.99", "WLR.PV.100", "WLR.PV.101", "WLR.PV.102", "WLR.PV.103", "WLR.PV.105", "WLR.PV.106", "WLR.PV.107", "WLR.PV.108", "WLR.PV.109")
only.guts.phyloseq <- prune_samples(only.gut.samples, all_samp_physeq) #Looks good. 254 samples!
dim(otu_table(only.guts.phyloseq)) #23641, 254
gut.otu.table.tp <- t(otu_table(only.guts.phyloseq))
dim(gut.otu.table.tp) #254 23641
class(gut.otu.table.tp)

# SOIL SAMPLES
only.soil.samples <-c("AN.S.1.2", "AN.S.3.2", "B2R.S.1.2", "B2R.S.3.2", "BCC.S.6.2", "BCC.S.10.2", "CTW.S.3.2", "CTW.S.5.2", "CWE.S.5.2", "CWE.S.11.2", "FBR.S.1.2", "FBR.S.2.2", "FW.S.1.2", "FW.S.6.2", "GCR.S.1.2", "GCR.S.3.2", "GEW.S.3.2", "GEW.S.5.2", "KP.S.4.2", "KP.S.6.2", "MR.S.1.2", "MR.S.4.2", "NW.S.2.2", "NW.S.6.2", "PME.S.6.2", "PME.S.9.2", "PMW.S.2.2", "PMW.S.3.2", "RCW.S.2.2", "RCW.S.5.2", "STR.S.1.2", "STR.S.3.2", "WLR.S.3.2", "WLR.S.6.2")
only.soil.phyloseq <- prune_samples(only.soil.samples, all_samp_physeq) #Looks good. 34 samples!
dim(otu_table(only.soil.phyloseq)) #23641, 34

# DROP GUT SAMPLES WITH LESS THAN 3500 READS AND SOIL SAMPLES WITH LESS THAN 26336 FOR RAREFYING
## GUT
gut.equal.to.or.greaterthan.3500 <-c("AN.PV.276", "AN.PV.277", "AN.PV.278", "AN.PV.279", "AN.PV.280", "AN.PV.281", "B2R.PV.309", "B2R.PV.310", "B2R.PV.311", "B2R.PV.312", "B2R.PV.318", "B2R.PV.319", "B2R.PV.321", "B2R.PV.322", "BCC.PD.285", "BCC.PD.288", "BCC.PD.291", "BCC.PD.292", "BCC.PD.295", "BCC.PD.296","BCC.PD.298", "BCC.PD.302", "BCC.PD.303", "BCC.PD.304", "BCC.PD.307", "BCC.PD.308", "BCC.PV.286", "BCC.PV.287", "BCC.PV.289", "BCC.PV.297","BCC.PV.300", "BCC.PV.301", "BR2.PD.324", "BR2.PD.325", "BR2.PD.326", "BR2.PD.327", "BR2.PD.328", "BR2.PD.329", "BR2.PD.330", "BR2.PD.332","BR2.PD.333", "BR2.PD.334", "BR2.PD.335", "BR2.PD.336", "CTW.PD.153", "CTW.PD.155", "CTW.PD.156", "CTW.PD.164", "CTW.PD.186", "CTW.PD.188","CTW.PD.191", "CTW.PD.192", "CTW.PV.157", "CTW.PV.158", "CTW.PV.162", "CTW.PV.183", "CTW.PV.184", "CTW.PV.185", "CTW.PV.193", "CWE.PD.14", "CWE.PD.16", "CWE.PD.17", "CWE.PD.2", "CWE.PD.3", "CWE.PD.4", "CWE.PD.5", "FBR.PD.123", "FBR.PD.124", "FBR.PD.125", "FBR.PD.126","FBR.PD.127", "FBR.PD.128", "FBR.PD.129", "FBR.PD.130", "FBR.PD.131", "FBR.PD.132", "FBR.PD.134", "FBR.PD.135", "FBR.PV.111", "FBR.PV.112","FBR.PV.113", "FBR.PV.114", "FBR.PV.115", "FBR.PV.116", "FBR.PV.120", "FBR.PV.121", "FBR.PV.122", "FW.PD.142", "FW.PD.146", "FW.PD.149", "FW.PD.165", "FW.PD.169", "FW.PD.173", "FW.PD.177", "FW.PD.178", "FW.PD.179", "FW.PD.180", "FW.PD.181", "FW.PV.138", "FW.PV.139", "FW.PV.141", "FW.PV.144", "FW.PV.168", "FW.PV.170", "FW.PV.171", "FW.PV.174", "FW.PV.175", "FW.PV.176", "FW.PV.182", "GCR.PD.45", "GCR.PD.46", "GCR.PD.47", "GCR.PD.48", "GCR.PD.49", "GCR.PD.52", "GEW.PD.235", "GEW.PD.236", "GEW.PD.237", "GEW.PD.241", "GEW.PD.242","GEW.PD.243", "GEW.PD.247", "GEW.PD.250", "GEW.PD.251", "GEW.PD.255", "GEW.PD.256", "GEW.PD.258", "GEW.PV.234", "GEW.PV.239", "GEW.PV.244","GEW.PV.245", "GEW.PV.246", "GEW.PV.249", "GEW.PV.252", "GEW.PV.253", "GEW.PV.254", "GEW.PV.259", "GEW.PV.260", "KP.PV.75", "KP.PV.78",  "KP.PV.79", "MR.PD.18", "MR.PD.19", "MR.PD.20", "MR.PD.23", "MR.PD.24", "MR.PD.25", "MR.PD.28", "MR.PD.29", "MR.PD.30",  "NW.PD.34", "NW.PD.35", "NW.PD.36", "NW.PD.37", "NW.PD.38", "NW.PD.41", "NW.PD.42", "NW.PD.43", "NW.PD.44", "PME.PV.224","PME.PV.225", "PME.PV.226", "PME.PV.228", "PME.PV.230", "PMW.PV.213", "PMW.PV.214", "PMW.PV.215", "PMW.PV.216", "PMW.PV.218", "PMW.PV.219", "PMW.PV.222", "RCW.PV.263", "RCW.PV.264", "RCW.PV.266", "RCW.PV.267", "RCW.PV.270", "RCW.PV.271", "RCW.PV.272", "RCW.PV.273", "RCW.PV.274","RCW.PV.275", "STR.PD.59", "STR.PD.60", "STR.PD.61", "STR.PD.62", "STR.PD.63", "STR.PD.64", "STR.PD.65", "STR.PD.66", "STR.PD.67", "STR.PD.68", "STR.PD.69", "STR.PD.71", "WLR.PV.102", "WLR.PV.106", "WLR.PV.107", "WLR.PV.109", "WLR.PV.98", "WLR.PV.99")
length(gut.equal.to.or.greaterthan.3500) #199
otu.gut.to.rarefy <- prune_samples(gut.equal.to.or.greaterthan.3500,only.guts.phyloseq)

## SOIL
soil.equaltoorgreaterthan26336 <-c("AN.S.1.2", "AN.S.3.2", "B2R.S.1.2", "B2R.S.3.2", "BCC.S.6.2", "BCC.S.10.2", "CTW.S.3.2", "CTW.S.5.2", "CWE.S.5.2", "CWE.S.11.2", "FBR.S.1.2", "FBR.S.2.2", "FW.S.1.2", "FW.S.6.2", "GCR.S.1.2", "GCR.S.3.2", "GEW.S.3.2", "GEW.S.5.2", "KP.S.4.2", "KP.S.6.2", "MR.S.1.2", "MR.S.4.2", "NW.S.2.2", "NW.S.6.2", "PME.S.6.2", "PME.S.9.2", "PMW.S.2.2", "PMW.S.3.2", "RCW.S.2.2", "RCW.S.5.2", "STR.S.1.2", "STR.S.3.2", "WLR.S.6.2")
length(soil.equaltoorgreaterthan26336) #33
otu.soil.to.rarefy <- prune_samples(soil.equaltoorgreaterthan26336, only.soil.phyloseq)

#################################################################################
#                 RAREFYING  MULTIPLE TIMES, JACCARD DIST, AND SHANNON
################################################################################

######################################### GUTS #########################################

# FUNCTION TO GET OTU OF PHYLOSEQ FOR RAREFYING (OR REALLY WORKING IN VEGAN IN GENERAL)
###### (switches rows and columns for vegan!)

psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

# GUT OTU TABLE OUT OF PHYLOSEQ
otu.gut.to.rarefy.veg <-psotu2veg(otu.gut.to.rarefy)
dim(otu.gut.to.rarefy.veg) #199 23641 (so columns are now ASVs)

# GETTING RID OF ASVs THAT ARE PRESENT IN NO GUT SAMPLES (to avoid memory problems with and to speed up rarefying)
otu.gut.to.rarefy.nozeros <- otu.gut.to.rarefy.veg[, colSums(otu.gut.to.rarefy.veg !=0) > 0]
dim(otu.gut.to.rarefy.nozeros) #199 rows and 1358 columns 
#check to make sure that you keep OTU names
colnames(otu.gut.to.rarefy.nozeros) 

# another way of doing this to check
otu.gut.to.rarefy.nozeros.sec.way <- otu.gut.to.rarefy.veg[,which(!apply(otu.gut.to.rarefy.veg, 2, FUN = function(x) {all(x == 0)}))]
dim(otu.gut.to.rarefy.nozeros.sec.way) #199 1358

set.seed(93)
# RAREFYING GUTS MULTIPLE TIMES AND THEN TAKING MEAN 
gut.r.3507.1k<-array(dim=c(199, 1358)) #199 rows, 1358 columns 
raw.rare<-list()
set.seed(93) 
for(j in 1:199){
  tempsamp<-array(dim=c(1000, 1358)) # rarefying 1000 times, columns match my number of columns 
  cat("\n",j,"of 199")
  for(i in 1:1000){ 
    tempsamp[i,]<-rrarefy(otu.gut.to.rarefy.nozeros[j,],3507) #3507 is where we want to rarefy to across gut samples
  }
  raw.rare[[j]]<-tempsamp
  gut.r.3507.1k[j,] <- apply(tempsamp,2,mean) # gets mean
}

dim(gut.r.3507.1k) #199 1358
sum(gut.r.3507.1k[,])/199 #equals 3507, so each sample was rarefied
class(gut.r.3507.1k)

# RE-APPLY ROW AND COLUMN NAMES
rownames(gut.r.3507.1k) <- rownames(otu.gut.to.rarefy.nozeros) 
colnames(gut.r.3507.1k) <- colnames(otu.gut.to.rarefy.nozeros) 

# STANDARDIZE BY PROPORTIONAL ABUNDANCE
gut.r.data <- decostand(gut.r.3507.1k, "total") # divided by row total, i.e. for each species
class(gut.r.data) # "matrix" "array"
dim(gut.r.data)
rownames(gut.r.data)

# SEPARATE P. VINDEX AND P. DIFFORMIS TO GET DIFFERENT DISTANCES
# P.VINDEX
PV.r.data.allcol <- gut.r.data[c(1:14, 27:32, 53:59, 79:87, 99:109, 128:141, 160:181, 194:199),]
dim(PV.r.data.allcol)  #89 1358
#filter out empty columns 
PV.r.data <- PV.r.data.allcol[, colSums(PV.r.data.allcol !=0) > 0]
dim(PV.r.data) #89, 855

#P. DIFFORMIS
PD.r.data.allcol <- gut.r.data[c(15:26, 33:52, 60:78, 88:98, 110:127, 142:159, 182:193),]
dim(PD.r.data.allcol) # 110 1358
#filter out empty columns
PD.r.data <- PD.r.data.allcol[, colSums(PD.r.data.allcol !=0) > 0]
dim(PD.r.data) #110, 907

#### GET SYMPATRY ONLY #####
row.names(gut.r.data)[7]
SYMP.r.data.allcol <- gut.r.data[c(7:59, 67:109,116:138), ]
dim(SYMP.r.data.allcol) #119, 1358
#filter out empty columns 
SYMP.r.data <- SYMP.r.data.allcol[, colSums(SYMP.r.data.allcol !=0) > 0]
dim(SYMP.r.data) #119, 943

#### GET ALLOPATRY ONLY 
ALLO.r.data.allcol <- gut.r.data[c(1:6, 60:66, 110:115, 139:199),]
dim(ALLO.r.data.allcol) #80, 1358
### filter out empty columns
ALLO.r.data <- ALLO.r.data.allcol[, colSums(ALLO.r.data.allcol !=0) > 0]
dim(ALLO.r.data) #80, 765

# GET JACCARD DISTANCE

# For all gut samples
gut.jacc.dist <- vegdist(gut.r.data, method = "jaccard")
str(gut.jacc.dist) #looks like we have sample names!

# For just P. vindex
PV.jacc.dist <- vegdist(PV.r.data, method = "jaccard")

# For just P. difformis
PD.jacc.dist <- vegdist(PD.r.data, method = "jaccard")
str(PD.jacc.dist)

# For just symp gut samples
SYMP.jacc.dist <- vegdist(SYMP.r.data, method = "jaccard")

# For just allo gut samples
ALLO.jacc.dist <- vegdist(ALLO.r.data, method = "jaccard")

# GET SHANNON DIVERSITY FROM RAREFIED GUTS (performed Aug. 31, )
soil_shann <- diversity(soil.r.26336.1k, index = "shannon", MARGIN = 1)


######################################### SOIL #########################################

# SOIL OTU TABLE OUT OF PHYLOSEQ
otu.soil.to.rarefy.veg <-psotu2veg(otu.soil.to.rarefy)
dim(otu.soil.to.rarefy.veg) # 33 23641 (so columns are now ASVs)

# GETTING RID OF ASVs THAT ARE PRESENT IN NO SOIL SAMPLES (to avoid memory problems and to speed up rarefying)
otu.soil.to.rarefy.nozeros <- otu.soil.to.rarefy.veg[, colSums(otu.soil.to.rarefy.veg !=0) > 0]
dim(otu.soil.to.rarefy.nozeros) #33 22365
#check to make sure that you keep OTU names
colnames(otu.soil.to.rarefy.nozeros) 

set.seed(93)

# RAREFYING SOIL MULTIPLE TIMES AND THEN TAKING MEAN 
soil.r.26336.1k<-array(dim=c(33, 22365)) #33 rows, 22365 columns in otu.soil.to.rarefy.nozeros
raw.rare<-list()
set.seed(93) 
for(j in 1:33){
  tempsamp<-array(dim=c(1000, 22365)) # rarefying 1000 times, columns match my number of columns 
  cat("\n",j,"of 33")
  for(i in 1:1000){ 
    tempsamp[i,]<-rrarefy(otu.soil.to.rarefy.nozeros[j,], 26336) #26336 is where we want to rarefy to. 
  }
  raw.rare[[j]]<-tempsamp
  soil.r.26336.1k[j,] <- apply(tempsamp,2,mean) #this step takes the mean. 
}

dim(soil.r.26336.1k) #33 22365
sum(soil.r.26336.1k[,])/33 #26336, means we rarefied correctly!

# RE-APPLY ROW AND COLUMN NAMES
rownames(soil.r.26336.1k) <- rownames(otu.soil.to.rarefy.nozeros) 
colnames(soil.r.26336.1k) <- colnames(otu.soil.to.rarefy.nozeros) 

# GET SHANNON DIVERSITY FROM RAREFIED SOIL
soil_shann <- diversity(soil.r.26336.1k, index = "shannon", MARGIN = 1)

#################################################################################
#                 CREATE NEW PHYLOSEQ OBJECTS WITH STUFF MADE ABOVE
################################################################################

###### 1. ALL GUT SAMPLES ######
# FIRST TRANSPOSE NEW RAREFIED GUT OTU TABLE SO TAXA = ROWS
gut.r.data.tp <-t(gut.r.data)
dim(gut.r.data.tp) # 1358  199

# MAKE TRANSPOSED OTU TABLE PHYLOSEQ OTU TABLE
gut_r_OTU = otu_table(gut.r.data.tp, taxa_are_rows = TRUE)

# READ IN NEW CSV FOR METADATA THAT HAS CORRECT SHANNON DATA COLLECTED ABOVE
metadata.shann = read.csv("phyloseq/all_samples/metadata_total_forR_csv.CSV", row.names=1) #need row names = 1 so that OTU names are consistent across objects
head(metadata.shann) 
dim(metadata.shann) #288, 20

gut_r_META = sample_data(metadata.shann)
#	(tree	was	already	imported	as	a	phyloseq	object)

#	check	that	your	OTU	names	are	consistent	across	objects
###yes, they are all the taxon ID, as opposed to the name.
taxa_names(TAX)
taxa_names(gut_r_OTU)
taxa_names(phy_tree) 

#	MERGE INTO ONE PHYLOSEQ OBJECT
gut.r.phyloseq<-	phyloseq(gut_r_OTU,	TAX,	gut_r_META,	phy_tree)
gut.r.phyloseq
dim(otu_table(gut.r.phyloseq))

#####
###### 2. ONLY P.VINDEX SAMPLES ######

# FIRST TRANSPOSE NEW RAREFIED GUT OTU TABLE SO TAXA = ROWS
PV.r.data.tp <-t(PV.r.data)
dim(PV.r.data.tp) # 855  89

# MAKE TRANSPOSED OTU TABLE PHYLOSEQ OTU TABLE
PV_r_OTU = otu_table(PV.r.data.tp, taxa_are_rows = TRUE)

#	check	that	your	OTU	names	are	consistent	across	objects
###yes, they are all the taxon ID, as opposed to the name.
taxa_names(TAX)
taxa_names(PV_r_OTU)
taxa_names(phy_tree)

#	MERGE INTO ONE PHYLOSEQ OBJECT
PV.r.phyloseq<-	phyloseq(PV_r_OTU,	TAX,	gut_r_META,	phy_tree)
PV.r.phyloseq #855 taxa and 89 samples. What we expect

#####
###### 3. ONLY P.DIFFORMIS SAMPLES ######

# FIRST TRANSPOSE NEW RAREFIED GUT OTU TABLE SO TAXA = ROWS
PD.r.data.tp <-t(PD.r.data)
dim(PD.r.data.tp) # 907 110

# MAKE TRANSPOSED OTU TABLE PHYLOSEQ OTU TABLE
PD_r_OTU = otu_table(PD.r.data.tp, taxa_are_rows = TRUE)

#	check	that	your	OTU	names	are	consistent	across	objects
###yes, they are all the taxon ID, as opposed to the name.
taxa_names(TAX)
taxa_names(PD_r_OTU)
taxa_names(phy_tree)

#	MERGE INTO ONE PHYLOSEQ OBJECT
PD.r.phyloseq<-	phyloseq(PD_r_OTU,	TAX,	gut_r_META,	phy_tree)
PD.r.phyloseq #907 taxa and 110 samples, WHAT WE expect


###### 4. ONLY SYMPATRIC SAMPLES ######
# FIRST TRANSPOSE NEW RAREFIED GUT OTU TABLE SO TAXA = ROWS
SYMP.r.data.tp <-t(SYMP.r.data)
dim(SYMP.r.data.tp) # 943  119

# MAKE TRANSPOSED OTU TABLE PHYLOSEQ OTU TABLE
SYMP_r_OTU = otu_table(SYMP.r.data.tp, taxa_are_rows = TRUE)

SYMP_r_META = sample_data(metadata.shann)
#	(tree	was	already	imported	as	a	phyloseq	object)

#	check	that	your	OTU	names	are	consistent	across	objects
###yes, they are all the taxon ID, as opposed to the name.
taxa_names(TAX)
taxa_names(SYMP_r_OTU)
taxa_names(phy_tree) 

#	MERGE INTO ONE PHYLOSEQ OBJECT
SYMP.r.phyloseq<-	phyloseq(SYMP_r_OTU,	TAX,	SYMP_r_META,	phy_tree)
SYMP.r.phyloseq
sample_data(SYMP.r.phyloseq) # looks good!

#################################################################################
#                      GET NEW, POST RAREFACTION UNIFRAC DISTANCES
################################################################################

# GUT UNIFRAC DISTANCES
guts.wUF.dist <- UniFrac(gut.r.phyloseq, weighted = TRUE)
guts.uwUF.dist <-UniFrac(gut.r.phyloseq, weighted = FALSE)

# P. VINDEX UNIFRAC
PV.wUF.dist <- UniFrac(PV.r.phyloseq, weighted = TRUE)
PV.uwUF.dist <-UniFrac(PV.r.phyloseq, weighted = FALSE)

# P. DIFFORMIS UNIFRAC
PD.wUF.dist <- UniFrac(PD.r.phyloseq, weighted = TRUE)
PD.uwUF.dist <-UniFrac(PD.r.phyloseq, weighted = FALSE)

# SYMPATRY ONLY UNIFRAC
SYMP.wUF.dist <- UniFrac(SYMP.r.phyloseq, weighted = TRUE)
SYMP.uwUF.dist <-UniFrac(SYMP.r.phyloseq, weighted = FALSE)

##### USEFUL FUNCTION (to make phyloseq objects compatible with vegan)
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

gut_r_vg_sampdat_catord_shavg <- read.csv(file = "gut_r_vg_sampdat_catord_shavg.csv", row.names = 1, header = TRUE)

Ord_Cattle_PresenceRank <- factor(gut_r_vg_sampdat_catord_shavg$CattlePresenceRank, levels = c("low", "medium", "high"), labels = c("low", "medium", "high"), ordered = TRUE)
length(Ord_Cattle_PresenceRank)

# Add in ordinal cattle variable because was messed up before
gut_r_veg_final_sampdat <- cbind(gut_r_vg_sampdat_catord_shavg[,c(1:8, 10:21)], Ord_Cattle_PresenceRank)

colnames(gut_r_veg_final_sampdat)

gut_r_veg_final_sampdat$Ord_Cattle_PresenceRank # gut_r_veg_final_sampdat has ordinal cattle variables 

#################################################################################
#                      SAVE OBJECTS AS R DATA FILES
################################################################################
##### re-saved June 24 when we made sure that cattle ordinal variables were correct
save(gut_r_veg_final_sampdat, all_samp_physeq_raw, psotu2veg, pssd2veg, gut.r.3507.1k, gut.r.data, gut.jacc.dist, soil.r.26336.1k, soil_shann, gut.r.phyloseq, guts.wUF.dist, guts.uwUF.dist, PV.r.data, PD.r.data, SYMP.r.data, ALLO.r.data, PV.jacc.dist, PD.jacc.dist, SYMP.jacc.dist, ALLO.jacc.dist, PV.r.phyloseq, PD.r.phyloseq, PV.wUF.dist, PV.uwUF.dist, PD.wUF.dist, PD.uwUF.dist, file = "ps_rare_diss_shann_.RData")