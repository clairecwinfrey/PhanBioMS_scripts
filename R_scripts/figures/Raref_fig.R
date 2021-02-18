######## DRAWING RAREFACTION CURVES IN R ###########

library(vegan)
load("ps_rare_diss_shann_.RData") # Created in PS_Rarefying_Dissim.R

setwd("/Users/clairewinfrey/Desktop/2019_biogeo_phan_GM/_All_pops")

length(rownames(sample_data(all_samp_physeq_raw)))

samp.col <- c(rep("blue", 9), rep("saddlebrown", 2), rep("blue", 12), rep("saddlebrown", 2), rep("gold", 12), rep("blue", 7), rep("saddlebrown", 2), rep("gold", 12), rep("gold", 12), 
              rep("blue", 11), rep("saddlebrown", 2), rep("gold", 8), rep("saddlebrown", 2), "black", rep("gold", 12),  rep("blue", 12), rep("saddlebrown", 2), rep("gold", 12), rep("blue", 12),
              rep("saddlebrown", 2), rep("gold", 8), rep("saddlebrown", 2), rep("gold", 12),rep("blue", 12), rep("saddlebrown", 2), rep("blue", 6), rep("saddlebrown", 2), rep("gold", 13), 
              rep("saddlebrown", 2), rep("gold", 11), rep("saddlebrown", 2), rep("black", 3), rep("blue", 7), rep("saddlebrown", 2), rep("blue", 12), rep("saddlebrown", 2), 
              rep("blue", 13), rep("saddlebrown", 2), rep("gold", 13), rep("saddlebrown", 2), rep("blue", 12), rep("saddlebrown", 2)) 

length(samp.col)

rare.plot <- rarecurve(t(otu_table(all_samp_physeq_raw)), step = 500, cex = 0.5, col = samp.col, label = FALSE, xlab = "Number of Sequences", ylab = "Number of ASVs")
