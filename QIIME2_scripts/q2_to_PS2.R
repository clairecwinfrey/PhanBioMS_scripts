####################################################################################################
# q2_to_PS2.R
# Takes files produced in 2.	q2_to_PS1.R and makes them phyloseq objects in R
# based on https://github.com/rwmurdoch/Manneheimia/blob/master/qiime2_to_physeq2.R
# Data: 2019-Biogeo_Phanaeus_GM
# April 26, 2020
####################################################################################################

library("ggplot2")
library("phyloseq")
library("ape")

#######################################################################################

otu.table.all = read.csv(file = "phyloseq/nonrarefied_OTU_matrix.csv", sep=",", row.names=1)
otu.table.all = as.matrix(otu.table.all)
head(otu.table.all)

taxonomy.all = read.csv("phyloseq/taxonomy_only.csv", sep=",", row.names = 1)
taxonomy.all = as.matrix(taxonomy.all)
head(taxonomy.all)

metadata.no.shan = read.csv("phyloseq/all_samples/metadata_nosoilshannon_forR_csv.csv", row.names=1) #need row names = 1 so that OTU names are consistent across objects
#for this, removed q2 "categories" label
#has all samples and all metadata, EXCEPT soil Shannon stuff
head(metadata.no.shan)

phy_tree = read_tree("phyloseq/tree.nwk") #Yay! It worked!

OTU = otu_table(otu.table.all, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy.all)
META = sample_data(metadata.no.shan)
#	(tree	was	already	imported	as	a	phyloseq	object)

#	CHECKING THAT OTU NAMES CONSISTENT
###yes, they are all the taxon ID, as opposed to the name.
taxa_names(TAX)
taxa_names(OTU)
taxa_names(phy_tree)

#	CHECK THAT SAMPLE NAMES ARE CONSISTENT
# Both use . not -.
sample_names(OTU)
sample_names(META)

#	MERGE INTO ONE PHYLOSEQ OBJECT
all_samp_physeq_raw	<-	phyloseq(OTU,	TAX,	META,	phy_tree)
dim(otu_table(all_samp_physeq_raw))
#still 23641   288,
######what phyloseq looks like!### numbers looooook goooood!!!

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 23641 taxa and 288 samples ]
#sample_data() Sample Data:       [ 287 samples by 14 sample variables ]
#tax_table()   Taxonomy Table:    [ 23641 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 23641 tips and 23640 internal nodes ]

sample_data(all_samp_physeq_raw)

save(otu.table.all, taxonomy.all, metadata.no.shan, phy_tree, OTU, TAX, META, all_samp_physeq_raw, file ="raw_phyloseq.dat.RData")
