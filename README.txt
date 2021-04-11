Biogeography of Phanaeus gut microbiome Project Scripts Overview
This is a description of how each file in scripts folder was used in our data analysis.

qiime_allsamples.sh
--This script was used in the qiime2 platform (version 2020.2). Specifically, it is the code for importing raw sequences, primer removal, dada2, getting rarefaction curves (just to try different depths NOT version used in paper), training the feature classifier, assigning taxonomy (SILVA 128 database), performing ASV fragment insertion on the SILVA 128 reference tree, filtering tree and ASV tables, and exporting tree and ASV tables.

q2_to_PS1.R
--This script is the first of two that takes qiime2 produced files and adapts them for use in R. 
Based on script here: https://github.com/rwmurdoch/Manneheimia/blob/master/qiime2_to_physeq1.R

q2_to_PS2.R
--This script takes files produced in q2_to_PS1.R and converts them into phyloseq objects
Based on script here: https://github.com/rwmurdoch/Manneheimia/blob/master/qiime2_to_physeq2.R

PS_Rarefying_Dissim.R
--This script finishes setting up phyloseq objects (based on various subsets of the data), rarefies the dataset, and gets dissimilarities 

Raref_fig.R
--This script creates rarefaction plot based on 3,500 rarefaction depth (for gut samples)

Taxonomic_barplots.R
--Script makes stacked barplots for relative abundance of different phyla and families in beetles’ guts (by Phanaeus species)

Boxplot_dist_comp_Nov1.R ---replace with Boxplot_comp_Dec27.R
--This script extracts dissimilarity comparisons for different combinations of P. vindex and P. difformis individuals in sympatry and in allopatry. Makes figure of boxplots

indicator_sp_analysis.R
--This script performs the species indicator analyses and obtains relative abundance for each identified bacterial ASV. First analysis is based on sympatry or allopatry and second analysis is based on the abundance of cattle present at each location

