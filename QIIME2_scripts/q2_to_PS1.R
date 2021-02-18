#q2_to_PS1.R

# This script is the first step at converting qiime2 files to an R-ready format.
# code adopted from https://github.com/rwmurdoch/Manneheimia/blob/master/qiime2_to_physeq1.R
# And based on this thread: https://forum.qiime2.org/t/converting-biom-files-with-taxonomic- info-for-import-in-r-with-phyloseq/2542/5

##### THINGS TO DO MANUALLY PRIOR TO RUNNING THIS SCRIPT ####
# In qiime_allsamples.sh, I exported an ASV table (r3500_final_feat_table_taxIDs.txt).
# Before running the following script, I changed #OTUID to OTUID in a text editor.
# r3500_final_feat_table_taxIDs.txt was then saved as r3500_OTU_table.txt, as used below.

otu <- read.table(file = "phyloseq/r3500_OTU_table.txt", header = TRUE)
head(otu)

# Before running the script below, I changed phyloseq/taxonomy_nospace.tsv to phyloseq/taxonomy_nospace.csv in Excel
tax <- read.csv(file="phyloseq/taxonomy_nospace.csv",header = TRUE)

head(tax)

# Because I filtered mitochondria, chloroplasts, etc. from my OTU table in qiime_allsamples.sh,
# and qiime2 doesn't filter out taxonomy, I had to merge these two files in R so that the length
# is right
merged_file <- merge(otu, tax, by.x = c("OTUID"), by.y = c("OTUID"))
head(merged_file)

write.table(merged_file, file = "combined_otu_tax_r3500.txt", sep = '\t', col.names = TRUE, row.names = FALSE)

# Now, open txt file that was just made in Excel and split into two files:
#1) for taxonomy only (only OTUID and taxonomic info columns)
######taxonomy_only_r3500.csv
#and 2)for the OTU matrix (containing only OTUID and abundances in each sample.)
###### OTU_matrix_r3500.csv
