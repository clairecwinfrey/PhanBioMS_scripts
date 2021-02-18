#!/bin/bash

####################################################################################################
# Bioinformatics processing of Phanaeus vindex and P. difformis gut microbiome data
# Beetles caught summer 2019, lab work and sequencing done fall 2019.
# March 19, 2020 (script worked on through April 2020)
# Claire Winfrey
####################################################################################################

# source activate qiime2-2020.2

####################################################################################################
#                                IMPORTING SEQUENCES
####################################################################################################
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path _manifest_files/All_samples_q2-manifest.csv \
    --output-path _All_pops/import_all_pops.qza \
    --input-format PairedEndFastqManifestPhred33 \

# Visualize sequence quality
qiime demux summarize \
--i-data _All_pops/import_all_pops.qza \
--o-visualization _All_pops/import_all_pops.qzv \

# results: 11,130,657 sequence counts! Median per sample: 37458.0. Range: 774-96086

####################################################################################################
#                                  PRIMER REMOVAL
####################################################################################################

# primers are 515F and 806R
qiime cutadapt trim-paired \
--i-demultiplexed-sequences import_all_pops.qza \
--p-cores 10 \
--p-front-f GTGYCAGCMGCCGCGGTAA \
--p-front-r GGACTACNVGGGTWTCTAAT \
--o-trimmed-sequences all-cutadapt.qza \
--p-error-rate 0.2 \
--verbose > all_cutadapt_stats.txt

###########Visualize
qiime demux summarize \
--i-data all-cutadapt.qza \
--o-visualization all-cutadapt.qzv \

####################################################################################################
#                          DADA2: QC AND FEATURE TABLE CONSTRUCTION
####################################################################################################
qiime dada2 denoise-paired \
--i-demultiplexed-seqs all-cutadapt.qza \
--p-trunc-len-f 185 \
--p-trunc-len-r 170 \
--p-n-threads 0 \
--o-table dada2_output/dada2_table.qza \
--o-representative-sequences dada2_output/dada2_rep_seq.qza \
--o-denoising-stats dada2_output/dada2_denoising-stats.qza \

# Get denoising stats
qiime tools export \
--input-path dada2_output/dada2_denoising-stats.qza \
--output-path exports/dada2_denoising-stats.txt \

##representatiave sequences: each feature in feature table will be represented by exactly one sequence (the joined-end sequences)

####################################################################################################
#                               VISUALIZING DADA2 RESULTS
####################################################################################################
# Create visualization of feature table and download sample frequency per sample csv
qiime feature-table summarize \
--i-table dada2_output/dada2_table.qza \
--o-visualization dada2_output/dada2_table.qzv \
####NEW RESULTS: #33,488 features and frequency 6,442,702

# Create visualization of feature table with metadata and download sample frequency per sample csv
qiime feature-table summarize \
--i-table dada2_output/dada2_table.qza \
--m-sample-metadata-file dada2_table_with_metadata.tsv \
--o-visualization dada2_output/metadata_table.qzv \
#had 33,488 features, with the total frequency of 6,442,702

# Representative sequences
qiime feature-table tabulate-seqs \
--i-data  dada2_output/dada2_rep_seq.qza \
--o-visualization dada2_output/dada2_rep_seq.qzv \

qiime metadata tabulate \
--m-input-file dada2_output/dada2_denoising-stats.qza \
--o-visualization dada2_output/dada2_denoising-stats.qzv \
#average percent of non-chimeric reads 51.39559028 (st dev 23.1853396)...

####################################################################################################
#                           OBTAINING & TRAINING THE FEATURE CLASSIFIER
####################################################################################################
####because qiime2 only has pre-trained classifer on my dataset for SILVA132

# Get scikit-learn 0.22.1. Can ONLY use these classifiers with scikit-learn 0.22.1
conda install --override-channels -c defaults scikit-learn=0.22.1

# Download and get Silva128 collection of reference files
wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_128_release.tgz
tar -xvf Silva_128_release.tgz
rm Silva_128_release.tgz

# Import 99% 16S only dataset
qiime tools import \
--type FeatureData[Sequence] \
--input-path /Users/brc_user/Desktop/CCW_Research/2019_biogeo_phan_GM/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/99/99_otus_16S.fasta \
--output-path /Users/brc_user/Desktop/CCW_Research/2019_biogeo_phan_GM/99_reps.qza

qiime tools import \
--type FeatureData[Taxonomy] \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path /Users/brc_user/Desktop/CCW_Research/2019_biogeo_phan_GM/SILVA_128_QIIME_release/taxonomy/16S_only/99/consensus_taxonomy_7_levels.txt \
--output-path /Users/brc_user/Desktop/CCW_Research/2019_biogeo_phan_GM/99_ref_taxonomy.qza

# Extract target region based on primer set
qiime feature-classifier extract-reads \
--i-sequences /Users/brc_user/Desktop/CCW_Research/2019_biogeo_phan_GM/99_reps.qza \
--p-f-primer GTGYCAGCMGCCGCGGTAA \
--p-r-primer GGACTACNVGGGTWTCTAAT \
--p-min-length 100 \
--p-max-length 400 \
--o-reads /Users/brc_user/Desktop/CCW_Research/2019_biogeo_phan_GM/99silva-ref-seqs.qza
--verbose \

# Train model
#note, used scikit-learn version 0.22.1
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads /Users/brc_user/Desktop/CCW_Research/2019_biogeo_phan_GM/99silva-ref-seqs.qza \
--i-reference-taxonomy /Users/brc_user/Desktop/CCW_Research/2019_biogeo_phan_GM/99_ref_taxonomy.qza \
--o-classifier /Users/brc_user/Desktop/CCW_Research/2019_biogeo_phan_GM/99silva_classifier.qza \
--verbose \

####################################################################################################
#                           TAXONOMIC CLASSIFICATION
####################################################################################################

qiime feature-classifier classify-sklearn \
--i-classifier /Users/brc_user/Desktop/CCW_Research/2019_biogeo_phan_GM/_All_pops/trained_silva_classifer/99silva_classifier.qza \
--p-n-jobs 6 \
--i-reads dada2_output/dada2_rep_seq.qza \
--o-classification seq_taxonomy.qza \
#TAXONOMY.QZA IS TYPE: FEATUREDATA[TAXONOMY] AND FORMAT:"TSVTaxonomyDirectoryFormat"

####################################################################################################
#                                  EXPORTING TAXONOMY
####################################################################################################
#SEE EXPORT QIIME2 TUTORIAL
#export taxonomy file as TSVTaxonomyDirectoryFormat (BUT WE'LL HAVE TO REEXPORT ONCE FILTERED)
qiime tools export \
--input-path seq_taxonomy.qza \
--output-path exports

####NEW RESULTS: #33,488
###33488 total lines of code (matches number of features of dada2_table.qzv); 32,534 assigned to at least kingdom,
#23,920 assigned to phylum; to species: 21,004 to family level; 18487 genera; 12,801 species!!!

####results: all 19,259 ind feat post dada2 were in it. Of these, 18445 were assigned to at least kingdom
#11,572 to phylum, class= 11,378, order=10,644; family=10,322; genus=8,484; 5142 to species
###this matches the number in the table immediately after dada2
###this also matches the number in the SEPP-filtered-table.qzv

#export feature table
#######################################################################################################
#                                  Q2 FRAGMENT INSERTION PLUGIN
######################################################################################################
##following steps from plugin tutorial: https://library.qiime2.org/plugins/q2-fragment-insertion/16/
#got SILVA 128 tree from: https://docs.qiime2.org/2020.2/data-resources/

qiime fragment-insertion sepp \
--i-representative-sequences dada2_output/dada2_rep_seq.qza \
--i-reference-database /Users/brc_user/Desktop/CCW_Research/2019_biogeo_phan_GM/sepp-refs-silva-128.qza \
--o-tree frag_insert_output/insertion-tree.qza \
--o-placements frag_insert_output/insertion-placements.qza \
--p-threads 11 \
--verbose \

# Filter feature-tables so that it only contains fragments that are in the insertion tree.
qiime fragment-insertion filter-features \
--i-table  dada2_output/dada2_table.qza \
--i-tree frag_insert_output/insertion-tree.qza \
--o-filtered-table frag_insert_output/sepp-filtered_table.qza \
--o-removed-table frag_insert_output/sepp-removed_table.qza \
--verbose \

# Visualize sepp filtered table
qiime feature-table summarize \
--i-table frag_insert_output/sepp-filtered_table.qza \
--o-visualization frag_insert_output/sepp-filtered_table.qzv \
###this shows that 33,488 features remained after SEPP (i.e. were incorporated into tree). Total frequency is 6,442,702!
###last time with wrong input file (rep seq from dada2 run on seqs where primers had not been trimmed off), this was 19,259 features across 2,188,804 total reads.

# Visualize sepp removed table
qiime feature-table summarize \
--i-table frag_insert_output/sepp-removed_table.qza \
--o-visualization frag_insert_output/sepp-removed_table.qzv \
###this shows that no features are in this table, i.e. they all made it into filtered table!

####################################################################################################
#                           FILTERING SEPP FILTERED FEATURE TABLES
####################################################################################################

# Removing unassigned at phylum level, chloroplasts, and mitochondria
qiime taxa filter-table \
--i-table frag_insert_output/sepp-filtered_table.qza \
--i-taxonomy seq_taxonomy.qza \
--p-include D_1 \
--p-exclude mitochondria,chloroplast \
--o-filtered-table filtered_post_sepp/sepp_table_with_phyla_no_mito-no-chloro.qza \
#from taxonomy.tsv, seems that 23920 have D_1 assigned, 31 have chloroplasts, and 240 have mitochondria.
#So, we expect 23920 minus 31 minus 240 = 23649 after filtering.

## Let's check out this table:
qiime feature-table summarize \
--i-table filtered_post_sepp/sepp_table_with_phyla_no_mito-no-chloro.qza \
--o-visualization filtered_post_sepp/sepp_table_with_phyla_no_mito-no-chloro.qzv \
#---23,649 features and 6,034,491 total frequency
#matches what we expect

qiime feature-table filter-features \
--i-table filtered_post_sepp/sepp_table_with_phyla_no_mito-no-chloro.qza \
--p-min-frequency 2 \
--o-filtered-table filtered_post_sepp/sepp_table_with_phyla_no_mito-no-chloro-nodoubl.qza \

# Make visualization on the final table version so we can explore rarefying
qiime feature-table summarize \
--i-table filtered_post_sepp/sepp_table_with_phyla_no_mito-no-chloro-nodoubl.qza \
--o-visualization filtered_post_sepp/sepp_table_with_phyla_no_mito-no-chloro-nodoubl.qzv \

####################################################################################################
#                           RAREFACTION CURVES
####################################################################################################

###############################################
#####make alpha rarefaction curves to see at what depth I should rarefy:

#by default will do observed_otus, shannon, and bc i have phylogeny, faith_pd
qiime diversity alpha-rarefaction \
--i-table filtered_post_sepp/sepp_table_with_phyla_no_mito-no-chloro-nodoubl.qza \
--i-phylogeny frag_insert_output/insertion-tree.qza \
--p-max-depth 5000 \
--p-metrics faith_pd \
--p-metrics shannon \
--p-metrics chao1 \
--p-metrics observed_otus \
--m-metadata-file all_samples_metadata.tsv \
--o-visualization rarefaction.qzv \

# let's look at greater sampling depth, so that we have a better idea of when soil samples peak.

qiime diversity alpha-rarefaction \
--i-table filtered_post_sepp/sepp_table_with_phyla_no_mito-no-chloro-nodoubl.qza \
--i-phylogeny frag_insert_output/insertion-tree.qza \
--p-max-depth 27000 \
--p-metrics faith_pd \
--p-metrics shannon \
--p-metrics chao1 \
--p-metrics observed_otus \
--m-metadata-file /Users/clairewinfrey/Desktop/2019_biogeo_phan_GM/_All_pops/qiime_metadata_8col.tsv \
--o-visualization rarefaction.27000.qzv \


###FAITH PD
#For range overlap looks good at about 3500. Not quite plateaued, but pretty flat.
#Soil depth curve was slowing but not flat at 3500. To plateau, looks like we'd have to go deeper than 5000
#Gut were flat around 1500,

####################################################################################################
#                                   RAREFYING
####################################################################################################
# We ended up rarefying to 3500 sequences (but tried 1500, 2000 too first (not shown in this script))

# Rarefy to 3500 seqs/sample:
qiime feature-table rarefy \
--i-table filtered_post_sepp/sepp_table_with_phyla_no_mito-no-chloro-nodoubl.qza \
--p-sampling-depth 3500 \
--o-rarefied-table rarefied/table_with_phyla_no_mito-no-chloro-nodoubl_r3500.qza \

# Take a look at rarefied table:
qiime feature-table summarize \
--i-table rarefied/table_with_phyla_no_mito-no-chloro-nodoubl_r3500.qza \
--o-visualization rarefied/table_with_phyla_no_mito-no-chloro-nodoubl_r3500.qzv \

######retained 16,743 features, 812,000 sequences (3500 seqs * 232 samples)

####################################################################################################
#                                    RAREFIED BARPLOTS
####################################################################################################

qiime taxa barplot \
--i-table rarefied/table_with_phyla_no_mito-no-chloro-nodoubl_r3500.qza \
--i-taxonomy seq_taxonomy.qza \
--m-metadata-file all_samples_metadata.tsv \
--o-visualization barplots/r3500-barplots.qzv \

# Preparing qiime2 files for import into phyloseq
# April 26, 2020
#following Robert Murdoch's instructions here: https://github.com/rwmurdoch/Manneheimia/blob/master/README.md
####################################################################################################

####################################################################################################
#                                  EXPORT ASV TABLE
####################################################################################################
qiime tools export \
--input-path rarefied/table_with_phyla_no_mito-no-chloro-nodoubl_r3500.qza \
--output-path phyloseq \

#####convert to tsv table
biom convert \
-i phyloseq/feature-table.biom \
-o phyloseq/r3500_OTU_table.txt \
--to-tsv \

####now, open up the otu_table.txt in a text editor and change #OTUID to OTUID (was #OTU ID, now is OTUID)

####################################################################################################
#                                  EXPORT TAXONOMY
####################################################################################################

#####Now, export taxonomy table
qiime tools export \
--input-path seq_taxonomy.qza \
--output-path phyloseq \
#is called taxonomy.tsv

####open this up in a text editor and change feature ID to OTUID

####################################################################################################
#                                     EXPORT TREE
####################################################################################################

###this tree is rooted
qiime tools export \
--input-path frag_insert_output/insertion-tree.qza \
--output-path phyloseq \

######now, we'll need to merge taxonomy and OTU tables in R and then output a
#merged file.

####### Next see R file: q2_to_PS1.R
