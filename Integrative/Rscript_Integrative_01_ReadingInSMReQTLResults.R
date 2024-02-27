#Comparing the eQTL/SMR results to our F0 and F2 DE Results
#Megan Hagenauer
#09-12-2023, additions/edits up through 02-05-2024
#updated to include results from F2 juvenile QTLS as additional validation

###################

sessionInfo()
# 
# R version 4.2.3 (2023-03-15)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Ventura 13.6
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] qqman_0.1.8                 plyr_1.8.8                  vcfR_1.14.0                 SummarizedExperiment_1.28.0
# [5] Biobase_2.58.0              GenomicRanges_1.50.2        GenomeInfoDb_1.34.9         IRanges_2.32.0             
# [9] S4Vectors_0.36.2            BiocGenerics_0.44.0         MatrixGenerics_1.10.0       matrixStats_1.0.0          
# 
# loaded via a namespace (and not attached):
#   [1] httr_1.4.6             pkgload_1.3.2          bit64_4.0.5            viridisLite_0.4.1     
# [5] splines_4.2.3          shiny_1.7.4            memuse_4.2-3           GenomeInfoDbData_1.2.9
# [9] remotes_2.4.2.1        sessioninfo_1.2.2      lattice_0.21-8         glue_1.6.2            
# [13] digest_0.6.31          promises_1.2.0.1       XVector_0.38.0         htmltools_0.5.5       
# [17] httpuv_1.6.9           Matrix_1.5-4.1         devtools_2.4.5         calibrate_1.7.7       
# [21] zlibbioc_1.44.0        purrr_1.0.1            xtable_1.8-4           processx_3.8.0        
# [25] later_1.3.0            mgcv_1.9-0             usethis_2.2.2          ellipsis_0.3.2        
# [29] cachem_1.0.8           cli_3.6.1              magrittr_2.0.3         crayon_1.5.2          
# [33] mime_0.12              memoise_2.0.1          ps_1.7.4               fs_1.6.1              
# [37] nlme_3.1-162           MASS_7.3-60            pkgbuild_1.4.2         vegan_2.6-4           
# [41] profvis_0.3.8          tools_4.2.3            data.table_1.14.8      prettyunits_1.1.1     
# [45] lifecycle_1.0.3        stringr_1.5.0          cluster_2.1.4          DelayedArray_0.24.0   
# [49] callr_3.7.3            compiler_4.2.3         rlang_1.1.1            grid_4.2.3            
# [53] RCurl_1.98-1.12        rstudioapi_0.14        htmlwidgets_1.6.2      miniUI_0.1.1.1        
# [57] bitops_1.0-7           R6_2.5.1               pinfsc50_1.2.0         fastmap_1.1.1         
# [61] gemma.R_1.3.2          bit_4.0.5              permute_0.9-7          ape_5.7-1             
# [65] stringi_1.7.12         parallel_4.2.3         Rcpp_1.0.10            vctrs_0.6.3           
# [69] urlchecker_1.0.1 

############

#Reading in relevant files:

setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/F2_eQTLs")

F0_Meta_F2_DEResults<-read.csv("TableS3_F0_Meta_F2_results_wLineageGenes_ForSuppl_20230419.csv", header=TRUE, stringsAsFactors = FALSE)
str(F0_Meta_F2_DEResults)
#Gene identifier column = ENSEMBL.ID..Rnor6.Ensembl.v.103.
tail(F0_Meta_F2_DEResults)
#There are two nonsense rows (all NA) tagged on the end
F0_Meta_F2_DEResults<-F0_Meta_F2_DEResults[-c(13787:13788),]

F2_eqtls<-read.delim("f2.eqtls_indep.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)

SMR<-read.delim("colocs.tsv", header=TRUE, sep="\t", stringsAsFactors = FALSE)

###################

#Reformatting the F2_eqtl results and SMR results & joining them together

str(F2_eqtls)
#gene identifier column = gene_id
#SNP identifier column= variant_id

str(SMR)
#gene identifier column = gene_id
#SNP identifier column= variant_id
max(table(SMR$variant_id))
#[1] 60
#Huh....I guess some of the variants show up more than once?
table(table(SMR$variant_id))
# 20   40   60 
# 5741 89    6 
#But most are just one per trait. Interesting.
#we need to split this by trait

table(SMR$trait)
# epm_boli              epm_distance_traveled epm_percent_time_open_arm       epm_time_immobile_s 
# 5937                      5937                      5937                      5937 
# latency_score_day6        latency_score_day7        lateral_loco_score        of_boli 
# 5937                      5937                      5937                      5937 
# of_distance_traveled of_percent_time_in_center      of_time_immobile_s      pca_index_day6only 
# 5937                      5937                      5937                      5937 
# pca_index_day7only      pca_index_days_6and7        prob_diff_day6        prob_diff_day7 
# 5937                      5937                      5937                      5937 
# rearing_loco_score        response_bias_day6        response_bias_day7          total_loco_score 
# 5937                      5937                      5937                      5937 

#Daniel included analyses for a bunch of behavioral variables that are actually subsets of a particular behavior (e.g. lateral+rearing=locoscore)
#I'm only going to grab the ones that are behaviors that we used in our F2 DE analysis, or analagous behaviors in the F2 juveniles:

SMR_Relevant<-list(SMR[SMR$trait=="total_loco_score",], 
                   SMR[SMR$trait=="epm_time_immobile_s",], 
                   SMR[SMR$trait=="epm_distance_traveled",], 
                   SMR[SMR$trait=="epm_percent_time_open_arm",], 
                   SMR[SMR$trait=="pca_index_days_6and7",],
                   SMR[SMR$trait=="of_distance_traveled",],
                   SMR[SMR$trait=="of_time_immobile_s",],
                   SMR[SMR$trait=="of_percent_time_in_center",]
) 

str(SMR_Relevant)
#List of 8

#are the variant ids in the same order for all of these traits?
sum(SMR_Relevant[[1]]$variant_id==SMR_Relevant[[2]]$variant_id)
#[1] 5937
sum(SMR_Relevant[[1]]$variant_id==SMR_Relevant[[3]]$variant_id)
sum(SMR_Relevant[[1]]$variant_id==SMR_Relevant[[4]]$variant_id)
sum(SMR_Relevant[[1]]$variant_id==SMR_Relevant[[5]]$variant_id)
sum(SMR_Relevant[[1]]$variant_id==SMR_Relevant[[6]]$variant_id)
sum(SMR_Relevant[[1]]$variant_id==SMR_Relevant[[7]]$variant_id)
sum(SMR_Relevant[[1]]$variant_id==SMR_Relevant[[8]]$variant_id)

#Nice - they are all in the same order (all sum to 5937)

#skip
# library(plyr)
# join_all(dfs=SMR_Relevant, by="variant_id", type="left")

colnames(SMR_Relevant[[1]])
#[1] "gene_id"    "gene_name"  "variant_id" "z_eQTL"     "trait"      "z_GWAS"     "T_SMR"      "p_SMR" 

for(i in c(1:8)){
  colnames(SMR_Relevant[[i]])[4:8]<-paste(SMR_Relevant[[i]]$trait[1], colnames(SMR_Relevant[[i]])[4:8], sep="_") 
}

SMR_Relevant_DF<-cbind.data.frame(SMR_Relevant[[1]][,c(-5)], 
                                  SMR_Relevant[[2]][,c(4,6:8)], 
                                  SMR_Relevant[[3]][,c(4,6:8)],
                                  SMR_Relevant[[4]][,c(4,6:8)],
                                  SMR_Relevant[[5]][,c(4,6:8)],
                                  SMR_Relevant[[6]][,c(4,6:8)],
                                  SMR_Relevant[[7]][,c(4,6:8)],
                                  SMR_Relevant[[8]][,c(4,6:8)])

colnames(SMR_Relevant_DF)
# [1] "gene_id"                          "gene_name"                        "variant_id"                       "total_loco_score_z_eQTL"         
# [5] "total_loco_score_z_GWAS"          "total_loco_score_T_SMR"           "total_loco_score_p_SMR"           "epm_time_immobile_s_z_eQTL"      
# [9] "epm_time_immobile_s_z_GWAS"       "epm_time_immobile_s_T_SMR"        "epm_time_immobile_s_p_SMR"        "epm_distance_traveled_z_eQTL"    
# [13] "epm_distance_traveled_z_GWAS"     "epm_distance_traveled_T_SMR"      "epm_distance_traveled_p_SMR"      "epm_percent_time_open_arm_z_eQTL"
# [17] "epm_percent_time_open_arm_z_GWAS" "epm_percent_time_open_arm_T_SMR"  "epm_percent_time_open_arm_p_SMR"  "pca_index_days_6and7_z_eQTL"     
# [21] "pca_index_days_6and7_z_GWAS"      "pca_index_days_6and7_T_SMR"       "pca_index_days_6and7_p_SMR"       "of_distance_traveled_z_eQTL"     
# [25] "of_distance_traveled_z_GWAS"      "of_distance_traveled_T_SMR"       "of_distance_traveled_p_SMR"       "of_time_immobile_s_z_eQTL"       
# [29] "of_time_immobile_s_z_GWAS"        "of_time_immobile_s_T_SMR"         "of_time_immobile_s_p_SMR"         "of_percent_time_in_center_z_eQTL"
# [33] "of_percent_time_in_center_z_GWAS" "of_percent_time_in_center_T_SMR"  "of_percent_time_in_center_p_SMR"    

str(SMR_Relevant_DF)
#'data.frame':	5937 obs. of  35 variables:

#to join them we need to make a combined variant id, gene id column
SMR_Relevant_DF$Gene_Variant_ID<-paste(SMR_Relevant_DF$gene_id, SMR_Relevant_DF$variant_id, sep="_")

str(F2_eqtls)
#'data.frame':	5937 obs. of  18 variables:
F2_eqtls$Gene_Variant_ID<-paste(F2_eqtls$gene_id, F2_eqtls$variant_id, sep="_")

#Just double-checking the numbers
sum(is.na(F2_eqtls$log2_aFC_alt_ref)==FALSE)
#[1] 5937
sum(is.na(F2_eqtls$log2_aFC_bLR_bHR)==FALSE)
#[1] 5703

library(plyr)
F2_eqtls_SMR<-join(F2_eqtls, SMR_Relevant_DF, by="Gene_Variant_ID", type="left", match="all")
str(F2_eqtls_SMR)
# 'data.frame':	5937 obs. of  54 variables:
# $ gene_id                         : chr  "ENSRNOG00000014330" "ENSRNOG00000015239" "ENSRNOG00000016381" "ENSRNOG00000013160" ...
# $ gene_name                       : chr  "Pcmt1" "Ginm1" "Ust" "Sash1" ...
# $ num_var                         : int  1872 2074 3342 3630 4222 4297 4297 3180 3180 2213 ...
# $ variant_id                      : chr  "chr1:1482318" "chr1:2828974" "chr1:3221846" "chr1:2774862" ...
# $ chrom                           : int  1 1 1 1 1 1 1 1 1 1 ...
# $ pos                             : int  1482318 2828974 3221846 2774862 4016846 4366305 3618606 4605341 5629945 6387063 ...
# $ ref                             : chr  "G" "G" "C" "G" ...
# $ alt                             : chr  "T" "A" "G" "A" ...
# $ af                              : num  0.0266 0.4119 0.3053 0.5922 0.6393 ...
# $ tss_distance                    : int  285300 -942993 -594099 71338 -253454 -354794 392905 47869 -976735 938105 ...
# $ pval_nominal                    : num  1.00e-07 9.45e-12 2.53e-04 5.77e-11 2.69e-35 ...
# $ slope                           : num  -0.613 -0.521 0.326 -0.502 0.814 ...
# $ slope_se                        : num  0.1114 0.0726 0.0876 0.0731 0.055 ...
# $ pval_beta                       : num  1.84e-05 2.95e-09 2.76e-02 2.12e-08 7.08e-31 ...
# $ rank                            : int  1 1 1 1 1 1 2 1 2 1 ...
# $ log2_aFC_alt_ref                : num  0.3959 -0.3171 0.0823 -0.1907 0.5983 ...
# $ log2_aFC_bLR_bHR                : num  0.3959 0.3171 0.0823 -0.1907 -0.5983 ...
# $ upregulated_in                  : chr  "bLR" "bLR" "bLR" "bHR" ...
# $ Gene_Variant_ID                 : chr  "ENSRNOG00000014330_chr1:1482318" "ENSRNOG00000015239_chr1:2828974" "ENSRNOG00000016381_chr1:3221846" "ENSRNOG00000013160_chr1:2774862" ...
# $ gene_id                         : chr  "ENSRNOG00000014330" "ENSRNOG00000015239" "ENSRNOG00000016381" "ENSRNOG00000013160" ...
# $ gene_name                       : chr  "Pcmt1" "Ginm1" "Ust" "Sash1" ...
# $ variant_id                      : chr  "chr1:1482318" "chr1:2828974" "chr1:3221846" "chr1:2774862" ...
# $ total_loco_score_z_eQTL         : num  -5.5 -7.18 3.72 -6.87 14.79 ...
# $ total_loco_score_z_GWAS         : num  -0.000944 0.323604 0.318613 -0.471269 0.682428 ...
# $ total_loco_score_T_SMR          : num  8.90e-07 1.05e-01 1.01e-01 2.21e-01 4.65e-01 ...
# $ total_loco_score_p_SMR          : num  0.999 0.746 0.751 0.638 0.495 ...
# $ epm_time_immobile_s_z_eQTL      : num  -5.5 -7.18 3.72 -6.87 14.79 ...
# $ epm_time_immobile_s_z_GWAS      : num  -0.497 -0.888 -0.854 0.841 -0.321 ...
# $ epm_time_immobile_s_T_SMR       : num  0.245 0.776 0.692 0.697 0.103 ...
# $ epm_time_immobile_s_p_SMR       : num  0.621 0.378 0.405 0.404 0.748 ...
# $ epm_distance_traveled_z_eQTL    : num  -5.5 -7.18 3.72 -6.87 14.79 ...
# $ epm_distance_traveled_z_GWAS    : num  0.977 0.709 0.827 -0.61 0.131 ...
# $ epm_distance_traveled_T_SMR     : num  0.9255 0.4974 0.651 0.3687 0.0171 ...
# $ epm_distance_traveled_p_SMR     : num  0.336 0.481 0.42 0.544 0.896 ...
# $ epm_percent_time_open_arm_z_eQTL: num  -5.5 -7.18 3.72 -6.87 14.79 ...
# $ epm_percent_time_open_arm_z_GWAS: num  1.454 0.734 0.132 -0.365 0.688 ...
# $ epm_percent_time_open_arm_T_SMR : num  1.976 0.5332 0.0174 0.1327 0.472 ...
# $ epm_percent_time_open_arm_p_SMR : num  0.16 0.465 0.895 0.716 0.492 ...
# $ pca_index_days_6and7_z_eQTL     : num  -5.5 -7.18 3.72 -6.87 14.79 ...
# $ pca_index_days_6and7_z_GWAS     : num  1.64369 0.2445 -0.00256 -0.42783 0.1521 ...
# $ pca_index_days_6and7_T_SMR      : num  2.48 5.97e-02 6.54e-06 1.82e-01 2.31e-02 ...
# $ pca_index_days_6and7_p_SMR      : num  0.115 0.807 0.998 0.669 0.879 ...
# $ of_distance_traveled_z_eQTL     : num  -5.5 -7.18 3.72 -6.87 14.79 ...
# $ of_distance_traveled_z_GWAS     : num  -1.01 1.096 0.741 -0.751 -0.383 ...
# $ of_distance_traveled_T_SMR      : num  0.987 1.174 0.528 0.558 0.147 ...
# $ of_distance_traveled_p_SMR      : num  0.32 0.279 0.468 0.455 0.702 ...
# $ of_time_immobile_s_z_eQTL       : num  -5.5 -7.18 3.72 -6.87 14.79 ...
# $ of_time_immobile_s_z_GWAS       : num  0.386 -1.301 -1.368 0.92 1.142 ...
# $ of_time_immobile_s_T_SMR        : num  0.148 1.638 1.648 0.832 1.295 ...
# $ of_time_immobile_s_p_SMR        : num  0.7 0.201 0.199 0.362 0.255 ...
# $ of_percent_time_in_center_z_eQTL: num  -5.5 -7.18 3.72 -6.87 14.79 ...
# $ of_percent_time_in_center_z_GWAS: num  1.433 1.172 0.31 -0.937 -1.353 ...
# $ of_percent_time_in_center_T_SMR : num  1.9233 1.337 0.0955 0.8615 1.8146 ...
# $ of_percent_time_in_center_p_SMR : num  0.165 0.248 0.757 0.353 0.178 ...

colnames(F2_eqtls_SMR)
# [1] "gene_id"                          "gene_name"                        "num_var"                          "variant_id"                      
# [5] "chrom"                            "pos"                              "ref"                              "alt"                             
# [9] "af"                               "tss_distance"                     "pval_nominal"                     "slope"                           
# [13] "slope_se"                         "pval_beta"                        "rank"                             "log2_aFC_alt_ref"                
# [17] "log2_aFC_bLR_bHR"                 "upregulated_in"                   "Gene_Variant_ID"                  "gene_id"                         
# [21] "gene_name"                        "variant_id"                       "total_loco_score_z_eQTL"          "total_loco_score_z_GWAS"         
# [25] "total_loco_score_T_SMR"           "total_loco_score_p_SMR"           "epm_time_immobile_s_z_eQTL"       "epm_time_immobile_s_z_GWAS"      
# [29] "epm_time_immobile_s_T_SMR"        "epm_time_immobile_s_p_SMR"        "epm_distance_traveled_z_eQTL"     "epm_distance_traveled_z_GWAS"    
# [33] "epm_distance_traveled_T_SMR"      "epm_distance_traveled_p_SMR"      "epm_percent_time_open_arm_z_eQTL" "epm_percent_time_open_arm_z_GWAS"
# [37] "epm_percent_time_open_arm_T_SMR"  "epm_percent_time_open_arm_p_SMR"  "pca_index_days_6and7_z_eQTL"      "pca_index_days_6and7_z_GWAS"     
# [41] "pca_index_days_6and7_T_SMR"       "pca_index_days_6and7_p_SMR"       "of_distance_traveled_z_eQTL"      "of_distance_traveled_z_GWAS"     
# [45] "of_distance_traveled_T_SMR"       "of_distance_traveled_p_SMR"       "of_time_immobile_s_z_eQTL"        "of_time_immobile_s_z_GWAS"       
# [49] "of_time_immobile_s_T_SMR"         "of_time_immobile_s_p_SMR"         "of_percent_time_in_center_z_eQTL" "of_percent_time_in_center_z_GWAS"
# [53] "of_percent_time_in_center_T_SMR"  "of_percent_time_in_center_p_SMR" 

#There are some redundant columns here (#20-22)
F2_eqtls_SMR<-F2_eqtls_SMR[,-c(20:22)]
colnames(F2_eqtls_SMR)

sum(is.na(F2_eqtls_SMR$log2_aFC_alt_ref)==FALSE)
#[1] 5937
#All eQTLs still present.

length(unique(F2_eqtls_SMR$variant_id))
#[1] 5836  # of variants (instead of gene/variant combinations)

#write.csv(F2_eqtls_SMR, "F2_eqtls_SMR_Joined.csv")
write.csv(F2_eqtls_SMR, "F2_eqtls_SMR_Joined_withJuvenileQTLs.csv")

#save.image("~/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/F2_eQTLs/Workspace_Joining_DEResults_eQTLs_20230926.RData")
save.image("~/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/F2_eQTLs/Workspace_Joining_DEResults_eQTLs_20240119.RData")

#########################