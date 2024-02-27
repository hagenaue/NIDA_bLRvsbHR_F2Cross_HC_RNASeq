############

#I'm going to add in one more piece of information to play with
#Upregulation or Downregulation in other rat models of internalizing behavior
#From the Birt et al. database
#The Garafola & Hen results are excluded because they are denoted with a direction of effect in the database
#This database only includes gene symbol (some of which is old annotation...)

OtherRatModels_Up<-read.csv("OtherRatModels_UpWInternalizing_NoHRLR.csv", header=TRUE, stringsAsFactors = FALSE)
str(OtherRatModels_Up)
# 'data.frame':	1187 obs. of  2 variables:
# $ GeneSymbol_Rat: chr  "Ak1" "Drg1" "Ghdc" "Tmem144" ...
# $ DatasetName   : chr  "Blaveri_2010_FlindersSensitiveLine_Upregulated_HC" "Blaveri_2010_FlindersSensitiveLine_Upregulated_HC" "Blaveri_2010_FlindersSensitiveLine_Upregulated_HC" "Blaveri_2010_FlindersSensitiveLine_Upregulated_HC" ...

max(table(OtherRatModels_Up$GeneSymbol_Rat))
#[1] 3
OtherRatModels_Up_Table<-table(OtherRatModels_Up$GeneSymbol_Rat)
str(OtherRatModels_Up_Table)
# 'table' int [1:1121(1d)] 1 1 1 2 2 1 1 1 2 1 ...
# - attr(*, "dimnames")=List of 1
# ..$ : chr [1:1121] "A2m" "Aars2" "Abat" "Abcb10" ...

OtherRatModels_Up_Table[OtherRatModels_Up_Table>1]
# Abcb10        Abcc9       Abhd13       Adam1a          Ak1         Aqp9        Arl16         Bmp3         Bphl 
# 2            2            2            2            2            2            2            2            2 
# Ccdc77       Chi3l1        Cisd2       Clec9a         Cmc1         Comt          Dcn         Drg1       Endod1 
# 2            2            2            2            2            2            2            2            2 
# Fxyd3        Grem2         Ist1        Larp6         Lcp1 LOC102546433       Lsm14b          Ltk        Magi1 
# 2            2            2            2            2            2            2            2            2 
# Mfap3l       Mfsd10        Mgst1        Mis12        Ophn1        P2rx4       Pik3r1         Pkia         Pnpo 
# 2            2            3            2            2            2            2            2            2 
# Ppfibp2        Ppp4c       Prss35         Ptgr1        Ptprj        Pvrl1        Rab3b          Rdx   RGD1305680 
# 2            2            2            2            2            2            2            2            2 
# Ripk2      S100a10      Sec23ip     Serpinb9      Slc24a3          Sst          St5        Strbp         Syt1 
# 2            2            2            2            2            2            2            2            2 
# Tgfbi        Thrsp         Thy1       Tmco5a      Tmem144        Tradd        Traf6         Tufm       Vps13c 
# 2            2            2            2            3            2            2            2            2 
# Zfp90 
# 2 

table(OtherRatModels_Up_Table)
# OtherRatModels_Up_Table
#   1    2    3 
# 1057   62    2 

#Actually in our dataset:
table(table(OtherRatModels_Up$GeneSymbol_Rat[OtherRatModels_Up$GeneSymbol_Rat%in%F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F0..Rnor6.Ensembl.v.103.]))
# 1   2   3 
# 765  55   2 

OtherRatModels_Down<-read.csv("OtherRatModels_UpWExternalizing_NoHRLR.csv", header=TRUE, stringsAsFactors = FALSE)
str(OtherRatModels_Down)
# 'data.frame':	1532 obs. of  2 variables:
# $ GeneSymbol_Rat: chr  "Acox3" "Apln" "Cav1" "Eif2b1" ...
# $ DatasetName   : chr  "Blaveri_2010_FlindersSensitiveLine_Downregulated_HC" "Blaveri_2010_FlindersSensitiveLine_Downregulated_HC" "Blaveri_2010_FlindersSensitiveLine_Downregulated_HC" "Blaveri_2010_FlindersSensitiveLine_Downregulated_HC"

OtherRatModels_Down_Table<-table(OtherRatModels_Down$GeneSymbol_Rat)
str(OtherRatModels_Down_Table)
# 'table' int [1:1437(1d)] 1 1 1 1 1 1 1 1 1 1 ...
# - attr(*, "dimnames")=List of 1
# ..$ : chr [1:1437] "Aak1" "Aarsd1" "Abca1" "Abca8a" ...

OtherRatModels_Down_Table[OtherRatModels_Down_Table>1]
# Adpgk      Agtpbp1        Akap5         Apln        Birc2         Blnk         C1qa           C3       Calml4 
# 2            2            2            2            2            2            2            2            2 
# Camkk2         Cav2        Ccnd1         Cd48         Cd74        Cep95        Clic2        Clic6       Col6a1 
# 2            2            2            3            2            2            4            2            2 
# Commd6        Creb1         Csad       Cyp4f4       Cyp4v3         Dab2       Deptor         Dpyd          Dut 
# 2            2            2            2            2            2            2            2            2 
# Eif2b1        Ephx2         Etfa         F11r          F2r           F5        Fbln2        Fcrl2        Frem3 
# 2            2            2            2            2            2            2            2            2 
# Fxyd5         Gldc          Gls        Golm1       Hapln1         Ilf3       Inpp5f        Itga7        Kalrn 
# 2            2            2            2            3            2            2            2            2 
# Krt8 LOC100174910 LOC100911253 LOC102552880          Lum         Mcat         Mdm4        Mfap3         Mfrp 
# 2            2            2            2            2            2            2            2            2 
# Mpeg1          Mx2        Nedd4         Nqo2        Otub1         Pak1         Perp      Phf20l1      Pik3c2a 
# 2            2            2            2            2            2            2            2            2 
# Prlr        Psme4         Pwp2         Pxdn   RGD1308772   RGD1309104   RGD1359508       Rnase4        Rsrp1 
# 2            2            2            2            2            2            2            2            2 
# RT1-Da      Rtn4ip1       Sez6l2      Slc31a1      Slc6a20      Sostdc1         Tcn2        Thoc1      Tmem100 
# 2            2            2            2            2            3            2            2            2 
# Tmem176a        Tmem2       Tspan7          Ttr      Tubgcp2        Vdac1       Vwa5b2        Yipf3       Zfp180 
# 2            2            2            2            2            2            2            2            2 

table(OtherRatModels_Down_Table)
# OtherRatModels_Down_Table
#   1    2    3    4 
# 1347   86    3    1 

#Actually in our dataset:
table(table(OtherRatModels_Down$GeneSymbol_Rat[OtherRatModels_Down$GeneSymbol_Rat%in%F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F0..Rnor6.Ensembl.v.103.]))
#   1    2    3 
# 1106   80    3

#Test code:
# i<-1
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME..Rnor6.Ensembl.v.88.[i]
# #[1] "Lrp11"
# OtherRatModels_Up_Table[names(OtherRatModels_Up_Table)=="Lrp11"][[1]]
# OtherRatModels_Down_Table[names(OtherRatModels_Down_Table)=="Lrp11"][[1]]
# #Error in OtherRatModels_Down_Table[names(OtherRatModels_Down_Table) ==  : 
#                                      #subscript out of bounds

F0_Meta_F2_DEResults_w_F2_eqtls_SMR$OtherInternalizingRatModels_Up<-rep(0,nrow(F0_Meta_F2_DEResults_w_F2_eqtls_SMR))
F0_Meta_F2_DEResults_w_F2_eqtls_SMR$OtherInternalizingRatModels_Down<-rep(0,nrow(F0_Meta_F2_DEResults_w_F2_eqtls_SMR))

for(i in c(1:nrow(F0_Meta_F2_DEResults_w_F2_eqtls_SMR))){
  GeneName<-F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F0..Rnor6.Ensembl.v.103.[i]
  if(GeneName%in%names(OtherRatModels_Up_Table)){
    F0_Meta_F2_DEResults_w_F2_eqtls_SMR$OtherInternalizingRatModels_Up[i]<-OtherRatModels_Up_Table[names(OtherRatModels_Up_Table)==GeneName][[1]]
  }else{}
  if(GeneName%in%names(OtherRatModels_Down_Table)){
    F0_Meta_F2_DEResults_w_F2_eqtls_SMR$OtherInternalizingRatModels_Down[i]<-OtherRatModels_Down_Table[names(OtherRatModels_Down_Table)==GeneName][[1]]
  }else{}
}

str(F0_Meta_F2_DEResults_w_F2_eqtls_SMR)
#'data.frame':	14353 obs. of  126 variables:

table(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$OtherInternalizingRatModels_Up)
# 0         1     2     3 
# 13457   825    68     3 

#Intersection of bLR up and bLR/bHR DE genes (candidates from the paper)
table(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$OtherInternalizingRatModels_Up,F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.)
#       0     1
# 0 12419  1038
# 1   717   108
# 2    44    24
# 3     1     2

#Definitely an enrichment

table(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$OtherInternalizingRatModels_Down)
#     0     1     2     3 
# 13085  1172    92     4 

#Intersection of bLR down and bLR/bHR DE genes (candidates from the paper)
table(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$OtherInternalizingRatModels_Down,F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.)
#       0     1
# 0 12096   989
# 1  1014   158
# 2    69    23
# 3     2     2

#Definitely an enrichment

############


#Processing the variables a little more before using them in correlations/PCA

colnames(F0_Meta_F2_DEResults_w_F2_eqtls_SMR)
# [1] "X"                                             "ENSEMBL.ID..Rnor6.Ensembl.v.103."              "GENENAME..Rnor6.Ensembl.v.88."                
# [4] "ENTREZID..Rnor6.Ensembl.v.88."                 "SEQCOORDSYSTEM..Rnor6.Ensembl.v.88."           "SEQNAME..Rnor6.Ensembl.v.88."                 
# [7] "SEQSTRAND..Rnor6.Ensembl.v.88."                "GENESEQSTART..Rnor6.Ensembl.v.88."             "GENESEQEND..Rnor6.Ensembl.v.88."              
# [10] "GENEBIOTYPE..Rnor6.Ensembl.v.88."              "GENENAME_F0..Rnor6.Ensembl.v.103."             "F0_Average.Expression..Log2.CPM."             
# [13] "F0_Log2FC_bLRvsbHR"                            "F0_Tstat_bLRvsbHR"                             "F0_Pvalue_bLRvsbHR"                           
# [16] "F0_FDR_bLRvsbHR"                               "MetaAnalysis_estimatedD_bLRvsbHR"              "MetaAnalysis_SE_bLRvsbHR"                     
# [19] "MetaAnalysis_Pvalue_bLRvsbHR"                  "MetaAnalysis_FDR_bLRvsbHR"                     "GENENAME_F2..Rnor6.Ensembl.v.103."            
# [22] "Included.in.F2.Candidate.Gene.Analysis."       "F2_Log2FC_LocoScore"                           "F2_Tstat_LocoScore"                           
# [25] "F2_Pvalue_LocoScore"                           "F2_FDR_LocoScore"                              "F2_Log2FC_EPM_DistanceTraveled"               
# [28] "F2_Tstat_EPM_DistanceTraveled"                 "F2_Pvalue_EPM_DistanceTraveled"                "F2_FDR_EPM_DistanceTraveled"                  
# [31] "F2_Log2FC_EPM_Time_Immobile"                   "F2_Tstat_EPM_Time_Immobile"                    "F2_Pvalue_EPM_Time_Immobile"                  
# [34] "F2_FDR_EPM_Time_Immobile"                      "F2_Log2FC_EPM_Percent_Time_Open_Arms"          "F2_Tstat_EPM_Percent_Time_Open_Arms"          
# [37] "F2_Pvalue_EPM_Percent_Time_Open_Arms"          "F2_FDR_EPM_Percent_Time_Open_Arms"             "F2_Log2FC_PCA_Index"                          
# [40] "F2_Tstat_PCA_Index"                            "F2_Pvalue_PCA_Index"                           "F2_FDR_PCA_Index"                             
# [43] "gene_name"                                     "num_var"                                       "variant_id"                                   
# [46] "chrom"                                         "pos"                                           "ref"                                          
# [49] "alt"                                           "af"                                            "tss_distance"                                 
# [52] "pval_nominal"                                  "slope"                                         "slope_se"                                     
# [55] "pval_beta"                                     "rank"                                          "log2_aFC_alt_ref"                             
# [58] "log2_aFC_bLR_bHR"                              "upregulated_in"                                "Gene_Variant_ID"                              
# [61] "total_loco_score_z_eQTL"                       "total_loco_score_z_GWAS"                       "total_loco_score_T_SMR"                       
# [64] "total_loco_score_p_SMR"                        "epm_time_immobile_s_z_eQTL"                    "epm_time_immobile_s_z_GWAS"                   
# [67] "epm_time_immobile_s_T_SMR"                     "epm_time_immobile_s_p_SMR"                     "epm_distance_traveled_z_eQTL"                 
# [70] "epm_distance_traveled_z_GWAS"                  "epm_distance_traveled_T_SMR"                   "epm_distance_traveled_p_SMR"                  
# [73] "epm_percent_time_open_arm_z_eQTL"              "epm_percent_time_open_arm_z_GWAS"              "epm_percent_time_open_arm_T_SMR"              
# [76] "epm_percent_time_open_arm_p_SMR"               "pca_index_days_6and7_z_eQTL"                   "pca_index_days_6and7_z_GWAS"                  
# [79] "pca_index_days_6and7_T_SMR"                    "pca_index_days_6and7_p_SMR"                    "of_distance_traveled_z_eQTL"                  
# [82] "of_distance_traveled_z_GWAS"                   "of_distance_traveled_T_SMR"                    "of_distance_traveled_p_SMR"                   
# [85] "of_time_immobile_s_z_eQTL"                     "of_time_immobile_s_z_GWAS"                     "of_time_immobile_s_T_SMR"                     
# [88] "of_time_immobile_s_p_SMR"                      "of_percent_time_in_center_z_eQTL"              "of_percent_time_in_center_z_GWAS"             
# [91] "of_percent_time_in_center_T_SMR"               "of_percent_time_in_center_p_SMR"               "CHROM"                                        
# [94] "POS"                                           "5739-JL-0021"                                  "5739-JL-0022"                                 
# [97] "5739-JL-0023"                                  "5739-JL-0024"                                  "5739-JL-0025"                                 
# [100] "5739-JL-0026"                                  "5739-JL-0027"                                  "5739-JL-0028"                                 
# [103] "5739-JL-0029"                                  "5739-JL-0030"                                  "5739-JL-0031"                                 
# [106] "5739-JL-0032"                                  "5739-JL-0033"                                  "5739-JL-0034"                                 
# [109] "5739-JL-0035"                                  "5739-JL-0036"                                  "5739-JL-0037"                                 
# [112] "5739-JL-0038"                                  "5739-JL-0039"                                  "5739-JL-0040"                                 
# [115] "Hs_FALSE"                                      "Hs_TRUE"                                       "Ht"                                           
# [118] "n_FALSE"                                       "n_TRUE"                                        "Gst"                                          
# [121] "Htmax"                                         "Gstmax"                                        "Gprimest"                                     
# [124] "total_loco_score_T_SMR_directionalv2"          "epm_time_immobile_s_T_SMR_directionalv2"       "epm_distance_traveled_T_SMR_directionalv2"    
# [127] "epm_percent_time_open_arm_T_SMR_directionalv2" "pca_index_days_6and7_T_SMR_directionalv2"      "of_distance_traveled_T_SMR_directionalv2"     
# [130] "of_time_immobile_s_T_SMR_directionalv2"        "of_percent_time_in_center_T_SMR_directionalv2" "OtherInternalizingRatModels_Up"               
# [133] "OtherInternalizingRatModels_Down"  

#Useful:
#To stick with relatively comparable values, we should try to grab the t-stats

#bHR/bLR genetic segregation:

#[123] "Gprimest" 
#No t-statistic available 

#bHR/bLR differential expression:

#[14] "F0_Tstat_bLRvsbHR"

#This could be used to calculate a Tstat:
# [17] "MetaAnalysis_estimatedD_bLRvsbHR"       
# [18] "MetaAnalysis_SE_bLRvsbHR"
F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_Tstat<-F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR/F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_SE_bLRvsbHR
hist(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_Tstat)
#Looks good: range-6 to 7
#So that would now be column:
#[134]

#Predicted bHR/bLR DE: This could be used to calculate a Tstat for aFC
# [53] "slope"                                  
# [54] "slope_se" 
# [58] "log2_aFC_bLR_bHR" (to give it correct directionality)
# [59] "upregulated_in" (maybe easier way to assign directionality?)
boxplot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$upregulated_in)
#yep, looks like "upregulated in" is the way to go
table(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$upregulated_in)
# bHR   bLR  equal 
# 2626  2897   224
#although some eQTLs are the same alleles in bHR/bLR - Daniel cut them out from the aFC predictions.
#not a negatable number
#why don't the NAs show up in the table?
sum(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$upregulated_in=="equal", na.rm=TRUE)
sum(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$upregulated_in=="bHR", na.rm=TRUE)
sum((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$upregulated_in=="bHR")&is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$upregulated_in)==FALSE)
#[1] 2626
F0_Meta_F2_DEResults_w_F2_eqtls_SMR$eQTL_Tstat<-abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope/F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope_se)
hist(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$eQTL_Tstat)
F0_Meta_F2_DEResults_w_F2_eqtls_SMR$eQTL_Tstat[(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$upregulated_in=="bHR")&is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$upregulated_in)==FALSE]<-(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$eQTL_Tstat[(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$upregulated_in=="bHR")&is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$upregulated_in)==FALSE])*(-1)
F0_Meta_F2_DEResults_w_F2_eqtls_SMR$eQTL_Tstat[(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$upregulated_in=="equal")&is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$upregulated_in)==FALSE]<-(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$eQTL_Tstat[(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$upregulated_in=="equal")&is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$upregulated_in)==FALSE])*(0)
hist(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$eQTL_Tstat)
#now the t-stats are directional (based on bHR/bLR genetics) and the "equal" variants are coded as 0
#So that would now be column:
#[135]

#F2 differential expression:
# [24] "F2_Tstat_LocoScore" 
# [28] "F2_Tstat_EPM_DistanceTraveled" 
# [32] "F2_Tstat_EPM_Time_Immobile"   
# [36] "F2_Tstat_EPM_Percent_Time_Open_Arms" 
# [40] "F2_Tstat_PCA_Index"

#F2 eQTL/QTL co-localization: SMR (non-directional)
# [63]  "total_loco_score_T_SMR"                       
# [67] "epm_time_immobile_s_T_SMR"                                  
# [71] "epm_distance_traveled_T_SMR"                                   
# [75] "epm_percent_time_open_arm_T_SMR"              
# [79] "pca_index_days_6and7_T_SMR"                 
# [83] "of_distance_traveled_T_SMR"                                   
# [87] "of_time_immobile_s_T_SMR"                     
# [91] "of_percent_time_in_center_T_SMR"                     

#F2 eQTL/QTL co-localization: SMR (directional)
# [124] "total_loco_score_T_SMR_directionalv2"          
# [125] "epm_time_immobile_s_T_SMR_directionalv2"       
# [126] "epm_distance_traveled_T_SMR_directionalv2"    
# [127] "epm_percent_time_open_arm_T_SMR_directionalv2" 
# [128] "pca_index_days_6and7_T_SMR_directionalv2"      
# [129] "of_distance_traveled_T_SMR_directionalv2"     
# [130] "of_time_immobile_s_T_SMR_directionalv2"       
# [131] "of_percent_time_in_center_T_SMR_directionalv2"

#Other models: (not t-stats)
# [132] "OtherInternalizingRatModels_Up"         
# [133] "OtherInternalizingRatModels_Down" 
#These variables are mostly 0's - let's combine this into a single variable:
F0_Meta_F2_DEResults_w_F2_eqtls_SMR$OtherInternalizingRatModels<-F0_Meta_F2_DEResults_w_F2_eqtls_SMR$OtherInternalizingRatModels_Up-F0_Meta_F2_DEResults_w_F2_eqtls_SMR$OtherInternalizingRatModels_Down

table(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$OtherInternalizingRatModels)
# -3    -2    -1     0     1     2     3 
# 4    82  1074 12403   731    56     3 
#It's still mostly 0's, but a little bit better.
#that would be column [136]

#write.csv(F0_Meta_F2_DEResults_w_F2_eqtls_SMR, "F0_Meta_F2_DEResults_w_F2_eqtls_SMR.csv")
write.csv(F0_Meta_F2_DEResults_w_F2_eqtls_SMR, "F0_Meta_F2_DEResults_w_F2_eqtls_SMR_wJuveniles.csv")
