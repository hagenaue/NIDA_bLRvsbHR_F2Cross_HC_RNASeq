#Just running some simple code joining our different U01 results for comparison with our previous meta-analysis:
#This is an updated version with DE results for PavCA, Time Immobile, and Distance Traveled (outlier removed)
#Messy: Some of the old code and output is mixed with the new code
#Megan Hagenauer
#2022-05-09

#Most recent workspace (maybe?):
load("~/Documents/Microarray Gen/HRLR/NIDA_U01/CompareF0F2Meta_Workspace_20220405_wGenetics.RData")
#Nope, that doesn't seem to be it, but I can't find a more recent one, so I guess I'll just read everything back in???

#Main Directory:
setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01")

#Grabbing F0 data joined with Late Gene RNA-Seq Meta-Analysis Results:
setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/f0data")

#F0 Limma Results come from TMM normalized data, with Ebayes corrected output.
#Model: Lineage + Sex + %rRNA + %Intergenic
#Model justification: %rRNA correlates strongly with PC1, %Intergenic correlates strongly with PC2 (and a little bit with PC1) but isn't also (significantly) correlated with Lineage (unlike %UTR or %Coding)
#Models using a Lineage*Sex interaction term didn't have any significant interaction effects... probably because our sample size is too small.

F0_Limma_Results_w_LateGeneMeta<-read.csv("F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen.csv", header=TRUE, stringsAsFactors = FALSE)

str(F0_Limma_Results_w_LateGeneMeta)
# 'data.frame':	13786 obs. of  36 variables:
#   $ X.1                            : int  1 2 3 4 5 6 7 8 9 10 ...
# $ ENSEMBL                        : chr  "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# $ AveExpr                        : num  6.55 9.32 3.84 4.34 4.3 ...
# $ Coef..Intercept.               : num  5.03 11.45 5.54 1.42 5.49 ...
# $ Coef.Lineage_AsFactorbLR       : num  0.12168 -0.05269 0.00477 0.09307 -0.03847 ...
# $ Coef.Sex_FactorFemale          : num  0.1659 -0.083 0.06 0.2199 -0.0687 ...
# $ Coef.RibosomePerc              : num  -16.5 43.9 38.3 32 26.5 ...
# $ Coef.Percent_Intergenic        : num  17.2 -30 -25.6 19.7 -17 ...
# $ t..Intercept.                  : num  9.22 24.12 11.84 2.86 14.75 ...
# $ t.Lineage_AsFactorbLR          : num  0.9776 -0.4866 0.0452 0.8088 -0.4557 ...
# $ t.Sex_FactorFemale             : num  1.393 -0.799 0.596 1.991 -0.847 ...
# $ t.RibosomePerc                 : num  -2.05 6.27 5.67 4.35 4.84 ...
# $ t.Percent_Intergenic           : num  3.44 -6.96 -5.9 4.41 -5 ...
# $ P.value..Intercept.            : num  2.15e-09 1.86e-18 1.44e-11 8.52e-03 1.33e-13 ...
# $ P.value.Lineage_AsFactorbLR    : num  0.338 0.631 0.964 0.427 0.653 ...
# $ P.value.Sex_FactorFemale       : num  0.1763 0.4322 0.5568 0.0578 0.4051 ...
# $ P.value.RibosomePerc           : num  5.08e-02 1.70e-06 7.48e-06 2.14e-04 6.12e-05 ...
# $ P.value.Percent_Intergenic     : num  2.10e-03 3.17e-07 4.15e-06 1.82e-04 4.02e-05 ...
# $ P.value.adj..Intercept.        : num  5.47e-09 5.39e-17 4.73e-11 1.15e-02 6.04e-13 ...
# $ P.value.adj.Lineage_AsFactorbLR: num  0.825 0.931 0.995 0.867 0.933 ...
# $ P.value.adj.Sex_FactorFemale   : num  0.588 0.786 0.849 0.429 0.767 ...
# $ P.value.adj.RibosomePerc       : num  6.82e-02 6.19e-06 2.23e-05 4.51e-04 1.46e-04 ...
# $ P.value.adj.Percent_Intergenic : num  4.69e-03 8.64e-06 3.54e-05 5.96e-04 1.79e-04 ...
# $ F                              : num  2553 6830 1243 1374 2385 ...
# $ F.p.value                      : num  3.73e-32 2.50e-37 2.23e-28 6.63e-29 8.49e-32 ...
# $ SYMBOL                         : chr  "Lrp11" "Pcmt1" "Nup43" "Lats1" ...
# $ X                              : int  6469 8404 8035 5999 5648 4514 9032 12112 13418 10708 ...
# $ estimate                       : num  0.1247 0.8168 0.5024 0.1994 -0.0682 ...
# $ SE                             : num  0.431 0.446 0.475 0.444 0.433 ...
# $ pval                           : num  0.7725 0.0671 0.2901 0.6535 0.8751 ...
# $ CI_lb                          : num  0.97 1.691 1.433 1.07 0.781 ...
# $ CI_ub                          : num  -0.7203 -0.0576 -0.4284 -0.6712 -0.9177 ...
# $ datasets                       : int  2 2 2 2 2 2 2 2 2 2 ...
# $ rawp                           : num  0.7725 0.0671 0.2901 0.6535 0.8751 ...
# $ BH                             : num  0.911 0.316 0.617 0.854 0.951 ...
# $ BY                             : num  1 1 1 1 1 ...

#Note - all of the columns following "SYMBOL" refer to the  bHR/bLR late generation RNA-Seq meta-analysis output from the Birt et al. Biological Psychiatry paper.

#Grabbing the F2 results - this .csv file contains the results from 4 models, since we were waffling about how to best control for noise.  I think M8 and M9 are probably the best. Model descriptions:

# M6_CovJustSexIntergenic: Locomotor Activity + Sex + % Intergenic. 
# M6 Justification: The most basic model. %Intergenic strongly correlates with PC1, which represents >60% of the variation in the dataset. This model most likely doesn't control strongly enough for sources of technical noise, I just wanted it as a comparison point.

# M9_CovSexIntergenicRrnaTechnical: Locomotor Activity + Sex + % Intergenic + % rRNA + Sequencing Batch + Dissector + STGTExperience
# M9 Justification: Controls for PC1 and PC2 as much as possible using %Intergenic and % rRNA, as well as three technical factors that we have strong reason to believe are important and that correlate with top PCs.
# This is probably a pretty good model, and strongly parallel to what was used for F0. May still have lingering noise that can get cleaned up.

# M8_CovSexRNAMetricsTechnical: Locomotor Activity + Sex + % Intergenic + % rRNA + % Intronic+ RNAConc + Sequencing Batch + Dissector + STGTExperience
# M8 Justification: Similar to M9, but includes two more RNA metrics that strongly correlated with top PCs (PC4-PC8): % Intronic+ RNAConc
# This is probably a pretty good model, but less parallel to what was used for F0.

#:M7_CovSexRNAMetricsUTRConfoundTechnical: Locomotor Activity + Sex + % Intergenic + % rRNA + % Intronic+ RNAConc + Sequencing Batch + Dissector + STGTExperience + % UTR
# M8 Justification: Similar to M8, but includes % UTR, which strongly correlates with PC5 in the F2 dataset and with PC2 in the F0 dataset, but was found to also correlate with HR/LR lineage in the F0 dataset. 
# This model makes me nervous, because it is possible that somehow the % UTR vs. HR/LR relationship is a genuine biological signal.

#None of these models produce any locomotor activity results that survive FDR correction, but peaking at the top results they do overlap previous HR/LR findings (so noisy, but not invalid)

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/reworktogethertuesday142022/reworktogethertuesday142022")

F2_Limma_Results_from4Models<-read.csv("F2_LocomotorActivityResults_FourModels_ForR.csv", header=TRUE, stringsAsFactors = FALSE)

str(F2_Limma_Results_from4Models)

colnames(F2_Limma_Results_from4Models )
#'data.frame':	14056 obs. of  158 variables:

# [1] "Symbol"                                                                 
# [2] "ENSEMBL"                                                                
# [3] "M6_CovJustSexIntergenic_AveExpr"                                        
# [4] "M6_CovJustSexIntergenic_Coef..Intercept."                               
# [5] "M6_CovJustSexIntergenic_Coef.Total_LocoScore"                           
# [6] "M6_CovJustSexIntergenic_Coef.Sex_AsFactorFemale"                        
# [7] "M6_CovJustSexIntergenic_Coef.Percent.Intergenic"                        
# [8] "M6_CovJustSexIntergenic_t..Intercept."                                  
# [9] "M6_CovJustSexIntergenic_t.Total_LocoScore"                              
# [10] "M6_CovJustSexIntergenic_t.Sex_AsFactorFemale"                           
# [11] "M6_CovJustSexIntergenic_t.Percent.Intergenic"                           
# [12] "M6_CovJustSexIntergenic_P.value..Intercept."                            
# [13] "M6_CovJustSexIntergenic_P.value.Total_LocoScore"                        
# [14] "M6_CovJustSexIntergenic_P.value.Sex_AsFactorFemale"                     
# [15] "M6_CovJustSexIntergenic_P.value.Percent.Intergenic"                     
# [16] "M6_CovJustSexIntergenic_P.value.adj..Intercept."                        
# [17] "M6_CovJustSexIntergenic_P.value.adj.Total_LocoScore"                    
# [18] "M6_CovJustSexIntergenic_P.value.adj.Sex_AsFactorFemale"                 
# [19] "M6_CovJustSexIntergenic_P.value.adj.Percent.Intergenic"                 
# [20] "M6_CovJustSexIntergenic_F"                                              
# [21] "M6_CovJustSexIntergenic_F.p.value"                                      
# [22] "M9_CovSexIntergenicRrnaTechnical_AveExpr"                               
# [23] "M9_CovSexIntergenicRrnaTechnical_Coef..Intercept."                      
# [24] "M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore"                  
# [25] "M9_CovSexIntergenicRrnaTechnical_Coef.Percent.Intergenic"               
# [26] "M9_CovSexIntergenicRrnaTechnical_Coef.RibosomePerc"                     
# [27] "M9_CovSexIntergenicRrnaTechnical_Coef.SequencingBatch2"                 
# [28] "M9_CovSexIntergenicRrnaTechnical_Coef.SequencingBatch3"                 
# [29] "M9_CovSexIntergenicRrnaTechnical_Coef.SexMale"                          
# [30] "M9_CovSexIntergenicRrnaTechnical_Coef.Dissector_AsNumeric"              
# [31] "M9_CovSexIntergenicRrnaTechnical_Coef.STGT_experienceTRUE"              
# [32] "M9_CovSexIntergenicRrnaTechnical_t..Intercept."                         
# [33] "M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore"                     
# [34] "M9_CovSexIntergenicRrnaTechnical_t.Percent.Intergenic"                  
# [35] "M9_CovSexIntergenicRrnaTechnical_t.RibosomePerc"                        
# [36] "M9_CovSexIntergenicRrnaTechnical_t.SequencingBatch2"                    
# [37] "M9_CovSexIntergenicRrnaTechnical_t.SequencingBatch3"                    
# [38] "M9_CovSexIntergenicRrnaTechnical_t.SexMale"                             
# [39] "M9_CovSexIntergenicRrnaTechnical_t.Dissector_AsNumeric"                 
# [40] "M9_CovSexIntergenicRrnaTechnical_t.STGT_experienceTRUE"                 
# [41] "M9_CovSexIntergenicRrnaTechnical_P.value..Intercept."                   
# [42] "M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore"               
# [43] "M9_CovSexIntergenicRrnaTechnical_P.value.Percent.Intergenic"            
# [44] "M9_CovSexIntergenicRrnaTechnical_P.value.RibosomePerc"                  
# [45] "M9_CovSexIntergenicRrnaTechnical_P.value.SequencingBatch2"              
# [46] "M9_CovSexIntergenicRrnaTechnical_P.value.SequencingBatch3"              
# [47] "M9_CovSexIntergenicRrnaTechnical_P.value.SexMale"                       
# [48] "M9_CovSexIntergenicRrnaTechnical_P.value.Dissector_AsNumeric"           
# [49] "M9_CovSexIntergenicRrnaTechnical_P.value.STGT_experienceTRUE"           
# [50] "M9_CovSexIntergenicRrnaTechnical_P.value.adj..Intercept."               
# [51] "M9_CovSexIntergenicRrnaTechnical_P.value.adj.Total_LocoScore"           
# [52] "M9_CovSexIntergenicRrnaTechnical_P.value.adj.Percent.Intergenic"        
# [53] "M9_CovSexIntergenicRrnaTechnical_P.value.adj.RibosomePerc"              
# [54] "M9_CovSexIntergenicRrnaTechnical_P.value.adj.SequencingBatch2"          
# [55] "M9_CovSexIntergenicRrnaTechnical_P.value.adj.SequencingBatch3"          
# [56] "M9_CovSexIntergenicRrnaTechnical_P.value.adj.SexMale"                   
# [57] "M9_CovSexIntergenicRrnaTechnical_P.value.adj.Dissector_AsNumeric"       
# [58] "M9_CovSexIntergenicRrnaTechnical_P.value.adj.STGT_experienceTRUE"       
# [59] "M9_CovSexIntergenicRrnaTechnical_F"                                     
# [60] "M9_CovSexIntergenicRrnaTechnical_F.p.value"                             
# [61] "M8_CovSexRNAMetricsTechnical_AveExpr"                                   
# [62] "M8_CovSexRNAMetricsTechnical_Coef..Intercept."                          
# [63] "M8_CovSexRNAMetricsTechnical_Coef.Total_LocoScore"                      
# [64] "M8_CovSexRNAMetricsTechnical_Coef.Percent.Intergenic"                   
# [65] "M8_CovSexRNAMetricsTechnical_Coef.RibosomePerc"                         
# [66] "M8_CovSexRNAMetricsTechnical_Coef.SequencingBatch2"                     
# [67] "M8_CovSexRNAMetricsTechnical_Coef.SequencingBatch3"                     
# [68] "M8_CovSexRNAMetricsTechnical_Coef.Percent.Intronic"                     
# [69] "M8_CovSexRNAMetricsTechnical_Coef.SexMale"                              
# [70] "M8_CovSexRNAMetricsTechnical_Coef.RNAconc"                              
# [71] "M8_CovSexRNAMetricsTechnical_Coef.Dissector_AsNumeric"                  
# [72] "M8_CovSexRNAMetricsTechnical_Coef.STGT_experienceTRUE"                  
# [73] "M8_CovSexRNAMetricsTechnical_t..Intercept."                             
# [74] "M8_CovSexRNAMetricsTechnical_t.Total_LocoScore"                         
# [75] "M8_CovSexRNAMetricsTechnical_t.Percent.Intergenic"                      
# [76] "M8_CovSexRNAMetricsTechnical_t.RibosomePerc"                            
# [77] "M8_CovSexRNAMetricsTechnical_t.SequencingBatch2"                        
# [78] "M8_CovSexRNAMetricsTechnical_t.SequencingBatch3"                        
# [79] "M8_CovSexRNAMetricsTechnical_t.Percent.Intronic"                        
# [80] "M8_CovSexRNAMetricsTechnical_t.SexMale"                                 
# [81] "M8_CovSexRNAMetricsTechnical_t.RNAconc"                                 
# [82] "M8_CovSexRNAMetricsTechnical_t.Dissector_AsNumeric"                     
# [83] "M8_CovSexRNAMetricsTechnical_t.STGT_experienceTRUE"                     
# [84] "M8_CovSexRNAMetricsTechnical_P.value..Intercept."                       
# [85] "M8_CovSexRNAMetricsTechnical_P.value.Total_LocoScore"                   
# [86] "M8_CovSexRNAMetricsTechnical_P.value.Percent.Intergenic"                
# [87] "M8_CovSexRNAMetricsTechnical_P.value.RibosomePerc"                      
# [88] "M8_CovSexRNAMetricsTechnical_P.value.SequencingBatch2"                  
# [89] "M8_CovSexRNAMetricsTechnical_P.value.SequencingBatch3"                  
# [90] "M8_CovSexRNAMetricsTechnical_P.value.Percent.Intronic"                  
# [91] "M8_CovSexRNAMetricsTechnical_P.value.SexMale"                           
# [92] "M8_CovSexRNAMetricsTechnical_P.value.RNAconc"                           
# [93] "M8_CovSexRNAMetricsTechnical_P.value.Dissector_AsNumeric"               
# [94] "M8_CovSexRNAMetricsTechnical_P.value.STGT_experienceTRUE"               
# [95] "M8_CovSexRNAMetricsTechnical_P.value.adj..Intercept."                   
# [96] "M8_CovSexRNAMetricsTechnical_P.value.adj.Total_LocoScore"               
# [97] "M8_CovSexRNAMetricsTechnical_P.value.adj.Percent.Intergenic"            
# [98] "M8_CovSexRNAMetricsTechnical_P.value.adj.RibosomePerc"                  
# [99] "M8_CovSexRNAMetricsTechnical_P.value.adj.SequencingBatch2"              
# [100] "M8_CovSexRNAMetricsTechnical_P.value.adj.SequencingBatch3"              
# [101] "M8_CovSexRNAMetricsTechnical_P.value.adj.Percent.Intronic"              
# [102] "M8_CovSexRNAMetricsTechnical_P.value.adj.SexMale"                       
# [103] "M8_CovSexRNAMetricsTechnical_P.value.adj.RNAconc"                       
# [104] "M8_CovSexRNAMetricsTechnical_P.value.adj.Dissector_AsNumeric"           
# [105] "M8_CovSexRNAMetricsTechnical_P.value.adj.STGT_experienceTRUE"           
# [106] "M8_CovSexRNAMetricsTechnical_F"                                         
# [107] "M8_CovSexRNAMetricsTechnical_F.p.value"                                 
# [108] "M7_CovSexRNAMetricsUTRConfoundTechnical_AveExpr"                        
# [109] "M7_CovSexRNAMetricsUTRConfoundTechnical_Coef..Intercept."               
# [110] "M7_CovSexRNAMetricsUTRConfoundTechnical_Coef.Total_LocoScore"           
# [111] "M7_CovSexRNAMetricsUTRConfoundTechnical_Coef.Percent.Intergenic"        
# [112] "M7_CovSexRNAMetricsUTRConfoundTechnical_Coef.RibosomePerc"              
# [113] "M7_CovSexRNAMetricsUTRConfoundTechnical_Coef.SequencingBatch2"          
# [114] "M7_CovSexRNAMetricsUTRConfoundTechnical_Coef.SequencingBatch3"          
# [115] "M7_CovSexRNAMetricsUTRConfoundTechnical_Coef.Percent.Intronic"          
# [116] "M7_CovSexRNAMetricsUTRConfoundTechnical_Coef.SexMale"                   
# [117] "M7_CovSexRNAMetricsUTRConfoundTechnical_Coef.Percent.UTR"               
# [118] "M7_CovSexRNAMetricsUTRConfoundTechnical_Coef.RNAconc"                   
# [119] "M7_CovSexRNAMetricsUTRConfoundTechnical_Coef.Dissector_AsNumeric"       
# [120] "M7_CovSexRNAMetricsUTRConfoundTechnical_Coef.STGT_experienceTRUE"       
# [121] "M7_CovSexRNAMetricsUTRConfoundTechnical_t..Intercept."                  
# [122] "M7_CovSexRNAMetricsUTRConfoundTechnical_t.Total_LocoScore"              
# [123] "M7_CovSexRNAMetricsUTRConfoundTechnical_t.Percent.Intergenic"           
# [124] "M7_CovSexRNAMetricsUTRConfoundTechnical_t.RibosomePerc"                 
# [125] "M7_CovSexRNAMetricsUTRConfoundTechnical_t.SequencingBatch2"             
# [126] "M7_CovSexRNAMetricsUTRConfoundTechnical_t.SequencingBatch3"             
# [127] "M7_CovSexRNAMetricsUTRConfoundTechnical_t.Percent.Intronic"             
# [128] "M7_CovSexRNAMetricsUTRConfoundTechnical_t.SexMale"                      
# [129] "M7_CovSexRNAMetricsUTRConfoundTechnical_t.Percent.UTR"                  
# [130] "M7_CovSexRNAMetricsUTRConfoundTechnical_t.RNAconc"                      
# [131] "M7_CovSexRNAMetricsUTRConfoundTechnical_t.Dissector_AsNumeric"          
# [132] "M7_CovSexRNAMetricsUTRConfoundTechnical_t.STGT_experienceTRUE"          
# [133] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value..Intercept."            
# [134] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.Total_LocoScore"        
# [135] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.Percent.Intergenic"     
# [136] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.RibosomePerc"           
# [137] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.SequencingBatch2"       
# [138] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.SequencingBatch3"       
# [139] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.Percent.Intronic"       
# [140] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.SexMale"                
# [141] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.Percent.UTR"            
# [142] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.RNAconc"                
# [143] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.Dissector_AsNumeric"    
# [144] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.STGT_experienceTRUE"    
# [145] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.adj..Intercept."        
# [146] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.adj.Total_LocoScore"    
# [147] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.adj.Percent.Intergenic" 
# [148] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.adj.RibosomePerc"       
# [149] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.adj.SequencingBatch2"   
# [150] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.adj.SequencingBatch3"   
# [151] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.adj.Percent.Intronic"   
# [152] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.adj.SexMale"            
# [153] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.adj.Percent.UTR"        
# [154] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.adj.RNAconc"            
# [155] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.adj.Dissector_AsNumeric"
# [156] "M7_CovSexRNAMetricsUTRConfoundTechnical_P.value.adj.STGT_experienceTRUE"
# [157] "M7_CovSexRNAMetricsUTRConfoundTechnical_F"                              
# [158] "M7_CovSexRNAMetricsUTRConfoundTechnical_F.p.value"   

#First, a quick check to see how similar the results from all of these models are:

#Correlation matrix of Locomotor score betas from the 4 models:

cor(as.matrix( cbind(F2_Limma_Results_from4Models$M6_CovJustSexIntergenic_Coef.Total_LocoScore, F2_Limma_Results_from4Models$M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore, F2_Limma_Results_from4Models$M8_CovSexRNAMetricsTechnical_Coef.Total_LocoScore, F2_Limma_Results_from4Models$M7_CovSexRNAMetricsUTRConfoundTechnical_Coef.Total_LocoScore)))

#          M6        M9        M8        M7    
# [,1]      [,2]      [,3]      [,4]
# [1,] 1.0000000 0.8885221 0.8154752 0.8112624
# [2,] 0.8885221 1.0000000 0.9580079 0.9498443
# [3,] 0.8154752 0.9580079 1.0000000 0.9784341
# [4,] 0.8112624 0.9498443 0.9784341 1.0000000

#So the results from the two models that I suspect are best (M9 and M8) are highly correlated (R=0.95), as are the results from M7.  The simplest model that doesn't control for suspected technical sources of noise is less correlated with the others.

#Here's the t-stat version:
cor(as.matrix( cbind(F2_Limma_Results_from4Models$M6_CovJustSexIntergenic_t.Total_LocoScore, F2_Limma_Results_from4Models$M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore, F2_Limma_Results_from4Models$M8_CovSexRNAMetricsTechnical_t.Total_LocoScore, F2_Limma_Results_from4Models$M7_CovSexRNAMetricsUTRConfoundTechnical_t.Total_LocoScore)))

#          M6        M9        M8        M7  
#         [,1]      [,2]      [,3]      [,4]
# [1,] 1.0000000 0.8396662 0.7583187 0.7807688
# [2,] 0.8396662 1.0000000 0.9235230 0.9171692
# [3,] 0.7583187 0.9235230 1.0000000 0.9694673
# [4,] 0.7807688 0.9171692 0.9694673 1.0000000

#Again, the results are pretty strongly correlated amongst all models that actually contain technical covariates


setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/reworktogethertuesday142022/reworktogethertuesday142022")

#Old version:

# F2_Limma_Results_More_20220121_Summarized<-read.csv("F2_Limma_Results_More_20220121_Summarized.csv", header=TRUE, stringsAsFactors = FALSE)
# str(F2_Limma_Results_More_20220121_Summarized)

F2_Limma_Results_More_20220121_Summarized<-read.csv("F2_Limma_Results_More_20220121_Summarized_NewVarAdded.csv", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results_More_20220121_Summarized)

# 'data.frame':	14056 obs. of  33 variables:
#   $ ENSEMBL                                                : chr  "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# $ Coef.Total_LocoScore_SimplerModel_wDissectionDay       : num  1.23e-05 -1.48e-06 -3.50e-05 -1.92e-05 -4.64e-06 -6.12e-05 1.05e-05 1.57e-05 5.78e-05 -5.50e-05 ...
# $ t.Total_LocoScore_SimplerModel_wDissectionDay          : num  0.503 -0.048 -0.851 -0.71 -0.135 ...
# $ P.value.Total_LocoScore_SimplerModel_wDissectionDay    : num  0.615 0.962 0.396 0.478 0.893 ...
# $ P.value.adj.Total_LocoScore_SimplerModel_wDissectionDay: num  0.999 0.999 0.999 0.999 0.999 ...
# $ Coef.EPM_Percent_Time_Open_Arm                         : num  -0.000409 -0.001465 -0.000238 -0.000266 0.000853 ...
# $ t.EPM_Percent_Time_Open_Arm                            : num  -0.577 -1.648 -0.203 -0.339 0.869 ...
# $ P.value.EPM_Percent_Time_Open_Arm                      : num  0.565 0.101 0.839 0.735 0.386 ...
# $ P.value.adj.EPM_Percent_Time_Open_Arm                  : num  0.994 0.994 0.997 0.996 0.994 ...
# $ Coef.EPM_DistanceTraveled                              : num  -1.01e-05 9.58e-06 -4.16e-06 -8.10e-06 -1.99e-06 -2.15e-05 1.47e-06 1.45e-05 5.79e-06 -1.16e-05 ...
# $ t.EPM_DistanceTraveled                                 : num  -0.654 0.4937 -0.1657 -0.4705 -0.0932 ...
# $ P.value.EPM_DistanceTraveled                           : num  0.514 0.622 0.869 0.638 0.926 ...
# $ P.value.adj.EPM_DistanceTraveled                       : num  1 1 1 1 1 ...
# $ Coef.LearningClassification_AsFactorGT                 : num  -0.001492 -0.00673 0.009725 0.024307 0.000185 ...
# $ Coef.LearningClassification_AsFactorIN                 : num  -0.00635 0.0242 0.03452 0.01126 0.03319 ...
# $ t.LearningClassification_AsFactorGT                    : num  -0.06379 -0.22906 0.24601 0.94721 0.00562 ...
# $ t.LearningClassification_AsFactorIN                    : num  -0.345 1.048 1.12 0.555 1.286 ...
# $ P.value.LearningClassification_AsFactorGT              : num  0.949 0.819 0.806 0.344 0.996 ...
# $ P.value.LearningClassification_AsFactorIN              : num  0.73 0.296 0.264 0.58 0.2 ...
# $ P.value.adj.LearningClassification_AsFactorGT          : num  1 1 1 1 1 ...
# $ P.value.adj.LearningClassification_AsFactorIN          : num  0.959 0.903 0.902 0.932 0.902 ...
# $ Coef.EPM_DistanceTraveled_noOutlier                    : num  -1.49e-05 -7.52e-06 -8.51e-06 -2.30e-05 -6.65e-06 -2.35e-05 -6.24e-06 2.52e-05 4.46e-05 -1.47e-05 ...
# $ t.EPM_DistanceTraveled_noOutlier                       : num  -0.869 -0.349 -0.302 -1.22 -0.28 ...
# $ P.value.EPM_DistanceTraveled_noOutlier                 : num  0.385 0.727 0.763 0.224 0.78 ...
# $ P.value.adj.EPM_DistanceTraveled_noOutlier             : num  1 1 1 1 1 ...
# $ Coef.PCA_Index_Days6and.7                              : num  0.00202 -0.00481 -0.02169 0.00367 0.00314 ...
# $ t.PCA_Index_Days6and.7                                 : num  0.132 -0.251 -0.875 0.224 0.152 ...
# $ P.value.PCA_Index_Days6and.7                           : num  0.895 0.802 0.383 0.823 0.88 ...
# $ P.value.adj.PCA_Index_Days6and.7                       : num  1 1 1 1 1 ...
# $ Coef.EPM_Time_Immobile                                 : num  -8.36e-05 -5.85e-05 6.00e-04 5.30e-05 8.89e-05 ...
# $ t.EPM_Time_Immobile                                    : num  -0.259 -0.143 1.111 0.149 0.196 ...
# $ P.value.EPM_Time_Immobile                              : num  0.796 0.886 0.268 0.882 0.845 ...
# $ P.value.adj.EPM_Time_Immobile                          : num  0.998 0.998 0.998 0.998 0.998 ...

library(plyr)

str(list(F0_Limma_Results_w_LateGeneMeta, F2_Limma_Results_from4Models, F2_Limma_Results_More_20220121_Summarized))

TempList<-list(F0_Limma_Results_w_LateGeneMeta, F2_Limma_Results_from4Models, F2_Limma_Results_More_20220121_Summarized)

F0_Meta_F2_4Models<-join_all(dfs=TempList, by="ENSEMBL", type="left")
#Note: I didn't update this object name to reflect the fact that I now also have output from (yet another) locomotor activity model, EPM, ST/GT

str(F0_Meta_F2_4Models)
#'data.frame':	13786 obs. of  225 variables:
#Well that's a marvelous beast, LOL.

#write.csv(F0_Meta_F2_4Models, "LimmaResults_F0_Meta_F2_5ModelsLocomotor_EPM_STGT.csv")

write.csv(F0_Meta_F2_4Models, "LimmaResults_F0_Meta_F2_5ModelsLocomotor_EPM_STGT_PavCA.csv")


#Correlation between T-stats in the different datasets:

cor(as.matrix(cbind(F0_Meta_F2_4Models$t.Lineage_AsFactorbLR, F0_Meta_F2_4Models$estimate, F0_Meta_F2_4Models$M6_CovJustSexIntergenic_t.Total_LocoScore, F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore, F0_Meta_F2_4Models$M8_CovSexRNAMetricsTechnical_t.Total_LocoScore, F0_Meta_F2_4Models$M7_CovSexRNAMetricsUTRConfoundTechnical_t.Total_LocoScore, F0_Meta_F2_4Models$t.Total_LocoScore_SimplerModel_wDissectionDay, F0_Meta_F2_4Models$t.EPM_Percent_Time_Open_Arm, F0_Meta_F2_4Models$t.EPM_DistanceTraveled, F0_Meta_F2_4Models$t.LearningClassification_AsFactorGT, F0_Meta_F2_4Models$t.LearningClassification_AsFactorIN,  F0_Meta_F2_4Models$t.EPM_DistanceTraveled_noOutlier, F0_Meta_F2_4Models$t.EPM_Time_Immobile, F0_Meta_F2_4Models$t.PCA_Index_Days6and.7)), use="pairwise.complete")

# [,1]        [,2]        [,3]       [,4]       [,5]       [,6]       [,7]        [,8]        [,9]
# [1,]  1.000000000  0.37571807 -0.25910977 -0.2212175 -0.2489206 -0.2176587 -0.1747805 -0.06101927 -0.22277762
# [2,]  0.375718070  1.00000000 -0.27316426 -0.1778227 -0.2073013 -0.2289370 -0.1584199 -0.04239608 -0.22422024
# [3,] -0.259109770 -0.27316426  1.00000000  0.8334972  0.7543923  0.7768770  0.7911415 -0.04183926  0.40851307
# [4,] -0.221217486 -0.17782270  0.83349718  1.0000000  0.9229983  0.9149396  0.9784730  0.23131285  0.36447488
# [5,] -0.248920584 -0.20730133  0.75439226  0.9229983  1.0000000  0.9697818  0.8801067  0.21124707  0.28440186
# [6,] -0.217658677 -0.22893700  0.77687701  0.9149396  0.9697818  1.0000000  0.8878354  0.16335389  0.33440947
# [7,] -0.174780515 -0.15841992  0.79114150  0.9784730  0.8801067  0.8878354  1.0000000  0.24070997  0.32928380
# [8,] -0.061019274 -0.04239608 -0.04183926  0.2313129  0.2112471  0.1633539  0.2407100  1.00000000  0.23948151
# [9,] -0.222777624 -0.22422024  0.40851307  0.3644749  0.2844019  0.3344095  0.3292838  0.23948151  1.00000000
# [10,] -0.021731100  0.04204025 -0.28158903 -0.3856298 -0.3484861 -0.3556784 -0.4126535 -0.07909066 -0.13743562
# [11,] -0.006838335  0.06905854 -0.04981391 -0.2223729 -0.3573524 -0.2527331 -0.2103689 -0.31652146  0.08376829
# [12,] -0.237829616 -0.24250213  0.34723461  0.4241617  0.3543640  0.3668509  0.3954705  0.45003539  0.88575437
# [13,]  0.283683698  0.10477430 -0.49116304 -0.5709958 -0.5337380 -0.4408689 -0.5176775 -0.25399701 -0.49827978
# [14,]  0.100984550 -0.04958306  0.21382213  0.2744848  0.3017771  0.3587527  0.3296361  0.02501600  0.06523872
# [,10]        [,11]       [,12]       [,13]       [,14]
# [1,] -0.02173110 -0.006838335 -0.23782962  0.28368370  0.10098455
# [2,]  0.04204025  0.069058538 -0.24250213  0.10477430 -0.04958306
# [3,] -0.28158903 -0.049813905  0.34723461 -0.49116304  0.21382213
# [4,] -0.38562982 -0.222372912  0.42416173 -0.57099582  0.27448478
# [5,] -0.34848605 -0.357352437  0.35436401 -0.53373802  0.30177713
# [6,] -0.35567837 -0.252733054  0.36685093 -0.44086892  0.35875265
# [7,] -0.41265345 -0.210368877  0.39547052 -0.51767753  0.32963612
# [8,] -0.07909066 -0.316521458  0.45003539 -0.25399701  0.02501600
# [9,] -0.13743562  0.083768288  0.88575437 -0.49827978  0.06523872
# [10,]  1.00000000  0.422223154 -0.21108223  0.26001711 -0.84568731
# [11,]  0.42222315  1.000000000 -0.08345006  0.11791258 -0.43048066
# [12,] -0.21108223 -0.083450065  1.00000000 -0.53756642  0.09649099
# [13,]  0.26001711  0.117912580 -0.53756642  1.00000000 -0.08125278
# [14,] -0.84568731 -0.430480657  0.09649099 -0.08125278  1.00000000


#Honestly, there's more signal there than I expected.
#Positive correlation between F0_Lineage_Tstats and LateGene RNA-Seq HR/LR Meta-analysis coefficents (R=0.35)
#Negative correlation between F0_Lineage_Tstats and F2_Locomotor_Tstats from all models. Interestingly, the relationship is actually strongest with the simplest model (?)
#Much weaker pattern for %TimeOpenArms_EPM & Distance_EPM: negative correlation with F0_Lineage, Positive correlation with locomotor score
#Much weaker pattern for STGT: GT has negative correlation with locomotor activity, positive correlation with IN (nice validation).

TempCorMatrix<-cor(as.matrix(cbind(F0_Meta_F2_4Models$t.Lineage_AsFactorbLR, F0_Meta_F2_4Models$estimate, F0_Meta_F2_4Models$M6_CovJustSexIntergenic_t.Total_LocoScore, F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore, F0_Meta_F2_4Models$M8_CovSexRNAMetricsTechnical_t.Total_LocoScore, F0_Meta_F2_4Models$M7_CovSexRNAMetricsUTRConfoundTechnical_t.Total_LocoScore, F0_Meta_F2_4Models$t.Total_LocoScore_SimplerModel_wDissectionDay, F0_Meta_F2_4Models$t.EPM_Percent_Time_Open_Arm, F0_Meta_F2_4Models$t.EPM_DistanceTraveled, F0_Meta_F2_4Models$t.LearningClassification_AsFactorGT, F0_Meta_F2_4Models$t.LearningClassification_AsFactorIN,  F0_Meta_F2_4Models$t.EPM_DistanceTraveled_noOutlier, F0_Meta_F2_4Models$t.EPM_Time_Immobile, F0_Meta_F2_4Models$t.PCA_Index_Days6and.7)), use="pairwise.complete")

#write.csv(TempCorMatrix, "CorMatrix_Tstats_LimmaResults_F0_Meta_F2.csv")
write.csv(TempCorMatrix, "CorMatrix_Tstats_LimmaResults_F0_Meta_F2_updated.csv")

#Non-parametric/Spearman version (with only the variables of biggest interest):
TempCorMatrix<-cor(as.matrix(cbind(F0_Meta_F2_4Models$t.Lineage_AsFactorbLR, F0_Meta_F2_4Models$estimate, F0_Meta_F2_4Models$t.EPM_Time_Immobile, F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore, F0_Meta_F2_4Models$t.EPM_DistanceTraveled_noOutlier, F0_Meta_F2_4Models$t.EPM_Percent_Time_Open_Arm, F0_Meta_F2_4Models$t.PCA_Index_Days6and.7)), use="pairwise.complete", method="spearman")

#write.csv(TempCorMatrix, "CorMatrix_Tstats_LimmaResults_F0_Meta_F2.csv")
write.csv(TempCorMatrix, "CorMatrix_Tstats_LimmaResults_F0_Meta_F2_updated_Spearman.csv")

cor.test(F0_Meta_F2_4Models$t.Lineage_AsFactorbLR, F0_Meta_F2_4Models$estimate, use="pairwise.complete", method="spearman")

# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$t.Lineage_AsFactorbLR and F0_Meta_F2_4Models$estimate
# S = 1.6093e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.308112

cor.test(F0_Meta_F2_4Models$t.Lineage_AsFactorbLR, F0_Meta_F2_4Models$t.EPM_Time_Immobile, use="pairwise.complete", method="spearman")
# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$t.Lineage_AsFactorbLR and F0_Meta_F2_4Models$t.EPM_Time_Immobile
# S = 2.7887e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.2950011 

cor.test(F0_Meta_F2_4Models$estimate, F0_Meta_F2_4Models$t.EPM_Time_Immobile, use="pairwise.complete", method="spearman")
# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$estimate and F0_Meta_F2_4Models$t.EPM_Time_Immobile
# S = 2.072e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.0798225 

cor.test(F0_Meta_F2_4Models$t.Lineage_AsFactorbLR, F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore, use="pairwise.complete", method="spearman")
# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$t.Lineage_AsFactorbLR and F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore
# S = 4.8563e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# -0.2276814 

cor.test(F0_Meta_F2_4Models$estimate, F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore, use="pairwise.complete", method="spearman")
# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$estimate and F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore
# S = 2.5791e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# -0.1453604 

cor.test(F0_Meta_F2_4Models$t.Lineage_AsFactorbLR, F0_Meta_F2_4Models$t.EPM_DistanceTraveled_noOutlier, use="pairwise.complete", method="spearman")
# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$t.Lineage_AsFactorbLR and F0_Meta_F2_4Models$t.EPM_DistanceTraveled_noOutlier
# S = 4.8992e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# -0.2385234 

cor.test(F0_Meta_F2_4Models$estimate, F0_Meta_F2_4Models$t.EPM_DistanceTraveled_noOutlier, use="pairwise.complete", method="spearman")
# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$estimate and F0_Meta_F2_4Models$t.EPM_DistanceTraveled_noOutlier
# S = 2.7823e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# -0.2355979 

cor.test(F0_Meta_F2_4Models$t.Lineage_AsFactorbLR, F0_Meta_F2_4Models$t.EPM_Percent_Time_Open_Arm, use="pairwise.complete", method="spearman")
# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$t.Lineage_AsFactorbLR and F0_Meta_F2_4Models$t.EPM_Percent_Time_Open_Arm
# S = 4.1596e+11, p-value = 2.549e-09
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# -0.05155981 

cor.test(F0_Meta_F2_4Models$estimate, F0_Meta_F2_4Models$t.EPM_Percent_Time_Open_Arm, use="pairwise.complete", method="spearman")
# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$estimate and F0_Meta_F2_4Models$t.EPM_Percent_Time_Open_Arm
# S = 2.3074e+11, p-value = 0.009428
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# -0.0246904 

cor.test(F0_Meta_F2_4Models$t.Lineage_AsFactorbLR, F0_Meta_F2_4Models$t.PCA_Index_Days6and.7, use="pairwise.complete", method="spearman")

# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$t.Lineage_AsFactorbLR and F0_Meta_F2_4Models$t.PCA_Index_Days6and.7
# S = 3.4732e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.1219682 

cor.test(F0_Meta_F2_4Models$estimate, F0_Meta_F2_4Models$t.PCA_Index_Days6and.7, use="pairwise.complete", method="spearman")
# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$estimate and F0_Meta_F2_4Models$t.PCA_Index_Days6and.7
# S = 2.3522e+11, p-value = 2.748e-06
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# -0.04457732 


#Correlation between Coef in the different datasets:

cor(as.matrix(cbind(F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR, F0_Meta_F2_4Models$estimate, F0_Meta_F2_4Models$M6_CovJustSexIntergenic_Coef.Total_LocoScore, F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore, F0_Meta_F2_4Models$M8_CovSexRNAMetricsTechnical_Coef.Total_LocoScore, F0_Meta_F2_4Models$M7_CovSexRNAMetricsUTRConfoundTechnical_Coef.Total_LocoScore, F0_Meta_F2_4Models$Coef.Total_LocoScore_SimplerModel_wDissectionDay, F0_Meta_F2_4Models$Coef.EPM_Percent_Time_Open_Arm, F0_Meta_F2_4Models$Coef.EPM_DistanceTraveled, F0_Meta_F2_4Models$Coef.LearningClassification_AsFactorGT, F0_Meta_F2_4Models$Coef.LearningClassification_AsFactorIN, F0_Meta_F2_4Models$Coef.EPM_DistanceTraveled_noOutlier, F0_Meta_F2_4Models$Coef.EPM_Time_Immobile, F0_Meta_F2_4Models$Coef.PCA_Index_Days6and.7)), use="pairwise.complete")
# [,1]        [,2]       [,3]       [,4]       [,5]       [,6]       [,7]        [,8]       [,9]
# [1,]  1.00000000  0.34545591 -0.2305792 -0.1958006 -0.2016754 -0.1785550 -0.1479704 -0.09041479 -0.2308912
# [2,]  0.34545591  1.00000000 -0.2080778 -0.1536555 -0.1696513 -0.1826409 -0.1436621 -0.05788521 -0.1792753
# [3,] -0.23057916 -0.20807780  1.0000000  0.8826453  0.8128090  0.8092741  0.8352606  0.11285811  0.4029076
# [4,] -0.19580056 -0.15365545  0.8826453  1.0000000  0.9578173  0.9470364  0.9799447  0.25373248  0.3750450
# [5,] -0.20167540 -0.16965133  0.8128090  0.9578173  1.0000000  0.9784620  0.9343210  0.25657312  0.3249628
# [6,] -0.17855500 -0.18264087  0.8092741  0.9470364  0.9784620  1.0000000  0.9381076  0.25023397  0.3509199
# [7,] -0.14797036 -0.14366211  0.8352606  0.9799447  0.9343210  0.9381076  1.0000000  0.25386703  0.3278227
# [8,] -0.09041479 -0.05788521  0.1128581  0.2537325  0.2565731  0.2502340  0.2538670  1.00000000  0.4379232
# [9,] -0.23089118 -0.17927526  0.4029076  0.3750450  0.3249628  0.3509199  0.3278227  0.43792324  1.0000000
# [10,]  0.06392852  0.03688858 -0.4066194 -0.3681757 -0.3284367 -0.3107149 -0.3436998 -0.11907187 -0.2228186
# [11,]  0.07437494  0.07158776 -0.1386096 -0.2070326 -0.2939380 -0.2317952 -0.1936002 -0.25358092 -0.1033762
# [12,] -0.26538972 -0.19532530  0.3783379  0.4012016  0.3609862  0.3696257  0.3568465  0.55687538  0.9151675
# [13,]  0.32222761  0.09653171 -0.5422123 -0.5192262 -0.4690136 -0.3934191 -0.4397952 -0.25512086 -0.5858664
# [14,] -0.05568478 -0.05227946  0.3848672  0.3576502  0.3545978  0.3796800  0.3691855  0.12480321  0.2334506
# [,10]       [,11]      [,12]       [,13]       [,14]
# [1,]  0.06392852  0.07437494 -0.2653897  0.32222761 -0.05568478
# [2,]  0.03688858  0.07158776 -0.1953253  0.09653171 -0.05227946
# [3,] -0.40661941 -0.13860964  0.3783379 -0.54221234  0.38486724
# [4,] -0.36817568 -0.20703258  0.4012016 -0.51922617  0.35765022
# [5,] -0.32843670 -0.29393805  0.3609862 -0.46901364  0.35459783
# [6,] -0.31071492 -0.23179520  0.3696257 -0.39341906  0.37968003
# [7,] -0.34369984 -0.19360020  0.3568465 -0.43979522  0.36918550
# [8,] -0.11907187 -0.25358092  0.5568754 -0.25512086  0.12480321
# [9,] -0.22281865 -0.10337620  0.9151675 -0.58586635  0.23345065
# [10,]  1.00000000  0.42072333 -0.2762660  0.32941623 -0.84517157
# [11,]  0.42072333  1.00000000 -0.1920410  0.06110788 -0.46397635
# [12,] -0.27626602 -0.19204099  1.0000000 -0.59871443  0.25874230
# [13,]  0.32941623  0.06110788 -0.5987144  1.00000000 -0.22969414
# [14,] -0.84517157 -0.46397635  0.2587423 -0.22969414  1.00000000

#I wish the meta-analysis had been run with log2FC instead of cohen's d so that it would be easier to directly compare effect sizes.
#Maybe re-run it before pulling together the final publication - or run a meta-analysis that combines all three available datasets.


TempCorMatrix<-cor(as.matrix(cbind(F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR, F0_Meta_F2_4Models$estimate, F0_Meta_F2_4Models$M6_CovJustSexIntergenic_Coef.Total_LocoScore, F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore, F0_Meta_F2_4Models$M8_CovSexRNAMetricsTechnical_Coef.Total_LocoScore, F0_Meta_F2_4Models$M7_CovSexRNAMetricsUTRConfoundTechnical_Coef.Total_LocoScore, F0_Meta_F2_4Models$Coef.Total_LocoScore_SimplerModel_wDissectionDay, F0_Meta_F2_4Models$Coef.EPM_Percent_Time_Open_Arm, F0_Meta_F2_4Models$Coef.EPM_DistanceTraveled, F0_Meta_F2_4Models$Coef.LearningClassification_AsFactorGT, F0_Meta_F2_4Models$Coef.LearningClassification_AsFactorIN, F0_Meta_F2_4Models$Coef.EPM_DistanceTraveled_noOutlier, F0_Meta_F2_4Models$Coef.EPM_Time_Immobile, F0_Meta_F2_4Models$Coef.PCA_Index_Days6and.7)), use="pairwise.complete")

#write.csv(TempCorMatrix, "CorMatrix_Coef_LimmaResults_F0_Meta_F2.csv")
write.csv(TempCorMatrix, "CorMatrix_Coef_LimmaResults_F0_Meta_F2_updated.csv")

#Using non-parametric/Spearman:
TempCorMatrix<-cor(as.matrix(cbind(F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR, F0_Meta_F2_4Models$estimate, F0_Meta_F2_4Models$Coef.EPM_Time_Immobile, F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore, F0_Meta_F2_4Models$Coef.EPM_DistanceTraveled_noOutlier, F0_Meta_F2_4Models$Coef.EPM_Percent_Time_Open_Arm,  F0_Meta_F2_4Models$Coef.PCA_Index_Days6and.7)), use="pairwise.complete", method="spearman")

#write.csv(TempCorMatrix, "CorMatrix_Coef_LimmaResults_F0_Meta_F2.csv")
write.csv(TempCorMatrix, "CorMatrix_Coef_LimmaResults_F0_Meta_F2_updated_Spearman.csv")


cor.test(F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR, F0_Meta_F2_4Models$estimate, use="pairwise.complete", method="spearman")
# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR and F0_Meta_F2_4Models$estimate
# S = 1.6332e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.29782 

cor.test(F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR, F0_Meta_F2_4Models$Coef.EPM_Time_Immobile, use="pairwise.complete", method="spearman")

# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR and F0_Meta_F2_4Models$Coef.EPM_Time_Immobile
# S = 2.7961e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.2931336 

cor.test(F0_Meta_F2_4Models$estimate, F0_Meta_F2_4Models$Coef.EPM_Time_Immobile, use="pairwise.complete", method="spearman")

# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$estimate and F0_Meta_F2_4Models$Coef.EPM_Time_Immobile
# S = 2.0883e+11, p-value = 2.11e-14
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.07261077 

cor.test(F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR, F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore, use="pairwise.complete", method="spearman")

# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR and F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore
# S = 4.8614e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# -0.2289803 

cor.test(F0_Meta_F2_4Models$estimate, F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore, use="pairwise.complete", method="spearman")

# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$estimate and F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore
# S = 2.5589e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# -0.1363942 

cor.test(F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR, F0_Meta_F2_4Models$Coef.EPM_DistanceTraveled_noOutlier, use="pairwise.complete", method="spearman")

# Spearman's rank correlation rho
# 
# data:  F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR and F0_Meta_F2_4Models$Coef.EPM_DistanceTraveled_noOutlier
# S = 4.8225e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# -0.2191319 


cor.test(F0_Meta_F2_4Models$estimate, F0_Meta_F2_4Models$Coef.EPM_DistanceTraveled_noOutlier, use="pairwise.complete", method="spearman")


cor.test(F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR, F0_Meta_F2_4Models$Coef.EPM_Percent_Time_Open_Arm, use="pairwise.complete", method="spearman")


cor.test(F0_Meta_F2_4Models$estimate, F0_Meta_F2_4Models$Coef.EPM_Percent_Time_Open_Arm, use="pairwise.complete", method="spearman")


cor.test(F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR, F0_Meta_F2_4Models$Coef.PCA_Index_Days6and.7, use="pairwise.complete", method="spearman")


cor.test(F0_Meta_F2_4Models$estimate, F0_Meta_F2_4Models$Coef.PCA_Index_Days6and.7, use="pairwise.complete", method="spearman")


#Outputting some parametric and non-parametric stats to go with that: 

summary.lm(lm(Coef.LearningClassification_AsFactorGT~estimate, data=F0_Meta_F2_4Models))
# Call:
#   lm(formula = Coef.LearningClassification_AsFactorGT ~ estimate, 
#      data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.91040 -0.02164 -0.00140  0.01938  0.65336 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.0032987  0.0004346  -7.590 3.45e-14 ***
#   estimate     0.0023142  0.0005963   3.881 0.000105 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.04563 on 11053 degrees of freedom
# (2731 observations deleted due to missingness)
# Multiple R-squared:  0.001361,	Adjusted R-squared:  0.00127 
# F-statistic: 15.06 on 1 and 11053 DF,  p-value: 0.0001047

summary.lm(lm(Coef.LearningClassification_AsFactorGT~Coef.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models))
# Call:
#   lm(formula = Coef.LearningClassification_AsFactorGT ~ Coef.Lineage_AsFactorbLR, 
#      data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.94516 -0.02254 -0.00169  0.02010  0.68364 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              -0.0029123  0.0004297  -6.778 1.27e-11 ***
#   Coef.Lineage_AsFactorbLR  0.0136283  0.0018422   7.398 1.47e-13 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.04952 on 13337 degrees of freedom
# (447 observations deleted due to missingness)
# Multiple R-squared:  0.004087,	Adjusted R-squared:  0.004012 
# F-statistic: 54.73 on 1 and 13337 DF,  p-value: 1.465e-13

summary.lm(lm(Coef.LearningClassification_AsFactorIN~estimate, data=F0_Meta_F2_4Models))

# Call:
#   lm(formula = Coef.LearningClassification_AsFactorIN ~ estimate, 
#      data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.63527 -0.02086  0.00548  0.02358  0.31701 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.0042993  0.0003770 -11.405  < 2e-16 ***
#   estimate     0.0039030  0.0005173   7.546 4.86e-14 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.03958 on 11053 degrees of freedom
# (2731 observations deleted due to missingness)
# Multiple R-squared:  0.005125,	Adjusted R-squared:  0.005035 
# F-statistic: 56.94 on 1 and 11053 DF,  p-value: 4.855e-14

summary.lm(lm(Coef.LearningClassification_AsFactorIN~Coef.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models))
# Call:
#   lm(formula = Coef.LearningClassification_AsFactorIN ~ Coef.Lineage_AsFactorbLR, 
#      data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.61760 -0.02127  0.00524  0.02374  0.38803 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              -0.0036931  0.0003629 -10.177   <2e-16 ***
#   Coef.Lineage_AsFactorbLR  0.0134009  0.0015559   8.613   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.04182 on 13337 degrees of freedom
# (447 observations deleted due to missingness)
# Multiple R-squared:  0.005532,	Adjusted R-squared:  0.005457 
# F-statistic: 74.19 on 1 and 13337 DF,  p-value: < 2.2e-16

summary.lm(lm(Coef.EPM_Percent_Time_Open_Arm~estimate, data=F0_Meta_F2_4Models))

# Call:
#   lm(formula = Coef.EPM_Percent_Time_Open_Arm ~ estimate, data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.032147 -0.000753 -0.000005  0.000730  0.020954 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.573e-04  1.470e-05  10.696  < 2e-16 ***
#   estimate    -1.230e-04  2.017e-05  -6.096 1.12e-09 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.001544 on 11053 degrees of freedom
# (2731 observations deleted due to missingness)
# Multiple R-squared:  0.003351,	Adjusted R-squared:  0.003261 
# F-statistic: 37.16 on 1 and 11053 DF,  p-value: 1.125e-09

summary.lm(lm(Coef.EPM_Percent_Time_Open_Arm~Coef.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models))
# Call:
#   lm(formula = Coef.EPM_Percent_Time_Open_Arm ~ Coef.Lineage_AsFactorbLR, 
#      data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.033025 -0.000778  0.000015  0.000762  0.022696 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               1.002e-04  1.463e-05    6.85 7.71e-12 ***
#   Coef.Lineage_AsFactorbLR -6.575e-04  6.271e-05  -10.48  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.001686 on 13337 degrees of freedom
# (447 observations deleted due to missingness)
# Multiple R-squared:  0.008175,	Adjusted R-squared:  0.0081 
# F-statistic: 109.9 on 1 and 13337 DF,  p-value: < 2.2e-16

summary.lm(lm(Coef.EPM_DistanceTraveled~estimate, data=F0_Meta_F2_4Models))

# Call:
#   lm(formula = Coef.EPM_DistanceTraveled ~ estimate, data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -8.612e-04 -1.318e-05  2.800e-07  1.283e-05  5.950e-04 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.153e-06  3.052e-07   10.33   <2e-16 ***
#   estimate    -8.023e-06  4.188e-07  -19.16   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 3.205e-05 on 11053 degrees of freedom
# (2731 observations deleted due to missingness)
# Multiple R-squared:  0.03214,	Adjusted R-squared:  0.03205 
# F-statistic:   367 on 1 and 11053 DF,  p-value: < 2.2e-16

summary.lm(lm(Coef.EPM_DistanceTraveled~Coef.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models))
# Call:
#   lm(formula = Coef.EPM_DistanceTraveled ~ Coef.Lineage_AsFactorbLR, 
#      data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -9.051e-04 -1.353e-05  1.500e-07  1.290e-05  4.904e-04 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               3.400e-06  2.939e-07   11.57   <2e-16 ***
#   Coef.Lineage_AsFactorbLR -3.453e-05  1.260e-06  -27.41   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 3.387e-05 on 13337 degrees of freedom
# (447 observations deleted due to missingness)
# Multiple R-squared:  0.05331,	Adjusted R-squared:  0.05324 
# F-statistic:   751 on 1 and 13337 DF,  p-value: < 2.2e-16


summary.lm(lm(Coef.EPM_DistanceTraveled_noOutlier~estimate, data=F0_Meta_F2_4Models))

# Call:
#   lm(formula = Coef.EPM_DistanceTraveled_noOutlier ~ estimate, 
#      data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -8.464e-04 -1.512e-05 -7.800e-07  1.360e-05  6.035e-04 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.325e-06  3.437e-07   9.674   <2e-16 ***
#   estimate    -9.875e-06  4.716e-07 -20.938   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.609e-05 on 11053 degrees of freedom
# (2731 observations deleted due to missingness)
# Multiple R-squared:  0.03815,	Adjusted R-squared:  0.03806 
# F-statistic: 438.4 on 1 and 11053 DF,  p-value: < 2.2e-16

sqrt(0.03815)
#[1] 0.1953202

summary.lm(lm(Coef.EPM_DistanceTraveled_noOutlier~Coef.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models))

# Call:
#   lm(formula = Coef.EPM_DistanceTraveled_noOutlier ~ Coef.Lineage_AsFactorbLR, 
#      data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -9.023e-04 -1.545e-05 -6.800e-07  1.377e-05  5.171e-04 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               3.168e-06  3.315e-07   9.558   <2e-16 ***
#   Coef.Lineage_AsFactorbLR -4.518e-05  1.421e-06 -31.789   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.82e-05 on 13337 degrees of freedom
# (447 observations deleted due to missingness)
# Multiple R-squared:  0.07043,	Adjusted R-squared:  0.07036 
# F-statistic:  1011 on 1 and 13337 DF,  p-value: < 2.2e-16

sqrt(0.07043)
#[1] 0.2653865

summary.lm(lm(Coef.EPM_Time_Immobile~estimate, data=F0_Meta_F2_4Models))

# Call:
#   lm(formula = Coef.EPM_Time_Immobile ~ estimate, data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0137940 -0.0003739 -0.0000238  0.0002874  0.0149527 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.159e-05  7.761e-06  -1.493    0.135    
# estimate     1.086e-04  1.065e-05  10.196   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0008149 on 11053 degrees of freedom
# (2731 observations deleted due to missingness)
# Multiple R-squared:  0.009318,	Adjusted R-squared:  0.009229 
# F-statistic:   104 on 1 and 11053 DF,  p-value: < 2.2e-16

sqrt(0.009318)
#[1] 0.09652979

summary.lm(lm(Coef.EPM_Time_Immobile~Coef.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models))

# Call:
#   lm(formula = Coef.EPM_Time_Immobile ~ Coef.Lineage_AsFactorbLR, 
#      data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0099343 -0.0003807 -0.0000367  0.0002706  0.0162063 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              -9.086e-06  7.284e-06  -1.247    0.212    
# Coef.Lineage_AsFactorbLR  1.228e-03  3.123e-05  39.309   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0008394 on 13337 degrees of freedom
# (447 observations deleted due to missingness)
# Multiple R-squared:  0.1038,	Adjusted R-squared:  0.1038 
# F-statistic:  1545 on 1 and 13337 DF,  p-value: < 2.2e-16

sqrt(0.1038)
#[1] 0.3221801


summary.lm(lm(Coef.PCA_Index_Days6and.7~estimate, data=F0_Meta_F2_4Models))
# Call:
#   lm(formula = Coef.PCA_Index_Days6and.7 ~ estimate, data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.33489 -0.01377  0.00007  0.01423  0.53230 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.0010646  0.0002898   3.674  0.00024 ***
#   estimate    -0.0021886  0.0003976  -5.504  3.8e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03043 on 11053 degrees of freedom
# (2731 observations deleted due to missingness)
# Multiple R-squared:  0.002733,	Adjusted R-squared:  0.002643 
# F-statistic: 30.29 on 1 and 11053 DF,  p-value: 3.799e-08

sqrt(0.002733)
#[1] 0.0522781


summary.lm(lm(Coef.PCA_Index_Days6and.7~Coef.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models))

# Call:
#   lm(formula = Coef.PCA_Index_Days6and.7 ~ Coef.Lineage_AsFactorbLR, 
#      data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.38937 -0.01441 -0.00016  0.01489  0.54972 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               0.0011533  0.0002848   4.049 5.17e-05 ***
#   Coef.Lineage_AsFactorbLR -0.0078654  0.0012212  -6.441 1.23e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03282 on 13337 degrees of freedom
# (447 observations deleted due to missingness)
# Multiple R-squared:  0.003101,	Adjusted R-squared:  0.003026 
# F-statistic: 41.48 on 1 and 13337 DF,  p-value: 1.229e-10

sqrt(0.003101)
#[1] 0.05568662

#plotting the relationship between LogFC for Total Loco Score and Lineage:

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01")

#Unfiltered: 

pdf("Scatterplot_LogFC_F2_Locomotor_vs_Meta_Lineage.pdf", width=5, height=6)
plot(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~estimate, data=F0_Meta_F2_4Models, xlab="bLR vs. bHR: Late Generation Meta-Analysis Estimated d", ylab="Locomotor Activity: F2 Log2FC")
TempLine<-lm(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~estimate, data=F0_Meta_F2_4Models)
abline(TempLine, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore ~ 
#        estimate, data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -5.408e-04 -2.022e-05  3.000e-06  2.382e-05  7.992e-04 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.567e-06  4.999e-07   3.134  0.00173 ** 
#   estimate    -1.121e-05  6.859e-07 -16.348  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 5.249e-05 on 11053 degrees of freedom
# (2731 observations deleted due to missingness)
# Multiple R-squared:  0.02361,	Adjusted R-squared:  0.02352 
# F-statistic: 267.3 on 1 and 11053 DF,  p-value: < 2.2e-16

F0_Meta_F2_4Models$MetaAnalysis_Tstat<-F0_Meta_F2_4Models$estimate/F0_Meta_F2_4Models$SE

pdf("Scatterplot_Tstat_F2_Locomotor_vs_Meta_Lineage.pdf", width=5, height=6)
plot(M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore~MetaAnalysis_Tstat, data=F0_Meta_F2_4Models, xlab="bLR vs. bHR: Late Generation Meta-Analysis T-stat", ylab="Locomotor Activity: F2 T-stat")
TempLine<-lm(M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore~MetaAnalysis_Tstat, data=F0_Meta_F2_4Models)
abline(TempLine, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore ~ 
#        MetaAnalysis_Tstat, data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.8940 -0.5806  0.0372  0.6013  4.0162 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.088207   0.008524   10.35   <2e-16 ***
#   MetaAnalysis_Tstat -0.106982   0.005701  -18.77   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.895 on 11053 degrees of freedom
# (2731 observations deleted due to missingness)
# Multiple R-squared:  0.03088,	Adjusted R-squared:  0.03079 
# F-statistic: 352.1 on 1 and 11053 DF,  p-value: < 2.2e-16


pdf("Scatterplot_LogFC_F2_Locomotor_vs_F0_Lineage.pdf", width=5, height=6)
plot(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~Coef.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models, xlab="bLR vs. bHR: F0 Log2FC ", ylab="Locomotor Activity: F2 Log2FC")
TempLine<-lm(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~Coef.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models)
abline(TempLine, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore ~ 
#        Coef.Lineage_AsFactorbLR, data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -1.080e-03 -2.117e-05  2.490e-06  2.394e-05  9.032e-04 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               2.035e-06  5.028e-07   4.047 5.21e-05 ***
#   Coef.Lineage_AsFactorbLR -4.971e-05  2.156e-06 -23.059  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 5.795e-05 on 13337 degrees of freedom
# (447 observations deleted due to missingness)
# Multiple R-squared:  0.03834,	Adjusted R-squared:  0.03827 
# F-statistic: 531.7 on 1 and 13337 DF,  p-value: < 2.2e-16

pdf("Scatterplot_Tstat_F2_Locomotor_vs_F0_Lineage.pdf", width=5, height=6)
plot(M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore~t.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models, xlab="bLR vs. bHR: F0 Tstat ", ylab="Locomotor Activity: F2 Tstat")
TempLine<-lm(M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore~t.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models)
abline(TempLine, lwd=3)
dev.off()

summary.lm(TempLine)

# Call:
#   lm(formula = M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore ~ 
#        t.Lineage_AsFactorbLR, data = F0_Meta_F2_4Models)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.6647 -0.5797  0.0268  0.5848  4.7542 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            0.094513   0.007812    12.1   <2e-16 ***
#   t.Lineage_AsFactorbLR -0.153421   0.005857   -26.2   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9004 on 13337 degrees of freedom
# (447 observations deleted due to missingness)
# Multiple R-squared:  0.04894,	Adjusted R-squared:  0.04887 
# F-statistic: 686.3 on 1 and 13337 DF,  p-value: < 2.2e-16


pdf("Scatterplot_Tstat_F2_Locomotor_vs_F0_Lineage_OutlierRemoved.pdf", width=5, height=6)
plot(M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore~t.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models[F0_Meta_F2_4Models$t.Lineage_AsFactorbLR<20,], xlab="bLR vs. bHR: F0 Tstat ", ylab="Locomotor Activity: F2 Tstat")
TempLine<-lm(M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore~t.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models[F0_Meta_F2_4Models$t.Lineage_AsFactorbLR<20,])
abline(TempLine, lwd=3)
dev.off()

summary.lm(TempLine)


#It would be nice to make a plot that contains all of these as layers (instead of individual plots for each)


pdf("Scatterplot_LogFC_F2_Locomotor_vs_Meta_Lineage_LayeredFilters.pdf", width=5, height=6)

plot(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~estimate, data=F0_Meta_F2_4Models, xlab="bLR vs. bHR: Late Generation Meta-Analysis Estimated d", ylab="Locomotor Activity: F2 Log2FC", col="darkgrey")
TempLine<-lm(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~estimate, data=F0_Meta_F2_4Models)
abline(TempLine, lwd=3, col="darkgrey")

points(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~estimate, data=F0_Meta_F2_4Models[((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE)),], col="black")
TempLine2<-lm(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~estimate, data=F0_Meta_F2_4Models[((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE)),])
abline(TempLine2, lwd=3, col="black")

points(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~estimate, data=F0_Meta_F2_4Models[((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) & (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE)),],  col="red")
TempLine3<-lm(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~estimate, data=F0_Meta_F2_4Models[((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) & (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE)),])
abline(TempLine3, lwd=3, col="red")

dev.off()

summary.lm(TempLine2)

# Call:
#   lm(formula = M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore ~ 
#        estimate, data = F0_Meta_F2_4Models[((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR < 
#                                                0.1 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR) == 
#                                                FALSE) | (F0_Meta_F2_4Models$BH < 0.1 & is.na(F0_Meta_F2_4Models$BH) == 
#                                                            FALSE)), ])
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -3.582e-04 -2.654e-05  3.660e-06  2.885e-05  6.549e-04 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -8.472e-07  2.493e-06  -0.340    0.734    
# estimate    -1.291e-05  1.399e-06  -9.225   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 7.154e-05 on 838 degrees of freedom
# (82 observations deleted due to missingness)
# Multiple R-squared:  0.09219,	Adjusted R-squared:  0.09111 
# F-statistic:  85.1 on 1 and 838 DF,  p-value: < 2.2e-16

summary.lm(TempLine3)

# Call:
#   lm(formula = M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore ~ 
#        estimate, data = F0_Meta_F2_4Models[((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR < 
#                                                0.1 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR) == 
#                                                FALSE) & (F0_Meta_F2_4Models$BH < 0.1 & is.na(F0_Meta_F2_4Models$BH) == 
#                                                            FALSE)), ])
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -2.612e-04 -9.175e-05 -1.227e-05  3.279e-05  6.157e-04 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  7.854e-06  1.749e-05   0.449  0.65504   
# estimate    -1.822e-05  6.022e-06  -3.025  0.00366 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0001375 on 60 degrees of freedom
# Multiple R-squared:  0.1323,	Adjusted R-squared:  0.1179 
# F-statistic: 9.151 on 1 and 60 DF,  p-value: 0.003656


#That's... kind of ugly. Let's try a simpler version (and fix the y-axis):

pdf("Scatterplot_LogFC_F2_Locomotor_vs_Meta_Lineage_LayeredFilters2.pdf", width=5, height=6)

plot(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~estimate, data=F0_Meta_F2_4Models, xlab="bLR vs. bHR: Late Generation Meta-Analysis Estimated d", ylab="Locomotor Activity: F2 Log2FC", col="black", ylim=c(-0.0008, 0.0008))
TempLine<-lm(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~estimate, data=F0_Meta_F2_4Models)
abline(TempLine, lwd=3, col="black")


points(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~estimate, data=F0_Meta_F2_4Models[((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) & (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE)),], col="red")
TempLine2<-lm(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~estimate, data=F0_Meta_F2_4Models[((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) & (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE)),])
abline(TempLine2, lwd=3, col="red")

dev.off()



pdf("Scatterplot_LogFC_F2_Locomotor_vs_F0_Lineage_LayeredFilters2.pdf", width=5, height=6)

plot(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~Coef.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models, xlab="bLR vs. bHR: F0 Log2FC", ylab="Locomotor Activity: F2 Log2FC", col="black")
TempLine<-lm(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~Coef.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models)
abline(TempLine, lwd=3, col="black")


points(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~Coef.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models[(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) & (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE),], col="red")

TempLine2<-lm(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~Coef.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models[(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) & (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE),])
abline(TempLine2, lwd=3, col="red")

dev.off()

# summary.lm(TempLine2)
# Call:
#   lm(formula = M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore ~ 
#        Coef.Lineage_AsFactorbLR, data = F0_Meta_F2_4Models[(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR < 
#                                                               0.1 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR) == 
#                                                               FALSE) & (F0_Meta_F2_4Models$BH < 0.1 & is.na(F0_Meta_F2_4Models$BH) == 
#                                                                           FALSE), ])
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -3.005e-04 -7.134e-05 -4.990e-06  3.471e-05  5.812e-04 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)               2.271e-06  1.808e-05   0.126   0.9005  
# Coef.Lineage_AsFactorbLR -3.113e-05  1.441e-05  -2.160   0.0348 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0001421 on 60 degrees of freedom
# Multiple R-squared:  0.07215,	Adjusted R-squared:  0.05668 
# F-statistic: 4.666 on 1 and 60 DF,  p-value: 0.03478

#1. This coding is ugly AF
#2. I should probably do a version that is more analagous to the way that we actually end up filtering results for overlap with QTLs (which has to be a little less stringent)




pdf("Scatterplot_LogFC_F2_Locomotor_vs_Meta_Lineage_LayeredFilters3.pdf", width=5, height=6)

plot(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~estimate, data=F0_Meta_F2_4Models, xlab="bLR vs. bHR: Late Generation Meta-Analysis Estimated d", ylab="Locomotor Activity: F2 Log2FC", col="black", ylim=c(-0.0008, 0.0008))
TempLine<-lm(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~estimate, data=F0_Meta_F2_4Models)
abline(TempLine, lwd=3, col="black")

#Subsetting data to require FDR<0.10 for bHR vs. bLR in either F0 or Late Gen Meta-Analysis or p<0.05 in both & a consistent direction of effect:
TempData<-F0_Meta_F2_4Models[((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE) | (F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR<0.05 & is.na(F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR)==FALSE & F0_Meta_F2_4Models$pval<0.05 & is.na(F0_Meta_F2_4Models$pval)==FALSE) & (F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR*F0_Meta_F2_4Models$estimate)>0),]

points(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~estimate, data=TempData, col="red")
TempLine2<-lm(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~estimate, data=TempData)
abline(TempLine2, lwd=3, col="red")

dev.off()

summary.lm(lm(TempLine2))
# Call:
#   lm(formula = TempLine2)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -5.356e-04 -2.704e-05  3.920e-06  3.055e-05  6.474e-04 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2.697e-07  2.451e-06  -0.110    0.912    
# estimate    -1.410e-05  1.442e-06  -9.783   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 7.591e-05 on 976 degrees of freedom
# (85 observations deleted due to missingness)
# Multiple R-squared:  0.0893,	Adjusted R-squared:  0.08837 
# F-statistic:  95.7 on 1 and 976 DF,  p-value: < 2.2e-16


pdf("Scatterplot_LogFC_F2_Locomotor_vs_F0_Lineage_LayeredFilters3.pdf", width=5, height=6)

plot(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~Coef.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models, xlab="bLR vs. bHR: F0 Log2FC", ylab="Locomotor Activity: F2 Log2FC", col="black", ylim=c(-0.0008, 0.0008))
TempLine<-lm(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~Coef.Lineage_AsFactorbLR, data=F0_Meta_F2_4Models)
abline(TempLine, lwd=3, col="black")

#Subsetting data to require FDR<0.10 for bHR vs. bLR in either F0 or Late Gen Meta-Analysis or p<0.05 in both & a consistent direction of effect:
TempData<-F0_Meta_F2_4Models[((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE) | (F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR<0.05 & is.na(F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR)==FALSE & F0_Meta_F2_4Models$pval<0.05 & is.na(F0_Meta_F2_4Models$pval)==FALSE) & (F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR*F0_Meta_F2_4Models$estimate)>0),]

points(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~Coef.Lineage_AsFactorbLR, data=TempData, col="red")
TempLine2<-lm(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore~Coef.Lineage_AsFactorbLR, data=TempData)
abline(TempLine2, lwd=3, col="red")

dev.off()

summary.lm(TempLine2)
# Call:
#   lm(formula = M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore ~ 
#        Coef.Lineage_AsFactorbLR, data = TempData)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -9.472e-04 -2.768e-05  1.830e-06  2.778e-05  8.333e-04 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              -1.056e-06  2.976e-06  -0.355    0.723    
# Coef.Lineage_AsFactorbLR -4.011e-05  4.961e-06  -8.084 1.73e-15 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 9.568e-05 on 1043 degrees of freedom
# (18 observations deleted due to missingness)
# Multiple R-squared:  0.05896,	Adjusted R-squared:  0.05806 
# F-statistic: 65.35 on 1 and 1043 DF,  p-value: 1.726e-15



#Which genes pop up in F0 and Meta?

#How many genes are in either dataset?
sum((is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (is.na(F0_Meta_F2_4Models$BH)==FALSE))
#[1] 13786

#How many genes are in both datasets?
sum((is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) & (is.na(F0_Meta_F2_4Models$BH)==FALSE))
#[1] 11175

#... and sig in the F0
sum((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) & (is.na(F0_Meta_F2_4Models$BH)==FALSE))
#[1] 137

#... and sig in the Late Gen Meta-Analysis
sum((is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) & (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE))
#[1] 767

#... sig in either (FDR<0.10):
sum((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE))
#[1] 922

#... sig in either (FDR<0.10) or p<0.05 in both
sum((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE) | (F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR<0.05 & is.na(F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR)==FALSE & F0_Meta_F2_4Models$pval<0.05 & is.na(F0_Meta_F2_4Models$pval)==FALSE))
#[1] 1080

1080/13786
#[1] 0.07834035

#... sig in either (FDR<0.10) or p<0.05 in both with consistent direction of effect:
sum((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE) | (F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR<0.05 & is.na(F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR)==FALSE & F0_Meta_F2_4Models$pval<0.05 & is.na(F0_Meta_F2_4Models$pval)==FALSE & (F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR*F0_Meta_F2_4Models$estimate)>0))
#[1] 1063

#Huda was asking for FDR calculations and volcano plots for just this subset:

F0_Meta_F2_4Models_LineageGenes<-F0_Meta_F2_4Models[((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE) | (F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR<0.05 & is.na(F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR)==FALSE & F0_Meta_F2_4Models$pval<0.05 & is.na(F0_Meta_F2_4Models$pval)==FALSE & (F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR*F0_Meta_F2_4Models$estimate)>0)),]

#outputting a correlation matrix just for these genes:
write.csv(cor(cbind(F0_Meta_F2_4Models_LineageGenes$Coef.Lineage_AsFactorbLR, F0_Meta_F2_4Models_LineageGenes$estimate, F0_Meta_F2_4Models_LineageGenes$M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore, F0_Meta_F2_4Models_LineageGenes$Coef.EPM_DistanceTraveled_noOutlier, F0_Meta_F2_4Models_LineageGenes$Coef.EPM_Time_Immobile, F0_Meta_F2_4Models_LineageGenes$Coef.EPM_Percent_Time_Open_Arm, F0_Meta_F2_4Models_LineageGenes$Coef.PCA_Index_Days6and.7), use="pairwise.complete.obs"), "CorMatrix_Log2FC_Behavior_OnlyHRLRGenes.csv")

#enrichment stats:
sum(F0_Meta_F2_4Models_LineageGenes$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore<0.05 & is.na(F0_Meta_F2_4Models_LineageGenes$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE)
#[1] 103

sum(F0_Meta_F2_4Models_LineageGenes$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore>0.05 & is.na(F0_Meta_F2_4Models_LineageGenes$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE)
#[1] 942

sum((is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (is.na(F0_Meta_F2_4Models$BH)==FALSE))
#[1] 13786

sum(((is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (is.na(F0_Meta_F2_4Models$BH)==FALSE)) & (F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore<0.05 & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE))
#[1] 504

sum(((is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (is.na(F0_Meta_F2_4Models$BH)==FALSE)) & (F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore>0.05 & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE))
#[1] 12835

CrossTable_NominalLocoScoreVsLineageGene<-rbind(c(103, 942), c((504-103), (12835-942)))
fisher.test(CrossTable_NominalLocoScoreVsLineageGene)

# Fisher's Exact Test for Count Data
# 
# data:  CrossTable_NominalLocoScoreVsLineageGene
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 2.559348 4.079357
# sample estimates:
# odds ratio 
# 3.242448 

sum(F0_Meta_F2_4Models_LineageGenes$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore<0.005 & is.na(F0_Meta_F2_4Models_LineageGenes$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE)
#[1] 25

sum(F0_Meta_F2_4Models_LineageGenes$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore>0.005 & is.na(F0_Meta_F2_4Models_LineageGenes$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE)
#[1] 1020

sum(((is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (is.na(F0_Meta_F2_4Models$BH)==FALSE)) & (F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore<0.005 & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE))
#[1] 61

sum(((is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (is.na(F0_Meta_F2_4Models$BH)==FALSE)) & (F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore>0.005 & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE))
#[1] 13278

CrossTable_Nominal005LocoScoreVsLineageGene<-rbind(c(25, 1020), c((61-25), (13278-1020)))
fisher.test(CrossTable_Nominal005LocoScoreVsLineageGene)

# Fisher's Exact Test for Count Data
# 
# data:  CrossTable_Nominal005LocoScoreVsLineageGene
# p-value = 9.618e-13
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 4.78011 14.35111
# sample estimates:
# odds ratio 
# 8.342773 

###########################################


#MH (2022-05-09): I didn't update any of the code or output after this point to reflect the inclusion of Time Immobile, Distance Traveled (outlier removed), or PavCA Index (Days 6 & 7) because the output after this point either A) Doesn't change, B) Changes, but isn't planned to be used in the paper (e.g., the manhattan plot for the p-values across all F2 RNA-Seq variables.

#Oh wait - correction: I jumped ahead and outputted the annotated version of the data.frame too.


###########################################


#... sig in both F0 and Late Gen Meta-Analysis (fdr<0.10)
F0_Meta_F2_4Models$Symbol[(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) & (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE)]

# [1] "Nkd2"     "Tnnt1"    "Pex11a"   "Plekhb1"  "Gal"      "Aldh3b1"  "Ctsf"     "Aldh1a1"  "Ifit1"    "Mrpl43"  
# [11] "Egfem1"   "Tmem144"  "Fcrl2"    "Dbt"      "Prss12"   "Tes"      "Asb15"    "Tmem176a" "Cd4"      "Rps4y2"  
# [21] "Rexo4"    "Sp3"      "Bloc1s6"  "Slx4ip"   "Pxmp4"    "Ggt7"     "Gss"      "Eif6"     "Mcmdc2"   "Nbn"     
# [31] "Hook1"    "Cyp2j4"   "C1qb"     "C1qc"     "C1qa"     "Smoc1"    "Utp20"    "Kif15"    "Exosc7"   "Oard1"   
# [41] "Nudt12"   "Alb"      "Ddc"      "Rnasel"   "Mettl16"  "Pnpo"     "Ghdc"     "Ephx2"    "Msra"     "Nudt18"  
# [51] "Rgcc"     "Gpr183"   "Masp1"    "Hbegf"    "Ednra"    "Txnl4b"   "Vps9d1"   "Ccdc167"  "Pkib"     "Rtn4ip1" 
# [61] "Grifin"   "Psmg3" 

#some familiar faces.

62/11175
#[1] 0.005548098


#using stricter FDR
F0_Meta_F2_4Models$Symbol[(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.05 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) & (F0_Meta_F2_4Models$BH<0.05 & is.na(F0_Meta_F2_4Models$BH)==FALSE)]

# [1] "Nkd2"     "Plekhb1"  "Gal"      "Aldh3b1"  "Aldh1a1"  "Ifit1"    "Mrpl43"   "Egfem1"   "Tmem144"  "Fcrl2"   
# [11] "Dbt"      "Prss12"   "Tes"      "Asb15"    "Tmem176a" "Rps4y2"   "Bloc1s6"  "Slx4ip"   "Pxmp4"    "Gss"     
# [21] "Mcmdc2"   "Cyp2j4"   "C1qb"     "C1qc"     "C1qa"     "Smoc1"    "Kif15"    "Exosc7"   "Nudt12"   "Alb"     
# [31] "Ddc"      "Ghdc"     "Rgcc"     "Txnl4b"   "Rtn4ip1"  "Grifin"   "Psmg3" 


########

#Now adding in locomotor activity too:


#How many genes are in either dataset & in the F2s:
sum(((is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (is.na(F0_Meta_F2_4Models$BH)==FALSE)) & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE)
#[1] 13339

#How many genes are in either dataset & in the F2s & have a nominal effect of locomotor activity (p<0.05):
sum(((is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (is.na(F0_Meta_F2_4Models$BH)==FALSE)) & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE & F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore<0.05)
#[1] 504

504/13339
#[1] 0.03778394

#How many genes are in either dataset & in the F2s & have a strong nominal effect of locomotor activity (p<0.005):
sum(((is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (is.na(F0_Meta_F2_4Models$BH)==FALSE)) & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE & F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore<0.005)
#[1] 61

61/13339
#[1] 0.004573056

#How many genes are in both datasets & in the F2s:
sum(((is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) & (is.na(F0_Meta_F2_4Models$BH)==FALSE)) & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE)
#[1] 11055

#How many genes are in both datasets & in the F2s & have a nominal effect of locomotor activity:
sum(((is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) & (is.na(F0_Meta_F2_4Models$BH)==FALSE)) & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE & F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore<0.05)
#[1] 384
384/11055

#[1] 0.03473541

#How many genes are in both datasets & in the F2s & have a strong nominal effect of locomotor activity:
sum(((is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) & (is.na(F0_Meta_F2_4Models$BH)==FALSE)) & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE & F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore<0.005)
#[1] 47

47/11055

#[1] 0.00425147

#... sig in either bHR/bLR dataset (FDR<0.10) & F2 locomotor p<0.05:

sum(((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE)) & (F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore<0.05 & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE))
#[1] 86

F0_Meta_F2_4Models$Symbol[((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE)) & (F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore<0.05 & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE)]

# [1] "Zfp551"     "Clasrp"     "Lsr"        NA           "Mfge8"      "Pex11a"     "Unc45a"     "Lipt2"     
# [9] "C2cd3"      "Ucp2"       "Plekhb1"    NA           "Kcnk4"      "Aldh1a1"    "Ifit2"      "Poc5"      
# [17] "RGD1564606" "RGD1359508" "Neurog2"    "Akr1b1"     "Rarres2"    "Gimap5"     "Tmem176a"   "Fam221a"   
# [25] "Abcg2"      "Tex261"     "Txnrd3"     "Fkbp4"      "Sp3"        "Ttc30a1"    "Commd9"     "Slx4ip"    
# [33] "Mmp24"      "Snhg11"     "Slc2a10"    "Pnisr"      "Serinc2"    "Pramef8"    "Sypl1"      "Ilvbl"     
# [41] "Akap8"      "Akap8l"     "Cyp4f4"     "Nudt4"      "Fzd6"       "Sun2"       "Robo3"      "Acad11"    
# [49] "Tmie"       "Cmc1"       NA           NA           "Rev1"       "Gls"        "Idh1"       "Art3"      
# [57] "Grxcr1"     "Mxd4"       "Ddc"        "Vstm2a"     "Cfh"        "Itln1"      "Wsb1"       "Kcnh6"     
# [65] "Hnrnpc"     "Prss55"     "Cdhr2"      "LOC690414"  "Bphl"       "Nqo2"       "Mkx"        "Fut10"     
# [73] "LOC302192"  NA           "Cyyr1"      "Cd200r1"    "Atp13a5"    "Dsc2"       "Camk4"      "Ist1"      
# [81] "Spg7"       "Vps9d1"     "Elfn1"      "Rasa4"      "Sfswap"     "Ift81"



#... sig in either (FDR<0.10) or p<0.05 in both and found in F2 (for comparison)
sum(((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE) | (F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR<0.05 & is.na(F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR)==FALSE & F0_Meta_F2_4Models$pval<0.05 & is.na(F0_Meta_F2_4Models$pval)==FALSE)) & (is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE))
#[1] 1062

#... sig in either (FDR<0.10) or p<0.05 in both & F2 locomotor p<0.05:
sum(((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE) | (F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR<0.05 & is.na(F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR)==FALSE & F0_Meta_F2_4Models$pval<0.05 & is.na(F0_Meta_F2_4Models$pval)==FALSE)) & (F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore<0.05 & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE))
#[1] 104

F0_Meta_F2_4Models$Symbol[((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE) | (F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR<0.05 & is.na(F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR)==FALSE & F0_Meta_F2_4Models$pval<0.05 & is.na(F0_Meta_F2_4Models$pval)==FALSE)) & (F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore<0.05 & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE)]

# [1] "Zfp551"     "Clasrp"     "Lsr"        NA           "Mfge8"      "Fanci"      "Pex11a"     "Wdr93"     
# [9] "Unc45a"     "Lipt2"      "C2cd3"      "Ucp2"       "Plekhb1"    "Mapk1ip1"   NA           "Kcnk4"     
# [17] "Aldh1a1"    "Ifit2"      "Poc5"       "RGD1564606" "Hmgcs1"     "RGD1359508" "Neurog2"    "Ppa2"      
# [25] "Akr1b1"     "Ptn"        "Rarres2"    "Gimap5"     "Tmem176a"   "Fam221a"    "Abcg2"      "Tex261"    
# [33] "Txnrd3"     "Fkbp4"      "Sp3"        "Ttc30a1"    "Commd9"     "Tmco5a"     "Slx4ip"     "Mmp24"     
# [41] "Snhg11"     "Slc2a10"    "Rgs20"      "Pnisr"      "Serinc2"    "Pramef8"    "Sypl1"      "Ilvbl"     
# [49] "Akap8"      "Akap8l"     "Cyp4f4"     "Nudt4"      "Rpl30"      "Fzd6"       "Sun2"       "Pus3"      
# [57] "Robo3"      "Plscr4"     "Acad11"     "Tmie"       "Cmc1"       NA           NA           "Rftn1"     
# [65] "Rev1"       "Mstn"       "Gls"        "Idh1"       "Art3"       "Grxcr1"     "Mxd4"       "Ddc"       
# [73] "Vstm2a"     "Cfh"        "Itln1"      "Wsb1"       "Kcnh6"      "Hnrnpc"     "Prss55"     "Cdhr2"     
# [81] "LOC690414"  "Bphl"       "Nqo2"       "Mkx"        "Fut10"      "LOC302192"  NA           "Cyyr1"     
# [89] "Cd200r1"    "Atp13a5"    "Ece2"       "Rock1"      "Dsc2"       "Camk4"      "Cndp1"      "Ist1"      
# [97] "Spg7"       "Vps9d1"     "Afg3l1"     "Elfn1"      "Rasa4"      "Sfswap"     "P2rx4"      "Ift81"


104/13339
#[1] 0.007796686

# % with either fdr<0.10 or p<0.05 in both bHR/bLR
1062/13339
#[1] 0.07961616

#% with nominal effect of locomotor activity
504/13339
#[1] 0.03778394

0.03778394*0.07961616
#[1] 0.003008212
#expected

#enrichment:
0.007796686/0.003008212
#[1] 2.591801

cbind(c(104,(1062-104)), c((504-104), (13339-504-1062+104)))

#             #bHR/bLR yes  #bHR/bLR no
#              [,1]         [,2]
# F2 Yes [1,]  104          400
# F2 No [2,]  958           11877

#double check:
104+400+958+11877
#[1] 13339

TempContingencyMatrix<-cbind(c(104,(1062-104)), c((504-104), (13339-504-1062+104)))

fisher.test(TempContingencyMatrix)
# Fisher's Exact Test for Count Data
# 
# data:  TempContingencyMatrix
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 2.546499 4.051404
# sample estimates:
# odds ratio 
# 3.22296 


#... sig in either (FDR<0.10) or p<0.05 in both & F2 locomotor p<0.005:
sum(((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE) | (F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR<0.05 & is.na(F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR)==FALSE & F0_Meta_F2_4Models$pval<0.05 & is.na(F0_Meta_F2_4Models$pval)==FALSE)) & (F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore<0.005 & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE))
#[1] 26

F0_Meta_F2_4Models$Symbol[((F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR)==FALSE) | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE) | (F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR<0.05 & is.na(F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR)==FALSE & F0_Meta_F2_4Models$pval<0.05 & is.na(F0_Meta_F2_4Models$pval)==FALSE)) & (F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore<0.005 & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE)]
# [1] NA           "Mfge8"      "Fanci"      "Wdr93"      "Ucp2"       NA           "RGD1359508" "Rarres2"   
# [9] "Gimap5"     "Tmem176a"   "Ttc30a1"    "Serinc2"    "Sypl1"      "Rpl30"      "Fzd6"       "Rftn1"     
# [17] "Mstn"       "Idh1"       "Art3"       "Bphl"       "Nqo2"       "Dsc2"       "Spg7"       "Vps9d1"    
# [25] "Afg3l1"     "Elfn1" 

26/13339
#[1] 0.001949172

# % with either fdr<0.10 or p<0.05 in both bHR/bLR
1062/13339
#[1] 0.07961616

#% with strong nominal effect (p<0.005) of locomotor activity
61/13339
#[1] 0.004573056

0.004573056*0.07961616
#[1] 0.0003640892
#expected

#enrichment:
0.001949172/0.0003640892
#[1] 5.353556


cbind(c(26,(1062-26)), c((61-26), (13339-61-1062+26)))

#             #bHR/bLR yes  #bHR/bLR no
#              [,1]         [,2]
# F2 Yes [1,]  26            35
# F2 No [2,]  1036           12242

#double check:
26+35+1036+12242
#[1] 13339

TempContingencyMatrix<-cbind(c(26,(1062-26)), c((61-26), (13339-61-1062+26)))

fisher.test(TempContingencyMatrix)

# Fisher's Exact Test for Count Data
# 
# data:  TempContingencyMatrix
# p-value = 1.604e-13
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 5.050814 15.069883
# sample estimates:
# odds ratio 
# 8.77503 


#Let's write that out - seems useful:

write.csv(F0_Meta_F2_4Models[(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE)) & F0_Meta_F2_4Models$M8_CovSexRNAMetricsTechnical_P.value.Total_LocoScore<0.05,], "F0_Meta_F2_4Models_wHRLRFDR10_andM8LocomotorP05.csv")

write.csv(F0_Meta_F2_4Models[(F0_Meta_F2_4Models$P.value.adj.Lineage_AsFactorbLR<0.10 | (F0_Meta_F2_4Models$BH<0.10 & is.na(F0_Meta_F2_4Models$BH)==FALSE)) & F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore<0.05,], "F0_Meta_F2_4Models_wHRLRFDR10_andM9LocomotorP05.csv")

setwd("~/Documents/Microarray Gen/HRLR/SecondaryAnalyses/RatGenomeDatabase_eQTL")


