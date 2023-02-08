#F2 HC RNA-Seq Dataset
#15_Limma Differential Expression Analyses
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.

#############################################################1/10/2022

#Differential Expression Analysis: Total LocoScore

#Model based on top correlates with PCs, stepwise regression for PCs

design<-model.matrix(~Total_LocoScore+Percent.Intergenic+RibosomePerc+SequencingBatch+Percent.Intronic+Sex+Percent.UTR+RNAconc+Dissector_AsNumeric+STGT_experience, data=F2_PCA_wMetaData)
design

v<-voom(F2_RawHitCount_dge_noLowHits_TMM, design, plot=TRUE)
v
# 
# An object of class "EList"

vfit<-lmFit(v,design)

str(vfit)
#This is the linear regression output (raw)
# 
# Formal class 'MArrayLM' [package "limma"] with 1 slot

efit<-eBayes(vfit)
#This is the linear regression output with an empirical Bayes correction

png("F2_TotalLocoScore_wPercent.UTR_20220110_MeanVarianceTrend.png")
plotSA(efit, main="F2_TotalLocoScore_wPercent.UTR_20220110: Mean-Variance Trend")
dev.off()
#Heteroskedasticity now gone! How pretty!

write.fit(efit, adjust="BH", file="F2_Limma_Results_TotalLocoScore_wPercent.UTR_20220110.txt")

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.05)
#How many results meet an FDR<0.05?

summary(dt) 
# (Intercept) Total_LocoScore Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 Percent.Intronic SexMale Percent.UTR RNAconc Dissector_AsNumeric STGT_experienceTRUE
# Down           265               0               6362         5842             2964             2424             3704     627        3823    1145                 962                 446
# NotSig        2191           14056               1692         3985             7627             8977             8142   12831        4665   12804               12598               13220
# Up           11600               0               6002         4229             3465             2655             2210     598        5568     107                 496                 390

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.10)
summary(dt)

#         (Intercept) Total_LocoScore Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 Percent.Intronic SexMale Percent.UTR RNAconc Dissector_AsNumeric STGT_experienceTRUE
# Down           349               0               6496         6230             3400             2863             4255     959        4057    1486                1505                 742
# NotSig        1892           14056               1445         3360             6604             7982             7157   12196        3992   12332               11762               12565
# Up           11815               0               6115         4466             4052             3211             2644     901        6007     238                 789                 749

#Conclusions:
#Our two RNA metric co-variates really matter - %Intergenic and %Ribosome and %UTR are related to the majority of genes in our dataset, %intronic is related to half the genes in the dataset
#RNAconc is related to a smaller number of genes (1600)
#Our technical co-variates are also well indicated - sequencing batch is related to at least half of the genes in our dataset (maybe more), Dissector is related to >2000 genes, and ST/GT experience and Sex are both related to around ~1500-1800 genes.
#... and locomotor score is related to *nothing*  :(

#########################

#Same differential expression model for total loco score, but without percent.UTR because it correlated with Lineage in the F0s

design<-model.matrix(~Total_LocoScore+Percent.Intergenic+RibosomePerc+SequencingBatch+Percent.Intronic+Sex+RNAconc+Dissector_AsNumeric+STGT_experience, data=F2_PCA_wMetaData)
design


v<-voom(F2_RawHitCount_dge_noLowHits_TMM, design, plot=TRUE)
v
# 
# An object of class "EList"

vfit<-lmFit(v,design)

str(vfit)
#This is the linear regression output (raw)
# 
# Formal class 'MArrayLM' [package "limma"] with 1 slot

efit<-eBayes(vfit)
#This is the linear regression output with an empirical Bayes correction

png("F2_TotalLocoScore_20220110_MeanVarianceTrend.png")
plotSA(efit, main="F2_TotalLocoScore_20220110: Mean-Variance Trend")
dev.off()
#Heteroskedasticity now gone! How pretty!

write.fit(efit, adjust="BH", file="F2_Limma_Results_TotalLocoScore_20220110.txt")

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.05)
#How many results meet an FDR<0.05?

summary(dt) 
#    (Intercept) Total_LocoScore Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 Percent.Intronic SexMale RNAconc Dissector_AsNumeric STGT_experienceTRUE
# Down           358               0               6322         5902             1535             2184             3729     596    1289                 777                 345
# NotSig        1132           14056               1781         3524            10626             9362             8128   12915   12371               12816               13224
# Up           12566               0               5953         4630             1895             2510             2199     545     396                 463                 487

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.10)

summary(dt) 
#         (Intercept) Total_LocoScore Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 Percent.Intronic SexMale RNAconc Dissector_AsNumeric STGT_experienceTRUE
# Down           417               0               6447         6223             2057             2668         4266        937    1811                1329                 601
# NotSig         953           14056               1523         2967             9536             8387         7144       12291   11554               11941               12620
# Up           12686               0               6086         4866             2463             3001         2646       828     691                 786                 835
# 


#Conclusions:
#Our two RNA metric co-variates really matter - %Intergenic and %Ribosome are related to the majority of genes in our dataset, %intronic is related to half the genes in the dataset
#RNAconc is doing a little less, but still realted to 2500 genes
#Our technical co-variates are also well indicated - sequencing batch is related to at least a third of the genes in our dataset (maybe more), Dissector is related to ~2000 genes, and ST/GT experience and Sex are both related to around ~1400-1700 genes.
#... and locomotor score is related to *nothing*  :(

#########################################This is the F2 Locoscore model we are using--M9  2/28/2022

# A model paralleling the F0 model - covariates that are strong correlates of PC1/PC2 in both datasets, and additional technical variables as co-variates that are strongly supported by theory.

design<-model.matrix(~Total_LocoScore+Percent.Intergenic+RibosomePerc+SequencingBatch+Sex+Dissector_AsNumeric+STGT_experience, data=F2_PCA_wMetaData)
design

v<-voom(F2_RawHitCount_dge_noLowHits_TMM, design, plot=TRUE)
v
# 
# An object of class "EList"

vfit<-lmFit(v,design)

str(vfit)
#This is the linear regression output (raw)
# 
# Formal class 'MArrayLM' [package "limma"] with 1 slot

efit<-eBayes(vfit)
#This is the linear regression output with an empirical Bayes correction

png("F2_TotalLocoScore_SimplerModel_20220110_MeanVarianceTrend.png")
plotSA(efit, main="F2_TotalLocoScore_SimplerModel_20220110: Mean-Variance Trend")
dev.off()
#Heteroskedasticity now gone! How pretty!

write.fit(efit, adjust="BH", file="F2_Limma_Results_TotalLocoScore_SimplerModel_20220110.txt")

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.05)
#How many results meet an FDR<0.05?

summary(dt) 
#      (Intercept) Total_LocoScore Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric STGT_experienceTRUE
# Down           322               0               6568         6112             3257             2790     581                1025                 453
# NotSig        1076           14056               1608         3386             8372             8396   12919               12260               13117
# Up           12658               0               5880         4558             2427             2870     556                 771                 486

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.10)

summary(dt) 
#         (Intercept) Total_LocoScore Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric STGT_experienceTRUE
# Down           374               0               6707         6432             3742             3298     894                1525                 769
# NotSig         898           14056               1343         2875             7340             7380   12312               11335               12433
# Up           12784               0               6006         4749             2974             3378     850                1196                 854
# 

#Conclusions: 
#Our RNA metric co-variates really matter - %Intergenic and %Ribosome are related to the majority of genes in our dataset
#Our technical co-variates are also well indicated - sequencing batch is related to at least a third of the genes in our dataset (maybe more), Dissector is related to ~2500 genes, and ST/GT experience and Sex are both related to around 1600 genes.
#... and locomotor score is related to *nothing*  :(

##############################Adding Dissection Date Collapsed

# A model paralleling the F0 model - covariates that are strong correlates of PC1/PC2 in both datasets, and additional technical variables as co-variates that are strongly supported by theory.

design<-model.matrix(~Total_LocoScore+Percent.Intergenic+RibosomePerc+SequencingBatch+Sex+Dissector_AsNumeric+STGT_experience+HPC_Dissection_Day, data=F2_PCA_wMetaData)
design

v<-voom(F2_RawHitCount_dge_noLowHits_TMM, design, plot=TRUE)
v
# 
# An object of class "EList"

vfit<-lmFit(v,design)

str(vfit)
#This is the linear regression output (raw)
# 
# Formal class 'MArrayLM' [package "limma"] with 1 slot

efit<-eBayes(vfit)
#This is the linear regression output with an empirical Bayes correction

png("F2_TotalLocoScore_SimplerModel_20220121_WithHPC_Dissection_Day_MeanVarianceTrend.png")
plotSA(efit, main="F2_TotalLocoScore_SimplerModel_20220121_WithHPC_Dissection_Day: Mean-Variance Trend")
dev.off()
#Heteroskedasticity now gone! How pretty!

write.fit(efit, adjust="BH", file="F2_Limma_Results_TotalLocoScore_SimplerModel_20220121_WithHPC_Dissection_Day.txt")

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.05)
#How many results meet an FDR<0.05?

summary(dt) 
#           (Intercept) Total_LocoScore Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric STGT_experienceTRUE HPC_Dissection_Day
# Down           262               0               6558         6121             1070             1434     565                 829                 342                428
# NotSig        1116           14056               1626         3370            11621            10810   12941               12612               13413              13483
# Up           12678               0               5872         4565             1365             1812     550                 615                 301                145

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.10)

summary(dt) 
#             (Intercept) Total_LocoScore Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric STGT_experienceTRUE HPC_Dissection_Day
# Down           312               0               6679         6439             1479             1858     886                1358                 624                856
# NotSig         940           14056               1391         2852            10818             9727   12317               11688               12825              12834
# Up           12804               0               5986         4765             1759             2471     853                1010                 607                366
# 


#####EPM variables now

# A model paralleling the F0 model - covariates that are strong correlates of PC1/PC2 in both datasets, and additional technical variables as co-variates that are strongly supported by theory.

design<-model.matrix(~EPM_Percent_Time_Open_Arm+Percent.Intergenic+RibosomePerc+SequencingBatch+Sex+Dissector_AsNumeric+STGT_experience, data=F2_PCA_wMetaData)
design

v<-voom(F2_RawHitCount_dge_noLowHits_TMM, design, plot=TRUE)
v
# 
# An object of class "EList"

vfit<-lmFit(v,design)

str(vfit)
#This is the linear regression output (raw)
# 
# Formal class 'MArrayLM' [package "limma"] with 1 slot

efit<-eBayes(vfit)
#This is the linear regression output with an empirical Bayes correction

png("F2_EPM_PercOA_SimplerModel_20220121_MeanVarianceTrend.png")
plotSA(efit, main="F2_EPM_PercOA_SimplerModel_20220121: Mean-Variance Trend")
dev.off()
#Heteroskedasticity now gone! How pretty!

write.fit(efit, adjust="BH", file="F2_Limma_Results_EPM_PercOA_SimplerModel_20220121.txt")

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.05)
#How many results meet an FDR<0.05?

summary(dt) 
#          (Intercept) EPM_Percent_Time_Open_Arm Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric STGT_experienceTRUE
# Down           335                         0               6561         6101             3304             2835     484                1026                 457
# NotSig        1069                     14056               1603         3393             8288             8330   13057               12261               13101
# Up           12652                         0               5892         4562             2464             2891     515                 769                 498

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.10)

summary(dt) 
#            (Intercept) EPM_Percent_Time_Open_Arm Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric STGT_experienceTRUE
# Down           378                         0               6706         6429             3817             3337     815                1527                 766
# NotSig         912                     14056               1332         2872             7220             7313   12409               11331               12431
# Up           12766                         0               6018         4755             3019             3406     832                1198                 859
# 



# A model paralleling the F0 model - covariates that are strong correlates of PC1/PC2 in both datasets, and additional technical variables as co-variates that are strongly supported by theory.

design<-model.matrix(~EPM_Time_Immobile+Percent.Intergenic+RibosomePerc+SequencingBatch+Sex+Dissector_AsNumeric+STGT_experience, data=F2_PCA_wMetaData)
design

v<-voom(F2_RawHitCount_dge_noLowHits_TMM, design, plot=TRUE)
v
# 
# An object of class "EList"

vfit<-lmFit(v,design)

str(vfit)
#This is the linear regression output (raw)
# 
# Formal class 'MArrayLM' [package "limma"] with 1 slot

efit<-eBayes(vfit)
#This is the linear regression output with an empirical Bayes correction

png("F2_EPM_Time_Immobile_SimplerModel_20220121_MeanVarianceTrend.png")
plotSA(efit, main="F2_EPM_Time_Immobile_SimplerModel_20220121: Mean-Variance Trend")
dev.off()
#Heteroskedasticity now gone! How pretty!

write.fit(efit, adjust="BH", file="F2_Limma_Results_EPM_Time_Immobile_SimplerModel_20220121.txt")

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.05)
#How many results meet an FDR<0.05?

summary(dt) 
#           (Intercept) EPM_Time_Immobile Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric STGT_experienceTRUE
# Down           310                 0               6571         6113             3342             2801     361                1050                 472
# NotSig        1096             14056               1605         3377             8222             8395   13294               12221               13090
# Up           12650                 0               5880         4566             2492             2860     401                 785                 494

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.10)


summary(dt) 
#            (Intercept) EPM_Time_Immobile Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric STGT_experienceTRUE
# Down           357                 0               6713         6433             3855             3301     640                1554                 777
# NotSig         924             14056               1338         2863             7163             7375   12699               11279               12431
# Up           12775                 0               6005         4760             3038             3380     717                1223                 848
# 



# A model paralleling the F0 model - covariates that are strong correlates of PC1/PC2 in both datasets, and additional technical variables as co-variates that are strongly supported by theory.

design<-model.matrix(~EPM_DistanceTraveled+Percent.Intergenic+RibosomePerc+SequencingBatch+Sex+Dissector_AsNumeric+STGT_experience, data=F2_PCA_wMetaData)
design

v<-voom(F2_RawHitCount_dge_noLowHits_TMM, design, plot=TRUE)
v
# 
# An object of class "EList"

vfit<-lmFit(v,design)

str(vfit)
#This is the linear regression output (raw)
# 
# Formal class 'MArrayLM' [package "limma"] with 1 slot

efit<-eBayes(vfit)
#This is the linear regression output with an empirical Bayes correction

png("F2_EPM_DistanceTraveled_SimplerModel_20220121_MeanVarianceTrend.png")
plotSA(efit, main="F2_EPM_DistanceTraveled_SimplerModel_20220121: Mean-Variance Trend")
dev.off()
#Heteroskedasticity now gone! How pretty!

write.fit(efit, adjust="BH", file="F2_Limma_Results_EPM_DistanceTraveled_SimplerModel_20220121.txt")

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.05)
#How many results meet an FDR<0.05?

summary(dt) 
#           (Intercept) EPM_DistanceTraveled Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric STGT_experienceTRUE
# Down           317                    0               6566         6115             3337             2810     378                1027                 447
# NotSig        1100                14056               1610         3406             8226             8372   13273               12249               13175
# Up           12639                    0               5880         4535             2493             2874     405                 780                 434

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.10)


summary(dt) 
#            (Intercept) EPM_DistanceTraveled Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric STGT_experienceTRUE
# Down           366                    0               6704         6433             3850             3318     677                1535                 733
# NotSig         919                14056               1347         2883             7159             7342   12711               11319               12551
# Up           12771                    0               6005         4740             3047             3396     668                1202                 772

####################

#After looking at the distributions for the data, we decided to come back and re-run distance traveled with the extreme outlier subject removed

hist(F2_PCA_wMetaData$EPM_DistanceTraveled)

summary(F2_PCA_wMetaData$EPM_DistanceTraveled)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1104    1736    2075    2075    2363    5209 

#Outlier has >5000, next highest is <3500

design<-model.matrix(~EPM_DistanceTraveled+Percent.Intergenic+RibosomePerc+SequencingBatch+Sex+Dissector_AsNumeric+STGT_experience, data=F2_PCA_wMetaData[F2_PCA_wMetaData$EPM_DistanceTraveled<5000,])
design

dim(design)
#[1] 244   9

v<-voom(F2_RawHitCount_dge_noLowHits_TMM[,F2_PCA_wMetaData$EPM_DistanceTraveled<5000], design, plot=TRUE)
v
# 
# An object of class "EList"

vfit<-lmFit(v,design)

str(vfit)
#This is the linear regression output (raw)
# 
# Formal class 'MArrayLM' [package "limma"] with 1 slot

efit<-eBayes(vfit)
#This is the linear regression output with an empirical Bayes correction

png("F2_EPM_DistanceTravelednoOutlier_SimplerModel_20220121_MeanVarianceTrend.png")
plotSA(efit, main="F2_EPM_DistanceTravelednoOutlier_SimplerModel_20220121: Mean-Variance Trend")
dev.off()
#Heteroskedasticity now gone! How pretty!

write.fit(efit, adjust="BH", file="F2_Limma_Results_EPM_DistanceTravelednoOutlier_SimplerModel_20220121.txt")

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.05)
#How many results meet an FDR<0.05?

summary(dt) 
# (Intercept) EPM_DistanceTraveled Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric STGT_experienceTRUE
# Down           320                    0               6563         6080             3329             2815     358                1023                 452
# NotSig        1111                14056               1610         3441             8232             8381   13315               12254               13156
# Up           12625                    0               5883         4535             2495             2860     383                 779                 448

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.10)

summary(dt) 
# (Intercept) EPM_DistanceTraveled Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric STGT_experienceTRUE
# Down           370                    0               6703         6417             3853             3299     609                1522                 736
# NotSig         931                14056               1348         2912             7150             7368   12815               11332               12535
# Up           12755                    0               6005         4727             3053             3389     632                1202                 785


#######ST/GT here

# A model paralleling the F0 model - covariates that are strong correlates of PC1/PC2 in both datasets, and additional technical variables as co-variates that are strongly supported by theory.

design<-model.matrix(~LearningClassification_AsFactor+Percent.Intergenic+RibosomePerc+SequencingBatch+Sex+Dissector_AsNumeric, data=F2_PCA_wMetaData)
design

v<-voom(F2_RawHitCount_dge_noLowHits_TMM, design, plot=TRUE)
v
# 
# An object of class "EList"

vfit<-lmFit(v,design)

str(vfit)
#This is the linear regression output (raw)
# 
# Formal class 'MArrayLM' [package "limma"] with 1 slot

efit<-eBayes(vfit)
#This is the linear regression output with an empirical Bayes correction

png("F2_LearningClassification_SimplerModel_20220121_MeanVarianceTrend.png")
plotSA(efit, main="F2_LearningClassification_SimplerModel_20220121: Mean-Variance Trend")
dev.off()
#Heteroskedasticity now gone! How pretty!

write.fit(efit, adjust="BH", file="F2_Limma_Results_LearningClassification_SimplerModel_20220121.txt")

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.05)
#How many results meet an FDR<0.05?

summary(dt) 
#             (Intercept) LearningClassification_AsFactor LearningClassification_AsFactorGT LearningClassification_AsFactorIN Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric
# Down           299                             466                                 0                                 0               6559         6074             3353             2785     491                 964
# NotSig        1059                           13138                             14056                             14056               1619         3430             8212             8450   13090               12341
# Up           12698                             452                                 0                                 0               5878         4552             2491             2821     475                 751

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.10)


summary(dt) 
#             (Intercept) LearningClassification_AsFactor LearningClassification_AsFactorGT LearningClassification_AsFactorIN Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric
# Down           356                             799                                 0                                 0               6697         6416             3866             3288     844                1469
# NotSig         871                           12505                             14056                             14055               1352         2878             7157             7423   12407               11397
# Up           12829                             752                                 0                                 1               6007         4762             3033             3345     805                1190
# 


# A model paralleling the F0 model - covariates that are strong correlates of PC1/PC2 in both datasets, and additional technical variables as co-variates that are strongly supported by theory.

design<-model.matrix(~EPM_Time_Immobile+Percent.Intergenic+RibosomePerc+SequencingBatch+Sex+Dissector_AsNumeric+STGT_experience, data=F2_PCA_wMetaData)
design

v<-voom(F2_RawHitCount_dge_noLowHits_TMM, design, plot=TRUE)
v
# 
# An object of class "EList"

vfit<-lmFit(v,design)

str(vfit)
#This is the linear regression output (raw)
# 
# Formal class 'MArrayLM' [package "limma"] with 1 slot

efit<-eBayes(vfit)
#This is the linear regression output with an empirical Bayes correction

png("F2_EPM_Time_Immobile_SimplerModel_20220121_MeanVarianceTrend.png")
plotSA(efit, main="F2_EPM_Time_Immobile_SimplerModel_20220121: Mean-Variance Trend")
dev.off()
#Heteroskedasticity now gone! How pretty!

write.fit(efit, adjust="BH", file="F2_Limma_Results_EPM_Time_Immobile_SimplerModel_20220121.txt")

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.05)
#How many results meet an FDR<0.05?

summary(dt) 
#           (Intercept) EPM_Time_Immobile Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric STGT_experienceTRUE
# Down           310                 0               6571         6113             3342             2801     361                1050                 472
# NotSig        1096             14056               1605         3377             8222             8395   13294               12221               13090
# Up           12650                 0               5880         4566             2492             2860     401                 785                 494

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.10)


summary(dt) 
#            (Intercept) EPM_Time_Immobile Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric STGT_experienceTRUE
# Down           357                 0               6713         6433             3855             3301     640                1554                 777
# NotSig         924             14056               1338         2863             7163             7375   12699               11279               12431
# Up           12775                 0               6005         4760             3038             3380     717                1223                 848
# 

################################

#I came back and ran a version using Days6 & 7 PavCA score after reviewing the distribution for the variable - there really doesn't seem to be any evidence of bimodality, so I'm not sure it makes sense to treat this as categorical instead of continuous for our purposes. We're probably just losing power.

# A model paralleling the F0 model - covariates that are strong correlates of PC1/PC2 in both datasets, and additional technical variables as co-variates that are strongly supported by theory.

#Note: Limma was unhappy that there wasn't PavCA data for the animals with no experience with the PavCA task, so we had to subset the data in order to use PavCA score as a predictor.
#This left a sample size of 205, but made it so that we could remove STGT experience from the model (one fewer predictor)

design<-model.matrix(~PCA_Index_Days6and.7+Percent.Intergenic+RibosomePerc+SequencingBatch+Sex+Dissector_AsNumeric, data=F2_PCA_wMetaData)
design
dim(design)
#[1] 205   8

dim(F2_RawHitCount_dge_noLowHits_TMM[,F2_PCA_wMetaData$STGT_experience==TRUE])
#[1] 14056   205

v<-voom(F2_RawHitCount_dge_noLowHits_TMM[,F2_PCA_wMetaData$STGT_experience==TRUE], design, plot=TRUE)
v
# 
# An object of class "EList"

vfit<-lmFit(v,design)

str(vfit)
#This is the linear regression output (raw)
# 
# Formal class 'MArrayLM' [package "limma"] with 1 slot

efit<-eBayes(vfit)
#This is the linear regression output with an empirical Bayes correction

png("F2_PavCADays67_SimplerModel_20220121_MeanVarianceTrend.png")
plotSA(efit, main="F2_PavCADays67_SimplerModel_20220121: Mean-Variance Trend")
dev.off()
#Heteroskedasticity now gone! How pretty!

write.fit(efit, adjust="BH", file="F2_Limma_Results_PavCADays67_SimplerModel_20220121.txt")

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.05)
#How many results meet an FDR<0.05?

summary(dt) 
# 
# (Intercept) PCA_Index_Days6and.7 Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric
# Down           221                    0               6340         5922             3942             3407     263                1599
# NotSig        1019                14056               1876         3021             6942             7269   13517               10982
# Up           12816                    0               5840         5113             3172             3380     276                1475

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.10)

summary(dt) 
# (Intercept) PCA_Index_Days6and.7 Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric
# Down           261                    0               6499         6165             4332             3882     473                2098
# NotSig         882                14056               1565         2553             6100             6277   13109                9842
# Up           12913                    0               5992         5338             3624             3897     474                2116


#######################

#Differential expression for Total LocoScore, using the simplest model (just Sex and the top correlate with PC1 as co-variates):


#How about locomotor activity, sex, Percent intergenic with Model 6?

colnames(F2_PCA_wMetaData)

design<-model.matrix(~Total_LocoScore +  Sex_AsFactor + Percent.Intergenic, data=F2_PCA_wMetaData)
design


v<-voom(F2_RawHitCount_dge_noLowHits_TMM, design, plot=TRUE)
v
# 
# An object of class "EList"

vfit<-lmFit(v,design)

str(vfit)
#This is the linear regression output (raw)
# 
# Formal class 'MArrayLM' [package "limma"] with 1 slot

efit<-eBayes(vfit)
#This is the linear regression output with an empirical Bayes correction

png("F2_LocoScore_Sex_Intergenic_M6_MeanVarianceTrend.png")
plotSA(efit, main="F2_LocoScore_Sex_Intergenic_M6: Mean-Variance Trend")
dev.off()
#Heteroskedasticity now gone! How pretty!

write.fit(efit, adjust="BH", file="F2_Limma_Results__LocoScore_Sex_Intergenic_M6.txt")

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.05)
#How many results meet an FDR<0.05?

summary(dt)
#         (Intercept) Total_LocoScore Sex_AsFactorFemale Percent.Intergenic
# Down          1005               0               2279               6498
# NotSig         581           14056               9015                758
# Up           12470               0               2762               6800


dt<-decideTests(efit, adjust.method = "BH", p.value = 0.10)
summary(dt)
# 
#           (Intercept) Total_LocoScore Sex_AsFactorFemale Percent.Intergenic
# Down          1048               0               2898               6549
# NotSig         492           14056               7872                642
# Up           12516               0               3286               6865
# 
# 

#Conclusions:
#Sex seems to matter more if %rRNA isn't included in the model (??)
#Locoscore still doesn't matter

###########                                                                           
# 

# Can our sample size handle these types of models?  Probably.
# Sample size = 245
# 
#E.g. a DE Model with DF=18 
# 245/18 = 13.6
# So basically 13 df to estimate each parameter - would mostly be able to reliably detect effect sizes *larger* than d=1 (or group differences bigger than 1 SD)
# 
# 
# E.g. a DE Model with DF=13 (once we include a variable of interest)
# 245/13 = 18.8
# So basically 19 df to estimate each parameter - would be able to reliably detect effect sizes a little smaller than d=1 (or group differences bigger than 1 SD)

##########################