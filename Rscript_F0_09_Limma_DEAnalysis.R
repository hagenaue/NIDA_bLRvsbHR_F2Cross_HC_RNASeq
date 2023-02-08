#F0 HC RNA-Seq Dataset
#09_Limma differential expression analysis
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-10, updated later for a few figures for the paper.

###########                                                                           

#Running differential expression analyses:

colnames(F0_PCA_wMetaData)


#Controlling for Sex and the strongest correlate with PC1 (%rRNA)

design<-model.matrix(~Lineage_AsFactor+Sex_Factor+RibosomePerc, data=F0_PCA_wMetaData)
design
# (Intercept) Lineage_AsFactorbLR Sex_FactorFemale RibosomePerc
# 1            1                   1                1     0.024283
# 2            1                   1                0     0.023922
# 3            1                   0                0     0.019931
# 4            1                   0                1     0.012275
# 5            1                   1                1     0.013630
# 6            1                   1                0     0.019710
# 7            1                   0                0     0.025424
# 8            1                   0                1     0.021875
# 9            1                   0                0     0.018035
# 11           1                   1                1     0.014526
# 12           1                   1                0     0.018377
# 13           1                   1                1     0.017040
# 14           1                   1                0     0.039761
# 15           1                   0                0     0.025320
# 16           1                   0                1     0.036275
# 17           1                   1                1     0.025660
# 18           1                   1                0     0.023043
# 19           1                   1                1     0.038016
# 20           1                   1                0     0.030922
# 21           1                   0                0     0.027106
# 22           1                   0                1     0.012322
# 23           1                   0                0     0.026196
# 24           1                   0                1     0.022625
# attr(,"assign")
# [1] 0 1 2 3
# attr(,"contrasts")
# attr(,"contrasts")$Lineage_AsFactor
# [1] "contr.treatment"
# 
# attr(,"contrasts")$Sex_Factor
# [1] "contr.treatment"
# 
# # attr(,"contrasts")$Sex_Factor
# # [1] "contr.treatment"          




v<-voom(F0_RawHitCount_dge_noLowHits_TMM, design, plot=TRUE)
v
# An object of class "EList"
# $targets
# group lib.size norm.factors
# SL483888     1 17064755     1.015457
# SL483889     1 21263036     1.120887
# SL483890     1 30995699     1.213766
# SL483891     1 29984350     1.303676
# SL483892     1 31849671     1.288386
# 18 more rows ...
# 
# $E
# SL483888 SL483889 SL483890 SL483891 SL483892 SL483893 SL483894 SL483895 SL483896 SL483898 SL483899 SL483900 SL483901 SL483902  SL483903 SL483904 SL483905 SL483906 SL483907  SL483908 SL483909 SL483910 SL483911
# ENSRNOG00000014303 6.687182 6.550088 6.718650 6.984316 7.026905 6.866643 6.679487 6.411326 6.303993 7.251160 6.472742 6.211799 6.756289 6.511979  5.722357 6.615202 6.048929 6.882992 6.705870  5.807843 7.032974 5.932293 6.578315
# ENSRNOG00000014330 9.124480 8.831335 8.534911 8.706212 8.638596 8.522793 9.355356 9.586258 9.475989 8.308396 9.550426 9.489858 9.917850 9.475019 10.483093 9.495841 9.967549 9.057023 9.666944 10.059180 8.745940 9.848881 9.625870
# ENSRNOG00000049505 3.685129 3.469307 3.228398 3.056758 3.208917 3.241441 3.643318 4.104367 3.892290 2.925101 4.045258 4.267532 4.394791 3.757179  4.721468 4.045206 4.044882 3.830099 4.083558  4.423524 3.664858 4.311358 4.285206
# ENSRNOG00000014916 4.556308 4.624523 4.586132 4.265039 4.318542 4.656890 4.369079 4.276153 3.595413 4.433725 4.227429 3.581513 5.140985 3.835807  4.350231 4.630802 3.814362 5.243522 4.707851  3.873263 4.030120 3.931813 4.816269
# ENSRNOG00000014996 4.335412 4.042995 3.933224 3.837014 3.695044 3.977092 4.299954 4.447362 4.416571 3.726791 4.326622 4.534718 4.678676 4.317370  4.879684 4.424156 4.477756 4.054333 4.557783  4.739315 3.906666 4.766084 4.595316
# 13781 more rows ...
# 
# $weights
# [,1]      [,2]     [,3]     [,4]     [,5]      [,6]      [,7]      [,8]      [,9]    [,10]     [,11]    [,12]     [,13]     [,14]     [,15]     [,16]     [,17]     [,18]     [,19]     [,20]    [,21]     [,22]     [,23]
# [1,] 10.273806 10.195769 9.979189 9.796795 9.527522 10.075205 10.332934 10.190486 10.108412 9.763635 10.088462 9.835132 10.403419 10.293990 10.410927 10.052548 10.269094 10.293165 10.354320 10.377940 9.694412 10.243956 10.254409
# [2,]  8.425706  8.153124 7.689892 8.020440 8.129542  8.205109  7.914864  7.993552  7.909902 8.297350  8.284019 8.246437  8.102689  7.836076  7.906167  8.061881  8.309885  7.812780  8.123043  7.964030 7.925853  7.716664  8.053374
# [3,]  7.996722  8.359318 9.193000 8.777389 8.675191  8.292697  8.698432  8.751476  8.770223 8.339650  8.142772 8.418212  8.325283  8.856220  8.804954  8.707830  8.049313  9.095411  8.359690  8.585835 8.960863  9.088082  8.625943
# [4,]  9.138633  9.269343 9.491928 9.398810 9.910331  9.268065  8.917118  9.230239  9.096165 9.596804  9.150884 9.633419  9.002622  9.074683  9.069806  9.775431  8.997081  9.953331  9.165549  8.782759 9.578717  9.293359  9.097367
# [5,]  8.590060  9.100734 9.918787 9.501882 9.424209  9.112253  9.354774  9.303924  9.555045 9.085047  9.000646 9.114567  8.785039  9.512491  9.099706  9.239722  8.830085  9.405095  8.976763  9.215199 9.680096  9.721790  9.169092
# 13781 more rows ...
# 
# $design
# (Intercept) Lineage_AsFactorbLR Sex_FactorFemale RibosomePerc
# 1           1                   1                1     0.024283
# 2           1                   1                0     0.023922
# 3           1                   0                0     0.019931
# 4           1                   0                1     0.012275
# 5           1                   1                1     0.013630
# 18 more rows ...

vfit<-lmFit(v,design)
#This is the linear regression output (raw)


str(vfit)
# Formal class 'MArrayLM' [package "limma"] with 1 slot
# ..@ .Data:List of 10
# .. ..$ : num [1:13786, 1:4] 6.76 8.43 2.96 3.4 3.77 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:13786] "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# .. .. .. ..$ : chr [1:4] "(Intercept)" "Lineage_AsFactorbLR" "Sex_FactorFemale" "RibosomePerc"
# .. ..$ : num [1:13786, 1:4] 0.244 0.275 0.262 0.253 0.254 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:13786] "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# .. .. .. ..$ : chr [1:4] "(Intercept)" "Lineage_AsFactorbLR" "Sex_FactorFemale" "RibosomePerc"
# .. ..$ : num [1:13786] 1.164 1.294 1.15 1.128 0.835 ...
# .. ..$ : num [1:13786] 19 19 19 19 19 19 19 19 19 19 ...
# .. ..$ : num [1:4, 1:4] 0.601 -0.053 -0.143 -19.798 -0.053 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:4] "(Intercept)" "Lineage_AsFactorbLR" "Sex_FactorFemale" "RibosomePerc"
# .. .. .. ..$ : chr [1:4] "(Intercept)" "Lineage_AsFactorbLR" "Sex_FactorFemale" "RibosomePerc"
# .. ..$ : int [1:4] 1 2 3 4
# .. ..$ : int 4
# .. ..$ : Named num [1:13786] 6.55 9.32 3.84 4.34 4.3 ...
# .. .. ..- attr(*, "names")= chr [1:13786] "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# .. ..$ : chr "ls"
# .. ..$ : num [1:23, 1:4] 1 1 1 1 1 1 1 1 1 1 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:23] "1" "2" "3" "4" ...
# .. .. .. ..$ : chr [1:4] "(Intercept)" "Lineage_AsFactorbLR" "Sex_FactorFemale" "RibosomePerc"
# .. .. ..- attr(*, "assign")= int [1:4] 0 1 2 3
# .. .. ..- attr(*, "contrasts")=List of 2
# .. .. .. ..$ Lineage_AsFactor: chr "contr.treatment"
# .. .. .. ..$ Sex_Factor      : chr "contr.treatment"
# ..$ names: chr [1:10] "coefficients" "stdev.unscaled" "sigma" "df.residual" ...

efit<-eBayes(vfit)
#This is the linear regression output with an empirical Bayes correction

png("F0_Sex_Lineage_rRNA_Model_MeanVarianceTrend.png")
plotSA(efit, main="F0_F0_Sex_Lineage_rRNA_Model: Mean-Variance Trend")
dev.off()
#Heteroskedasticity now gone! How pretty!

write.fit(efit, adjust="BH", file="F0_Limma_Results_Sex_Lineage_rRNA_.txt")

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.05)
#How many results meet an FDR<0.05?

summary(dt)
#              (Intercept) Lineage_AsFactorbLR Sex_FactorFemale RibosomePerc
# Down           436                  67               10         4756
# NotSig        1088               13666            13772         4637
# Up           12262                  53                4         4393



dt<-decideTests(efit, adjust.method = "BH", p.value = 0.1)
summary(dt)
#            (Intercept) Lineage_AsFactorbLR Sex_FactorFemale RibosomePerc
# Down           510                 113               12         5139
# NotSig         902               13568            13770         3872
# Up           12374                 105                4         4775

##################

#Model: Controlling for Sex and the strongest correlates with PC1 and PC2 (%rRNA, %UTR). Consider results carefully: %UTR also correlates with lineage, so we could be controlling for a confound but we could also be eliminating true biological signal.

design<-model.matrix(~Lineage_AsFactor+Sex_Factor+RibosomePerc+Percent_UTR, data=F0_PCA_wMetaData)
design
#   (Intercept) Lineage_AsFactorbLR Sex_FactorFemale RibosomePerc Percent_UTR
# 1            1                   1                1     0.024283    0.386023
# 2            1                   1                0     0.023922    0.375106
# 3            1                   0                0     0.019931    0.379926
# 4            1                   0                1     0.012275    0.348627
# 5            1                   1                1     0.013630    0.356408
# 6            1                   1                0     0.019710    0.387615
# 7            1                   0                0     0.025424    0.360804
# 8            1                   0                1     0.021875    0.354907
# 9            1                   0                0     0.018035    0.339328
# 11           1                   1                1     0.014526    0.369650
# 12           1                   1                0     0.018377    0.361717
# 13           1                   1                1     0.017040    0.334694
# 14           1                   1                0     0.039761    0.381602
# 15           1                   0                0     0.025320    0.340916
# 16           1                   0                1     0.036275    0.334584
# 17           1                   1                1     0.025660    0.394791
# 18           1                   1                0     0.023043    0.333451
# 19           1                   1                1     0.038016    0.410534
# 20           1                   1                0     0.030922    0.375220
# 21           1                   0                0     0.027106    0.329144
# 22           1                   0                1     0.012322    0.342515
# 23           1                   0                0     0.026196    0.336564
# 24           1                   0                1     0.022625    0.366782
# attr(,"assign")
# [1] 0 1 2 3 4
# attr(,"contrasts")
# attr(,"contrasts")$Lineage_AsFactor
# [1] "contr.treatment"
# 
# attr(,"contrasts")$Sex_Factor



v<-voom(F0_RawHitCount_dge_noLowHits_TMM, design, plot=TRUE)
v


vfit<-lmFit(v,design)

str(vfit)
#This is the linear regression output (raw)

#

efit<-eBayes(vfit)
#This is the linear regression output with an empirical Bayes correction

png("F0_Sex_Lineage_rRNA_UTR_MeanVarianceTrend.png")
plotSA(efit, main="F0_Sex_Lineage_rRNA_UTR: Mean-Variance Trend")
dev.off()
#Heteroskedasticity now gone! How pretty!

write.fit(efit, adjust="BH", file="F0_Limma_Results_Sex_Lineage_rRNA_UTR.txt")

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.05)
#How many results meet an FDR<0.05?

summary(dt)
#            (Intercept) Lineage_AsFactorbLR Sex_FactorFemale RibosomePerc Percent_UTR
# Down           739                  22               10         5639        3232
# NotSig        5339               13723            13771         3424        7189
# Up            7708                  41                5         4723        3365


dt<-decideTests(efit, adjust.method = "BH", p.value = 0.1)
#How many results meet an FDR<0.1?

summary(dt)
#              (Intercept) Lineage_AsFactorbLR Sex_FactorFemale RibosomePerc Percent_UTR
# Down           979                  30               12         5935        3955
# NotSig        4586               13703            13768         2858        5893
# Up            8221                  53                6         4993        3938


#####################

#Model: Controlling for Sex, the top correlate with PC1 (%rRNA), and a strong correlate with both PC1 & PC2 that was the top correlate with PC1 in the F2 data (% intergenic, to make the analyses as parallel as possible)

design<-model.matrix(~Lineage_AsFactor+Sex_Factor+RibosomePerc+Percent_Intergenic, data=F0_PCA_wMetaData)
design

v<-voom(F0_RawHitCount_dge_noLowHits_TMM, design, plot=TRUE)
v


vfit<-lmFit(v,design)

str(vfit)
#This is the linear regression output (raw)

#

efit<-eBayes(vfit)
#This is the linear regression output with an empirical Bayes correction

png("F0_Sex_Lineage_rRNA_Intergenic_MeanVarianceTrend.png")
plotSA(efit, main="F0_Sex_Lineage_rRNA_Intergenic: Mean-Variance Trend")
dev.off()
#Heteroskedasticity now gone! How pretty!

write.fit(efit, adjust="BH", file="F0_Limma_Results_Sex_Lineage_rRNA_Intergenic.txt")

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.05)
#How many results meet an FDR<0.05?

summary(dt)
#           (Intercept) Lineage_AsFactorbLR Sex_FactorFemale RibosomePerc Percent_Intergenic
# Down           296                  60               15         5155               4305
# NotSig        2778               13642            13763         3818               5010
# Up           10712                  84                8         4813               4471

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.1)
#How many results meet an FDR<0.1?

summary(dt)
#              (Intercept) Lineage_AsFactorbLR Sex_FactorFemale RibosomePerc Percent_Intergenic
# Down           378                  86               22         5474               4658
# NotSig        2341               13569            13743         3164               4183
# Up           11067                 131               21         5148               4945



################################################
#Using sex as an interaction term
design<-model.matrix(~Lineage_AsFactor*Sex_Factor+RibosomePerc, data=F0_PCA_wMetaData)
design


v<-voom(F0_RawHitCount_dge_noLowHits_TMM, design, plot=TRUE)
v

vfit<-lmFit(v,design)

str(vfit)
#This is the linear regression output (raw)

#

efit<-eBayes(vfit)
#This is the linear regression output with an empirical Bayes correction


png("F0_LineageSex_Interaction_RibosomePerc_MeanVarianceTrend_SexInteraction.png")
plotSA(efit, main="F0_LineageSex_Interaction_RibosomePerc: Mean-Variance Trend")
dev.off()
#Heteroskedasticity now gone! How pretty!

write.fit(efit, adjust="BH", file="F0_Limma_Results_LineageSex_Interaction_RibosomePerc.txt")

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.05)
#How many results meet an FDR<0.05?

summary(dt)
# (Intercept) Lineage_AsFactorbLR Sex_FactorFemale RibosomePerc Lineage_AsFactorbLR:Sex_FactorFemale
# Down           440                   8                9         4726                                    0
# NotSig        1116               13769            13774         4698                                13786
# Up           12230                   9                3         4362                                    0

dt<-decideTests(efit, adjust.method = "BH", p.value = 0.1)

#How many results meet an FDR<0.1?

summary(dt)

#          (Intercept) Lineage_AsFactorbLR Sex_FactorFemale RibosomePerc Lineage_AsFactorbLR:Sex_FactorFemale
# Down           503                  11               10         5124                                    0
# NotSig         944               13762            13773         3906                                13786
# Up           12339                  13                3         4756                                    0



####################################