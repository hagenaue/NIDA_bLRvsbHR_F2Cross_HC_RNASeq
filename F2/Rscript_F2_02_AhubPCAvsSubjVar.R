#F2 HC RNA-Seq Dataset
#02_Comparing the PCA output from Ahub to our MetaData to identify potential sources of large scale noise
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.



#Taking a quick peak at the top principal components of variation (from Ahub) and other technical and categorical variables in the RNA-Seq data
#Which co-variates seem to be driving them? (or at least strongly related)
#These co-variates are things we should consider including in our differential expression model


colnames(F2_PCA_wMetaData)
# [1] "row_name"                        "PC1"                             "PC2"                             "PC3"                            
# [5] "PC4"                             "PC5"                             "PC6"                             "RNAconc"                        
# [9] "PF.Reads"                        "RibosomePerc"                    "mRNAperc"                        "RIN"                            
# [13] "DV200"                           "HPC_Dissection_Date"             "HPC_Dissection_Day"              "HPC_Side_Used"                  
# [17] "Dissector"                       "Age.days."                       "Family"                          "Sex"                            
# [21] "Rat"                             "Lateral_LocoScore"               "Rearing_LocoScore"               "Total_LocoScore"                
# [25] "EPM_Percent_Time_Open_Arm"       "EPM_Boli"                        "EPM_DistanceTraveled"            "EPM_Time_Immobile"              
# [29] "LearningClassification"          "PCA_Index_Days6and.7"            "Response_Bias_Day7"              "Prob_Diff_Day7"                 
# [33] "Latency_Score_Day7"              "PCA_Index_Day7only"              "Response_Bias_Day6"              "Prob_Diff_Day6"                 
# [37] "Latency_Score_Day6"              "PCA_Index_Day6only"              "EPM_RawDataFile"                 "EPM_Trial"                      
# [41] "EPM_TestDate"                    "EPM_TestDay"                     "EPM_TestOrderByDay"              "Loco_TestDate"                  
# [45] "Loco_TestDay"                    "Loco_Batch"                      "Loco_Rack"                       "Loco_Box"                       
# [49] "DOB"                             "DOD"                             "Sex_AsFactor"                    "Sex_AsNumeric"                  
# [53] "Dissector_AsFactor"              "Dissector_AsNumeric"             "HPC_Side_AsFactor"               "HPC_Side_AsNumeric"             
# [57] "HPC_Dissection_Date_AsFactor"    "LearningClassification_AsFactor" "STGT_experience"                 "STGT_experience_AsNumeric"  

F2_PCA_MetaData_CorMatrix<-cor(as.matrix(F2_PCA_wMetaData[,c(2:13,15,18, 22:28, 30:38, 42:43, 45, 52, 54, 56, 61)]), use="pairwise.complete.obs")

write.csv(F2_PCA_MetaData_CorMatrix, "F2_PCA_MetaData_CorMatrix.csv")

#cannot have ST/GT experience and PCA variables in correlation matrix at same time

F2_PCA_MetaData_CorMatrix<-cor(as.matrix(F2_PCA_wMetaData[,c(2:13,15,18, 22:28, 42:43, 45, 52, 54, 56, 60, 61)]), use="pairwise.complete.obs")

write.csv(F2_PCA_MetaData_CorMatrix, "F2_PCA_MetaData_CorMatrix_withoutPCAvar.csv")

F2_PCA_MetaData_CorMatrix_noDistanceOutlier<-cor(as.matrix(F2_PCA_wMetaData[F2_PCA_wMetaData$EPM_DistanceTraveled<5000,c(2:13,15,18, 22:28, 30:38, 42:43, 45, 52, 54, 56, 61)]), use="pairwise.complete.obs")

write.csv(F2_PCA_MetaData_CorMatrix_noDistanceOutlier, "F2_PCA_MetaData_CorMatrix_noDistanceOutlier.csv")



#Checking whether there are hidden batch variables and/or whether any of our existing batch variables can be collapsed into more concise versions:


pdf("F2_PC1vsSeqID_Numeric.pdf", width=6, height=6)
plot(PC1~SeqID_Numeric, data=F2_PCA_wMetaData, pch=18, main="F2", xlab="SeqID_Numeric")
BestFitLine<-lm(PC1~SeqID_Numeric, data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()

#Looks like there are three clear sequencing batches, and the third one has different PC1 than the others.

pdf("F2_PC2vsSeqID_Numeric.pdf", width=6, height=6)
plot(PC2~SeqID_Numeric, data=F2_PCA_wMetaData, pch=18, main="F2", xlab="SeqID_Numeric")
BestFitLine<-lm(PC2~SeqID_Numeric, data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()
#Same sequencing batch looks weird here.

#What is the relationship with dissection day?
pdf("F2_DissectionDayvsSeqID_Numeric.pdf", width=6, height=6)
plot(HPC_Dissection_Day~SeqID_Numeric, data=F2_PCA_wMetaData, pch=18, main="F2", xlab="SeqID_Numeric")
BestFitLine<-lm(HPC_Dissection_Day~SeqID_Numeric, data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()
#Huh. That's fascinating. The dissection batches tightly overlap with the sequencing batches.

#Looks like we need a new variable: Sequencing Batch

F2_PCA_wMetaData$SequencingBatch<-F2_PCA_wMetaData$SeqID_Numeric

F2_PCA_wMetaData$SequencingBatch[F2_PCA_wMetaData$SeqID_Numeric<470200]<-1
F2_PCA_wMetaData$SequencingBatch[F2_PCA_wMetaData$SeqID_Numeric>470200 & F2_PCA_wMetaData$SeqID_Numeric<470700]<-2
F2_PCA_wMetaData$SequencingBatch[F2_PCA_wMetaData$SeqID_Numeric>470700]<-3

table(F2_PCA_wMetaData$SequencingBatch)
# 1  2  3 
# 95 95 60

F2_PCA_wMetaData$SequencingBatch<-as.factor(as.character(F2_PCA_wMetaData$SequencingBatch))

pdf("F2_PC1vsSequencingBatch.pdf", height=6, width=4)
boxplot(PC1~SequencingBatch, data=F2_PCA_wMetaData, main="F2", ylab="PC1")
stripchart(PC1~SequencingBatch, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

anova(lm(PC1~SequencingBatch, data=F2_PCA_wMetaData))
# Analysis of Variance Table
# 
# Response: PC1
#                    Df Sum Sq  Mean Sq F value    Pr(>F)    
# SequencingBatch   2 0.1546 0.077298  22.584 9.826e-10 ***
#   Residuals       247 0.8454 0.003423                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


pdf("F2_PC2vsSequencingBatch.pdf", height=6, width=4)
boxplot(PC2~SequencingBatch, data=F2_PCA_wMetaData, main="F2", ylab="PC2")
stripchart(PC2~SequencingBatch, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

anova(lm(PC2~SequencingBatch, data=F2_PCA_wMetaData))
# Analysis of Variance Table
# 
# Response: PC2
# Df Sum Sq  Mean Sq F value    Pr(>F)    
# SequencingBatch   2 0.1178 0.058898   16.49 1.896e-07 ***
#   Residuals       247 0.8822 0.003572                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

pdf("F2_PC3vsSequencingBatch.pdf", height=6, width=4)
boxplot(PC3~SequencingBatch, data=F2_PCA_wMetaData, main="F2", ylab="PC3")
stripchart(PC3~SequencingBatch, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

anova(lm(PC3~SequencingBatch, data=F2_PCA_wMetaData))

# Analysis of Variance Table
# 
# Response: PC3
# Df  Sum Sq  Mean Sq F value    Pr(>F)    
# SequencingBatch   2 0.31924 0.159621  57.915 < 2.2e-16 ***
#   Residuals       247 0.68076 0.002756                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Absolutely *stunning* artifact. - Megan, come back to this and use it as an example 

pdf("F2_PC2vsRIN.pdf", width=6, height=6)
plot(PC2~RIN, data=F2_PCA_wMetaData, pch=18, main="F2", xlab="RIN")
BestFitLine<-lm(PC2~RIN, data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()

#Looks like we may have two outliers in RIN - 6.5 and 7.5, otherwise there is a very strong relationship with PC2, which is surprising because all of the other RIN values are >8.3

pdf("F2_PC2vsRIN_noOutliers.pdf", width=6, height=6)
plot(PC2~RIN, data=F2_PCA_wMetaData[F2_PCA_wMetaData$RIN>8,], pch=18, main="F2", xlab="RIN")
BestFitLine<-lm(PC2~RIN, data=F2_PCA_wMetaData[F2_PCA_wMetaData$RIN>8,])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(BestFitLine)

# Call:
#   lm(formula = PC2 ~ RIN, data = F2_PCA_wMetaData[F2_PCA_wMetaData$RIN > 
#                                                     8, ])
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.143162 -0.040007 -0.005504  0.034454  0.172943 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.36499    0.14992  -9.105   <2e-16 ***
#   RIN          0.15123    0.01661   9.107   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.05511 on 246 degrees of freedom
# Multiple R-squared:  0.2521,	Adjusted R-squared:  0.2491 
# F-statistic: 82.94 on 1 and 246 DF,  p-value: < 2.2e-16


#Absolutely gorgeous artifact.

#I'm going to guess that sequencing ID and RIN go partially together.

pdf("F2_RINvsSequencingBatch_noOutliers.pdf", height=6, width=4)
boxplot(RIN~SequencingBatch, data=F2_PCA_wMetaData[F2_PCA_wMetaData$RIN>8,], main="F2", ylab="RIN")
stripchart(RIN~SequencingBatch, data=F2_PCA_wMetaData[F2_PCA_wMetaData$RIN>8,], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

#Nope - Pretty unconvincing.

levels(F2_PCA_wMetaData$Dissector_AsFactor)<-list(One="EHB",Two="KA")
levels(F2_PCA_wMetaData$Dissector_AsFactor)

levels(F2_PCA_wMetaData$Dissector)<-list(One="EHB",Two="KA")
levels(F2_PCA_wMetaData$Dissector)
# $One
# [1] "EHB"
# 
# $Two
# [1] "KA"

str(F2_PCA_wMetaData)


pdf("F2_PC3vsFamily_Boxplot.pdf", height=6, width=6)
boxplot(PC3~Family, data=F2_PCA_wMetaData, main="F2", ylab="PC3")
stripchart(PC3~Family, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

anova(lm(PC3~Family, data=F2_PCA_wMetaData))
# Analysis of Variance Table
# 
# Response: PC3
# Df  Sum Sq   Mean Sq F value    Pr(>F)    
# Family     11 0.12889 0.0117170  3.2012 0.0004445 ***
#   Residuals 238 0.87111 0.0036601                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


pdf("F2_PC3vsHPC_Dissection_Date_Boxplot.pdf", height=6, width=6)
boxplot(PC3~HPC_Dissection_Date, data=F2_PCA_wMetaData, main="F2", ylab="PC3")
stripchart(PC3~HPC_Dissection_Date, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()


pdf("F2_PC3vsHPC_Dissection_Date_Boxplot_Width12.pdf", height=6, width=12)
boxplot(PC3~HPC_Dissection_Date, data=F2_PCA_wMetaData, main="F2", ylab="PC3")
stripchart(PC3~HPC_Dissection_Date, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()


anova(lm(PC3~HPC_Dissection_Date, data=F2_PCA_wMetaData))
# Analysis of Variance Table
# 
#Response: PC3
#                        Df  Sum Sq  Mean Sq F value    Pr(>F)    
# HPC_Dissection_Date   9 0.42657 0.047397  19.837 < 2.2e-16 ***
#   Residuals           240 0.57343 0.002389                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


pdf("F2_PC1vsHPC_Dissection_Date_Boxplot_Width12.pdf", height=6, width=12)
boxplot(PC1~HPC_Dissection_Date, data=F2_PCA_wMetaData, main="F2", ylab="PC1")
stripchart(PC1~HPC_Dissection_Date, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

anova(lm(PC1~HPC_Dissection_Date, data=F2_PCA_wMetaData))

# Analysis of Variance Table

# Response: PC1
#                      Df  Sum Sq   Mean Sq F value    Pr(>F)    
# HPC_Dissection_Date   9 0.17506 0.0194512  5.6589 4.191e-07 ***
#   Residuals           240 0.82494 0.0034372                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


pdf("F2_PC2vsHPC_Dissection_Date_Boxplot_Width12.pdf", height=6, width=12)
boxplot(PC2~HPC_Dissection_Date, data=F2_PCA_wMetaData, main="F2", ylab="PC2")
stripchart(PC2~HPC_Dissection_Date, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

anova(lm(PC2~HPC_Dissection_Date, data=F2_PCA_wMetaData))
# Analysis of Variance Table
# 
# Response: PC2
# Df  Sum Sq   Mean Sq F value    Pr(>F)    
# HPC_Dissection_Date   9 0.12234 0.0135933  3.7172 0.0002183 ***
#   Residuals           240 0.87766 0.0036569                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


pdf("F2_PC3vsDissectionDay.pdf", width=6, height=6)
plot(PC3~HPC_Dissection_Day, data=F2_PCA_wMetaData, pch=18, main="F2", xlab="Dissection Day")
BestFitLine<-lm(PC3~HPC_Dissection_Day, data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()


summary.lm(lm(PC3~HPC_Dissection_Day, data=F2_PCA_wMetaData))
# Call:
#   lm(formula = PC3 ~ HPC_Dissection_Day, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.144008 -0.041771 -0.003346  0.045104  0.245925 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)         0.018000   0.009604   1.874   0.0621 .
# HPC_Dissection_Day -0.002327   0.001130  -2.060   0.0405 *
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.06296 on 248 degrees of freedom
# Multiple R-squared:  0.01682,	Adjusted R-squared:  0.01285 
# F-statistic: 4.242 on 1 and 248 DF,  p-value: 0.04047


#Can we come up with a more efficient variable for dissection batch?

table(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor)

F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Releveled<-relevel(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor, ref="12/3/2020")


PCA_ByDissectionBatch<-cbind(tapply(F2_PCA_wMetaData$PC1, F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Releveled, mean),
                             tapply(F2_PCA_wMetaData$PC2, F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Releveled, mean),
                             tapply(F2_PCA_wMetaData$PC3, F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Releveled, mean))

write.csv(PCA_ByDissectionBatch, "PCA_ByDissectionBatch.csv")

F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed<-as.character(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Releveled)

F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed[F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed=="12/1/2020"]<-"Batch1"
F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed[F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed=="12/2/2020"]<-"Batch2"
F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed[F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed=="12/4/2020"]<-"Batch2"
F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed[F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed=="12/3/2020"]<-"Batch3"
F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed[F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed=="12/7/2020"]<-"Batch4"
F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed[F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed=="12/8/2020"]<-"Batch5"
F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed[F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed=="12/9/2020"]<-"Batch5"
F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed[F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed=="12/10/2020"]<-"Batch5"
F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed[F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed=="12/11/2020"]<-"Batch6"
F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed[F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed=="12/14/2020"]<-"Batch6"


F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed<-as.factor(F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed)
F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed<-relevel(F2_PCA_wMetaData$HPC_Dissection_Date_Collapsed, ref="Batch5")


pdf("F2_PC1vsHPC_Dissection_Batch_Collapsed_Boxplot_Width12.pdf", height=6, width=12)
boxplot(PC1~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData, main="F2", ylab="PC1")
stripchart(PC1~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(PC1~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData))

# Call:
#   lm(formula = PC1 ~ HPC_Dissection_Date_Collapsed, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.091965 -0.050759 -0.007492  0.046719  0.158487 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         -0.008807   0.005873  -1.500  0.13499    
# HPC_Dissection_Date_CollapsedBatch1 -0.055770   0.018655  -2.990  0.00308 ** 
# HPC_Dissection_Date_CollapsedBatch2 -0.005128   0.010624  -0.483  0.62974    
# HPC_Dissection_Date_CollapsedBatch3 -0.015344   0.016758  -0.916  0.36077    
# HPC_Dissection_Date_CollapsedBatch4  0.026599   0.012386   2.148  0.03274 *  
# HPC_Dissection_Date_CollapsedBatch6  0.047776   0.010041   4.758 3.35e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.05873 on 244 degrees of freedom
# Multiple R-squared:  0.1585,	Adjusted R-squared:  0.1412 
# F-statistic: 9.191 on 5 and 244 DF,  p-value: 4.987e-08

summary.lm(lm(PC1~HPC_Dissection_Date_AsFactor_Releveled, data=F2_PCA_wMetaData))
# Call:
#   lm(formula = PC1 ~ HPC_Dissection_Date_AsFactor_Releveled, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.09999 -0.05161 -0.00669  0.04675  0.14624 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                      -0.024151   0.015669  -1.541 0.124551    
# HPC_Dissection_Date_AsFactor_Releveled12/1/2020  -0.040426   0.023622  -1.711 0.088302 .  
# HPC_Dissection_Date_AsFactor_Releveled12/10/2020  0.027595   0.018466   1.494 0.136393    
# HPC_Dissection_Date_AsFactor_Releveled12/11/2020  0.071942   0.018786   3.829 0.000164 ***
# HPC_Dissection_Date_AsFactor_Releveled12/14/2020  0.049005   0.020430   2.399 0.017218 *  
# HPC_Dissection_Date_AsFactor_Releveled12/2/2020   0.009940   0.021456   0.463 0.643595    
# HPC_Dissection_Date_AsFactor_Releveled12/4/2020   0.010374   0.019191   0.541 0.589305    
# HPC_Dissection_Date_AsFactor_Releveled12/7/2020   0.041943   0.019080   2.198 0.028883 *  
# HPC_Dissection_Date_AsFactor_Releveled12/8/2020   0.013231   0.018618   0.711 0.477964    
# HPC_Dissection_Date_AsFactor_Releveled12/9/2020   0.003037   0.018976   0.160 0.872965    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.05863 on 240 degrees of freedom
# Multiple R-squared:  0.1751,	Adjusted R-squared:  0.1441 
# F-statistic: 5.659 on 9 and 240 DF,  p-value: 4.191e-07

#Based on Adj.R.squared, the collapsed variable performs basically as well as the full one.


anova(lm(PC1~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData))

# Analysis of Variance Table
# 
# Response: PC1
# Df  Sum Sq  Mean Sq F value    Pr(>F)    
# HPC_Dissection_Date_Collapsed   5 0.15848 0.031697  9.1906 4.987e-08 ***
#   Residuals                     244 0.84152 0.003449                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


anova(lm(F2_PCA_wMetaData$PC1~F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor))
# Analysis of Variance Table
# 
# Response: F2_PCA_wMetaData$PC1
# Df  Sum Sq   Mean Sq F value    Pr(>F)    
# F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor   9 0.17506 0.0194512  5.6589 4.191e-07 ***
#   Residuals                                     240 0.82494 0.0034372                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#Date collapsed performs better


anova(lm(F2_PCA_wMetaData$PC2~F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor))

# Analysis of Variance Table
# 
# Response: F2_PCA_wMetaData$PC2
# Df  Sum Sq   Mean Sq F value    Pr(>F)    
# F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor   9 0.12234 0.0135933  3.7172 0.0002183 ***
#   Residuals                                     240 0.87766 0.0036569                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



pdf("F2_PC2vsHPC_Dissection_Batch_Collapsed_Boxplot_Width12.pdf", height=6, width=12)
boxplot(PC2~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData, main="F2", ylab="PC2")
stripchart(PC2~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(PC2~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData))

# Call:
#   lm(formula = PC2 ~ HPC_Dissection_Date_Collapsed, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.145663 -0.040177 -0.007128  0.038713  0.178034 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         -0.008997   0.006056  -1.486  0.13865    
# HPC_Dissection_Date_CollapsedBatch1 -0.001615   0.019238  -0.084  0.93315    
# HPC_Dissection_Date_CollapsedBatch2 -0.005514   0.010956  -0.503  0.61519    
# HPC_Dissection_Date_CollapsedBatch3  0.045012   0.017281   2.605  0.00976 ** 
# HPC_Dissection_Date_CollapsedBatch4 -0.010732   0.012773  -0.840  0.40162    
# HPC_Dissection_Date_CollapsedBatch6  0.042130   0.010354   4.069 6.38e-05 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.06056 on 244 degrees of freedom
# Multiple R-squared:  0.1051,	Adjusted R-squared:  0.08679 
# F-statistic: 5.733 on 5 and 244 DF,  p-value: 5.059e-05


summary.lm(lm(PC2~HPC_Dissection_Date_AsFactor_Releveled, data=F2_PCA_wMetaData))

# Call:
#   lm(formula = PC2 ~ HPC_Dissection_Date_AsFactor_Releveled, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.138751 -0.043057 -0.006422  0.038218  0.174550 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                                       0.036015   0.016162   2.228  0.02678 * 
# HPC_Dissection_Date_AsFactor_Releveled12/1/2020  -0.046627   0.024365  -1.914  0.05685 . 
# HPC_Dissection_Date_AsFactor_Releveled12/10/2020 -0.037199   0.019047  -1.953  0.05198 . 
# HPC_Dissection_Date_AsFactor_Releveled12/11/2020 -0.009794   0.019377  -0.505  0.61370   
# HPC_Dissection_Date_AsFactor_Releveled12/14/2020  0.008178   0.021073   0.388  0.69830   
# HPC_Dissection_Date_AsFactor_Releveled12/2/2020  -0.069302   0.022131  -3.131  0.00196 **
# HPC_Dissection_Date_AsFactor_Releveled12/4/2020  -0.039798   0.019794  -2.011  0.04549 * 
# HPC_Dissection_Date_AsFactor_Releveled12/7/2020  -0.055744   0.019680  -2.832  0.00501 **
# HPC_Dissection_Date_AsFactor_Releveled12/8/2020  -0.045814   0.019203  -2.386  0.01782 * 
# HPC_Dissection_Date_AsFactor_Releveled12/9/2020  -0.053479   0.019573  -2.732  0.00676 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.06047 on 240 degrees of freedom
# Multiple R-squared:  0.1223,	Adjusted R-squared:  0.08943 
# F-statistic: 3.717 on 9 and 240 DF,  p-value: 0.0002183

#Based on adj.r.squared again, the collapsed variable performs just as well.


anova(lm(PC2~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData))

# Analysis of Variance Table
# 
# Response: PC2
# Df  Sum Sq   Mean Sq F value    Pr(>F)    
# HPC_Dissection_Date_Collapsed   5 0.10513 0.0210262  5.7331 5.059e-05 ***
#   Residuals                     244 0.89487 0.0036675                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


pdf("F2_PC3vsHPC_Dissection_Batch_Collapsed_Boxplot_Width12.pdf", height=6, width=12)
boxplot(PC3~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData, main="F2", ylab="PC3")
stripchart(PC3~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

anova(lm(PC3~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData))

# Analysis of Variance Table
# 
# Response: PC3
# Df  Sum Sq  Mean Sq F value    Pr(>F)    
# HPC_Dissection_Date_Collapsed   5 0.41124 0.082248  34.086 < 2.2e-16 ***
#   Residuals                     244 0.58876 0.002413                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

summary.lm(lm(PC3~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData))

# Call:
#   lm(formula = PC3 ~ HPC_Dissection_Date_Collapsed, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.127088 -0.032272  0.001482  0.033843  0.189061 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         -0.041052   0.004912  -8.357 4.91e-15 ***
# HPC_Dissection_Date_CollapsedBatch1  0.124038   0.015604   7.949 6.97e-14 ***
# HPC_Dissection_Date_CollapsedBatch2  0.078897   0.008886   8.878  < 2e-16 ***
# HPC_Dissection_Date_CollapsedBatch3  0.020761   0.014017   1.481  0.13986    
# HPC_Dissection_Date_CollapsedBatch4  0.027678   0.010360   2.672  0.00806 ** 
# HPC_Dissection_Date_CollapsedBatch6  0.083340   0.008398   9.923  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.04912 on 244 degrees of freedom
# Multiple R-squared:  0.4112,	Adjusted R-squared:  0.3992 
# F-statistic: 34.09 on 5 and 244 DF,  p-value: < 2.2e-16

summary.lm(lm(PC3~HPC_Dissection_Date_AsFactor_Releveled, data=F2_PCA_wMetaData))

# Call:
#   lm(formula = PC3 ~ HPC_Dissection_Date_AsFactor_Releveled, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.121274 -0.031428  0.000552  0.029888  0.179034 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                      -0.020290   0.013064  -1.553 0.121703    
# HPC_Dissection_Date_AsFactor_Releveled12/1/2020   0.103277   0.019694   5.244 3.45e-07 ***
# HPC_Dissection_Date_AsFactor_Releveled12/10/2020 -0.014539   0.015396  -0.944 0.345938    
# HPC_Dissection_Date_AsFactor_Releveled12/11/2020  0.056312   0.015663   3.595 0.000393 ***
# HPC_Dissection_Date_AsFactor_Releveled12/14/2020  0.072605   0.017033   4.263 2.91e-05 ***
# HPC_Dissection_Date_AsFactor_Releveled12/2/2020   0.068310   0.017888   3.819 0.000171 ***
# HPC_Dissection_Date_AsFactor_Releveled12/4/2020   0.052322   0.016000   3.270 0.001233 ** 
# HPC_Dissection_Date_AsFactor_Releveled12/7/2020   0.006917   0.015908   0.435 0.664098    
# HPC_Dissection_Date_AsFactor_Releveled12/8/2020  -0.014239   0.015522  -0.917 0.359877    
# HPC_Dissection_Date_AsFactor_Releveled12/9/2020  -0.035620   0.015821  -2.251 0.025264 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.04888 on 240 degrees of freedom
# Multiple R-squared:  0.4266,	Adjusted R-squared:  0.4051 
# F-statistic: 19.84 on 9 and 240 DF,  p-value: < 2.2e-16


anova(lm(F2_PCA_wMetaData$PC3~F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor))

# Analysis of Variance Table
# 
# Response: F2_PCA_wMetaData$PC3
# Df  Sum Sq  Mean Sq F value    Pr(>F)    
# F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor   9 0.42657 0.047397  19.837 < 2.2e-16 ***
#   Residuals                                     240 0.57343 0.002389                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#Again, based on adj.r.squared, the collapsed variable performs almost as well as the full variable.



#Looks like Age as a variable is just standing in for Sex:
pdf("F2_Age.days.vsSex.pdf", height=6, width=6)
boxplot(Age.days.~Sex, data=F2_PCA_wMetaData, main="F2", ylab="Age (days)")
stripchart(Age.days.~Sex, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(Age.days.~Sex, data=F2_PCA_wMetaData))
# Call:
#   lm(formula = Age.days. ~ Sex, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -7.528 -1.528  0.736  1.736  3.736 
# 
# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    128.2640     0.2129  602.39   <2e-16 ***
#   SexMale      -7.7360     0.3011  -25.69   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 2.381 on 248 degrees of freedom
# Multiple R-squared:  0.7269,	Adjusted R-squared:  0.7258 
# F-statistic:   660 on 1 and 248 DF,  p-value: < 2.2e-16


anova(lm(F2_PCA_wMetaData$PC1~F2_PCA_wMetaData$Family))
# Analysis of Variance Table
# 
# Response: F2_PCA_wMetaData$PC1
# Df  Sum Sq   Mean Sq F value    Pr(>F)    
# F2_PCA_wMetaData$Family  11 0.12445 0.0113134  3.0753 0.0007034 ***
#   Residuals               238 0.87555 0.0036788                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


anova(lm(F2_PCA_wMetaData$PC2~F2_PCA_wMetaData$Family))
# Analysis of Variance Table
# 
# Response: F2_PCA_wMetaData$PC2
# Df  Sum Sq   Mean Sq F value   Pr(>F)   
# F2_PCA_wMetaData$Family  11 0.11991 0.0109006  2.9478 0.001116 **
#   Residuals               238 0.88009 0.0036979                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

anova(lm(F2_PCA_wMetaData$PC3~F2_PCA_wMetaData$Family))
# Analysis of Variance Table
# 
# Response: F2_PCA_wMetaData$PC3
# Df  Sum Sq   Mean Sq F value    Pr(>F)    
# F2_PCA_wMetaData$Family  11 0.12889 0.0117170  3.2012 0.0004445 ***
#   Residuals               238 0.87111 0.0036601                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



#Checking sac date real quick:

table(F2_PCA_wMetaData$DOD, F2_PCA_wMetaData$Sex)

#             Female Male
# 12/10/2013      0   28
# 12/11/2013      0   55
# 12/13/2013      9    0
# 12/17/2013     59    0
# 12/18/2013     57    0
# 12/4/2013       0    4
# 12/5/2013       0   10
# 12/6/2013       0    9
# 12/9/2013       0   19

pdf("F2_PC1vsDODSex.pdf", height=6, width=20)
boxplot(PC1~DOD+Sex, data=F2_PCA_wMetaData, main="F2", ylab="PC1")
stripchart(PC1~DOD+Sex, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(PC1~DOD, data=F2_PCA_wMetaData))

# Call:
#   lm(formula = PC1 ~ DOD, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.121482 -0.054257 -0.001125  0.048149  0.163059 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   -0.015718   0.011673  -1.347  0.17938   
# DOD12/11/2013  0.005810   0.014339   0.405  0.68570   
# DOD12/13/2013  0.009600   0.023668   0.406  0.68540   
# DOD12/17/2013  0.036428   0.014175   2.570  0.01077 * 
# DOD12/18/2013  0.002339   0.014254   0.164  0.86981   
# DOD12/4/2013   0.029637   0.033016   0.898  0.37026   
# DOD12/5/2013   0.032746   0.022754   1.439  0.15142   
# DOD12/6/2013   0.071130   0.023668   3.005  0.00293 **
# DOD12/9/2013   0.008154   0.018359   0.444  0.65735   
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.06177 on 241 degrees of freedom
# Multiple R-squared:  0.08056,	Adjusted R-squared:  0.05004 
# F-statistic: 2.639 on 8 and 241 DF,  p-value: 0.008631

pdf("F2_PC2vsDODSex.pdf", height=6, width=20)
boxplot(PC2~DOD+Sex, data=F2_PCA_wMetaData, main="F2", ylab="PC2")
stripchart(PC2~DOD+Sex, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(PC2~DOD, data=F2_PCA_wMetaData))


# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.135919 -0.037808 -0.005634  0.031596  0.167233 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   -0.013323   0.011519  -1.157  0.24858   
# DOD12/11/2013 -0.004706   0.014150  -0.333  0.73976   
# DOD12/13/2013 -0.013688   0.023355  -0.586  0.55838   
# DOD12/17/2013  0.039677   0.013988   2.837  0.00495 **
#   DOD12/18/2013  0.005743   0.014066   0.408  0.68345   
# DOD12/4/2013   0.011069   0.032580   0.340  0.73435   
# DOD12/5/2013   0.003963   0.022454   0.176  0.86006   
# DOD12/6/2013   0.070908   0.023355   3.036  0.00266 **
#   DOD12/9/2013   0.016965   0.018117   0.936  0.34998   
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.06095 on 241 degrees of freedom
# Multiple R-squared:  0.1047,	Adjusted R-squared:  0.07494 
# F-statistic: 3.521 on 8 and 241 DF,  p-value: 0.0007155

pdf("F2_PC3vsDODSex.pdf", height=6, width=20)
boxplot(PC3~DOD+Sex, data=F2_PCA_wMetaData, main="F2", ylab="PC3")
stripchart(PC3~DOD+Sex, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(PC3~DOD, data=F2_PCA_wMetaData))


# Call:
#   lm(formula = PC3 ~ DOD, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.149976 -0.040648 -0.004066  0.038605  0.191560 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)   -0.015193   0.011681  -1.301   0.1946  
# DOD12/11/2013  0.005012   0.014350   0.349   0.7272  
# DOD12/13/2013  0.046717   0.023685   1.972   0.0497 *
#   DOD12/17/2013  0.027552   0.014185   1.942   0.0533 .
# DOD12/18/2013  0.018220   0.014265   1.277   0.2027  
# DOD12/4/2013   0.033016   0.033039   0.999   0.3186  
# DOD12/5/2013   0.029153   0.022771   1.280   0.2017  
# DOD12/6/2013   0.054983   0.023685   2.321   0.0211 *
#   DOD12/9/2013  -0.025281   0.018372  -1.376   0.1701  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.06181 on 241 degrees of freedom
# Multiple R-squared:  0.07924,	Adjusted R-squared:  0.04867 
# F-statistic: 2.592 on 8 and 241 DF,  p-value: 0.009819


#Sac date... might actually matter. But is hierarchially nested in Sex, so makes modeling much more complicated in Limma.  Going to skip for now, but have it in the back of our mind if we still need to clean up a lot of noise.
