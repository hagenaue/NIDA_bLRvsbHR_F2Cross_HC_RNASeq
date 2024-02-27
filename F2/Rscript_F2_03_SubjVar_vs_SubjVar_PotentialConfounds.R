#F2 HC RNA-Seq Dataset
#03_Comparing Subject variables to Subject variables and screening for confounds
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.


colnames(F2_PCA_wMetaData)


pdf("F2_Total_LocoScore_vs_HPC_Dissection_Date_AsFactor.pdf", height=6, width=8)
boxplot(Total_LocoScore~HPC_Dissection_Date_AsFactor, data=F2_PCA_wMetaData, main="F2", ylab="Total_LocoScore")
stripchart(Total_LocoScore~HPC_Dissection_Date_AsFactor, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(Total_LocoScore~HPC_Dissection_Date_AsFactor, data=F2_PCA_wMetaData))

# Call:
#   lm(formula = Total_LocoScore ~ HPC_Dissection_Date_AsFactor, 
#      data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -628.91 -224.63  -53.63  193.17  852.59 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                              904.91      87.12  10.387  < 2e-16 ***
#   HPC_Dissection_Date_AsFactor12/10/2020  -444.50     100.23  -4.435 1.42e-05 ***
#   HPC_Dissection_Date_AsFactor12/11/2020  -347.84     101.85  -3.415 0.000751 ***
#   HPC_Dissection_Date_AsFactor12/14/2020  -456.28     109.48  -4.168 4.32e-05 ***
#   HPC_Dissection_Date_AsFactor12/2/2020   -428.03     113.18  -3.782 0.000197 ***
#   HPC_Dissection_Date_AsFactor12/3/2020   -396.77     116.42  -3.408 0.000770 ***
#   HPC_Dissection_Date_AsFactor12/4/2020   -291.52     102.82  -2.835 0.004979 ** 
#   HPC_Dissection_Date_AsFactor12/7/2020   -299.08     102.32  -2.923 0.003805 ** 
#   HPC_Dissection_Date_AsFactor12/8/2020   -427.76     100.23  -4.268 2.86e-05 ***
#   HPC_Dissection_Date_AsFactor12/9/2020   -432.44     101.85  -4.246 3.14e-05 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 289 on 235 degrees of freedom
# Multiple R-squared:  0.1115,	Adjusted R-squared:  0.0775 
# F-statistic: 3.278 on 9 and 235 DF,  p-value: 0.0008841


pdf("F2_Total_LocoScore_vs_HPC_Dissection_Date_Collapsed.pdf", height=6, width=8)
boxplot(Total_LocoScore~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData, main="F2", ylab="Total_LocoScore")
stripchart(Total_LocoScore~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(Total_LocoScore~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData))
# 
# Call:
#   lm(formula = Total_LocoScore ~ HPC_Dissection_Date_Collapsed, 
#      data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -628.91 -234.91  -47.91  195.09  843.09 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           469.91      29.19  16.100  < 2e-16 ***
#   HPC_Dissection_Date_CollapsedBatch1   435.00      91.88   4.735 3.76e-06 ***
#   HPC_Dissection_Date_CollapsedBatch2    93.84      52.43   1.790   0.0748 .  
# HPC_Dissection_Date_CollapsedBatch3    38.23      82.55   0.463   0.6437    
# HPC_Dissection_Date_CollapsedBatch4   135.92      61.08   2.225   0.0270 *  
#   HPC_Dissection_Date_CollapsedBatch6    45.11      50.55   0.892   0.3731    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 288.9 on 239 degrees of freedom
# Multiple R-squared:  0.09651,	Adjusted R-squared:  0.07761 
# F-statistic: 5.106 on 5 and 239 DF,  p-value: 0.0001813


pdf("F2_Total_LocoScore_vs_HPC_Dissector_DissectionDate.pdf", height=6, width=8)
boxplot(Total_LocoScore~HPC_Dissector_DissectionDate, data=F2_PCA_wMetaData, main="F2", ylab="Total_LocoScore")
stripchart(Total_LocoScore~HPC_Dissector_DissectionDate, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(Total_LocoScore~HPC_Dissector_DissectionDate, data=F2_PCA_wMetaData))
# 
# Call:
#   lm(formula = Total_LocoScore ~ HPC_Dissection_Date_Collapsed, 
#      data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -628.91 -234.91  -47.91  195.09  843.09 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           469.91      29.19  16.100  < 2e-16 ***
#   HPC_Dissection_Date_CollapsedBatch1   435.00      91.88   4.735 3.76e-06 ***
#   HPC_Dissection_Date_CollapsedBatch2    93.84      52.43   1.790   0.0748 .  
# HPC_Dissection_Date_CollapsedBatch3    38.23      82.55   0.463   0.6437    
# HPC_Dissection_Date_CollapsedBatch4   135.92      61.08   2.225   0.0270 *  
#   HPC_Dissection_Date_CollapsedBatch6    45.11      50.55   0.892   0.3731    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 288.9 on 239 degrees of freedom
# Multiple R-squared:  0.09651,	Adjusted R-squared:  0.07761 
# F-statistic: 5.106 on 5 and 239 DF,  p-value: 0.0001813


pdf("F2_Total_LocoScore_vs_HPC_Dissector_DissectionDate_Collapsed.pdf", height=6, width=8)
boxplot(Total_LocoScore~HPC_Dissector_DissectionDate_Collapsed, data=F2_PCA_wMetaData, main="F2", ylab="Total_LocoScore")
stripchart(Total_LocoScore~HPC_Dissector_DissectionDate_Collapsed, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(Total_LocoScore~HPC_Dissector_DissectionDate_Collapsed, data=F2_PCA_wMetaData))
#
# Call:
#   lm(formula = Total_LocoScore ~ HPC_Dissector_DissectionDate_Collapsed, 
#      data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -628.91 -217.67  -50.06  190.86  872.33 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                     440.67      46.40   9.498  < 2e-16 ***
#   HPC_Dissector_DissectionDate_CollapsedBatch1     66.26      92.79   0.714   0.4759    
# HPC_Dissector_DissectionDate_CollapsedBatch10    33.18      66.04   0.502   0.6159    
# HPC_Dissector_DissectionDate_CollapsedBatch11    90.40      69.11   1.308   0.1921    
# HPC_Dissector_DissectionDate_CollapsedBatch2    214.70      98.92   2.171   0.0310 *  
#   HPC_Dissector_DissectionDate_CollapsedBatch4     72.33     112.46   0.643   0.5207    
# HPC_Dissector_DissectionDate_CollapsedBatch5     19.11     107.15   0.178   0.8586    
# HPC_Dissector_DissectionDate_CollapsedBatch6    464.24      98.92   4.693 4.58e-06 ***
#   HPC_Dissector_DissectionDate_CollapsedBatch7     36.21      86.02   0.421   0.6742    
# HPC_Dissector_DissectionDate_CollapsedBatch8     67.48      90.27   0.747   0.4555    
# HPC_Dissector_DissectionDate_CollapsedBatch9    148.30      60.89   2.436   0.0156 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 289.7 on 234 degrees of freedom
# Multiple R-squared:  0.1105,	Adjusted R-squared:  0.07248 
# F-statistic: 2.907 on 10 and 234 DF,  p-value: 0.0019
# 

pdf("F2_Total_LocoScore_vs_SequencingBatch.pdf", height=6, width=8)
boxplot(Total_LocoScore~SequencingBatch, data=F2_PCA_wMetaData, main="F2", ylab="Total_LocoScore")
stripchart(Total_LocoScore~SequencingBatch, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(Total_LocoScore~SequencingBatch, data=F2_PCA_wMetaData))
# 
# 
# Call:
#   lm(formula = Total_LocoScore ~ SequencingBatch, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -480.73 -243.28  -50.73  208.72  864.12 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        596.88      30.51  19.565  < 2e-16 ***
#   SequencingBatch2  -119.61      43.26  -2.765  0.00613 ** 
#   SequencingBatch3   -76.15      50.09  -1.520  0.12978    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 297.3 on 242 degrees of freedom
# Multiple R-squared:  0.03115,	Adjusted R-squared:  0.02314 
# F-statistic:  3.89 on 2 and 242 DF,  p-value: 0.02173


pdf("F2_Percent.Intronic_vs_Total_LocoScore.pdf", height=6, width=6)
plot(Percent.Intronic~Total_LocoScore, data=F2_PCA_wMetaData, pch=18, main="F2", ylab="Percent.Intronic", xlab="Total_LocoScore")
BestFitLine<-lm(Percent.Intronic~Total_LocoScore, data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(Percent.Intronic~Total_LocoScore, data=F2_PCA_wMetaData))
# 
# Call:
#   lm(formula = Percent.Intronic ~ Total_LocoScore, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0058472 -0.0023458 -0.0003845  0.0015158  0.0159112 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      3.743e-02  4.375e-04  85.549   <2e-16 ***
#   Total_LocoScore -1.690e-06  7.146e-07  -2.365   0.0188 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.003358 on 243 degrees of freedom
# Multiple R-squared:  0.0225,	Adjusted R-squared:  0.01847 
# F-statistic: 5.593 on 1 and 243 DF,  p-value: 0.01882


pdf("F2_DV200_vs_EPM_Percent_Time_Open_Arm.pdf", height=6, width=6)
plot(DV200~EPM_Percent_Time_Open_Arm, data=F2_PCA_wMetaData, pch=18, main="F2", ylab="DV200", xlab="EPM_Percent_Time_Open_Arm")
BestFitLine<-lm(DV200~EPM_Percent_Time_Open_Arm, data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(DV200~EPM_Percent_Time_Open_Arm, data=F2_PCA_wMetaData))

# 
# Call:
#   lm(formula = DV200 ~ EPM_Percent_Time_Open_Arm, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.33525 -0.14114  0.03373  0.17134  0.76766 
# 
# Coefficients:
#   Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)               96.843668   0.027654 3502.024   <2e-16 ***
#   EPM_Percent_Time_Open_Arm -0.003200   0.001608   -1.989   0.0478 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.2633 on 243 degrees of freedom
# Multiple R-squared:  0.01603,	Adjusted R-squared:  0.01198 
# F-statistic: 3.958 on 1 and 243 DF,  p-value: 0.04777
# 

pdf("F2_EPM_Percent_Time_Open_Arm_vs_Age.days.pdf", height=6, width=6)
plot(EPM_Percent_Time_Open_Arm~Age.days., data=F2_PCA_wMetaData, pch=18, main="F2", ylab="EPM_Percent_Time_Open_Arm", xlab="Age (days)")
BestFitLine<-lm(EPM_Percent_Time_Open_Arm~Age.days., data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(EPM_Percent_Time_Open_Arm~Age.days., data=F2_PCA_wMetaData))
# 
# 
# Call:
#   lm(formula = EPM_Percent_Time_Open_Arm ~ Age.days., data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -16.945  -8.075  -2.037   6.681  29.443 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -71.8495    17.6582  -4.069 6.39e-05 ***
#   Age.days.     0.6872     0.1418   4.845 2.26e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 10.03 on 243 degrees of freedom
# Multiple R-squared:  0.08809,	Adjusted R-squared:  0.08434 
# F-statistic: 23.47 on 1 and 243 DF,  p-value: 2.257e-06


pdf("F2_EPM_DistanceTraveled_vs_Age.days.pdf", height=6, width=6)
plot(EPM_DistanceTraveled~Age.days., data=F2_PCA_wMetaData, pch=18, main="F2", ylab="EPM_DistanceTraveled", xlab="Age (days)")
BestFitLine<-lm(EPM_DistanceTraveled~Age.days., data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(EPM_DistanceTraveled~Age.days., data=F2_PCA_wMetaData))

# 
# Call:
#   lm(formula = EPM_DistanceTraveled ~ Age.days., data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -872.5 -321.3  -46.8  237.4 3211.4 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1936.239    813.843  -2.379   0.0181 *  
#   Age.days.      32.243      6.537   4.932 1.51e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 462.2 on 243 degrees of freedom
# Multiple R-squared:  0.09101,	Adjusted R-squared:  0.08727 
# F-statistic: 24.33 on 1 and 243 DF,  p-value: 1.507e-06


pdf("F2_EPM_Time_Immobile_vs_Age.days.pdf", height=6, width=6)
plot(EPM_Time_Immobile~Age.days., data=F2_PCA_wMetaData, pch=18, main="F2", ylab="EPM_Time_Immobile", xlab="Age (days)")
BestFitLine<-lm(EPM_Time_Immobile~Age.days., data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(EPM_Time_Immobile~Age.days., data=F2_PCA_wMetaData))

# 
# Call:
#   lm(formula = EPM_Time_Immobile ~ Age.days., data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -55.158 -15.150  -1.839  13.880  65.801 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 462.5872    39.3168  11.766   <2e-16 ***
#   Age.days.    -2.9199     0.3158  -9.246   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 22.33 on 243 degrees of freedom
# Multiple R-squared:  0.2602,	Adjusted R-squared:  0.2572 
# F-statistic: 85.49 on 1 and 243 DF,  p-value: < 2.2e-16


pdf("F2_EPM_Percent_Time_Open_Arm_vs_DOD.pdf", height=6, width=6)
boxplot(EPM_Percent_Time_Open_Arm~DOD, data=F2_PCA_wMetaData, pch=18, main="F2", ylab="EPM_Percent_Time_Open_Arm", xlab="DOD")
stripchart(EPM_Percent_Time_Open_Arm~DOD, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()


summary.lm(lm(EPM_Percent_Time_Open_Arm~DOD, data=F2_PCA_wMetaData))
# 
# Call:
#   lm(formula = EPM_Percent_Time_Open_Arm ~ DOD, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -18.586  -8.020  -1.808   6.882  29.442 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    11.3018     1.9007   5.946 9.82e-09 ***
#   DOD12/11/2013   0.6862     2.3349   0.294  0.76910    
# DOD12/13/2013   5.2560     3.8538   1.364  0.17392    
# DOD12/17/2013   3.1998     2.3144   1.383  0.16812    
# DOD12/18/2013   7.3146     2.3279   3.142  0.00189 ** 
#   DOD12/4/2013   -8.1785     6.1098  -1.339  0.18200    
# DOD12/5/2013   -0.9618     3.7051  -0.260  0.79541    
# DOD12/6/2013   -2.0707     3.8538  -0.537  0.59156    
# DOD12/9/2013   -3.1036     3.0924  -1.004  0.31659    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 10.06 on 236 degrees of freedom
# Multiple R-squared:  0.1093,	Adjusted R-squared:  0.07908 
# F-statistic: 3.619 on 8 and 236 DF,  p-value: 0.0005449


pdf("F2_EPM_DistanceTraveled_vs_DOD.pdf", height=6, width=6)
boxplot(EPM_DistanceTraveled~DOD, data=F2_PCA_wMetaData, pch=18, main="F2", ylab="EPM_DistanceTraveled", xlab="DOD")
stripchart(EPM_DistanceTraveled~DOD, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(EPM_DistanceTraveled~DOD, data=F2_PCA_wMetaData))
# 
# Call:
#   lm(formula = EPM_DistanceTraveled ~ DOD, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -855.47 -312.81  -46.83  266.99 3071.51 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    2137.31      86.86  24.607   <2e-16 ***
#   DOD12/11/2013  -178.09     106.70  -1.669   0.0964 .  
# DOD12/13/2013   -89.84     176.11  -0.510   0.6104    
# DOD12/17/2013   -29.41     105.77  -0.278   0.7812    
# DOD12/18/2013   179.89     106.38   1.691   0.0921 .  
# DOD12/4/2013   -584.68     279.21  -2.094   0.0373 *  
#   DOD12/5/2013   -299.22     169.32  -1.767   0.0785 .  
# DOD12/6/2013   -343.56     176.11  -1.951   0.0523 .  
# DOD12/9/2013   -300.87     141.32  -2.129   0.0343 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 459.6 on 236 degrees of freedom
# Multiple R-squared:  0.1271,	Adjusted R-squared:  0.09753 
# F-statistic: 4.296 on 8 and 236 DF,  p-value: 7.62e-05


pdf("F2_EPM_Time_Immobile_vs_DOD.pdf", height=6, width=6)
boxplot(EPM_Time_Immobile~DOD, data=F2_PCA_wMetaData, pch=18, main="F2", ylab="EPM_Time_Immobile", xlab="DOD")
stripchart(EPM_Time_Immobile~DOD, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(EPM_Time_Immobile~DOD, data=F2_PCA_wMetaData))

# 
# Call:
#   lm(formula = EPM_Time_Immobile ~ DOD, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -46.314 -15.177  -1.057  13.572  72.013 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   107.6321     4.1242  26.097  < 2e-16 ***
#   DOD12/11/2013  -0.5651     5.0664  -0.112 0.911292    
# DOD12/13/2013 -12.4821     8.3623  -1.493 0.136858    
# DOD12/17/2013 -21.9063     5.0220  -4.362 1.93e-05 ***
#   DOD12/18/2013 -21.7534     5.0512  -4.307 2.43e-05 ***
#   DOD12/4/2013   46.2879    13.2576   3.491 0.000573 ***
#   DOD12/5/2013   16.1619     8.0396   2.010 0.045540 *  
#   DOD12/6/2013   16.0690     8.3623   1.922 0.055860 .  
# DOD12/9/2013    8.6261     6.7101   1.286 0.199862    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 21.82 on 236 degrees of freedom
# Multiple R-squared:  0.3138,	Adjusted R-squared:  0.2905 
# F-statistic: 13.49 on 8 and 236 DF,  p-value: 4.765e-16



str(F2_PCA_wMetaData)

#######Learning Classification here with HPC_Dissection_Date variables--this did not work yet

pdf("F2_LearningClassification_AsFactor_vs_HPC_Dissection_Date_AsFactor.pdf", height=6, width=8)
boxplot(LearningClassification_AsFactor~HPC_Dissection_Date_AsFactor, data=F2_PCA_wMetaData, main="F2", ylab="LearningClassification")
stripchart(LearningClassification_AsFactor~HPC_Dissection_Date_AsFactor, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(LearningClassification~HPC_Dissection_Date_AsFactor, data=F2_PCA_wMetaData))
anova(lm(PC1~SequencingBatch, data=F2_PCA_wMetaData))


pdf("F2_Total_LocoScore_vs_HPC_Dissection_Date_Collapsed.pdf", height=6, width=8)
boxplot(Total_LocoScore~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData, main="F2", ylab="Total_LocoScore")
stripchart(Total_LocoScore~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(Total_LocoScore~HPC_Dissection_Date_Collapsed, data=F2_PCA_wMetaData))



pdf("F2_Total_LocoScore_vs_HPC_Dissector_DissectionDate.pdf", height=6, width=8)
boxplot(Total_LocoScore~HPC_Dissector_DissectionDate, data=F2_PCA_wMetaData, main="F2", ylab="Total_LocoScore")
stripchart(Total_LocoScore~HPC_Dissector_DissectionDate, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(Total_LocoScore~HPC_Dissector_DissectionDate, data=F2_PCA_wMetaData))


pdf("F2_Total_LocoScore_vs_HPC_Dissector_DissectionDate_Collapsed.pdf", height=6, width=8)
boxplot(Total_LocoScore~HPC_Dissector_DissectionDate_Collapsed, data=F2_PCA_wMetaData, main="F2", ylab="Total_LocoScore")
stripchart(Total_LocoScore~HPC_Dissector_DissectionDate_Collapsed, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(Total_LocoScore~HPC_Dissector_DissectionDate_Collapsed, data=F2_PCA_wMetaData))
#



