#F2 HC RNA-Seq Dataset
#13_Screening for Noise and Potential Confounds: Comparing the PCA Output to All Variables
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.


#Examining the relationship between variables of interest, potential co-variates, and the top PCs
#Goal 1: Identify potential confounds/mediating variables (variables with relationships with our variables of interest)
#Goal 2: Identify variables that might be contributing large-scale noise (correlated with top PCs)


F2_PCA_MetaData_CorMatrix<-cor(as.matrix(cbind(F2_PCA_wMetaData[,c(2:13,15,18, 22:28, 42:43, 45, 52, 54, 56, 60, 61, 65)],RNAMetricsForFinal_F2s[,c(4:8,10)])), use="pairwise.complete.obs")

write.csv(F2_PCA_MetaData_CorMatrix, "F2_PCA_MetaData_CorMatrix_withoutPCAvar_NoOutliers_moreRNAmetrics.csv")


colnames(F2_PCA_wMetaData)

SubjectFactorVariables<-as.matrix(as.data.frame(F2_PCA_wMetaData[,c(14,19,20,29,41,44,46:48,50,59,62,64,69,70)],stringsAsFactors = TRUE))

str(as.data.frame(F2_PCA_wMetaData[,c(14,19,20,29,41,44,46:48,50,59,62,64,69,70)],stringsAsFactors = TRUE))
# 'data.frame':	245 obs. of  15 variables:
#   $ HPC_Dissection_Date                   : chr  "12/1/2020" "12/1/2020" "12/1/2020" "12/1/2020" ...
# $ Family                                : chr  "4" "9" "9" "9" ...
# $ Sex                                   : chr  "Male" "Male" "Male" "Male" ...
# $ LearningClassification                : chr  "IN" "ST" "ST" "IN" ...
# $ EPM_TestDate                          : chr  "10/21/2013" "10/25/2013" "10/22/2013" "10/25/2013" ...
# $ Loco_TestDate                         : chr  "9/24/2013" "9/25/2013" "9/25/2013" "9/25/2013" ...
# $ Loco_Batch                            : int  1 3 3 3 1 2 1 3 2 1 ...
# $ Loco_Rack                             : int  2 1 1 1 1 1 1 2 1 1 ...
# $ Loco_Box                              : int  10 13 17 15 15 4 3 15 7 18 ...
# $ DOD                                   : chr  "12/10/2013" "12/11/2013" "12/11/2013" "12/11/2013" ...
# $ STGT_experience                       : Factor w/ 2 levels "FALSE","TRUE": 2 2 2 2 2 2 2 2 2 2 ...
# $ SequencingBatch                       : Factor w/ 3 levels "1","2","3": 1 1 1 1 1 1 1 1 1 1 ...
# $ HPC_Dissection_Date_Collapsed         : Factor w/ 6 levels "Batch5","Batch1",..: 2 2 2 2 2 2 2 2 2 2 ...
# $ HPC_Dissector_DissectionDate          : chr  "KA_12/1/2020" "KA_12/1/2020" "KA_12/1/2020" "KA_12/1/2020" ...
# $ HPC_Dissector_DissectionDate_Collapsed: Factor w/ 11 levels "Batch3","Batch1",..: 8 8 8 8 8 8 8 8 8 8 ...

colnames(SubjectContinuousVariables)

SubjectContinuousVariables<-as.matrix(F2_PCA_wMetaData[,c(2:13,15,18, 22:28, 42:43, 45, 46:48, 52, 54, 56, 61, 65:68,74:77,80:84)])

str(SubjectFactorVariables)
# chr [1:245, 1:12] "12/1/2020" "12/1/2020" "12/1/2020" "12/1/2020" "12/1/2020" "12/1/2020" "12/1/2020" "12/1/2020" "12/1/2020" "12/1/2020" "12/2/2020" "12/2/2020" "12/2/2020" "12/2/2020" "12/2/2020" "12/2/2020" "12/2/2020" "12/2/2020" "12/2/2020" ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:245] "63" "181" "185" "183" ...
# ..$ : chr [1:12] "HPC_Dissection_Date" "Family" "Sex" "LearningClassification" ...

head(SubjectFactorVariables)
HPC_Dissection_Date Family Sex      LearningClassification EPM_TestDate Loco_TestDate Loco_Batch Loco_Rack Loco_Box DOD          STGT_experience SequencingBatch HPC_Dissection_Date_Collapsed HPC_Dissector_DissectionDate
# 63  "12/1/2020"         "4"    "Male"   "IN"                   "10/21/2013" "9/24/2013"   "1"        "2"       "10"     "12/10/2013" "TRUE"          "1"             "Batch1"                      "KA_12/1/2020"              
# 181 "12/1/2020"         "9"    "Male"   "ST"                   "10/25/2013" "9/25/2013"   "3"        "1"       "13"     "12/11/2013" "TRUE"          "1"             "Batch1"                      "KA_12/1/2020"              
# 185 "12/1/2020"         "9"    "Male"   "ST"                   "10/22/2013" "9/25/2013"   "3"        "1"       "17"     "12/11/2013" "TRUE"          "1"             "Batch1"                      "KA_12/1/2020"              
# 183 "12/1/2020"         "9"    "Male"   "IN"                   "10/25/2013" "9/25/2013"   "3"        "1"       "15"     "12/11/2013" "TRUE"          "1"             "Batch1"                      "KA_12/1/2020"              
# 83  "12/1/2020"         "4"    "Male"   "ST"                   "10/21/2013" "9/24/2013"   "1"        "1"       "15"     "12/10/2013" "TRUE"          "1"             "Batch1"                      "KA_12/1/2020"              
# 5   "12/1/2020"         "1"    "Female" "ST"                   "10/30/2013" "9/26/2013"   "2"        "1"       " 4"     "12/18/2013" "TRUE"          "1"             "Batch1"                      "KA_12/1/2020"              
# HPC_Dissector_DissectionDate_Collapsed
# 63  "Batch6"                              
# 181 "Batch6"                              
# 185 "Batch6"                              
# 183 "Batch6"                              
# 83  "Batch6"                              
# 5   "Batch6"                              

#Outputting histograms for each continuous variable:
for (i in 1:length(SubjectContinuousVariables[1,])){
  png(paste(paste("Histogram of", colnames(SubjectContinuousVariables)[i], sep="  "), "png", sep=".")) 
  hist(SubjectContinuousVariables[, i], col=i+1)
  dev.off() 
}

#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables:
for (i in 1:length(SubjectContinuousVariables[1,])){
  for(j in 1:length(SubjectContinuousVariables[1,])){
    png(paste("14", paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), "png", sep=".")) 
    plot(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j], main=paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), xlab=colnames(SubjectContinuousVariables)[j], ylab=colnames(SubjectContinuousVariables)[i])
    RegressionLine<-lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j])
    abline(RegressionLine, col=2)
    mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4)))
    dev.off() 
  } 
}

#Using boxplots to visually examine the relationships between the continuous subject variables and categorical subject variables:
for (i in 1:length(SubjectContinuousVariables[1,])){
  for(j in 1:length(SubjectFactorVariables[1,])){
    png(paste("14", paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), "png", sep=".")) 
    boxplot(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j], main=paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), xlab=colnames(SubjectFactorVariables)[j], ylab=colnames(SubjectContinuousVariables)[i])
    mtext(paste("p-value=", round(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4)))
    dev.off() 
  } 
}

#Creating a text file of contingency tables to visually examine the relationships between categorical subject variables:
CrossTabsIV<-file("Cross Tabs Between Subject Factors.txt")
out<-c(
  capture.output(
    for (i in 1:length(SubjectFactorVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        ContingencyTable<-table(SubjectFactorVariables[,i],SubjectFactorVariables[,j])
        print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(ContingencyTable)
        print(paste("p-value=", chisq.test(ContingencyTable)$p.value)) 
      } 
    }
  )
)
cat(out, file="Cross Tabs Between Subject Factors.txt", sep="\n", append=TRUE)
close(CrossTabsIV)
rm(out)

library(car)

StatisticalRelationshipsIV<-file("Statistical Relationships between Subject Variables.txt")
out<-c(
  
  #Using linear regression to examine the statistical relationships between the continuous subject variables:
  
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "))
        print(summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j])))
      } 
    }
  ),
  
  #Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
  
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))) 
      } 
    }
  ),
  
  #Using chi-square to examine the statistical relationships between the categorical subject variables:
  
  capture.output(
    for (i in 1:length(SubjectFactorVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(chisq.test(ContingencyTable)) 
      } 
    }
  )
  
)
cat(out, file="Statistical Relationships between Subject Variables.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIV)
rm(out)

#Flagging variables that are collinear with other subject variables:

FlaggedRelationshipsBetweenIV<-file("Flagged Relationships Between Subject Variables.txt")
out<-c(
  
  #Using linear regression to examine the statistical relationships between the continuous subject variables:
  
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        if(summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j]))$coefficient[8]<0.05){
          print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], "p-value=", summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j]))$coefficient[8], sep="  "))}else{}
      } 
    }
  ),
  
  #Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
  
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        if(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
          print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], sep="  ")) 
        }else{} 
      } 
    }
  ),
  
  #Using chi-square to examine the statistical relationships between the categorical subject variables:
  
  capture.output(
    for (i in 1:length(SubjectFactorVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        ContingencyTable<-table(SubjectFactorVariables[,i], SubjectFactorVariables[,j])
        if(chisq.test(ContingencyTable)$p.value<0.05){
          print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", chisq.test(ContingencyTable)$p.value, sep="  "))
        }else{}
      }
    }
  )
)
cat(out, file="Flagged Relationships Between Subject Variables.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIV)
rm(out)

########2/24/2022 Elaine creating higher quality plots

colnames(F2_PCA_wMetaData)

pdf("F2_PC1vsRibosomePerc.pdf", width=5, height=5)
plot(PC1~RibosomePerc, data=F2_PCA_wMetaData, pch=18, ylab="PC1", xlab="%rRNA")
BestFitLine<-lm(PC1~RibosomePerc, data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()

#stats for figure legend:
summary.lm(lm(PC1~RibosomePerc, data=F2_PCA_wMetaData))
# Call:
#   lm(formula = PC1 ~ RibosomePerc, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -64.868 -14.686  -0.873  15.140  58.871 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -100.961      3.891  -25.95   <2e-16 ***
#   RibosomePerc 8893.192    319.447   27.84   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 22.06 on 243 degrees of freedom
# Multiple R-squared:  0.7613,	Adjusted R-squared:  0.7603 
# F-statistic:   775 on 1 and 243 DF,  p-value: < 2.2e-16


pdf("F2_PC1vsPercentIntergenic.pdf", width=5, height=5)
plot(PC1~Percent.Intergenic, data=F2_PCA_wMetaData, pch=18, ylab="PC1", xlab="%Intergenic")
BestFitLine<-lm(PC1~Percent.Intergenic, data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(PC1~Percent.Intergenic, data=F2_PCA_wMetaData))
# Call:
#   lm(formula = PC1 ~ Percent.Intergenic, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -32.213  -9.456  -0.042   9.502  47.256 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          336.218      7.339   45.81   <2e-16 ***
#   Percent.Intergenic -2514.858     54.456  -46.18   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 14.44 on 243 degrees of freedom
# Multiple R-squared:  0.8977,	Adjusted R-squared:  0.8973 
# F-statistic:  2133 on 1 and 243 DF,  p-value: < 2.2e-16


pdf("F2_PC2vsRibosomePerc.pdf", width=5, height=5)
plot(PC2~RibosomePerc, data=F2_PCA_wMetaData, pch=18, ylab="PC2", xlab="%rRNA")
BestFitLine<-lm(PC2~RibosomePerc, data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(PC2~RibosomePerc, data=F2_PCA_wMetaData))
# Call:
#   lm(formula = PC2 ~ RibosomePerc, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -37.084  -7.634   1.364   9.168  36.355 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)     5.637      2.350   2.398   0.0172 *
#   RibosomePerc -496.506    192.975  -2.573   0.0107 *
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 13.33 on 243 degrees of freedom
# Multiple R-squared:  0.02652,	Adjusted R-squared:  0.02251 
# F-statistic:  6.62 on 1 and 243 DF,  p-value: 0.01068

sqrt(0.02652)

pdf("F2_PC2vsPercentIntergenic.pdf", width=5, height=5)
plot(PC2~Percent.Intergenic, data=F2_PCA_wMetaData, pch=18, ylab="PC2", xlab="%Intergenic")
BestFitLine<-lm(PC2~Percent.Intergenic, data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(PC2~Percent.Intergenic, data=F2_PCA_wMetaData))

# Call:
#   lm(formula = PC2 ~ Percent.Intergenic, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -35.690  -8.343   1.062   8.552  39.224 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)           5.381      6.855   0.785    0.433
# Percent.Intergenic  -40.249     50.868  -0.791    0.430
# 
# Residual standard error: 13.49 on 243 degrees of freedom
# Multiple R-squared:  0.00257,	Adjusted R-squared:  -0.001535 
# F-statistic: 0.6261 on 1 and 243 DF,  p-value: 0.4296


pdf("F2_PC3vsSex.pdf", width=5, height=5)
boxplot(PC3~Sex, data=F2_PCA_wMetaData, xlab="Sex", ylab="PC3")
stripchart(PC3~Sex, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(PC3~Sex, data=F2_PCA_wMetaData))
# Call:
#   lm(formula = PC3 ~ Sex, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -16.3977  -2.1893   0.2403   2.6234   9.9046 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -10.5762     0.3857  -27.42   <2e-16 ***
#   SexMale      21.2392     0.5466   38.85   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 4.278 on 243 degrees of freedom
# Multiple R-squared:  0.8614,	Adjusted R-squared:  0.8608 
# F-statistic:  1510 on 1 and 243 DF,  p-value: < 2.2e-16


pdf("F2_PercentIntronic_vs_Total_Locoscore.pdf", width=5, height=5)
plot(Percent.Intronic~Total_LocoScore, data=F2_PCA_wMetaData, pch=18, ylab="%Intronic", main="F2", xlab="Total LocoScore")
BestFitLine<-lm(Percent.Intronic~Total_LocoScore, data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(Percent.Intronic~Total_LocoScore, data=F2_PCA_wMetaData))

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


pdf("F2_DV200_vs_EPM_PercTimeOA.pdf", width=5, height=5)
plot(DV200~EPM_Percent_Time_Open_Arm, data=F2_PCA_wMetaData, pch=18, ylab="DV200", main="F2", xlab="EPM: %Time in Open Arms")
BestFitLine<-lm(DV200~EPM_Percent_Time_Open_Arm, data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(DV200~EPM_Percent_Time_Open_Arm, data=F2_PCA_wMetaData))

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



