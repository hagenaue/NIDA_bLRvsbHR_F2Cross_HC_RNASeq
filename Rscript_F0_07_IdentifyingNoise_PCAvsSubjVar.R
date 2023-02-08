#F0 HC RNA-Seq Dataset
#07_Identifying sources of noise: PCA vs. Subject Variables
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-10, updated later for a few figures for the paper.



#Code that automatically loops through all variables (including the PCA output) and screens for relationships.
#This helps identify potential confounds or mediating variables (variables with relationships with lineage) 
#It also helps identify potential sources of large-scale noise (variables with relationships with the top PCs)

str(F0_PCA_wMetaData)

colnames(F0_PCA_wMetaData)

SubjectContinuousVariables<-as.matrix(F0_PCA_wMetaData[,c(5:23,26)])

SubjectFactorVariables<-as.matrix(F0_PCA_wMetaData[,c(24,25)])


str(SubjectFactorVariables)
#  chr [1:23, 1:2] "Female" "Male" "Male" "Female" "Female" "Male" "Male" "Female" "Male" "Female" "Male" "Female" "Male" "Male" "Female" "Female" "Male" "Female" "Male" "Male" "Female" "Male" "Female" "bLR" "bLR" "bHR" "bHR" "bLR" "bLR" "bHR" ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:23] "1" "2" "3" "4" ...
# ..$ : chr [1:2] "Sex_Factor" "Lineage_AsFactor"

F0_PCA_cormatrix<-cor(SubjectContinuousVariables)
write.csv(F0_PCA_cormatrix,"F0_PCA_cormatrix.csv")


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

colnames(F0_PCA_wMetaData)

pdf("F0_PC1vsRibosomePerc.pdf", width=5, height=5)
plot(PC1~RibosomePerc, data=F0_PCA_wMetaData, pch=18, ylab="PC1", xlab="%rRNA")
BestFitLine<-lm(PC1~RibosomePerc, data=F0_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()

pdf("F0_PC1vsPercentIntergenic.pdf", width=5, height=5)
plot(PC1~Percent_Intergenic, data=F0_PCA_wMetaData, pch=18, ylab="PC1", xlab="%Intergenic")
BestFitLine<-lm(PC1~Percent_Intergenic, data=F0_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()


pdf("F0_PC2vsRibosomePerc.pdf", width=5, height=5)
plot(PC2~RibosomePerc, data=F0_PCA_wMetaData, pch=18, ylab="PC2", xlab="%rRNA")
BestFitLine<-lm(PC2~RibosomePerc, data=F0_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()

pdf("F0_PC2vsPercentIntergenic.pdf", width=5, height=5)
plot(PC2~Percent_Intergenic, data=F0_PCA_wMetaData, pch=18, ylab="PC2", xlab="%Intergenic")
BestFitLine<-lm(PC2~Percent_Intergenic, data=F0_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()


pdf("F0_PC3vsSex.pdf", width=5, height=5)
boxplot(PC3~Sex, data=F0_PCA_wMetaData, xlab="Sex", ylab="PC3")
stripchart(PC3~Sex, data=F0_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()








####################


#Examining some of the relationships more closely in an effort to decide on best co-variates for our model:




#Double-checking whether Lineage is related to any of our potential co-variates 
#... which may provide interesting biological insight in their own right.


pdf("F0_RNAConc_vs_Lineage.pdf", height=6, width=4)
boxplot(Conc~Lineage_AsFactor, data=F0_PCA_wMetaData, main="F0", ylab="RNA Concentration (ng/uL)")
stripchart(Conc~Lineage_AsFactor, data=F0_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(Conc~Lineage_AsFactor, data=F0_PCA_wMetaData))

# Call:
#   lm(formula = Conc ~ Lineage_AsFactor, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -29.854 -21.725  -7.174   5.301  85.299 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           57.060      9.451   6.038 5.43e-06 ***
#   Lineage_AsFactorbLR   -9.951     13.084  -0.761    0.455    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 31.34 on 21 degrees of freedom
# Multiple R-squared:  0.02681,	Adjusted R-squared:  -0.01953 
# F-statistic: 0.5785 on 1 and 21 DF,  p-value: 0.4554


pdf("F0_mRNAPercent_vs_Lineage.pdf", height=6, width=4)
boxplot(mRNAperc~Lineage_AsFactor, data=F0_PCA_wMetaData, main="F0", ylab="mRNA %")
stripchart(mRNAperc~Lineage_AsFactor, data=F0_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(mRNAperc~Lineage_AsFactor, data=F0_PCA_wMetaData))

# Call:
#   lm(formula = mRNAperc ~ Lineage_AsFactor, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.036625 -0.006242  0.003964  0.007647  0.028815 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.844659   0.005104 165.498   <2e-16 ***
#   Lineage_AsFactorbLR -0.012073   0.007066  -1.709    0.102    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.01693 on 21 degrees of freedom
# Multiple R-squared:  0.1221,	Adjusted R-squared:  0.08026 
# F-statistic:  2.92 on 1 and 21 DF,  p-value: 0.1022


pdf("F0_RibosomePerc_vs_Lineage.pdf", height=6, width=4)
boxplot(RibosomePerc~Lineage_AsFactor, data=F0_PCA_wMetaData, main="F0", ylab="rRNA %")
stripchart(RibosomePerc~Lineage_AsFactor, data=F0_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(RibosomePerc~Lineage_AsFactor, data=F0_PCA_wMetaData))

# Call:
#   lm(formula = RibosomePerc ~ Lineage_AsFactor, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0104442 -0.0050758 -0.0001522  0.0033205  0.0156868 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.022489   0.002343   9.597 3.96e-09 ***
#   Lineage_AsFactorbLR 0.001585   0.003244   0.488     0.63    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.007772 on 21 degrees of freedom
# Multiple R-squared:  0.01123,	Adjusted R-squared:  -0.03585 
# F-statistic: 0.2386 on 1 and 21 DF,  p-value: 0.6303



pdf("F0_Reads_vs_Lineage.pdf", height=6, width=4)
boxplot(PF_Reads~Lineage_AsFactor, data=F0_PCA_wMetaData, main="F0", ylab="PF.Reads")
stripchart(PF_Reads~Lineage_AsFactor, data=F0_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()


summary.lm(lm(PF_Reads~Lineage_AsFactor, data=F0_PCA_wMetaData))

# Call:
#   lm(formula = PF_Reads ~ Lineage_AsFactor, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -5780703 -2239428   215387  2225513  6177761 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         33403351    1026906  32.528   <2e-16 ***
#   Lineage_AsFactorbLR -3021536    1421687  -2.125   0.0456 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 3406000 on 21 degrees of freedom
# Multiple R-squared:  0.177,	Adjusted R-squared:  0.1378 
# F-statistic: 4.517 on 1 and 21 DF,  p-value: 0.04559


#Interesting - bLRs have fewer total reads (post-filter/QC) than bHRs
#This is useful because all of the RNA-Seq data is typically normalized to total reads (e.g., CPM - counts per million)
#So... if you have bLRs looking like they have elevated levels of something - e.g.C1qc - that is relative to their total reads. If their total reads are less, if you ran ISH, you might not see the same result.
#This means that interpreting anything down-regulated in bLRs is easy
#Interpreting upregulation in bLRs is trickier, probably best to prioritize the most extreme results - less extreme may not validate using methods without the same normalization


pdf("F0_Age_vs_Lineage.pdf", height=6, width=4)
boxplot(Age_days~Lineage_AsFactor, data=F0_PCA_wMetaData, main="F0", ylab="Age (days)")
stripchart(Age_days~Lineage_AsFactor, data=F0_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()


summary.lm(lm(Age_days~Lineage_AsFactor, data=F0_PCA_wMetaData))
# Call:
#   lm(formula = Age_days ~ Lineage_AsFactor, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -24.17 -12.17 -10.17  14.09  16.83 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         175.9091     4.4277  39.730   <2e-16 ***
#   Lineage_AsFactorbLR  -0.7424     6.1298  -0.121    0.905    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 14.68 on 21 degrees of freedom
# Multiple R-squared:  0.000698,	Adjusted R-squared:  -0.04689 
# F-statistic: 0.01467 on 1 and 21 DF,  p-value: 0.9048


#######################################################

#Checking whether sex is related to any of our potential co-variates
#... which may provide biological insight in its own right (but also could easily be technical in origin)

pdf("F0_RNAConc_vs_Sex.pdf", height=6, width=4)
boxplot(Conc~Sex_Factor, data=F0_PCA_wMetaData, main="F0", ylab="RNA Concentration (ng/uL)")
stripchart(Conc~Sex_Factor, data=F0_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(Conc~Sex_Factor, data=F0_PCA_wMetaData))

# Call:
#   lm(formula = Conc ~ Sex_Factor, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -37.842 -22.741  -4.891   8.007  74.582 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        43.754      8.807   4.968 6.47e-05 ***
#   Sex_FactorFemale   16.965     12.735   1.332    0.197    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 30.51 on 21 degrees of freedom
# Multiple R-squared:  0.07792,	Adjusted R-squared:  0.03401 
# F-statistic: 1.775 on 1 and 21 DF,  p-value: 0.1971


pdf("F0_PercentRNA_vs_Sex.pdf", height=6, width=4)
boxplot(mRNAperc~Sex_Factor, data=F0_PCA_wMetaData, main="F0", ylab="Percent mRNA")
stripchart(mRNAperc~Sex_Factor, data=F0_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(mRNAperc~Sex_Factor, data=F0_PCA_wMetaData))
# 
# Call:
#   lm(formula = mRNAperc ~ Sex_Factor, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.044116 -0.006885  0.003393  0.011750  0.022718 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.836786   0.005191 161.188   <2e-16 ***
#   Sex_FactorFemale 0.003291   0.007507   0.438    0.666    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.01798 on 21 degrees of freedom
# Multiple R-squared:  0.009072,	Adjusted R-squared:  -0.03812 
# F-statistic: 0.1923 on 1 and 21 DF,  p-value: 0.6655

pdf("F0_RibosomePerc_vs_Sex.pdf", height=6, width=4)
boxplot(RibosomePerc~Sex_Factor, data=F0_PCA_wMetaData, main="F0", ylab="Percent Ribosome")
stripchart(RibosomePerc~Sex_Factor, data=F0_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(RibosomePerc~Sex_Factor, data=F0_PCA_wMetaData))

# Call:
#   lm(formula = RibosomePerc ~ Sex_Factor, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0094093 -0.0057688  0.0001907  0.0024462  0.0163317 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       0.024812   0.002206   11.25 2.39e-10 ***
#   Sex_FactorFemale -0.003128   0.003190   -0.98    0.338    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.007643 on 21 degrees of freedom
# Multiple R-squared:  0.04377,	Adjusted R-squared:  -0.001762 
# F-statistic: 0.9613 on 1 and 21 DF,  p-value: 0.338

pdf("F0_Reads_vs_Sex.pdf", height=6, width=4)
boxplot(PF_Reads~Sex_Factor, data=F0_PCA_wMetaData, main="F0", ylab="PF.Reads")
stripchart(PF_Reads~Sex_Factor, data=F0_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(PF_Reads~Sex_Factor, data=F0_PCA_wMetaData))

# Call:
#   lm(formula = PF_Reads ~ Sex_Factor, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -7729484 -2298671  -767971  2606908  8215938 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      31365174    1073943  29.206   <2e-16 ***
#   Sex_FactorFemale   965424    1552919   0.622    0.541    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 3720000 on 21 degrees of freedom
# Multiple R-squared:  0.01807,	Adjusted R-squared:  -0.02869 
# F-statistic: 0.3865 on 1 and 21 DF,  p-value: 0.5408

##################################

#Exploring the relationship between the top PCs and some of the co-variates more carefully:

summary.lm(lm(PC1~RibosomePerc, data=F0_PCA_wMetaData))

# Call:
#   lm(formula = PC1 ~ RibosomePerc, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -60.295 -27.595   0.723  22.232  47.369 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    133.60      20.29   6.584 1.61e-06 ***
#   RibosomePerc -5729.85     828.77  -6.914 7.86e-07 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 29.68 on 21 degrees of freedom
# Multiple R-squared:  0.6948,	Adjusted R-squared:  0.6802 
# F-statistic:  47.8 on 1 and 21 DF,  p-value: 7.86e-07

summary.lm(lm(PC1~RibosomePerc+Conc, data=F0_PCA_wMetaData))
# Call:
#   lm(formula = PC1 ~ RibosomePerc + Conc, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -47.838 -15.941  -0.414  20.988  39.041 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     51.6614    32.8839   1.571 0.131864    
# RibosomePerc -3778.1196   972.7084  -3.884 0.000922 ***
#   Conc             0.7024     0.2393   2.935 0.008182 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 25.43 on 20 degrees of freedom
# Multiple R-squared:  0.7867,	Adjusted R-squared:  0.7653 
# F-statistic: 36.88 on 2 and 20 DF,  p-value: 1.952e-07

summary.lm(lm(PC1~RibosomePerc+Percent_Intergenic, data=F0_PCA_wMetaData))

# 
# Call:
#   lm(formula = PC1 ~ RibosomePerc + Percent_Intergenic, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -47.873 -13.402  -1.273  18.308  26.354 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          -45.08      37.39  -1.206    0.242    
# RibosomePerc       -5670.98     557.91 -10.165 2.40e-09 ***
#   Percent_Intergenic  1727.79     336.53   5.134 5.05e-05 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 19.98 on 20 degrees of freedom
# Multiple R-squared:  0.8683,	Adjusted R-squared:  0.8552 
# F-statistic: 65.94 on 2 and 20 DF,  p-value: 1.568e-09

summary.lm(lm(PC2~Percent_UTR, data=F0_PCA_wMetaData))
# 
# Call:
#   lm(formula = PC2 ~ Percent_UTR, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -18.796  -5.594  -1.127   7.747  13.775 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   386.94      30.69   12.61 2.92e-11 ***
#   Percent_UTR -1072.13      84.88  -12.63 2.81e-11 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 9.066 on 21 degrees of freedom
# Multiple R-squared:  0.8837,	Adjusted R-squared:  0.8781 
# F-statistic: 159.5 on 1 and 21 DF,  p-value: 2.814e-11

summary.lm(lm(Percent_UTR~Lineage_AsFactor, data=F0_PCA_wMetaData))
# 
# Call:
#   lm(formula = Percent_UTR ~ Lineage_AsFactor, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.038783 -0.011254  0.000073  0.013019  0.038300 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.348554   0.005955  58.534  < 2e-16 ***
#   Lineage_AsFactorbLR 0.023680   0.008244   2.872  0.00912 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.01975 on 21 degrees of freedom
# Multiple R-squared:  0.2821,	Adjusted R-squared:  0.2479 
# F-statistic: 8.251 on 1 and 21 DF,  p-value: 0.009116

#Percent UTR is closely related to PC2, but is also related to Lineage. :(  That makes it debatable to use as a co-variate.

#Pecent intergenic had the strongest relationship with PC1 in the F2 dataset. How about here?

summary.lm(lm(PC2~Percent_Intergenic, data=F0_PCA_wMetaData))
# 
# Call:
#   lm(formula = PC2 ~ Percent_Intergenic, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -36.096  -9.262   6.164  12.626  28.713 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          155.55      31.18   4.988 6.16e-05 ***
#   Percent_Intergenic -1515.81     301.70  -5.024 5.66e-05 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 17.91 on 21 degrees of freedom
# Multiple R-squared:  0.5459,	Adjusted R-squared:  0.5243 
# F-statistic: 25.24 on 1 and 21 DF,  p-value: 5.66e-05
# 

summary.lm(lm(Percent_Intergenic~Lineage_AsFactor, data=F0_PCA_wMetaData))
# 
# 
# Call:
#   lm(formula = Percent_Intergenic ~ Lineage_AsFactor, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.017295 -0.008933 -0.002262  0.003323  0.036935 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.098341   0.003687  26.673   <2e-16 ***
#   Lineage_AsFactorbLR 0.008197   0.005104   1.606    0.123    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.01223 on 21 degrees of freedom
# Multiple R-squared:  0.1094,	Adjusted R-squared:  0.06697 
# F-statistic: 2.579 on 1 and 21 DF,  p-value: 0.1232


#... and percent intergenic is less related to bHR/bLR lineage. Better co-variate?

