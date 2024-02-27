#F0 HC RNA-Seq Dataset
#03_Relationships Between MetaData Variables
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-10, updated later for a few figures for the paper.



pdf("F0_LibrarySize_vs_Lineage.pdf", height=6, width=4)
boxplot(LibrarySize~Lineage_AsFactor, data=F0_PCA_wMetaData, main="F0", ylab="Library Size")
stripchart(LibrarySize~Lineage_AsFactor, data=F0_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(LibrarySize~Lineage_AsFactor, data=F0_PCA_wMetaData))

# 
# Call:
#   lm(formula = LibrarySize ~ Lineage_AsFactor, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -4301271 -1242693  -330577  1757161  4539869 
# 
# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         23903180     707556  33.783  < 2e-16 ***
#   Lineage_AsFactorbLR -2777010     979567  -2.835  0.00992 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 2347000 on 21 degrees of freedom
# Multiple R-squared:  0.2768,	Adjusted R-squared:  0.2423 
# F-statistic: 8.037 on 1 and 21 DF,  p-value: 0.00992

#Similar to pf.reads - library sizes varies with Lineage.
#May be a confound worth considering... or a genuine biological signal??

pdf("F0_LibrarySize_vs_Sex.pdf", height=6, width=4)
boxplot(LibrarySize~Sex_Factor, data=F0_PCA_wMetaData, main="F0", ylab="Library Size")
stripchart(LibrarySize~Sex_Factor, data=F0_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(LibrarySize~Sex_Factor, data=F0_PCA_wMetaData))

# 
# Call:
#   lm(formula = LibrarySize ~ Sex_Factor, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -5972900 -2023103  -347151  1892918  6303613 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      22139436     790364  28.012   <2e-16 ***
#   Sex_FactorFemale   658364    1142864   0.576    0.571    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 2738000 on 21 degrees of freedom
# Multiple R-squared:  0.01556,	Adjusted R-squared:  -0.03132 
# F-statistic: 0.3319 on 1 and 21 DF,  p-value: 0.5707

####################################################1/27/2022 Elaine

colnames(F0_PCA_wMetaData)


pdf("F0_Percent_UTR_vs_Lineage.pdf", height=6, width=4)
boxplot(Percent_UTR~Lineage_AsFactor, data=F0_PCA_wMetaData, main="F0", ylab="%UTR")
stripchart(Percent_UTR~Lineage_AsFactor, data=F0_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(Percent_UTR~Lineage_AsFactor, data=F0_PCA_wMetaData))

# Call:
#   lm(formula = Percent_UTR ~ Lineage_AsFactor, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.038783 -0.011254  0.000073  0.013019  0.038300 
# 
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           0.348554   0.005955  58.534  < 2e-16 ***
#   Lineage_AsFactorbLR 0.023680   0.008244   2.872  0.00912 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.01975 on 21 degrees of freedom
# Multiple R-squared:  0.2821,	Adjusted R-squared:  0.2479 
# F-statistic: 8.251 on 1 and 21 DF,  p-value: 0.009116



pdf("F0_Percent_Coding_vs_Lineage.pdf", height=6, width=4)
boxplot(Percent_Coding~Lineage_AsFactor, data=F0_PCA_wMetaData, main="F0", ylab="%Coding")
stripchart(Percent_Coding~Lineage_AsFactor, data=F0_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(Percent_Coding~Lineage_AsFactor, data=F0_PCA_wMetaData))


# Call:
#   lm(formula = Percent_Coding ~ Lineage_AsFactor, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.07492 -0.01934  0.00552  0.02077  0.06636 
# 
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           0.49610    0.01039  47.727   <2e-16 ***
#   Lineage_AsFactorbLR -0.03575    0.01439  -2.484   0.0215 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.03448 on 21 degrees of freedom
# Multiple R-squared:  0.2272,	Adjusted R-squared:  0.1904 
# F-statistic: 6.172 on 1 and 21 DF,  p-value: 0.02149


pdf("F0_RIN_vs_Lineage.pdf", height=6, width=4)
boxplot(RIN~Lineage_AsFactor, data=F0_PCA_wMetaData, main="F0", ylab="RIN")
stripchart(RIN~Lineage_AsFactor, data=F0_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()

summary.lm(lm(RIN~Lineage_AsFactor, data=F0_PCA_wMetaData))

# Call:
#   lm(formula = RIN ~ Lineage_AsFactor, data = F0_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.6250 -0.1409  0.0750  0.1920  0.4750 
# 
# Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            9.59091    0.08658 110.775   <2e-16 ***
#   Lineage_AsFactorbLR -0.26591    0.11986  -2.218   0.0377 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.2872 on 21 degrees of freedom
# Multiple R-squared:  0.1899,	Adjusted R-squared:  0.1513 
# F-statistic: 4.921 on 1 and 21 DF,  p-value: 0.03768


pdf("F0_PF_Reads_vs_Lineage.pdf", height=6, width=4)
boxplot(PF_Reads~Lineage_AsFactor, data=F0_PCA_wMetaData, main="F0", ylab="PF_Reads")
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
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            33403351    1026906  32.528   <2e-16 ***
#   Lineage_AsFactorbLR -3021536    1421687  -2.125   0.0456 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 3406000 on 21 degrees of freedom
# Multiple R-squared:  0.177,	Adjusted R-squared:  0.1378 
# F-statistic: 4.517 on 1 and 21 DF,  p-value: 0.04559

######################################################
