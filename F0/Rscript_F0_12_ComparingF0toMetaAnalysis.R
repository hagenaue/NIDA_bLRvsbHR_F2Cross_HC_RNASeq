#F0 HC RNA-Seq Dataset
#12_Comparing F0 Results to the Results of our late generation bHR/bLR Meta-Analysis (Birt et al. 2021)
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-10, updated later for a few figures for the paper.


#Comparing our results with results from the Late Generation bHR/bLR meta-analysis of hippocampal RNA-Seq data (from Birt Hagenauer et al. 2021)

MetaAnalysisOutput_JustLateGen<-read.csv("metaOutput_JustLateGen_FDR.csv", header=TRUE, stringsAsFactors = FALSE)
str(MetaAnalysisOutput_JustLateGen)

# 'data.frame':	13319 obs. of  11 variables:
#   $ X         : int  19 20 21 22 23 24 25 26 27 28 ...
# $ GeneSymbol: chr  "A2m" "A3galt2" "Aaas" "Aacs" ...
# $ estimate  : num  0.482 -0.486 0.817 -1.048 -0.436 ...
# $ SE        : num  0.438 0.441 0.443 0.467 0.457 ...
# $ pval      : num  0.2716 0.2713 0.0649 0.0249 0.3399 ...
# $ CI_lb     : num  -0.3773 -1.3509 -0.0504 -1.9638 -1.3312 ...
# $ CI_ub     : num  1.341 0.38 1.685 -0.132 0.459 ...
# $ datasets  : int  2 2 2 2 2 2 2 2 2 2 ...
# $ rawp      : num  0.2716 0.2713 0.0649 0.0249 0.3399 ...
# $ BH        : num  0.599 0.599 0.312 0.197 0.659 ...
# $ BY        : num  1 1 1 1 1 1 1 1 1 1 ..


colnames(MetaAnalysisOutput_JustLateGen)[2]<-"SYMBOL"


F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen<-join(F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped, MetaAnalysisOutput_JustLateGen, by="SYMBOL", type="left")

str(F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen)

#re-run on 11/23/2021
#'data.frame':	13860 obs. of  36 variables:



plot(Coef.Lineage_AsFactorbLR~estimate, data=F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen)
#Ah - the estimates are reversed (bLR = reference group) - I just double-checked C1qa and Bmp4 in the output

F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$estimate<-(F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$estimate)*-1

F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$CI_lb<-(F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$CI_lb)*-1

F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$CI_ub<-(F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$CI_ub)*-1


write.csv(F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen, "F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen.csv")

pdf("U01RNASeq_Vs_LateGeneMetaAnalysis_bLRvsbHR.pdf", width=6, height=6)
plot(Coef.Lineage_AsFactorbLR~estimate, data=F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen, ylab="U01 RNA-Seq: bLR vs. bHR Log2FC", xlab="Meta-Analysis of F37 & F43 RNA-Seq: bLR vs. bHR Estimated Effect Size", ylim=c(-5,5))
BestFitLine<-lm(Coef.Lineage_AsFactorbLR~estimate, data=F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen)
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(BestFitLine)
# 
# Call:
#   lm(formula = Coef.Lineage_AsFactorbLR ~ estimate, data = F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.9166 -0.0970 -0.0061  0.0890  2.7139 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.012034   0.001782   6.753 1.53e-11 ***
#   estimate    0.095302   0.002449  38.911  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1881 on 11173 degrees of freedom
# (2611 observations deleted due to missingness)
# Multiple R-squared:  0.1193,	Adjusted R-squared:  0.1193 
# F-statistic:  1514 on 1 and 11173 DF,  p-value: < 2.2e-16



pdf("U01RNASeq_Vs_LateGeneMetaAnalysis_bLRvsbHR_Tstat.pdf", width=6, height=6)
plot(t.Lineage_AsFactorbLR~estimate, data=F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen, ylab="U01 RNA-Seq: bLR vs. bHR Tstat", xlab="Meta-Analysis of F37 & F43 RNA-Seq: bLR vs. bHR Estimated Effect Size", ylim=c(-12,12))
BestFitLine<-lm(t.Lineage_AsFactorbLR~estimate, data=F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen)
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(BestFitLine)
# 

# Call:
#   lm(formula = t.Lineage_AsFactorbLR ~ estimate, data = F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -6.1636 -0.7653 -0.0254  0.7533  7.8661 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.07011    0.01122   6.246 4.36e-10 ***
#   estimate     0.66108    0.01543  42.854  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 1.185 on 11173 degrees of freedom
# (2611 observations deleted due to missingness)
# Multiple R-squared:  0.1412,	Adjusted R-squared:  0.1411 
# F-statistic:  1836 on 1 and 11173 DF,  p-value: < 2.2e-16


pdf("U01RNASeq_Vs_LateGeneMetaAnalysis_bLRvsbHR_FDR10Meta.pdf", width=6, height=6)
plot(Coef.Lineage_AsFactorbLR~estimate, data=F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen[F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH<0.1,], ylab="U01 RNA-Seq: bLR vs. bHR Log2FC", xlab="Meta-Analysis of F37 & F43 RNA-Seq: bLR vs. bHR Estimated Effect Size", ylim=c(-5,5))
BestFitLine<-lm(Coef.Lineage_AsFactorbLR~estimate, data=F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen[F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH<0.1,])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(BestFitLine)
# 
# Call:
#   lm(formula = Coef.Lineage_AsFactorbLR ~ estimate, data = F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen[F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH < 
#                                                                                                                                                   0.1, ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.8064 -0.1320  0.0024  0.1372  2.6590 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.012763   0.012827  -0.995     0.32    
# estimate     0.124421   0.006933  17.946   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.3519 on 765 degrees of freedom
# (2611 observations deleted due to missingness)
# Multiple R-squared:  0.2963,	Adjusted R-squared:  0.2953 
# F-statistic: 322.1 on 1 and 765 DF,  p-value: < 2.2e-16


pdf("U01RNASeq_Vs_LateGeneMetaAnalysis_bLRvsbHR_Tstat_FDR10Meta.pdf", width=6, height=6)
plot(t.Lineage_AsFactorbLR~estimate, data=F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen[F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH<0.1,], ylab="U01 RNA-Seq: bLR vs. bHR Tstat", xlab="Meta-Analysis of F37 & F43 RNA-Seq: bLR vs. bHR Estimated Effect Size", ylim=c(-12,12))
BestFitLine<-lm(t.Lineage_AsFactorbLR~estimate, data=F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen[F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH<0.1,])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(BestFitLine)
# 
# Call:
#   lm(formula = t.Lineage_AsFactorbLR ~ estimate, data = F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen[F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH < 
#                                                                                                                                                0.1, ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.9095 -0.8927 -0.0432  0.8915  7.2547 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.01649    0.05278  -0.313    0.755    
# estimate     0.78954    0.02853  27.678   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 1.448 on 765 degrees of freedom
# (2611 observations deleted due to missingness)
# Multiple R-squared:  0.5003,	Adjusted R-squared:  0.4997 
# F-statistic: 766.1 on 1 and 765 DF,  p-value: < 2.2e-16

