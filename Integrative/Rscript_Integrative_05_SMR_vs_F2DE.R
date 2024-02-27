
############################################

#Moving on to the SMR Co-Localization Analyses:

####################################

#What about the correlation between F2 Log2FC and F2 SMR?

plot(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Log2FC_LocoScore)~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR, xlab="SMR (T): Total LocoScore", ylab="Abs(F2 Log2FC Locoscore)")
#pretty rough
plot(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore)~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR, xlab="SMR (T): Total LocoScore", ylab="Abs(T) for F2 Log2FC Locoscore")
#pretty rough, but there's a few genes that stand out

cor(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore), F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR, use="pairwise.complete.obs")
#[1] 0.2038793 (Pearson)
cor(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore), F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR, use="pairwise.complete.obs", method="spearman")
#[1] 0.1187908 (Spearman)

#I dug around, and apparently we can assign the SMR output a directionality (since both the QTL and the eQTL results have a direction)

#Note: this is assuming that the QTL has the same phenotype associated with "reference" as bHR/bLR - worth double-checking...
#I should be able to get that information from the GWAS Z-score (as long as the effects were relatively normally distributed...)
#and from the eQTL

#So if GWAS Z-score is positive for locoscore, that means more locoscore with the alt allele
#So if GWAS Z-score is negative for locoscore, that means less locoscore with the alt allele
#If eQTL t-stat is positive, that means more gene expression with the alt allele
#If eQTL t-stat is negative, that means less gene expression with the alt allele

#So GWAS Z>0 and eQTL T>0, then more locoscore is associated with more gene expression (or hypothetically + Log2FC)
#So GWAS Z<0 and eQTL T<0, then less locoscore is associated with less gene expression (or hypothetically also + Log2FC)
#GWAS Z>0 and eQTL T<0, then more locoscore is associated with less gene expression (or hypothetically - Log2FC)
#GWAS Z<0 and eQTL T>0, then less locoscore is associated with more gene expression (or hypothetically - Log2FC)

#Mathematically, those relationships are pretty handy
F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2<-abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR)
F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[which((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref*F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_GWAS)<0)]<-F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[which((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref*F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_GWAS)<0)]*(-1)


#Making the same plot again using this new version with predicted direction:
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2)
#ooh - much nicer.


#Making a version that uses our additional criteria for filtering and labeling:

CriteriaForPlot<-(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1 & ((bHRGenes==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR<0)|(bLRGenes==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR>0)))
sum(CriteriaForPlot, na.rm=TRUE)
#[1] 492

pdf("F2DE_Locoscore_vs_PredictedF2DE_LocoScore_AlleVariants_bHRbLRgenes_NoColor_topeGenesBlack.pdf", width=5, height=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2, col="grey", ylab="F2 DE: LocoScore T-stat", xlab="Predicted F2 DE: SMR T-stat")
Trendline1<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2)
abline(Trendline1, col="grey", lwd=3)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
Trendline2<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
abline(Trendline2, col="black", lwd=3)
dev.off()

summary.lm(Trendline1)

# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore ~ 
#        F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.4128 -0.6137  0.0121  0.6252  4.4688 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                               0.08253    0.01268   6.508 8.28e-11 ***
#   F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2  0.10245    0.00408  25.111  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.9614 on 5745 degrees of freedom
# (8606 observations deleted due to missingness)
# Multiple R-squared:  0.0989,	Adjusted R-squared:  0.09874 
# F-statistic: 630.5 on 1 and 5745 DF,  p-value: < 2.2e-16

pdf("F2DE_LocoScore_vs_PredictedF2DE_Locoscore_onlybHRbLRgenes_p05topeGenes.pdf", width=5, height=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)], ylab="F2 DE: LocoScore T-stat", xlab="Predicted F2 DE: SMR T-stat")
Trendline2<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[which(CriteriaForPlot==TRUE & bHRGenes==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_LocoScore<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_p_SMR<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2>0)]
       ~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[which(CriteriaForPlot==TRUE & bHRGenes==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_LocoScore<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_p_SMR<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2>0)], col="green4", pch=16)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[which(CriteriaForPlot==TRUE & bLRGenes==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_LocoScore<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_p_SMR<0.05  & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2<0)]
       ~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[which(CriteriaForPlot==TRUE & bLRGenes==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_LocoScore<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_p_SMR<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2<0)], col="red", pch=16)
abline(Trendline2, col="black", lwd=3)
dev.off()

summary.lm(Trendline2)

# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[which(CriteriaForPlot == 
#                                                                               TRUE)] ~ F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[which(CriteriaForPlot == 
#                                                                                                                                                                         TRUE)])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.9600 -0.6613 -0.0068  0.7213  3.4544 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                                                                                              -0.029193   0.049538  -0.589    0.556
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[which(CriteriaForPlot == TRUE)]  0.143857   0.009531  15.094   <2e-16
# 
# (Intercept)                                                                                                 
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[which(CriteriaForPlot == TRUE)] ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 1.097 on 490 degrees of freedom
# Multiple R-squared:  0.3174,	Adjusted R-squared:  0.316 
# F-statistic: 227.8 on 1 and 490 DF,  p-value: < 2.2e-16

sqrt(0.3174)
#[1] 0.5633826

#For labeling:
cbind.data.frame(Gene=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F2..Rnor6.Ensembl.v.103.[which(CriteriaForPlot==TRUE & bHRGenes==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_LocoScore<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_p_SMR<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2>0)],
                 DETstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[which(CriteriaForPlot==TRUE & bHRGenes==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_LocoScore<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_p_SMR<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2>0)],
                 SMRTstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[which(CriteriaForPlot==TRUE & bHRGenes==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_LocoScore<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_p_SMR<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2>0)])

#             Gene  DETstat  SMRTstat
# 1          Zfp551 2.421338  8.881697
# 2  AABR07071904.1 3.867393 36.444578
# 3           Mfge8 2.916811 16.891925
# 4           Lipt2 2.251408 14.777376
# 5           C2cd3 2.461613 21.405944
# 6            Ucp2 3.235008 21.737626
# 7             Sp3 2.348308  8.119622
# 8         Ttc30a1 3.909581  9.624517
# 9          Tmco5a 1.989535  4.260536
# 10          Cndp1 2.701293  6.831436
# 11           Spg7 3.789752 13.428439
# 12           Spg7 3.789752  7.542803
# 13         Afg3l1 3.484937 17.073236

cbind.data.frame(Gene=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F2..Rnor6.Ensembl.v.103.[which(CriteriaForPlot==TRUE & bLRGenes==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_LocoScore<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_p_SMR<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2<0)],
                 DETstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[which(CriteriaForPlot==TRUE & bLRGenes==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_LocoScore<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_p_SMR<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2<0)],
                 SMRTstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[which(CriteriaForPlot==TRUE & bLRGenes==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_LocoScore<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_p_SMR<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2<0)])

#       Gene   DETstat   SMRTstat
# 1      Lsr -2.434178 -17.249650
# 2    Fanci -3.007313 -17.944510
# 3   Pex11a -2.757793 -18.661885
# 4    Wdr93 -3.129942 -17.101626
# 5   Unc45a -2.175854 -14.773957
# 6  Plekhb1 -2.807592 -24.329036
# 7  Rarres2 -3.126275  -6.162816
# 8   Gimap5 -3.282441  -4.998288
# 9     Fzd6 -3.633133 -17.942523
# 10    Idh1 -4.006376  -9.857362
# 11    Ist1 -2.017224 -13.666533
# 12  Vps9d1 -3.650486 -15.552153


######################################

#What do these predicted vs. real DE plots look like for some of the other behaviors?

F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2<-abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR)
F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2[which((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref*F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_z_GWAS)<0)]<-F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2[which((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref*F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_z_GWAS)<0)]*(-1)




pdf("F2DE_vs_PredictedF2DE_TimeImmobile_AlleVariants_bHRbLRgenes_NoColor_TopeGenesBlack.pdf", width=5, height=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2, col="grey", ylab="F2 DE: EPM Time Immobile T-stat", xlab="Predicted F2 DE: SMR T-stat")
Trendline1<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2)
abline(Trendline1, col="grey", lwd=3)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
Trendline2<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
abline(Trendline2, col="black", lwd=3)
dev.off()

summary.lm(Trendline1)

# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile ~ 
#        F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.1094 -0.6789  0.0163  0.6268  3.5298 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                 -0.066560   0.012745  -5.223 1.83e-07 ***
#   F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2  0.159440   0.007657  20.823  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.9661 on 5745 degrees of freedom
# (8606 observations deleted due to missingness)
# Multiple R-squared:  0.07018,	Adjusted R-squared:  0.07002 
# F-statistic: 433.6 on 1 and 5745 DF,  p-value: < 2.2e-16

sqrt(0.07018)
#[1] 0.2649151

summary.lm(Trendline2)

# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(CriteriaForPlot == 
#                                                                                       TRUE)] ~ F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2[which(CriteriaForPlot == 
#                                                                                                                                                                                    TRUE)])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.0968 -0.6305  0.0324  0.6636  3.3302 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                                                  0.06064    0.04448   1.363    0.173    
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2[which(CriteriaForPlot == TRUE)]  0.23516    0.02167  10.853   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.9866 on 490 degrees of freedom
# Multiple R-squared:  0.1938,	Adjusted R-squared:  0.1921 
# F-statistic: 117.8 on 1 and 490 DF,  p-value: < 2.2e-16

sqrt(0.1938)
#[1] 0.4402272

#With filtering:
pdf("F2DE_vs_PredictedF2DE_TimeImmobile_onlybHRbLRgenes_p05topgenes_TopeGenes.pdf", width=5, height=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)], ylab="F2 DE: EPM Time Immobile T-stat", xlab="Predicted F2 DE: SMR T-stat")
Trendline2<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
abline(Trendline2, col="black", lwd=3)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2>0)]
       ~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2>0)], col="red", pch=16)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2<0)]
       ~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2<0)], col="green4", pch=16)
dev.off()

#for labeling:
cbind.data.frame(Gene=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F2..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2>0)],
                 DETstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2>0)],
                 SMRTstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2>0)])

#     Gene  DETstat SMRTstat
# 1     Lsr 2.164527 5.079367
# 2 Plekhb1 3.473916 7.698405
# 3    <NA> 2.654667 4.377452
# 4 Tmem144 2.415786 4.569390

#the NA was the Rn60_1_2216.1 from the old v88 annotation

cbind.data.frame(Gene=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F2..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2<0)],
                 DETstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2<0)],
                 SMRTstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_time_immobile_s_T_SMR_directionalv2<0)])
# Gene   DETstat   SMRTstat
# 1    <NA> -2.327407 -11.114464
# 2   C2cd3 -2.199529  -5.969622
# 3 Ccdc167 -2.169617  -3.937098

#The NA is AABR07071904.1 from the v88 annotation


#################

#Predicted vs. Real F2 DE for Distance Traveled:

F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2<-abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR)
F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2[which((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref*F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_z_GWAS)<0)]<-F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2[which((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref*F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_z_GWAS)<0)]*(-1)



pdf("F2DE_vs_PredictedF2DE_EPMDistance__AlleVariants_bHRbLRgenes_NoColor_TopeGenesBlack.pdf", width=5, height=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2, col="grey", ylab="F2 DE: EPM Distance Traveled T-stat", xlab="Predicted F2 DE: SMR T-stat")
Trendline1<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2)
abline(Trendline1, col="grey", lwd=3)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
Trendline2<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
abline(Trendline2, col="black", lwd=3)
dev.off()

summary.lm(Trendline1)

# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled ~ 
#        F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.7293 -0.5759 -0.0053  0.5621  3.0284 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                    0.09541    0.01148    8.31   <2e-16 ***
#   F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2  0.14413    0.00596   24.18   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.8703 on 5745 degrees of freedom
# (8606 observations deleted due to missingness)
# Multiple R-squared:  0.0924,	Adjusted R-squared:  0.09224 
# F-statistic: 584.9 on 1 and 5745 DF,  p-value: < 2.2e-16

sqrt(0.0924)
#[1] 0.3039737

summary.lm(Trendline2)

# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(CriteriaForPlot == 
#                                                                                          TRUE)] ~ F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2[which(CriteriaForPlot == 
#                                                                                                                                                                                         TRUE)])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.76695 -0.61020  0.02357  0.63706  2.86522 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                                                   -0.005146   0.042634  -0.121    0.904    
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2[which(CriteriaForPlot == TRUE)]  0.249275   0.016302  15.291   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.9455 on 490 degrees of freedom
# Multiple R-squared:  0.323,	Adjusted R-squared:  0.3217 
# F-statistic: 233.8 on 1 and 490 DF,  p-value: < 2.2e-16

sqrt(0.323)
#[1] 0.5683309


#With filtering:
pdf("F2DE_vs_PredictedF2DE_EPMDistance_onlybHRbLRgenes_p05topgenes_TopeGenes.pdf", width=5, height=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)], ylab="F2 DE: EPM Distance Traveled T-stat", xlab="Predicted F2 DE: SMR T-stat")
Trendline2<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
abline(Trendline2, col="black", lwd=3)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2>0)]
       ~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2>0)], col="green4", pch=16)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2<0)]
       ~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2<0)], col="red", pch=16)
dev.off()

#for labeling:
cbind.data.frame(Gene=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F2..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2>0)],
                 DETstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2>0)],
                 SMRTstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2>0)])

#       Gene  DETstat  SMRTstat
# 1     <NA> 4.044018 21.676743
# 2    Mfge8 3.075677  4.315368
# 3    Fcrl2 2.352862  5.391312
# 4  Ndufaf6 2.194078  5.150067
# 5  Ankrd54 3.179116  4.868909
# 6     Maff 3.411116  4.206738
# 7     Sun2 2.566107  4.269776
# 8    Oard1 2.774671  4.580773
# 9   Loxhd1 2.716807  3.961857
# 10  Afg3l1 2.060752  4.479271

#NA is AABR07071904.1 in the v88 annotation


cbind.data.frame(Gene=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F2..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2<0)],
                 DETstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2<0)],
                 SMRTstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_distance_traveled_T_SMR_directionalv2<0)])

#     Gene   DETstat  SMRTstat
# 1   Fanci -1.972399 -4.538167
# 2   Wdr93 -3.298761 -4.328929
# 3  Unc45a -2.215841 -4.534881
# 4 Tmem144 -2.835746 -7.624189
# 5    Fzd6 -2.977735 -6.401529



#################################

#F2 predicted DE vs. Real F2 DE: EPM Percent time in the open arms

F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2<-abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR)
F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2[which((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref*F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_z_GWAS)<0)]<-F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2[which((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref*F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_z_GWAS)<0)]*(-1)



pdf("F2DE_vs_PredictedF2DE_EPMOpenArms__AlleVariants_bHRbLRgenes_NoColor_topeGenesBlack.pdf", width=5, height=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2, col="grey", ylab="F2 DE: EPM %Time Open Arms T-stat", xlab="Predicted F2 DE: SMR T-stat")
Trendline1<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2)
abline(Trendline1, col="grey", lwd=3)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
Trendline2<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
abline(Trendline2, col="black", lwd=3)
dev.off()

summary.lm(Trendline1)
# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms ~ 
#        F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.3547 -0.6782 -0.0044  0.6706  3.6765 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                                                                       0.062499   0.013019   4.801 1.62e-06
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2 0.149815   0.008247  18.165  < 2e-16
# 
# (Intercept)                                                                       ***
#   F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9869 on 5745 degrees of freedom
# (8606 observations deleted due to missingness)
# Multiple R-squared:  0.05432,	Adjusted R-squared:  0.05415 
# F-statistic:   330 on 1 and 5745 DF,  p-value: < 2.2e-16
sqrt(0.05432)
#[1] 0.2330665

summary.lm(Trendline2)
# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms[which(CriteriaForPlot == 
#                                                                                                TRUE)] ~ F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2[which(CriteriaForPlot == 
#                                                                                                                                                                                                   TRUE)])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.3639 -0.7537  0.0601  0.6743  3.5973 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                                                       -0.0006804  0.0444796  -0.015    0.988    
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2[which(CriteriaForPlot == TRUE)]  0.1924786  0.0245864   7.829 3.07e-14 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.9864 on 490 degrees of freedom
# Multiple R-squared:  0.1112,	Adjusted R-squared:  0.1094 
# F-statistic: 61.29 on 1 and 490 DF,  p-value: 3.07e-14

sqrt(0.1112)
#[1] 0.3334666


#With filtering:
pdf("F2DE_vs_PredictedF2DE_EPMOpenArms_onlybHRbLRgenes_p05topgenes_topeGenes.pdf", width=5, height=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)], ylab="F2 DE: EPM %Time Open Arms T-stat", xlab="Predicted F2 DE: SMR T-stat")
Trendline2<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
abline(Trendline2, col="black", lwd=3)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2>0)]
       ~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2>0)], col="green4", pch=16)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2<0)]
       ~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2<0)], col="red", pch=16)
dev.off()

#for labeling:
cbind.data.frame(Gene=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F2..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2>0)],
                 DETstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2>0)],
                 SMRTstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2>0)])

# Gene  DETstat SMRTstat
# 1 Maff 1.988845 4.567456
# 2 Sun2 3.262713 4.674993

cbind.data.frame(Gene=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F2..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2<0)],
                 DETstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2<0)],
                 SMRTstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$epm_percent_time_open_arm_T_SMR_directionalv2<0)])

# Gene   DETstat  SMRTstat
# 1 Hba-a3 -2.019930 -5.387071
# 2   Ghdc -2.098144 -5.923704




##########

F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2<-abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR)
F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2[which((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref*F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_z_GWAS)<0)]<-F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2[which((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref*F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_z_GWAS)<0)]*(-1)



pdf("F2DE_vs_PredictedF2DE_PCAIndex__AlleVariants_bHRbLRgenes_NoColor_topeGenesBlack.pdf", width=5, height=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2, col="grey", ylab="F2 DE: PCA Index T-stat", xlab="Predicted F2 DE: SMR T-stat")
Trendline1<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2)
abline(Trendline1, col="grey", lwd=3)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
Trendline2<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
abline(Trendline2, col="black", lwd=3)
dev.off()

summary.lm(Trendline1)

# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index ~ 
#        F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.4377 -0.5858  0.0326  0.6049  3.4135 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                  0.006245   0.011925   0.524    0.601    
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2 0.164762   0.007486  22.010   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.904 on 5745 degrees of freedom
# (8606 observations deleted due to missingness)
# Multiple R-squared:  0.07777,	Adjusted R-squared:  0.07761 
# F-statistic: 484.4 on 1 and 5745 DF,  p-value: < 2.2e-16

sqrt(0.07777)
#[1] 0.2788727

summary.lm(Trendline2)

# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index[which(CriteriaForPlot == 
#                                                                               TRUE)] ~ F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2[which(CriteriaForPlot == 
#                                                                                                                                                                             TRUE)])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.90925 -0.59302  0.05953  0.58600  2.22816 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                                                   0.02209    0.03940   0.561    0.575    
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2[which(CriteriaForPlot == TRUE)]  0.27083    0.02159  12.545   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.8735 on 490 degrees of freedom
# Multiple R-squared:  0.2431,	Adjusted R-squared:  0.2416 
# F-statistic: 157.4 on 1 and 490 DF,  p-value: < 2.2e-16

sqrt(0.2431)
#[1] 0.4930517


#With filtering:
pdf("F2DE_vs_PredictedF2DE_PCAIndex_onlybHRbLRgenes_p05topgenes_topeGenes.pdf", width=5, height=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)], ylab="F2 DE: PCA Index T-stat", xlab="Predicted F2 DE: SMR T-stat")
Trendline2<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
abline(Trendline2, col="black", lwd=3)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_PCA_Index<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2>0)]
       ~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_PCA_Index<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2>0)], col="green4", pch=16)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_PCA_Index<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2<0)]
       ~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_PCA_Index<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2<0)], col="red", pch=16)
dev.off()

#for labeling:
cbind.data.frame(Gene=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F2..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_PCA_Index<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2>0)],
                 DETstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_PCA_Index<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2>0)],
                 SMRTstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_PCA_Index<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2>0)])

#   Gene  DETstat SMRTstat
# 1 Mcee 2.314172 9.266280
# 2 Mcee 2.314172 7.461211

cbind.data.frame(Gene=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F2..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_PCA_Index<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2<0)],
                 DETstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_PCA_Index<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2<0)],
                 SMRTstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_PCA_Index<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_PCA_Index<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pca_index_days_6and7_T_SMR_directionalv2<0)])

# Gene   DETstat  SMRTstat
# 1  Fanci -2.172143 -4.006590
# 2 Pex11a -2.192097 -4.131513
# 3   Fzd6 -2.965447 -7.592240
# 4   Pdxp -2.031475 -8.517536
# 5   Cmc1 -2.037153 -4.167888
# 6   <NA> -2.091921 -5.230065

# NA is AABR07022144.1 in the v88 annotation


########################

#I should probably make some plots and run some stats for comparisons with the juvenile F2 SMR results too.

F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2<-abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR)
F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2[which((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref*F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_z_GWAS)<0)]<-F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2[which((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref*F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_z_GWAS)<0)]*(-1)

F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2<-abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR)
F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2[which((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref*F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_z_GWAS)<0)]<-F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2[which((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref*F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_z_GWAS)<0)]*(-1)

F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2<-abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR)
F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2[which((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref*F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_z_GWAS)<0)]<-F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2[which((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref*F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_z_GWAS)<0)]*(-1)



pdf("F2DE_vs_PredictedF2DE_EPMDistance_vs_OFDistance_AlleVariants_bHRbLRgenes_NoColor_topeGenesBlack.pdf", width=5, height=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2, col="grey", ylab="F2 DE: EPM Distance Traveled T-stat", xlab="Predicted F2 DE: OF Distance Traveled SMR T-stat")
Trendline1<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2)
abline(Trendline1, col="grey", lwd=3)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
Trendline2<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
abline(Trendline2, col="black", lwd=3)
dev.off()

summary.lm(Trendline1)

# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled ~ 
#        F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.1426 -0.5921 -0.0210  0.5785  3.2674 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                  0.096367   0.011961   8.057 9.47e-16 ***
#   F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2 0.050361   0.005359   9.398  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.9066 on 5744 degrees of freedom
# (8607 observations deleted due to missingness)
# Multiple R-squared:  0.01514,	Adjusted R-squared:  0.01497 
# F-statistic: 88.32 on 1 and 5744 DF,  p-value: < 2.2e-16

sqrt(0.01514)
#[1] 0.1230447

summary.lm(Trendline2)

# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(CriteriaForPlot == 
#                                                                                          TRUE)] ~ F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2[which(CriteriaForPlot == 
#                                                                                                                                                                                        TRUE)])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.0942 -0.7138 -0.0325  0.6698  3.2126 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                                                  -0.000612   0.048587  -0.013     0.99    
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2[which(CriteriaForPlot == TRUE)]  0.115520   0.014003   8.249 1.48e-15 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 1.077 on 490 degrees of freedom
# Multiple R-squared:  0.1219,	Adjusted R-squared:  0.1202 
# F-statistic: 68.05 on 1 and 490 DF,  p-value: 1.478e-15

sqrt(0.1219)
#[1] 0.3491418

#With filtering:
pdf("F2DE_vs_PredictedF2DE_EPMDistance_vs_OFDistance_onlybHRbLRgenes_p05topgenes_topeGenes.pdf", width=5, height=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)], ylab="F2 DE: EPM Distance Traveled T-stat", xlab="Predicted F2 DE: OF Distance Traveled SMR T-stat")
Trendline2<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
abline(Trendline2, col="black", lwd=3)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2>0)]
       ~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2>0)], col="green4", pch=16)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2<0)]
       ~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2<0)], col="red", pch=16)
dev.off()

#for labeling:
cbind.data.frame(Gene=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F2..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2>0)],
                 DETstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2>0)],
                 SMRTstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2>0)])

# Gene  DETstat SMRTstat
# 1  <NA> 4.044018 14.91432
# 2 Mfge8 3.075677 15.56612

#NA is AABR07071904.1 in the v88 annotation


cbind.data.frame(Gene=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F2..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2<0)],
                 DETstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2<0)],
                 SMRTstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_DistanceTraveled<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_distance_traveled_T_SMR_directionalv2<0)])

# Gene   DETstat   SMRTstat
# 1  Fanci -1.972399 -16.178725
# 2  Wdr93 -3.298761 -15.744017
# 3 Unc45a -2.215841 -11.011351
# 4   <NA> -2.460040  -9.277316
# 5   <NA> -2.681771  -5.685532

#The first NA is Rn60_1_2216.1, the second NA is AABR07043510.1 in the v88 annotation


#####


pdf("F2DE_vs_PredictedF2DE_EPMTimeImmobile_vs_OFTimeImmobile_AlleVariants_bHRbLRgenes_NoColor_topeGenesBlack.pdf", width=5, height=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2, col="grey", ylab="F2 DE: EPM Time Immobile T-stat", xlab="Predicted F2 DE: OF Time Immobile SMR T-stat")
Trendline1<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2)
abline(Trendline1, col="grey", lwd=3)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
Trendline2<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
abline(Trendline2, col="black", lwd=3)
dev.off()

summary.lm(Trendline1)

# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile ~ 
#        F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.8363 -0.7065  0.0090  0.6617  3.6288 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                -0.06988    0.01315  -5.315 1.11e-07 ***
#   F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2  0.05505    0.00709   7.766 9.56e-15 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.9967 on 5744 degrees of freedom
# (8607 observations deleted due to missingness)
# Multiple R-squared:  0.01039,	Adjusted R-squared:  0.01022 
# F-statistic:  60.3 on 1 and 5744 DF,  p-value: 9.56e-15

sqrt(0.01039)
#[1] 0.1019313

summary.lm(Trendline2)

# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(CriteriaForPlot == 
#                                                                                       TRUE)] ~ F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2[which(CriteriaForPlot == 
#                                                                                                                                                                                   TRUE)])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.6478 -0.7265  0.0120  0.7091  3.0509 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                                                 0.06123    0.04783   1.280    0.201    
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2[which(CriteriaForPlot == TRUE)]  0.11473    0.01920   5.975 4.44e-09 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 1.061 on 490 degrees of freedom
# Multiple R-squared:  0.0679,	Adjusted R-squared:  0.066 
# F-statistic:  35.7 on 1 and 490 DF,  p-value: 4.441e-09

sqrt(0.0679)
#[1] 0.2605763


#With filtering:
pdf("F2DE_vs_PredictedF2DE_EPMTimeImmobile_vs_OFTimeImmobile_onlybHRbLRgenes_p05topgenes_topeGenes.pdf", width=5, height=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)], ylab="F2 DE: EPM Time Immobile T-stat", xlab="Predicted F2 DE: OF Time Immobile SMR T-stat")
Trendline2<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
abline(Trendline2, col="black", lwd=3)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2>0)]
       ~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2>0)], col="red", pch=16)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2<0)]
       ~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2<0)], col="green4", pch=16)
dev.off()

#for labeling:
cbind.data.frame(Gene=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F2..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2>0)],
                 DETstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2>0)],
                 SMRTstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2>0)])

#     Gene  DETstat  SMRTstat
# 1   Wdr93 3.615647  7.163010
# 2 Plekhb1 3.473916 12.454720
# 3    Fzd6 2.890888  4.678009

cbind.data.frame(Gene=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$GENENAME_F2..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2<0)],
                 DETstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2<0)],
                 SMRTstat=F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Time_Immobile<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_time_immobile_s_T_SMR_directionalv2<0)])
# Gene   DETstat   SMRTstat
# 1    <NA> -2.327407 -11.358801
# 2   Mfge8 -2.136213  -7.125957
# 3   C2cd3 -2.199529 -12.326647
# 4 Ankrd54 -2.181199  -3.880168

#The NA is AABR07071904.1 in the v88 annotation



############



pdf("F2DE_vs_PredictedF2DE_EPMOpenArms_vs_OFTimeCenter_AlleVariants_bHRbLRgenes_NoColor_topeGenesBlack.pdf", width=5, height=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2, col="grey", ylab="F2 DE: EPM %Time Open Arms T-stat", xlab="Predicted F2 DE: OF %Time in Center SMR T-stat")
Trendline1<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2)
abline(Trendline1, col="grey", lwd=3)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
Trendline2<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
abline(Trendline2, col="black", lwd=3)
dev.off()

summary.lm(Trendline1)
# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms ~ 
#        F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.9268 -0.7127 -0.0068  0.6894  4.0079 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                       0.063697   0.013390   4.757 2.01e-06 ***
#   F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2 0.002281   0.007893   0.289    0.773    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 1.015 on 5744 degrees of freedom
# (8607 observations deleted due to missingness)
# Multiple R-squared:  1.453e-05,	Adjusted R-squared:  -0.0001596 
# F-statistic: 0.08349 on 1 and 5744 DF,  p-value: 0.7726


sqrt(1.453e-05)
#[1] 0.003811824
#basically no relationship

summary.lm(Trendline2)

# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms[which(CriteriaForPlot ==
#                                                                                                TRUE)] ~ F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2[which(CriteriaForPlot ==
#                                                                                                                                                                                                   TRUE)])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max
# -3.6678 -0.7829  0.0764  0.6797  3.5849
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                                                                                                       0.004634   0.047042   0.099   0.9216
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2[which(CriteriaForPlot == TRUE)] 0.036854   0.021917   1.682   0.0933 .
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 1.043 on 490 degrees of freedom
# Multiple R-squared:  0.005737,	Adjusted R-squared:  0.003708
# F-statistic: 2.828 on 1 and 490 DF,  p-value: 0.0933


sqrt(0.005737)
#[1] 0.07574299
#Basically no relationship

#With filtering:
pdf("F2DE_vs_PredictedF2DE_EPMOpenArms_vs_OFTimeCenter_onlybHRbLRgenes_p05topgenes_topeGenes.pdf", width=5, height=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)], ylab="F2 DE: EPM %Time Open Arms T-stat", xlab="Predicted F2 DE: OF %Time in Center SMR T-stat")
Trendline2<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms[which(CriteriaForPlot==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2[which(CriteriaForPlot==TRUE)])
abline(Trendline2, col="black", lwd=3)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2>0)]
       ~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2>0)], col="green4", pch=16)
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2<0)]
       ~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_p_SMR<0.05 & CriteriaForPlot==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_EPM_Percent_Time_Open_Arms<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$of_percent_time_in_center_T_SMR_directionalv2<0)], col="red", pch=16)
dev.off()
#No genes to label


########################

#Daniel was asking how much the correlation between SMR Tstat and F2 DE Tstat is driven by DE Genes being more likely to be strong eQTLs?

#...which is basically like asking about the strength of the QTL prediction on its own vs. the eQTL

#We might need to switch to non-directional stats for this to be tackled easily statistically
#First, the simple non-directional version, no filtering:
summary.lm(lm(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore)~abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2)))

# Call:
#   lm(formula = abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore) ~ 
#        abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2))
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.7210 -0.4723 -0.1326  0.3229  3.6579 
# 
# Coefficients:
#   Estimate Std. Error
# (Intercept)                                                                   0.707777   0.009798
# abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2) 0.049753   0.003152
# t value Pr(>|t|)    
# (Intercept)                                                                     72.24   <2e-16 ***
#   abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2)   15.79   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6223 on 5745 degrees of freedom
# (8606 observations deleted due to missingness)
# Multiple R-squared:  0.04157,	Adjusted R-squared:  0.0414 
# F-statistic: 249.2 on 1 and 5745 DF,  p-value: < 2.2e-16

plot(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore)~abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2))
#Without directionality, that is certainly uglier lol

#The simple non-directional version with filtering for bHR/bLR genes:
summary.lm(lm(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])~abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])))

# Call:
#   lm(formula = abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 
#                                                                             0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis. == 
#                                                                             1]) ~ abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 
#                                                                                                                                                                  0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis. == 
#                                                                                                                                                                  1]))
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.7933 -0.5727 -0.1469  0.4150  3.3946 
# 
# Coefficients:
#   Estimate
# (Intercept)                                                                                                                                                                                                           0.795242
# abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis. == 1]) 0.073459
# Std. Error
# (Intercept)                                                                                                                                                                                                             0.041917
# abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis. == 1])   0.008151
# t value
# (Intercept)                                                                                                                                                                                                            18.972
# abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis. == 1])   9.012
# Pr(>|t|)
# (Intercept)                                                                                                                                                                                                             <2e-16
# abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis. == 1])   <2e-16
# 
# (Intercept)                                                                                                                                                                                                           ***
#   abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis. == 1]) ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7883 on 518 degrees of freedom
# (391 observations deleted due to missingness)
# Multiple R-squared:  0.1355,	Adjusted R-squared:  0.1339 
# F-statistic: 81.22 on 1 and 518 DF,  p-value: < 2.2e-16

plot(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])~abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_T_SMR_directionalv2[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1]))


#Non-directional relationship with eQTLs (using Z-score), no filtering:
summary.lm(lm(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore)~abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_eQTL)))

# Call:
#   lm(formula = abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore) ~ 
#        abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_eQTL))
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.1941 -0.4849 -0.1352  0.3281  3.6739 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                                                      0.636015   0.017271   36.83   <2e-16
# abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_eQTL) 0.023162   0.002245   10.32   <2e-16
# 
# (Intercept)                                                      ***
#   abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_eQTL) ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6298 on 5745 degrees of freedom
# (8606 observations deleted due to missingness)
# Multiple R-squared:  0.01819,	Adjusted R-squared:  0.01802 
# F-statistic: 106.4 on 1 and 5745 DF,  p-value: < 2.2e-16
sqrt(0.01819)
#0.1348703

##Non-directional relationship with eQTLs (using Z-score), with filtering:
summary.lm(lm(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])~abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_eQTL[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])))

# Residual standard error: 0.8439 on 518 degrees of freedom
# (391 observations deleted due to missingness)
# Multiple R-squared:  0.009389,	Adjusted R-squared:  0.007477 
# F-statistic:  4.91 on 1 and 518 DF,  p-value: 0.02714
sqrt(0.009389)
#[1] 0.09689685

plot(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])~abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_eQTL[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1]))
#no high leverage points - so lack of relationship not due to high leverage

#Non-directional relationship with eQTLs (using aFC), no filtering:
summary.lm(lm(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore)~abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref)))

# Call:
#   lm(formula = abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore) ~ 
#        abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref))
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.6492 -0.4810 -0.1332  0.3251  3.3175 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                               0.739603   0.009833   75.22   <2e-16 ***
#   abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref) 0.161762   0.016148   10.02   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6301 on 5745 degrees of freedom
# (8606 observations deleted due to missingness)
# Multiple R-squared:  0.01717,	Adjusted R-squared:  0.017 
# F-statistic: 100.4 on 1 and 5745 DF,  p-value: < 2.2e-16

##Non-directional relationship with eQTLs (using aFC), with filtering:
summary.lm(lm(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])~abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])))

# Residual standard error: 0.8451 on 518 degrees of freedom
# (391 observations deleted due to missingness)
# Multiple R-squared:  0.006589,	Adjusted R-squared:  0.004671 
# F-statistic: 3.436 on 1 and 518 DF,  p-value: 0.06437

sqrt(0.006589)
#[1] 0.08117266

plot(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])~abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1]))
#Lots of high leverage points with this version of the analysis, but no visible relationship whatsoever.


#Non-directional relationship with QTLs (using Z-score), no filtering:
summary.lm(lm(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore)~abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_GWAS)))

# Call:
#   lm(formula = abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore) ~
#        abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_GWAS))
# 
# Residuals:
#   Min      1Q  Median      3Q     Max
# -1.3552 -0.4732 -0.1351  0.3255  3.6379
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                                                      0.660054   0.012977   50.86   <2e-16
# abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_GWAS) 0.120549   0.009129   13.21   <2e-16
# 
# (Intercept)                                                      ***
#   abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_GWAS) ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.6262 on 5745 degrees of freedom
# (8606 observations deleted due to missingness)
# Multiple R-squared:  0.02946,	Adjusted R-squared:  0.02929 
# F-statistic: 174.4 on 1 and 5745 DF,  p-value: < 2.2e-16

##Non-directional relationship with QTLs (using Z-score), with filtering:
summary.lm(lm(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])~abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_GWAS[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])))

# Residual standard error: 0.8026 on 518 degrees of freedom
# (391 observations deleted due to missingness)
# Multiple R-squared:  0.104,	Adjusted R-squared:  0.1022 
# F-statistic: 60.11 on 1 and 518 DF,  p-value: 4.795e-14

sqrt(0.104)
#[1] 0.3224903

#How much is that just driven by AABR...
plot(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])~abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_GWAS[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1]))
#nope, there's actually much fewer high leverage points.

#Combining QTLs and eQTLs w/o official SMR...

summary.lm(lm(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])~abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_GWAS[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])+abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_eQTL[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])))
# Pr(>|t|)
# (Intercept)                                                                                                                                                                                              2.24e-11
# abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_GWAS[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis. == 1]) 1.34e-13
# abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_eQTL[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis. == 1])   0.0853
# 
# (Intercept)                                                                                                                                                                                              ***
#   abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_GWAS[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis. == 1]) ***
#   abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_eQTL[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis. == 1]) .  
# Residual standard error: 0.8011 on 517 degrees of freedom
# (391 observations deleted due to missingness)
# Multiple R-squared:  0.1091,	Adjusted R-squared:  0.1056 
# F-statistic: 31.65 on 2 and 517 DF,  p-value: 1.076e-13

summary.lm(lm(abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Tstat_LocoScore[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])~abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_GWAS[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])+abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 &F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1])))

# abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_GWAS[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis. == 1]) 3.22e-14
# abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis. == 1])           0.038
# 
# (Intercept)                                                                                                                                                                                              ***
#   abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$total_loco_score_z_GWAS[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis. == 1]) ***
#   abs(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis. == 1])        *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8 on 517 degrees of freedom
# (391 observations deleted due to missingness)
# Multiple R-squared:  0.1114,	Adjusted R-squared:  0.108 
# F-statistic: 32.41 on 2 and 517 DF,  p-value: 5.492e-14


#So the SMR calculation does have added value
