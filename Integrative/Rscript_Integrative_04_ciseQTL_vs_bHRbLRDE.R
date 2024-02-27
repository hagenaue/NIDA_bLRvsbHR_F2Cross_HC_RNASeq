
#Can the cis-eQTL (aFC) database predict bHR/bLR differential expression?

#First, double-checking how our results from the whole dataset play out in this data-frame (with some eGenes represented more than once & limited to genes expressed in F2)
#Meta-Analysis vs. F0
cor(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR,F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR, use="pairwise.complete.obs")
#[1] 0.3658721 (Pearson)
cor(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR,F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR, use="pairwise.complete.obs", method="spearman")
#[1] 0.3089781 (Spearman)

#Filtering by bHR/bLR segregation:
#Since hypothetically there should only be large aFC differentiating bHR/bLR if the genotypes segregate for that eVariant
#This would hypothetically represent the bHR/bLR DE that is driven by genetics

#Previously, I had found that variants with complete segregation of bHR/bLR into reference vs. heterozygote (0/0 vs. 0/1, chr1_166318198) have Gprimest of 0.27
#Whereas complete segregation (0/0 vs. 1/1) is Gprimest of 1
#I decided to use 0.27 as my cut-off
sum(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27, na.rm=TRUE)
#[1] 2500
sum(is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest)==FALSE)
#[1] 5747
2500/5747
#[1] 0.4350096
#So that's still 43% of the eVariants. Definitely a liberal criteria.

#First, double-checking how our results from the whole dataset play out in this subset:
#Meta-Analysis vs. F0
cor(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27],F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27], use="pairwise.complete.obs")
#[1] 0.5484089 (Pearson)
#(higher than without filtering)
cor(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27],F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27], use="pairwise.complete.obs", method="spearman")
#[1] 0.5605639
#(higher than without filtering)

#F0 vs. aFC
cor(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27],F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27], use="pairwise.complete.obs")
#[1] 0.7668082 (Pearson)
cor(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27],F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27], use="pairwise.complete.obs", method="spearman")
#[1] 0.6303265 (Spearman)
#sooooo pretty

pdf("Plot_F0_Log2FC_vs_aFC_Gst027.pdf", height=5, width=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27], ylab="F0: bLR vs. bHR Log2FC", xlab="Allelic Fold Change (aFC)", main="bLR/bHR segregated variants (F0 Gst>0.27)")
dev.off()
#So pretty...

#Publication quality figure:

library("scales")

pdf("Plot_F0_Log2FC_vs_aFC_Gst027_Nocolor.pdf", height=5, width=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27], ylab="F0: bLR vs. bHR Log2FC", xlab="Allelic Fold Change (aFC): Predicted bLR vs. bHR Log2FC", main="bLR/bHR Segregated Variants (F0 GST>0.27)", col=alpha("dimgrey", 0.5))
TrendLine<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27])
abline(TrendLine, col=1, lwd=3)
dev.off()


summary.lm(TrendLine)
# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 
#                                                                         0.27] ~ F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 
#                                                                                                                                        0.27])
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.77925 -0.11070 -0.00742  0.10350  3.15031 
# 
# Coefficients:
#   Estimate
# (Intercept)                                                                                               0.004515
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27] 0.528367
# Std. Error
# (Intercept)                                                                                                 0.004494
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27]   0.008849
# t value
# (Intercept)                                                                                                 1.005
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27]  59.708
# Pr(>|t|)
# (Intercept)                                                                                                  0.315
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27]   <2e-16
# 
# (Intercept)                                                                                                  
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27] ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2247 on 2498 degrees of freedom
# (8608 observations deleted due to missingness)
# Multiple R-squared:  0.588,	Adjusted R-squared:  0.5878 
# F-statistic:  3565 on 1 and 2498 DF,  p-value: < 2.2e-16

##############

#Congruence: 

#How many are there?

#bHR/bLR DE Genes (From the F2 candidate analysis)
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1)]))
#[1] 1045

#bHR/bLR DE eGenes:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$pval_beta)==FALSE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1)]))
#[1] 654

#Segregated bHR/bLR DE eGenes
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1)]))
#[1] 480


##################


#If I'm going to divide them up into bHR-like and bLR-like, I would need to refer to the meta-analysis results too:

bLRGenes<-((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_FDR_bLRvsbHR<0.10)
           |(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_FDR_bLRvsbHR<0.10)
           |(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Pvalue_bLRvsbHR<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_Pvalue_bLRvsbHR<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR>0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR>0))
sum(bLRGenes, na.rm=TRUE)
#[1] 689
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bLRGenes==TRUE)]))
#[1] 619
#This is the total number of bLR-genes - not just the ones that are also represented in the F2 dataset (608)

#Number of bLR genes that are eGenes:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bLRGenes==TRUE & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope)==FALSE)]))
#[1] 380

#Number of bLR gene cis-eQTLs:
length(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bLRGenes==TRUE & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope)==FALSE)])
#[1] 450

#Number of bLR genes that are bHR/bLR segregated eGenes:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bLRGenes & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope)==FALSE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)]))
#[1] 264

#Number of bLR gene cis-eQTLs that are bHR/bLR segregated:
length(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bLRGenes & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope)==FALSE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)])
#[1] 278

#Number of bLR genes that are eGenes with bHR/bLR segregated eVariants and bHR/bLR differential expression matching predictions:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR>0 & bLRGenes==TRUE)]))
#[1] 250

#Number of bLR gene cis-eQTLs with bHR/bLR segregated eVariants and bHR/bLR differential expression matching predictions:
length(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR>0 & bLRGenes==TRUE)])
#[1] 262


bHRGenes<-((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_FDR_bLRvsbHR<0.10)
           |(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_FDR_bLRvsbHR<0.10)
           |(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Pvalue_bLRvsbHR<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_Pvalue_bLRvsbHR<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR<0 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR<0))
sum(bHRGenes, na.rm=TRUE)
#[1] 501
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bHRGenes==TRUE)]))
#[1] 444
#(this is the number of bHR genes total - the number that also have F2 data is 437)

#Number of bHR genes that are eGenes:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bHRGenes==TRUE & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope)==FALSE)]))
#[1] 274

#Number of bHR gene cis-eQTLs:
length(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bHRGenes==TRUE & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope)==FALSE)])
#[1] 331

#Number of bHR genes that are bHR/bLR segregated eGenes:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bHRGenes==TRUE & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope)==FALSE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)]))
#[1] 216

#Number of bHR gene cis-eQTLs that are bHR/bLR segregated:
length(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bHRGenes==TRUE & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope)==FALSE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)])
#[1] 242

#Number of bHR genes that are eGenes with bHR/bLR segregated eVariants and differential expression matching predictions:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR<0 & bHRGenes==TRUE)]))
#[1] 206

#Number of bHR gene cis-eQTLs with bHR/bLR segregated eVariants and differential expression matching predictions:
length(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR<0 & bHRGenes==TRUE)])
#[1] 230

#Unique bHR/bLR genes that are eGenes:
380+274
#[1] 654
654/1450
#[1] 0.4510345
#A higher percentage than the rate of eGenes in the whole dataset (which is a little more than 1/3)

#Unique bHR/bLR genes that are eGenes with bHR/bLR segregated variants:
264+216
#[1] 480
480/654
#[1] 0.733945
#Definitely a higher percentage than the rate of bHR/bLR segregated variants in the whole cis-eQTL database (which is around 43%)

#Unique bHR/bLR genes that are eGenes with bHR/bLR segregated variants and bHR/bLR differential expression matching predictions:
250+206
#[1] 456 (wow - that's a pretty high proportion. I guess that makes sense.)
456/480
#[1] 0.95

#bHR/bLR gene cis-eQTLs with bHR/bLR segregated variants:
278+242
#[1] 520

#bHR/bLR gene cis-eQTLS with bHR/bLR segregated variants and bHR/bLR differential expression matching predictions:
262+230
#[1] 492 (wow - that's a pretty high proportion. I guess that makes sense)


bLRLikeBehaviorGenes<-((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_LocoScore<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Log2FC_LocoScore<0)
                       |(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Log2FC_EPM_Time_Immobile>0)
                       |(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Log2FC_EPM_DistanceTraveled<0)
                       |(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Log2FC_EPM_Percent_Time_Open_Arms<0)
                       |(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_PCA_Index <0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Log2FC_PCA_Index<0))
sum(bLRLikeBehaviorGenes, na.rm=TRUE)
#[1] 980
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bLRLikeBehaviorGenes==TRUE)]))
#[1] 924 


bHRLikeBehaviorGenes<-((F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_LocoScore<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Log2FC_LocoScore>0)
                       |(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Log2FC_EPM_Time_Immobile<0)
                       |(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Log2FC_EPM_DistanceTraveled>0)
                       |(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Log2FC_EPM_Percent_Time_Open_Arms>0)
                       |(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Pvalue_PCA_Index<0.05 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F2_Log2FC_PCA_Index>0))
sum(bHRLikeBehaviorGenes, na.rm=TRUE) 
#[1] 1140 
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bHRLikeBehaviorGenes==TRUE)]))
#[1] 1075 

sum(bLRGenes & bLRLikeBehaviorGenes, na.rm=TRUE)
#[1] 127
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bLRGenes==TRUE & bLRLikeBehaviorGenes==TRUE)]))
#[1] 111

sum(bHRGenes & bHRLikeBehaviorGenes, na.rm=TRUE)
#[1] 94
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bHRGenes==TRUE & bHRLikeBehaviorGenes==TRUE)]))
#[1] 81 


############

#Going a little bit further: Examining the overlap between our top candidate genes (bHR/bLR Genes that are DE with bHR/bLR-like F2 Behavior) and the cis-eQTL database

#bLR genes + bLR-like behavior genes:

#Number of bLR genes + bLR-like behavior genes that are eGenes:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bLRGenes==TRUE & bLRLikeBehaviorGenes==TRUE & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope)==FALSE)]))
#[1] 75
75/111
#[1] 0.6756757
#That's definitely a higher rate than the rest of the dataset

#Number of  bLR genes + bLR-like behavior genes cis-eQTLs:
length(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bLRGenes==TRUE & bLRLikeBehaviorGenes==TRUE & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope)==FALSE)])
#[1] 91

#Number of bLR genes + bLR-like behavior genes that are bHR/bLR segregated eGenes:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bLRGenes==TRUE & bLRLikeBehaviorGenes==TRUE & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope)==FALSE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)]))
#[1] 57
57/75
#[1] 0.76 - a much higher rate than amongst the eGenes as a whole (43%)

#Number of bLR genes + bLR-like behavior gene cis-eQTLs that are bHR/bLR segregated:
length(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bLRGenes==TRUE & bLRLikeBehaviorGenes==TRUE & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope)==FALSE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)])
#[1] 62

#Number of bLR genes + bLR-like behavior genes that are eGenes with bHR/bLR segregated eVariants and bHR/bLR differential expression matching predictions:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bLRGenes==TRUE & bLRLikeBehaviorGenes==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR>0)]))
#[1] 56
56/57
#[1] 0.9824561 - Nice

#Number of bLR gene cis-eQTLs with bHR/bLR segregated eVariants and bHR/bLR differential expression matching predictions:
length(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bLRGenes==TRUE & bLRLikeBehaviorGenes==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR>0)])
#[1] 61


#Time for bHR genes + bHR-like behavior genes:

#Number of bHR + bHR-like behavior genes that are eGenes:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bHRGenes==TRUE & bHRLikeBehaviorGenes==TRUE & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope)==FALSE)]))
#[1] 57
57/81
#[1] 0.7037037 definitely an enrichment over the full dataset

#Number of bHR + bHR-like behavior gene cis-eQTLs:
length(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bHRGenes==TRUE & bHRLikeBehaviorGenes==TRUE & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope)==FALSE)])
#[1] 70

#Number of bHR + bHR-like behavior genes that are bHR/bLR segregated eGenes:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bHRGenes==TRUE & bHRLikeBehaviorGenes==TRUE & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope)==FALSE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)]))
#[1] 47
47/57
#[1] 0.8245614 - again, much higher rate than in the general dataset.

#Number of bHR + bHR-like behavior gene cis-eQTLs that are bHR/bLR segregated:
length(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bHRGenes==TRUE & bHRLikeBehaviorGenes==TRUE & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$slope)==FALSE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)])
#[1] 51

#Number of bHR + bHR-like behavior genes that are eGenes with bHR/bLR segregated eVariants and differential expression matching predictions:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bHRGenes==TRUE & bHRLikeBehaviorGenes==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR<0)]))
#[1] 44
44/47
#[1] 0.9361702 - nice, the congruence is really strong.

#Number of bHR + bHR-like behavior gene cis-eQTLs with bHR/bLR segregated eVariants and differential expression matching predictions:
length(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bHRGenes==TRUE & bHRLikeBehaviorGenes==TRUE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR<0)])
#[1] 47

#Unique bHR/bLR+bHR/bLR behavior genes that are eGenes:
75+57
#[1] 132
132/(81+111)
#[1] 0.6875
#A higher percentage than the rate of eGenes in the whole dataset (which is a little more than 1/3)

#Unique bHR/bLR genes that are eGenes with bHR/bLR segregated variants:
57+47
#[1] 104
104/132
#[1] 0.7878788
#Definitely a higher percentage than the rate of bHR/bLR segregated variants in the whole cis-eQTL database (which is around 43%)

#Unique bHR/bLR genes that are eGenes with bHR/bLR segregated variants and bHR/bLR differential expression matching predictions:
56+44
#[1] 100
100/104
#[1] 0.9615385 - really pretty.

#For honing in within the SMR results, we'll need the number of eQTLs:

#bHR/bLR gene cis-eQTLs with bHR/bLR segregated variants:
62+51
#[1] 113

#bHR/bLR gene cis-eQTLS with bHR/bLR segregated variants and bHR/bLR differential expression matching predictions:
61+47
#[1] 108 

#Do we want to mark these off in the figure some way?


##################

pdf("Plot_F0_Log2FC_vs_aFC_Gst027_color_bHRbLRDEcandidateGenes.pdf", height=5, width=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)], ylab="F0: bLR vs. bHR Log2FC", xlab="Allelic Fold Change (aFC): Predicted bLR vs. bHR Log2FC", main="bLR/bHR Segregated Variants (F0 GST>0.27)", col=alpha("darkgrey", 0.5))
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR>0 & bLRGenes==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR>0 & bLRGenes==TRUE)], col=alpha("red", 0.25))
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR<0 & bHRGenes==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR<0 & bHRGenes==TRUE)], col=alpha("darkgreen", 0.25))
TrendLine<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)])
abline(TrendLine, col=1, lwd=3)
dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$F0_Log2FC_bLRvsbHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 
#                                                                               0.27)] ~ F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 
#                                                                                                                                                     0.27)])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.77925 -0.11070 -0.00742  0.10350  3.15031 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                                                      0.004515   0.004494   1.005    0.315    
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27)] 0.528367   0.008849  59.708   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.2247 on 2498 degrees of freedom
# Multiple R-squared:  0.588,	Adjusted R-squared:  0.5878 
# F-statistic:  3565 on 1 and 2498 DF,  p-value: < 2.2e-16


sqrt(0.588)
#0.7668116


###################


#Meta-Analysis vs. aFC
cor(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)], F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)], use="pairwise.complete.obs")
#[1] 0.5231619 (Pearson)

cor(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)], F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)], use="pairwise.complete.obs", method="spearman")
#[1] 0.6062834 (Spearman)

#So roughly equivalent to the relationship between F0 and aFC when rank-based

#No color:
pdf("Plot_MetaAnalysisD_vs_aFC_Gst027.pdf", height=5, width=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)], ylab="Meta-Analysis: bLR vs. bHR Estimated d", xlab="Allelic Fold Change (aFC)", main="bLR/bHR segregated variants (F0 Gst>0.27)")
dev.off()


pdf("Plot_MetaAnalysisD_vs_aFC_Gst027_color_bHRbLRDEcandidateGenes.pdf", height=5, width=5)
plot(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)], ylab="Meta-Analysis: bLR vs. bHR Estimated d", xlab="Allelic Fold Change (aFC): Predicted bLR vs. bHR Log2FC", main="bLR/bHR Segregated Variants (F0 GST>0.27)", col=alpha("darkgrey", 0.5))
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR>0 & bLRGenes==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR>0 & bLRGenes==TRUE)], col=alpha("red", 0.25))
points(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR<0 & bHRGenes==TRUE)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27 & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR<0 & bHRGenes==TRUE)], col=alpha("darkgreen", 0.25))
TrendLine<-lm(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)]~F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)])
abline(TrendLine, col=1, lwd=3)
dev.off()

summary.lm(TrendLine)
# Call:
#   lm(formula = F0_Meta_F2_DEResults_w_F2_eqtls_SMR$MetaAnalysis_estimatedD_bLRvsbHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 
#                                                                                             0.27)] ~ F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 
#                                                                                                                                                                   0.27)])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.0210 -0.5142 -0.0027  0.5213  7.2410 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                                                                       0.02678    0.01951   1.373     0.17    
# F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_bLR_bHR[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest > 0.27)]  1.16926    0.04145  28.211   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.8969 on 2112 degrees of freedom
# (386 observations deleted due to missingness)
# Multiple R-squared:  0.2737,	Adjusted R-squared:  0.2734 
# F-statistic: 795.9 on 1 and 2112 DF,  p-value: < 2.2e-16

sqrt(0.2737)
#[1] 0.5231635

