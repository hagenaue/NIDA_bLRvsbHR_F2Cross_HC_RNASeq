########################

#Let's grab some F2 genotypes for QTL and eQTL plots:

setwd("/Users/hagenaue/Downloads/genotype_data/F2_genotypes_VCF")
#Based on the file size, we should be able to read these files in individually
list.files()

ExtractGT<-function(filename, VariantsOfInterest){
  vcf <- read.vcfR(filename, verbose = FALSE )
  gt_V2 <- extract.gt(vcf, element = c('GT'), as.numeric = FALSE)
  temp<-strsplit(row.names(gt_V2), "_")
  temp2<-do.call(rbind.data.frame, lapply(temp, as.character))
  rm(temp)
  colnames(temp2)[1]<-"CHROM"
  colnames(temp2)[2]<-"POS"
  temp2$POS<-as.numeric(temp2$POS)
  temp2$variant_id<-paste(temp2$CHROM, temp2$POS, sep=":")
  temp3<-cbind.data.frame(vcf@fix[,4], vcf@fix[,5])
  colnames(temp3)<-c("VCF_REF", "VCF_ALT")
  temp3$VCF_variant_id<-paste(vcf@fix[,1], vcf@fix[,2], sep=":")
  rm(vcf)
  gt_V2<-cbind.data.frame(temp2, temp3, gt_V2)
  rm(temp2, temp3)
  return(gt_V2[which(gt_V2$variant_id%in%VariantsOfInterest),])
  rm(gt_V2)
}

#I double-checked how things are working using Chr 19, since it is small:
Chr19_gt_F2<-ExtractGT("chr19_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr19_gt_F2)
#reading in that file took ~2.5 minutes
#subjects seem to be indicated with #s 87-852
#that doesn't match Genotyping id or RNA-Seq ID
#It does match rfid in the raw_phenotypes_N539.csv file in the Chitre et al. release
#That file has behavioral data.
#But the gene expression data ("ResidualsMatrix_wMean.csv" in the WGCNA folder) only has RNA-Seq sequencing id
#Looks like rfid is Rat" in "HRLR_F2_pca_matrix_RNAconc.txt"
#...which I do not appear to have on my computer. :(
#Also, the reference vs. alt allele isn't defined (so we don't know what 0 and 1 are)
#Although I suppose we could figure that out from the freq info?

#Let's double-check this
setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/F2_eQTLs")

F2_Phenotypes_forGT<-read.csv("raw_phenotypes_N539.csv", header=TRUE, stringsAsFactors = FALSE)

F2_Phenotypes_forGT$rfid
#this only includes #310-852
colnames(Chr19_gt_F2)[-c(1:6)])
#this includes #s 87-852
length(colnames(Chr19_gt_F2)[-c(1:6)])
#[1] 629
#And many more rats than we have in the F2 sample.
#Presumably some don't have adult behavioral data?
sum(is.na(F2_Phenotypes_forGT$total_loco_score)==FALSE)
#[1] 323
#That is a larger population than what we have RNA-Sequencing data on
#So the SMR analysis and DE results aren't completely circular...
sum(colnames(Chr19_gt_F2)%in%F2_Phenotypes_forGT$rfid[is.na(F2_Phenotypes_forGT$total_loco_score)==FALSE])
#[1] 323
#And vice versa:
sum(F2_Phenotypes_forGT$rfid[is.na(F2_Phenotypes_forGT$total_loco_score)==FALSE]%in%colnames(Chr19_gt_F2))
#[1] 323
#Great.

#Testing out whether I can visually recreate a QTL:
F2_Phenotypes_forGT_wAdultData<-F2_Phenotypes_forGT[is.na(F2_Phenotypes_forGT$total_loco_score)==FALSE,]
F2_Phenotypes_forGT_wAdultData$Rat<-as.character(F2_Phenotypes_forGT_wAdultData$rfid)
temp<-Chr19_gt_F2[,-c(1:6)]
temp2<-temp[,colnames(temp)%in%F2_Phenotypes_forGT_wAdultData$Rat]
sum(colnames(temp2)==F2_Phenotypes_forGT_wAdultData$Rat)
#[1] 323 - they are in the same order
cbind(colnames(temp2),F2_Phenotypes_forGT_wAdultData$Rat)
#Sanity check - looks good, they are in the same order

pdf("Boxplot_QTL_Afg3l1.pdf", height=5, width=4)
boxplot(F2_Phenotypes_forGT_wAdultData$total_loco_score~as.character(temp2[Chr19_gt_F2$variant_id=="chr19:56219141",]), xlab="chr19:56219141 (Afg3l1)", ylab="Total LocoScore", main="QTL", col=c("red", "gold", "green4"))
stripchart(F2_Phenotypes_forGT_wAdultData$total_loco_score~as.character(temp2[Chr19_gt_F2$variant_id=="chr19:56219141",]), vertical = TRUE,  method = "jitter", add = TRUE, cex=1, pch=1, col = 'black')
dev.off()

summary.lm(lm(F2_Phenotypes_forGT_wAdultData$total_loco_score~as.character(temp2[Chr19_gt_F2$variant_id=="chr19:56219141",])))
# Call:
#   lm(formula = F2_Phenotypes_forGT_wAdultData$total_loco_score ~ 
#        as.character(temp2[Chr19_gt_F2$variant_id == "chr19:56219141", 
#        ]))
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -534.13 -212.53  -38.04  170.47  925.98 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                            461.04      24.36  18.924  < 2e-16 ***
#   as.character(temp2[Chr19_gt_F2$variant_id == "chr19:56219141", ])0/1    73.98      33.33   2.220 0.027138 *  
#   as.character(temp2[Chr19_gt_F2$variant_id == "chr19:56219141", ])1/1   182.09      47.52   3.832 0.000153 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 276.7 on 320 degrees of freedom
# Multiple R-squared:  0.04591,	Adjusted R-squared:  0.03995 
# F-statistic:   7.7 on 2 and 320 DF,  p-value: 0.0005421

#That's without all of the normalization that was done to the behavioral data, so not bad.
#The 0/0 in this case clearly has more frequency than the 0/1
#Locomotor activity is higher in 1/1 than in 0/0 (positive slope)

Chr19_gt_F2$VCF_REF[Chr19_gt_F2$variant_id=="chr19:56219141"]
#[1] "C"
Chr19_gt_F2$VCF_ALT[Chr19_gt_F2$variant_id=="chr19:56219141"]
#[1] "T"

#Let's try another in Spg7:

pdf("Boxplot_QTL_Spg7.pdf", height=5, width=4)
boxplot(F2_Phenotypes_forGT_wAdultData$total_loco_score~as.character(temp2[Chr19_gt_F2$variant_id=="chr19:55989619",]), xlab="chr19:55989619 (Spg7)", ylab="Total LocoScore", main="QTL", col=c("red", "gold", "green4"))
stripchart(F2_Phenotypes_forGT_wAdultData$total_loco_score~as.character(temp2[Chr19_gt_F2$variant_id=="chr19:55989619",]), vertical = TRUE,  method = "jitter", add = TRUE, cex=1, pch=1, col = 'black')
dev.off()
#In this case, 0/0 clearly has *less* frequency than 1/1 (which means that it seems like the alternate allele)
#The slope is negative

Chr19_gt_F2$VCF_REF[Chr19_gt_F2$variant_id=="chr19:55989619"]
#[1] "C"
Chr19_gt_F2$VCF_ALT[Chr19_gt_F2$variant_id=="chr19:55989619"]
#[1] "T"


#Should we try another?

pdf("Boxplot_QTL_Vps9d1.pdf", height=5, width=4)
boxplot(F2_Phenotypes_forGT_wAdultData$total_loco_score~as.character(temp2[Chr19_gt_F2$variant_id=="chr19:55921978",]), xlab="chr19:55921978 (Vps9d1)", ylab="Total LocoScore", main="QTL", col=c("red", "gold", "green4"))
stripchart(F2_Phenotypes_forGT_wAdultData$total_loco_score~as.character(temp2[Chr19_gt_F2$variant_id=="chr19:55921978",]), vertical = TRUE,  method = "jitter", add = TRUE, cex=1, pch=1, col = 'black')
dev.off()
#In this case, 0/0 clearly has *less* frequency than 1/1 (which means that it seems like the alternate allele)
#The slope is negative

Chr19_gt_F2$VCF_REF[Chr19_gt_F2$variant_id=="chr19:55921978"]
#[1] "A"
Chr19_gt_F2$VCF_ALT[Chr19_gt_F2$variant_id=="chr19:55921978"]
#[1] "G"

#A few notes:
#1) A2 is the "reference" allele in the Chitre et al. results
#2) The Chitre et al. results define the reference allele in a manner so that Freq is almost always <0.5
#3) The reference allele in the Chitre results does not necessarily match the reference in the VCF
#4) Ref is the reference allele in the eQTL results.
#5) The reference allele in the eQTL results can have an af ranging from 0-1 and is not always the same as in Chitre et al.
#6) The reference allele in the eQTL results appears to be the same as the reference allele in the VCF (Daniel confirmed this)
#7) The definition of reference allele in the Z scored GWAS results is likely to be the same as in the eQTL analysis. (Daniel confirmed this)


Chr1_gt_F2<-ExtractGT("chr01_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr1_gt_F2)

Chr2_gt_F2<-ExtractGT("chr02_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr2_gt_F2)

Chr3_gt_F2<-ExtractGT("chr03_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr3_gt_V2)

Chr4_gt_F2<-ExtractGT("chr04_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr4_gt_F2)

Chr5_gt_F2<-ExtractGT("chr05_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr5_gt_F2)

Chr6_gt_F2<-ExtractGT("chr06_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr6_gt_F2)

Chr7_gt_F2<-ExtractGT("chr07_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr7_gt_F2)

Chr8_gt_F2<-ExtractGT("chr08_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr8_gt_F2)

Chr9_gt_F2<-ExtractGT("chr09_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr9_gt_F2)

Chr10_gt_F2<-ExtractGT("chr10_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr10_gt_F2)

Chr11_gt_F2<-ExtractGT("chr11_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr11_gt_F2)

Chr12_gt_F2<-ExtractGT("chr12_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr12_gt_F2)

Chr13_gt_F2<-ExtractGT("chr13_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr13_gt_F2)

Chr14_gt_F2<-ExtractGT("chr14_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr14_gt_F2)

Chr15_gt_F2<-ExtractGT("chr15_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr15_gt_F2)

Chr16_gt_F2<-ExtractGT("chr16_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr16_gt_F2)

Chr17_gt_F2<-ExtractGT("chr17_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr17_gt_F2)

Chr18_gt_F2<-ExtractGT("chr18_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr18_gt_F2)

Chr20_gt_F2<-ExtractGT("chr20_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr20_gt_F2)



setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/F2_eQTLs/VcfR_F0_Genetics")

AllAutosomes_gt_F2<-rbind.data.frame(Chr1_gt_F2, 
                                     Chr2_gt_F2, 
                                     Chr3_gt_F2, 
                                     Chr4_gt_F2, 
                                     Chr5_gt_F2,
                                     Chr6_gt_F2,
                                     Chr7_gt_F2,
                                     Chr8_gt_F2,
                                     Chr9_gt_F2,
                                     Chr10_gt_F2,
                                     Chr11_gt_F2,
                                     Chr12_gt_F2,
                                     Chr13_gt_F2,
                                     Chr14_gt_F2,
                                     Chr15_gt_F2,
                                     Chr16_gt_F2,
                                     Chr17_gt_F2,
                                     Chr18_gt_F2,
                                     Chr19_gt_F2,
                                     Chr20_gt_F2
)

str(AllAutosomes_gt_F2)
#'data.frame':	5653 obs. of  32 variables:

write.csv(AllAutosomes_gt_F2,"AllAutosomes_gt_F2.csv")

save.image("~/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/F2_eQTLs/VcfR_F0_Genetics/Workspace_VcfR_F0_Genetics.RData")

rm(Chr1_gt_F2, 
   Chr2_gt_F2, 
   Chr3_gt_F2, 
   Chr4_gt_F2, 
   Chr5_gt_F2,
   Chr6_gt_F2,
   Chr7_gt_F2,
   Chr8_gt_F2,
   Chr9_gt_F2,
   Chr10_gt_F2,
   Chr11_gt_F2,
   Chr12_gt_F2,
   Chr13_gt_F2,
   Chr14_gt_F2,
   Chr15_gt_F2,
   Chr16_gt_F2,
   Chr17_gt_F2,
   Chr18_gt_F2,
   Chr19_gt_F2,
   Chr20_gt_F2
)
