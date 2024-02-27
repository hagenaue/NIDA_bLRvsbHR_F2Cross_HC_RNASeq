
#Reading in the F0 genotyping results:
#Downloaded from the Chitre et al. supplemental data and results
#A longer version of this code that contains all the troubleshooting is found in "Code_F0_Genetics_Vcf.R"

setwd("/Users/hagenaue/Downloads/genotype_data/F0_genotypes")
list.files()

#These are the bHR Sample #s:
PhenotypeCluster<-colnames(gt)%in%c("5739-JL-0033", "5739-JL-0030", "5739-JL-0036", "5739-JL-0035","5739-JL-0032","5739-JL-0027","5739-JL-0029","5739-JL-0022","5739-JL-0028", "5739-JL-0023")
PhenotypeCluster

library(vcfR)

# filename<-"chr1_flt.vcf.gz"
# VariantsOfInterest<-F2_eqtls_SMR$variant_id

ExtractGTandDiff<-function(filename, VariantsOfInterest){
  vcf <- read.vcfR(filename, verbose = FALSE )
  gt_V2 <- extract.gt(vcf, element = c('GT'), as.numeric = FALSE)
  temp<-strsplit(row.names(gt_V2), "_")
  temp2<-do.call(rbind.data.frame, lapply(temp, as.character))
  rm(temp)
  colnames(temp2)[1]<-"CHROM"
  colnames(temp2)[2]<-"POS"
  temp2$POS<-as.numeric(temp2$POS)
  temp2$variant_id<-paste(temp2$CHROM, temp2$POS, sep=":")
  gt_V2<-cbind.data.frame(temp2, gt_V2)
  rm(temp2)
  myDiff <- genetic_diff(vcf, pops = as.factor(PhenotypeCluster), method = 'nei')
  rm(vcf)
  myDiff$variant_id<-paste(myDiff$CHROM, myDiff$POS, sep=":")
  gt_myDiff<-cbind.data.frame(gt_V2, myDiff[,c(-1,-2,-12)])
  rm(gt_V2, myDiff)
  return(gt_myDiff[which(gt_myDiff$variant_id%in%VariantsOfInterest),])
  rm(gt_myDiff)
}

Chr1_gt_V2<-ExtractGTandDiff("chr1_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr1_gt_V2)

Chr2_gt_V2<-ExtractGTandDiff("chr2_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr2_gt_V2)

Chr3_gt_V2<-ExtractGTandDiff("chr3_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr3_gt_V2)

Chr4_gt_V2<-ExtractGTandDiff("chr4_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr4_gt_V2)

Chr5_gt_V2<-ExtractGTandDiff("chr5_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr5_gt_V2)

Chr6_gt_V2<-ExtractGTandDiff("chr6_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr6_gt_V2)

Chr7_gt_V2<-ExtractGTandDiff("chr7_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr7_gt_V2)

Chr8_gt_V2<-ExtractGTandDiff("chr8_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr8_gt_V2)

Chr9_gt_V2<-ExtractGTandDiff("chr9_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr9_gt_V2)

Chr10_gt_V2<-ExtractGTandDiff("chr10_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr10_gt_V2)

Chr11_gt_V2<-ExtractGTandDiff("chr11_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr11_gt_V2)

Chr12_gt_V2<-ExtractGTandDiff("chr12_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr12_gt_V2)

Chr13_gt_V2<-ExtractGTandDiff("chr13_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr13_gt_V2)

Chr14_gt_V2<-ExtractGTandDiff("chr14_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr14_gt_V2)

Chr15_gt_V2<-ExtractGTandDiff("chr15_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr15_gt_V2)

Chr16_gt_V2<-ExtractGTandDiff("chr16_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr16_gt_V2)

Chr17_gt_V2<-ExtractGTandDiff("chr17_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr17_gt_V2)

Chr18_gt_V2<-ExtractGTandDiff("chr18_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr18_gt_V2)

Chr19_gt_V2<-ExtractGTandDiff("chr19_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr19_gt_V2)

Chr20_gt_V2<-ExtractGTandDiff("chr20_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(Chr20_gt_V2)

#I don't think there are any X chromosome results in the aFC data -  should I bother with it?

ChrX_gt_V2<-ExtractGTandDiff("chrX_flt.vcf.gz",F2_eqtls_SMR$variant_id) 
str(ChrX_gt_V2)
#yeah, there's nothing there.


setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/F2_eQTLs/VcfR_F0_Genetics")

AllAutosomes_gt_V2<-rbind.data.frame(Chr1_gt_V2, 
                                     Chr2_gt_V2, 
                                     Chr3_gt_V2, 
                                     Chr4_gt_V2, 
                                     Chr5_gt_V2,
                                     Chr6_gt_V2,
                                     Chr7_gt_V2,
                                     Chr8_gt_V2,
                                     Chr9_gt_V2,
                                     Chr10_gt_V2,
                                     Chr11_gt_V2,
                                     Chr12_gt_V2,
                                     Chr13_gt_V2,
                                     Chr14_gt_V2,
                                     Chr15_gt_V2,
                                     Chr16_gt_V2,
                                     Chr17_gt_V2,
                                     Chr18_gt_V2,
                                     Chr19_gt_V2,
                                     Chr20_gt_V2
)

str(AllAutosomes_gt_V2)
#'data.frame':	5653 obs. of  32 variables:
#'note: this is slightly different from the number of unique eVariants (5836) - some must not have survived QC in the F0 data (since it is a much smaller dataset)

hist(AllAutosomes_gt_V2$Gprimest)
summary(AllAutosomes_gt_V2$Gprimest)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.06304 0.20879 0.30889 0.54535 1.00000
#So around half the (HC expressed) eVariants are at least partially segregated
#and around 1/4 of the (Hc expressed) eVariants are very segregated.

write.csv(AllAutosomes_gt_V2,"AllAutosomes_gt_V2.csv")

save.image("~/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/F2_eQTLs/VcfR_F0_Genetics/Workspace_VcfR_F0_Genetics.RData")

rm(Chr1_gt_V2, 
   Chr2_gt_V2, 
   Chr3_gt_V2, 
   Chr4_gt_V2, 
   Chr5_gt_V2,
   Chr6_gt_V2,
   Chr7_gt_V2,
   Chr8_gt_V2,
   Chr9_gt_V2,
   Chr10_gt_V2,
   Chr11_gt_V2,
   Chr12_gt_V2,
   Chr13_gt_V2,
   Chr14_gt_V2,
   Chr15_gt_V2,
   Chr16_gt_V2,
   Chr17_gt_V2,
   Chr18_gt_V2,
   Chr19_gt_V2,
   Chr20_gt_V2,
   ChrX_gt_V2
)

#########################

#Joining the F0 genotype results to the eQTL and SMR results:

#double-checking the formatting:

str(AllAutosomes_gt_V2)
#'data.frame':	5653 obs. of  32 variables:

colnames(AllAutosomes_gt_V2)
#"variant_id"
head(AllAutosomes_gt_V2$variant_id)
#[1] "chr1:1482318"

str(F2_eqtls_SMR)
#'data.frame':	5937 obs. of  51 variables

head(F2_eqtls_SMR$variant_id)
#[1] "chr1:1482318"

#Yep, same formatting.

F2_eqtls_SMR_GT<-join(F2_eqtls_SMR, AllAutosomes_gt_V2, by="variant_id", type="left")

str(F2_eqtls_SMR_GT)
#'data.frame':	5937 obs. of  82 variables:
#Excellent.

sum(is.na(F2_eqtls_SMR_GT$log2_aFC_alt_ref)==FALSE)
#[1] 5937
#All eQTLs still present.

length(unique(F2_eqtls_SMR_GT$variant_id))
#[1] 5836  # of variants (instead of gene/variant combinations)

#Previously, I had found that variants with complete segregation of bHR/bLR into reference vs. heterozygote (0/0 vs. 0/1, chr1_166318198) have Gprimest of 0.27
#Whereas complete segregation (0/0 vs. 1/1) is Gprimest of 1
#I decided to use 0.27 as my cut-off

#Partially segregated:

sum(F2_eqtls_SMR_GT$Gprimest>0.27, na.rm=TRUE)
#[1] 2503
2503/length(F2_eqtls_SMR_GT$Gprimest)
#[1] 0.4215934

#unique variants:
length(unique(F2_eqtls_SMR_GT$variant_id[which(F2_eqtls_SMR_GT$Gprimest>0.27)]))
#[1] 2452
2452/5836
#[1] 0.4201508

#unique eGenes:
length(unique(F2_eqtls_SMR_GT$gene_id[which(F2_eqtls_SMR_GT$Gprimest>0.27)]))
#[1] 2398

#unique eQTLs:
length(unique(F2_eqtls_SMR_GT$Gene_Variant_ID[which(F2_eqtls_SMR_GT$Gprimest>0.27)]))
#[1] 2503

#Fully segregated:

sum(F2_eqtls_SMR_GT$Gprimest==1, na.rm=TRUE)
#[1] 23
23/length(F2_eqtls_SMR_GT$Gprimest)
#[1] 0.00387401

#unique eVariants:
length(unique(F2_eqtls_SMR_GT$variant_id[which(F2_eqtls_SMR_GT$Gprimest==1)]))
#[1] 23
23/5836
#[1] 0.003941056

sum(F2_eqtls_SMR_GT$Gprimest==0, na.rm=TRUE)
#[1] 210
210/length(F2_eqtls_SMR_GT$Gprimest)
#[1] 0.0353714

#unique eVariants:
length(unique(F2_eqtls_SMR_GT$variant_id[which(F2_eqtls_SMR_GT$Gprimest==0)]))
#[1] 205
205/5836
#[1] 0.0351268


#########################