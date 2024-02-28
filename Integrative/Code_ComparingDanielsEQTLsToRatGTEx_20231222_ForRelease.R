#Comparing Daniel Munro's hippocampal eQTLs from our F2 data to RatGTEx eQTLs
#Megan Hagenauer
#2023-10-25, updated 2023-12-22

###############

sessionInfo()
# R version 4.2.3 (2023-03-15)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Ventura 13.6.3
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] plyr_1.8.8                  VennDiagram_1.7.3           futile.logger_1.4.3         scales_1.2.1                SummarizedExperiment_1.28.0
# [6] Biobase_2.58.0              GenomicRanges_1.50.2        GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.2           
# [11] BiocGenerics_0.44.0         MatrixGenerics_1.10.0       matrixStats_1.0.0          
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.10            lattice_0.21-8         prettyunits_1.1.1      ps_1.7.4               digest_0.6.31          mime_0.12              R6_2.5.1              
# [8] futile.options_1.0.1   httr_1.4.6             zlibbioc_1.44.0        rlang_1.1.1            rstudioapi_0.14        data.table_1.14.8      miniUI_0.1.1.1        
# [15] callr_3.7.3            urlchecker_1.0.1       Matrix_1.5-4.1         devtools_2.4.5         stringr_1.5.0          htmlwidgets_1.6.2      RCurl_1.98-1.12       
# [22] bit_4.0.5              munsell_0.5.0          shiny_1.7.4            DelayedArray_0.24.0    compiler_4.2.3         httpuv_1.6.9           pkgbuild_1.4.2        
# [29] htmltools_0.5.5        GenomeInfoDbData_1.2.9 crayon_1.5.2           later_1.3.0            bitops_1.0-7           xtable_1.8-4           lifecycle_1.0.3       
# [36] formatR_1.14           magrittr_2.0.3         cli_3.6.1              stringi_1.7.12         cachem_1.0.8           farver_2.1.1           XVector_0.38.0        
# [43] fs_1.6.1               promises_1.2.0.1       remotes_2.4.2.1        ellipsis_0.3.2         vctrs_0.6.3            lambda.r_1.2.4         gemma.R_1.3.2         
# [50] tools_4.2.3            bit64_4.0.5            glue_1.6.2             purrr_1.0.1            processx_3.8.0         pkgload_1.3.2          fastmap_1.1.1         
# [57] colorspace_2.1-0       sessioninfo_1.2.2      memoise_2.0.1          profvis_0.3.8          usethis_2.2.2  

##################

#RatGTEx eQTLs:

#RatGTEx eQTLs downloaded from:
#https://ratgtex.org/data/eqtl/rn6.eqtls_indep.txt
#on Oct 25 2023

#here's the description of the RatGTEx project:
#https://ratgtex.org

#Note: 

#All of the brain eQTLs in this database come from heterogeneous stock rats
#There are many brain regions available
#These are all cis-eQTLs (all studies underpowered to detect trans eQTLs)

#The "Brain hemisphere" eQTLs come from the largest sample size (n=339)
###Unpublished:  Francesca Telese, UCSD, "Pilot: Creating the dataset for TWAS in HS rats"

#The basolateral amygdala, prefrontal cortex, and nucleus accumbens core come a separate sample (n=185-191)
###Unpublished: Suzanne Mitchell and Robert Hitzemann, OHSU, "U01DA046077: Identification of Genetic Features of Delay Discounting Using a Heterogeneous Stock Rat Model"

#The other regions have pretty dinky sample sizes (n=75-81) and come from a separate sample.
###Publication: Munro et al. 2022: "The regulatory landscape of multiple brain regions in outbred heterogeneous stock rats"

setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/F2_eQTLs")

RatGTEx<-read.delim("RatGTEX_rn6eqtls_indep_Download20231025.txt", sep="\t", header=TRUE)

str(RatGTEx)

# 'data.frame':	71343 obs. of  17 variables:
# $ tissue      : chr  "Adipose" "Adipose" "Adipose" "Adipose" ...
# $ gene_id     : chr  "ENSRNOG00000047964" "ENSRNOG00000050370" "ENSRNOG00000032365" "ENSRNOG00000032365" ...
# $ gene_name   : chr  "AABR07000089.1" "Vom2r6" "Vom2r5" "Vom2r5" ...
# $ num_var     : int  17 17 20 20 38 74 77 82 92 97 ...
# $ variant_id  : chr  "chr1:1521157" "chr1:1690027" "chr1:1757115" "chr1:1762774" ...
# $ chrom       : int  1 1 1 1 1 1 1 1 1 1 ...
# $ pos         : int  1521157 1690027 1757115 1762774 2015941 2361293 1771476 2359016 2222622 2336604 ...
# $ ref         : chr  "C" "G" "G" "G" ...
# $ alt         : chr  "G" "A" "A" "A" ...
# $ af          : num  0.595 0.56 0.661 0.944 0.678 ...
# $ tss_distance: int  821745 -934382 -943598 -949257 835010 -968964 -329290 -884170 605238 633909 ...
# $ pval_nominal: num  2.59e-06 1.17e-05 1.33e-21 3.65e-05 1.09e-67 ...
# $ slope       : num  -0.355 -0.276 -0.496 0.456 -0.87 ...
# $ slope_se    : num  0.0744 0.0622 0.0489 0.109 0.0406 ...
# $ pval_beta   : num  7.10e-05 2.27e-04 1.08e-18 9.59e-04 4.68e-62 ...
# $ rank        : int  1 1 1 2 1 1 1 1 1 1 ...
# $ log2_aFC    : num  -0.7624 -0.581 -0.5965 0.0712 -0.5668 ...

#I'm going to assume that their ref/alt won't necessarily match our own.
#Actually, Daniel says that they should match.

table(RatGTEx$tissue)
# Adipose     BLA   Brain     Eye      IL     LHb   Liver    NAcc   NAcc2     OFC      PL     PL2 
# 11809    6064   13293     567    3847    3281    9461    3366    5164    3798    3906    6787


##########################

#Reading in Daniel's database of conditionally-independent HC eQTLs (created using the F2 data)

F2_eQTLs<-read.delim("f2.eqtls_indep.txt", sep="\t", header=TRUE)

str(F2_eQTLs)
# 'data.frame':	5937 obs. of  18 variables:
# $ gene_id         : chr  "ENSRNOG00000014330" "ENSRNOG00000015239" "ENSRNOG00000016381" "ENSRNOG00000013160" ...
# $ gene_name       : chr  "Pcmt1" "Ginm1" "Ust" "Sash1" ...
# $ num_var         : int  1872 2074 3342 3630 4222 4297 4297 3180 3180 2213 ...
# $ variant_id      : chr  "chr1:1482318" "chr1:2828974" "chr1:3221846" "chr1:2774862" ...
# $ chrom           : int  1 1 1 1 1 1 1 1 1 1 ...
# $ pos             : int  1482318 2828974 3221846 2774862 4016846 4366305 3618606 4605341 5629945 6387063 ...
# $ ref             : chr  "G" "G" "C" "G" ...
# $ alt             : chr  "T" "A" "G" "A" ...
# $ af              : num  0.0266 0.4119 0.3053 0.5922 0.6393 ...
# $ tss_distance    : int  285300 -942993 -594099 71338 -253454 -354794 392905 47869 -976735 938105 ...
# $ pval_nominal    : num  1.00e-07 9.45e-12 2.53e-04 5.77e-11 2.69e-35 ...
# $ slope           : num  -0.613 -0.521 0.326 -0.502 0.814 ...
# $ slope_se        : num  0.1114 0.0726 0.0876 0.0731 0.055 ...
# $ pval_beta       : num  1.84e-05 2.95e-09 2.76e-02 2.12e-08 7.08e-31 ...
# $ rank            : int  1 1 1 1 1 1 2 1 2 1 ...
# $ log2_aFC_alt_ref: num  0.3959 -0.3171 0.0823 -0.1907 0.5983 ...
# $ log2_aFC_bLR_bHR: num  0.3959 0.3171 0.0823 -0.1907 -0.5983 ...
# $ upregulated_in  : chr  "bLR" "bLR" "bLR" "bHR" ...

#Outputting some basic characteristics for the paper:

pdf("Histogram_HCeQTLs_TssDistance.pdf", height=4, width=5)
hist(F2_eQTLs$tss_distance, breaks=100, col="grey")
dev.off()

table(table(F2_eQTLs$gene_id))
# 1    2    3    4 
# 4798  523   27    3

table(table(F2_eQTLs$variant_id))

# 1    2    3 
# 5741   89    6 

plot(slope~log2_aFC_alt_ref, data=F2_eQTLs)

length(unique(F2_eQTLs$variant_id))
#[1] 5836


#####################

#Examining overlap between our conditionally-independent HC eQTLs and Rat GTEx

sum(F2_eQTLs$variant_id%in%RatGTEx$variant_id)
#[1] 115
#hmmmm...
#Well, this may be a short-lived endeavor, LOL
#As far as I can tell, it doesn't seem to be a formatting issue
#Talking to Daniel, it seems like this lack of overlap is likely to be due the uncertainty when choosing eVariants introduced by LD blocks

#Maybe instead we could peek at whether the same genes seem to be regulated by eQTLs?

sum(unique(F2_eQTLs$gene_id)%in%unique(RatGTEx$gene_id))
#[1] 4786
#That's pretty decent overlap, but could just be due to random chance.

length(unique(RatGTEx$gene_id))
#[1] 16555

length(unique(RatGTEx$variant_id))
#[1] 61146

#For interpretation, it would be helpful to know how many of the RatGTEx eGenes aren't even detectably expressed in the rat HC (i.e. our full list, not our eGene list):

Genes_inF2Data<-read.csv("Genes_inF2ExpressionDataSet.csv", header=TRUE, stringsAsFactors = FALSE)
str(Genes_inF2Data)
# 'data.frame':	14056 obs. of  1 variable:
# $ gene_id: chr  "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...

sum(unique(RatGTEx$gene_id)%in%unique(Genes_inF2Data$gene_id))
#[1] 11365

#Overlap by region:

#Brain Hemisphere

sum(F2_eQTLs$gene_id%in%RatGTEx$gene_id[RatGTEx$tissue=="Brain"])
#[1] 4339

sum(unique(RatGTEx$gene_id[RatGTEx$tissue=="Brain"])%in%unique(Genes_inF2Data$gene_id))
#[1] 8222

sum(unique(F2_eQTLs$gene_id)%in%unique(RatGTEx$gene_id[RatGTEx$tissue=="Brain"]))
#[1] 3854

length(unique(F2_eQTLs$gene_id))
#[1] 5351

3854/5351
#[1] 0.7202392 
#72% of the eGenes that Daniel identified are already known to be brain eGenes.

length(unique(RatGTEx$gene_id[RatGTEx$tissue=="Brain"]))
#[1] 10295
#But a pretty large number of genes that are expressed in the brain have been identified as eGenes...
#So not sure how much that buys us.

3854/10295
#[1] 0.3743565

library("VennDiagram")   

pdf("Venn_RatGTEx_Brain_Vs_OurHC_eGenes.pdf", height=2, width=2)
grid.newpage()                                        
draw.pairwise.venn(5351,10295, 3854,col=c("grey", "grey"), fill=c("lightgrey", "lightgrey"))
dev.off()

#BLA

sum(unique(F2_eQTLs$gene_id)%in%unique(RatGTEx$gene_id[RatGTEx$tissue=="BLA"]))
#[1] 2495

2495/5351
#[1] 0.466268

length(unique(RatGTEx$gene_id[RatGTEx$tissue=="BLA"]))
#[1] 5436

sum(unique(RatGTEx$gene_id[RatGTEx$tissue=="BLA"])%in%unique(Genes_inF2Data$gene_id))
#[1] 4628

2495/5436
#[1] 0.4589772

pdf("Venn_RatGTEx_BLA_Vs_OurHC_eGenes.pdf", height=2, width=2)
grid.newpage()                                        
draw.pairwise.venn(5351,5436,2495,col=c("grey", "grey"), fill=c("lightgrey", "lightgrey"))
dev.off()

sum(unique(F2_eQTLs$gene_id)%in%unique(RatGTEx$gene_id[RatGTEx$tissue=="PL"]))
#[1] 1756

1756/5351
#[1] 0.328163

length(unique(RatGTEx$gene_id[RatGTEx$tissue=="PL"]))
#[1] 3755

sum(unique(RatGTEx$gene_id[RatGTEx$tissue=="PL"])%in%unique(Genes_inF2Data$gene_id))
#[1] 3180

1756/3755
#[1] 0.4676431

pdf("Venn_RatGTEx_PL_Vs_OurHC_eGenes.pdf", height=2, width=2)
grid.newpage()                                        
draw.pairwise.venn(5351,3755, 1756,col=c("grey", "grey"), fill=c("lightgrey", "lightgrey"))
dev.off()

sum(unique(F2_eQTLs$gene_id)%in%unique(RatGTEx$gene_id[RatGTEx$tissue=="PL2"]))
#[1] 2649

2649/5351
#[1] 0.4950477

length(unique(RatGTEx$gene_id[RatGTEx$tissue=="PL2"]))
#[1] 5985

sum(unique(RatGTEx$gene_id[RatGTEx$tissue=="PL2"])%in%unique(Genes_inF2Data$gene_id))
#[1] 5052

2649/5985
#[1] 0.4426065

pdf("Venn_RatGTEx_PL2_Vs_OurHC_eGenes.pdf", height=2, width=2)
grid.newpage()                                        
draw.pairwise.venn(5351,5985, 2649,col=c("grey", "grey"), fill=c("lightgrey", "lightgrey"))
dev.off()

sum(unique(F2_eQTLs$gene_id)%in%unique(RatGTEx$gene_id[RatGTEx$tissue=="OFC"]))
#[1] 1749

1749/5351
#[1] 0.3268548

length(unique(RatGTEx$gene_id[RatGTEx$tissue=="OFC"]))
#[1] 3638

sum(unique(RatGTEx$gene_id[RatGTEx$tissue=="OFC"])%in%unique(Genes_inF2Data$gene_id))
#[1] 3116

1749/3638
#[1] 0.4807587

pdf("Venn_RatGTEx_OFC_Vs_OurHC_eGenes.pdf", height=2, width=2)
grid.newpage()                                        
draw.pairwise.venn(5351,3638, 1749,col=c("grey", "grey"), fill=c("lightgrey", "lightgrey"))
dev.off()

sum(unique(F2_eQTLs$gene_id)%in%unique(RatGTEx$gene_id[RatGTEx$tissue=="IL"]))
#[1] 1739

1739/5351
#[1] 0.324986

length(unique(RatGTEx$gene_id[RatGTEx$tissue=="IL"]))
#[1] 3701

sum(unique(RatGTEx$gene_id[RatGTEx$tissue=="IL"])%in%unique(Genes_inF2Data$gene_id))
#[1] 3147

1739/3701
#[1] 0.469873

pdf("Venn_RatGTEx_OFC_Vs_OurHC_eGenes.pdf", height=2, width=2)
grid.newpage()                                        
draw.pairwise.venn(5351,3701, 1739,col=c("grey", "grey"), fill=c("lightgrey", "lightgrey"))
dev.off()


#non-cortical regions:

sum(unique(F2_eQTLs$gene_id)%in%unique(RatGTEx$gene_id[RatGTEx$tissue=="LHb"]))
#[1] 1512

1512/5351
#[1] 0.282564

length(unique(RatGTEx$gene_id[RatGTEx$tissue=="LHb"]))
#[1] 3162

sum(unique(RatGTEx$gene_id[RatGTEx$tissue=="LHb"])%in%unique(Genes_inF2Data$gene_id))
#[1] 2688

1512/3162
#[1] 0.4781784

pdf("Venn_RatGTEx_LHb_Vs_OurHC_eGenes.pdf", height=2, width=2)
grid.newpage()                                        
draw.pairwise.venn(5351,3162,1512,col=c("grey", "grey"), fill=c("lightgrey", "lightgrey"))
dev.off()

sum(unique(F2_eQTLs$gene_id)%in%unique(RatGTEx$gene_id[RatGTEx$tissue=="NAcc"]))
#[1] 1528

1528/5351
#[1] 0.2855541

length(unique(RatGTEx$gene_id[RatGTEx$tissue=="NAcc"]))
#[1] 3260

sum(unique(RatGTEx$gene_id[RatGTEx$tissue=="NAcc"])%in%unique(Genes_inF2Data$gene_id))
#[1] 2774

1528/3260
#[1] 0.4687117

pdf("Venn_RatGTEx_NAcc_Vs_OurHC_eGenes.pdf", height=2, width=2)
grid.newpage()                                        
draw.pairwise.venn(5351,3260,1528,col=c("grey", "grey"), fill=c("lightgrey", "lightgrey"))
dev.off()

sum(unique(F2_eQTLs$gene_id)%in%unique(RatGTEx$gene_id[RatGTEx$tissue=="NAcc2"]))
#[1] 2171

2171/5351
#[1] 0.4057186

length(unique(RatGTEx$gene_id[RatGTEx$tissue=="NAcc2"]))
#[1] 4717

sum(unique(RatGTEx$gene_id[RatGTEx$tissue=="NAcc2"])%in%unique(Genes_inF2Data$gene_id))
#[1] 3983

2171/4717
#[1] 0.4602502

pdf("Venn_RatGTEx_NAcc2_Vs_OurHC_eGenes.pdf", height=2, width=2)
grid.newpage()                                        
draw.pairwise.venn(5351,4717,2171,col=c("grey", "grey"), fill=c("lightgrey", "lightgrey"))
dev.off()

#Non brain tissue:

sum(unique(F2_eQTLs$gene_id)%in%unique(RatGTEx$gene_id[RatGTEx$tissue=="Adipose"]))
#[1] 3233

3233/5351
#[1] 0.6041861

length(unique(RatGTEx$gene_id[RatGTEx$tissue=="Adipose"]))
#[1] 9750

sum(unique(RatGTEx$gene_id[RatGTEx$tissue=="Adipose"])%in%unique(Genes_inF2Data$gene_id))
#[1] 7216

3233/9750
#[1] 0.3315897

pdf("Venn_RatGTEx_Adipose_Vs_OurHC_eGenes.pdf", height=2, width=2)
grid.newpage()                                        
draw.pairwise.venn(5351,9750,3233,col=c("grey", "grey"), fill=c("lightgrey", "lightgrey"))
dev.off()

sum(unique(F2_eQTLs$gene_id)%in%unique(RatGTEx$gene_id[RatGTEx$tissue=="Liver"]))
#[1] 2792

2792/5351
#[1] 0.6041861

length(unique(RatGTEx$gene_id[RatGTEx$tissue=="Liver"]))
#[1] 7818

sum(unique(RatGTEx$gene_id[RatGTEx$tissue=="Liver"])%in%unique(Genes_inF2Data$gene_id))
#[1] 6238

2792/7818
#[1] 0.3571246

pdf("Venn_RatGTEx_Liver_Vs_OurHC_eGenes.pdf", height=2, width=2)
grid.newpage()                                        
draw.pairwise.venn(5351,7818,2792,col=c("grey", "grey"), fill=c("lightgrey", "lightgrey"))
dev.off()

###############

#I wonder how frequently genes are showing up in RatGTEx?

#Hmm... I realized that is a complicated question to ask, because some genes will show up repeatedly due to being in LD with a lot of other genes/variants
#So it seems like the best way to ask the question would be to ask how many tissues each gene shows up in

RatGTEx_GeneFrequency<-data.frame(gene_id=rep("NA",length(unique(RatGTEx$gene_id))), gene_name=rep("NA",length(unique(RatGTEx$gene_id))), NumberOfTissues=rep(0,length(unique(RatGTEx$gene_id))), Tissues=rep("NA",length(unique(RatGTEx$gene_id))))
str(RatGTEx_GeneFrequency)

for(i in c(1:length(unique(RatGTEx$gene_id)))){
  print(i)
  TempGene<-names(table(unique(RatGTEx$gene_id)))[i]
  RatGTEx_GeneFrequency[i,1]<-TempGene
  RatGTEx_GeneFrequency[i,2]<-RatGTEx$gene_name[RatGTEx$gene_id==TempGene][1]
  RatGTEx_GeneFrequency[i,3]<-length(unique(RatGTEx$tissue[RatGTEx$gene_id==TempGene]))
  #RatGTEx_GeneFrequency[i,4]<-unique(RatGTEx$tissue[RatGTEx$gene_id==names(table(unique(RatGTEx$gene_id)))[i]])
}
  

write.csv(RatGTEx_GeneFrequency, "RatGTEx_GeneFrequency.csv")

table(RatGTEx_GeneFrequency$NumberOfTissues)
# 1    2    3    4    5    6    7    8    9   10   11   12 
# 5009 2921 2057 1442 1045  907  713  600  617  547  555  142 

str(F2_eQTLs)

library(plyr)
F2_eQTLs_vs_RatGTEx_GeneFrequency<-join(F2_eQTLs, RatGTEx_GeneFrequency, by="gene_id", type="left" )
write.csv(F2_eQTLs_vs_RatGTEx_GeneFrequency, "F2_eQTLs_vs_RatGTEx_GeneFrequency.csv")

pdf("Boxplot_AbsLog2AFC_vs_GTExNumberOfTissues.pdf", height=12, width=4)
boxplot(abs(F2_eQTLs_vs_RatGTEx_GeneFrequency$log2_aFC_alt_ref)~F2_eQTLs_vs_RatGTEx_GeneFrequency$NumberOfTissues)
dev.off()

pdf("Boxplot_AbsSlope_vs_GTExNumberOfTissues.pdf", height=12, width=4)
boxplot(abs(F2_eQTLs_vs_RatGTEx_GeneFrequency$slope)~F2_eQTLs_vs_RatGTEx_GeneFrequency$NumberOfTissues)
dev.off()
#That has a pretty clear increasing pattern after 7

pdf("Boxplot_Log10p_vs_GTExNumberOfTissues.pdf", height=12, width=4)
boxplot(log10(F2_eQTLs_vs_RatGTEx_GeneFrequency$pval_beta)~F2_eQTLs_vs_RatGTEx_GeneFrequency$NumberOfTissues)
dev.off()
#That has a pretty clear correlation too

F2_eQTLs_GeneIDs<-data.frame(gene_id=unique(F2_eQTLs$gene_id), HC_eQTL=rep(1, length(unique(F2_eQTLs$gene_id))))


F2_eQTLs_vs_RatGTEx_GeneFrequencyFull<-join(F2_eQTLs_GeneIDs, RatGTEx_GeneFrequency, by="gene_id", type="full" )
str(F2_eQTLs_vs_RatGTEx_GeneFrequencyFull)
#'data.frame':	17120 obs. of  6 variables

F2_eQTLs_vs_RatGTEx_GeneFrequencyFull$NumberOfTissues[is.na(F2_eQTLs_vs_RatGTEx_GeneFrequencyFull$NumberOfTissues)]<-0
F2_eQTLs_vs_RatGTEx_GeneFrequencyFull$HC_eQTL[is.na(F2_eQTLs_vs_RatGTEx_GeneFrequencyFull$HC_eQTL)]<-0
F2_eQTLs_vs_RatGTEx_GeneFrequencyFull$NumberOfTissuesWHC<-F2_eQTLs_vs_RatGTEx_GeneFrequencyFull$NumberOfTissues+F2_eQTLs_vs_RatGTEx_GeneFrequencyFull$HC_eQTL

table(F2_eQTLs_vs_RatGTEx_GeneFrequencyFull$NumberOfTissues)
# 0    1    2    3    4    5    6    7    8    9   10   11   12 
# 565 5009 2921 2057 1442 1045  907  713  600  617  547  555  142 

table(F2_eQTLs_vs_RatGTEx_GeneFrequencyFull$NumberOfTissuesWHC)
#   1    2    3    4    5    6    7    8    9   10   11   12   13 
# 5011 2855 2096 1546 1119  934  760  653  565  573  510  404   94 

write.csv(F2_eQTLs_vs_RatGTEx_GeneFrequencyFull, "F2_eQTLs_vs_RatGTEx_GeneFrequencyFull.csv")

table(F2_eQTLs_vs_RatGTEx_GeneFrequencyFull$HC_eQTL, F2_eQTLs_vs_RatGTEx_GeneFrequencyFull$NumberOfTissues)

#     0    1    2    3    4    5    6    7    8    9   10   11   12
# 0    0 4446 2292 1467  956  633  522  375  315  280  236  199   48
# 1  565  563  629  590  486  412  385  338  285  337  311  356   94


#That's interesting - some of our top genes are definitely eGenes in a variety of tissues.


#################
#Maybe the magnitude could be a source of validation?
#Like, are we identifying similar genes as *strongly* regulated by genetics.

library(plyr)

F2_eQTLs_toJoin<-F2_eQTLs
colnames(F2_eQTLs_toJoin)[-1]<-paste("F2", colnames(F2_eQTLs_toJoin)[-1], sep="_")

F2_eQTLs_Vs_RatGTEx_byEnsembl<-join(F2_eQTLs_toJoin, RatGTEx, by="gene_id", match="all", type="inner")
str(F2_eQTLs_Vs_RatGTEx_byEnsembl)
#'data.frame':	34742 obs. of  34 variables:

F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl<-F2_eQTLs_Vs_RatGTEx_byEnsembl[F2_eQTLs_Vs_RatGTEx_byEnsembl$tissue=="Brain",]
str(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl)
#'data.frame':	5874 obs. of  34 variables:

sum(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$variant_id==F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_variant_id)
#[1] 2
#lololol

plot(abs(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_log2_aFC_alt_ref)~abs(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$log2_aFC))
#Not terrible
cor(abs(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_log2_aFC_alt_ref),abs(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$log2_aFC), use="pairwise.complete.obs")
#[1] 0.6905708

cor(abs(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_log2_aFC_alt_ref),abs(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$log2_aFC), use="pairwise.complete.obs", method="spearman")
#[1] 0.4550737
#so...is that actually evidence supporting our eQTLs?
#Or does that mean that the genes that are highly variable in one dataset are highly variable in another?

#I wonder if the distance from transcription start site is similar for the variants
#(like if they are vaguely in the same neighborhood)
#And or the bp

plot(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_pos~F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$pos)
#Strongly linear, but that doesn't mean very much because they were already joined by gene and these are cis-Eqtls

plot(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_tss_distance~F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$tss_distance)
#Basically no relationship that I can see - plus sign shape.
#If anything, I would say that the tss_distances that are most extreme in one dataset are not often present in the other.
#Which I guess would make sense becuase those are the tss_distances that the the most likely to be randomly assigned eVariants within LD blocks

cor(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_tss_distance,  F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$tss_distance)
#[1] 0.05427136
#Maybe we could filter by that though?

hist(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_tss_distance)
hist(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$tss_distance)

head(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_tss_distance)
F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_tss_distance_rounded<-round(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_tss_distance, digits=-5)
head(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_tss_distance_rounded)
F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_tss_distance_rounded[c(1:100)]

F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$tss_distance_rounded<-round(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$tss_distance, digits=-5)

sum(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$tss_distance_rounded==F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_tss_distance_rounded)
#[1] 570
#Ouch.

F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl_SimilarTss<-F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl[F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$tss_distance_rounded==F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_tss_distance_rounded,]

plot(abs(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl_SimilarTss$F2_log2_aFC_alt_ref)~abs(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl_SimilarTss$log2_aFC))
#Well, there's a positive correlation, but it seems driven by a few extreme effects (aFC>1)

cor(abs(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl_SimilarTss$F2_log2_aFC_alt_ref),abs(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl_SimilarTss$log2_aFC))
#[1] 0.7884664

cor(abs(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl_SimilarTss$F2_log2_aFC_alt_ref),abs(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl_SimilarTss$log2_aFC), method="spearman")
#[1] 0.5310822

plot(abs(F2_log2_aFC_alt_ref)~abs(log2_aFC), data=F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl_SimilarTss[F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl_SimilarTss$F2_rank==1,])
#not a big difference

write.csv(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl_SimilarTss, "F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl_SimilarTss.csv")


plot(log10(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_pval_nominal)~log10(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$pval_nominal))
#Definitely... not strongly correlated, especially since we already filtered for significance.

cor(log10(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$F2_pval_nominal),log10(F2_eQTLs_Vs_RatGTEx_Brain_byEnsembl$pval_nominal), use="pairwise.complete.obs")
#[1] 0.1451754


#I'm going to try to increase the # of comparisons by including other brain-related tissue.
table(RatGTEx$tissue)

length(unique(RatGTEx$gene_id[RatGTEx$tissue!="Eye"& RatGTEx$tissue!="Liver" & RatGTEx$tissue!="Adipose"]))
#[1] 13097
#that's ~3,000 more than what we had with just brain

F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl<-F2_eQTLs_Vs_RatGTEx_byEnsembl[F2_eQTLs_Vs_RatGTEx_byEnsembl$tissue!="Eye"& F2_eQTLs_Vs_RatGTEx_byEnsembl$tissue!="Liver" & F2_eQTLs_Vs_RatGTEx_byEnsembl$tissue!="Adipose",]
str(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl)
#'data.frame':	25844 obs. of  34 variables:
#'Lots of repetition now.

sum(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$variant_id==F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$F2_variant_id, na.rm=TRUE)
#[1] 9
#lololol
F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$gene_name[which(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$variant_id==F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$F2_variant_id)]
#[1] "Thbs2" "Myh7b" "Ptpro" "Galt"  "Cox16" "Brcc3" "Stap2" "Iqce"  "Pomk"  

plot(abs(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$F2_log2_aFC_alt_ref)~abs(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$log2_aFC))
#Messier - many more points now.
cor(abs(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$F2_log2_aFC_alt_ref),abs(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$log2_aFC), use="pairwise.complete.obs")
#[1] 0.6729271

cor(abs(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$F2_log2_aFC_alt_ref),abs(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$log2_aFC), use="pairwise.complete.obs", method="spearman")
#[1] 0.5105281

pdf("Plot_RatGTEx_BraineQTL_aFC_vs_F2_HCeQTL_aFC.pdf", width=5, height=5)
plot(abs(log2_aFC)~abs(F2_log2_aFC_alt_ref), data=F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl, ylab="RatGTEx Brain eQTLs: aFC", xlab="F2 Hippocampal eQTLs: aFC")
TrendLine<-lm(abs(log2_aFC)~abs(F2_log2_aFC_alt_ref), data=F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl)
abline(TrendLine, col="red", lwd=3)
dev.off()

#so...is that actually evidence supporting our eQTLs?
#Or does that mean that the genes that are highly variable in one dataset are highly variable in another?

#I wonder if the distance from transcription start site is similar for the variants
#(like if they are vaguely in the same neighborhood)
#And or the bp

plot(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$F2_tss_distance~F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$tss_distance)
#Basically no relationship that I can see - plus sign shape.
#If anything, I would say that the tss_distances that are most extreme in one dataset are not often present in the other.
cor(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$F2_tss_distance,  F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$tss_distance, use="pairwise.complete.obs")
#[1] 0.04928481
#Maybe we could filter by that though?

F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$F2_tss_distance_rounded<-round(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$F2_tss_distance, digits=-5)

F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$tss_distance_rounded<-round(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$tss_distance, digits=-5)

sum(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$tss_distance_rounded==F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$F2_tss_distance_rounded, na.rm=TRUE)
#[1] 2738
#Better

F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_SimilarTss<-F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl[F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$tss_distance_rounded==F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$F2_tss_distance_rounded,]

plot(abs(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_SimilarTss$F2_log2_aFC_alt_ref)~abs(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_SimilarTss$log2_aFC))
#Well, there's a positive correlation
#Still driven by  extreme effects (aFC>1), but there are many more points now.

pdf("Plot_RatGTEx_BraineQTL_SimilarTss_aFC_vs_F2_HCeQTL_aFC.pdf", width=5, height=5)
plot(abs(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_SimilarTss$log2_aFC)~abs(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_SimilarTss$F2_log2_aFC_alt_ref), ylab="Nearby RatGTEx Brain eQTLs: aFC", xlab="F2 Hippocampal eQTLs: aFC")
TrendLine<-lm(abs(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_SimilarTss$log2_aFC)~abs(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_SimilarTss$F2_log2_aFC_alt_ref))
abline(TrendLine, col="red", lwd=3)
dev.off()

cor(abs(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_SimilarTss$F2_log2_aFC_alt_ref),abs(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_SimilarTss$log2_aFC), use="pairwise.complete.obs")
#[1] 0.7749707

cor(abs(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_SimilarTss$F2_log2_aFC_alt_ref),abs(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_SimilarTss$log2_aFC), use="pairwise.complete.obs", method="spearman")
#[1] 0.6202615

plot(abs(F2_log2_aFC_alt_ref)~abs(log2_aFC), data=F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_SimilarTss[F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_SimilarTss$F2_rank==1,])
#that is a little bit prettier
summary.lm(lm(abs(F2_log2_aFC_alt_ref)~abs(log2_aFC), data=F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_SimilarTss[F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_SimilarTss$F2_rank==1,]))
# Multiple R-squared:  0.6423,	Adjusted R-squared:  0.6422 
# F-statistic:  4363 on 1 and 2430 DF,  p-value: < 2.2e-16
sqrt(0.6423)
#[1] 0.8014362
#Not a huge jump though.

write.csv(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_SimilarTss, "F2_eQTLs_Vs_RatGTEx_AllBrainRegions_byEnsembl_SimilarTss.csv")

#So a few comments:
#1) eGenes that are strongly related to genetic variation in our study also tend to be strongly related to genetic variation in other eQTL studies
#2) Some of the eGenes that are strong eQTLs and differentially expressed in our study *are not* eGenes with similar tss variants in RatGTEx:
#This is particularly notable for some of our strongest eQTLs: 
#e.g., Fcrl2, AABR07071904.1, Rpl30, Fanci, Ucp2, Plekhb1, Tmem144, Afg3l1, Fzd6, Spg7, Idh1, Sun2, C2cd3, Bphl, Trhr,  Mfge8, 
#Also: Lsr,  Lipt2,  Vps9d1, Ist1, Sp3, Ttc30a1, Rarres2, Gimap5
#Or have much weaker eQTLs (and/or only one of them) in RatGTEx: Bmp4, Tmco5a
#Wdr93, Pex11a, Nqo2, Mcee have eVariants with similar tss and similar strength aFC

#Quite notably, these eGenes in our F2s aren't found *at all* as eGenes in the brain in RatGTEx:
F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$F2_gene_name[
F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl$F2_gene_name%in%c("Fcrl2", "AABR07071904.1", "Rpl30", "Fanci", "Ucp2", "Plekhb1", "Tmem144", "Afg3l1", "Fzd6", "Spg7", "Idh1", "Sun2", "C2cd3", "Bphl", "Trhr",  "Mfge8", "Lsr",  "Lipt2",  "Vps9d1", "Ist1", "Sp3", "Ttc30a1", "Rarres2", "Gimap5")]
#Only "Sp3", "Ttc30a1", "Idh1", "Spg7" aren't in the larger database, all of the other genes are in there just with loci located farther away from the tss site we identified.
#I wonder how far away the others are from the variants we identified?

write.csv(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl, "F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl.csv")

#Interesting. 
#So Fcrl2, Nqo2, Trhr, Fzd6, Wdr93, Tmco5a, Pex11a, Bphl are associated with huge aFC for other variants. (esp. Fcrl2!)
#Fanci, Tmem144, Mfge8, Afg3l1, Rpl30 are respectable
#The size of the effects for Ucp2 are kind of wussy, but the p-value for one of them is huge
#Same goes for Plekhb1
#AABR07071904.1 is super wussy

############

#Just one more check: How much are the extreme aFCs driven by low-level expressed genes?

F2_Residuals<-read.csv("ResidualsMatrix_wMean.csv", header=TRUE, stringsAsFactors = FALSE, row.names=1)
str(F2_Residuals)
#'data.frame':	14056 obs. of  245 variables:
#'
F2_MeanExpression<-apply(F2_Residuals, 1, mean)
head(F2_MeanExpression)

sum(names(F2_MeanExpression)%in%F2_eQTLs_More$gene_id)
#[1] 13547

F2_MeanExpression_DF<-data.frame(gene_id=names(F2_MeanExpression), F2_MeanExpression)

F2_MeanExpression_DF$F2_SD_Expression<-apply(F2_Residuals, 1, sd)

F2_eQTLs_vs_MeanExpression<-join(F2_eQTLs, F2_MeanExpression_DF, by="gene_id")
str(F2_eQTLs_vs_MeanExpression)
plot(abs(F2_eQTLs_vs_MeanExpression$log2_aFC_alt_ref)~F2_eQTLs_vs_MeanExpression$F2_MeanExpression)
#Definitely negatively correlated

pdf("Plot_aFC_vs_MeanExpression.pdf", height=5, width=5)
plot(abs(F2_eQTLs_vs_MeanExpression$log2_aFC_alt_ref)~F2_eQTLs_vs_MeanExpression$F2_MeanExpression, xlab="F2 Mean Expression", ylab="Abs(Log2 aFC)")
dev.off()

cor(abs(F2_eQTLs_vs_MeanExpression$log2_aFC_alt_ref), F2_eQTLs_vs_MeanExpression$F2_MeanExpression)
#[1] -0.3914107
cor(abs(F2_eQTLs_vs_MeanExpression$log2_aFC_alt_ref), F2_eQTLs_vs_MeanExpression$F2_MeanExpression, method="spearman")
#[1] -0.4941592

pdf("Plot_aFC_vs_SD_Expression.pdf", height=5, width=5)
plot(abs(F2_eQTLs_vs_MeanExpression$log2_aFC_alt_ref)~F2_eQTLs_vs_MeanExpression$F2_SD_Expression, xlab="F2 Expression StDev", ylab="Abs(Log2 aFC)")
dev.off()

plot(abs(F2_eQTLs_vs_MeanExpression$log2_aFC_alt_ref)~F2_eQTLs_vs_MeanExpression$F2_SD_Expression)
#positively correlated, although that makes sense if we're assuming that genetic variation is the main source of variation for these eGenes
cor(abs(F2_eQTLs_vs_MeanExpression$log2_aFC_alt_ref), F2_eQTLs_vs_MeanExpression$F2_SD_Expression)
#[1] 0.772334
cor(abs(F2_eQTLs_vs_MeanExpression$log2_aFC_alt_ref), F2_eQTLs_vs_MeanExpression$F2_SD_Expression, method="spearman")
#[1] 0.6147486

summary.lm(lm(abs(F2_eQTLs_vs_MeanExpression$log2_aFC_alt_ref)~F2_eQTLs_vs_MeanExpression$F2_MeanExpression+F2_eQTLs_vs_MeanExpression$F2_SD_Expression))
# Call:
#   lm(formula = abs(F2_eQTLs_vs_MeanExpression$log2_aFC_alt_ref) ~ 
#        F2_eQTLs_vs_MeanExpression$F2_MeanExpression + F2_eQTLs_vs_MeanExpression$F2_SD_Expression)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.7622 -0.1178 -0.0164  0.0839  6.0665 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                  -0.216944   0.015536 -13.964  < 2e-16 ***
#   F2_eQTLs_vs_MeanExpression$F2_MeanExpression  0.013092   0.002396   5.464 4.84e-08 ***
#   F2_eQTLs_vs_MeanExpression$F2_SD_Expression   2.084347   0.025692  81.129  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3371 on 5934 degrees of freedom
# Multiple R-squared:  0.5985,	Adjusted R-squared:  0.5984 
# F-statistic:  4423 on 2 and 5934 DF,  p-value: < 2.2e-16

#That's not much stronger than the prediction from just SD_Expression
#Over 12 units of Log2FC mean expression you would only expect
0.13*12
#[1] 1.56 difference in Log2FC
#O.k., maybe that is a lot.
#Whereas SD expression really only ranges up to 3, so at max 6 difference in Log2FC


#Is it only the aFC that shows this pattern?

plot(log10(F2_eQTLs_vs_MeanExpression$pval_nominal)~F2_eQTLs_vs_MeanExpression$F2_MeanExpression)
#not much relationship, if anything, the strongest p-values are from the middle-expressed genes
plot(log10(F2_eQTLs_vs_MeanExpression$pval_beta)~F2_eQTLs_vs_MeanExpression$F2_MeanExpression)
#not much relationship, if anything, the strongest p-values are from the middle-expressed genes

plot(log10(F2_eQTLs_vs_MeanExpression$pval_nominal)~F2_eQTLs_vs_MeanExpression$F2_SD_Expression)
plot(log10(F2_eQTLs_vs_MeanExpression$pval_beta)~F2_eQTLs_vs_MeanExpression$F2_SD_Expression)
#No relationship

#Next question: 
#Does low level gene expression mediate the ability of the strength of our eGenes (aFC) to predict RatGTEx eGenes?

F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression<-join(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl, F2_MeanExpression_DF, by="gene_id", type="inner", match="all")
dim(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression)
#[1] 25844    38

colnames(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression)

summary.lm(lm(abs(F2_log2_aFC_alt_ref)~abs(log2_aFC), data=F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression))
# Call:
#   lm(formula = abs(F2_log2_aFC_alt_ref) ~ abs(log2_aFC), data = F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.6801 -0.1626 -0.0523  0.0864  6.4749 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.107835   0.003546   30.41   <2e-16 ***
#   abs(log2_aFC) 0.690501   0.004722  146.24   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4666 on 25842 degrees of freedom
# Multiple R-squared:  0.4528,	Adjusted R-squared:  0.4528 
# F-statistic: 2.139e+04 on 1 and 25842 DF,  p-value: < 2.2e-16

summary.lm(lm(abs(F2_log2_aFC_alt_ref)~abs(log2_aFC)+F2_MeanExpression, data=F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression))
# Call:
#   lm(formula = abs(F2_log2_aFC_alt_ref) ~ abs(log2_aFC) + F2_MeanExpression, 
#      data = F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.1896 -0.1673 -0.0366  0.1031  6.2071 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        0.365743   0.008539   42.83   <2e-16 ***
#   abs(log2_aFC)      0.613619   0.005177  118.54   <2e-16 ***
#   F2_MeanExpression -0.048740   0.001474  -33.06   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4571 on 25841 degrees of freedom
# Multiple R-squared:  0.475,	Adjusted R-squared:  0.475 
# F-statistic: 1.169e+04 on 2 and 25841 DF,  p-value: < 2.2e-16

#Nice. Abs(log2FC) is almost as predictive after controlling for MeanExpression

summary.lm(lm(abs(F2_log2_aFC_alt_ref)~abs(log2_aFC)+F2_MeanExpression+F2_SD_Expression, data=F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression))
# Call:
#   lm(formula = abs(F2_log2_aFC_alt_ref) ~ abs(log2_aFC) + F2_MeanExpression + 
#        F2_SD_Expression, data = F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.8937 -0.1248 -0.0214  0.0845  5.9140 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       -0.150014   0.007662 -19.580  < 2e-16 ***
#   abs(log2_aFC)      0.160398   0.005266  30.460  < 2e-16 ***
#   F2_MeanExpression  0.005894   0.001211   4.867 1.14e-06 ***
#   F2_SD_Expression   1.815562   0.013730 132.236  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.353 on 25840 degrees of freedom
# Multiple R-squared:  0.6869,	Adjusted R-squared:  0.6869 
# F-statistic: 1.89e+04 on 3 and 25840 DF,  p-value: < 2.2e-16

#And Abs(log2FC) is also still predictive (although less so) after controlling for SD Expression
#Although the interpretation for this finding is up for debate...

hist(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression$F2_MeanExpression)
summary(F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression$F2_MeanExpression)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.593   3.234   4.755   4.611   6.090  11.403 

pdf("Plot_RatGTEx_BraineQTL_aFC_vs_F2_HCeQTL_aFC_ColorLowExpression.pdf", width=5, height=5)
plot(abs(log2_aFC)~abs(F2_log2_aFC_alt_ref), data=F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression, ylab="RatGTEx Brain eQTLs: aFC", xlab="F2 Hippocampal eQTLs: aFC")
points(abs(log2_aFC)~abs(F2_log2_aFC_alt_ref), data=F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression[F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression$F2_MeanExpression<2,], col="blue")
points(abs(log2_aFC)~abs(F2_log2_aFC_alt_ref), data=F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression[F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression$F2_MeanExpression>2,], col="black")
TrendLine<-lm(abs(log2_aFC)~abs(F2_log2_aFC_alt_ref), data=F2_eQTLs_Vs_RatGTEx_Brain2_byEnsembl_wMeanExpression)
abline(TrendLine, col="red", lwd=3)
dev.off()

F2_eQTLs_vs_RatGTEx_GeneFrequencyFull_wMeanExpression<-join(F2_eQTLs_vs_RatGTEx_GeneFrequencyFull, F2_MeanExpression_DF, by="gene_id", type="left")
dim(F2_eQTLs_vs_RatGTEx_GeneFrequencyFull_wMeanExpression)
#[1] 17120     8

boxplot(F2_eQTLs_vs_RatGTEx_GeneFrequencyFull_wMeanExpression$F2_MeanExpression~F2_eQTLs_vs_RatGTEx_GeneFrequencyFull_wMeanExpression$NumberOfTissues)
plot(F2_eQTLs_vs_RatGTEx_GeneFrequencyFull_wMeanExpression$F2_MeanExpression~F2_eQTLs_vs_RatGTEx_GeneFrequencyFull_wMeanExpression$NumberOfTissues)
#Wow - there's really no relationship between F2 Mean Expression in the HC and the number of tissues that have an eGene in RatGTEx
#I was really expecting that genes with higher expression levels would be more likely to show up repeatedly (just because they would be more likely to be expressed everywhere)

##################################################

#Ok, trying again - Daniel sent us the document with all of the nominally significant cis-eQTLs for each gene (not just the FDR-significant, conditionally independent ones)

F2_eQTL_AllSig<-read.delim("f2.cis_qtl_signif.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_eQTL_AllSig)
#'data.frame':	7930908 obs. of  11 variables:

sum(F2_eQTL_AllSig$variant_id%in%RatGTEx$variant_id)
#[1] 119414
#There we go.

colnames(F2_eQTL_AllSig)
#Oh weird -for some reason ENSEMBL ID is listed as "phenotype_id" - probably because the eQTL and QTL software are the same maybe?
length(unique(F2_eQTL_AllSig$phenotype_id))
#[1] 5351
length(unique(F2_eQTL_AllSig$variant_id))

F2_eQTL_AllSig$eGene_eVariant<-paste(F2_eQTL_AllSig$phenotype_id, F2_eQTL_AllSig$variant_id, sep="_")

head(F2_eQTL_AllSig$eGene_eVariant)

#################

#Peeking at whether the coding variants identified in the Chitre et al. paper are in this database:

Chitre_PutativeCausalCodingVariants<-c("chr1:94610123","chr7:83403144")

write.csv(F2_eQTL_AllSig[F2_eQTL_AllSig$variant_id%in%Chitre_PutativeCausalCodingVariants,], "F2_eQTL_AllSig_ChitreCodingVariants.csv")

#It looks like the Plekhf1 coding variant is in LD with the ENSRNOG00000052237 eVariant, and also serves an eVariant for ENSRNOG00000052237
#It is just a little bit less significant (e-39 vs. e-53) - I'm not sure how meaningful that is in a sample size this big.

write.csv(RatGTEx[RatGTEx$variant_id%in%Chitre_PutativeCausalCodingVariants,], "RatGTEx_ChitreCodingVariants.csv")
#these variants aren't represented as eQTLs in RatGTEx

#Distance to Plekhf1 start:
94836296-94609377
#226919
#Distance to Plekhf1 stop:
94836296-94610854
#225442

#visualizing the eVariants for ENSRNOG00000052237:

str(F2_eQTL_AllSig[F2_eQTL_AllSig$phenotype_id=="ENSRNOG00000052237",])
#'data.frame':	2633 obs. of  11 variables:

pdf("Plot_AABR07071904_eQTL_Block_wPlekhf1missense.pdf", height=5, width=8)
plot(slope~tss_distance, data=F2_eQTL_AllSig[F2_eQTL_AllSig$phenotype_id=="ENSRNOG00000052237",], main="AABR07071904.1 (ENSRNOG00000052237) eQTL")
abline(h=0, col="grey")#no slope
abline(v=167111, col="red")#our top eVariant
abline(v=0, col="grey")#tss for ENSRNOG00000052237
abline(v=-226173, col="blue")#Plekhf1 missense coding variant
abline(v=-226919, col="blue", lty=2)#Plekhf1 start site
abline(v=-225442, col="blue", lty=2)#Plekhf1 stop site
dev.off()


#Rnor6 Ensembl coordinates for ENSRNOG00000052237: 94,836,296 - 94,840,706
#Rnor6 Ensembl coordinates for Plekhf1: 94,609,377 - 94,610,854 (-)

#Daniel is going to go back and confirm this using more traditional methods of identifying LD blocks.


#################

#Comparing with RatGTEx

RatGTEx$eGene_eVariant<-paste(RatGTEx$gene_id, RatGTEx$variant_id, sep="_")

sum(F2_eQTL_AllSig$eGene_eVariant%in%RatGTEx$eGene_eVariant)
#[1] 10321
#Much better

F2_eQTL_AllSig_toJoin<-F2_eQTL_AllSig
colnames(F2_eQTL_AllSig_toJoin)[-11]<-paste("F2", colnames(F2_eQTL_AllSig_toJoin)[-11], sep="_")

F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant<-join(F2_eQTL_AllSig_toJoin, RatGTEx, by="eGene_eVariant", match="all", type="inner")
str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant)
#'data.frame':	10685 obs. of  28 variables:

write.csv(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant, "F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant.csv")

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$variant_id))
#[1] 9843
length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id))
#[1] 2885
#So approximately half of the eGenes in our F2 eQTL database.

#There isn't information about the alt/ref alleles in this output file
#So I guess we'll take it on faith that both databases are using the same alt/ref?


#I should probably output one with our top e-genes to take a peek at alternative explanations for the effects we are seeing:

TopEGenesOfInterest<-c("ENSRNOG00000021053",
                       "ENSRNOG00000052237",
                       "ENSRNOG00000016327",
                       "ENSRNOG00000016689",
                       "ENSRNOG00000012357",
                       "ENSRNOG00000015003",
                       "ENSRNOG00000026514",
                       "ENSRNOG00000017510",
                       "ENSRNOG00000016906",
                       "ENSRNOG00000017608",
                       "ENSRNOG00000018627",
                       "ENSRNOG00000017854",
                       "ENSRNOG00000004660",
                       "ENSRNOG00000028904",
                       "ENSRNOG00000026994",
                       "ENSRNOG00000038330",
                       "ENSRNOG00000010081",
                       "ENSRNOG00000016164",
                       "ENSRNOG00000017820",
                       "ENSRNOG00000017577",
                       "ENSRNOG00000015144")


write.csv(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id%in%TopEGenesOfInterest,], "F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant_JustTopGenesOfInterest.csv")

#I'm curious as to whether the top eGenes we identified tend to have larger slopes in Rat GTEx than other eQTLs on average
hist(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id%in%TopEGenesOfInterest], breaks=12, xlim=c(-5,5))
#probably easier to visualize in absolute terms:
hist(abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id%in%TopEGenesOfInterest]), breaks=10, xlim=c(0,5))

summary(abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id%in%TopEGenesOfInterest]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2001  0.6084  0.9613  0.9383  1.2733  1.7252

hist(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5], breaks=100, xlim=c(-5,5))
#abs:
hist(abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5]), breaks=100, xlim=c(0,5))

summary(abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5]), breaks=100, xlim=c(0,5))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.05595 0.37505 0.56190 0.66776 0.82174 4.93677 

#Yep, the distribution for our top eGenes is definitely right shifted.
#Cool.


#########################

plot(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_slope~F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)
#There's an extreme outlier in the RatGTEx data that is making it hard to visualize

summary(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -600132.0      -0.6      -0.2     -56.0       0.5     245.3

plot(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope>-6e05,])
#Apparently there are a few other extreme outliers in the RatGTEx dataset too

plot(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<10,])
#Well, it's not exactly beautiful, but there's a positive correlation.

summary.lm(lm(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<10,]))
# Call:
#   lm(formula = slope ~ F2_slope, data = F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope) < 
#                                                                                    10, ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -8.2892 -0.2774  0.0281  0.2993  9.5566 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.054097   0.006679   -8.10  6.1e-16 ***
#   F2_slope     0.711672   0.009329   76.29  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6879 on 10647 degrees of freedom
# Multiple R-squared:  0.3534,	Adjusted R-squared:  0.3534 
# F-statistic:  5820 on 1 and 10647 DF,  p-value: < 2.2e-16

sqrt(0.3534)
#[1] 0.5944746

#... and not a wussy correlation at that.

#I should probably divide these figures/stats up by RatGTEx tissue type:

table(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue)
# Adipose   BLA   Brain     Eye      IL     LHb   Liver    NAcc   NAcc2     OFC      PL     PL2 
# 1221    1063    1855     118     729     626     997     609     853     749     757    1108 


#Maybe first take a peek at all Brain tissue types?

#There are a few outliers

hist(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Adipose" & F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Liver" &  F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Eye"], breaks=100)

str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Adipose" & F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Liver" &  F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Eye" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)>5,])
#'data.frame':	46 obs. of  28 variables:

unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Adipose" & F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Liver" &  F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Eye"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)>5], breaks=100)
# [1] "ENSRNOG00000019968" "ENSRNOG00000048043" "ENSRNOG00000033837" "ENSRNOG00000007542"
# [5] "ENSRNOG00000007338" "ENSRNOG00000018109" "ENSRNOG00000003147" "ENSRNOG00000046493"
# [9] "ENSRNOG00000051344"

str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Adipose" & F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Liver" &  F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Eye" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,])
#'data.frame':	8303 obs. of  28 variables:


library(scales)
pdf("Scatterplot_F2_eQTLs_AllSig_Vs_RatGTEX_AllBrainTissueTypes.pdf", height=5, width=5)
plot(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Adipose" & F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Liver" &  F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Eye"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,], col = alpha(1, 0.2), xlab="F2 Hippocampal eQTLs: Slope", ylab="Rat GTEx All Brain Tissue eQTLs: Slope", main="All Significant F2 Hippocampal eQTLs")
TrendLine<-lm(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Adipose" & F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Liver" &  F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Eye"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,])
abline(TrendLine, col=2, lwd=3)
dev.off()

pdf("Scatterplot_F2_eQTLs_AllSig_Vs_RatGTEX_AllBrainTissueTypes_Lighter.pdf", height=7, width=4)
plot(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Adipose" & F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Liver" &  F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Eye"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,], col = alpha(1, 0.05), xlab="F2 Hippocampal eQTLs: Slope", ylab="Rat GTEx All Brain Tissue eQTLs: Slope", main="All Significant F2 Hippocampal eQTLs")
TrendLine<-lm(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Adipose" & F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Liver" &  F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue!="Eye"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = slope ~ F2_slope, data = F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue != 
#                                                                                    "Adipose" & F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue != 
#                                                                                    "Liver" & F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue != 
#                                                                                    "Eye" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope) < 
#                                                                                 5, ])
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -5.5668 -0.2299  0.0417  0.2828  4.9208 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.068253   0.007045  -9.688   <2e-16 ***
#   F2_slope     0.811077   0.009672  83.860   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6406 on 8301 degrees of freedom
# Multiple R-squared:  0.4586,	Adjusted R-squared:  0.4586 
# F-statistic:  7032 on 1 and 8301 DF,  p-value: < 2.2e-16

sqrt(0.4586)
#[1] 0.6772001



#Brain Hemisphere (this is an entire hemisphere and has the biggest sample size in RatGTEx)


str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Brain",])
#'data.frame':	1855 obs. of  28 variables:
#'
length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Brain"]))
#[1] 1556

#Zoomed in on this tissue, there seems to be another extreme outlier:

str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Brain" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)>5,])
#'data.frame':	34 obs. of  28 variables:
#Huh, that's actually kind-of a lot of data points. Are they all the same eGenes?
length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Brain" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)>5]))
#[1] 5
#Ah - so I guess maybe the points are right on top of each other.

str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Brain" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,])
#'data.frame':	1821 obs. of  28 variables:

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Brain" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5]))
#[1] 1552
#So lots of genes represented in what's left.

plot(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Brain" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,])
#ooh - that's pretty.

library(scales)
pdf("Scatterplot_F2_eQTLs_AllSig_Vs_RatGTEX_Brain.pdf", height=5, width=5)
plot(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Brain" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,], col = alpha(1, 0.4), xlab="F2 Hippocampal eQTLs: Slope", ylab="Rat GTEx Brain eQTLs: Slope", main="All Significant F2 Hippocampal eQTLs")
TrendLine<-lm(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Brain" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)
# Call:
#   lm(formula = slope ~ F2_slope, data = F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue == 
#                                                                                    "Brain" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope) < 
#                                                                                    5, ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -5.5894 -0.2469  0.1071  0.3681  5.1613 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.16652    0.01843  -9.035   <2e-16 ***
#   F2_slope     0.90164    0.02640  34.157   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7848 on 1819 degrees of freedom
# Multiple R-squared:  0.3908,	Adjusted R-squared:  0.3904 
# F-statistic:  1167 on 1 and 1819 DF,  p-value: < 2.2e-16

sqrt(0.3908)
#[1] 0.62514

#It is probably worthwhile to run a non-parametric version of the stats that includes the extreme outliers.

cor(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Brain"], F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Brain"], method="spearman")
#[1] 0.6280302

#BLA


str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="BLA",])
#'data.frame':	1063 obs. of  28 variables:

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="BLA"]))
#[1] 969

#BLA has a few outlier points too:
summary(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="BLA"])
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -3.94755 -0.58835 -0.26779 -0.01578  0.52732  5.97247 

#Outlier
str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="BLA" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)>5,])
#'data.frame':	8 obs. of  28 variables:

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="BLA" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)>5]))
#[1] 1

#What's left:
str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="BLA" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,])
#'data.frame':	1055 obs. of  28 variables:

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="BLA" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5]))
#[1] 968
#968 genes

pdf("Scatterplot_F2_eQTLs_AllSig_Vs_RatGTEX_BLA.pdf", height=5, width=5)
plot(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="BLA" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,], col = alpha(1, 0.4), xlab="F2 Hippocampal eQTLs: Slope", ylab="Rat GTEx BLA eQTLs: Slope", main="All Significant F2 Hippocampal eQTLs")
TrendLine<-lm(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="BLA" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = slope ~ F2_slope, data = F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue == 
#                                                                                    "BLA" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope) < 
#                                                                                    5, ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -3.10617 -0.21664  0.01554  0.20110  2.55552 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.02451    0.01714   -1.43    0.153    
# F2_slope     0.81038    0.02345   34.56   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5556 on 1053 degrees of freedom
# Multiple R-squared:  0.5315,	Adjusted R-squared:  0.531 
# F-statistic:  1194 on 1 and 1053 DF,  p-value: < 2.2e-16
sqrt(0.5315)
#[1] 0.7290405

#It is probably worthwhile to run a non-parametric version of the stats that includes the extreme outliers.
cor(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="BLA"], F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="BLA"], method="spearman")
#[1] 0.7522498

#PL

str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="PL",])
#'data.frame':	757 obs. of  28 variables:

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="PL"]))
#[1] 739

pdf("Scatterplot_F2_eQTLs_AllSig_Vs_RatGTEX_PL.pdf", height=5, width=5)
plot(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="PL",], col = alpha(1, 0.4), xlab="F2 Hippocampal eQTLs: Slope", ylab="Rat GTEx PL eQTLs: Slope", main="All Significant F2 Hippocampal eQTLs")
TrendLine<-lm(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="PL",])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = slope ~ F2_slope, data = F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue == 
#                                                                                    "PL", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.0623 -0.2865 -0.0140  0.2554  3.2421 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.02365    0.02127  -1.112    0.267    
# F2_slope     0.77817    0.02940  26.465   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5838 on 755 degrees of freedom
# Multiple R-squared:  0.4812,	Adjusted R-squared:  0.4806 
# F-statistic: 700.4 on 1 and 755 DF,  p-value: < 2.2e-16

sqrt(0.4812)
#[1] 0.6936858

#It is probably worthwhile to run a non-parametric version of the stats that includes the extreme outliers.
cor(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="PL"], F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="PL"], method="spearman")
#[1] 0.6944644

#PL2

str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="PL2",])
#'data.frame':	1108 obs. of  28 variables:

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="PL2"]))
#[1] 956

#This one has a few outliers too:

summary(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="PL2"])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -8.9345 -0.6576 -0.3234 -0.2228  0.4989  3.9606 

#Outlier
str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="PL2" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)>5,])
#'data.frame':	1 obs. of  28 variables:

#What's left
str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="PL2" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,])
#'data.frame':	1107 obs. of  28 variables:

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="PL2" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5]))
#[1] 956
#956 unique genes

pdf("Scatterplot_F2_eQTLs_AllSig_Vs_RatGTEX_PL2.pdf", height=5, width=5)
plot(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="PL2"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,], col = alpha(1, 0.4), xlab="F2 Hippocampal eQTLs: Slope", ylab="Rat GTEx PL2 eQTLs: Slope", main="All Significant F2 Hippocampal eQTLs")
TrendLine<-lm(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="PL2"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)
# Call:
#   lm(formula = slope ~ F2_slope, data = F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue == 
#                                                                                    "PL2" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope) < 
#                                                                                    5, ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.5785 -0.2013  0.0994  0.3401  3.1983 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.12276    0.02238  -5.486  5.1e-08 ***
#   F2_slope     0.99373    0.02953  33.649  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7389 on 1105 degrees of freedom
# Multiple R-squared:  0.5061,	Adjusted R-squared:  0.5057 
# F-statistic:  1132 on 1 and 1105 DF,  p-value: < 2.2e-16

sqrt(0.5061)
#[1] 0.7114071

#It is probably worthwhile to run a non-parametric version of the stats that includes the extreme outliers.
cor(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="PL2"], F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="PL2"], method="spearman")
#[1] 0.7692336



#IL
#No outliers

str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="IL",])
#'data.frame':	729 obs. of  28 variables:

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="IL"]))
#[1] 713

pdf("Scatterplot_F2_eQTLs_AllSig_Vs_RatGTEX_IL.pdf", height=5, width=5)
plot(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="IL",], col = alpha(1, 0.4), xlab="F2 Hippocampal eQTLs: Slope", ylab="Rat GTEx IL eQTLs: Slope", main="All Significant F2 Hippocampal eQTLs")
TrendLine<-lm(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="IL",])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = slope ~ F2_slope, data = F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue == 
#                                                                                    "IL", ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.04098 -0.27392 -0.00898  0.27711  2.56728 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.01570    0.02123   -0.74     0.46    
# F2_slope     0.77288    0.02902   26.63   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5727 on 727 degrees of freedom
# Multiple R-squared:  0.4938,	Adjusted R-squared:  0.4932 
# F-statistic: 709.3 on 1 and 727 DF,  p-value: < 2.2e-16

sqrt(0.4938)
#[1] 0.702709

#It is probably worthwhile to run a non-parametric version of the stats:
cor(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="IL"], F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="IL"], method="spearman")
#[1] 0.6917697


#OFC 
#No outliers

str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="OFC" ,])
#'data.frame':	749 obs. of  28 variables:

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="OFC"]))
#[1] 733
#unique genes

pdf("Scatterplot_F2_eQTLs_AllSig_Vs_RatGTEX_OFC.pdf", height=5, width=5)
plot(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="OFC",], col = alpha(1, 0.4), xlab="F2 Hippocampal eQTLs: Slope", ylab="Rat GTEx OFC eQTLs: Slope", main="All Significant F2 Hippocampal eQTLs")
TrendLine<-lm(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="OFC",])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)

# Call:
#   lm(formula = slope ~ F2_slope, data = F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue == 
#                                                                                    "OFC", ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.54201 -0.27724 -0.00806  0.27279  2.70095 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.002418   0.019844  -0.122    0.903    
# F2_slope     0.778499   0.027293  28.523   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5422 on 747 degrees of freedom
# Multiple R-squared:  0.5213,	Adjusted R-squared:  0.5207 
# F-statistic: 813.6 on 1 and 747 DF,  p-value: < 2.2e-16

sqrt(0.5213)
#[1] 0.7220111

#It is probably worthwhile to run a non-parametric version of the stats:
cor(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="OFC"], F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="OFC"], method="spearman")
#[1] 0.7044997

##########################

#Presumably the eQTLs would be more different in the other (non-cortical) brain regions?

table(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue)

Adipose   BLA   Brain     Eye     IL     LHb   Liver    NAcc   NAcc2     OFC      PL     PL2 
1221    1063    1855     118     729     626     997     609     853     749     757    1108 

#LHb
#This one has outlier datapoints
#The slopes are generally smaller - I wonder if I should use a different definition of outlier

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="LHb"]))
#[1] 617
#unique genes

str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="LHb"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)>5 ,])
#'data.frame':	2 obs. of  28 variables:
#outliers

#What's left
str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="LHb"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5 ,])
#'data.frame':	624 obs. of  28 variables:

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="LHb"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5]))
#[1] 616
#unique genes

pdf("Scatterplot_F2_eQTLs_AllSig_Vs_RatGTEX_LHb.pdf", height=5, width=5)
plot(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="LHb" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,], col = alpha(1, 0.4), xlab="F2 Hippocampal eQTLs: Slope", ylab="Rat GTEx LHb eQTLs: Slope", main="All Significant F2 Hippocampal eQTLs")
TrendLine<-lm(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="LHb"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)
# Call:
#   lm(formula = slope ~ F2_slope, data = F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue == 
#                                                                                    "LHb" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope) < 
#                                                                                    5, ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.8266 -0.2934 -0.0248  0.2559  4.6090 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.0003926  0.0227434   0.017    0.986    
# F2_slope    0.7017496  0.0306172  22.920   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5657 on 622 degrees of freedom
# Multiple R-squared:  0.4579,	Adjusted R-squared:  0.457 
# F-statistic: 525.3 on 1 and 622 DF,  p-value: < 2.2e-16

sqrt(0.4579)
#[1] 0.6766831

cor(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="LHb"], F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="LHb"], method="spearman")
#[1] 0.6678689


#Nucleus Accumbens:

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="NAcc"]))
#[1] 604
#unique genes

str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="NAcc"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)>5 ,])
#'data.frame':	1 obs. of  28 variables:
#outliers

#What's left
str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="NAcc"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5 ,])
#'data.frame':	608 obs. of  28 variables:

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="NAcc"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5]))
#[1] 603
#unique genes

pdf("Scatterplot_F2_eQTLs_AllSig_Vs_RatGTEX_NAcc.pdf", height=5, width=5)
plot(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="NAcc" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,], col = alpha(1, 0.4), xlab="F2 Hippocampal eQTLs: Slope", ylab="Rat GTEx NAcc eQTLs: Slope", main="All Significant F2 Hippocampal eQTLs")
TrendLine<-lm(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="NAcc"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)
# Call:
#   lm(formula = slope ~ F2_slope, data = F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue == 
#                                                                                    "NAcc" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope) < 
#                                                                                    5, ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.22114 -0.31484  0.01365  0.31677  2.46454 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.06124    0.02379  -2.574   0.0103 *  
#   F2_slope     0.73712    0.03295  22.373   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5862 on 606 degrees of freedom
# Multiple R-squared:  0.4524,	Adjusted R-squared:  0.4515 
# F-statistic: 500.6 on 1 and 606 DF,  p-value: < 2.2e-16

sqrt(0.452)
#[1] 0.6723095

cor(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="NAcc"], F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="NAcc"], method="spearman")
#[1] 0.6484169

#NAcc2

str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="NAcc2"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)>5 ,])
#'data.frame':	0 obs. of  28 variables:
#no outliers

#What's left
str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="NAcc2"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5 ,])
#'data.frame':	853 obs. of  28 variables:

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="NAcc2"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5]))
#[1] 804
#unique genes

pdf("Scatterplot_F2_eQTLs_AllSig_Vs_RatGTEX_NAcc2.pdf", height=5, width=5)
plot(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="NAcc2" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,], col = alpha(1, 0.4), xlab="F2 Hippocampal eQTLs: Slope", ylab="Rat GTEx NAcc2 eQTLs: Slope", main="All Significant F2 Hippocampal eQTLs")
TrendLine<-lm(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="NAcc2"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)
# Call:
#   lm(formula = slope ~ F2_slope, data = F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue == 
#                                                                                    "NAcc2" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope) < 
#                                                                                    5, ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.98062 -0.16552  0.01238  0.18844  1.37149 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.02188    0.01391  -1.573    0.116    
# F2_slope     0.60633    0.01863  32.548   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4061 on 851 degrees of freedom
# Multiple R-squared:  0.5545,	Adjusted R-squared:  0.554 
# F-statistic:  1059 on 1 and 851 DF,  p-value: < 2.2e-16

sqrt(0.5545)
#[1] 0.7446476

#Non-parametric:
cor(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="NAcc2"], F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="NAcc2"], method="spearman")
#[1] 0.7405411


#So basically as strong as the cortical tissue, especially for the NAcc derived from the larger sample size (n=188)


#I wonder if the comparison is weaker with tissue that is very different from the brain (although the n's are also super different, so this comparison may be difficult to interpret)

#Liver

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Liver"]))
#[1] 915
#unique genes

str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Liver"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)>5 ,])
#'data.frame':	5 obs. of  28 variables:
#outliers

#What's left
str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Liver"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5 ,])
#'data.frame':	992 obs. of  28 variables:

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Liver" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5]))
#[1] 912
#unique genes

pdf("Scatterplot_F2_eQTLs_AllSig_Vs_RatGTEX_Liver.pdf", height=5, width=5)
plot(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Liver" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,], col = alpha(1, 0.4), xlab="F2 Hippocampal eQTLs: Slope", ylab="Rat GTEx Liver eQTLs: Slope", main="All Significant F2 Hippocampal eQTLs")
TrendLine<-lm(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Liver"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)
# Call:
#   lm(formula = slope ~ F2_slope, data = F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue == 
#                                                                                    "Liver" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope) < 
#                                                                                    5, ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.15408 -0.29772 -0.03768  0.30395  2.68336 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.04083    0.01504  -2.714  0.00676 ** 
#   F2_slope     0.22265    0.02276   9.782  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4734 on 990 degrees of freedom
# Multiple R-squared:  0.08813,	Adjusted R-squared:  0.08721 
# F-statistic: 95.69 on 1 and 990 DF,  p-value: < 2.2e-16

sqrt(0.08813)
#[1] 0.296867

#Non-parametric:
cor(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Liver"], F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Liver"], method="spearman")
#[1] 0.3061357

#Adipose:

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Adipose"]))
#[1] 1101
#unique genes

str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Adipose"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)>5 ,])
#'data.frame':	4 obs. of  28 variables:
#outliers

#What's left
str(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Adipose"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5 ,])
#'data.frame':	1217 obs. of  28 variables:

length(unique(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_phenotype_id[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Adipose" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5]))
#[1] 1097
#unique genes

pdf("Scatterplot_F2_eQTLs_AllSig_Vs_RatGTEX_Adipose.pdf", height=5, width=5)
plot(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Adipose" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,], col = alpha(1, 0.4), xlab="F2 Hippocampal eQTLs: Slope", ylab="Rat GTEx Adipose eQTLs: Slope", main="All Significant F2 Hippocampal eQTLs")
TrendLine<-lm(slope~F2_slope, data=F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Adipose"& abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope)<5,])
abline(TrendLine, col=2, lwd=3)
dev.off()

summary.lm(TrendLine)
# Call:
#   lm(formula = slope ~ F2_slope, data = F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue == 
#                                                                                    "Adipose" & abs(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope) < 
#                                                                                    5, ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.4561 -0.3226 -0.0186  0.3484  1.9077 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.01511    0.01380  -1.095    0.274    
# F2_slope     0.27399    0.02074  13.213   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4804 on 1215 degrees of freedom
# Multiple R-squared:  0.1256,	Adjusted R-squared:  0.1249 
# F-statistic: 174.6 on 1 and 1215 DF,  p-value: < 2.2e-16

sqrt(0.1256)
#[1] 0.3544009

cor(F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Adipose"], F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$F2_slope[F2_eQTLs_AllSig_Vs_RatGTEx_byGeneVariant$tissue=="Adipose"], method="spearman")
#[1] 0.3478177

#Well that's super fun. I'll need to be careful not to over-interpret the results though, especially since both of these tissue come from the same cohort/study.




