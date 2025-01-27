#Making some quick sex difference volcano plots
#Megan Hagenauer
#Dec 12 2024

setwd("~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/Manuscript/NextSubmission/FrontiersInMolNeuro/Revisions/NewSexDiff_SupplTable")

SexDiffLineageLoco_DE<-read.csv("ToMakeSupplTable_JustLineageLocomotorSexDiff.csv", header=TRUE, stringsAsFactors = FALSE)

str(SexDiffLineageLoco_DE)
#I should reverse the direction of effect in the F2 output to match the F0 (i.e., make males the reference instead of females)

colnames(SexDiffLineageLoco_DE)

SexDiffLineageLoco_DE$F2_Log2FC_FvsM<-SexDiffLineageLoco_DE$F2_Log2FC_MvsF*(-1)

SexDiffLineageLoco_DE$F2_Tstat_FvsM<-SexDiffLineageLoco_DE$F2_Tstat_MvsF*(-1)

#Here's the volcano plot code from the F0s to parallel for the F2s
colnames(SexDiffLineageLoco_DE)

colnames(SexDiffLineageLoco_DE)[1]<-"ENSEMBL"



#Volcano plot for Sex in the F2s

tiff("VolcanoPlot_F2Sex_noInteraction_wGold_FDR_0point1.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(SexDiffLineageLoco_DE, plot(F2_Log2FC_FvsM, -log10(F2_Pvalue_MvsF), pch=19, main="Sex", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-6,6), ylim=c(0,12), xlab="Log2FC", ylab="-log10(p-value)"))
with(subset(SexDiffLineageLoco_DE, abs(F2_Log2FC_FvsM)>1), points(F2_Log2FC_FvsM, -log10(F2_Pvalue_MvsF), pch=19, col="red", cex=0.6))
with(subset(SexDiffLineageLoco_DE, F2_FDR_MvsF<0.1), points(F2_Log2FC_FvsM, -log10(F2_Pvalue_MvsF), pch=19, col="green", cex=0.6))
with(subset(SexDiffLineageLoco_DE, abs(F2_Log2FC_FvsM)>1 & F2_FDR_MvsF<0.1), points(F2_Log2FC_FvsM, -log10(F2_Pvalue_MvsF), pch=19, col="gold", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.1", "both"), col=c("red", "green", "gold"), pch=19, cex=1)
dev.off()


#Volcano plot for Log2FC in the F2s

#I had to remove the xlim from this code because the Log2FC is on a completely different scale for LocoScore (xlim=c(-6,6),)
#Likewise, it really doesn't make sense to color code Log2FC>1 on this scale.

tiff("VolcanoPlot_F2LocoScore_noInteraction_wGold_FDR_0point1.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(SexDiffLineageLoco_DE, plot(F2_Log2FC_LocoScore, -log10(F2_Pvalue_LocoScore), pch=19, main="LocoScore", cex.lab=1.8, cex.main=2, cex=0.6,  ylim=c(0,12), xlab="Log2FC", ylab="-log10(p-value)"))
# with(subset(SexDiffLineageLoco_DE, abs(F2_Log2FC_LocoScore)>1), points(F2_Log2FC_LocoScore, -log10(F2_Pvalue_LocoScore), pch=19, col="red", cex=0.6))
with(subset(SexDiffLineageLoco_DE, F2_FDR_LocoScore<0.1), points(F2_Log2FC_LocoScore, -log10(F2_Pvalue_LocoScore), pch=19, col="green", cex=0.6))
# with(subset(SexDiffLineageLoco_DE, abs(F2_Log2FC_LocoScore)>1 & F2_FDR_LocoScore<0.1), points(F2_Log2FC_LocoScore, -log10(F2_Pvalue_LocoScore), pch=19, col="gold", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.1", "both"), col=c("red", "green", "gold"), pch=19, cex=1)
dev.off()

#How does a unit of Log2FC compare to bHR vs. bLR?

setwd("~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/F2_eQTLs/Input")

Behavior<-read.csv("Akil_HRLR_F0_F1_F2_PhenotypeData_ForR.csv", header = TRUE, stringsAsFactors = FALSE)

str(Behavior)

setwd("~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/Manuscript/NextSubmission/FrontiersInMolNeuro/Revisions/NewSexDiff_SupplTable")

Behavior$Sex[Behavior$Sex=="Female"]<-"F"

Behavior$Sex[Behavior$Sex=="Male"]<-"M"

Behavior$Phenotype[Behavior$Phenotype=="I1"]<-NA

Behavior$Phenotype[Behavior$Phenotype=="I2"]<-NA

Behavior$Phenotype[is.na(Behavior$Phenotype)]<-""


Behavior$Sex_AsFactor<-as.factor(Behavior$Sex)
levels(Behavior$Sex_AsFactor)
#[1] "F" "M"

#Setting the reference levels for factors:
Behavior$Sex_AsFactor<-relevel(Behavior$Sex_AsFactor, ref="M")
levels(Behavior$Sex_AsFactor)
#[1] "M" "F"


Behavior$LearningClassification_AsFactor<-as.factor(Behavior$LearningClassification)
levels(Behavior$LearningClassification_AsFactor)
#[1] ""   "GT" "IN" "ST"

#Setting the reference levels for factors:
Behavior$LearningClassification_AsFactor<-relevel(Behavior$LearningClassification_AsFactor, ref="ST")
levels(Behavior$LearningClassification_AsFactor)
#[1] "ST" ""   "GT" "IN"

Behavior$GenPheno<-paste(Behavior$Study_Generation,Behavior$Phenotype,sep = " ")

str(Behavior)

Behavior$SexGenPheno<-paste(Behavior$Sex, Behavior$GenPheno, sep="_")

tapply(Behavior$Total_LocoScore[Behavior$RNASeq=="Elaine"], INDEX=Behavior$SexGenPheno[Behavior$RNASeq=="Elaine"], mean)
# F_F0 HR    F_F0 LR      F_F2     M_F0 HR    M_F0 LR      M_F2  
# 2446.33333   40.50000  574.66400 2033.16667   52.66667  500.00800

#So mean differences between F0 bHR and bLR:
#Females:
2446.33333-40.50000
#[1] 2405.833

#Males:
2033.16667-52.66667
#[1] 1980.5

#Average across sex:
(2405.833+1980.5)/2
#[1] 2193.167


#So bHR/bLR differences are similar to a difference of 2193 LocoScore
#So hypothetically we could multiply the LocoScore betas for the differential expression output by 2193 to make them more interpretable/comparable with phenotype

hist(SexDiffLineageLoco_DE$F2_Log2FC_LocoScore)

SexDiffLineageLoco_DE$F2_Log2FC_LocoScore_UnitsHRLR<-SexDiffLineageLoco_DE$F2_Log2FC_LocoScore*2193.167

hist(SexDiffLineageLoco_DE$F2_Log2FC_LocoScore_UnitsHRLR)
#That's actually pretty cool.

tapply(Behavior$Total_LocoScore[Behavior$RNASeq=="Elaine"], INDEX=Behavior$SexGenPheno[Behavior$RNASeq=="Elaine"], max)
# F_F0 HR F_F0 LR   F_F2  M_F0 HR M_F0 LR   M_F2  
# 2949      69    1461    2320     101    1319 

tapply(Behavior$Total_LocoScore[Behavior$RNASeq=="Elaine"], INDEX=Behavior$SexGenPheno[Behavior$RNASeq=="Elaine"], min)
# F_F0 HR F_F0 LR   F_F2  M_F0 HR M_F0 LR   M_F2  
# 2152      23      75    1533      30      40 

#Compared to the range of the F2s:

#females:
1461-75
#[1] 1386

#males:
1319-40
#[1] 1279

tapply(Behavior$Total_LocoScore[Behavior$RNASeq=="Elaine"], INDEX=Behavior$SexGenPheno[Behavior$RNASeq=="Elaine"], sd)
# F_F0 HR   F_F0 LR     F_F2    M_F0 HR   M_F0 LR     M_F2  
# 278.42390  17.82975 326.14097 301.60532  27.30324 272.39371 

#Effect size (d) for bHR/bLR differences in locomotor activity:

#In females:
2405.833/sqrt(((278.42390^2)+(17.82975^2))/2)
#[1] 12.1951

#In males:
1980.5/sqrt(((301.60532^2)+(27.30324^2))/2)
#[1] 9.248655

#Size of the range in the F2s in units of standard deviation:

#females:
1386/326.14097
#[1] 4.249696

#males:
1279/272.39371
#[1] 4.695409

#... and that's the range.  Whereas the "average" deviation from the mean is the stdev - so in general, about half of the subjects are going to fall within +/-1 stdev units.  


#I'm going to make a volcano plot using the bHR/bLR version of LocoScore:

tiff("VolcanoPlot_F2LocoScoreUnitsHRLR_noInteraction_wGold_FDR_0point1.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(SexDiffLineageLoco_DE, plot(F2_Log2FC_LocoScore_UnitsHRLR, -log10(F2_Pvalue_LocoScore), pch=19, main="LocoScore", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-6,6), ylim=c(0,12), xlab="Log2FC", ylab="-log10(p-value)"))
with(subset(SexDiffLineageLoco_DE, abs(F2_Log2FC_LocoScore_UnitsHRLR)>1), points(F2_Log2FC_LocoScore_UnitsHRLR, -log10(F2_Pvalue_LocoScore), pch=19, col="red", cex=0.6))
with(subset(SexDiffLineageLoco_DE, F2_FDR_LocoScore<0.1), points(F2_Log2FC_LocoScore_UnitsHRLR, -log10(F2_Pvalue_LocoScore), pch=19, col="green", cex=0.6))
with(subset(SexDiffLineageLoco_DE, abs(F2_Log2FC_LocoScore_UnitsHRLR)>1 & F2_FDR_LocoScore<0.1), points(F2_Log2FC_LocoScore_UnitsHRLR, -log10(F2_Pvalue_LocoScore), pch=19, col="gold", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.1", "both"), col=c("red", "green", "gold"), pch=19, cex=1)
dev.off()


#Another interesting question might be how well sex differences in gene expression in the F0s predicted sex differences in gene expression in the F2s (i.e., how well predictions hold up when testing *the same* concept instead of related concepts of phenotype and locoscore)



cor.test(SexDiffLineageLoco_DE$F0_Log2FC_FvsM, SexDiffLineageLoco_DE$F2_Log2FC_FvsM, use="pairwise.complete", method="spearman")

# Spearman's rank correlation rho
# 
# data:  SexDiffLineageLoco_DE$F0_Log2FC_FvsM and SexDiffLineageLoco_DE$F2_Log2FC_FvsM
# S = 3.6673e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.07290722

#Wow, that's actually pretty weak.
#I guess because the effect sizes are smaller, so the power to detect the effects in the smaller sample is pretty wussy.

#for comparison:
#bHR/bLR F0 vs. Meta-analysis estimate was rho=0.29782
#bHR/bLR F0 vs. F2 Locoscore was rho= -0.2289803 
#bHR/bLR F0 vs. F2 time immobile was rho=0.2931336 

#Let's visualize it:
plot(SexDiffLineageLoco_DE$F2_Log2FC_FvsM~SexDiffLineageLoco_DE$F0_Log2FC_FvsM, xlab="Log2FC: F0 F vs. M", ylab="Log2FC: F2: F vs. M")
#LOLOLOL - yeah, there's some super extreme values there, which are definitely consistent between datasets.

#what about everyone in between?

plot(SexDiffLineageLoco_DE$F2_Log2FC_FvsM~SexDiffLineageLoco_DE$F0_Log2FC_FvsM, xlim=c(-3,3), ylim=c(-1.5,1.5), xlab="Log2FC: F0 F vs. M", ylab="Log2FC: F2: F vs. M")
#Not much relationship - it mostly just looks like a ball.

#that means that the pearson correlation is probably much higher than spearman:
cor.test(SexDiffLineageLoco_DE$F0_Log2FC_FvsM, SexDiffLineageLoco_DE$F2_Log2FC_FvsM, use="pairwise.complete", method="pearson")

#yep
#Pearson's product-moment correlation
# data:  SexDiffLineageLoco_DE$F0_Log2FC_FvsM and SexDiffLineageLoco_DE$F2_Log2FC_FvsM
# t = 121.14, df = 13337, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.7156138 0.7317763
# sample estimates:
#       cor 
# 0.7237943 


#Let's make an RRHO too:

library(RRHO)

TempDF<-SexDiffLineageLoco_DE[is.na(SexDiffLineageLoco_DE$F0_Tstat_FvsM)==FALSE & is.na(SexDiffLineageLoco_DE$F2_Tstat_FvsM)==FALSE, ]

str(TempDF)
#'data.frame':	13339 obs. of  24 variables:

list1<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$F0_Tstat_FvsM)

list2<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$F2_Tstat_FvsM)

getwd()
#[1] "/Users/hagenaue/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/Manuscript/NextSubmission/FrontiersInMolNeuro/Revisions/NewSexDiff_SupplTable"

RRHO(list1, list2, labels=c("F0: F vs. M", "F2: F vs. M"), plots=TRUE, alternative="two.sided", outputdir="/Users/hagenaue/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/Manuscript/NextSubmission/FrontiersInMolNeuro/Revisions/NewSexDiff_SupplTable", BY=TRUE, log10.ind=TRUE)

#So there is convergence at the extremes, but not much else.



#Another quick question:

#Do the significant sex difference genes overlap with the bHR/bLR genes?

SexDiffLineageLoco_DE$ENSEMBL[which((SexDiffLineageLoco_DE$F0_FDR_FvsM<0.10 | SexDiffLineageLoco_DE$F2_FDR_MvsF<0.1) & SexDiffLineageLoco_DE$F0_FDR_bLRvsbHR<0.1)]
# [1] "ENSRNOG00000020456" "ENSRNOG00000019708" "ENSRNOG00000020906"
# [4] "ENSRNOG00000012278" "ENSRNOG00000016294" "ENSRNOG00000017063"
# [7] "ENSRNOG00000012749" "ENSRNOG00000033564" "ENSRNOG00000008138"
# [10] "ENSRNOG00000013720" "ENSRNOG00000004327" "ENSRNOG00000009068"
# [13] "ENSRNOG00000034191" "ENSRNOG00000005365" "ENSRNOG00000002764"
# [16] "ENSRNOG00000026783" "ENSRNOG00000020751" "ENSRNOG00000018070"
# [19] "ENSRNOG00000018938" "ENSRNOG00000050040"

#Yep, there's definitely some (which makes sense, given how many there are)

#I should probably pull in the whole list for comparison:

str(HPC_DEResults)

#Double-checking that the order is the same:
sum(HPC_DEResults$ENSEMBL==SexDiffLineageLoco_DE$ENSEMBL)

#Apparently there are two blank rows that read in accidentally with the HPC_DEResults

HPC_DEResults[c(13787:13788),]
# X ENSEMBL GENENAME..Rnor6.Ensembl.v.88. ENTREZID..Rnor6.Ensembl.v.88.
# 13787 NA                                                                  NA
# 13788 NA                                                                  NA

HPC_DEResults2<-HPC_DEResults[-c(13787:13788),]

sum(HPC_DEResults2$ENSEMBL==SexDiffLineageLoco_DE$ENSEMBL)
#[1] 13786
sum(HPC_DEResults2$ENSEMBL!=SexDiffLineageLoco_DE$ENSEMBL)
#[1] 0

#Good.

SexDiffLineageLoco_DE$ENSEMBL[which((SexDiffLineageLoco_DE$F0_FDR_FvsM<0.10 | SexDiffLineageLoco_DE$F2_FDR_MvsF<0.1) & HPC_DEResults2$bHR_bLR_DEG==1)]

# [1] "ENSRNOG00000013351" "ENSRNOG00000017826" "ENSRNOG00000018808"
# [4] "ENSRNOG00000014166" "ENSRNOG00000016055" "ENSRNOG00000015413"
# [7] "ENSRNOG00000023521" "ENSRNOG00000021053" "ENSRNOG00000021067"
# [10] "ENSRNOG00000021102" "ENSRNOG00000009932" "ENSRNOG00000018421"
# [13] "ENSRNOG00000017854" "ENSRNOG00000013290" "ENSRNOG00000020456"
# [16] "ENSRNOG00000012806" "ENSRNOG00000024144" "ENSRNOG00000020533"
# [19] "ENSRNOG00000019708" "ENSRNOG00000020299" "ENSRNOG00000020906"
# [22] "ENSRNOG00000021164" "ENSRNOG00000012278" "ENSRNOG00000009953"
# [25] "ENSRNOG00000009822" "ENSRNOG00000016550" "ENSRNOG00000020027"
# [28] "ENSRNOG00000021157" "ENSRNOG00000011559" "ENSRNOG00000013409"
# [31] "ENSRNOG00000010972" "ENSRNOG00000015425" "ENSRNOG00000010466"
# [34] "ENSRNOG00000046905" "ENSRNOG00000031173" "ENSRNOG00000006228"
# [37] "ENSRNOG00000005615" "ENSRNOG00000006003" "ENSRNOG00000057556"
# [40] "ENSRNOG00000016294" "ENSRNOG00000006444" "ENSRNOG00000017063"
# [43] "ENSRNOG00000017172" "ENSRNOG00000010897" "ENSRNOG00000012172"
# [46] "ENSRNOG00000014204" "ENSRNOG00000021234" "ENSRNOG00000021269"
# [49] "ENSRNOG00000008124" "ENSRNOG00000015668" "ENSRNOG00000014869"
# [52] "ENSRNOG00000007060" "ENSRNOG00000010381" "ENSRNOG00000013642"
# [55] "ENSRNOG00000000145" "ENSRNOG00000019609" "ENSRNOG00000000122"
# [58] "ENSRNOG00000033984" "ENSRNOG00000012749" "ENSRNOG00000009183"
# [61] "ENSRNOG00000021526" "ENSRNOG00000022694" "ENSRNOG00000011550"
# [64] "ENSRNOG00000019981" "ENSRNOG00000010208" "ENSRNOG00000039666"
# [67] "ENSRNOG00000005387" "ENSRNOG00000003313" "ENSRNOG00000030418"
# [70] "ENSRNOG00000029885" "ENSRNOG00000005124" "ENSRNOG00000047088"
# [73] "ENSRNOG00000004198" "ENSRNOG00000059714" "ENSRNOG00000011489"
# [76] "ENSRNOG00000033564" "ENSRNOG00000030273" "ENSRNOG00000007545"
# [79] "ENSRNOG00000027124" "ENSRNOG00000005034" "ENSRNOG00000005931"
# [82] "ENSRNOG00000004403" "ENSRNOG00000030334" "ENSRNOG00000012295"
# [85] "ENSRNOG00000012886" "ENSRNOG00000009592" "ENSRNOG00000009756"
# [88] "ENSRNOG00000004659" "ENSRNOG00000059605" "ENSRNOG00000026643"
# [91] "ENSRNOG00000011582" "ENSRNOG00000008904" "ENSRNOG00000010688"
# [94] "ENSRNOG00000008138" "ENSRNOG00000010524" "ENSRNOG00000007442"
# [97] "ENSRNOG00000008487" "ENSRNOG00000020721" "ENSRNOG00000055524"
# [100] "ENSRNOG00000048411" "ENSRNOG00000013852" "ENSRNOG00000017258"
# [103] "ENSRNOG00000016567" "ENSRNOG00000002212" "ENSRNOG00000002866"
# [106] "ENSRNOG00000002643" "ENSRNOG00000003924" "ENSRNOG00000003069"
# [109] "ENSRNOG00000013720" "ENSRNOG00000004327" "ENSRNOG00000009642"
# [112] "ENSRNOG00000009935" "ENSRNOG00000003453" "ENSRNOG00000000033"
# [115] "ENSRNOG00000051548" "ENSRNOG00000009068" "ENSRNOG00000003687"
# [118] "ENSRNOG00000034191" "ENSRNOG00000034015" "ENSRNOG00000006019"
# [121] "ENSRNOG00000003654" "ENSRNOG00000019713" "ENSRNOG00000019729"
# [124] "ENSRNOG00000020140" "ENSRNOG00000003269" "ENSRNOG00000005365"
# [127] "ENSRNOG00000002810" "ENSRNOG00000012155" "ENSRNOG00000002764"
# [130] "ENSRNOG00000026783" "ENSRNOG00000059618" "ENSRNOG00000008689"
# [133] "ENSRNOG00000015941" "ENSRNOG00000018247" "ENSRNOG00000020751"
# [136] "ENSRNOG00000020792" "ENSRNOG00000009474" "ENSRNOG00000004604"
# [139] "ENSRNOG00000048712" "ENSRNOG00000036664" "ENSRNOG00000011704"
# [142] "ENSRNOG00000011621" "ENSRNOG00000018651" "ENSRNOG00000011947"
# [145] "ENSRNOG00000018070" "ENSRNOG00000017511" "ENSRNOG00000018938"
# [148] "ENSRNOG00000050040" "ENSRNOG00000016585" "ENSRNOG00000049949"
# [151] "ENSRNOG00000052111" "ENSRNOG00000022500" "ENSRNOG00000012654"
# [154] "ENSRNOG00000037957" "ENSRNOG00000024066" "ENSRNOG00000027469"
# [157] "ENSRNOG00000020084" "ENSRNOG00000018796" "ENSRNOG00000016971"
# [160] "ENSRNOG00000017745" "ENSRNOG00000012428" "ENSRNOG00000042274"
# [163] "ENSRNOG00000018784" "ENSRNOG00000015150" "ENSRNOG00000026994"
# [166] "ENSRNOG00000019022" "ENSRNOG00000000609" "ENSRNOG00000050655"
# [169] "ENSRNOG00000043436" "ENSRNOG00000000975" "ENSRNOG00000000906"
# [172] "ENSRNOG00000047329" "ENSRNOG00000001241" "ENSRNOG00000001317"
# [175] "ENSRNOG00000001426" "ENSRNOG00000000931" "ENSRNOG00000001173"

sum(HPC_DEResults2$bHR_bLR_DEG==1)
#[1] 1045

175/1045
#[1] 0.1674641

175/1679 
#[1] 0.1042287

#Useful.
#I should add that to our suppl table.

#Visually double-checking that the order of the formatted supplement matches the current dataframe:
head(SexDiffLineageLoco_DE$ENSEMBL)
tail(SexDiffLineageLoco_DE$ENSEMBL)
#Yep, looks good.

write.csv(data.frame(HPC_DEResults2$ENSEMBL, HPC_DEResults2$bHR_bLR_DEG), "bHR_bLR_DEG.csv")

save.image("~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/HRLR_NONHPC/bHR_bLR_NonHPC_Meta_AndSexDiff.RData")
