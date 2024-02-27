
################


#Making a manhattan plot of the SMR results that we used in our analysis

#https://www.r-graph-gallery.com/101_Manhattan_plot.html

# Oh - I can see why Daniel did Bonferonni now as the correction.
# The manhattan plots are p-value
# But the FDR correction will have a slightly different p-value for each variable that meets threshold for "significant"

#Alright, it's actually simpler than expected:
#For adults, only locoscore and distance traveled survive FDR correction.
#And only one gene survives for FDR for distance traveled (AABR07071904.1)
#For LocoScore: FDR<0.05=p<0.000119 and FDR<0.10=p<0.000546
#For juveniles, only of distance traveled survives FDR correction.
#no genes with FDR<0.05, FDR<0.10=p<9.05e-05

#So we could really just make the manhattan plots focus on the SMR results for a single behavior for each age group and have everything important represented

library("scales")

library(qqman)
#Citation appreciated but not required:
# Turner, (2018). qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots. Journal of Open Source Software, 3(25), 731, https://doi.org/10.21105/joss.00731.

SMR_OnlyRelevant_ForPlotting<-SMR[which(SMR$trait=="total_loco_score"|SMR$trait=="epm_time_immobile_s"|SMR$trait=="epm_distance_traveled"|SMR$trait=="epm_percent_time_open_arm"|SMR$trait=="pca_index_days_6and7"),]
colnames(SMR_OnlyRelevant_ForPlotting)
head(SMR_OnlyRelevant_ForPlotting$variant_id)
TempSMR<-do.call(rbind, strsplit(SMR_OnlyRelevant_ForPlotting$variant_id, split=":", fixed=TRUE))
head(TempSMR)
TempSMR2<-do.call(rbind, strsplit(TempSMR[,1], split="chr", fixed=TRUE))
head(TempSMR2)
SMR_OnlyRelevant_ForPlotting$CHR<-as.numeric(TempSMR2[,2])
SMR_OnlyRelevant_ForPlotting$LOC<-as.numeric(TempSMR[,2])
table(SMR_OnlyRelevant_ForPlotting$CHR)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20 
# 3885 1855 2160 1585 1710 1455 1750 1730 1350 2700 1050 1025  745 1000 1180 1050  790  955  780  930 
#Only autosomes. Perfect.
rm(TempSMR, TempSMR2)
colnames(SMR_OnlyRelevant_ForPlotting)

setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/F2_eQTLs")

pdf("ManhattanPlot_SMR_Pval.pdf", width=16, height=5)
manhattan(SMR_OnlyRelevant_ForPlotting, chr="CHR", bp="LOC", snp="variant_id", p="p_SMR", genomewideline=FALSE, suggestiveline=FALSE)
dev.off()

#Booger - the variant id denotes all of the SMR results for a variant, whether they are the strong result or not.
#I'll need to make a new column that combines variant id and trait.
SMR_OnlyRelevant_ForPlotting$variant_trait_forplotting<-paste(SMR_OnlyRelevant_ForPlotting$variant_id, SMR_OnlyRelevant_ForPlotting$trait, sep="_")


#just SMR results for LocoScore marked off with a simple FDR<0.05 and FDR<0.10:

SMR_OnlyLocoScore_forPlotting<-SMR[which(SMR$trait=="total_loco_score"),]
colnames(SMR_OnlyLocoScore_forPlotting)
head(SMR_OnlyLocoScore_forPlotting$variant_id)
TempSMR<-do.call(rbind, strsplit(SMR_OnlyLocoScore_forPlotting$variant_id, split=":", fixed=TRUE))
head(TempSMR)
TempSMR2<-do.call(rbind, strsplit(TempSMR[,1], split="chr", fixed=TRUE))
head(TempSMR2)
SMR_OnlyLocoScore_forPlotting$CHR<-as.numeric(TempSMR2[,2])
SMR_OnlyLocoScore_forPlotting$LOC<-as.numeric(TempSMR[,2])
table(SMR_OnlyLocoScore_forPlotting$CHR)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 777 371 432 317 342 291 350 346 270 540 210 205 149 200 236 210 158 191 156 186 

SMR_OnlyLocoScore_forPlotting$variant_trait_forplotting<-paste(SMR_OnlyLocoScore_forPlotting$variant_id, SMR_OnlyLocoScore_forPlotting$trait, sep="_")

pdf("ManhattanPlot_SMR_LocoScore_Pval_ConvergentColored_FDR05andFDR10_forAllCriteria.pdf", width=10, height=5)
manhattan(SMR_OnlyLocoScore_forPlotting, chr="CHR", bp="LOC", snp="variant_trait_forplotting", p="p_SMR", genomewideline=-log10(0.0001186), suggestiveline=-log10(0.000546), highlight=c(
  "chr1:89443807_total_loco_score",
  "chr1:95003407_total_loco_score",
  "chr1:141037045_total_loco_score",
  "chr1:141210512_total_loco_score",
  "chr1:141294154_total_loco_score",
  "chr1:141312442_total_loco_score",
  "chr1:141313585_total_loco_score",
  "chr1:164561826_total_loco_score",
  "chr1:165794789_total_loco_score",
  "chr1:166318198_total_loco_score",
  "chr1:166392350_total_loco_score",
  "chr7:77982760_total_loco_score",
  "chr19:41897300_total_loco_score",
  "chr19:55921978_total_loco_score",
  "chr19:55989619_total_loco_score",
  "chr19:56219141_total_loco_score"))
dev.off()

#No Mcee for this set of results - It's significant SMR and DE results don't line up  "chr1:125425048_total_loco_score"
#I'm pretty sure that is a false negative - it's just barely over the threshold - but it would be too complicated to argue for it in the space we have in the manuscript.

######################

###ooooh - I just realized that we have SMR results based on the juvenile QTLs too
# let's make a Manhattan plot for those - might serve as good additional evidence, even though the behavior is a little different.

table(SMR$trait)
SMR_OnlyRelevant_ForPlotting_Juvenile<-SMR[which(SMR$trait=="of_distance_traveled"|SMR$trait=="of_percent_time_in_center"|SMR$trait=="of_time_immobile_s"),]
colnames(SMR_OnlyRelevant_ForPlotting_Juvenile)
head(SMR_OnlyRelevant_ForPlotting_Juvenile$variant_id)
TempSMR<-do.call(rbind, strsplit(SMR_OnlyRelevant_ForPlotting_Juvenile$variant_id, split=":", fixed=TRUE))
head(TempSMR)
TempSMR2<-do.call(rbind, strsplit(TempSMR[,1], split="chr", fixed=TRUE))
head(TempSMR2)
SMR_OnlyRelevant_ForPlotting_Juvenile$CHR<-as.numeric(TempSMR2[,2])
SMR_OnlyRelevant_ForPlotting_Juvenile$LOC<-as.numeric(TempSMR[,2])
table(SMR_OnlyRelevant_ForPlotting_Juvenile$CHR)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20 
# 2331 1113 1296  951 1026  873 1050 1038  810 1620  630  615  447  600  708  630  474  573  468  558
#Only autosomes. Perfect.
rm(TempSMR, TempSMR2)
colnames(SMR_OnlyRelevant_ForPlotting_Juvenile)

plot.new()
hist(SMR_OnlyRelevant_ForPlotting_Juvenile$p_SMR)

sum(SMR_OnlyRelevant_ForPlotting_Juvenile$p_SMR=="inf")
#[1] NA
sum(is.na(SMR_OnlyRelevant_ForPlotting_Juvenile$p_SMR))
#[1] 3

SMR_OnlyRelevant_ForPlotting_Juvenile<-SMR_OnlyRelevant_ForPlotting_Juvenile[is.na(SMR_OnlyRelevant_ForPlotting_Juvenile$p_SMR)==FALSE,]

SMR_OnlyRelevant_ForPlotting_Juvenile$variant_trait_forplotting<-paste(SMR_OnlyRelevant_ForPlotting_Juvenile$variant_id, SMR_OnlyRelevant_ForPlotting_Juvenile$trait, sep="_")

pdf("ManhattanPlot_SMR_Pval_Juvenile.pdf", width=16, height=5)
manhattan(SMR_OnlyRelevant_ForPlotting_Juvenile, chr="CHR", bp="LOC", snp="variant_id", p="p_SMR", genomewideline=FALSE, suggestiveline=FALSE)
dev.off()



#just SMR results for OF Distance Traveled marked off with a simple FDR<0.05 and FDR<0.10:

table(SMR$trait)
SMR_OnlyOFDistance_ForPlotting_Juvenile<-SMR[which(SMR$trait=="of_distance_traveled"),]
colnames(SMR_OnlyOFDistance_ForPlotting_Juvenile)
head(SMR_OnlyOFDistance_ForPlotting_Juvenile$variant_id)
TempSMR<-do.call(rbind, strsplit(SMR_OnlyOFDistance_ForPlotting_Juvenile$variant_id, split=":", fixed=TRUE))
head(TempSMR)
TempSMR2<-do.call(rbind, strsplit(TempSMR[,1], split="chr", fixed=TRUE))
head(TempSMR2)
SMR_OnlyOFDistance_ForPlotting_Juvenile$CHR<-as.numeric(TempSMR2[,2])
SMR_OnlyOFDistance_ForPlotting_Juvenile$LOC<-as.numeric(TempSMR[,2])
table(SMR_OnlyOFDistance_ForPlotting_Juvenile$CHR)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
# 777 371 432 317 342 291 350 346 270 540 210 205 149 200 236 210 158 191 156 186

rm(TempSMR, TempSMR2)
colnames(SMR_OnlyOFDistance_ForPlotting_Juvenile)

sum(SMR_OnlyOFDistance_ForPlotting_Juvenile$p_SMR=="inf")
#[1] NA
sum(is.na(SMR_OnlyOFDistance_ForPlotting_Juvenile$p_SMR))
#[1] 1

SMR_OnlyOFDistance_ForPlotting_Juvenile<-SMR_OnlyOFDistance_ForPlotting_Juvenile[is.na(SMR_OnlyOFDistance_ForPlotting_Juvenile$p_SMR)==FALSE,]

SMR_OnlyOFDistance_ForPlotting_Juvenile$variant_trait_forplotting<-paste(SMR_OnlyOFDistance_ForPlotting_Juvenile$variant_id, SMR_OnlyOFDistance_ForPlotting_Juvenile$trait, sep="_")

#There isn't an FDR<0.05, so I just put that line outside the plot 

pdf("ManhattanPlot_SMR_LocoScore_Pval_Juvenile_Colored_FDR05_FDR10.pdf", width=10, height=5)
manhattan(SMR_OnlyRelevant_ForPlotting_Juvenile, chr="CHR", bp="LOC", snp="variant_trait_forplotting", p="p_SMR", genomewideline=-log10(0.000008421762), suggestiveline=-log10(9.05e-05),  highlight=c(
  "chr1:166318198_of_distance_traveled",
  "chr1:166392350_of_distance_traveled",
  "chr1:141294154_of_distance_traveled", 
  "chr1:141037045_of_distance_traveled",
  "chr1:141312442_of_distance_traveled",
  "chr1:164561826_of_distance_traveled",
  "chr1:141313585_of_distance_traveled"))
dev.off()

#These genes don't survive cut-off anymore for the juveniles:
#AABR07071904.1, C2cd3, Mcee:
#"chr1:95003407_of_distance_traveled",
#"chr1:165794789_of_distance_traveled",
#"chr1:125425048_of_distance_traveled"
#"chr1:166153166_of_distance_traveled"
