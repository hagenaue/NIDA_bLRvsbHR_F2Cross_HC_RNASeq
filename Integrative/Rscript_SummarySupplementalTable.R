#Quickly making a supplemental table summarizing all of our main results from our current NIDA_U01 output
#2023-02-02
#Megan Hagenauer

setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/DEResults_F0F2_SummaryTables")

F0_Meta_F2_results<-read.csv("F0_Meta_F2_results_ForSupplTable.csv", header=TRUE, stringsAsFactors = FALSE)

str(F0_Meta_F2_results)
colnames(F0_Meta_F2_results)

F0_Meta_F2_results_toJoin<-F0_Meta_F2_results[,c(-21, -25, -29, -33, -37)]

dim(F0_Meta_F2_results_toJoin)
#[1] 14148    41

length(unique(F0_Meta_F2_results_toJoin$ENSEMBL))
#[1] 13786

table(table(F0_Meta_F2_results_toJoin$ENSEMBL))
# 1     2     3     4     5     6     7     8    11 
# 13507   236    26    10     1     2     1     2     1

dim(unique(F0_Meta_F2_results_toJoin[,-1]))
#[1] 14148    40
#grrr... I wonder how much of this was just created after we added the v88 annotation.

colnames(F0_Meta_F2_results_toJoin)

dim(unique(F0_Meta_F2_results_toJoin[,c(-1,-c(33:41))]))
#[1] 13786    31

#Yep. Alright - the next question is whether they were run through the analysis that way previously...
#Yes - but with the notable point that this is the full dataset for the F0s - with the F2 and bHR/bLR meta-analysis data joined to it. 
#The F2 data actually included 14,056 transcripts that survived filtering

#Our notes say later 13,339 ENSEMBL IDs present in both F0 and F2 - maybe that is what we should include in the supplement?

UniqueENSEMBL<-data.frame("ENSEMBL"=(unique(F0_Meta_F2_results_toJoin$ENSEMBL)), stringsAsFactors=FALSE)

library(plyr)

F0_Meta_F2_results_toJoin_UniqueENSEMBL<-join(x=UniqueENSEMBL, y=F0_Meta_F2_results_toJoin, by="ENSEMBL", type="left", match="first")
dim(F0_Meta_F2_results_toJoin_UniqueENSEMBL)
#[1] 13786    41



setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/reworktogethertuesday142022/reworktogethertuesday142022")

F0_Meta_F2_LineageGenes<-read.csv("F0_Meta_F2_4Models_LineageGenes_Updated.csv", header=TRUE, stringsAsFactors = FALSE)

str(F0_Meta_F2_LineageGenes)
colnames(F0_Meta_F2_LineageGenes)

F0_Meta_F2_LineageGenes_toJoin<-F0_Meta_F2_LineageGenes[,c(3,227,228,229,232,233,234)]

setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/DEResults_F0F2_SummaryTables")

#double-checking
sum(F0_Meta_F2_LineageGenes_toJoin$ENSEMBL%in%F0_Meta_F2_results_toJoin_UniqueENSEMBL$ENSEMBL)
#[1] 1063

F0_Meta_F2_results_wLineageGenes_ForSuppl<-join(F0_Meta_F2_results_toJoin_UniqueENSEMBL, F0_Meta_F2_LineageGenes_toJoin, by="ENSEMBL")

setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/DEResults_F0F2_SummaryTables")

write.csv(F0_Meta_F2_results_wLineageGenes_ForSuppl, "F0_Meta_F2_results_wLineageGenes_ForSuppl.csv")

#Double-checking a few things for the legend:

sum(F0_Meta_F2_results_wLineageGenes_ForSuppl$GENENAME==F0_Meta_F2_results_wLineageGenes_ForSuppl$SYMBOL.1)
#[1] 13786

sum(F0_Meta_F2_results_wLineageGenes_ForSuppl$SYMBOL%in%F0_Meta_F2_results_wLineageGenes_ForSuppl$SYMBOL.1)
#[1] 12275
sum(is.na(F0_Meta_F2_results_wLineageGenes_ForSuppl$SYMBOL))
#[1] 789
12275+789
#13064
sum(duplicated(F0_Meta_F2_results_wLineageGenes_ForSuppl$SYMBOL))
#[1] 818
sum(duplicated(F0_Meta_F2_results_wLineageGenes_ForSuppl$SYMBOL)&is.na(F0_Meta_F2_results_wLineageGenes_ForSuppl$SYMBOL)==FALSE)
#[1] 30
F0_Meta_F2_results_wLineageGenes_ForSuppl$SYMBOL[(duplicated(F0_Meta_F2_results_wLineageGenes_ForSuppl$SYMBOL)&is.na(F0_Meta_F2_results_wLineageGenes_ForSuppl$SYMBOL)==FALSE)]

#So approximately 700 of the meta-analysis gene symbols aren't in the F0/F2 dataframe

F0_Meta_F2_results_wLineageGenes_ForSuppl$SYMBOL[(F0_Meta_F2_results_wLineageGenes_ForSuppl$SYMBOL!=F0_Meta_F2_results_wLineageGenes_ForSuppl$SYMBOL.1) & is.na(F0_Meta_F2_results_wLineageGenes_ForSuppl$SYMBOL)==FALSE]

#Ah - it's equivalent to whatever the annotation was that was used for F2:
sum(F0_Meta_F2_results_wLineageGenes_ForSuppl$SYMBOL%in%F0_Meta_F2_results_wLineageGenes_ForSuppl$Symbol)
#[1] 13489
#Yep - that, plus some NAs.


#I should probably look up what that annotation was
#Rnor6 Ensemblv103
