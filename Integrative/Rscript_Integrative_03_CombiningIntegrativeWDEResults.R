#########################

#Combining the genetic results with the DE results:

colnames(F0_Meta_F2_DEResults)
str(F0_Meta_F2_DEResults)
#'data.frame':	13786 obs. of  42 variables:

max(table(F0_Meta_F2_DEResults$ENSEMBL.ID..Rnor6.Ensembl.v.103.))
#[1] 1
length(unique(F0_Meta_F2_DEResults$ENSEMBL.ID..Rnor6.Ensembl.v.103.))
#[1] 13786

str(F2_eqtls_SMR_GT)
#'data.frame':	5937 obs. of  82 variables:

F2_eqtls_SMR_GT_toJoin<-F2_eqtls_SMR_GT
colnames(F2_eqtls_SMR_GT_toJoin)[1]<-"ENSEMBL.ID..Rnor6.Ensembl.v.103." 

length(unique(F2_eqtls_SMR_GT_toJoin$ENSEMBL.ID..Rnor6.Ensembl.v.103.))
#[1] 5351
5937-5351
#[1] 586- Ensembl genes with multiple eVariants

str(F0_Meta_F2_DEResults)
#'data.frame':	13786 obs. of  42 variables:

colnames(F0_Meta_F2_DEResults)
length(unique(F0_Meta_F2_DEResults$ENSEMBL.ID..Rnor6.Ensembl.v.103.))
#[1] 13786
sum(is.na(F0_Meta_F2_DEResults$F0_Log2FC_bLRvsbHR)==FALSE)
#[1] 13786
sum(is.na(F0_Meta_F2_DEResults$F2_Log2FC_LocoScore)==FALSE)
#[1] 13339

#It looks like this data.frame doesn't include any of the F2 genes that weren't also found in the F0s... maybe?
#Although "F0_Meta_F2_4Models_withRnor6v88Coord_JustGenes.csv" has 14,148 rows, all of which have F0 data too.
#Ok, going back to the code, it looks like some of those rows were created when adding v88 annotation (for previous version of paper).

#O.k. Looking back at the F0 and F2 results in their original form:
#the F0 data included 13,786 Ensembl transcripts that survived filtering (which matches above)
#the F2 data included 14,056 Ensembl transcripts that survived filtering
#Our notes say later 13,339 ENSEMBL IDs present in both F0 and F2 
#I just double-checked this using the F0 and F2 results and confirmed it is true.
#So this is the full dataset for the F0s - with the F2 and bHR/bLR meta-analysis data joined to it. 

sum(F2_eqtls_SMR_GT$gene_id%in% F0_Meta_F2_DEResults$ENSEMBL.ID..Rnor6.Ensembl.v.103.)
#[1] 5747

F0_Meta_F2_DEResults_w_F2_eqtls_SMR<-join(F0_Meta_F2_DEResults, F2_eqtls_SMR_GT_toJoin, by="ENSEMBL.ID..Rnor6.Ensembl.v.103.", type="left", match="all")

str(F0_Meta_F2_DEResults_w_F2_eqtls_SMR)
#'data.frame':	14353 obs. of  123 variables:
#So fewer genes in the DE dataset have multiple eVariants - less than the eQTL results
14353-13788
#[1] 565
sum(is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref)==FALSE)
#[1] 5747

#write.csv(F0_Meta_F2_DEResults_w_F2_eqtls_SMR, "F0_Meta_F2_DEResults_w_F2_eqtls_SMR_GT.csv")
write.csv(F0_Meta_F2_DEResults_w_F2_eqtls_SMR, "F0_Meta_F2_DEResults_w_F2_eqtls_SMR_GT_wJuvenileSMR.csv")

####################


#Documenting and visualizing data subsetting by bHR/bLR DEGs and bHR/bLR segregation:

length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.))
#[1] 13786

#bHR/bLR DEGs:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1)]))
#[1] 1045

#bHR/bLR DEGs that are eGenes:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1 & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref)==FALSE)]))
#[1] 654

#bHR/bLR DEGs that are bHR/bLR segregated eGenes:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Included.in.F2.Candidate.Gene.Analysis.==1 & is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref)==FALSE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)]))
#[1] 480

#Segregated eGenes in the DE dataset:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref)==FALSE & F0_Meta_F2_DEResults_w_F2_eqtls_SMR$Gprimest>0.27)]))
#[1] 2395

#eGenes in the DE dataset:
length(unique(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(is.na(F0_Meta_F2_DEResults_w_F2_eqtls_SMR$log2_aFC_alt_ref)==FALSE)]))
#[1] 5180

#It might help to have some sort of visualization for this for the paper, since we're narrowing our scope again:
library(VennDiagram)
grid.newpage() 
draw.pairwise.venn(2395,1045,480,col=c("darkgrey", "darkgrey"), fill=c("grey", "grey"))

