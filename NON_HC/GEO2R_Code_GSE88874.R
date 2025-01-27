#Preparing the amygdala bHR/bLR data from GSE88874 for a meta-analysis
#Megan Hagenauer 
#12-02-2024

################################################################

#Code from GEO2R for differential expression analysis of GSE88874
#This is a differential expression analysis using just the adult bLR mother/bLRs and adult bHR mother/bHRs
# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0

#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE88874", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL19519", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX11111XXXXX",
               "XXXXX00000")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

exprs(gset) <- normalizeBetweenArrays(exprs(gset)) # normalize data

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Adult bHR","Adult bLR"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GB_ACC"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
# Please note that the code provided to generate graphs serves as a guidance to
# the users. It does not replicate the exact GEO2R web display due to multitude
# of graphical options.
# 
# The following will produce basic volcano plot using limma function:
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE88874", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE88874", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 5, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=5", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE88874")


#############

#Reading the GEO2R results into R to prepare them for the meta-analysis input.
#This code is adapted from Brain Data Alchemy Project (v.2023 and 2024)

GSE88874_DE<-read.delim("GSE88874.top.table.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE)

str(GSE88874_DE)
# 'data.frame':	22197 obs. of  7 variables:
# $ ID       : chr  "NM_001012119" "XM_001074233" "XM_001066392" "XM_001065805" ...
# $ adj.P.Val: num  4.8e-05 4.8e-05 4.8e-05 4.8e-05 4.8e-05 4.8e-05 4.8e-05 4.8e-05 4.8e-05 4.8e-05 ...
# $ P.Value  : num  6.00e-09 7.82e-09 1.02e-08 1.10e-08 1.38e-08 1.77e-08 2.16e-08 2.22e-08 2.50e-08 2.92e-08 ...
# $ t        : num  14.7 14.4 -14.1 -14 -13.7 ...
# $ B        : num  10.9 10.7 10.5 10.4 10.2 ...
# $ logFC    : num  2.26 2.13 -1.94 -1.92 -1.98 ...
# $ GB_ACC   : chr  "NM_001012119" "XM_001074233" "XM_001066392" "XM_001065805" ...

colnames(GSE88874_DE)[1]<-"ACCNUM"

table(table(GSE88874_DE[,1]))
# 1 
# 22197 
#each row has a separate accession number

#Adding annotation:

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Rn.eg.db")

library(org.Rn.eg.db)

UniKeys <- keys(org.Rn.eg.db, keytype="ACCNUM")

# temp <- select(org.Rn.eg.db,
#                keys = UniKeys,
#                columns = c("SYMBOL", "ENTREZID", "ENSEMBL"),
#                keytype = "ACCNUM")

temp2 <- select(org.Rn.eg.db,
               keys = UniKeys,
               columns = c("SYMBOL","ENTREZID"),
               keytype = "ACCNUM")


# check for duplicate genes
sum(duplicated(GSE88874_DE$ACCNUM))
#[1] 0

#Since there are none, we'll just proceed
annotateddata <- join(temp, GSE88874_DE, by="ACCNUM", type="inner")

nrow(annotateddata)
#[1] 8252
#That's certainly better than what was available on Gemma... I wonder if I should alert them to the bug.
sum(!duplicated(annotateddata$ACCNUM))
#[1] 8234
sum(!duplicated(annotateddata$ENTREZID))
#[1] 8117
sum(!duplicated(annotateddata$SYMBOL))
#[1] 8117
8252-8117
#[1] 135
#Only 135 duplicates (accession numbers that match to more than one Entrez ID) - not bad.
#But some of those are also entrez ids mapping to more than one ensembl id - we should probably add that annotation back in later.
 
write.csv(annotateddata, "GSE88874_DE_annotated.csv")

#Trying again:

annotateddata <- join(temp2, GSE88874_DE, by="ACCNUM", type="inner")

nrow(annotateddata)
#[1] 8234
sum(!duplicated(annotateddata$ACCNUM))
#[1] 8234
sum(!duplicated(annotateddata$ENTREZID))
#[1] 8117
8234-8117
#[1] 117 accession numbers that map to the same Entrez ID as another accession number

write.csv(annotateddata, "GSE88874_DE_annotated_noMultimapped.csv")


#########################

#Preparing the results to be fed into a meta-analysis:

#Making the column names match up with previous coding:

DE_Results<-annotateddata
colnames(DE_Results)
# [1] "ACCNUM"    "SYMBOL"    "ENTREZID"  "adj.P.Val" "P.Value"   "t"        
# [7] "B"         "logFC"     "GB_ACC" 

colnames(DE_Results)[3]<-"NCBIid"
colnames(DE_Results)[2]<-"GeneSymbol"


#Using the Brain Data Alchemy Function (v.2024):

FilteringDEResults_GoodAnnotation<-function(DE_Results){
  
  print("# of rows in results")
  print(nrow(DE_Results))
  
  print("# of rows with missing NCBI annotation:")
  print(sum(DE_Results$NCBIid==""|DE_Results$NCBIid=="null"))
  
  print("# of rows with NA NCBI annotation:")
  print(sum(is.na(DE_Results$NCBIid)))
  
  print("# of rows with missing Gene Symbol annotation:")
  print(sum(DE_Results$GeneSymbol==""|DE_Results$GeneSymbol=="null"))
  
  print("# of rows mapped to multiple NCBI_IDs:")
  print(length(grep('\\|', DE_Results$NCBIid)))
  
  print("# of rows mapped to multiple Gene Symbols:")
  print(length(grep('\\|', DE_Results$GeneSymbol)))
  
  #I only want the subset of data which contains rows that do not contain an NCBI EntrezID of ""
  DE_Results_NoNA<-DE_Results[(DE_Results$NCBIid==""|DE_Results$NCBIid=="null")==FALSE & is.na(DE_Results$NCBIid)==FALSE,]
  
  #I also only want the subset of data that is annotated with a single gene (not ambiguously mapped to more than one gene)
  if(length(grep('\\|', DE_Results_NoNA$NCBIid))==0){
    DE_Results_GoodAnnotation<<-DE_Results_NoNA
  }else{
    #I only want rows annotated with a single Gene Symbol (no pipe):
    DE_Results_GoodAnnotation<<-DE_Results_NoNA[-(grep('\\|', DE_Results_NoNA$NCBIid)),]
  }
  #I used a double arrow in that conditional to place DE_Results_GoodAnnotation back out in the environment outside the function 
  
  print("# of rows with good annotation")
  print(nrow(DE_Results_GoodAnnotation))
  
  #For record keeping (sometimes useful for troubleshooting later)
  write.csv(DE_Results_GoodAnnotation, "DE_Results_GoodAnnotation.csv")
  
  rm(DE_Results_NoNA, DE_Results)
  
  print("Outputted object: DE_Results_GoodAnnotation")
}

FilteringDEResults_GoodAnnotation(DE_Results)

# [1] "# of rows in results"
# [1] 8234
# [1] "# of rows with missing NCBI annotation:"
# [1] 0
# [1] "# of rows with NA NCBI annotation:"
# [1] 0
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] 0
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 0
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 0
# [1] "# of rows with good annotation"
# [1] 8234
# [1] "Outputted object: DE_Results_GoodAnnotation"


##########

colnames(DE_Results_GoodAnnotation)
# [1] "ACCNUM"     "GeneSymbol" "NCBIid"     "adj.P.Val"  "P.Value"    "t"    
# [7] "B"          "logFC"      "GB_ACC"  

NamesOfFoldChangeColumns<-c("logFC")
NamesOfTstatColumns<-c("t")
ComparisonsOfInterest<-c("bHR_vs_bLR")
#Note that the reference level appears reversed from our standard ref=bHR
#I double-checked by plotting NPY in GEO2R - NPY is strongly upregulated in bLRs in the amygdala
#We'll need to eventually fix that to run comparisons.

setwd("~/University of Michigan Dropbox/Megan Hagenauer/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/HR_LR_Amygdala_Meta/GSE88874/GEO2R")

library(tidyverse)
library(dplyr)

#The Brain Data Alchemy (v.2024) function CollapsingDEResults_OneResultPerGene seems to not be working anymore due to issues with the select() function
#I reverted back to the simpler code from v.2023

#This code will be dataset specific
#We need to extract the Log2FC ("Coef") and T-statistic ("t.") columns for the statistical contrasts relevant to our meta-analysis and place them into their own matrix

FoldChanges<-cbind(DE_Results_GoodAnnotation$logFC)

Tstats<-cbind(DE_Results_GoodAnnotation$t)

#Making the row names for the Log2FC and Tstat matrices the Entrez ID gene annotation:
row.names(FoldChanges)<-DE_Results_GoodAnnotation$NCBIid

row.names(Tstats)<-DE_Results_GoodAnnotation$NCBIid

#Let's rename our columns to something nicer describing the effect of interest:
#Note - we later discovered that this name needs to include the dataset identifier (GSEID#) for later joining and plotting purposes

ComparisonsOfInterest<-c("GSE88874_bHR_vs_bLR")

colnames(FoldChanges)<-ComparisonsOfInterest
colnames(Tstats)<-ComparisonsOfInterest

ExtractingDEResults<-function(GSE_ID, FoldChanges, Tstats){
  
  #We calculate the standard error by dividing the log2FC by the tstat
 StandardErrors<-FoldChanges/Tstats
  str(StandardErrors)
  
  #For running our meta-analysis, we are actually going to need the sampling variance instead of the standard error
  #The sampling variance is just the standard error squared.
  
  SamplingVars<-(StandardErrors)^2
  str(SamplingVars)
  
  TempMasterResults<-list(Log2FC=FoldChanges, Tstat=Tstats, SE=StandardErrors, SV=SamplingVars)
  
  assign(paste("DEResults", GSE_ID, sep="_"), TempMasterResults, envir = as.environment(1))
  
  print(paste("Output: Named DEResults", GSE_ID, sep="_"))
  
  rm(TempMasterResults, SamplingVars, StandardErrors, FoldChanges, Tstats)
  
}

ExtractingDEResults("GSE88874", FoldChanges, Tstats)
  
# num [1:8234, 1] 0.316 0.238 0.201 0.148 0.165 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:8234] "24153" "24157" "24158" "24159" ...
# ..$ : chr "GSE88874_bHR_vs_bLR"
# num [1:8234, 1] 0.0998 0.0568 0.0405 0.0219 0.0271 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:8234] "24153" "24157" "24158" "24159" ...
# ..$ : chr "GSE88874_bHR_vs_bLR"
# [1] "Output: Named DEResults_GSE88874"

#Flipping the direction of effect:

names(DEResults_GSE88874)
#[1] "Log2FC" "Tstat"  "SE"     "SV"  

str(DEResults_GSE88874[[1]])
# num [1:8234, 1] 0.1793 0.0468 -0.1668 0.7749 0.2533 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:8234] "24153" "24157" "24158" "24159" ...
# ..$ : chr "GSE88874_bHR_vs_bLR"

head(DEResults_GSE88874[[1]][,1])
# 24153      24157      24158      24159      24161      24162 
# 0.1792924  0.0467865 -0.1667984  0.7748948  0.2533114  0.6314374 

temp<-DEResults_GSE88874[[1]][,1]*(-1)

head(temp)
# 24153      24157      24158      24159      24161      24162 
# -0.1792924 -0.0467865  0.1667984 -0.7748948 -0.2533114 -0.6314374

DEResults_GSE88874[[1]][,1]<-temp

head(DEResults_GSE88874[[1]][,1])
# 24153      24157      24158      24159      24161      24162 
# -0.1792924 -0.0467865  0.1667984 -0.7748948 -0.2533114 -0.6314374 


head(DEResults_GSE88874[[2]][,1])
# 24153      24157      24158      24159      24161      24162 
# 0.5674450  0.1962427 -0.8287590  5.2397009  1.5374113  2.9588062 

temp<-DEResults_GSE88874[[2]][,1]*(-1)

head(temp)
# 24153      24157      24158      24159      24161      24162 
# -0.5674450 -0.1962427  0.8287590 -5.2397009 -1.5374113 -2.9588062 

DEResults_GSE88874[[2]][,1]<-temp

head(DEResults_GSE88874[[2]][,1])
# 24153      24157      24158      24159      24161      24162 
# -0.5674450 -0.1962427  0.8287590 -5.2397009 -1.5374113 -2.9588062 

#I guess I probably need to fix the columnnames too...

colnames(DEResults_GSE88874[[1]])
#[1] "GSE88874_bHR_vs_bLR"

#Yep
colnames(DEResults_GSE88874[[1]])<-"GSE88874_bLR_vs_bHR"
colnames(DEResults_GSE88874[[2]])<-"GSE88874_bLR_vs_bHR"
colnames(DEResults_GSE88874[[3]])<-"GSE88874_bLR_vs_bHR"
colnames(DEResults_GSE88874[[4]])<-"GSE88874_bLR_vs_bHR"


table(table(row.names(DEResults_GSE88874[[1]])))

# 1    2    3    4    5    9 
# 8024   78   12    1    1    1 

