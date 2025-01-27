#Preparing the bHR/bLR data from GSE86893 for a meta-analysis
#Megan Hagenauer 
#12-02-2024

#Amygdala
#Differential expression results downloaded from Gemma on 12-02-24
#For convenience, I converted the file to a .csv and removed the extra header information

setwd("~/University of Michigan Dropbox/Megan Hagenauer/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/HR_LR_Amygdala_Meta/GSE86893/13862_GSE86893_diffExpAnalysis_131475")

list.files()

GSE86893_DE<-read.csv("resultset_ID490846.data_Amygdala.csv", header=TRUE, stringsAsFactors = FALSE)

str(GSE86893_DE)
# 'data.frame':	20058 obs. of  7 variables:
# $ Element_Name           : int  259167 245961 64455 501688 364475 29705 266711 689499 680424 361493 ...
# $ Gene_Symbol            : chr  "Art5" "Dmgdh" "Rasd1" "Lsm11" ...
# $ Gene_Name              : chr  "ADP-ribosyltransferase 5" "dimethylglycine dehydrogenase" "ras related dexamethasone induced 1" "LSM11, U7 small nuclear RNA associated" ...
# $ NCBI_ID                : int  259167 245961 64455 501688 364475 29705 266711 689499 680424 361493 ...
# $ FoldChange_resistant.to: num  -1.743 -1.25 0.439 0.333 -0.227 ...
# $ Tstat_resistant.to     : num  -24.26 -10.8 10.09 7.65 -7.34 ...
# $ PValue_resistant.to    : num  3.23e-10 7.84e-07 1.46e-06 1.74e-05 2.48e-05 ...

#Looks like we'll need to reverse the direction of effect for all of these too.

##############

#Renaming things to recycle Brain Data Alchemy Project code:

DE_Results<-GSE86893_DE

colnames(DE_Results)
# [1] "Element_Name"            "Gene_Symbol"             "Gene_Name"              
# [4] "NCBI_ID"                 "FoldChange_resistant.to" "Tstat_resistant.to"     
# [7] "PValue_resistant.to"    

#Looking at Gemma's experimental design tab, apparently "resistant to" is there word for bLR and "sensitive to" is their word for bHR
#They have the reference set to "sensitive to" or bHR

colnames(DE_Results)[2]<-"GeneSymbol"
colnames(DE_Results)[4]<-"NCBIid"

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
# [1] 20058
# [1] "# of rows with missing NCBI annotation:"
# [1] NA
# [1] "# of rows with NA NCBI annotation:"
# [1] 2171
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] 2171
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 0
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 0
# [1] "# of rows with good annotation"
# [1] 17887
# [1] "Outputted object: DE_Results_GoodAnnotation"


colnames(DE_Results)

#We need to extract the Log2FC ("Coef") and T-statistic ("t.") columns for the statistical contrasts relevant to our meta-analysis and place them into their own matrix

FoldChanges<-cbind(DE_Results_GoodAnnotation$FoldChange_resistant.to)

Tstats<-cbind(DE_Results_GoodAnnotation$Tstat_resistant.to)

#Making the row names for the Log2FC and Tstat matrices the Entrez ID gene annotation:
row.names(FoldChanges)<-DE_Results_GoodAnnotation$NCBIid

row.names(Tstats)<-DE_Results_GoodAnnotation$NCBIid

#Let's rename our columns to something nicer describing the effect of interest:
#Note - we later discovered that this name needs to include the dataset identifier (GSEID#) for later joining and plotting purposes

ComparisonsOfInterest<-c("GSE86893_bLR_vs_bHR")

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

ExtractingDEResults("GSE86893", FoldChanges, Tstats)

# num [1:17887, 1] 0.0719 0.1158 0.0435 0.0435 0.0309 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:17887] "259167" "245961" "64455" "501688" ...
# ..$ : chr "GSE86893_bLR_vs_bHR"
# num [1:17887, 1] 0.005163 0.013399 0.001894 0.001894 0.000952 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:17887] "259167" "245961" "64455" "501688" ...
# ..$ : chr "GSE86893_bLR_vs_bHR"
# [1] "Output: Named DEResults_GSE86893"

table(table(names(DEResults_GSE86893[[1]][,1])))
# 1 
# 17887 
