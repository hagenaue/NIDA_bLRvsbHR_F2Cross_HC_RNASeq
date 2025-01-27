#Code for running a (very simple) bLR vs. bHR meta-analysis of public Non-HPC data
#GSE88874 (Amygdala), GSE86893 (Amygdala, Dorsal Raphe), John Stead's F4's (Cortex, Hypothalamus)
#Megan Hagenauer
#12-09-2024

#This code is adapted from the 2024 version of the Brain Data Alchemy Project

library(plyr)

install.packages("metafor")

library(metafor)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("multtest")

library(multtest)

##############

DEResults_GSE86893_AMY<-DEResults_GSE86893

colnames(DEResults_GSE86893_AMY[[1]])<-"GSE86893_AMY_bLR_vs_bHR"
colnames(DEResults_GSE86893_AMY[[2]])<-"GSE86893_AMY_bLR_vs_bHR"
colnames(DEResults_GSE86893_AMY[[3]])<-"GSE86893_AMY_bLR_vs_bHR"
colnames(DEResults_GSE86893_AMY[[4]])<-"GSE86893_AMY_bLR_vs_bHR"

#I had to come back and fix this one - the direction of effect was reversed:
str(DEResults_GSE86893_AMY)
head(DEResults_GSE86893_AMY[[1]])
# GSE86893_AMY_bLR_vs_bHR
# 259167                  1.7431
# 245961                  1.2498
# 64455                  -0.4393
# 501688                 -0.3329
# 364475                  0.2265
# 29705                   0.2528
temp<-DEResults_GSE86893_AMY[[1]]*-1
head(temp)
# GSE86893_AMY_bLR_vs_bHR
# 259167                 -1.7431
# 245961                 -1.2498
# 64455                   0.4393
# 501688                  0.3329
# 364475                 -0.2265
# 29705                  -0.2528
DEResults_GSE86893_AMY[[1]]<-temp
head(DEResults_GSE86893_AMY[[1]])
# GSE86893_AMY_bLR_vs_bHR
# 259167                 -1.7431
# 245961                 -1.2498
# 64455                   0.4393
# 501688                  0.3329
# 364475                 -0.2265
# 29705                  -0.2528

head(DEResults_GSE86893_AMY[[2]])
# GSE86893_AMY_bLR_vs_bHR
# 259167                 24.2593
# 245961                 10.7969
# 64455                 -10.0943
# 501688                 -7.6499
# 364475                  7.3408
# 29705                   6.8921
temp<-DEResults_GSE86893_AMY[[2]]*-1
head(temp)
# GSE86893_AMY_bLR_vs_bHR
# 259167                -24.2593
# 245961                -10.7969
# 64455                  10.0943
# 501688                  7.6499
# 364475                 -7.3408
# 29705                  -6.8921

DEResults_GSE86893_AMY[[2]]<-temp
head(DEResults_GSE86893_AMY[[2]])
# GSE86893_AMY_bLR_vs_bHR
# 259167                -24.2593
# 245961                -10.7969
# 64455                  10.0943
# 501688                  7.6499
# 364475                 -7.3408
# 29705                  -6.8921

colnames(DEResults_GSE88874_AMY[[1]])<-"GSE88874_AMY_bLR_vs_bHR"
colnames(DEResults_GSE88874_AMY[[2]])<-"GSE88874_AMY_bLR_vs_bHR"
colnames(DEResults_GSE88874_AMY[[3]])<-"GSE88874_AMY_bLR_vs_bHR"
colnames(DEResults_GSE88874_AMY[[4]])<-"GSE88874_AMY_bLR_vs_bHR"

##############

#Aligning datasets:

#A function for aligning all of our rat differential expression results from different datasets into a single data frame for Log2FCs and sampling variances (SVs):

ListOfRatDEResults<-list(DEResults_F4_CTX, DEResults_F4_HYP, DEResults_GSE86893_Raphe, DEResults_GSE86893_AMY, DEResults_GSE88874_AMY)

AligningRatDatasets<-function(ListOfRatDEResults){
  
  #Making an empty list to hold our results:
  Rat_MetaAnalysis_FoldChange_Dfs<-list()
  
  #Looping over all of the rat differential expression results:
  for(i in c(1:length(ListOfRatDEResults))){
    
    #Placing each of the log2FC results for each dataset into a single list
    #Each element in the list is formatted so that the rownames are Rat Entrez Gene ID and then there are columns containing the Log2FC for the differential expression results:
    Rat_MetaAnalysis_FoldChange_Dfs[[i]]<-data.frame(Rat_EntrezGene.ID=row.names(ListOfRatDEResults[[i]][[1]]),ListOfRatDEResults[[i]][[1]], stringsAsFactors=FALSE)
  }
  
  #Letting the user know the structure of the set of Log2FC dataframes that we are starting out with:
  print("Rat_MetaAnalysis_FoldChange_Dfs:")
  print(str(Rat_MetaAnalysis_FoldChange_Dfs))
  
  #Running "join all" on the list to align all of the results by Entrez Gene ID and make them a single data frame:
  Rat_MetaAnalysis_FoldChanges<<-join_all(Rat_MetaAnalysis_FoldChange_Dfs, by="Rat_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  #Letting the user know the structure of the dataframe that we just created:
  print("Rat_MetaAnalysis_FoldChanges:")
  print(str(Rat_MetaAnalysis_FoldChanges))
  
  #Doing the same steps for the sampling variances:
  
  Rat_MetaAnalysis_SV_Dfs<-list()
  
  for(i in c(1:length(ListOfRatDEResults))){
    Rat_MetaAnalysis_SV_Dfs[[i]]<-data.frame(Rat_EntrezGene.ID=row.names(ListOfRatDEResults[[i]][[4]]),ListOfRatDEResults[[i]][[4]], stringsAsFactors=FALSE)
  }
  
  print("Rat_MetaAnalysis_SV_Dfs:")
  print(str(Rat_MetaAnalysis_SV_Dfs))
  
  Rat_MetaAnalysis_SV<<-join_all(Rat_MetaAnalysis_SV_Dfs, by="Rat_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  
  print("Rat_MetaAnalysis_SV:")
  print(str(Rat_MetaAnalysis_SV))
  
  #Cleaning up our environment to remove unneeded objects:
  rm(Rat_MetaAnalysis_SV_Dfs, Rat_MetaAnalysis_FoldChange_Dfs)
}

AligningRatDatasets(ListOfRatDEResults)

# # [1] "Rat_MetaAnalysis_FoldChange_Dfs:"
# List of 5
# $ :'data.frame':	10069 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID: chr [1:10069] "100036582" "100036765" "100125364" "100125370" ...
# ..$ F4_CTX_bLR_vs_bHR: num [1:10069] 0.01111 -0.01187 -0.01497 -0.01287 -0.00556 ...
# $ :'data.frame':	10069 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID: chr [1:10069] "100036582" "100036765" "100125364" "100125370" ...
# ..$ F4_HYP_bLR_vs_bHR: num [1:10069] 0.0588 0.0519 0.0535 0.0207 -0.0502 ...
# $ :'data.frame':	18068 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID        : chr [1:18068] "108352578" "108350052" "103692946" "171409" ...
# ..$ GSE86893_Raphe_bLR_vs_bHR: num [1:18068] -2.04 1.96 -4.52 4.07 -4.41 ...
# $ :'data.frame':	17887 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID      : chr [1:17887] "259167" "245961" "64455" "501688" ...
# ..$ GSE86893_AMY_bLR_vs_bHR: num [1:17887] -1.743 -1.25 0.439 0.333 -0.227 ...
# $ :'data.frame':	8252 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID      : chr [1:8252] "24153" "24157" "24158" "24159" ...
# ..$ GSE88874_AMY_bLR_vs_bHR: num [1:8252] -0.1793 -0.0468 0.1668 -0.7749 -0.2533 ...
# NULL
# [1] "Rat_MetaAnalysis_FoldChanges:"
# 'data.frame':	20173 obs. of  6 variables:
#   $ Rat_EntrezGene.ID        : chr  "100036582" "100036765" "100125364" "100125370" ...
# $ F4_CTX_bLR_vs_bHR        : num  0.01111 -0.01187 -0.01497 -0.01287 -0.00556 ...
# $ F4_HYP_bLR_vs_bHR        : num  0.0588 0.0519 0.0535 0.0207 -0.0502 ...
# $ GSE86893_Raphe_bLR_vs_bHR: num  NA -0.0129 -0.0812 -0.0041 NA ...
# $ GSE86893_AMY_bLR_vs_bHR  : num  NA -0.0284 -0.0578 0.1051 NA ...
# $ GSE88874_AMY_bLR_vs_bHR  : num  NA -0.417 NA NA NA ...
# NULL
# [1] "Rat_MetaAnalysis_SV_Dfs:"
# List of 5
# $ :'data.frame':	10069 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID: chr [1:10069] "100036582" "100036765" "100125364" "100125370" ...
# ..$ F4_CTX_bLR_vs_bHR: num [1:10069] 0.00349 0.00257 0.00525 0.00461 0.00227 ...
# $ :'data.frame':	10069 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID: chr [1:10069] "100036582" "100036765" "100125364" "100125370" ...
# ..$ F4_HYP_bLR_vs_bHR: num [1:10069] 0.00385 0.00283 0.00245 0.01332 0.00211 ...
# $ :'data.frame':	18068 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID        : chr [1:18068] "108352578" "108350052" "103692946" "171409" ...
# ..$ GSE86893_Raphe_bLR_vs_bHR: num [1:18068] 0.0238 0.0261 0.1794 0.2143 0.2784 ...
# $ :'data.frame':	17887 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID      : chr [1:17887] "259167" "245961" "64455" "501688" ...
# ..$ GSE86893_AMY_bLR_vs_bHR: num [1:17887] 0.005163 0.013399 0.001894 0.001894 0.000952 ...
# $ :'data.frame':	8252 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID      : chr [1:8252] "24153" "24157" "24158" "24159" ...
# ..$ GSE88874_AMY_bLR_vs_bHR: num [1:8252] 0.0998 0.0568 0.0405 0.0219 0.0271 ...
# NULL
# [1] "Rat_MetaAnalysis_SV:"
# 'data.frame':	20173 obs. of  6 variables:
#   $ Rat_EntrezGene.ID        : chr  "100036582" "100036765" "100125364" "100125370" ...
# $ F4_CTX_bLR_vs_bHR        : num  0.00349 0.00257 0.00525 0.00461 0.00227 ...
# $ F4_HYP_bLR_vs_bHR        : num  0.00385 0.00283 0.00245 0.01332 0.00211 ...
# $ GSE86893_Raphe_bLR_vs_bHR: num  NA 0.00545 0.00639 0.00812 NA ...
# $ GSE86893_AMY_bLR_vs_bHR  : num  NA 0.00738 0.00267 0.00412 NA ...
# $ GSE88874_AMY_bLR_vs_bHR  : num  NA 0.0838 NA NA NA ...
# NULL


##############

#I tweaked the next bit of code a bit, because we only have rat datasets (no mouse datasets) in the meta-analysis, so no need for orthologs
#But we probably do want some extra gene annotation eventually.
#And the column names & numbers aren't going to match the meta-analysis code... so I'll need to tweak that next too...

MetaAnalysis_FoldChanges<-Rat_MetaAnalysis_FoldChanges
  
MetaAnalysis_SV<-Rat_MetaAnalysis_SV

##############

#Should we see if there is any correlation between the results in the studies?
#Probably not - these are pretty dinky little datasets

FoldChange_CorMatrix<-cor(as.matrix(MetaAnalysis_FoldChanges[,-1]), method="spearman", use="pairwise.complete.obs")

write.csv(FoldChange_CorMatrix, "FoldChange_CorMatrix.csv")

#Yep, nada

##############

#Adapting the meta-analysis function to fit our purposes...
#I needed to change all of the references to the annotation columns to fit the fact that we currently only have rat Entrez ID

str(MetaAnalysis_FoldChanges)
# 'data.frame':	20173 obs. of  6 variables:
#   $ Rat_EntrezGene.ID        : chr  "100036582" "100036765" "100125364" "100125370" ...
# $ F4_CTX_bLR_vs_bHR        : num  0.01111 -0.01187 -0.01497 -0.01287 -0.00556 ...
# $ F4_HYP_bLR_vs_bHR        : num  0.0588 0.0519 0.0535 0.0207 -0.0502 ...
# $ GSE86893_Raphe_bLR_vs_bHR: num  NA -0.0129 -0.0812 -0.0041 NA ...
# $ GSE86893_AMY_bLR_vs_bHR  : num  NA -0.0284 -0.0578 0.1051 NA ...
# $ GSE88874_AMY_bLR_vs_bHR  : num  NA -0.417 NA NA NA ...

RunBasicMetaAnalysis<-function(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV){
  
  #The function first provides information about how many of the statistical contrasts have NA values as their differential expression results for each gene:
  MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges[,-1], 1, function(y) sum(is.na(y)))
  
  print("Table of # of NAs per Row (Gene):")
  print(table(MetaAnalysis_FoldChanges_NAsPerRow))
  
  #Then any row (gene) that has too many NAs is removed from the analysis:
  MetaAnalysis_FoldChanges_ForMeta<<-MetaAnalysis_FoldChanges[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  MetaAnalysis_SV_ForMeta<<-MetaAnalysis_SV[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  
  print("MetaAnalysis_FoldChanges_ForMeta:")
  print(str(MetaAnalysis_FoldChanges_ForMeta))
  
  #I'm going to make an empty matrix to store the results of my meta-analysis:
  metaOutput<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), 6)
  
  #And then run a loop that run's a meta-analysis on the differential expression results (i.e., the columns that aren't annotation) for each gene (row):
  for(i in c(1:nrow(MetaAnalysis_FoldChanges_ForMeta))){
    
    #When pulling out the log2FC values and sampling variances (SV) for each gene, we use the function as.numeric to make sure they are in numeric matrix format because this is the required input format for the meta-analysis function that we will use:
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-1])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-1])
    
    #I added a function tryCatch that double-checks that the meta-analysis function (rma) doesn't produce errors (which breaks the loop):
    skip_to_next <- FALSE
    tryCatch(TempMeta<-rma(effect, var), error = function(e) {skip_to_next <<- TRUE})
    
    #If everything looks good, we move on to running the meta-analysis using a model that treats the variation in Log2FC across studies as random effects:
    if(skip_to_next){}else{
      TempMeta<-rma(effect, var)
      metaOutput[i, 1]<-TempMeta$b #gives estimate Log2FC
      metaOutput[i, 2]<-TempMeta$se #gives standard error
      metaOutput[i, 3]<-TempMeta$pval #gives pval
      metaOutput[i, 4]<-TempMeta$ci.lb #gives confidence interval lower bound
      metaOutput[i, 5]<-TempMeta$ci.ub #gives confidence interval upper bound
      metaOutput[i, 6]<-NumberOfComparisons-sum(is.na(effect))#Number of comparisons with data
      rm(TempMeta)
    }
    rm(effect, var)
  }
  
  #Naming the columns in our output:
  colnames(metaOutput)<-c("Log2FC_estimate", "SE", "pval", "CI_lb", "CI_ub", "Number_Of_Comparisons")
  
  #The row names for our output are the rat entrez ids: 
  row.names(metaOutput)<-MetaAnalysis_FoldChanges_ForMeta[,1]
  
  #We return this output back into our global environment
  metaOutput<<-metaOutput
  MetaAnalysis_Annotation<<-MetaAnalysis_FoldChanges_ForMeta[,1]
  return(metaOutput)
  return(MetaAnalysis_Annotation)
  
  #... and provide the user with an update about the newly created object:
  
  print("metaOutput:")
  print(str(metaOutput))
  
  print("Top of metaOutput:")
  print(head(metaOutput))
  
  print("Bottom of metaOutput")
  print(tail(metaOutput))
  
}

NumberOfComparisons<-5
CutOffForNAs<-3

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)

# [1] "Table of # of NAs per Row (Gene):"
# MetaAnalysis_FoldChanges_NAsPerRow
# 0    1    2    3    4 
# 5251 4111 2147 7033 1631 
# [1] "MetaAnalysis_FoldChanges_ForMeta:"
# 'data.frame':	11509 obs. of  6 variables:
#   $ Rat_EntrezGene.ID        : chr  "100036765" "100125364" "100125370" "100125373" ...
# $ F4_CTX_bLR_vs_bHR        : num  -0.0119 -0.015 -0.0129 -0.0251 0.0494 ...
# $ F4_HYP_bLR_vs_bHR        : num  0.0519 0.0535 0.0207 -0.0215 -0.0862 ...
# $ GSE86893_Raphe_bLR_vs_bHR: num  -0.0129 -0.0812 -0.0041 -0.1385 0.0094 ...
# $ GSE86893_AMY_bLR_vs_bHR  : num  -0.0284 -0.0578 0.1051 -0.1363 0.1393 ...
# $ GSE88874_AMY_bLR_vs_bHR  : num  -0.417 NA NA NA NA ...
# NULL
# There were 50 or more warnings (use warnings() to see the first 50)


#######

#FDR correction and adding annotation:

#I'm going to need to adapt this Brain Data Alchemy function too...
#I'm going to grab the annotation that we already extracted from org.Rn.eg.db 
#EntrezVsGeneSymbol
str(EntrezVsGeneSymbol)
# 'data.frame':	47813 obs. of  2 variables:
# $ EntrezGeneID: chr  "100034253" "100036765" "100049583" "100101342" ...
# $ GeneSymbol  : chr  "Gnl3l" "Ccdc92" "Trex1" "Vom2r-ps11" ...

colnames(EntrezVsGeneSymbol)[1]<-"ENTREZID"

str(MetaAnalysis_Annotation)
#chr [1:11509] "100036765" "100125364" "100125370" "100125373" ...

temp<-data.frame(ENTREZID=MetaAnalysis_Annotation)
str(temp)
# 'data.frame':	11509 obs. of  1 variable:
#   $ ENTREZID: chr  "100036765" "100125364" "100125370" "100125373" ...
MetaAnalysis_Annotation<-temp

colnames(metaOutput)
#[1] "Log2FC_estimate"       "SE"                   
# [3] "pval"                  "CI_lb"                
# [5] "CI_ub"                 "Number_Of_Comparisons"

FalseDiscoveryCorrection<-function(metaOutput, EntrezVsGeneSymbol, MetaAnalysis_Annotation){
  
  #This calculates the false discovery rate, or q-value, for each of our p-values using the Benjamini-Hochberg procedure:
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,3], proc=c("BH"))
  
  #Then we put those results back into the order of our orginal output:
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  #And bind the false discovery rate (FDR) to the rest of the meta-analysis output:
  metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])
  
  #And name that column FDR:
  colnames(metaOutputFDR)[7]<-"FDR"
  
  #These results are returned to our global environment:
  metaOutputFDR<<-metaOutputFDR
  
  #We let the user know the basic structure of the meta-analysis output with FDR added to it (just to make sure everything still looks good)
  print("metaOutputFDR:")
  print(str(metaOutputFDR))
  
  #Then we make a dataframe that adds the annotation to that output:
  TempDF<-cbind.data.frame(metaOutputFDR, MetaAnalysis_Annotation)
  #And then adds even more detailed gene annotation:
  
  #First the detailed annotation for the rat genes:
  TempDF3<-join(TempDF, EntrezVsGeneSymbol, by="ENTREZID", type="left", match="first")
  
  #This is renamed and returned to our global environment:
  metaOutputFDR_annotated<-TempDF3
  metaOutputFDR_annotated<<-metaOutputFDR_annotated
  
  #And written out into our working directory:
  write.csv(metaOutputFDR_annotated, "metaOutputFDR_annotated.csv")
  
  #Then we make a version of the output in order by p-value:
  metaOutputFDR_OrderbyPval<<-metaOutputFDR_annotated[order(metaOutputFDR_annotated[,4]),]
  
  #Let's write out a version of the output in order by p-value:
  write.csv(metaOutputFDR_OrderbyPval, "metaOutputFDR_orderedByPval.csv")
  
  #And give the user some information about their results:
  
  print("Do we have any genes that are statistically significant following loose false discovery rate correction (FDR<0.10)?")
  print(sum(metaOutputFDR_annotated[,8]<0.10, na.rm=TRUE))
  
  print("Do we have any genes that are statistically significant following traditional false discovery rate correction (FDR<0.05)?")
  print(sum(metaOutputFDR_annotated[,8]<0.05, na.rm=TRUE))
  
  print("What are the top results?")
  print(head(metaOutputFDR_annotated[order(metaOutputFDR_annotated[,4]),]))
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
} 

FalseDiscoveryCorrection(metaOutput, EntrezVsGeneSymbol, MetaAnalysis_Annotation)

# [1] "metaOutputFDR:"
# num [1:11509, 1:7] 0.00234 -0.01583 0.03614 -0.08597 0.01119 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:11509] "100036765" "100125364" "100125370" "100125373" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
# NULL
# [1] "Do we have any genes that are statistically significant following loose false discovery rate correction (FDR<0.10)?"
# [1] 36
# [1] "Do we have any genes that are statistically significant following traditional false discovery rate correction (FDR<0.05)?"
# [1] 20
# [1] "What are the top results?"
# ENTREZID Log2FC_estimate         SE         pval      CI_lb
# X500282   500282     -0.08860318 0.01485010 2.423904e-09 -0.1177088
# X58966     58966      0.26177757 0.04475902 4.957297e-09  0.1740515
# X314637   314637      0.25493207 0.04557914 2.229549e-08  0.1655986
# X291978   291978      0.20320426 0.04298672 2.277037e-06  0.1189518
# X24231     24231      0.34438017 0.07440458 3.683577e-06  0.1985499
# X367153   367153     -0.22229081 0.04856485 4.712520e-06 -0.3174762
# CI_ub Number_Of_Comparisons          FDR GeneSymbol
# X500282 -0.05949752                     5 2.789671e-05      Arl8b
# X58966   0.34950364                     5 2.852676e-05      Ramp2
# X314637  0.34426554                     3 8.553295e-05    Slc39a3
# X291978  0.28745668                     4 6.551606e-03       Dus2
# X24231   0.49021046                     3 8.478857e-03         C2
# X367153 -0.12710545                     3 9.039398e-03      Cep70

################

#Making a few forest plots:

#I'll probably need to adapt the Brain Data Alchemy Code for this too...

MakeForestPlots<-function(metaOutputFDR_annotated, EntrezIDAsCharacter, species){
  
  #I originally wrote this function using only mouse Entrez IDs as input 
  #but then I realized that the function didn't work for genes that were only found in rats
  #so now the function allows either rat or mouse Entrez ids as input, and includes a conditional (if/else) statement
  
  if(species=="Mouse"){

    
  }else if(species=="Rat"){
    
    #This set of code does all of the same processes as above, but interpreting the EntrezID as a Rat Entrez ID:
    
    RatGeneSymbol<-metaOutputFDR_annotated$GeneSymbol[which(metaOutputFDR_annotated$ENTREZID==EntrezIDAsCharacter)][1]
    
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[which(MetaAnalysis_FoldChanges_ForMeta$Rat_EntrezGene.ID==EntrezIDAsCharacter),-1][1,])
    
    var<-as.numeric(MetaAnalysis_SV_ForMeta[which(MetaAnalysis_FoldChanges_ForMeta$Rat_EntrezGene.ID==EntrezIDAsCharacter),-1][1,]) 
    
  }else{
    
    print("Please use either 'Mouse' or 'Rat' to indicate whether you are using annotation for mouse or rat genes")
    
  }
  
  #This code makes the Forest Plot
  
  #First it opens up a .pdf file to output the plot into:
  #It automatically names that file with the mouse and rat ene symbols
  pdf(paste("ForestPlot_", "Rat", RatGeneSymbol, EntrezIDAsCharacter,".pdf", sep="_"), height=5, width=8)
  
  #This code makes the forest plot:
  #Note that the x-axis limits are currently set to -3 to 3
  #This may be too big or too small for visualizing the results for some genes.
  forest.rma(rma(effect, var), slab=colnames(MetaAnalysis_FoldChanges_ForMeta)[-1],  xlim=c(-3, 3))
  
  #This code labels the forest plot with the mouse and rat gene symbols:
  mtext(paste("Rat", EntrezIDAsCharacter, RatGeneSymbol, sep="_"), line=-1.5, cex=2)
  
  #This closes the connection to the .pdf file, finishing the plot
  dev.off()
  
}


setwd("~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/HRLR_NONHPC/Forestplots")

#Making forest plots for some of the top genes:

#Maybe I should just loop it

for(i in c(1:length(metaOutputFDR_annotated$ENTREZID[which(metaOutputFDR_annotated$FDR<0.10)]))){
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter=metaOutputFDR_OrderbyPval$ENTREZID[i], species="Rat")
}
# #Were clearly still having issues with data from genes with more than one NM
# 
# for(i in c(1:length(metaOutputFDR_annotated$ENTREZID[which(metaOutputFDR_annotated$FDR<0.10)]))){
#    print(
#   metaOutputFDR_annotated$ENTREZID[which(metaOutputFDR_annotated$ENTREZID==metaOutputFDR_OrderbyPval$ENTREZID[i])])
# }


setwd("~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/HRLR_NONHPC/ForestPlotsHPC")

#Making forest plots for some of the top genes from the Hippocampus that were also implicated by genetics:

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="54315", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="25277", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="307833", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="64355", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="85249", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="24479", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="353231", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="361239", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="361239", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="64471", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="308759", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="361968", species="Rat")
#Tmem144 apparently isn't in the meta-analysis output 

#Gimap5
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="246774", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="499189", species="Rat")
#Wdr93 not found



#Genes suggested previously by genetic analyses:

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="83610", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="308543", species="Rat")

#Double-checking genes that were mentioned in the original papers:

#Itgb3bp - listed as down-regulated in the bLR amygdala in Cohen 2017 (GSE86893)
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="362548", species="Rat")
#... and it is (now fixed) and down-regulated in the bLR amygdala in the forest plot for that dataset too.

#Cd74 is listed as upregulated in the bLR raphe in Cohen 2017 (GSE86893)
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="25599", species="Rat")
#upregulated in the bLR raphe in the forest plot

#Trhr
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="25570", species="Rat")
#doesn't seem to show much

#Bmp4
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="25296", species="Rat")
#doesn't seem to show much

#Cacng3
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="140724", species="Rat")
#doesn't seem to show much

#Cxcr4
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="60628", species="Rat")
#up in HRs in the cortex (down in bLRs)... barely... matching John's report.

#Sstr3 (Smstr28)
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="171044", species="Rat")
#Ok, that one matches John's earlier report - up in HR or down in bLR in the Cortex

#Slc38a2
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="29642", species="Rat")
#That's not significant like it was in John's earlier report, but still up in HR or down in LR in the HYP in both

#Slc12a2
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="83629", species="Rat")
#Yep, still up in HR or down in LR in the HYP in John's earlier report and now

#The Cohen 2015 does not discuss bHR/bLR differential expression, so we can't pull out results from our analysis for comparison

###################


#Interesting: Many of the top genes from the HPC replicate in one dataset (GSE88874) but not the other (GSE86893) - I wonder if the GSE88874 dataset included more amygdala tissue that is HPC-like (e.g., basolateral) vs. striatal-like (CEA, MEA).

#Cohen 2015 (GSE88874) cites this paper for their amygdala dissection, but it really doesn't say anything besides the fact that they used hole punches:
#https://www.sciencedirect.com/science/article/pii/S0006899313011086?via=ihub

#Cohen 2017 (GSE86893) doesn't describe their dissection at all.

#But it is worth noting that GSE88874 was collected as part of an experiment that also include HPC and GSE86893 was collected as part of an experiment that included dorsal raphe, so I could imagine the sections used for punching being more rostral for GSE88874 and more caudal for GSE86893. 

#################

#Comparison to bHR/bLR hippocampal results:

setwd("~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/HRLR_NONHPC")

HPC_DEResults<-read.csv("TableS3_F0_Meta_F2_results_ForSuppl_20240223_DEResults.csv", header=TRUE, stringsAsFactors = FALSE)

str(HPC_DEResults)
#'data.frame':	13788 obs. of  37 variables:

#for joining to amygdala meta-analysis results:
colnames(HPC_DEResults)[2]
#[1] "ENSEMBL.ID..Rnor6.Ensembl.v.103."

colnames(metaOutputFDR_annotated)

UniKeys2 <- keys(org.Rn.eg.db, keytype="ENTREZID")

MoreAnnotation <- select(org.Rn.eg.db,
               keys = UniKeys2,
              columns = c("ENSEMBL"),
              keytype = c("ENTREZID"))

dim(metaOutputFDR_annotated)
#[1] 11509     9

metaOutputFDR_moreannotated<-join(metaOutputFDR_annotated, MoreAnnotation, by="ENTREZID", type="left")

dim(metaOutputFDR_moreannotated)
#[1] 11551    10
#So there are a few ENSEMBL ids that are represented multiple times (because they map to more than one ENTREZ id).

colnames(HPC_DEResults)[2]<-"ENSEMBL"

HPC_DEResults_vs_nonHPC_DEResults<-join(HPC_DEResults, metaOutputFDR_moreannotated, by="ENSEMBL", type="left", match="all")

str(HPC_DEResults_vs_nonHPC_DEResults)
#'data.frame':	13972 obs. of  46 variables:
13972-13788
#[1] 184
#So there are a few ENSEMBL ids that are represented multiple times 

write.csv(HPC_DEResults_vs_nonHPC_DEResults, "HPC_DEResults_vs_NonHPC_DEResults.csv")

pdf("Scatterplot_F0_vs_nonHPC.pdf", height=5, width=4)
plot(HPC_DEResults_vs_nonHPC_DEResults$F0_Log2FC_bLRvsbHR~HPC_DEResults_vs_nonHPC_DEResults$Log2FC_estimate, ylim=c(-4,4), xlab="nonHPC Log2FC", ylab="F0 HPC Log2FC")
dev.off()
#That's a pretty wussy positive correlation, but a few genes stand out.

cor(HPC_DEResults_vs_nonHPC_DEResults$F0_Log2FC_bLRvsbHR, HPC_DEResults_vs_nonHPC_DEResults$Log2FC_estimate, method="spearman", use="pairwise.complete.obs")
#[1] 0.1494256

#Looks like we used RRHOs ranked by tstat for the paper - I should probably do something parallel here:

head((HPC_DEResults_vs_nonHPC_DEResults$Log2FC_estimate/HPC_DEResults_vs_nonHPC_DEResults$SE))

HPC_DEResults_vs_nonHPC_DEResults$nonHPC_Meta_Tstat<-(HPC_DEResults_vs_nonHPC_DEResults$Log2FC_estimate/HPC_DEResults_vs_nonHPC_DEResults$SE)

pdf("Scatterplot_F0_vs_nonHPC_tstat.pdf", height=5, width=4)
plot(HPC_DEResults_vs_nonHPC_DEResults$F0_Tstat_bLRvsbHR~HPC_DEResults_vs_nonHPC_DEResults$nonHPC_Meta_Tstat, xlab="nonHPC Tstat", ylab="F0 HPC Tstat")
dev.off()
#Not horrible - a slight positive correlation is actually visible

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RRHO")

library(RRHO)

TempDF<-HPC_DEResults_vs_nonHPC_DEResults[is.na(HPC_DEResults_vs_nonHPC_DEResults$F0_Tstat_bLRvsbHR)==FALSE & is.na(HPC_DEResults_vs_nonHPC_DEResults$nonHPC_Meta_Tstat)==FALSE, ]

str(TempDF)
#'data.frame':	9322 obs. of  47 variables:

table(table(TempDF$ENSEMBL))
# 1    2    3    4    5    6    9 
# 8989  132   13    1    1    2    1 
#RRHO is unhappy that there is redundancy/multimapped genes

TempDF$Row<-c(1:nrow(TempDF))

table(table(TempDF$Row))

list1<-data.frame(ENSEMBL=TempDF$Row, Metric=TempDF$F0_Tstat_bLRvsbHR)

list2<-data.frame(ENSEMBL=TempDF$Row, Metric=TempDF$nonHPC_Meta_Tstat)

setwd("~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/HRLR_NONHPC/rrho")

RRHO(list1, list2, labels=c("F0 HPC: bLR vs. bHR", "Non-HPC: bLR vs. bHR"), plots=TRUE, alternative="two.sided", outputdir="~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/HRLR_NONHPC/rrho",  BY=TRUE, log10.ind=TRUE)

cor.test(HPC_DEResults_vs_nonHPC_DEResults$F0_Tstat_bLRvsbHR, HPC_DEResults_vs_nonHPC_DEResults$nonHPC_Meta_Tstat, method="spearman", use="pairwise.complete.obs")
#	Spearman's rank correlation rho
# data:  HPC_DEResults_vs_nonHPC_DEResults$F0_Tstat_bLRvsbHR and HPC_DEResults_vs_nonHPC_DEResults$nonHPC_Meta_Tstat
# S = 1.1399e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.1556875

cor.test(HPC_DEResults_vs_nonHPC_DEResults$F0_Log2FC_bLRvsbHR, HPC_DEResults_vs_nonHPC_DEResults$Log2FC_estimate, method="spearman", use="pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  HPC_DEResults_vs_nonHPC_DEResults$F0_Log2FC_bLRvsbHR and HPC_DEResults_vs_nonHPC_DEResults$Log2FC_estimate
# S = 1.1484e+11, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1494256 

table(TempDF$bHR_bLR_DEG==TRUE)
# FALSE  TRUE 
# 8561   761

list1<-data.frame(ENSEMBL=TempDF$Row[TempDF$bHR_bLR_DEG==TRUE], Metric=TempDF$F0_Tstat_bLRvsbHR[TempDF$bHR_bLR_DEG==TRUE])

list2<-data.frame(ENSEMBL=TempDF$Row[TempDF$bHR_bLR_DEG==TRUE], Metric=TempDF$nonHPC_Meta_Tstat[TempDF$bHR_bLR_DEG==TRUE])

setwd("~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/HRLR_NONHPC/rrhoHPCDEG")

RRHO(list1, list2, labels=c("F0 HPC: bLR vs. bHR", "Non-HPC: bLR vs. bHR"), plots=TRUE, alternative="two.sided", outputdir="~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/HRLR_NONHPC/rrhoHPCDEG",  BY=TRUE, log10.ind=TRUE)


cor.test(HPC_DEResults_vs_nonHPC_DEResults$F0_Tstat_bLRvsbHR[HPC_DEResults_vs_nonHPC_DEResults$bHR_bLR_DEG==TRUE], HPC_DEResults_vs_nonHPC_DEResults$nonHPC_Meta_Tstat[HPC_DEResults_vs_nonHPC_DEResults$bHR_bLR_DEG==TRUE], method="spearman", use="pairwise.complete.obs")
#	Spearman's rank correlation rho
# data:  HPC_DEResults_vs_nonHPC_DEResults$F0_Tstat_bLRvsbHR[HPC_DEResults_vs_nonHPC_DEResults$bHR_bLR_DEG == TRUE] and HPC_DEResults_vs_nonHPC_DEResults$nonHPC_Meta_Tstat[HPC_DEResults_vs_nonHPC_DEResults$bHR_bLR_DEG == TRUE]
# S = 49294352, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.3288877 

cor.test(HPC_DEResults_vs_nonHPC_DEResults$F0_Log2FC_bLRvsbHR[HPC_DEResults_vs_nonHPC_DEResults$bHR_bLR_DEG==TRUE], HPC_DEResults_vs_nonHPC_DEResults$Log2FC_estimate[HPC_DEResults_vs_nonHPC_DEResults$bHR_bLR_DEG==TRUE], method="spearman", use="pairwise.complete.obs")
# Spearman's rank correlation rho
# 
# data:  HPC_DEResults_vs_nonHPC_DEResults$F0_Log2FC_bLRvsbHR[HPC_DEResults_vs_nonHPC_DEResults$bHR_bLR_DEG == TRUE] and HPC_DEResults_vs_nonHPC_DEResults$Log2FC_estimate[HPC_DEResults_vs_nonHPC_DEResults$bHR_bLR_DEG == TRUE]
# S = 51674543, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2964829

HPC_DEResults_vs_nonHPC_DEResults$GENENAME_F2..Rnor6.Ensembl.v.103.[which( HPC_DEResults_vs_nonHPC_DEResults$F0_Log2FC_bLRvsbHR>0.5 &  HPC_DEResults_vs_nonHPC_DEResults$Log2FC_estimate>0.5)]
# [1] "Tnnt1"        "Cym"          "Ptgds"        "Fhad1"       
# [5] "LOC100362027" "LOC100362027" "Ghdc"         "Grifin"      
# [9] NA  

HPC_DEResults_vs_nonHPC_DEResults$GENENAME_F2..Rnor6.Ensembl.v.103.[which( HPC_DEResults_vs_nonHPC_DEResults$F0_Log2FC_bLRvsbHR<(-0.5) &  HPC_DEResults_vs_nonHPC_DEResults$Log2FC_estimate<(-0.5))]
#[1] "Fam111a" "Prss12"  NA        "Cyp4f5" 

#Interesting.

pdf("Scatterplot_HPC_Meta_vs_nonHPC.pdf", height=5, width=4)
plot(HPC_DEResults_vs_nonHPC_DEResults$MetaAnalysis_estimatedD_bLRvsbHR~HPC_DEResults_vs_nonHPC_DEResults$Log2FC_estimate, xlab="nonHPC Log2FC", ylab="Meta HPC d")
dev.off()
#Again, basically a blob with a few exceptions

pdf("Scatterplot_HPC_Meta_vs_nonHPC_tstat.pdf", height=5, width=4)
plot(HPC_DEResults_vs_nonHPC_DEResults$MetaAnalysis_estimatedD_bLRvsbHR~HPC_DEResults_vs_nonHPC_DEResults$nonHPC_Meta_Tstat, xlab="nonHPC tstat", ylab="Meta HPC d")
dev.off()

cor.test(HPC_DEResults_vs_nonHPC_DEResults$MetaAnalysis_estimatedD_bLRvsbHR, HPC_DEResults_vs_nonHPC_DEResults$Log2FC_estimate, method="spearman", use="pairwise.complete.obs")
# Spearman's rank correlation rho
# 
# data:  HPC_DEResults_vs_nonHPC_DEResults$MetaAnalysis_estimatedD_bLRvsbHR and HPC_DEResults_vs_nonHPC_DEResults$Log2FC_estimate
# S = 9.1315e+10, p-value = 3.278e-09
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.06462958

HPC_DEResults_vs_nonHPC_DEResults$GENENAME_F2..Rnor6.Ensembl.v.103.[which( HPC_DEResults_vs_nonHPC_DEResults$MetaAnalysis_estimatedD_bLRvsbHR>2 &  HPC_DEResults_vs_nonHPC_DEResults$Log2FC_estimate>0.5)]
#[1] "Tnnt1"  "Ghdc"   "Rpl17"  "Grifin"

HPC_DEResults_vs_nonHPC_DEResults$GENENAME_F2..Rnor6.Ensembl.v.103.[which( HPC_DEResults_vs_nonHPC_DEResults$MetaAnalysis_estimatedD_bLRvsbHR<(-2) &  HPC_DEResults_vs_nonHPC_DEResults$Log2FC_estimate<(-0.5))]
#[1] "Prss12"

HPC_DEResults_vs_nonHPC_DEResults$GENENAME_F2..Rnor6.Ensembl.v.103.[which((HPC_DEResults_vs_nonHPC_DEResults$F0_FDR_bLRvsbHR<0.1 | HPC_DEResults_vs_nonHPC_DEResults$MetaAnalysis_FDR_bLRvsbHR<0.1) & HPC_DEResults_vs_nonHPC_DEResults$pval<0.05)]
# [1] "Stxbp5"    "Vip"       "Syt5"      "Fkrp"      "Calm3"    
# [6] "Kif22"     "Vkorc1"    "Sipa1"     "Fam111a"   "Aldh1a1"  
# [11] "Ldb1"      "Golph3l"   "Cym"       "Cnn3"      "Prss12"   
# [16] "Trmt10a"   "Chpf2"     "Rarres2"   "Tmem176a"  "Abcg2"    
# [21] "Ctnna2"    "Card9"     "Tor1b"     "Hat1"      "Spi1"     
# [26] "Chgb"      "Acss2"     "Mmp24"     "C1qc"      "C1qa"     
# [31] "Fhad1"     "Stag2"     "Rbmx"      "Iah1"      "Nudt4"    
# [36] "Pdxp"      "Ilf3"      "Abhd14a"   "Gls"       "Cab39"    
# [41] "Fcgr3a"    "Rgs7"      "Capn2"     "Aida"      "Rpusd1"   
# [46] "Adam19"    "G3bp1"     "Trappc1"   "Wsb1"      "Wnk4"     
# [51] "Arf2"      "Cd300lf"   "Fbxo34"    "Nudt18"    "Cdyl"     
# [56] "LOC302192" "Dzip3"     "Setd6"     "Man2b1"    "Cbfb"     
# [61] "Dus2"      "Zfp90"     "RT1-T24-4" "Rtn4ip1"   "Psmg3"    
# [66] "Eif2b1"    "Ift81" 

HPC_DEResults_vs_nonHPC_DEResults$GENENAME_F2..Rnor6.Ensembl.v.103.[which((HPC_DEResults_vs_nonHPC_DEResults$bHR_bLR_DEG==TRUE) & HPC_DEResults_vs_nonHPC_DEResults$pval<0.05)]
# [1] "Stxbp5"    "Vip"       "Syt5"      "Fkrp"      "Calm3"    
# [6] "Kif22"     "Vkorc1"    "Sipa1"     "Fam111a"   "Aldh1a1"  
# [11] "Ldb1"      "Golph3l"   "Cym"       "Cnn3"      "Prss12"   
# [16] "Trmt10a"   "Chpf2"     "Rarres2"   "Tmem176a"  "Abcg2"    
# [21] "Ctnna2"    "Card9"     "Tor1b"     "Rc3h2"     "Hat1"     
# [26] "Spi1"      "Chgb"      "Acss2"     "Mmp24"     "C1qc"     
# [31] "C1qa"      "Fhad1"     "Stag2"     "Rbmx"      "Iah1"     
# [36] "Snx13"     "Esyt2"     "Tmem196"   "Nudt4"     "Pdxp"     
# [41] "Ilf3"      "Abhd14a"   "Gls"       "Cab39"     "Fcgr3a"   
# [46] "Rgs7"      "Capn2"     "Aida"      "Rpusd1"    "Adam19"   
# [51] "G3bp1"     "Trappc1"   "Wsb1"      "Naglu"     "Wnk4"     
# [56] "Arf2"      "Cd300lf"   "Fbxo34"    "Nudt18"    "Idnk"     
# [61] "Cdyl"      "LOC302192" "Setd4"     "Dzip3"     "Setd6"    
# [66] "Man2b1"    "Cbfb"      "Dus2"      "Zfp90"     "RT1-T24-4"
# [71] "Rtn4ip1"   "Psmg3"     "Eif2b1"    "P2rx4"     "Ift81"

#More strict version:
HPC_DEResults_vs_nonHPC_DEResults$GENENAME_F2..Rnor6.Ensembl.v.103.[which((HPC_DEResults_vs_nonHPC_DEResults$F0_FDR_bLRvsbHR<0.1 | HPC_DEResults_vs_nonHPC_DEResults$MetaAnalysis_FDR_bLRvsbHR<0.1) & HPC_DEResults_vs_nonHPC_DEResults$FDR<0.1)]
#[1] "Syt5"    "Fam111a" "Aldh1a1" "Hat1"    "Stag2"   "Wnk4"   
#[7] "Dus2"

HPC_DEResults_vs_nonHPC_DEResults$GENENAME_F2..Rnor6.Ensembl.v.103.[which((HPC_DEResults_vs_nonHPC_DEResults$bHR_bLR_DEG==TRUE) & HPC_DEResults_vs_nonHPC_DEResults$FDR<0.1)]
#[1] "Syt5"    "Fam111a" "Aldh1a1" "Hat1"    "Stag2"   "Wnk4"   
#[7] "Dus2" 

#"C1qa" 
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="298566", species="Rat")

#"C1qc" 362634
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="362634", species="Rat")
#more convincing

#Ghdc
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="303542", species="Rat")

#Grifin
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="117130", species="Rat")
#That one needs a bigger xlim

#Zfp90
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="498945", species="Rat")
#Driven by one dataset

#Rarres2 
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="297073", species="Rat")
#reasonably convincing

#Gls 24398
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="24398", species="Rat")
#mostly driven by results from one study
#Syt5
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="54309", species="Rat")

#Fam111a
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="499322", species="Rat")
#That one needs a larger xlim!

#"Aldh1a1"
#24188
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="24188", species="Rat")
#Pretty

# "Hat1"
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="296501", species="Rat")
#Driven by one study

#"Stag2"
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="313304", species="Rat")
#reasonably compelling

#"Wnk4"
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="287715", species="Rat")
#driven by one dataset

#"Dus2"
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="291978", species="Rat")
#reasonably compelling

save.image("~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/HRLR_NONHPC/ForestPlotsHPC/bHR_bLR_NonHPC_Meta.RData")

sessionInfo()

# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.1.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/Detroit
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets 
# [7] methods   base     
# 
# other attached packages:
#   [1] multtest_2.62.0              metafor_4.6-0               
# [3] numDeriv_2016.8-1.1          metadat_1.2-0               
# [5] Matrix_1.7-1                 limma_3.62.1                
# [7] rae230arnentrezgprobe_25.0.0 rae230arnentrezgcdf_25.0.0  
# [9] rae230arnentrezg.db_25.0.0   plyr_1.8.9                  
# [11] org.Rn.eg.db_3.20.0          AnnotationDbi_1.68.0        
# [13] IRanges_2.40.0               S4Vectors_0.44.0            
# [15] affy_1.84.0                  Biobase_2.66.0              
# [17] BiocGenerics_0.52.0          BiocManager_1.30.25         
# 
# loaded via a namespace (and not attached):
#   [1] RSQLite_2.3.9           lattice_0.22-6         
# [3] grid_4.4.2              fastmap_1.2.0          
# [5] blob_1.2.4              jsonlite_1.8.9         
# [7] GenomeInfoDb_1.42.1     DBI_1.2.3              
# [9] survival_3.7-0          httr_1.4.7             
# [11] UCSC.utils_1.2.0        preprocessCore_1.68.0  
# [13] Biostrings_2.74.0       cli_3.6.3              
# [15] rlang_1.1.4             crayon_1.5.3           
# [17] XVector_0.46.0          bit64_4.5.2            
# [19] splines_4.4.2           cachem_1.1.0           
# [21] tools_4.4.2             memoise_2.0.1          
# [23] mathjaxr_1.6-0          GenomeInfoDbData_1.2.13
# [25] vctrs_0.6.5             R6_2.5.1               
# [27] png_0.1-8               zlibbioc_1.52.0        
# [29] KEGGREST_1.46.0         bit_4.5.0.1            
# [31] MASS_7.3-61             pkgconfig_2.0.3        
# [33] affyio_1.76.0           Rcpp_1.0.13-1          
# [35] statmod_1.5.0           rstudioapi_0.17.1      
# [37] nlme_3.1-166            compiler_4.4.2

str(HPC_DEResults)
colnames(HPC_DEResults)

colnames(GSE88874_DE)
head(GSE88874_DE)
str(annotateddata)
colnames(annotateddata)
#[1] "ACCNUM"    "SYMBOL"    "ENTREZID"  "ENSEMBL"   "adj.P.Val"
#[6] "P.Value"   "t"         "B"         "logFC"     "GB_ACC"  

annotateddata[annotateddata$SYMBOL=="Ucp2",]
# ACCNUM SYMBOL ENTREZID            ENSEMBL       adj.P.Val  P.Value
# 1502 NM_019354   Ucp2    54315 ENSRNOG00000017854 0.0023207 0.000157
#         t        B    logFC    GB_ACC
# 1502 5.45556 1.017057 1.047481 NM_019354

annotateddata[annotateddata$SYMBOL=="Mfge8",]
# ACCNUM SYMBOL ENTREZID            ENSEMBL adj.P.Val
# 582 NM_001040186  Mfge8    25277 ENSRNOG00000017510 0.0013944
# 583    NM_012811  Mfge8    25277 ENSRNOG00000017510 0.0015340
# P.Value        t        B     logFC       GB_ACC
# 582 7.39e-05 5.938772 1.780857 1.0193712 NM_001040186
# 583 8.49e-05 5.847984 1.639947 0.9711744    NM_012811

#Yep, apparently those results weren't typos - they actually are huge log2fc

#Double checking direction of effect one more time:
annotateddata[annotateddata$SYMBOL=="Npy",]
# ACCNUM SYMBOL ENTREZID            ENSEMBL adj.P.Val  P.Value
# 224 NM_012614    Npy    24604 ENSRNOG00000046449 0.0033939 0.000274
#         t         B      logFC    GB_ACC
# 224 -5.108545 0.4478018 -0.6882057 NM_012614

DEResults_GSE88874[[1]][names(DEResults_GSE88874[[1]][,1])=="24604",1]
#[1] 0.6882057
#Good, this matches the direction of effect on GEO

#Another double-check:
#In the McCoy paper they report lower levels of Grin2B in bLRs
DEResults_GSE88874[[1]][names(DEResults_GSE88874[[1]][,1])=="24410",1]
#[1] -2.086542
#yep, looks good.

hist(annotateddata$logFC)
#It looks like the log2fc and tstats in this dataset in general are huge. Lots of genes with +/-1 log2FC
hist(annotateddata$t)
#And huge t-stats too for a dinky little dataset lol
#lots over abs(5)

sum(is.na(HPC_DEResults$ENSEMBL))
#[1] 0
sum(is.na(annotateddata$ENSEMBL))
#[1] 218

GSE88874_DE_Annotated<-annotateddata[is.na(annotateddata$ENSEMBL)==FALSE,]
str(GSE88874_DE_Annotated)
#'data.frame':	8034 obs. of  10 variables:
sum(is.na(GSE88874_DE_Annotated$ENSEMBL))
#[1] 0

HPC_DEResults_vs_GSE88874<-join(HPC_DEResults, GSE88874_DE_Annotated, by="ENSEMBL", type="left")

cor.test(HPC_DEResults_vs_GSE88874$F0_Tstat_bLRvsbHR, HPC_DEResults_vs_GSE88874$t, method="spearman", use="pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  HPC_DEResults_vs_GSE88874$F0_Tstat_bLRvsbHR and HPC_DEResults_vs_GSE88874$t
# S = 3.7402e+10, p-value = 0.0125
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.03221027

cor.test(HPC_DEResults_vs_GSE88874$F0_Log2FC_bLRvsbHR, HPC_DEResults_vs_GSE88874$logFC, method="spearman", use="pairwise.complete.obs")
#	Spearman's rank correlation rho
# data:  HPC_DEResults_vs_GSE88874$F0_Log2FC_bLRvsbHR and HPC_DEResults_vs_GSE88874$logFC
# S = 3.8191e+10, p-value = 2.796e-05
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.05399892

plot(HPC_DEResults_vs_GSE88874$F0_Tstat_bLRvsbHR~ HPC_DEResults_vs_GSE88874$t)

plot(HPC_DEResults_vs_GSE88874$F0_Log2FC_bLRvsbHR~ HPC_DEResults_vs_GSE88874$logFC)

#Interesting - it is definitely a negative correlation, even though many of our top genes are super in the same direction.



colnames(GSE86893_DE)
head(GSE86893_DE)

#Itgb3bp - listed as down-regulated in the bLR amygdala in Cohen 2017 (GSE86893)
GSE86893_DE[GSE86893_DE$Gene_Symbol=="Itgb3bp",]
# Element_Name Gene_Symbol                               Gene_Name
# 2255       362548     Itgb3bp integrin subunit beta 3 binding protein
# NCBI_ID FoldChange_resistant.to Tstat_resistant.to
# 2255  362548                 -0.2421            -1.9605
# PValue_resistant.to
# 2255              0.0784

#Looks good.

# colnames(HPC_DEResults)
# HPC_DEResults$
# 
#   MoreAnnotation <- select(org.Rn.eg.db,
#                            keys = UniKeys2,
#                            columns = c("ENSEMBL"),
#                            keytype = c("ENTREZID"))
#   
#   
# DEResults_GSE86893_AMY[[1]],DEResults_GSE86893_AMY[[2]], DEResults_GSE86893_AMY[[3]], 
# 
# head(DEResults_GSE86893_AMY[[2]])
# DEResults_GSE88874_AMY


#Making an RRHO plot for just the AMY meta-analysis:

setwd("~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/HR_LR_Amygdala_Meta")

list.files()

HPC_DEResults_vs_AMYMeta<-read.csv("HPC_DEResults_vs_AMY_DEResults.csv", header=TRUE, stringsAsFactors = FALSE)

str(HPC_DEResults_vs_AMYMeta)

HPC_DEResults_vs_AMYMeta$AMY_Meta_Tstat<-HPC_DEResults_vs_AMYMeta$Log2FC_estimate/HPC_DEResults_vs_AMYMeta$SE
  
library(RRHO)

TempDF2<-HPC_DEResults_vs_AMYMeta[is.na(HPC_DEResults_vs_AMYMeta$F0_Tstat_bLRvsbHR)==FALSE & is.na(HPC_DEResults_vs_AMYMeta$AMY_Meta_Tstat)==FALSE, ]

str(TempDF2)
#'data.frame':	5968 obs. of  48 variables:

table(table(TempDF2$ENSEMBL))
# 1    2    3    4    5    9 
# 5730   95   10    1    1    1 
#RRHO is unhappy that there is redundancy/multimapped genes

TempDF2$Row<-c(1:nrow(TempDF2))

table(table(TempDF2$Row))

list1<-data.frame(ENSEMBL=TempDF2$Row, Metric=TempDF2$F0_Tstat_bLRvsbHR)

list2<-data.frame(ENSEMBL=TempDF2$Row, Metric=TempDF2$AMY_Meta_Tstat)

setwd("~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/HR_LR_Amygdala_Meta/rrho")

RRHO(list1, list2, labels=c("F0 HPC: bLR vs. bHR", "AMY: bLR vs. bHR"), plots=TRUE, alternative="two.sided", outputdir="~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/HR_LR_Amygdala_Meta/rrho",  BY=TRUE, log10.ind=TRUE)

cor.test(TempDF2$F0_Tstat_bLRvsbHR, TempDF2$AMY_Meta_Tstat, method="spearman", use="pairwise.complete.obs")
#	Spearman's rank correlation rho
# data:  TempDF2$F0_Tstat_bLRvsbHR and TempDF2$AMY_Meta_Tstat
# S = 3.0935e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.1267907

cor.test(TempDF2$F0_Log2FC_bLRvsbHR, TempDF2$Log2FC_estimate, method="spearman", use="pairwise.complete.obs")
#	Spearman's rank correlation rho
# data:  TempDF2$F0_Log2FC_bLRvsbHR and TempDF2$Log2FC_estimate
# S = 3.2077e+10, p-value = 2.459e-13
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.09457585

setwd("~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/HR_LR_Amygdala_Meta/rrhoHPCDEG")

length(TempDF2$Row[TempDF2$bHR_bLR_DEG==TRUE])

list1<-data.frame(ENSEMBL=TempDF2$Row[TempDF2$bHR_bLR_DEG==TRUE], Metric=TempDF2$F0_Tstat_bLRvsbHR[TempDF2$bHR_bLR_DEG==TRUE])

list2<-data.frame(ENSEMBL=TempDF2$Row[TempDF2$bHR_bLR_DEG==TRUE], Metric=TempDF2$AMY_Meta_Tstat[TempDF2$bHR_bLR_DEG==TRUE])

RRHO(list1, list2, labels=c("F0 HPC: bLR vs. bHR", "AMY: bLR vs. bHR"), plots=TRUE, alternative="two.sided", outputdir="~/University of Michigan Dropbox/Megan Hagenauer/Laptop/Animal/HRLR/NIDA_U01/HR_LR_Amygdala_Meta/rrhoHPCDEG", BY=TRUE, log10.ind=TRUE)

cor.test(TempDF2$F0_Tstat_bLRvsbHR[TempDF2$bHR_bLR_DEG==TRUE], TempDF2$AMY_Meta_Tstat[TempDF2$bHR_bLR_DEG==TRUE], method="spearman", use="pairwise.complete.obs")
#	Spearman's rank correlation rho
# data:  TempDF2$F0_Tstat_bLRvsbHR[TempDF2$bHR_bLR_DEG == TRUE] and TempDF2$AMY_Meta_Tstat[TempDF2$bHR_bLR_DEG == TRUE]
# S = 18926483, p-value = 1.713e-10
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.270751

cor.test(TempDF2$F0_Log2FC_bLRvsbHR[TempDF2$bHR_bLR_DEG==TRUE], TempDF2$Log2FC_estimate[TempDF2$bHR_bLR_DEG==TRUE], method="spearman", use="pairwise.complete.obs")
# Spearman's rank correlation rho
# data:  TempDF2$F0_Log2FC_bLRvsbHR[TempDF2$bHR_bLR_DEG == TRUE] and TempDF2$Log2FC_estimate[TempDF2$bHR_bLR_DEG == TRUE]
# S = 21380688, p-value = 3.968e-05
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#      rho 
# 0.176189 

TempDF2$GENENAME..Rnor6.Ensembl.v.88.[TempDF2$FDR<0.1]
# [1] "Syt5"         "Apoc1"        "Lgals7"       "Prmt3"       
# [5] "Ubfd1"        "Asrgl1"       "Kazald1"      "Ablim1"      
# [9] "Ndufs6"       "Cym"          "Cyr61"        "Cav2"        
# [13] "Cav1"         "Egr4"         "Ccdc174"      "Surf4"       
# [17] "Plpp7"        "Chrna1"       "Hdc"          "Wwp1"        
# [21] "Gabrr2"       "Mtap"         "Tekt2"        "A3galt2"     
# [25] "Pla2g2c"      "Efnb1"        "Plp1"         "RGD1307315"  
# [29] "Tmem196"      "Csrp2"        "LOC100359574" "Nt5e"        
# [33] "RGD1310507"   "Nradd"        "Entpd3"       "Rcan2"       
# [37] "Il1rl1"       "Yes1"         "Bmp3"         "Pf4"         
# [41] "Pds5a"        "Hba-a2"       "Rpl30"        "Aurkb"       
# [45] "Dnajc7"       "Tcam1"        "Rnase4"       "Nefm"        
# [49] "Elf1"         "Gadd45g"      "Pip4k2a"      "Eif4ebp1"    
# [53] "Cd47"         "Dgcr2"        "Grp"          "Rrad"        
# [57] "Farsa"        "Inpp4b"       "Ggnbp1"  

TempDF2$GENENAME..Rnor6.Ensembl.v.88.[TempDF2$FDR<0.05]
# [1] "Apoc1"        "Lgals7"       "Prmt3"        "Ubfd1"       
# [5] "Asrgl1"       "Kazald1"      "Ablim1"       "Ndufs6"      
# [9] "Cyr61"        "Cav1"         "Egr4"         "Ccdc174"     
# [13] "Plpp7"        "Chrna1"       "Hdc"          "Wwp1"        
# [17] "Gabrr2"       "Tekt2"        "A3galt2"      "Efnb1"       
# [21] "RGD1307315"   "Tmem196"      "Csrp2"        "LOC100359574"
# [25] "Nt5e"         "RGD1310507"   "Nradd"        "Entpd3"      
# [29] "Rcan2"        "Il1rl1"       "Yes1"         "Bmp3"        
# [33] "Pf4"          "Pds5a"        "Rpl30"        "Dnajc7"      
# [37] "Tcam1"        "Rnase4"       "Nefm"         "Elf1"        
# [41] "Gadd45g"      "Pip4k2a"      "Eif4ebp1"     "Cd47"        
# [45] "Grp"          "Farsa"        "Inpp4b"  

TempDF2$GENENAME..Rnor6.Ensembl.v.88.[TempDF2$bHR_bLR_DEG & TempDF2$FDR<0.1]
#[1] "Syt5"    "Cym"     "Cav2"    "Cav1"    "Efnb1"   "Tmem196"

TempDF2$GENENAME..Rnor6.Ensembl.v.88.[TempDF2$bHR_bLR_DEG & TempDF2$FDR<0.05]
#[1] "Cav1"    "Efnb1"   "Tmem196"

TempDF2$GENENAME..Rnor6.Ensembl.v.88.[TempDF2$bHR_bLR_DEG & TempDF2$pval<0.05]
# [1] "Syt5"       "Fkrp"       "Hsd3b7"     "Itgad"      "Gal"       
# [6] "Chka"       "Aldh1a1"    "Crhbp"      "Cartpt"     "Golph3l"   
# [11] "Cym"        "Cxxc4"      "Trmt10a"    "Cav2"       "Cav1"      
# [16] "Akr1b1"     "Tmem176a"   "Gadd45a"    "Atp2b2"     "Nek6"      
# [21] "Hat1"       "P2rx3"      "Nr1h3"      "Gss"        "Cyp2j4"    
# [26] "C1qc"       "Fhad1"      "Efnb1"      "Apln"       "Tmem196"   
# [31] "Nudt4"      "Ilf3"       "Htr3a"      "Trip10"     "Ica1l"     
# [36] "Nefh"       "Ugp2"       "RGD1563962" "Fcgr3a"     "Kcnip1"    
# [41] "Kcnh6"      "Idnk"       "Zdhhc2"     "LOC681180"  "Cd200r1"   
# [46] "Rpl17"      "Cndp1"      "Maf"        "Spg7"       "RT1-T24-4" 
# [51] "RT1-Db1"    "Pknox1"     "Grifin"     "P2rx4"  
