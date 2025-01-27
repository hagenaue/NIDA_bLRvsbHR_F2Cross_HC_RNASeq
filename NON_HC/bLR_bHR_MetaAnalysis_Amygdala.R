#Code for running a (very simple) bLR vs. bHR meta-analysis of public amygdala data
#GSE88874 and GSE86893
#Megan Hagenauer
#12-02-2024

#This code is adapted from the 2024 version of the Brain Data Alchemy Project

library(plyr)

install.packages("metafor")

library(metafor)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("multtest")

library(multtest)

##############

#Aligning datasets:

#A function for aligning all of our rat differential expression results from different datasets into a single data frame for Log2FCs and sampling variances (SVs):

ListOfRatDEResults<-list(DEResults_GSE86893, DEResults_GSE88874)

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

# [1] "Rat_MetaAnalysis_FoldChange_Dfs:"
# List of 2
# $ :'data.frame':	17887 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID  : chr [1:17887] "259167" "245961" "64455" "501688" ...
# ..$ GSE86893_bLR_vs_bHR: num [1:17887] -1.743 -1.25 0.439 0.333 -0.227 ...
# $ :'data.frame':	8234 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID  : chr [1:8234] "24153" "24157" "24158" "24159" ...
# ..$ GSE88874_bLR_vs_bHR: num [1:8234] -0.1793 -0.0468 0.1668 -0.7749 -0.2533 ...
# NULL
# [1] "Rat_MetaAnalysis_FoldChanges:"
# 'data.frame':	19103 obs. of  3 variables:
#   $ Rat_EntrezGene.ID  : chr  "259167" "245961" "64455" "501688" ...
# $ GSE86893_bLR_vs_bHR: num  -1.743 -1.25 0.439 0.333 -0.227 ...
# $ GSE88874_bLR_vs_bHR: num  0.2003 0.0882 NA NA NA ...
# NULL
# [1] "Rat_MetaAnalysis_SV_Dfs:"
# List of 2
# $ :'data.frame':	17887 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID  : chr [1:17887] "259167" "245961" "64455" "501688" ...
# ..$ GSE86893_bLR_vs_bHR: num [1:17887] 0.005163 0.013399 0.001894 0.001894 0.000952 ...
# $ :'data.frame':	8234 obs. of  2 variables:
#   ..$ Rat_EntrezGene.ID  : chr [1:8234] "24153" "24157" "24158" "24159" ...
# ..$ GSE88874_bLR_vs_bHR: num [1:8234] 0.0998 0.0568 0.0405 0.0219 0.0271 ...
# NULL
# [1] "Rat_MetaAnalysis_SV:"
# 'data.frame':	19103 obs. of  3 variables:
#   $ Rat_EntrezGene.ID  : chr  "259167" "245961" "64455" "501688" ...
# $ GSE86893_bLR_vs_bHR: num  0.005163 0.013399 0.001894 0.001894 0.000952 ...
# $ GSE88874_bLR_vs_bHR: num  0.072 0.102 NA NA NA ...
# NULL

##############

#I tweaked the next bit of code a bit, because we only have rat datasets (no mouse datasets) in the meta-analysis, so no need for orthologs
#But we probably do want some extra gene annotation eventually.
#And the column names & numbers aren't going to match the meta-analysis code... so I'll need to tweak that next too...

MetaAnalysis_FoldChanges<-Rat_MetaAnalysis_FoldChanges
  
MetaAnalysis_SV<-Rat_MetaAnalysis_SV

##############

#Should we see if there is any correlation between the results in the two studies?
#Probably not - these are pretty dinky little datasets

plot(MetaAnalysis_FoldChanges$GSE86893_bLR_vs_bHR~MetaAnalysis_FoldChanges$GSE88874_bLR_vs_bHR)

#Yep, nada

cor(MetaAnalysis_FoldChanges$GSE86893_bLR_vs_bHR, MetaAnalysis_FoldChanges$GSE88874_bLR_vs_bHR, method="spearman", use="pairwise.complete.obs")
#[1] -0.0629709

##############

#Adapting the meta-analysis function to fit our purposes...
#I needed to change all of the references to the annotation columns to fit the fact that we currently only have rat Entrez ID

str(MetaAnalysis_FoldChanges)
# 'data.frame':	19103 obs. of  3 variables:
# $ Rat_EntrezGene.ID  : chr  "259167" "245961" "64455" "501688" ...
# $ GSE86893_bLR_vs_bHR: num  -1.743 -1.25 0.439 0.333 -0.227 ...
# $ GSE88874_bLR_vs_bHR: num  0.2003 0.0882 NA NA NA ...

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

NumberOfComparisons<-2
CutOffForNAs<-1

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)

# [1] "Table of # of NAs per Row (Gene):"
# MetaAnalysis_FoldChanges_NAsPerRow
# 0     1 
# 7133 11970 
# [1] "MetaAnalysis_FoldChanges_ForMeta:"
# 'data.frame':	7133 obs. of  3 variables:
# $ Rat_EntrezGene.ID  : chr  "259167" "245961" "29705" "266711" ...
# $ GSE86893_bLR_vs_bHR: num  -1.743 -1.25 -0.253 0.185 -3.192 ...
# $ GSE88874_bLR_vs_bHR: num  0.2003 0.0882 -0.9439 -1.0449 0.2622 ...
# NULL


#######

#FDR correction and adding annotation:

#I'm going to need to adapt this Brain Data Alchemy function too...
#I'm going to grab the annotation that we already extracted from org.Rn.eg.db for GSE88874
#That annotation is in the object annotateddata

str(annotateddata)

# 'data.frame':	8234 obs. of  9 variables:
# $ ACCNUM   : chr  "NM_012488" "NM_012489" "NM_016986" "NM_016987" ...
# $ SYMBOL   : chr  "A2m" "Acaa1a" "Acadm" "Acly" ...
# $ ENTREZID : chr  "24153" "24157" "24158" "24159" ...
# $ adj.P.Val: num  0.71995 0.91079 0.58895 0.00294 0.29196 ...
# $ P.Value  : num  0.581 0.848 0.424 0.000222 0.151 0.0122 0.609 0.56 0.0568 0.00163 ...
# $ t        : num  0.567 0.196 -0.829 5.24 1.537 ...
# $ B        : num  -6.693 -6.843 -6.505 0.665 -5.709 ...
# $ logFC    : num  0.1793 0.0468 -0.1668 0.7749 0.2533 ...
# $ GB_ACC   : chr  "NM_012488" "NM_012489" "NM_016986" "NM_016987" ...

str(MetaAnalysis_Annotation)
#chr [1:7133] "259167" "245961" "29705" "266711" "361493" "24252" "297822" ...

temp<-data.frame(ENTREZID=MetaAnalysis_Annotation)
str(temp)
# 'data.frame':	7150 obs. of  1 variable:
#   $ ENTREZID: chr  "259167" "245961" "29705" "266711" ...
MetaAnalysis_Annotation<-temp

colnames(metaOutput)
# [1] "Log2FC_estimate"       "SE"                    "pval"                 
# [4] "CI_lb"                 "CI_ub"                 "Number_Of_Comparisons"

FalseDiscoveryCorrection<-function(metaOutput, annotateddata, MetaAnalysis_Annotation){
  
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
  TempDF3<-join(TempDF, annotateddata[,c(2:3)], by="ENTREZID", type="left", match="first")
  
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

FalseDiscoveryCorrection(metaOutput, annotateddata, MetaAnalysis_Annotation)

# # [1] "metaOutputFDR:"
# num [1:7133, 1:7] -0.789 -0.614 -0.581 -0.399 -1.439 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:7133] "259167" "245961" "29705" "266711" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
# NULL
# [1] "Do we have any genes that are statistically significant following loose false discovery rate correction (FDR<0.10)?"
# [1] 86
# [1] "Do we have any genes that are statistically significant following traditional false discovery rate correction (FDR<0.05)?"
# [1] 67
# [1] "What are the top results?"
# ENTREZID Log2FC_estimate         SE         pval      CI_lb      CI_ub
# X64640     64640       1.0130727 0.14982626 1.364296e-11  0.7194186  1.3067267
# X24252     24252       0.4329265 0.07016338 6.817618e-10  0.2954088  0.5704442
# X297930   297930      -0.3243864 0.05291787 8.788137e-10 -0.4281035 -0.2206692
# X297804   297804      -0.3234212 0.05557449 5.898892e-09 -0.4323452 -0.2144972
# X307989   307989      -0.5660702 0.10164755 2.562771e-08 -0.7652958 -0.3668447
# X25292     25292       1.0548990 0.19106037 3.365267e-08  0.6804276  1.4293704
# Number_Of_Comparisons          FDR SYMBOL
# X64640                      2 9.731523e-08  Rpl30
# X24252                      2 2.089526e-06  Cebpa
# X297930                     2 2.089526e-06   Wwp1
# X297804                     2 1.051920e-05  Plag1
# X307989                     2 3.656050e-05 Ablim1
# X25292                      2 4.000742e-05  Apoc1

#Looks like we have some familiar faces there. :)

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
    
    RatGeneSymbol<-metaOutputFDR_annotated$SYMBOL[which(metaOutputFDR_annotated$ENTREZID==EntrezIDAsCharacter)][1]
    
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

setwd("~/University of Michigan Dropbox/Megan Hagenauer/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/HR_LR_Amygdala_Meta/ForestPlots_topAMYGenes")

#Making forest plots for some of the top genes:

# MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="64640", species="Rat")
# 
# MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="24252", species="Rat")
# 
# MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="297930", species="Rat")
# 
# MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="297804", species="Rat")
# 
# MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="307989", species="Rat")
# 
# MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="25292", species="Rat")
# 
# MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="25404", species="Rat")
# 
# #A few more of the top hits that might be of theoretical interest:
# 
# MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="25404", species="Rat")
# 
# MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="58959", species="Rat")
# 
# MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="24604", species="Rat")
# 
# MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="25555", species="Rat")
# 
# MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="24517", species="Rat")
# 
# MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="29245", species="Rat")

#Maybe I should just loop it

for(i in c(1:length(metaOutputFDR_annotated$ENTREZID[which(metaOutputFDR_annotated$FDR<0.10)]))){
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter=metaOutputFDR_OrderbyPval$ENTREZID[i], species="Rat")
}
#Were clearly still having issues with data from genes with more than one NM

for(i in c(1:length(metaOutputFDR_annotated$ENTREZID[which(metaOutputFDR_annotated$FDR<0.10)]))){
   print(
  metaOutputFDR_annotated$ENTREZID[which(metaOutputFDR_annotated$ENTREZID==metaOutputFDR_OrderbyPval$ENTREZID[i])])
}

# [1] "64640"
# [1] "24252"
# [1] "297930"
# [1] "297804"
# [1] "307989"
# [1] "25292"
# [1] "25404" "25404"
# [1] "64352"
# [1] "114491"
# [1] "293997"
# [1] "291005"
# [1] "362793"
# [1] "296635"
# [1] "85424"
# [1] "116568"
# [1] "24588"
# [1] "140666"
# [1] "56759"
# [1] "65130"
# [1] "500750"
# [1] "365900"
# [1] "171101"
# [1] "300092"
# [1] "303589"
# [1] "303536"
# [1] "29695"
# [1] "305343"
# [1] "317241"
# [1] "25129"
# [1] "29364"
# [1] "24884"
# [1] "58813"
# [1] "293454"
# [1] "116677"
# [1] "360918"
# [1] "171553"
# [1] "89820"
# [1] "246307"
# [1] "25334"
# [1] "246143"
# [1] "29478"
# [1] "25667"
# [1] "316077"
# [1] "297458"
# [1] "29317"
# [1] "170945"
# [1] "29518"
# [1] "25556"
# [1] "116636"
# [1] "116723"
# [1] "288917"
# [1] "63878"
# [1] "59305"
# [1] "499943"
# [1] "83476"
# [1] "25378"
# [1] "29636"
# [1] "498918"
# [1] "315963"
# [1] "25186"
# [1] "498065"
# [1] "81710"
# [1] "79557"
# [1] "29701"
# [1] "24443"
# [1] "298532"
# [1] "116699"
# [1] "310782"
# [1] "64463"
# [1] "680018"
# [1] "54309"
# [1] "24943"
# [1] "56825"
# [1] "363425"
# [1] "114592"
# [1] "304407"
# [1] "298227"
# [1] "83521"
# [1] "494520"
# [1] "29132"
# [1] "304346"
# [1] "289453"
# [1] "360742"
# [1] "25120"
# [1] "29387"
# [1] "25632"

#Only "25404" (Cav1) is a duplicate.

setwd("~/University of Michigan Dropbox/Megan Hagenauer/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/HR_LR_Amygdala_Meta/ForestPlots_topHCGenes")


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

#Genes suggested previously by genetic analyses:

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="83610", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="308543", species="Rat")

#Double-checking genes that were mentioned in the original papers:

#Itgb3bp - listed as down-regulated in the amygdala in Cohen 2017 (GSE86893)
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="362548", species="Rat")
#... and it is down-regulated in the forest plot for that dataset too. Good.

#The Cohen 2015 does not discuss bHR/bLR differential expression, so we can't pull out results from our analysis for comparison

###################


#Interesting: Many of the top genes from the HPC replicate in one dataset (GSE88874) but not the other (GSE86893) - I wonder if the GSE88874 dataset included more amygdala tissue that is HPC-like (e.g., basolateral) vs. striatal-like (CEA, MEA).

#Cohen 2015 (GSE88874) cites this paper for their amygdala dissection, but it really doesn't say anything besides the fact that they used hole punches:
#https://www.sciencedirect.com/science/article/pii/S0006899313011086?via=ihub

#Cohen 2017 (GSE86893) doesn't describe their dissection at all.

#But it is worth noting that GSE88874 was collected as part of an experiment that also include HPC and GSE86893 was collected as part of an experiment that included dorsal raphe, so I could imagine the sections used for punching being more rostral for GSE88874 and more caudal for GSE86893. 

#################

#Comparison to bHR/bLR hippocampal results:

setwd("~/University of Michigan Dropbox/Megan Hagenauer/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/HR_LR_Amygdala_Meta")

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
#[1] 7133    9

metaOutputFDR_moreannotated<-join(metaOutputFDR_annotated, MoreAnnotation, by="ENTREZID", type="left")

dim(metaOutputFDR_moreannotated)
#[1] 7150   10
#So there are a few ENSEMBL ids that are represented multiple times (because they map to more than one ENTREZ id).

colnames(HPC_DEResults)[2]<-"ENSEMBL"

HPC_DEResults_vs_AMY_DEResults<-join(HPC_DEResults, metaOutputFDR_moreannotated, by="ENSEMBL", type="left", match="all")

str(HPC_DEResults_vs_AMY_DEResults)
#'data.frame':	13918 obs. of  46 variables:
13918-13788
#[1] 130
#So there are a few ENSEMBL ids that are represented multiple times 

write.csv(HPC_DEResults_vs_AMY_DEResults, "HPC_DEResults_vs_AMY_DEResults.csv")

pdf("Scatterplot_F0_vs_AMY.pdf", height=5, width=4)
plot(HPC_DEResults_vs_AMY_DEResults$F0_Log2FC_bLRvsbHR~HPC_DEResults_vs_AMY_DEResults$Log2FC_estimate, ylim=c(-4,4), xlab="AMY Log2FC", ylab="F0 HPC Log2FC")
dev.off()
#That's a pretty wussy positive correlation, but a few genes stand out.

cor(HPC_DEResults_vs_AMY_DEResults$F0_Log2FC_bLRvsbHR, HPC_DEResults_vs_AMY_DEResults$Log2FC_estimate, method="spearman", use="pairwise.complete.obs")
#[1] 0.09457585

HPC_DEResults_vs_AMY_DEResults$GENENAME_F2..Rnor6.Ensembl.v.103.[which( HPC_DEResults_vs_AMY_DEResults$F0_Log2FC_bLRvsbHR>0.5 &  HPC_DEResults_vs_AMY_DEResults$Log2FC_estimate>0.5)]
# [1] "Slc22a6"      "Cym"          "Rarres2"      "Ptgds"       
# [5] NA             "A3galt2"      "C1qb"         "Fhad1"       
# [9] "Slc16a8"      "Igfbp6"       "Kif15"        NA            
# [13] "LOC100362027" NA             NA             "Tcam1"       
# [17] "Slc16a3"      "Cd200r1"      "Ttr"          "Grifin"      
# [21] NA 

HPC_DEResults_vs_AMY_DEResults$GENENAME_F2..Rnor6.Ensembl.v.103.[which( HPC_DEResults_vs_AMY_DEResults$F0_Log2FC_bLRvsbHR<(-0.5) &  HPC_DEResults_vs_AMY_DEResults$Log2FC_estimate<(-0.5))]
#[1] "Plagl1"  "Gal"     "Prss12"  "Tes"     NA        "Tmem196" "Cacng5" 

#Interesting.

pdf("Scatterplot_HPC_Meta_vs_AMY.pdf", height=5, width=4)
plot(HPC_DEResults_vs_AMY_DEResults$MetaAnalysis_estimatedD_bLRvsbHR~HPC_DEResults_vs_AMY_DEResults$Log2FC_estimate, xlab="AMY Log2FC", ylab="Meta HPC d")
dev.off()
#Again, basically a blob with a few exceptions

cor(HPC_DEResults_vs_AMY_DEResults$MetaAnalysis_estimatedD_bLRvsbHR, HPC_DEResults_vs_AMY_DEResults$Log2FC_estimate, method="spearman", use="pairwise.complete.obs")
#[1] -0.03315768

HPC_DEResults_vs_AMY_DEResults$GENENAME_F2..Rnor6.Ensembl.v.103.[which( HPC_DEResults_vs_AMY_DEResults$MetaAnalysis_estimatedD_bLRvsbHR>2 &  HPC_DEResults_vs_AMY_DEResults$Log2FC_estimate>0.5)]
#[1] "Tmem176a" "P2rx3"    "C1qb"     "Kif15"    "Rpl17"    "Grifin"   

HPC_DEResults_vs_AMY_DEResults$GENENAME_F2..Rnor6.Ensembl.v.103.[which( HPC_DEResults_vs_AMY_DEResults$MetaAnalysis_estimatedD_bLRvsbHR<(-2) &  HPC_DEResults_vs_AMY_DEResults$Log2FC_estimate<(-0.5))]
#[1] "Prss12" "Tes" 

HPC_DEResults_vs_AMY_DEResults$GENENAME_F2..Rnor6.Ensembl.v.103.[which((HPC_DEResults_vs_AMY_DEResults$F0_FDR_bLRvsbHR<0.1 | HPC_DEResults_vs_AMY_DEResults$MetaAnalysis_FDR_bLRvsbHR<0.1) & HPC_DEResults_vs_AMY_DEResults$pval<0.05)]
# [1] "Syt5"      "Fkrp"      "Hsd3b7"    "Itgam"     "Gal"       "Chka"     
# [7] "Aldh1a1"   "Crhbp"     "Cartpt"    "Golph3l"   "Cym"       "Cxxc4"    
# [13] "Trmt10a"   "Cav2"      "Cav1"      "Akr1b1"    "Tmem176a"  "Atp2b2"   
# [19] "Nek6"      "Hat1"      NA          "P2rx3"     "Nr1h3"     "Gss"      
# [25] "Cyp2j4"    "C1qc"      "Fhad1"     "Efnb1"     "Apln"      "Nudt4"    
# [31] "Ilf3"      "Htr3a"     "Trip10"    "Ica1l"     "Nefh"      "Ugp2"     
# [37] "Rabif"     "Fcgr3a"    "Kcnip1"    "Kcnh6"     "LOC302192" "Cd200r1"  
# [43] "Rpl17"     "Spg7"      "RT1-T24-4" "RT1-Db1"   "Pknox1"    "Grifin"   

#Well, there are some familiar faces - Cav1&Cav2, C1qc, Spg7
#Along with some functionally interesting additional findings - Crhbp, Htr3a, Cartpt, P2rx3

#More strict version:
HPC_DEResults_vs_AMY_DEResults$GENENAME_F2..Rnor6.Ensembl.v.103.[which((HPC_DEResults_vs_AMY_DEResults$F0_FDR_bLRvsbHR<0.1 | HPC_DEResults_vs_AMY_DEResults$MetaAnalysis_FDR_bLRvsbHR<0.1) & HPC_DEResults_vs_AMY_DEResults$FDR<0.1)]
#[1] "Syt5"  "Cym"   "Cav2"  "Cav1"  NA      "Efnb1"

HPC_DEResults_vs_AMY_DEResults$ENTREZID[which((HPC_DEResults_vs_AMY_DEResults$F0_FDR_bLRvsbHR<0.1 | HPC_DEResults_vs_AMY_DEResults$MetaAnalysis_FDR_bLRvsbHR<0.1) & HPC_DEResults_vs_AMY_DEResults$FDR<0.1)]

#[1] "54309"  "56825"  "363425" "25404"  "79557"  "25186" 

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="54309", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="56825", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="25404", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="363425", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="79557", species="Rat")

MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="25186", species="Rat")

#"C1qc" 362634
MakeForestPlots(metaOutputFDR_annotated, EntrezIDAsCharacter="362634", species="Rat")