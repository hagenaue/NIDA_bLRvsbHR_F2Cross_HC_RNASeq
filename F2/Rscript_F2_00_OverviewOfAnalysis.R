#NIDA U01 project: Genetics and Transcriptomics of bHR/bLR F0-F1-F2cross
#This code focuses on the analysis of the hippocampal RNA-Seq data from the F2 male and female rats in relationship to locomotor activity, behavior on the elevated plus maze (EPM), and sign-tracking/goal-tracking

#Elaine Hebda-Bauer and Megan Hagenauer
#Nov 2021-Jan 2022

########################################

#Project Background:

#This study was performed using hippocampal tissue punches dissected by from F2 male and female bHR/bLR rats under the supervision of Elaine Hebda-Bauer

#The RNA processing and sequencing was organized by Elaine Hebda-Bauer and processed in tandem at the discovery life sciences Core facility with the F2 hippocampal RNA-Seq study (but in a separate sequencing batch)

#Analysis version information:

#prior to 11/2021: Differential expression analysis performed automatically on Ahub. Because it was automatic, there was not enough QC, no evaluation of noise/confounds, and ended with awkward model construction. 
#12/2021: Performed analysis in R, but missed a few QC steps and potential technical co-variates
#1/2022: Redid analysis in R, more QC steps and potential co-variates - performed quickly/slap-dash to meet abstract deadline
#1/10/2022:  Re-Did analysis carefully, removing two samples with mislabeling, using Megan's PCA code, etc

############


#To do list:

#Package and version information needs updating
#Fill in section of methods with potential confounding relationships (variables of interest vs. all other variables)
#Fill in stats for relationship between PC1 and %intergenic and PC3 and Sex
#Define R cut-off for eliminating variables that were redundant with each other:
# - %intergenic vs. %mRNA = R=-0.97 - choose %intergenic (correlates with PC1)
# - % coding vs. %MRNA = R=0.85 - choose % coding (correlated with PC5)
# - dissection day vs. sequencing ID = R=0.90 - choose dissection day (more correlated with PC1, PC2, and PC6)
# - library size vs. pf.reads = 0.99 - choose library size.
# - Sex vs. Age = 0.85
#Maybe re-run stepwise with the redundant variables removed--COMPLETED 1/21/2022
#Add results for EPM and STGT--COMPLETED 1/21/2022
#Compare the EPM and STGT results to F0/meta-analysis
#Make basic graphs of top bHR/bLR genes - clean data first? 
#Make pretty ones for some of our top F0/F2 combined genes.
#Make Volcano plots
#Output data cleaned of confounds

############
