#NIDA U01 project: Genetics and Transcriptomics of bHR/bLR F0-F1-F2cross
#This code focuses on the analysis of the hippocampal RNA-Seq data from the F0 male and female bHR/bLR rats

#Elaine Hebda-Bauer and Megan Hagenauer
#Nov 2021-Jan 2022

########################################

#Project Background:

#This study was performed using whole hippocampal tissue dissected by from F0 male and female bHR/bLR rats by  Peter Blandino, but the RNA processing and sequencing was organized by Elaine Hebda-Bauer and processed in tandem at the discovery life sciences Core facility with the F2 hippocampal RNA-Seq study (but in a separate sequencing batch)

#Previously, there was also a smaller hippocampal RNA-Seq study performed on F0 (F37) bHR/bLR males by Peter Blandino. That data was analyzed as part of the Birt Hagenauer et al. 2021 paper and is not included in the current analysis.

#Analysis version information:

#prior to 11/2021: Differential expression analysis performed automatically on Ahub. Because it was automatic, there was not enough QC, no evaluation of noise/confounds, and ended with awkward model construction. 
#11/2021: Began re-analysis in R, but decided to add TMM normalize, analyze F2s, and come back to it.
#12/2021: Performed analysis in R, but missed a few QC steps and potential technical co-variates
#1/2022: Redid analysis in R, more QC steps and potential co-variates - performed quickly/slap-dash to meet abstract deadline
#1/10/2022:  Re-Did analysis carefully, removing one low RIN sample, using Megan's PCA code, etc

########################################

#To do list:
#Package and version information needs updating
#Find out how low the RIN value is for sl483897 in the F0 dataset, output range of RIN for remaining samples.
#Double-check remaining sample size.
#Output full stats to back up relationship between Lineage and some of the co-variates of interest.
#Output differential expression results from a model with no co-variates.
#Output a quick formal comparison of the results from the models examined (correlation matrix)
#Make basic graphs of top bHR/bLR genes based on the previous meta-analysis (now that the data is QC-ed).  
#Make pretty ones for some of our top F0/F2 combined genes.
#Do we need to re-do the F0 Volcano plots??
#Do we also need to re-do the F0 vs. Meta-analysis plots?

########################################
