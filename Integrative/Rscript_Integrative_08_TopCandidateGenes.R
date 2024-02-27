

#So we have multiple noisy pieces of data
#0) F0 genetic results: large effect sizes, but lots of irrelevant genes
#1) F0 & meta-analysis DE: large effect sizes, but lots of irrelevant genes along for the ride. 
#2) allelic FoldChange: large effect sizes, genetic, but lots of irrelevant genes along for the ride, maybe missing some effects due to underpowered eQTL analysis.
#3) F2 DE: small effect sizes, noisy, but actually specific to behavior and less LD, not necessarily genetic effects
#4) F2 SMR: larger effect sizes, genetic, but based on broad QTLs and not filtered for LD (so also noisy)
#5) Birt et al. database of other rat models: very noisy, but probably has different LD than our model, but could suggest generalization

#What I would like: A gene that...
#1) Is differentially expressed in bLR/bHR in a manner that is likely to be due to genetics (congruent aFC, and the variants that drive the aFC are segregated)
#2) Is differentially expressed in the F2s in a manner that is likely to be due to genetics (SMR (but which)?)
#3) Kudos if in Birt et al. database


#The coding for this is *a lot* of conditionals
#To make it easier to code/read
AllResults<-F0_Meta_F2_DEResults_w_F2_eqtls_SMR
str(AllResults)

#Useful from earlier:
bHRGenes
bHRLikeBehaviorGenes
bLRGenes
bLRLikeBehaviorGenes

#Double-checking the number of genes representing cis-eQTLs that are segregated in bHR/bLR:
length(unique(AllResults$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(is.na(AllResults$log2_aFC_bLR_bHR)==FALSE & (is.na(AllResults$Gprimest)==FALSE & AllResults$Gprimest>0.27))]))
#[1] 2395


#As shown earlier, only 108 results (cis-eQTLs) meet the initial filtering criteria based on differential expression in bHR/bLR and with bHR/bLR-like behavior that matches predictions based on segregated bHR/bLR eVariants (aFC):
AllResults$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which((is.na(AllResults$Included.in.F2.Candidate.Gene.Analysis.)==FALSE & AllResults$Included.in.F2.Candidate.Gene.Analysis.=="1") & 
                                                    (is.na(AllResults$log2_aFC_bLR_bHR)==FALSE & ((bHRGenes==TRUE & bHRLikeBehaviorGenes==TRUE & AllResults$log2_aFC_bLR_bHR<0)|(bLRGenes==TRUE & bLRLikeBehaviorGenes==TRUE & AllResults$log2_aFC_bLR_bHR>0)) ) & #because all of the aFC are nominally significant (p<0.05)... not sure why some are so small
                                                    (is.na(AllResults$Gprimest)==FALSE & AllResults$Gprimest>0.27) &
                                                    (is.na(AllResults$F2_Log2FC_LocoScore)==FALSE) &
                                                    (is.na(AllResults$total_loco_score_p_SMR)==FALSE))]



# Daniel used a Bonferonni correction for the SMR results, but actually used FDR for the cis-eQTLs. Not sure why he switched it up.
# And FDR is used throughout the rest of the paper.
# Ah - I think maybe he switched it up because an FDR correction wasn't part of the official output from the SMR analysis & he just wanted a simple threshold

#So let's add an FDR corrected p-value to the SMR results.

library(multtest)

colnames(AllResults)

SMROutput<-AllResults[,c(64,68,72,76,80,84,88,92)]

FalseDiscoveryCorrection<-function(SMROutput,i){
  
  tempPvalAdj<-mt.rawp2adjp(SMROutput[,i], proc=c("BH"))
  
  PvalAdj<-tempPvalAdj$adjp[order(tempPvalAdj$index),]
  
  SMROutputFDR<-cbind(SMROutputFDR, PvalAdj[,2])
  
  colnames(SMROutputFDR)[8+i]<-paste("FDR_", colnames(SMROutput[i]), sep="")
  
  SMROutputFDR<<-SMROutputFDR
  
  print(sum(PvalAdj[,2]<0.05, na.rm=TRUE))
  print(sum(PvalAdj[,2]<0.10, na.rm=TRUE))
  
  rm(tempPvalAdj, PvalAdj)
  
}


SMROutputFDR<-SMROutput

for(i in c(1:ncol(SMROutput))){
  FalseDiscoveryCorrection(SMROutput,i)
}

# [1] 38
# [1] 79
# [1] 0
# [1] 0
# [1] 1
# [1] 1
# [1] 0
# [1] 0
# [1] 0
# [1] 0
# [1] 0
# [1] 13
# [1] 0
# [1] 0
# [1] 0
# [1] 0

colnames(SMROutputFDR)

AllResults<-cbind(AllResults, SMROutputFDR[,c(9:16)])

AllResults$GENENAME..Rnor6.Ensembl.v.88.[which(AllResults$FDR_total_loco_score_p_SMR<0.05)]
# [1] "Gltscr2"        "Kptn"           "Atp5sl"         "Selenov"        "Lsr"            "LOC690000"      "Pop4"           "AC120712.2"    
# [9] "AABR07071904.1" "Siglec5"        "Ptpn5"          "AABR07003304.2" "Gas2"           "Mcee"           "Apba2"          "Mfge8"         
# [17] "Fanci"          "Pex11a"         "Wdr93"          "AABR07004404.1" "Unc45a"         "Hddc3"          "Lipt2"          "Pgm2l1"        
# [25] "AABR07004881.1" "C2cd3"          "Ucp2"           "Coa4"           "Plekhb1"        "Fam168a"        "Tmem9b"         "Fzd6"          
# [33] "Sybu"           "Kcnv1"          "Ctu2"           "Vps9d1"         "Afg3l1"         "RGD1559896" 

AllResults$GENENAME..Rnor6.Ensembl.v.88.[which(AllResults$FDR_total_loco_score_p_SMR<0.10)]
# [1] "Pdcd6"          "Nkd2"           "Rps6ka2"        "Fam120b"        "LOC103689961"   "Gltscr2"        "Napa"           "Kptn"          
# [9] "Calm2"          "Rsph6a"         "Atp5sl"         "RGD1560854"     "Selenov"        "Eid2"           "Mrps12"         "Lsr"           
# [17] "LOC690000"      "Pop4"           "AC120712.2"     "AABR07071904.1" "Siglec5"        "Dzf17"          "LOC102546648"   "Ptpn5"         
# [25] "Zdhhc13"        "AABR07003304.2" "Prmt3"          "Gas2"           "Mtmr10"         "Mcee"           "Mcee"           "Apba2"         
# [33] "Fam174b"        "Slco3a1"        "Mrps11"         "Aen"            "Mfge8"          "Fanci"          "Polg"           "Pex11a"        
# [41] "Wdr93"          "AABR07004404.1" "Vps33b"         "Unc45a"         "Hddc3"          "Lipt2"          "Pgm2l1"         "AABR07004881.1"
# [49] "P4ha3"          "C2cd3"          "Ucp2"           "Dnajb13"        "Coa4"           "Plekhb1"        "Fam168a"        "Tmem9b"        
# [57] "Nrip3"          "Adm"            "Pdk1"           "Map3k20"        "Rgs22"          "Fzd6"           "Slc25a32"       "Nudcd1"        
# [65] "Ebag9"          "Sybu"           "Kcnv1"          "Atp7b"          "LOC102556347"   "Ist1"           "Txnl4b"         "Klhdc4"        
# [73] "Zfpm1"          "Ctu2"           "Spg7"           "Vps9d1"         "Afg3l1"         "RGD1559896"     "Gnpat"


AllResults$GENENAME..Rnor6.Ensembl.v.88.[which(AllResults$FDR_epm_distance_traveled_p_SMR<0.05)]
#[1] "AABR07071904.1"

AllResults$GENENAME..Rnor6.Ensembl.v.88.[which(AllResults$FDR_of_distance_traveled_p_SMR<0.10)]
# [1] "Gas2"           "Apba2"          "Mfge8"          "Fanci"          "Pex11a"         "Wdr93"          "Arpin"          "AABR07004404.1"
# [9] "Hddc3"          "Lipt2"          "Pgm2l1"         "Ucp2"           "Plekhb1"   



CandidateGene_FDR05<-(is.na(AllResults$Included.in.F2.Candidate.Gene.Analysis.)==FALSE & AllResults$Included.in.F2.Candidate.Gene.Analysis.=="1") & 
  (is.na(AllResults$log2_aFC_bLR_bHR)==FALSE & ((bHRGenes==TRUE & AllResults$log2_aFC_bLR_bHR<0)|(bLRGenes==TRUE & AllResults$log2_aFC_bLR_bHR>0)) ) & #because all of the aFC are nominally significant (p<0.05)... not sure why some are so small
  (is.na(AllResults$Gprimest)==FALSE & AllResults$Gprimest>0.27) &
  (is.na(AllResults$F2_Log2FC_LocoScore)==FALSE) &
  (is.na(AllResults$total_loco_score_p_SMR)==FALSE) &
  (
    (AllResults$F2_Pvalue_LocoScore<0.05 & AllResults$FDR_total_loco_score_p_SMR<0.05 & (AllResults$F2_Log2FC_LocoScore*AllResults$total_loco_score_T_SMR_directionalv2)>0)|
      (AllResults$F2_Pvalue_EPM_DistanceTraveled<0.05 & ((AllResults$FDR_epm_distance_traveled_p_SMR<0.05 & (AllResults$F2_Log2FC_EPM_DistanceTraveled*AllResults$epm_distance_traveled_T_SMR_directionalv2)>0)|(AllResults$FDR_of_distance_traveled_p_SMR<0.05 & (AllResults$F2_Log2FC_EPM_DistanceTraveled*AllResults$of_distance_traveled_T_SMR_directionalv2)>0)))|
      (AllResults$F2_Pvalue_EPM_Time_Immobile<0.05 & ((AllResults$FDR_epm_time_immobile_s_p_SMR<0.05 & (AllResults$F2_Log2FC_EPM_Time_Immobile*AllResults$epm_time_immobile_s_T_SMR_directionalv2)>0)|(AllResults$FDR_of_time_immobile_s_p_SMR<0.05 & (AllResults$F2_Log2FC_EPM_Time_Immobile*AllResults$of_time_immobile_s_T_SMR_directionalv2)>0)))|
      (AllResults$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & ((AllResults$FDR_epm_percent_time_open_arm_p_SMR<0.05 & (AllResults$F2_Log2FC_EPM_Percent_Time_Open_Arms*AllResults$epm_percent_time_open_arm_T_SMR_directionalv2)>0)|(AllResults$FDR_of_percent_time_in_center_p_SMR<0.05 & (AllResults$F2_Log2FC_EPM_Percent_Time_Open_Arms*AllResults$of_percent_time_in_center_T_SMR_directionalv2)>0)))|
      (AllResults$F2_Pvalue_PCA_Index<0.05 & (AllResults$FDR_pca_index_days_6and7_p_SMR<0.05 & (AllResults$F2_Log2FC_PCA_Index*AllResults$pca_index_days_6and7_T_SMR_directionalv2)>0))
  )


head(CandidateGene_FDR05)
sum(CandidateGene_FDR05, na.rm=TRUE)
#[1] 14


cbind(AllResults$ENSEMBL.ID..Rnor6.Ensembl.v.103.[CandidateGene_FDR05], AllResults$GENENAME..Rnor6.Ensembl.v.88.[CandidateGene_FDR05])
#     [,1]                 [,2]            
# [1,] "ENSRNOG00000021053" "Lsr"           
# [2,] "ENSRNOG00000052237" "AABR07071904.1"
# [3,] "ENSRNOG00000017510" "Mfge8"         
# [4,] "ENSRNOG00000016689" "Fanci"         
# [5,] "ENSRNOG00000015003" "Pex11a"        
# [6,] "ENSRNOG00000026514" "Wdr93"         
# [7,] "ENSRNOG00000012357" "Unc45a"        
# [8,] "ENSRNOG00000016906" "Lipt2"         
# [9,] "ENSRNOG00000017608" "C2cd3"         
# [10,] "ENSRNOG00000017854" "Ucp2"          
# [11,] "ENSRNOG00000018627" "Plekhb1"       
# [12,] "ENSRNOG00000004660" "Fzd6"          
# [13,] "ENSRNOG00000028904" "Vps9d1"        
# [14,] "ENSRNOG00000026994" "Afg3l1" 


#FDR10 version:

CandidateGene_FDR10<-(is.na(AllResults$Included.in.F2.Candidate.Gene.Analysis.)==FALSE & AllResults$Included.in.F2.Candidate.Gene.Analysis.=="1") & 
  (is.na(AllResults$log2_aFC_bLR_bHR)==FALSE & ((bHRGenes==TRUE & AllResults$log2_aFC_bLR_bHR<0)|(bLRGenes==TRUE & AllResults$log2_aFC_bLR_bHR>0)) ) & #because all of the aFC are nominally significant (p<0.05)... not sure why some are so small
  (is.na(AllResults$Gprimest)==FALSE & AllResults$Gprimest>0.27) &
  (is.na(AllResults$F2_Log2FC_LocoScore)==FALSE) &
  (is.na(AllResults$total_loco_score_p_SMR)==FALSE) &
  (
    (AllResults$F2_Pvalue_LocoScore<0.05 & AllResults$FDR_total_loco_score_p_SMR<0.10 & (AllResults$F2_Log2FC_LocoScore*AllResults$total_loco_score_T_SMR_directionalv2)>0)|
      (AllResults$F2_Pvalue_EPM_DistanceTraveled<0.05 & ((AllResults$FDR_epm_distance_traveled_p_SMR<0.10 & (AllResults$F2_Log2FC_EPM_DistanceTraveled*AllResults$epm_distance_traveled_T_SMR_directionalv2)>0)|(AllResults$FDR_of_distance_traveled_p_SMR<0.10 & (AllResults$F2_Log2FC_EPM_DistanceTraveled*AllResults$of_distance_traveled_T_SMR_directionalv2)>0)))|
      (AllResults$F2_Pvalue_EPM_Time_Immobile<0.05 & ((AllResults$FDR_epm_time_immobile_s_p_SMR<0.10 & (AllResults$F2_Log2FC_EPM_Time_Immobile*AllResults$epm_time_immobile_s_T_SMR_directionalv2)>0)|(AllResults$FDR_of_time_immobile_s_p_SMR<0.10 & (AllResults$F2_Log2FC_EPM_Time_Immobile*AllResults$of_time_immobile_s_T_SMR_directionalv2)>0)))|
      (AllResults$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & ((AllResults$FDR_epm_percent_time_open_arm_p_SMR<0.10 & (AllResults$F2_Log2FC_EPM_Percent_Time_Open_Arms*AllResults$epm_percent_time_open_arm_T_SMR_directionalv2)>0)|(AllResults$FDR_of_percent_time_in_center_p_SMR<0.10 & (AllResults$F2_Log2FC_EPM_Percent_Time_Open_Arms*AllResults$of_percent_time_in_center_T_SMR_directionalv2)>0)))|
      (AllResults$F2_Pvalue_PCA_Index<0.05 & (AllResults$FDR_pca_index_days_6and7_p_SMR<0.10 & (AllResults$F2_Log2FC_PCA_Index*AllResults$pca_index_days_6and7_T_SMR_directionalv2)>0))
  )


head(CandidateGene_FDR10)
sum(CandidateGene_FDR10, na.rm=TRUE)
#[1] 16

cbind(AllResults$ENSEMBL.ID..Rnor6.Ensembl.v.103.[CandidateGene_FDR10], AllResults$GENENAME..Rnor6.Ensembl.v.88.[CandidateGene_FDR10])
# [1,] "ENSRNOG00000021053" "Lsr"           
# [2,] "ENSRNOG00000052237" "AABR07071904.1"
# [3,] "ENSRNOG00000017510" "Mfge8"         
# [4,] "ENSRNOG00000016689" "Fanci"         
# [5,] "ENSRNOG00000015003" "Pex11a"        
# [6,] "ENSRNOG00000026514" "Wdr93"         
# [7,] "ENSRNOG00000012357" "Unc45a"        
# [8,] "ENSRNOG00000016906" "Lipt2"         
# [9,] "ENSRNOG00000017608" "C2cd3"         
# [10,] "ENSRNOG00000017854" "Ucp2"          
# [11,] "ENSRNOG00000018627" "Plekhb1"       
# [12,] "ENSRNOG00000004660" "Fzd6"          
# [13,] "ENSRNOG00000015144" "Ist1"          
# [14,] "ENSRNOG00000015150" "Spg7"          
# [15,] "ENSRNOG00000028904" "Vps9d1"        
# [16,] "ENSRNOG00000026994" "Afg3l1" 


write.csv(AllResults, "AllResults_20240129.csv")


#For adults, only locoscore and distance traveled survive FDR correction.
#And only one gene survives for FDR for distance traveled (AABR07071904.1)
#For LocoScore: FDR<0.05=p<0.000119 and FDR<0.10=p<0.000546
#For juveniles, only of distance traveled survives FDR correction.
#no genes with FDR<0.05, FDR<0.10=p<9.05e-05


######################