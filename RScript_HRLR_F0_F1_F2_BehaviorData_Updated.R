
###bHR/bLR F0 F1 F2 Behavioral Data Analysis
#Beginning 2/9/2022


setwd("//n05-corea-cifs.umhs.med.umich.edu/Home1/Users/hebda/My Documents/HRLR_F0_F1_F2_BehaviorData")

getwd()
#[1] "\\\\n05-corea-cifs.umhs.med.umich.edu/Home1/Users/hebda/My Documents/HRLR_F0_F1_F2_BehaviorData"


Behavior<-read.csv("Akil_HRLR_F0_F1_F2_PhenotypeData_ForR.csv", header = TRUE, stringsAsFactors = FALSE)

str(Behavior)
#'data.frame':	849 obs. of  38 variables:
# $ RatNumForNIHstudy        : int  1 2 3 4 5 6 7 8 9 10 ...
# $ RatUniqueID              : int  13148 12941 12932 13155 13177 12965 12960 13176 13204 12972 ...
# $ WGS                      : chr  "Jun" "Jun" "Jun" "Jun" ...
# $ LibraryID                : chr  "SL483888" "" "" "SL483889" ...
# $ RNASeq                   : chr  "Elaine" "Peter" "" "Elaine" ...
# $ Adult_Behavioral_Data    : chr  "Yes" "Yes" "Yes" "Yes" ...
# $ HaveST_GT_Data           : chr  "No" "No" "No" "No" ...
# $ DOB                      : chr  "11/13/2012" "11/11/2012" "11/11/2012" "11/13/2012" ...
# $ DOD                      : chr  "5/22/2013" "4/25/2013" "5/22/2013" "4/25/2013" ...
# $ Age_days                 : int  189 164 191 162 188 163 190 164 192 163 ...
# $ Mother                   : int  NA NA NA NA NA NA NA NA NA NA ...
# $ Father                   : int  NA NA NA NA NA NA NA NA NA NA ...
# $ On_Pedigree              : chr  "Yes" "Yes" "Yes" "Yes" ...
# $ Phenotype                : chr  "LR" "HR" "HR" "LR" ...
# $ Family                   : int  1 1 1 1 2 2 2 2 3 3 ...
# $ Study_Generation         : chr  "F0" "F0" "F0" "F0" ...
# $ Sex                      : chr  "Female" "Male" "Female" "Male" ...
# $ Lateral_LocoScore        : int  12 747 833 12 14 715 890 22 17 827 ...
# $ Rearing_LocoScore        : int  11 1606 1448 18 6 1480 2059 22 21 1674 ...
# $ Total_LocoScore          : int  23 2353 2267 30 20 2195 2949 44 38 2501 ...
# $ EPM_Percent_Time_Open_Arm: num  2.62 32.16 64.36 5.16 11.05 ...
# $ EPM_Boli                 : int  1 0 0 4 6 0 0 5 1 0 ...
# $ EPM_Distance_Traveled    : num  1730 4212 2675 4126 1700 ...
# $ EPM_Time_Immobile        : num  29.6 32.7 34.9 4.4 20.2 ...
# $ OF_Percent_Time_In_Center: num  NA NA NA NA NA NA NA NA NA NA ...
# $ OF_Boli                  : int  NA NA NA NA NA NA NA NA NA NA ...
# $ OF_Distance_Traveled     : num  NA NA NA NA NA NA NA NA NA NA ...
# $ OF_Time_Immobile         : num  NA NA NA NA NA NA NA NA NA NA ...
# $ LearningClassification   : chr  "" "" "" "" ...
# $ PCA_Index_Days6and.7     : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Response_Bias_Day7       : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Prob_Diff_Day7           : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Latency_Score_Day7       : num  NA NA NA NA NA NA NA NA NA NA ...
# $ PCA_Index_Day7only       : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Response_Bias_Day6       : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Prob_Diff_Day6           : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Latency_Score_Day6       : num  NA NA NA NA NA NA NA NA NA NA ...
# $ PCA_Index_Day6only       : num  NA NA NA NA NA NA NA NA NA NA ...

table(Behavior$EPM_Boli)





Behavior$Sex[Behavior$Sex=="Female"]<-"F"

Behavior$Sex[Behavior$Sex=="Male"]<-"M"

table(Behavior$Study_Generation)

Behavior$Phenotype[Behavior$Phenotype=="I1"]<-NA

Behavior$Phenotype[Behavior$Phenotype=="I2"]<-NA

Behavior$Phenotype[is.na(Behavior$Phenotype)]<-""


Behavior$Sex_AsFactor<-as.factor(Behavior$Sex)
levels(Behavior$Sex_AsFactor)
#[1] "F" "M"

#Setting the reference levels for factors:
Behavior$Sex_AsFactor<-relevel(Behavior$Sex_AsFactor, ref="M")
levels(Behavior$Sex_AsFactor)
#[1] "M" "F"


Behavior$LearningClassification_AsFactor<-as.factor(Behavior$LearningClassification)
levels(Behavior$LearningClassification_AsFactor)
#[1] ""   "GT" "IN" "ST"

#Setting the reference levels for factors:
Behavior$LearningClassification_AsFactor<-relevel(Behavior$LearningClassification_AsFactor, ref="ST")
levels(Behavior$LearningClassification_AsFactor)
#[1] "ST" ""   "GT" "IN"



str(Behavior)
# 'data.frame':	849 obs. of  41 variables:
#   $ RatNumForNIHstudy              : int  1 2 3 4 5 6 7 8 9 10 ...
# $ RatUniqueID                    : int  13148 12941 12932 13155 13177 12965 12960 13176 13204 12972 ...
# $ WGS                            : chr  "Jun" "Jun" "Jun" "Jun" ...
# $ LibraryID                      : chr  "SL483888" "" "" "SL483889" ...
# $ RNASeq                         : chr  "Elaine" "Peter" "" "Elaine" ...
# $ Adult_Behavioral_Data          : chr  "Yes" "Yes" "Yes" "Yes" ...
# $ HaveST_GT_Data                 : chr  "No" "No" "No" "No" ...
# $ DOB                            : chr  "11/13/2012" "11/11/2012" "11/11/2012" "11/13/2012" ...
# $ DOD                            : chr  "5/22/2013" "4/25/2013" "5/22/2013" "4/25/2013" ...
# $ Age_days                       : int  189 164 191 162 188 163 190 164 192 163 ...
# $ Mother                         : int  NA NA NA NA NA NA NA NA NA NA ...
# $ Father                         : int  NA NA NA NA NA NA NA NA NA NA ...
# $ On_Pedigree                    : chr  "Yes" "Yes" "Yes" "Yes" ...
# $ Phenotype                      : chr  "LR" "HR" "HR" "LR" ...
# $ Family                         : int  1 1 1 1 2 2 2 2 3 3 ...
# $ Study_Generation               : chr  "F0" "F0" "F0" "F0" ...
# $ Sex                            : chr  "F" "M" "F" "M" ...
# $ Lateral_LocoScore              : int  12 747 833 12 14 715 890 22 17 827 ...
# $ Rearing_LocoScore              : int  11 1606 1448 18 6 1480 2059 22 21 1674 ...
# $ Total_LocoScore                : int  23 2353 2267 30 20 2195 2949 44 38 2501 ...
# $ EPM_Percent_Time_Open_Arm      : num  2.62 32.16 64.36 5.16 11.05 ...
# $ EPM_Boli                       : int  1 0 0 4 6 0 0 5 1 0 ...
# $ EPM_Distance_Traveled          : num  1730 4212 2675 4126 1700 ...
# $ EPM_Time_Immobile              : num  29.6 32.7 34.9 4.4 20.2 ...
# $ OF_Percent_Time_In_Center      : num  NA NA NA NA NA NA NA NA NA NA ...
# $ OF_Boli                        : int  NA NA NA NA NA NA NA NA NA NA ...
# $ OF_Distance_Traveled           : num  NA NA NA NA NA NA NA NA NA NA ...
# $ OF_Time_Immobile               : num  NA NA NA NA NA NA NA NA NA NA ...
# $ LearningClassification         : chr  "" "" "" "" ...
# $ PCA_Index_Days6and.7           : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Response_Bias_Day7             : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Prob_Diff_Day7                 : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Latency_Score_Day7             : num  NA NA NA NA NA NA NA NA NA NA ...
# $ PCA_Index_Day7only             : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Response_Bias_Day6             : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Prob_Diff_Day6                 : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Latency_Score_Day6             : num  NA NA NA NA NA NA NA NA NA NA ...
# $ PCA_Index_Day6only             : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Sex_AsFactor                   : Factor w/ 2 levels "M","F": 2 1 2 1 2 1 2 1 2 1 ...
# $ GenPheno                       : chr  "F0 LR" "F0 HR" "F0 HR" "F0 LR" ...
# $ LearningClassification_AsFactor: Factor w/ 4 levels "ST","","GT","IN": 2 2 2 2 2 2 2 2 2 2 ...

Behavior$GenPheno<-paste(Behavior$Study_Generation,Behavior$Phenotype,sep = " ")

#########################################These work with Sex as Character

pdf("TotalLocoScore_All_SexGen.pdf", width=9.5, height=6)
boxplot(Total_LocoScore~Sex*GenPheno, data=Behavior, ylab="LocoScore", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(Total_LocoScore~Sex*GenPheno, data=Behavior, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()

#col=c("green4","green2", "grey48", "grey","firebrick4", "firebrick2", "gold3", "gold")


pdf("TotalLocoScore_Abe_SexGen.pdf", width=9.5, height=6)
boxplot(Total_LocoScore~Sex*GenPheno, data=Behavior[Behavior$gDNAtoPalmerLab=="Yes",], ylab="LocoScore", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(Total_LocoScore~Sex*GenPheno, data=Behavior[Behavior$gDNAtoPalmerLab=="Yes",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()


pdf("TotalLocoScore_Abe_SexGen_STGT.pdf", width=9.5, height=6)
boxplot(Total_LocoScore~Sex*GenPheno, data=Behavior[Behavior$gDNAtoPalmerLab=="Yes"&Behavior$HaveST_GT_Data=="Yes",], ylab="LocoScore", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(Total_LocoScore~Sex*GenPheno, data=Behavior[Behavior$gDNAtoPalmerLab=="Yes"&Behavior$HaveST_GT_Data=="Yes",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()


Behavior_Abe_STGT<-Behavior[Behavior$gDNAtoPalmerLab=="Yes"&Behavior$HaveST_GT_Data=="Yes",]

table(Behavior_Abe_STGT$LearningClassification, Behavior_Abe_STGT$Sex, Behavior_Abe_STGT$GenPheno)

###########################################All Rats here


pdf("TotalLocoScore_All_SexGen_MaleRef.pdf", width=9.5, height=6)
boxplot(Total_LocoScore~Sex_AsFactor*GenPheno, data=Behavior, ylab="LocoScore", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(Total_LocoScore~Sex_AsFactor*GenPheno, data=Behavior, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()


#########EPM_Percent_Time_Open_Arm--3 with 87, 80, and 75% (the next highest is 65%); only the 87% rat (#120) was used for Abe's WGS, none for RNASeq
pdf("EPM_Percent_Time_Open_Arm_All_SexGen_MaleRef.pdf", width=9.5, height=6)
boxplot(EPM_Percent_Time_Open_Arm~Sex_AsFactor*GenPheno, data=Behavior, ylab="Percent Time in Open Arms", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart( EPM_Percent_Time_Open_Arm~Sex_AsFactor*GenPheno, data=Behavior, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()


######Checkout 5 M-F1 outliers For EPM distance traveled--None are in the RNAseq or Palmer's WGS samples
pdf("EPM_Distance_Traveled_All_SexGen_MaleRef.pdf", width=9.5, height=6)
boxplot(EPM_Distance_Traveled~Sex_AsFactor*GenPheno, data=Behavior, ylab="Distance Traveled in EPM", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart( EPM_Distance_Traveled~Sex_AsFactor*GenPheno, data=Behavior, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

pdf("EPM_Distance_Traveled_AllnoOutlier_SexGen_MaleRef.pdf", width=9.5, height=6)
boxplot(EPM_Distance_Traveled~Sex_AsFactor*GenPheno, data=Behavior[Behavior$EPM_Distance_Traveled<5000,], ylab="Distance Traveled in EPM", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart( EPM_Distance_Traveled~Sex_AsFactor*GenPheno, data=Behavior[Behavior$EPM_Distance_Traveled<5000,], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()


pdf("EPM_Time_Immobile_All_SexGen_MaleRef.pdf", width=9.5, height=6)
boxplot(EPM_Time_Immobile~Sex_AsFactor*GenPheno, data=Behavior, ylab="Time Immobile in EPM", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(EPM_Time_Immobile~Sex_AsFactor*GenPheno, data=Behavior, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()


pdf("PCA_Index_Days6and.7_All_SexGen_MaleRef.pdf", width=9.5, height=6)
boxplot(PCA_Index_Days6and.7~Sex_AsFactor*GenPheno, data=Behavior, ylab="PCA_Index_Days6and.7", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(PCA_Index_Days6and.7~Sex_AsFactor*GenPheno, data=Behavior, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()


pdf("PCA_Index_Days6and.7_All_SexGen_MaleRef_LearnClass_Yaxis_Title.pdf", width=9.5, height=6)
boxplot(PCA_Index_Days6and.7~Sex_AsFactor*GenPheno, data=Behavior, ylab="Learning Classification", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(PCA_Index_Days6and.7~Sex_AsFactor*GenPheno, data=Behavior, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()


#########Have GT, IN, and ST only here as desired, but would also like to have ST as the reference level
pdf("Total_LocoScore_LearnClass_F2.pdf", width=9.5, height=8)
boxplot(Total_LocoScore~LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], ylab="LocoScore", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(Total_LocoScore~LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()

########Have ST as the reference level here, but now the samples without ST/GT data show up here as a group
pdf("Total_LocoScore_LearnClass_AsFactor_F2.pdf", width=9.5, height=8)
boxplot(Total_LocoScore~LearningClassification_AsFactor, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], ylab="LocoScore", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(Total_LocoScore~LearningClassification_AsFactor, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()

####################Rats for RNASeq--24 F0s and 250 F2s


pdf("TotalLocoScore_RNAseq_SexGen_MaleRef.pdf", width=9.5, height=6)
boxplot(Total_LocoScore~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], ylab="LocoScore", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(Total_LocoScore~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()


pdf("EPM_Percent_Time_Open_Arm_RNAseq_SexGen_MaleRef.pdf", width=9.5, height=6)
boxplot(EPM_Percent_Time_Open_Arm~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], ylab="Percent Time in Open Arms", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart( EPM_Percent_Time_Open_Arm~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()


pdf("EPM_Distance_Traveled_RNAseq_SexGen_MaleRef.pdf", width=9.5, height=6)
boxplot(EPM_Distance_Traveled~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], ylab="Distance Traveled in EPM", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart( EPM_Distance_Traveled~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()


pdf("EPM_Distance_Traveled_noOutlier_RNAseq_SexGen_MaleRef.pdf", width=9.5, height=6)
boxplot(EPM_Distance_Traveled~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine" & Behavior$EPM_Distance_Traveled<5000,], ylab="Distance Traveled in EPM", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart( EPM_Distance_Traveled~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine" & Behavior$EPM_Distance_Traveled<5000,], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()



pdf("EPM_Time_Immobile__RNAseq_SexGen_MaleRef.pdf", width=9.5, height=6)
boxplot(EPM_Time_Immobile~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], ylab="Time Immobile in EPM", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(EPM_Time_Immobile~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()


pdf("PCA_Index_Days6and.7__RNAseq_SexGen_MaleRef.pdf", width=9.5, height=6)
boxplot(PCA_Index_Days6and.7~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], ylab="PCA_Index_Days6and.7", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(PCA_Index_Days6and.7~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()


pdf("PCA_Index_Days6and.7_All_SexGen_MaleRef_LearnClass_Yaxis_Title.pdf", width=9.5, height=6)
boxplot(PCA_Index_Days6and.7~Sex_AsFactor*GenPheno, data=Behavior, ylab="Learning Classification", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(PCA_Index_Days6and.7~Sex_AsFactor*GenPheno, data=Behavior, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()


#########Have GT, IN, and ST only here as desired, but would also like to have ST as the reference level
#This generates same boxplot as above since only F2 rats had ST/GT data and all were used for RNASeq
pdf("Total_LocoScore_LearnClass_F2.pdf", width=9.5, height=8)
boxplot(Total_LocoScore~LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], ylab="LocoScore", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(Total_LocoScore~LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

########Have ST as the reference level here, but now the samples without ST/GT data show up here as a group
#This generates same boxplot as above since only F2 rats had ST/GT data and all were used for RNASeq
pdf("Total_LocoScore_LearnClass_AsFactor_F2.pdf", width=9.5, height=8)
boxplot(Total_LocoScore~LearningClassification_AsFactor, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], ylab="LocoScore", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(Total_LocoScore~LearningClassification_AsFactor, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()


####################2/10/2022  Rats for RNASeq--24 F0s and 250 F2s--These are color boxplots

#col=c("green4","green2", "grey48", "grey","firebrick4", "firebrick2", "gold3", "gold")

#col=c("green4","green2","firebrick4", "firebrick2", "grey57", "grey80")


pdf("TotalLocoScore_RNAseq_SexGen_MaleRef_Color.pdf", width=9.5, height=6)
boxplot(Total_LocoScore~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], ylab="LocoScore", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("green4","green2","firebrick4", "firebrick2", "grey57", "grey80"))
stripchart(Total_LocoScore~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(Total_LocoScore~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], contrasts=list(GenPheno=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: Total_LocoScore
#                         Sum Sq  Df  F value  Pr(>F)    
# (Intercept)           46717433   1 538.7853 < 2e-16 ***
# Sex_AsFactor            331419   1   3.8222 0.05162 .  
# GenPheno              36897051   2 212.7641 < 2e-16 ***
# Sex_AsFactor:GenPheno   358060   2   2.0647 0.12887    
# Residuals             23237962 268                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


pdf("EPM_Percent_Time_Open_Arm_RNAseq_SexGen_MaleRef_Color.pdf", width=9.5, height=6)
boxplot(EPM_Percent_Time_Open_Arm~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], ylab="Percent Time in Open Arms", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("green4","green2","firebrick4", "firebrick2", "grey57", "grey80"))
stripchart( EPM_Percent_Time_Open_Arm~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(EPM_Percent_Time_Open_Arm~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], contrasts=list(GenPheno=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: EPM_Percent_Time_Open_Arm
#                          Sum Sq  Df  F value    Pr(>F)    
# (Intercept)             14198.9   1 129.2723 < 2.2e-16 ***
#   Sex_AsFactor            832.8   1   7.5818  0.006299 ** 
#   GenPheno               3105.6   2  14.1372 1.456e-06 ***
#   Sex_AsFactor:GenPheno  1078.2   2   4.9082  0.008063 ** 
#   Residuals             29436.4 268                       
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


pdf("EPM_Distance_Traveled_RNAseq_SexGen_MaleRef_Color.pdf", width=9.5, height=6)
boxplot(EPM_Distance_Traveled~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], ylab="Distance Traveled in EPM", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("green4","green2","firebrick4", "firebrick2", "grey57", "grey80"))
stripchart( EPM_Distance_Traveled~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(EPM_Distance_Traveled~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], contrasts=list(GenPheno=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: EPM_Distance_Traveled
#                             Sum Sq  Df  F value    Pr(>F)    
# (Intercept)              277178343   1 1121.509 < 2.2e-16 ***
#   Sex_AsFactor            2554138   1   10.335  0.001466 ** 
#   GenPheno               18082139   2   36.582 8.977e-15 ***
#   Sex_AsFactor:GenPheno   6068022   2   12.276 7.922e-06 ***
#   Residuals              66235572 268                       
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Behavior[Behavior$EPM_Distance_Traveled<5000,]

pdf("EPM_Distance_Traveled_noOutlier_RNAseq_SexGen_MaleRef_Color.pdf", width=9.5, height=6)
boxplot(EPM_Distance_Traveled~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine" & Behavior$EPM_Distance_Traveled<5000,], ylab="Distance Traveled in EPM", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("green4","green2","firebrick4", "firebrick2", "grey57", "grey80"))
stripchart( EPM_Distance_Traveled~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine" & Behavior$EPM_Distance_Traveled<5000,], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(EPM_Distance_Traveled~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine"& Behavior$EPM_Distance_Traveled<5000,], contrasts=list(GenPheno=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: EPM_Distance_Traveled
# Sum Sq  Df  F value    Pr(>F)    
# (Intercept)           276095185   1 1327.020 < 2.2e-16 ***
#   Sex_AsFactor            2453347   1   11.792 0.0006895 ***
#   GenPheno               18268156   2   43.902 < 2.2e-16 ***
#   Sex_AsFactor:GenPheno   6368999   2   15.306 5.093e-07 ***
#   Residuals              55551083 267                       
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


pdf("EPM_Time_Immobile__RNAseq_SexGen_MaleRef_Color.pdf", width=9.5, height=6)
boxplot(EPM_Time_Immobile~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], ylab="Time Immobile in EPM", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("green4","green2","firebrick4", "firebrick2", "grey57", "grey80"))
stripchart(EPM_Time_Immobile~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(EPM_Time_Immobile~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], contrasts=list(GenPheno=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: EPM_Time_Immobile
# Sum Sq  Df  F value    Pr(>F)    
# (Intercept)           307656   1 442.8507 < 2.2e-16 ***
#   Sex_AsFactor            6768   1   9.7418  0.001998 ** 
#   GenPheno               54589   2  39.2889 1.088e-15 ***
#   Sex_AsFactor:GenPheno   5047   2   3.6327  0.027756 *  
#   Residuals             186184 268                       
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


pdf("PCA_Index_Days6and.7__RNAseq_SexGen_MaleRef_Color.pdf", width=9.5, height=6)
boxplot(PCA_Index_Days6and.7~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], ylab="PCA_Index_Days6and.7", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("grey57", "grey80"))
stripchart(PCA_Index_Days6and.7~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(PCA_Index_Days6and.7~Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine",], contrasts=list(Sex_AsFactor=contr.sum)), type=3)
# 
# Anova Table (Type III tests)
# 
# Response: PCA_Index_Days6and.7
#              Sum Sq  Df F value    Pr(>F)    
# (Intercept)   4.828   1  18.164 3.078e-05 ***
#   Sex_AsFactor 13.223   1  49.743 2.579e-11 ***
#   Residuals    55.024 207                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


pdf("PCA_Index_Days6and.7_All_SexGen_MaleRef_LearnClass_Yaxis_Title_Color.pdf", width=9.5, height=6)
boxplot(PCA_Index_Days6and.7~Sex_AsFactor*GenPheno, data=Behavior, ylab="Learning Classification", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("grey57", "grey80"))
stripchart(PCA_Index_Days6and.7~Sex_AsFactor*GenPheno, data=Behavior, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()


#########Have GT, IN, and ST only here as desired, but would also like to have ST as the reference level----Sex added here
#This generates same boxplot as above since only F2 rats had ST/GT data and all were used for RNASeq
pdf("Total_LocoScore_LearnClass_BySex_F2.pdf", width=9.5, height=8)
boxplot(Total_LocoScore~Sex_AsFactor*LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], ylab="LocoScore", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(Total_LocoScore~Sex_AsFactor*LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()

#########Have GT, IN, and ST only here as desired, but would also like to have ST as the reference level----Sex and Color added here
#This generates same boxplot as above since only F2 rats had ST/GT data and all were used for RNASeq

pdf("Total_LocoScore_LearnClass_BySex_F2_Color.pdf", width=9.5, height=8)
boxplot(Total_LocoScore~Sex_AsFactor*LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], ylab="LocoScore", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("grey57", "grey80", "grey57", "grey80", "grey57", "grey80"))
stripchart(Total_LocoScore~Sex_AsFactor*LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

BehaviorForPCA<-Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",]

table(BehaviorForPCA$Sex, BehaviorForPCA$LearningClassification)

#   GT IN ST
# F 11 30 60
# M 29 58 21

fisher.test(table(BehaviorForPCA$Sex, BehaviorForPCA$LearningClassification))

# Fisher's Exact Test for Count Data
# 
# data:  table(BehaviorForPCA$Sex, BehaviorForPCA$LearningClassification)
# p-value = 1.254e-08
# alternative hypothesis: two.sided

table(BehaviorForPCA$Sex[BehaviorForPCA$LearningClassification%in%c("ST", "GT")], BehaviorForPCA$LearningClassification[BehaviorForPCA$LearningClassification%in%c("ST", "GT")])
#   GT ST
# F 11 60
# M 29 21

fisher.test(table(BehaviorForPCA$Sex[BehaviorForPCA$LearningClassification%in%c("ST", "GT")], BehaviorForPCA$LearningClassification[BehaviorForPCA$LearningClassification%in%c("ST", "GT")]))

# Fisher's Exact Test for Count Data
# 
# data:  table(BehaviorForPCA$Sex[BehaviorForPCA$LearningClassification %in% c("ST", "GT")], BehaviorForPCA$LearningClassification[BehaviorForPCA$LearningClassification %in% c("ST", "GT")])
# p-value = 1.472e-06
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.05121211 0.33569458
# sample estimates:
# odds ratio 
#  0.1354166 


11/60
#[1] 0.1833333

29/21
#[1] 1.380952

0.1833333/1.380952
#[1] 0.1327586

BehaviorForPCA$LearningClassification_AsFactor<-droplevels(BehaviorForPCA$LearningClassification_AsFactor)

library(car)

Anova(lm(Total_LocoScore~LearningClassification_AsFactor*Sex_AsFactor, data=BehaviorForPCA, contrasts=list(LearningClassification_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: Total_LocoScore
# Sum Sq  Df  F value    Pr(>F)    
# (Intercept)                                  40763976   1 514.3942 < 2.2e-16 ***
#   LearningClassification_AsFactor               3832498   2  24.1808 3.806e-10 ***
#   Sex_AsFactor                                    94577   1   1.1935    0.2759    
# LearningClassification_AsFactor:Sex_AsFactor   247412   2   1.5610    0.2124    
# Residuals                                    16087054 203                       
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


pdf("EPM_Percent_Time_Open_Arm_LearnClass_BySex_F2_Color.pdf", width=9.5, height=8)
boxplot(EPM_Percent_Time_Open_Arm~Sex_AsFactor*LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], ylab="EPM: % Time Open Arm", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("grey57", "grey80", "grey57", "grey80", "grey57", "grey80"))
stripchart(EPM_Percent_Time_Open_Arm~Sex_AsFactor*LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(EPM_Percent_Time_Open_Arm~LearningClassification_AsFactor*Sex_AsFactor, data=BehaviorForPCA, contrasts=list(LearningClassification_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: EPM_Percent_Time_Open_Arm
# Sum Sq  Df  F value  Pr(>F)    
# (Intercept)                                  28059.2   1 282.6061 < 2e-16 ***
#   LearningClassification_AsFactor                490.9   2   2.4720 0.08695 .  
# Sex_AsFactor                                   455.5   1   4.5877 0.03339 *  
#   LearningClassification_AsFactor:Sex_AsFactor   278.2   2   1.4008 0.24878    
# Residuals                                    20155.4 203                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


pdf("EPM_Distance_Traveled_LearnClass_BySex_F2_Color.pdf", width=9.5, height=8)
boxplot(EPM_Distance_Traveled~Sex_AsFactor*LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], ylab="EPM: Distance Traveled", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("grey57", "grey80", "grey57", "grey80", "grey57", "grey80"))
stripchart(EPM_Distance_Traveled~Sex_AsFactor*LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(EPM_Distance_Traveled~LearningClassification_AsFactor*Sex_AsFactor, data=BehaviorForPCA, contrasts=list(LearningClassification_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: EPM_Distance_Traveled
# Sum Sq  Df   F value    Pr(>F)    
# (Intercept)                                  649931211   1 3056.7211 < 2.2e-16 ***
#   LearningClassification_AsFactor                2444181   2    5.7477  0.003732 ** 
#   Sex_AsFactor                                    495042   1    2.3283  0.128601    
# LearningClassification_AsFactor:Sex_AsFactor   1108196   2    2.6060  0.076298 .  
# Residuals                                     43162603 203                        
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


pdf("EPM_Distance_Traveled_noOutlier_LearnClass_BySex_F2_Color.pdf", width=9.5, height=8)
boxplot(EPM_Distance_Traveled~Sex_AsFactor*LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2" & Behavior$EPM_Distance_Traveled<5000,], ylab="EPM: Distance Traveled", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("grey57", "grey80", "grey57", "grey80", "grey57", "grey80"))
stripchart(EPM_Distance_Traveled~Sex_AsFactor*LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2"& Behavior$EPM_Distance_Traveled<5000,], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(EPM_Distance_Traveled~LearningClassification_AsFactor*Sex_AsFactor, data=BehaviorForPCA[BehaviorForPCA$EPM_Distance_Traveled<5000,], contrasts=list(LearningClassification_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: EPM_Distance_Traveled
# Sum Sq  Df   F value    Pr(>F)    
# (Intercept)                                  634660224   1 3978.5148 < 2.2e-16 ***
#   LearningClassification_AsFactor                3196936   2   10.0204   7.092e-05 ***
#   Sex_AsFactor                                    879827   1    5.5154   0.01982 *  
#   LearningClassification_AsFactor:Sex_AsFactor    757553   2    2.3744   0.09566 .  
# Residuals                                     32223423 202                        
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



pdf("EPM_Time_Immobile_LearnClass_BySex_F2_Color.pdf", width=9.5, height=8)
boxplot(EPM_Time_Immobile~Sex_AsFactor*LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], ylab="EPM: Time Immobile", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("grey57", "grey80", "grey57", "grey80", "grey57", "grey80"))
stripchart(EPM_Time_Immobile~Sex_AsFactor*LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(EPM_Time_Immobile~LearningClassification_AsFactor*Sex_AsFactor, data=BehaviorForPCA, contrasts=list(LearningClassification_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: EPM_Time_Immobile
# Sum Sq  Df   F value    Pr(>F)    
# (Intercept)                                  1461012   1 3287.0258 < 2.2e-16 ***
#   LearningClassification_AsFactor                 6400   2    7.1991 0.0009536 ***
#   Sex_AsFactor                                   16318   1   36.7120 6.552e-09 ***
#   LearningClassification_AsFactor:Sex_AsFactor     715   2    0.8047 0.4486462    
# Residuals                                      90229 203                        
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

pdf("EPM_Boli_LearnClass_BySex_F2_Color.pdf", width=9.5, height=8)
boxplot(EPM_Boli~Sex_AsFactor*LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], ylab="EPM: Boli", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("grey57", "grey80", "grey57", "grey80", "grey57", "grey80"))
stripchart(EPM_Boli~Sex_AsFactor*LearningClassification, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(EPM_Boli~LearningClassification_AsFactor*Sex_AsFactor, data=BehaviorForPCA, contrasts=list(LearningClassification_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: EPM_Boli
# Sum Sq  Df F value    Pr(>F)    
# (Intercept)                                   57.61   1 28.2189 2.861e-07 ***
#   LearningClassification_AsFactor                7.45   2  1.8249  0.163891    
# Sex_AsFactor                                  19.05   1  9.3333  0.002555 ** 
#   LearningClassification_AsFactor:Sex_AsFactor   3.77   2  0.9238  0.398688    
# Residuals                                    410.36 201                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


########Have ST as the reference level here, but now the samples without ST/GT data show up here as a group--DID NOT REDO THIS GRAPH HERE
#This generates same boxplot as above since only F2 rats had ST/GT data and all were used for RNASeq

pdf("Total_LocoScore_LearnClass_AsFactor_F2.pdf_Color", width=9.5, height=8)
boxplot(Total_LocoScore~LearningClassification_AsFactor, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], ylab="LocoScore", xlab="", cex.lab=1.5, cex.axis=1.2)
stripchart(Total_LocoScore~LearningClassification_AsFactor, data=Behavior[Behavior$HaveST_GT_Data=="Yes"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()


######################

#Making a quick correlation matrix:

colnames(Behavior)

table(Behavior$Phenotype)
# 
#      HR  LR 
# 801  24  24

table(Behavior$GenPheno)

# F0 HR F0 LR   F1    F2  
# 24    24   261   540 


Behavior_HRLR_forRNASeq<-Behavior[Behavior$RNASeq=="Elaine" & Behavior$Phenotype!="",]
dim(Behavior_HRLR_forRNASeq)
#[1] 24 41
#Note - one of these samples was removed from the actual RNA-Seq analysis due to low RIN

#Making HRs the reference group in a numeric version of the phenotype variable:
Behavior_HRLR_forRNASeq$Phenotype_AsNumeric<-rep(0, 24)

Behavior_HRLR_forRNASeq$Phenotype_AsNumeric[Behavior_HRLR_forRNASeq$Phenotype=="LR"]<-1

write.csv(cor(as.matrix(Behavior_HRLR_forRNASeq[,c(18:24, 42)]), use="pairwise.complete.obs"), "F0_CorrelationMatrix_Behavior_RNASeqAnimals.csv")

Behavior_F2_forRNASeq<-Behavior[Behavior$RNASeq=="Elaine" &Behavior$Study_Generation=="F2",]
dim(Behavior_F2_forRNASeq)
#[1] 250  43

write.csv(cor(as.matrix(Behavior_F2_forRNASeq[,c(18:24, 30:38)]), use="pairwise.complete.obs"), "F2_CorrelationMatrix_Behavior_RNASeqAnimals_All250.csv")

write.csv(cor(as.matrix(Behavior_F2_forRNASeq[Behavior_F2_forRNASeq$EPM_Distance_Traveled<5000,c(18:24, 30:38)]), use="pairwise.complete.obs"), "F2_CorrelationMatrix_Behavior_RNASeqAnimals_249noDistanceTravelOutlier.csv")



#################4/15/2022 Make scatterplot of locosore vs EPM%time in OA for F2s############


pdf("F2_Total_LocoScoreVsEPM_Percent_Time_Open_Arm.pdf", width=6, height=6)
plot(Total_LocoScore~EPM_Percent_Time_Open_Arm, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="LocoScore", pch=18, xlab="EPM %Time in Open Arms")
BestFitLine<-lm(Total_LocoScore~EPM_Percent_Time_Open_Arm, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(Total_LocoScore~EPM_Percent_Time_Open_Arm, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))
# Call:
#   lm(formula = Total_LocoScore ~ EPM_Percent_Time_Open_Arm, data = Behavior[Behavior$RNASeq == 
#                                                                               "Elaine" & Behavior$Study_Generation == "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -533.95 -228.56  -66.93  191.98  954.65 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                420.654     30.020  14.013  < 2e-16 ***
#   EPM_Percent_Time_Open_Arm    8.519      1.738   4.901 1.72e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 289.1 on 248 degrees of freedom
# Multiple R-squared:  0.0883,	Adjusted R-squared:  0.08462 
# F-statistic: 24.02 on 1 and 248 DF,  p-value: 1.724e-06

sqrt( 0.0883)
#[1] 0.2971532


pdf("F2_Total_LocoScoreVsEPM_Percent_Time_Open_Arm_ColSex.pdf", width=6, height=6)
plot(Total_LocoScore~EPM_Percent_Time_Open_Arm, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="LocoScore", pch=18, xlab="EPM %Time in Open Arms",col=as.numeric(Behavior$Sex_AsFactor[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"]))
BestFitLine<-lm(Total_LocoScore~EPM_Percent_Time_Open_Arm, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=1, lwd=3)
dev.off()

summary.lm(lm(Total_LocoScore~EPM_Percent_Time_Open_Arm*Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))

# Call:
#   lm(formula = Total_LocoScore ~ EPM_Percent_Time_Open_Arm * Sex_AsFactor, 
#      data = Behavior[Behavior$RNASeq == "Elaine" & Behavior$Study_Generation == 
#                        "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -538.67 -226.52  -70.26  194.55  938.27 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                             409.0497    39.2517  10.421  < 2e-16 ***
#   EPM_Percent_Time_Open_Arm                 8.4167     2.7262   3.087  0.00225 ** 
#   Sex_AsFactorF                            33.6443    61.9293   0.543  0.58744    
# EPM_Percent_Time_Open_Arm:Sex_AsFactorF  -0.4604     3.6506  -0.126  0.89974    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 290 on 246 degrees of freedom
# Multiple R-squared:  0.09028,	Adjusted R-squared:  0.07919 
# F-statistic: 8.138 on 3 and 246 DF,  p-value: 3.458e-05


colnames(Behavior)

pdf("F2_Total_LocoScoreVsEPM_Time_Immobile.pdf", width=6, height=6)
plot(Total_LocoScore~EPM_Time_Immobile, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="LocoScore", pch=18, xlab="EPM Time Immobile")
BestFitLine<-lm(Total_LocoScore~EPM_Time_Immobile, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(Total_LocoScore~EPM_Time_Immobile, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))
#
# Call:
#   lm(formula = Total_LocoScore ~ EPM_Time_Immobile, data = Behavior[Behavior$RNASeq == 
#                                                                       "Elaine" & Behavior$Study_Generation == "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -647.73 -200.78  -38.26  180.46  831.77 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       1052.2523    67.7067  15.541  < 2e-16 ***
#   EPM_Time_Immobile   -5.2023     0.6618  -7.861 1.17e-13 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 270.9 on 248 degrees of freedom
# Multiple R-squared:  0.1995,	Adjusted R-squared:  0.1962 
# F-statistic:  61.8 on 1 and 248 DF,  p-value: 1.165e-13

sqrt(0.1995)
#[1] 0.4466542


pdf("F2_Total_LocoScoreVsEPM_EPM_Time_Immobile_ColSex.pdf", width=6, height=6)
plot(Total_LocoScore~EPM_Time_Immobile, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="LocoScore", pch=18, xlab="EPM Time Immobile",col=as.numeric(Behavior$Sex_AsFactor[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"]))
BestFitLine<-lm(Total_LocoScore~EPM_Time_Immobile, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=1, lwd=3)
dev.off()

summary.lm(lm(Total_LocoScore~EPM_Time_Immobile*Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))

# 
# Call:
#   lm(formula = Total_LocoScore ~ EPM_Time_Immobile * Sex_AsFactor, 
#      data = Behavior[Behavior$RNASeq == "Elaine" & Behavior$Study_Generation == 
#                        "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -681.21 -196.67  -47.66  163.34  847.95 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                     1093.0424   105.8939  10.322  < 2e-16 ***
#   EPM_Time_Immobile                 -5.3054     0.9226  -5.751 2.62e-08 ***
#   Sex_AsFactorF                    105.9621   157.2823   0.674    0.501    
# EPM_Time_Immobile:Sex_AsFactorF   -1.9393     1.6106  -1.204    0.230    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 269.1 on 246 degrees of freedom
# Multiple R-squared:  0.2165,	Adjusted R-squared:  0.207 
# F-statistic: 22.66 on 3 and 246 DF,  p-value: 5.449e-13










#########################5/6/2022 Elaine#######################
colnames(Behavior)

pdf("F2_EPM_Percent_Time_Open_ArmVsEPM_Distance_Traveled.pdf", width=6, height=6)
plot(EPM_Percent_Time_Open_Arm~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM % Time Open Arms", pch=18, xlab="EPM Distance Traveled")
BestFitLine<-lm(EPM_Percent_Time_Open_Arm~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(EPM_Percent_Time_Open_Arm~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))

# Call:
#   lm(formula = EPM_Percent_Time_Open_Arm ~ EPM_Distance_Traveled, 
#      data = Behavior[Behavior$RNASeq == "Elaine" & Behavior$Study_Generation == 
#                        "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -32.375  -6.045  -1.015   5.230  27.001 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           -8.478041   2.554434  -3.319  0.00104 ** 
#   EPM_Distance_Traveled  0.010659   0.001196   8.915  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 9.192 on 248 degrees of freedom
# Multiple R-squared:  0.2427,	Adjusted R-squared:  0.2396 
# F-statistic: 79.48 on 1 and 248 DF,  p-value: < 2.2e-16

sqrt( 0.2427)
#[1] 0.4926459

pdf("F2_EPM_Percent_Time_Open_ArmVsEPM_Distance_Traveled_noOutlier.pdf", width=6, height=6)
plot(EPM_Percent_Time_Open_Arm~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,], ylab="EPM % Time Open Arms", pch=18, xlab="EPM Distance Traveled")
BestFitLine<-lm(EPM_Percent_Time_Open_Arm~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(EPM_Percent_Time_Open_Arm~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,]))
# Call:
#   lm(formula = EPM_Percent_Time_Open_Arm ~ EPM_Distance_Traveled, 
#      data = Behavior[Behavior$RNASeq == "Elaine" & Behavior$Study_Generation == 
#                        "F2" & Behavior$EPM_Distance_Traveled < 5000, ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -27.2969  -5.7035  -0.9877   5.1423  26.7039 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           -12.616037   2.690401  -4.689 4.54e-06 ***
#   EPM_Distance_Traveled   0.012723   0.001272  10.003  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 8.928 on 247 degrees of freedom
# Multiple R-squared:  0.2883,	Adjusted R-squared:  0.2854 
# F-statistic: 100.1 on 1 and 247 DF,  p-value: < 2.2e-16

sqrt(0.2883)
#[1] 0.5369358

pdf("F2_EPM_Percent_Time_Open_ArmVsEPM_Distance_Traveled_ColSex.pdf", width=6, height=6)
plot(EPM_Percent_Time_Open_Arm~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM % Time Open Arms", pch=18, xlab="EPM Distance Traveled",col=as.numeric(Behavior$Sex_AsFactor[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"]))
BestFitLine<-lm(EPM_Percent_Time_Open_Arm~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=1, lwd=3)
dev.off()

summary.lm(lm(EPM_Percent_Time_Open_Arm~EPM_Distance_Traveled*Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))
# Call:
#   lm(formula = EPM_Percent_Time_Open_Arm ~ EPM_Distance_Traveled * 
#        Sex_AsFactor, data = Behavior[Behavior$RNASeq == "Elaine" & 
#                                        Behavior$Study_Generation == "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -20.992  -6.278  -1.068   5.035  29.129 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         -4.105204   3.162391  -1.298   0.1955    
# EPM_Distance_Traveled                0.007635   0.001566   4.876 1.94e-06 ***
#   Sex_AsFactorF                       -7.981763   5.349135  -1.492   0.1369    
# EPM_Distance_Traveled:Sex_AsFactorF  0.005355   0.002478   2.161   0.0316 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 9.003 on 246 degrees of freedom
# Multiple R-squared:  0.2792,	Adjusted R-squared:  0.2704 
# F-statistic: 31.77 on 3 and 246 DF,  p-value: < 2.2e-16


pdf("F2_EPM_Percent_Time_Open_ArmVsEPM_Distance_TravelednoOutlier_ColSex.pdf", width=6, height=6)
plot(EPM_Percent_Time_Open_Arm~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,], ylab="EPM % Time Open Arms", pch=18, xlab="EPM Distance Traveled",col=as.numeric(Behavior$Sex_AsFactor[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"]))
BestFitLine<-lm(EPM_Percent_Time_Open_Arm~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,])
abline(BestFitLine, col=1, lwd=3)
dev.off()

summary.lm(lm(EPM_Percent_Time_Open_Arm~EPM_Distance_Traveled*Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,]))

# Call:
#   lm(formula = EPM_Percent_Time_Open_Arm ~ EPM_Distance_Traveled * 
#        Sex_AsFactor, data = Behavior[Behavior$RNASeq == "Elaine" & 
#                                        Behavior$Study_Generation == "F2" & Behavior$EPM_Distance_Traveled < 
#                                        5000, ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -23.5957  -5.5445  -0.8155   4.9631  28.2761 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         -9.867076   3.700390  -2.666  0.00818 ** 
#   EPM_Distance_Traveled                0.010713   0.001875   5.713 3.21e-08 ***
#   Sex_AsFactorF                       -2.219891   5.636162  -0.394  0.69402    
# EPM_Distance_Traveled:Sex_AsFactorF  0.002277   0.002664   0.855  0.39351    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 8.872 on 245 degrees of freedom
# Multiple R-squared:  0.3029,	Adjusted R-squared:  0.2944 
# F-statistic: 35.49 on 3 and 245 DF,  p-value: < 2.2e-16



pdf("F2_EPM_Percent_Time_Open_ArmVsEPM_Time_Immobile.pdf", width=6, height=6)
plot(EPM_Percent_Time_Open_Arm~EPM_Time_Immobile, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM % Time Open Arms", pch=18, xlab="EPM Time Immobile")
BestFitLine<-lm(EPM_Percent_Time_Open_Arm~EPM_Time_Immobile, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(EPM_Percent_Time_Open_Arm~EPM_Time_Immobile, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))
#
# Call:
#   lm(formula = EPM_Percent_Time_Open_Arm ~ EPM_Time_Immobile, data = Behavior[Behavior$RNASeq == 
#                                                                                 "Elaine" & Behavior$Study_Generation == "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -17.668  -6.776  -1.465   5.257  26.779 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       31.86531    2.35489  13.532  < 2e-16 ***
#   EPM_Time_Immobile -0.18356    0.02302  -7.975 5.61e-14 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 9.423 on 248 degrees of freedom
# Multiple R-squared:  0.2041,	Adjusted R-squared:  0.2009 
# F-statistic:  63.6 on 1 and 248 DF,  p-value: 5.613e-14

sqrt(0.2041)
#[1] 0.4517743

pdf("F2_EPM_Percent_Time_Open_ArmVsEPM_Time_Immobile_ColSex.pdf", width=6, height=6)
plot(EPM_Percent_Time_Open_Arm~EPM_Time_Immobile, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM % Time Open Arms", pch=18, xlab="EPM Time Immobile",col=as.numeric(Behavior$Sex_AsFactor[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"]))
BestFitLine<-lm(EPM_Percent_Time_Open_Arm~EPM_Time_Immobile, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=1, lwd=3)
dev.off()

summary.lm(lm(EPM_Percent_Time_Open_Arm~EPM_Time_Immobile*Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))
# Call:
# lm(formula = EPM_Percent_Time_Open_Arm ~ EPM_Time_Immobile * 
#      Sex_AsFactor, data = Behavior[Behavior$RNASeq == "Elaine" & 
#                                      Behavior$Study_Generation == "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -18.803  -6.554  -1.573   5.251  27.770 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                     26.92044    3.70021   7.275 4.61e-12 ***
#   EPM_Time_Immobile               -0.14416    0.03224  -4.472 1.19e-05 ***
#   Sex_AsFactorF                    8.85207    5.49585   1.611    0.109    
# EPM_Time_Immobile:Sex_AsFactorF -0.07847    0.05628  -1.394    0.164    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 9.403 on 246 degrees of freedom
# Multiple R-squared:  0.2138,	Adjusted R-squared:  0.2042 
# F-statistic:  22.3 on 3 and 246 DF,  p-value: 8.319e-13






#########################
pdf("F2_Total_LocoScoreVsEPM_DistanceTraveled.pdf", width=6, height=6)
plot(Total_LocoScore~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="LocoScore", pch=18, xlab="EPM Distance Traveled")
BestFitLine<-lm(Total_LocoScore~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(Total_LocoScore~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))

# Call:
#   lm(formula = Total_LocoScore ~ EPM_Distance_Traveled, data = Behavior[Behavior$RNASeq == 
#                                                                           "Elaine" & Behavior$Study_Generation == "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -733.56 -205.48  -30.41  164.25  792.30 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           -33.34628   75.47503  -0.442    0.659    
# EPM_Distance_Traveled   0.27433    0.03533   7.765 2.15e-13 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 271.6 on 248 degrees of freedom
# Multiple R-squared:  0.1956,	Adjusted R-squared:  0.1923 
# F-statistic: 60.29 on 1 and 248 DF,  p-value: 2.148e-13

sqrt(0.1956)
#[1] 0.4422669

pdf("F2_Total_LocoScoreVsEPM_DistanceTravelednoOutlier.pdf", width=6, height=6)
plot(Total_LocoScore~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,], ylab="LocoScore", pch=18, xlab="EPM Distance Traveled")
BestFitLine<-lm(Total_LocoScore~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(Total_LocoScore~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,]))

#Call:
 # lm(formula = Total_LocoScore ~ EPM_Distance_Traveled, data = Behavior[Behavior$RNASeq == 
                                                                          # "Elaine" & Behavior$Study_Generation == "F2" & Behavior$EPM_Distance_Traveled < 
                                                                          # 5000, ])

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -554.99 -191.64  -37.81  176.16  790.58 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           -127.10664   80.53511  -1.578    0.116    
# EPM_Distance_Traveled    0.32109    0.03808   8.433 2.85e-15 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 267.3 on 247 degrees of freedom
# Multiple R-squared:  0.2235,	Adjusted R-squared:  0.2204 
# F-statistic: 71.11 on 1 and 247 DF,  p-value: 2.846e-15

sqrt(0.2235)
#[1] 0.4727579

pdf("F2_Total_LocoScoreVsEPM_DistanceTraveled_ColSex.pdf", width=6, height=6)
plot(Total_LocoScore~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="LocoScore", pch=18, xlab="EPM Distance Traveled",col=as.numeric(Behavior$Sex_AsFactor[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"]))
BestFitLine<-lm(Total_LocoScore~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=1, lwd=3)
dev.off()

summary.lm(lm(Total_LocoScore~EPM_Distance_Traveled*Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))

# Call:
#   lm(formula = Total_LocoScore ~ EPM_Distance_Traveled * Sex_AsFactor, 
#      data = Behavior[Behavior$RNASeq == "Elaine" & Behavior$Study_Generation == 
#                        "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -573.41 -206.87  -43.36  167.77  799.04 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           58.80335   95.29299   0.617    0.538    
# EPM_Distance_Traveled                  0.22589    0.04718   4.788 2.91e-06 ***
#   Sex_AsFactorF                       -242.64372  161.18661  -1.505    0.134    
# EPM_Distance_Traveled:Sex_AsFactorF    0.11773    0.07466   1.577    0.116    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 271.3 on 246 degrees of freedom
# Multiple R-squared:  0.2037,	Adjusted R-squared:  0.194 
# F-statistic: 20.98 on 3 and 246 DF,  p-value: 3.914e-12


pdf("F2_Total_LocoScoreVsEPM_DistanceTravelednoOutlier_ColSex.pdf", width=6, height=6)
plot(Total_LocoScore~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,], ylab="LocoScore", pch=18, xlab="EPM Distance Traveled",col=as.numeric(Behavior$Sex_AsFactor[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"]))
BestFitLine<-lm(Total_LocoScore~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,])
abline(BestFitLine, col=1, lwd=3)
dev.off()

summary.lm(lm(Total_LocoScore~EPM_Distance_Traveled*Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,]))

# Call:
#   lm(formula = Total_LocoScore ~ EPM_Distance_Traveled * Sex_AsFactor, 
#      data = Behavior[Behavior$RNASeq == "Elaine" & Behavior$Study_Generation == 
#                        "F2" & Behavior$EPM_Distance_Traveled < 5000, ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -554.24 -198.97  -43.69  166.62  784.78 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         -98.58280  111.84239  -0.881    0.379    
# EPM_Distance_Traveled                 0.30996    0.05668   5.469 1.12e-07 ***
#   Sex_AsFactorF                       -85.25757  170.35012  -0.500    0.617    
# EPM_Distance_Traveled:Sex_AsFactorF   0.03365    0.08052   0.418    0.676    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 268.2 on 245 degrees of freedom
# Multiple R-squared:  0.2247,	Adjusted R-squared:  0.2152 
# F-statistic: 23.67 on 3 and 245 DF,  p-value: 1.737e-13


pdf("F2_EPMTimeImmobileVsEPM_DistanceTraveled.pdf", width=6, height=6)
plot(EPM_Time_Immobile~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM Time Immobile", pch=18, xlab="EPM Distance Traveled")
BestFitLine<-lm(EPM_Time_Immobile~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(EPM_Time_Immobile~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))

# Call:
#   lm(formula = EPM_Time_Immobile ~ EPM_Distance_Traveled, data = Behavior[Behavior$RNASeq == 
#                                                                             "Elaine" & Behavior$Study_Generation == "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -40.228 -10.911  -1.561   9.115  99.530 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           182.855552   4.719715   38.74   <2e-16 ***
#   EPM_Distance_Traveled  -0.040319   0.002209  -18.25   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 16.98 on 248 degrees of freedom
# Multiple R-squared:  0.5732,	Adjusted R-squared:  0.5715 
# F-statistic: 333.1 on 1 and 248 DF,  p-value: < 2.2e-16

sqrt(0.5732)
#[1] 0.7570997


pdf("F2_EPMTimeImmobileVsEPM_DistanceTravelednoOutlier.pdf", width=6, height=6)
plot(EPM_Time_Immobile~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,], ylab="EPM Time Immobile", pch=18, xlab="EPM Distance Traveled")
BestFitLine<-lm(EPM_Time_Immobile~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(EPM_Time_Immobile~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,]))

# Call:
#   lm(formula = EPM_Time_Immobile ~ EPM_Distance_Traveled, data = Behavior[Behavior$RNASeq == 
#                                                                             "Elaine" & Behavior$Study_Generation == "F2" & Behavior$EPM_Distance_Traveled < 
#                                                                             5000, ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -41.468 -10.775  -0.779  10.189  44.910 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           195.576930   4.680681   41.78   <2e-16 ***
#   EPM_Distance_Traveled  -0.046665   0.002213  -21.09   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 15.53 on 247 degrees of freedom
# Multiple R-squared:  0.6429,	Adjusted R-squared:  0.6414 
# F-statistic: 444.6 on 1 and 247 DF,  p-value: < 2.2e-16

sqrt(0.6429)
#[1] 0.8018105


pdf("F2_EPMTimeImmobileVsEPM_DistanceTraveled_ColSex.pdf", width=6, height=6)
plot(EPM_Time_Immobile~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM Time Immobile", pch=18, xlab="EPM Distance Traveled",col=as.numeric(Behavior$Sex_AsFactor[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"]))
BestFitLine<-lm(EPM_Time_Immobile~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=1, lwd=3)
dev.off()

summary.lm(lm(EPM_Time_Immobile~EPM_Distance_Traveled*Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))

# Call:
#   lm(formula = EPM_Time_Immobile ~ EPM_Distance_Traveled * Sex_AsFactor, 
#      data = Behavior[Behavior$RNASeq == "Elaine" & Behavior$Study_Generation == 
#                        "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -36.244  -8.797  -0.956   8.166  82.417 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         184.868586   5.277030  35.033   <2e-16 ***
#   EPM_Distance_Traveled                -0.037420   0.002613 -14.323   <2e-16 ***
#   Sex_AsFactorF                       -24.549414   8.926015  -2.750   0.0064 ** 
#   EPM_Distance_Traveled:Sex_AsFactorF   0.003834   0.004134   0.927   0.3547    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 15.02 on 246 degrees of freedom
# Multiple R-squared:  0.6687,	Adjusted R-squared:  0.6646 
# F-statistic: 165.5 on 3 and 246 DF,  p-value: < 2.2e-16

pdf("F2_EPMTimeImmobileVsEPM_DistanceTravelednoOutlier_ColSex.pdf", width=6, height=6)
plot(EPM_Time_Immobile~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,], ylab="EPM Time Immobile", pch=18, xlab="EPM Distance Traveled",col=as.numeric(Behavior$Sex_AsFactor[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"]))
BestFitLine<-lm(EPM_Time_Immobile~EPM_Distance_Traveled, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,])
abline(BestFitLine, col=1, lwd=3)
dev.off()

summary.lm(lm(EPM_Time_Immobile~EPM_Distance_Traveled*Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,]))

# Call:
#   lm(formula = EPM_Time_Immobile ~ EPM_Distance_Traveled * Sex_AsFactor, 
#      data = Behavior[Behavior$RNASeq == "Elaine" & Behavior$Study_Generation == 
#                        "F2" & Behavior$EPM_Distance_Traveled < 5000, ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -36.318  -8.700  -0.895   7.966  36.978 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         207.489955   5.678275  36.541  < 2e-16 ***
#   EPM_Distance_Traveled                -0.049505   0.002878 -17.203  < 2e-16 ***
#   Sex_AsFactorF                       -47.170784   8.648732  -5.454  1.2e-07 ***
#   EPM_Distance_Traveled:Sex_AsFactorF   0.015918   0.004088   3.894 0.000127 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 13.61 on 245 degrees of freedom
# Multiple R-squared:  0.7279,	Adjusted R-squared:  0.7246 
# F-statistic: 218.5 on 3 and 245 DF,  p-value: < 2.2e-16


pdf("F2_Total_LocoScoreVsBoli.pdf", width=6, height=6)
plot(Total_LocoScore~EPM_Boli, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="LocoScore", pch=18, xlab="EPM Boli")
BestFitLine<-lm(Total_LocoScore~EPM_Boli, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(Total_LocoScore~EPM_Boli, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))

# Call:
#   lm(formula = Total_LocoScore ~ EPM_Boli, data = Behavior[Behavior$RNASeq == 
#                                                              "Elaine" & Behavior$Study_Generation == "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -502.25 -243.85  -54.25  212.75  893.75 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   567.25      20.68  27.432  < 2e-16 ***
#   EPM_Boli      -40.90      11.91  -3.435 0.000696 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 296.7 on 245 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.04594,	Adjusted R-squared:  0.04205 
# F-statistic:  11.8 on 1 and 245 DF,  p-value: 0.0006963

sqrt(0.04594)
#[1] 0.2143362

pdf("F2_Total_LocoScoreVsBoli_ColSex.pdf", width=6, height=6)
plot(Total_LocoScore~EPM_Boli, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="LocoScore", pch=18, xlab="EPM Boli", col=as.numeric(Behavior$Sex_AsFactor[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"]))
BestFitLine<-lm(Total_LocoScore~EPM_Boli, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=1, lwd=3)
dev.off()

summary.lm(lm(Total_LocoScore~EPM_Boli*Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))

# Call:
#   lm(formula = Total_LocoScore ~ EPM_Boli * Sex_AsFactor, data = Behavior[Behavior$RNASeq == 
#                                                                             "Elaine" & Behavior$Study_Generation == "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -483.39 -252.89  -47.39  202.11  878.61 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             546.712     31.275  17.481  < 2e-16 ***
#   EPM_Boli                -37.658     13.365  -2.818  0.00524 ** 
#   Sex_AsFactorF            35.675     41.858   0.852  0.39490    
# EPM_Boli:Sex_AsFactorF    4.926     40.262   0.122  0.90272    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 297.4 on 243 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.04939,	Adjusted R-squared:  0.03765 
# F-statistic: 4.208 on 3 and 243 DF,  p-value: 0.006327


pdf("F2_EPM_Percent_Time_Open_ArmVsBoli.pdf", width=6, height=6)
plot(EPM_Percent_Time_Open_Arm~EPM_Boli, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM: % Time Open Arm", pch=18, xlab="EPM Boli")
BestFitLine<-lm(EPM_Percent_Time_Open_Arm~EPM_Boli, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(EPM_Percent_Time_Open_Arm~EPM_Boli, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))

# Call:
#   lm(formula = EPM_Percent_Time_Open_Arm ~ EPM_Boli, data = Behavior[Behavior$RNASeq == 
#                                                                        "Elaine" & Behavior$Study_Generation == "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -14.752  -8.353  -1.963   7.388  31.677 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  14.7524     0.7204  20.477  < 2e-16 ***
#   EPM_Boli     -1.3996     0.4148  -3.374 0.000861 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 10.34 on 245 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.0444,	Adjusted R-squared:  0.0405 
# F-statistic: 11.38 on 1 and 245 DF,  p-value: 0.0008614


pdf("F2_EPM_Percent_Time_Open_ArmVsBoli_ColSex.pdf", width=6, height=6)
plot(EPM_Percent_Time_Open_Arm~EPM_Boli, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM: % Time Open Arm", pch=18, xlab="EPM Boli", col=as.numeric(Behavior$Sex_AsFactor[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"]))
BestFitLine<-lm(EPM_Percent_Time_Open_Arm~EPM_Boli, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=1, lwd=3)
dev.off()

summary.lm(lm(EPM_Percent_Time_Open_Arm~EPM_Boli*Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))

# Call:
#   lm(formula = EPM_Percent_Time_Open_Arm ~ EPM_Boli * Sex_AsFactor, 
#      data = Behavior[Behavior$RNASeq == "Elaine" & Behavior$Study_Generation == 
#                        "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -16.772  -8.012  -2.042   6.268  29.423 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             12.0070     1.0636  11.289  < 2e-16 ***
#   EPM_Boli                -0.9347     0.4545  -2.056 0.040808 *  
#   Sex_AsFactorF            4.7954     1.4235   3.369 0.000878 ***
#   EPM_Boli:Sex_AsFactorF   0.2885     1.3692   0.211 0.833280    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 10.11 on 243 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.09274,	Adjusted R-squared:  0.08154 
# F-statistic:  8.28 on 3 and 243 DF,  p-value: 2.888e-05


pdf("F2_EPM_Percent_Time_Open_ArmVsPavCA.pdf", width=6, height=6)
plot(EPM_Percent_Time_Open_Arm~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM: % Time Open Arm", pch=18, xlab="PavCA Index: Days 6 & 7")
BestFitLine<-lm(EPM_Percent_Time_Open_Arm~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(EPM_Percent_Time_Open_Arm~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))

# Call:
#   lm(formula = EPM_Percent_Time_Open_Arm ~ PCA_Index_Days6and.7, 
#      data = Behavior[Behavior$RNASeq == "Elaine" & Behavior$Study_Generation == 
#                        "F2", ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -16.8792  -7.5112  -0.9176   7.4509  24.5504 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           13.0725     0.7112  18.380  < 2e-16 ***
#   PCA_Index_Days6and.7   5.0362     1.2071   4.172 4.44e-05 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 9.972 on 207 degrees of freedom
# (41 observations deleted due to missingness)
# Multiple R-squared:  0.07757,	Adjusted R-squared:  0.07312 
#F-statistic: 17.41 on 1 and 207 DF,  p-value: 4.436e-05

sqrt( 0.07757)
#[1] 0.2785139

pdf("F2_EPM_Percent_Time_Open_ArmVsPavCA_ColSex.pdf", width=6, height=6)
plot(EPM_Percent_Time_Open_Arm~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM: % Time Open Arm", pch=18, xlab="PavCA Index: Days 6 & 7", col=as.numeric(Behavior$Sex_AsFactor[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"]))
BestFitLine<-lm(EPM_Percent_Time_Open_Arm~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=1, lwd=3)
dev.off()

summary.lm(lm(EPM_Percent_Time_Open_Arm~PCA_Index_Days6and.7*Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))
# Call:
#   lm(formula = EPM_Percent_Time_Open_Arm ~ PCA_Index_Days6and.7 * 
#        Sex_AsFactor, data = Behavior[Behavior$RNASeq == "Elaine" & 
#                                        Behavior$Study_Generation == "F2", ])
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -17.7356  -8.3011  -0.6977   6.8253  26.5764 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         11.5801     0.9683  11.959   <2e-16 ***
#   PCA_Index_Days6and.7                 3.9830     1.8050   2.207   0.0284 *  
#   Sex_AsFactorF                        3.6833     1.5946   2.310   0.0219 *  
#   PCA_Index_Days6and.7:Sex_AsFactorF  -0.6984     2.6770  -0.261   0.7944    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 9.888 on 205 degrees of freedom
# (41 observations deleted due to missingness)
# Multiple R-squared:  0.1018,	Adjusted R-squared:  0.08864 
#F-statistic: 7.744 on 3 and 205 DF,  p-value: 6.339e-05


pdf("F2_EPM_Distance_Traveled_ArmVsPavCA.pdf", width=6, height=6)
plot(EPM_Distance_Traveled~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM: Distance Traveled", pch=18, xlab="PavCA Index: Days 6 & 7")
BestFitLine<-lm(EPM_Distance_Traveled~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(EPM_Distance_Traveled~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))

# Call:
#   lm(formula = EPM_Distance_Traveled ~ PCA_Index_Days6and.7, data = Behavior[Behavior$RNASeq == 
#                                                                                "Elaine" & Behavior$Study_Generation == "F2", ])
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -812.4 -300.7  -64.5  242.6 3371.5 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           2064.21      33.23   62.12  < 2e-16 ***
#   PCA_Index_Days6and.7   257.14      56.39    4.56 8.76e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 465.9 on 207 degrees of freedom
# (41 observations deleted due to missingness)
# Multiple R-squared:  0.09128,	Adjusted R-squared:  0.08689 
# F-statistic: 20.79 on 1 and 207 DF,  p-value: 8.755e-06


pdf("F2_EPM_Distance_Traveled_Open_ArmVsPavCA_ColSex.pdf", width=6, height=6)
plot(EPM_Distance_Traveled~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM: Distance Traveled", pch=18, xlab="PavCA Index: Days 6 & 7", col=as.numeric(Behavior$Sex_AsFactor[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"]))
BestFitLine<-lm(EPM_Distance_Traveled~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=1, lwd=3)
dev.off()

summary.lm(lm(EPM_Distance_Traveled~PCA_Index_Days6and.7*Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))

# Call:
#   lm(formula = EPM_Distance_Traveled ~ PCA_Index_Days6and.7 * Sex_AsFactor, 
#      data = Behavior[Behavior$RNASeq == "Elaine" & Behavior$Study_Generation == 
#                        "F2", ])
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -869.6 -308.8  -59.2  248.5 3361.6 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         1995.90      45.31  44.050   <2e-16 ***
#   PCA_Index_Days6and.7                 168.57      84.46   1.996   0.0473 *  
#   Sex_AsFactorF                        141.20      74.62   1.892   0.0599 .  
# PCA_Index_Days6and.7:Sex_AsFactorF    65.58     125.27   0.523   0.6012    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 462.7 on 205 degrees of freedom
# (41 observations deleted due to missingness)
# Multiple R-squared:  0.1123,	Adjusted R-squared:  0.09931 
#F-statistic: 8.645 on 3 and 205 DF,  p-value: 1.984e-05



pdf("F2_EPM_Distance_Traveled_noOutierVsPavCA.pdf", width=6, height=6)
plot(EPM_Distance_Traveled~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,], ylab="EPM: Distance Traveled", pch=18, xlab="PavCA Index: Days 6 & 7")
BestFitLine<-lm(EPM_Distance_Traveled~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(EPM_Distance_Traveled~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,]))
# Call:
#   lm(formula = EPM_Distance_Traveled ~ PCA_Index_Days6and.7, data = Behavior[Behavior$RNASeq == 
#                                                                                "Elaine" & Behavior$Study_Generation == "F2" & Behavior$EPM_Distance_Traveled < 
#                                                                                5000, ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -823.84 -273.23  -38.29  269.50 1240.29 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           2040.31      28.83  70.781  < 2e-16 ***
#   PCA_Index_Days6and.7   308.86      49.07   6.294 1.83e-09 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 402.2 on 206 degrees of freedom
# (41 observations deleted due to missingness)
# Multiple R-squared:  0.1613,	Adjusted R-squared:  0.1572 
# F-statistic: 39.62 on 1 and 206 DF,  p-value: 1.826e-09

sqrt(0.1613)
#[1] 0.4016217

pdf("F2_EPM_Distance_Traveled_noOutlierVsPavCA_ColSex.pdf", width=6, height=6)
plot(EPM_Distance_Traveled~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,], ylab="EPM: Distance Traveled", pch=18, xlab="PavCA Index: Days 6 & 7", col=as.numeric(Behavior$Sex_AsFactor[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"]))
BestFitLine<-lm(EPM_Distance_Traveled~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(EPM_Distance_Traveled~PCA_Index_Days6and.7*Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"&Behavior$EPM_Distance_Traveled<5000,]))
# Call:
#   lm(formula = EPM_Distance_Traveled ~ PCA_Index_Days6and.7 * Sex_AsFactor, 
#      data = Behavior[Behavior$RNASeq == "Elaine" & Behavior$Study_Generation == 
#                        "F2" & Behavior$EPM_Distance_Traveled < 5000, ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -869.62 -305.35  -42.39  249.11 1270.95 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         1972.82      39.02  50.554  < 2e-16 ***
#   PCA_Index_Days6and.7                 258.90      73.33   3.531 0.000512 ***
#   Sex_AsFactorF                        164.28      64.17   2.560 0.011185 *  
#   PCA_Index_Days6and.7:Sex_AsFactorF   -24.76     108.15  -0.229 0.819109    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 397.6 on 204 degrees of freedom
# (41 observations deleted due to missingness)
# Multiple R-squared:  0.1886,	Adjusted R-squared:  0.1767 
# F-statistic: 15.81 on 3 and 204 DF,  p-value: 2.799e-09





pdf("F2_Total_LocoScore_VsPavCA.pdf", width=6, height=6)
plot(Total_LocoScore~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="LocoScore", pch=18, xlab="PavCA Index: Days 6 & 7")
BestFitLine<-lm(Total_LocoScore~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(Total_LocoScore~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))
# 
# Call:
#   lm(formula = Total_LocoScore ~ PCA_Index_Days6and.7, data = Behavior[Behavior$RNASeq == 
#                                                                          "Elaine" & Behavior$Study_Generation == "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -577.62 -212.39  -19.86  195.28 1046.82 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            506.40      19.84  25.528  < 2e-16 ***
#   PCA_Index_Days6and.7   249.36      33.67   7.407 3.22e-12 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 278.1 on 207 degrees of freedom
# (41 observations deleted due to missingness)
# Multiple R-squared:  0.2095,	Adjusted R-squared:  0.2057 
# F-statistic: 54.86 on 1 and 207 DF,  p-value: 3.224e-12

sqrt(0.2095)
#[1] 0.4577117

pdf("F2_Total_LocoScoreVsPavCA_ColSex.pdf", width=6, height=6)
plot(Total_LocoScore~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="LocoScore", pch=18, xlab="PavCA Index: Days 6 & 7", col=as.numeric(Behavior$Sex_AsFactor[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"]))
BestFitLine<-lm(Total_LocoScore~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=1, lwd=3)
dev.off()

summary.lm(lm(Total_LocoScore~PCA_Index_Days6and.7*Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))

# Call:
#   lm(formula = Total_LocoScore ~ PCA_Index_Days6and.7 * Sex_AsFactor, 
#      data = Behavior[Behavior$RNASeq == "Elaine" & Behavior$Study_Generation == 
#                        "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -550.30 -199.65  -32.97  169.57 1073.16 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          536.35      27.19  19.726  < 2e-16 ***
#   PCA_Index_Days6and.7                 295.97      50.68   5.840 2.03e-08 ***
#   Sex_AsFactorF                        -56.64      44.78  -1.265    0.207    
# PCA_Index_Days6and.7:Sex_AsFactorF   -47.54      75.17  -0.632    0.528    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 277.7 on 205 degrees of freedom
# (41 observations deleted due to missingness)
# Multiple R-squared:  0.2198,	Adjusted R-squared:  0.2084 
# F-statistic: 19.25 on 3 and 205 DF,  p-value: 4.895e-11




pdf("F2_EPM_Time_Immobile_VsPavCA.pdf", width=6, height=6)
plot(EPM_Time_Immobile~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM Time Immobile", pch=18, xlab="PavCA Index: Days 6 & 7")
BestFitLine<-lm(EPM_Time_Immobile~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(lm(EPM_Time_Immobile~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))
# 
# Call:
#   lm(formula = EPM_Time_Immobile ~ PCA_Index_Days6and.7, data = Behavior[Behavior$RNASeq == 
#                                                                            "Elaine" & Behavior$Study_Generation == "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -51.665 -16.191  -2.156  13.058  67.027 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           100.377      1.611  62.312  < 2e-16 ***
#   PCA_Index_Days6and.7  -18.527      2.734  -6.777 1.26e-10 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 22.59 on 207 degrees of freedom
# (41 observations deleted due to missingness)
# Multiple R-squared:  0.1816,	Adjusted R-squared:  0.1776 
# F-statistic: 45.92 on 1 and 207 DF,  p-value: 1.255e-10

sqrt(0.1816)
#[1] 0.4261455

pdf("F2_EPM_Time_ImmobileVsPavCA_ColSex.pdf", width=6, height=6)
plot(EPM_Time_Immobile~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM Time Immobile", pch=18, xlab="PavCA Index: Days 6 & 7", col=as.numeric(Behavior$Sex_AsFactor[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"]))
BestFitLine<-lm(EPM_Time_Immobile~PCA_Index_Days6and.7, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",])
abline(BestFitLine, col=1, lwd=3)
dev.off()

summary.lm(lm(EPM_Time_Immobile~PCA_Index_Days6and.7*Sex_AsFactor, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",]))

# 
# Call:
#   lm(formula = EPM_Time_Immobile ~ PCA_Index_Days6and.7 * Sex_AsFactor, 
#      data = Behavior[Behavior$RNASeq == "Elaine" & Behavior$Study_Generation == 
#                        "F2", ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -47.601 -12.667  -1.786  11.954  62.294 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         108.308      2.053  52.744  < 2e-16 ***
#   PCA_Index_Days6and.7                -13.221      3.828  -3.454 0.000671 ***
#   Sex_AsFactorF                       -19.770      3.382  -5.846 1.96e-08 ***
#   PCA_Index_Days6and.7:Sex_AsFactorF    4.415      5.677   0.778 0.437692    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 20.97 on 205 degrees of freedom
# (41 observations deleted due to missingness)
# Multiple R-squared:  0.3013,	Adjusted R-squared:  0.2911 
# F-statistic: 29.47 on 3 and 205 DF,  p-value: 6.929e-16






#General notes:
#The interaction term with sex is never significant, so the relationship between each of the variables seems to be similar in males and females (even when there is a sex effect in the variables themselves) - odd exception: % Time Open Arms vs. Boli.
#Controlling for sex effects never makes any of the relationships between variables disappear, although it does weaken some of the relationships with PavCA and EPM Boli


#################################



table(Behavior$Family[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"], Behavior$Sex[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2"])

#     F  M
# 1   7 10
# 2  10 11
# 3  10 12
# 4  13 13
# 5   7 12
# 6  10  9
# 7   9 12
# 8   9 11
# 9  13 14
# 10  8  8
# 11 21  8
# 12  8  5

Behavior$Family_AsFactor<-as.factor(Behavior$Family)

pdf("Total_LocoScore_Family_BySex_F2_Color.pdf", width=16, height=8)
boxplot(Total_LocoScore~Sex_AsFactor*Family, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="LocoScore", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("grey57", "grey80", "grey57", "grey80", "grey57", "grey80"))
stripchart(Total_LocoScore~Sex_AsFactor*Family, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(Total_LocoScore~Sex_AsFactor+Family_AsFactor, data=Behavior, contrasts=list(Family_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: Total_LocoScore
# Sum Sq  Df   F value    Pr(>F)    
# (Intercept)     247909065   1 1317.0149 < 2.2e-16 ***
#   Sex_AsFactor      2373795   1   12.6108 0.0004127 ***
#   Family_AsFactor   1574629  11    0.7605 0.6798795    
# Residuals       116517825 619                        
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 

pdf("EPM_Percent_Time_Open_Arm_Family_BySex_F2_Color.pdf", width=16, height=8)
boxplot(EPM_Percent_Time_Open_Arm~Sex_AsFactor*Family, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM: % Time Open Arm", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("grey57", "grey80", "grey57", "grey80", "grey57", "grey80"))
stripchart(EPM_Percent_Time_Open_Arm~Sex_AsFactor*Family, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(EPM_Percent_Time_Open_Arm~Sex_AsFactor+Family_AsFactor, data=Behavior, contrasts=list(Family_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: EPM_Percent_Time_Open_Arm
# Sum Sq  Df  F value  Pr(>F)    
# (Intercept)     149887   1 839.7813 < 2e-16 ***
#   Sex_AsFactor       145   1   0.8105 0.36833    
# Family_AsFactor   3273  11   1.6670 0.07703 .  
# Residuals       110481 619                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


pdf("EPM_Distance_Traveled_Family_BySex_F2_Color.pdf", width=16, height=8)
boxplot(EPM_Distance_Traveled~Sex_AsFactor*Family, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM: Distance Traveled", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("grey57", "grey80", "grey57", "grey80", "grey57", "grey80"))
stripchart(EPM_Distance_Traveled~Sex_AsFactor*Family, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(EPM_Distance_Traveled~Sex_AsFactor+Family_AsFactor, data=Behavior, contrasts=list(Family_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: EPM_Distance_Traveled
# Sum Sq  Df  F value    Pr(>F)    
# (Intercept)     3729407984   1 758.3546 < 2.2e-16 ***
#   Sex_AsFactor       9993793   1   2.0322    0.1545    
# Family_AsFactor  282046104  11   5.2139 6.992e-08 ***
#   Residuals       3044095030 619                       
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



pdf("EPM_Time_Immobile_Family_BySex_F2_Color.pdf", width=16, height=8)
boxplot(EPM_Time_Immobile~Sex_AsFactor*Family, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM: Time Immobile", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("grey57", "grey80", "grey57", "grey80", "grey57", "grey80"))
stripchart(EPM_Time_Immobile~Sex_AsFactor*Family, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(EPM_Time_Immobile~Sex_AsFactor+Family_AsFactor, data=Behavior, contrasts=list(Family_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: EPM_Time_Immobile
# Sum Sq  Df  F value    Pr(>F)    
# (Intercept)     3455054   1 2555.873 < 2.2e-16 ***
#   Sex_AsFactor      33364   1   24.681 8.764e-07 ***
#   Family_AsFactor   14484  11    0.974    0.4688    
# Residuals        836770 619                       
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


pdf("EPM_Boli_Family_BySex_F2_Color.pdf", width=16, height=8)
boxplot(EPM_Boli~Sex_AsFactor*Family, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="EPM: Boli", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("grey57", "grey80", "grey57", "grey80", "grey57", "grey80"))
stripchart(EPM_Boli~Sex_AsFactor*Family, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(EPM_Boli~Sex_AsFactor+Family_AsFactor, data=Behavior, contrasts=list(Family_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: EPM_Boli
# Sum Sq  Df  F value Pr(>F)    
# (Intercept)      315.50   1 152.7791 <2e-16 ***
#   Sex_AsFactor     155.39   1  75.2496 <2e-16 ***
#   Family_AsFactor   27.32  11   1.2028 0.2814    
# Residuals       1270.01 615                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

pdf("PavCA_Family_BySex_F2_Color.pdf", width=16, height=8)
boxplot(PCA_Index_Days6and.7~Sex_AsFactor*Family, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], ylab="PavCA Index: Days 6 & 7", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("grey57", "grey80", "grey57", "grey80", "grey57", "grey80"))
stripchart(PCA_Index_Days6and.7~Sex_AsFactor*Family, data=Behavior[Behavior$RNASeq=="Elaine"&Behavior$Study_Generation=="F2",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
dev.off()

Anova(lm(PCA_Index_Days6and.7~Sex_AsFactor+Family_AsFactor, data=Behavior, contrasts=list(Family_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: PCA_Index_Days6and.7
# Sum Sq  Df F value    Pr(>F)    
# (Intercept)      4.433   1 17.7810 3.783e-05 ***
#   Sex_AsFactor     9.613   1 38.5583 3.114e-09 ***
#   Family_AsFactor  6.160  11  2.2462   0.01362 *  
#   Residuals       48.865 196                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



table(Behavior$Family, Behavior$Mother)

# 1  3  5  7  9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 53 54 65 66 76 78 87 90 99 104 113 115 124 127 132 133 142 145 151 164 167 176 180 189 191 198 203 208 209 217 220 231 232 240 242 255 256 266 268 275 278 285 286 295 298
# 1  12 12  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 11 11  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  12  12   0   0   0   0   0   0   0   0   0   0   0   0   0   0
# 2   0  0 12 12  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 11  11   0   0   0   0   0   0   0   0   0  11  12   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
# 3   0  0  0  0 12 12  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 12 12  0  0  0   0  12  12   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
# 4   0  0  0  0  0  0 12 12  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  12  12   0   0   0   0  11  12   0   0   0   0   0   0
# 5   0  0  0  0  0  0  0  0 12 12  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0  12  10   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  12  12   0   0   0   0
# 6   0  0  0  0  0  0  0  0  0  0 12  6  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  12  11   0   0   0   0   0   0  12  12   0   0   0   0   0   0   0   0
# 7   0  0  0  0  0  0  0  0  0  0  0  0 11 12  0  0  0  0  0  0  0  0  0  0  0  0  8 12  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   9  11   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
# 8   0  0  0  0  0  0  0  0  0  0  0  0  0  0 12  8  0  0  0  0  0  0  0  0 11 12  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  12  12   0   0   0   0   0   0   0   0   0   0
# 9   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  6 11  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   6   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  12  12   0   0
# 10  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 12 12  0  0  0  0  0  0  0  0  0  0  0  0  0   0   0   0  12  12   0   0   0   0   0   0   0   0   0   0   0  12  12   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
# 11  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 12 10  0  0  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0  12  12   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
# 12  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 12  5  0  0  0  0  0  0  0  0  0   0   0   0   0   0  12  12   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  12  12
# 
# 304 306
# 1    0   0
# 2    0   0
# 3    0   0
# 4    0   0
# 5    0   0
# 6    0   0
# 7    0   0
# 8    0   0
# 9    0   0
# 10   0   0
# 11  12  12
# 12   0   0


table(Behavior$Family, Behavior$Father)

# 2  4  6  8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 49 52 57 63 70 73 82 83 93 98 105 109 120 122 130 131 137 139 147 149 158 161 170 171 183 186 195 196 204 207 211 213 224 227 236 238 248 249 259 261 271 272 282 288 293
# 1  12 12  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 12 12  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  11  11   0   0   0   0   0   0   0   0   0   0   0   0   0
# 2   0  0 12 12  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 11 12   0   0   0   0   0   0   0   0   0   0  11  11   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
# 3   0  0  0  0 12 12  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 12 12  0  0  0  0  12  12   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
# 4   0  0  0  0  0  0 12 12  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  12  11   0   0   0   0  12  12   0   0   0   0   0
# 5   0  0  0  0  0  0  0  0 12 12  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0  12  12   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  10  12   0   0   0
# 6   0  0  0  0  0  0  0  0  0  0 12  6  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  12  12   0   0   0   0   0   0  12  11   0   0   0   0   0   0   0
# 7   0  0  0  0  0  0  0  0  0  0  0  0 11 12  0  0  0  0  0  0  0  0  0  0  0  0 11  9  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   8  12   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
# 8   0  0  0  0  0  0  0  0  0  0  0  0  0  0 12  8  0  0  0  0  0  0  0  0 12 12  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  11  12   0   0   0   0   0   0   0   0   0
# 9   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  6 11  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0  12  12   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   6   0   0
# 10  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 12 12  0  0  0  0  0  0  0  0  0  0  0  0  0  0   0   0  12  12   0   0   0   0   0   0   0   0   0   0   0   0  12  12   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
# 11  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 12 10  0  0  0  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0  12  12   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
# 12  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 12  5  0  0  0  0  0  0  0  0  0  0   0   0   0   0  12  12   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  12  12
# 
# 300 301
# 1    0   0
# 2    0   0
# 3    0   0
# 4    0   0
# 5    0   0
# 6    0   0
# 7    0   0
# 8    0   0
# 9    0   0
# 10   0   0
# 11  12  12
# 12   0   0


cbind(Behavior$Family[Behavior$Family=="1"], Behavior$Mother[Behavior$Family=="1"], Behavior$Father[Behavior$Family=="1"])

#The parents of the F2s are all siblings from the same F1 family (HRXLR cross) but might be from different litters from the same parents???
#So asking about "family" means asking about influence of grandparents, whereas asking about "Mother" (or Father) means asking about litter.
#And asking about litter in this case, means partially asking about selective breeding, because the most extreme males and females were chosen to be bred together from each family.


#Conclusion either mother or father can designate a litter for the F2 generation.

Behavior$Mother_AsFactor<-as.factor(Behavior$Mother)

Anova(lm(Total_LocoScore~Sex_AsFactor+Mother_AsFactor, data=Behavior, contrasts=list(Mother_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: Total_LocoScore
# Sum Sq  Df   F value    Pr(>F)    
# (Intercept)     199322090   1 3679.0089 < 2.2e-16 ***
#   Sex_AsFactor      1818819   1   33.5710 1.178e-08 ***
#   Mother_AsFactor  11780197  51    4.2634 < 2.2e-16 ***
#   Residuals        28768625 531                        
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Anova(lm(EPM_Percent_Time_Open_Arm~Sex_AsFactor+Mother_AsFactor, data=Behavior, contrasts=list(Mother_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: EPM_Percent_Time_Open_Arm
# Sum Sq  Df  F value    Pr(>F)    
# (Intercept)     125516   1 826.4729 < 2.2e-16 ***
#   Sex_AsFactor       752   1   4.9508  0.026497 *  
#   Mother_AsFactor  12094  51   1.5615  0.009659 ** 
#   Residuals        80642 531                       
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Anova(lm(EPM_Distance_Traveled~Sex_AsFactor+Mother_AsFactor, data=Behavior, contrasts=list(Mother_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: EPM_Distance_Traveled
# Sum Sq  Df  F value    Pr(>F)    
# (Intercept)     3299646625   1 696.4698 < 2.2e-16 ***
#   Sex_AsFactor      10944897   1   2.3102    0.1291    
# Mother_AsFactor  736886805  51   3.0498 1.366e-10 ***
#   Residuals       2515704564 531                       
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Anova(lm(EPM_Time_Immobile~Sex_AsFactor+Mother_AsFactor, data=Behavior, contrasts=list(Mother_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: EPM_Time_Immobile
# Sum Sq  Df  F value    Pr(>F)    
# (Intercept)     3156129   1 7298.799 < 2.2e-16 ***
#   Sex_AsFactor      18191   1   42.069  2.02e-10 ***
#   Mother_AsFactor  470586  51   21.339 < 2.2e-16 ***
#   Residuals        229614 531                       
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Anova(lm(EPM_Boli~Sex_AsFactor+Mother_AsFactor, data=Behavior, contrasts=list(Mother_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)
# Anova Table (Type III tests)
# 
# Response: EPM_Boli
# Sum Sq  Df  F value Pr(>F)    
# (Intercept)     234.53   1 127.4851 <2e-16 ***
#   Sex_AsFactor    153.70   1  83.5492 <2e-16 ***
#   Mother_AsFactor  93.13  51   0.9927 0.4914    
# Residuals       969.49 527                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


Anova(lm(PCA_Index_Days6and.7~Sex_AsFactor+Mother_AsFactor, data=Behavior, contrasts=list(Mother_AsFactor=contr.sum, Sex_AsFactor=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: PCA_Index_Days6and.7
# Sum Sq  Df F value   Pr(>F)    
# (Intercept)      4.056   1 16.7369 6.47e-05 ***
#   Sex_AsFactor     7.626   1 31.4641 7.55e-08 ***
#   Mother_AsFactor 11.400  27  1.7421  0.01792 *  
#   Residuals       43.625 180                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



