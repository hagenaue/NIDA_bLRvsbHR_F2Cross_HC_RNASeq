#F2 HC RNA-Seq Dataset
#19_Looking at the dissection variable in a little more depth
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.

#########################################2/11/2022

levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Releveled)
#[1] "12/3/2020"  "12/1/2020"  "12/10/2020" "12/11/2020" "12/14/2020" "12/2/2020"  "12/4/2020"  "12/7/2020"  "12/8/2020"  "12/9/2020" 

#F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Releveled<-relevel(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor, ref="12/3/2020")


F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted<-factor(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor, levels=c("12/1/2020","12/2/2020","12/3/2020","12/4/2020","12/7/2020", "12/8/2020","12/9/2020","12/10/2020", "12/11/2020", "12/14/2020"))

levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)
#[1] "12/1/2020"  "12/2/2020"  "12/3/2020"  "12/4/2020"  "12/7/2020"  "12/8/2020"  "12/9/2020"  "12/10/2020" "12/11/2020" "12/14/2020"

# 
# pdf("TotalLocoScore_RNAseq_SexGen_MaleRef_Color.pdf", width=9.5, height=6)
# boxplot(Total_LocoScore~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], ylab="LocoScore", xlab="", cex.lab=1.5, cex.axis=1.2, col=c("green4","green2","firebrick4", "firebrick2", "grey57", "grey80"))
# stripchart(Total_LocoScore~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')
# 
# dev.off()

pdf("RNAconc_RNAseq_HPC_Dissection_Date_byDissector.pdf", width=9.5, height=6)
boxplot(RNAconc~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, ylab="RNA Conc (ng/ul)", xlab="", las=2,cex.lab=1.5, cex.axis=1.2, col=c("grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57"))
stripchart(RNAconc~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()


pdf("RibosomePerc_RNAseq_HPC_Dissection_Date_byDissector.pdf", width=9.5, height=6)
boxplot(RibosomePerc~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, ylab="% rRNA", xlab="", las=2, cex.lab=1.5, cex.axis=1.2, col=c("grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57"))
stripchart(RibosomePerc~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()


pdf("Percent_Intergenic_RNAseq_HPC_Dissection_Date_byDissector.pdf", width=9.5, height=6)
boxplot(Percent.Intergenic~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, ylab="% Intergenic", xlab="", las=2, cex.lab=1.5, cex.axis=1.2, col=c("grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57"))
stripchart(Percent.Intergenic~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()

####shorten HPC_Dissection_Date_Sorted values

levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)[levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)=="12/1/2020"]<-"1st"
levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)[levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)=="12/2/2020"]<-"2nd"
levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)[levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)=="12/3/2020"]<-"3rd"
levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)[levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)=="12/4/2020"]<-"4th"
levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)[levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)=="12/7/2020"]<-"5th"
levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)[levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)=="12/8/2020"]<-"6th"
levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)[levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)=="12/9/2020"]<-"7th"
levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)[levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)=="12/10/2020"]<-"8th"
levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)[levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)=="12/11/2020"]<-"9th"
levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)[levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)=="12/14/2020"]<-"10th"

levels(F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Sorted)
#[1] "1st"  "2nd"  "3rd"  "4th"  "5th"  "6th"  "7th"  "8th"  "9th"  "10th"

str(F2_PCA_wMetaData)


###Re-do3 boxplots from above to reflect shortened HPC_Dissection_Date_AsFactor_Sorted

pdf("RNAconc_RNAseq_HPC_Dissection_Date_byDissector.pdf", width=9.5, height=6)
boxplot(RNAconc~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, ylab="RNA Conc (ng/ul)", xlab="", las=2,cex.lab=1.5, cex.axis=1.2, col=c("grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57"))
stripchart(RNAconc~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()

pdf("RibosomePerc_RNAseq_HPC_Dissection_Date_byDissector.pdf", width=9.5, height=6)
boxplot(RibosomePerc~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, ylab="% rRNA", xlab="", las=2, cex.lab=1.5, cex.axis=1.2, ylim=c(0, 0.03), col=c("grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57"))
stripchart(RibosomePerc~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()


pdf("Percent_Intergenic_RNAseq_HPC_Dissection_Date_byDissector.pdf", width=9.5, height=6)
boxplot(Percent.Intergenic~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, ylab="% Intergenic", xlab="", las=2, cex.lab=1.5, cex.axis=1.2, col=c("grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57"))
stripchart(Percent.Intergenic~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()

##ylim begins at 0 here
pdf("Percent_Intergenic_RNAseq_HPC_Dissection_Date_byDissector_ylim0.pdf", width=9.5, height=6)
boxplot(Percent.Intergenic~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, ylab="% Intergenic", xlab="", las=2, cex.lab=1.5, cex.axis=1.2, ylim=c(0, 0.2), col=c("grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57"))
stripchart(Percent.Intergenic~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()

#######Boxplots without Y axis label

pdf("RNAconc_RNAseq_HPC_Dissection_Date_byDissector_No_YaxisLabel.pdf", width=9.5, height=6)
boxplot(RNAconc~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, ylab="", xlab="", las=2,cex.lab=1.5, cex.axis=1.2, col=c("grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57"))
stripchart(RNAconc~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()


pdf("RibosomePerc_RNAseq_HPC_Dissection_Date_byDissector_No_YaxisLabel.pdf", width=9.5, height=6)
boxplot(RibosomePerc~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, ylab="", xlab="", las=2, cex.lab=1.5, cex.axis=1.2, ylim=c(0, 0.03), col=c("grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57"))
stripchart(RibosomePerc~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()


pdf("Percent_Intergenic_RNAseq_HPC_Dissection_Date_byDissector_No_YaxisLabel.pdf", width=9.5, height=6)
boxplot(Percent.Intergenic~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, ylab="", xlab="", las=2, cex.lab=1.5, cex.axis=1.2, col=c("grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57"))
stripchart(Percent.Intergenic~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()

##ylim begins at 0 here
pdf("Percent_Intergenic_RNAseq_HPC_Dissection_Date_byDissector_ylim0_No_YaxisLabel.pdf", width=9.5, height=6)
boxplot(Percent.Intergenic~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, ylab="", xlab="", las=2, cex.lab=1.5, cex.axis=1.2, ylim=c(0, 0.2), col=c("grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57","grey80","grey57"))
stripchart(Percent.Intergenic~F2_PCA_wMetaData$Dissector_AsFactor*HPC_Dissection_Date_AsFactor_Sorted, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=0.5, col='black')

dev.off()







