#F2 HC RNA-Seq Dataset
#20_Making violin plots to illustrate the distributions for the covariates/RNA metrics in the F0 and F2 datasets
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.


write.csv(F2_PCA_wMetaData,"F2_PCA_wMetaData_ViolinPlots")

####Added F0 RNAmetrics to file and re-named it: F0_F2_PCA_wMetaData_ViolinPlots

F0_F2_RNAmetric_ViolinPlots<-read.delim("F0_F2_PCA_wMetaData_ViolinPlots2.txt", sep='\t', header=TRUE, stringsAsFactors = FALSE)

str(F0_F2_RNAmetric_ViolinPlots)
# 'data.frame':	268 obs. of  15 variables:
#   $ row_name          : chr  "sl469967" "sl469968" "sl469969" "sl469970" ...
# $ Gen               : chr  "F2" "F2" "F2" "F2" ...
# $ RNAconc           : num  127.2 98.5 143.9 191.4 125.3 ...
# $ PF.Reads          : int  30800397 27366626 29238544 38111688 30064848 26047278 26871577 27028205 25643674 28446948 ...
# $ RIN               : num  8.8 8.8 9 9 9 9 8.9 9.1 9.1 9.1 ...
# $ DV200             : num  96.4 96.6 96.8 96.6 96.9 96.8 96.5 96.8 96.6 96.8 ...
# $ LibrarySize       : int  21826046 19217774 20269455 26635890 21240639 18267068 18970374 18962467 17916830 19853360 ...
# $ Percent.Ribosomal : num  0.0096 0.0101 0.0086 0.0082 0.0088 0.0072 0.0084 0.0112 0.0075 0.0072 ...
# $ Percent.Coding    : num  0.476 0.466 0.452 0.466 0.459 ...
# $ Percent.UTR       : num  0.327 0.335 0.33 0.326 0.336 ...
# $ Percent.Intronic  : num  0.0316 0.032 0.0383 0.041 0.0345 ...
# $ Percent.Intergenic: num  0.156 0.158 0.171 0.159 0.161 ...
# $ Percent.mRNA      : num  0.802 0.8 0.782 0.792 0.796 ...
# $ PCT_USABLE_BASES  : num  0.784 0.782 0.764 0.775 0.779 ...
# $ MEDIAN_CV_COVERAGE: num  0.445 0.432 0.424 0.422 0.423 ...

F0_F2_RNAmetric_ViolinPlots$Gen<-as.factor(F0_F2_RNAmetric_ViolinPlots$Gen)

str(F0_F2_RNAmetric_ViolinPlots)
# 
# 'data.frame':	268 obs. of  15 variables:
#   $ row_name          : chr  "sl469967" "sl469968" "sl469969" "sl469970" ...
# $ Gen               : Factor w/ 2 levels "F0","F2": 2 2 2 2 2 2 2 2 2 2 ...
# $ RNAconc           : num  127.2 98.5 143.9 191.4 125.3 ...
# $ PF.Reads          : int  30800397 27366626 29238544 38111688 30064848 26047278 26871577 27028205 25643674 28446948 ...
# $ RIN               : num  8.8 8.8 9 9 9 9 8.9 9.1 9.1 9.1 ...
# $ DV200             : num  96.4 96.6 96.8 96.6 96.9 96.8 96.5 96.8 96.6 96.8 ...
# $ LibrarySize       : int  21826046 19217774 20269455 26635890 21240639 18267068 18970374 18962467 17916830 19853360 ...
# $ Percent.Ribosomal : num  0.0096 0.0101 0.0086 0.0082 0.0088 0.0072 0.0084 0.0112 0.0075 0.0072 ...
# $ Percent.Coding    : num  0.476 0.466 0.452 0.466 0.459 ...
# $ Percent.UTR       : num  0.327 0.335 0.33 0.326 0.336 ...
# $ Percent.Intronic  : num  0.0316 0.032 0.0383 0.041 0.0345 ...
# $ Percent.Intergenic: num  0.156 0.158 0.171 0.159 0.161 ...
# $ Percent.mRNA      : num  0.802 0.8 0.782 0.792 0.796 ...
# $ PCT_USABLE_BASES  : num  0.784 0.782 0.764 0.775 0.779 ...
# $ MEDIAN_CV_COVERAGE: num  0.445 0.432 0.424 0.422 0.423 ...

library(ggplot2)
# Basic violin plot
RNAconc <- ggplot(F0_F2_RNAmetric_ViolinPlots$RNAconc), aes(y=RNAconc)) + 
  geom_violin()

#Error: `data` must be a data frame, or other object coercible by `fortify()`, not a numeric vector.

#all of these did not work####
# RNAconc <- ggplot(data=data.frame(F0_F2_RNAmetric_ViolinPlots$RNAconc), aes(y=RNAconc)) + 
#   geom_violin()
# 
# RNAconc <- ggplot(F0_F2_RNAmetric_ViolinPlots, aes(x=Gen, y=RNAconc)) + 
#   geom_violin()
# 
# pdf("RNAconc_F0_F2_Violin.pdf", height=6, width=6)
# p <- ggplot(F0_F2_RNAmetric_ViolinPlots, aes(x=Gen, y=RNAconc)) + 
#   geom_violin()
# dev.off()
# 
# pdf("RNAconc_F0_F2_Violin.pdf", height=6, width=6)
# p <- ggplot(data=F0_F2_RNAmetric_ViolinPlots, aes(x=Gen, y=RNAconc)) + 
#   geom_violin()
# dev.off()
# 
# pdf("RNAconc_F0_F2_Violin.pdf", height=6, width=6)
# p <- ggplot(data=F0_F2_RNAmetric_ViolinPlots$RNAconc, aes(x=Gen, y=RNAconc)) + 
#   geom_violin()
# dev.off()
# 
# pdf("RNAconc_F0_F2_Violin.pdf", height=6, width=6)
# p <- F0_F2_RNAmetric_ViolinPlots$RNAconc (aes(x=Gen, y=RNAconc) + 
#   geom_violin()
# dev.off()
# 
# ggplot(data=F0_F2_RNAmetric_ViolinPlots$RNAconc, aes(x=Gen, y=RNAconc)) +
#   geom_violin()
# 
# 
# ggplot(data=F0_F2_RNAmetric_ViolinPlots, aes(x=Gen, y=RNAconc)) +
#   geom_violin()
# 
# pdf("RNAconc_F0_F2_Violin.pdf", height=6, width=6)
# ggplot(F0_F2_RNAmetric_ViolinPlots, aes(x=Gen, y=RNAconc)) +
#   geom_violin()
# 
# pdf("RNAconc_F0_F2_Violin.pdf", height=6, width=6)
# ggplot(data=F0_F2_RNAmetric_ViolinPlots, aes(x=Gen, y=RNAconc)) +
#   geom_violin()


####This works######
pdf("RNAconc_F0_F2_Violin.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), RNAconc))
p + geom_violin()
dev.off()

###DOn't like the way this one looks--horizontal lines are too wide###################
pdf("RNAconc_F0_F2_Violin_quantiles.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), RNAconc))
p + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
dev.off()

################This one works: median and quartiles ##############################
pdf("RNAconc_F0_F2_Violin_median_quartile.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), RNAconc))
p + geom_violin() + stat_summary(fun=median, geom = "point", size=2, color="black") + geom_boxplot(width=0.1)
dev.off()

################This one works: mean ##############################
pdf("RNAconc_F0_F2_Violin_mean.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), RNAconc))
p + geom_violin() + stat_summary(fun=mean, geom="point", shape=23, size=2)
dev.off()                                 

################This one is better: median and thinner quartile##############################
pdf("RNAconc_F0_F2_Violin_median_quartile_thinner.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), RNAconc))
p + geom_violin() + stat_summary(fun=median, geom = "point", size=2, color="black") + geom_boxplot(width=0.05)
dev.off()


################This one works: with jitter but many points outside violin##############################
pdf("RNAconc_F0_F2_Violin_wJitter.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), RNAconc))
p + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2))
dev.off()   


######################this one works: reverse color_ BG white, violins grey   median and quartile########################
pdf("RNAconc_F0_F2_Violin_median_quartile_thinner_BGwhite_violinGrey_Trim.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), RNAconc))
p + geom_violin(fill='#A4A4A4') + stat_summary(fun=median, geom = "point", size=2, color="black") + geom_boxplot(width=0.05) + theme_minimal()
dev.off()

######################this one works: reverse color_ BG white, violins grey  median and quartile########################
pdf("RNAconc_F0_F2_Violin_median_quartile_thinner_BGwhite_violinGrey_NoTrim.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), RNAconc))
p + geom_violin(trim=FALSE, fill='#A4A4A4') + stat_summary(fun=median, geom = "point", size=2, color="black") + geom_boxplot(width=0.05) + theme_minimal()
dev.off()


#####3/29/2022 More violin plots use no trim version################################

pdf("LibrarySize_F0_F2_Violin_median_quartile_thinner_BGwhite_violinGrey_NoTrim.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), LibrarySize))
p + geom_violin(trim=FALSE, fill='#A4A4A4') + stat_summary(fun=median, geom = "point", size=2, color="black") + geom_boxplot(width=0.05) + theme_minimal()
dev.off()

pdf("RIN_F0_F2_Violin_median_quartile_thinner_BGwhite_violinGrey_NoTrim.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), RIN))
p + geom_violin(trim=FALSE, fill='#A4A4A4') + stat_summary(fun=median, geom = "point", size=2, color="black") + geom_boxplot(width=0.05) + theme_minimal()
dev.off()

pdf("DV200_F0_F2_Violin_median_quartile_thinner_BGwhite_violinGrey_NoTrim.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), DV200))
p + geom_violin(trim=FALSE, fill='#A4A4A4') + stat_summary(fun=median, geom = "point", size=2, color="black") + geom_boxplot(width=0.05) + theme_minimal()
dev.off()

pdf("Percent.mRNA_F0_F2_Violin_median_quartile_thinner_BGwhite_violinGrey_NoTrim.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), Percent.mRNA))
p + geom_violin(trim=FALSE, fill='#A4A4A4') + stat_summary(fun=median, geom = "point", size=2, color="black") + geom_boxplot(width=0.05) + theme_minimal()
dev.off()

pdf("Percent.Ribosomal_F0_F2_Violin_median_quartile_thinner_BGwhite_violinGrey_NoTrim.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), Percent.Ribosomal))
p + geom_violin(trim=FALSE, fill='#A4A4A4') + stat_summary(fun=median, geom = "point", size=2, color="black") + geom_boxplot(width=0.05) + theme_minimal()
dev.off()

pdf("Percent.Coding_F0_F2_Violin_median_quartile_thinner_BGwhite_violinGrey_NoTrim.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), Percent.Coding))
p + geom_violin(trim=FALSE, fill='#A4A4A4') + stat_summary(fun=median, geom = "point", size=2, color="black") + geom_boxplot(width=0.05) + theme_minimal()
dev.off()

pdf("Percent.UTR_F0_F2_Violin_median_quartile_thinner_BGwhite_violinGrey_NoTrim.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), Percent.UTR))
p + geom_violin(trim=FALSE, fill='#A4A4A4') + stat_summary(fun=median, geom = "point", size=2, color="black") + geom_boxplot(width=0.05) + theme_minimal()
dev.off()

pdf("Percent.Intronic_F0_F2_Violin_median_quartile_thinner_BGwhite_violinGrey_NoTrim.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), Percent.Intronic))
p + geom_violin(trim=FALSE, fill='#A4A4A4') + stat_summary(fun=median, geom = "point", size=2, color="black") + geom_boxplot(width=0.05) + theme_minimal()
dev.off()

pdf("RNAconc_F0_F2_Violin_median_quartile_thinner_BGwhite_violinGrey_NoTrim.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), RNAconc))
p + geom_violin(trim=FALSE, fill='#A4A4A4') + stat_summary(fun=median, geom = "point", size=2, color="black") + geom_boxplot(width=0.05) + theme_minimal()
dev.off()

pdf("Percent.Intergenic_F0_F2_Violin_median_quartile_thinner_BGwhite_violinGrey_NoTrim.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), Percent.Intergenic))
p + geom_violin(trim=FALSE, fill='#A4A4A4') + stat_summary(fun=median, geom = "point", size=2, color="black") + geom_boxplot(width=0.05) + theme_minimal()
dev.off()

pdf("MEDIAN_CV_COVERAGE_F0_F2_Violin_median_quartile_thinner_BGwhite_violinGrey_NoTrim.pdf", height=6, width=6)
p<-ggplot(F0_F2_RNAmetric_ViolinPlots, aes(factor(Gen), MEDIAN_CV_COVERAGE))
p + geom_violin(trim=FALSE, fill='#A4A4A4') + stat_summary(fun=median, geom = "point", size=2, color="black") + geom_boxplot(width=0.05) + theme_minimal()
dev.off()