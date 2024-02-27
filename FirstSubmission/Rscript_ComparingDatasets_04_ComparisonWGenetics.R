#Comparing RNA-Seq Results across datasets: F0 vs. F2 vs. bHR/bLR Meta-analysis
#04_Comparing our differential expression results to the Goldman and Chitre genetic (QTL) results
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-27, updated later for a few figures for the paper.


#What also was a qtl in the Goldman study?

OldGoldmanQTLresults<-read.csv("EntrezID_GeneSymbol_GeneLocations_Rnor6_Rnor5_wQTLresults.csv", header=TRUE, stringsAsFactors = FALSE)

colnames(OldGoldmanQTLresults)

colnames(OldGoldmanQTLresults)[3]<-"SYMBOL"

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01")

colnames(F0_Meta_F2_4Models)

sum(is.na(OldGoldmanQTLresults$SYMBOL))
#[1] 0

head(OldGoldmanQTLresults$SYMBOL)

F0_Meta_F2_4Models_OldGoldmanQTLresults<-join(F0_Meta_F2_4Models, OldGoldmanQTLresults, by="SYMBOL", type="left")

#Note - joining by symbol is really not ideal. 
#We should get an up-to-date version of the QTLs from Abe's lab 
#and up-to-date ENSEMBL based annotation w/ gene locations for whichever genome build Abe is using.
#In the meantime: 
#I used Abe's regional association plots, and trying to encompass the full region with LOD>4, then added an 1 MB bin to each side (Rnor6)
#For any gene with both HR/LR associations and F2_Behavioral associations, I filled in any missing coordinates using NCBI's Genome

write.csv(F0_Meta_F2_4Models_OldGoldmanQTLresults, "F0_Meta_F2_LocomotorEPMSTGT_OldGoldmanQTLresults.csv")


#But good enough for preliminary results.


#Coming back to this, let's run a better analysis:

#First, reading in annotation for the genes (Rnor6 coordinates) derived from the same package/source as our earlier output for the F0s/F2s (from Elaine)

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01")

MoreAnnotation<-read.csv("EnsemblVsGeneSymbol_wCHRLOC.csv", header=TRUE, stringsAsFactors = FALSE)

str(MoreAnnotation)
# $ X        : int  1 2 3 4 5 6 7 8 9 10 ...
# $ ENSEMBL  : chr  "ENSRNOG00000017701" "ENSRNOG00000028896" "ENSRNOG00000028896" "ENSRNOG00000032908" ...
# $ SYMBOL   : chr  "Asip" "A2m" "A2m" "Acaa1a" ...
# $ CHRLOCCHR: chr  "3" "4" "4" "8" ...
# $ CHRLOC   : int  150492009 154423167 154309425 128027879 -260124417 -88392247 -49836064 -49836801 -99120378 -49836685 ...
# $ CHRLOCEND: int  150579870 154473038 154359138 128036471 -260148589 -88442845 -49851714 -49851648 -99121144 -49851660 ...

#Note that there are often more than one ChrLOC for a particular gene:
table(table(MoreAnnotation$ENSEMBL))
#     1     2     3     4     5     6     7     8     9    11 
# 19649  1052   186    74    13     8     2     9     3     2 
#But that is certainly the minority of genes. 

str(F0_Meta_F2_4Models)
#'data.frame':	13786 obs. of  213 variables:
#'
F0_Meta_F2_4Models_withRnor6Coord<-join(F0_Meta_F2_4Models, MoreAnnotation, by="ENSEMBL", type="left")

str(F0_Meta_F2_4Models_withRnor6Coord)
#'data.frame':	14449 obs. of  218 variables:
14449-13786
#[1] 663
#So 663 rows represent Ensembl genes that are represented in other rows.

#write.csv(F0_Meta_F2_4Models_withRnor6Coord, "F0_Meta_F2_4Models_withRnor6Coord.csv")
write.csv(F0_Meta_F2_4Models_withRnor6Coord, "F0_Meta_F2_4Models_withRnor6Coord_Updated.csv")


#Hmm.... I've discovered a problem with this annotation that was also present in the old annotation:
#We have several top genes that have ENSEMBL annotation (but no gene symbol) that *are not included* in this annotation even though they are present in the Rnor6 in NCBI's Genome Data Viewer: https://www.ncbi.nlm.nih.gov/genome/gdv/browser/genome/?id=GCF_000001895.5
#... and some of these are top hits.
#Why?  I can go back an add in the annotation for the top hits by hand again, but I was assuming that they were missing previously because my annotation was old... now I'm not sure why they are missing. Elaine's version of the annotation package was downloaded recently. (???)  Maybe the package is out of date????
#.... or maybe something about how we outputted it made it so that nothing was outputtted that didn't have a gene symbol? Let's check:

sum(is.na(MoreAnnotation$SYMBOL))
#[1] 0
#Ah. :(
#Yep - just took a quick peek by hand in Excel, and it definitely looks like nothing outputted that lacked a gene symbol - I'm guessing that is because of how we outputted the annotation from the package (versus a deficiency in the package itself) - we'll have to check.

#Time to try another package:

#https://bioconductor.org/packages/devel/bioc/vignettes/ensembldb/inst/doc/ensembldb.html

library(ensembldb)

#library(EnsDb.Rnor.v103)

# 
# DB <- ensDbFromGRanges(Y, path = tempdir(), version = 103,
#                        organism = "Rattus_norvegicus")


library(AnnotationHub)

## Load the annotation resource.
ah <- AnnotationHub()

## Query for all available EnsDb databases
query(ah, "EnsDb")

ahDb <- query(ah, pattern = c("Rattus Norvegicus", "EnsDb", 103))

ahDb
# AnnotationHub with 0 records
# # snapshotDate(): 2017-04-25 

ahDb <- query(ah, pattern = c("Rattus Norvegicus", "EnsDb", 99))

ahDb
# AnnotationHub with 0 records
# # snapshotDate(): 2017-04-25 

ahDb <- query(ah, pattern = c("Rattus Norvegicus", "EnsDb"))
ahDb
#AnnotationHub with 2 records
# snapshotDate(): 2017-04-25 
# $dataprovider: Ensembl
# $species: Rattus norvegicus
# $rdataclass: EnsDb
# additional mcols(): taxonomyid, genome, description, coordinate_1_based, maintainer,
#   rdatadateadded, preparerclass, tags, rdatapath, sourceurl, sourcetype 
# retrieve records with, e.g., 'object[["AH53239"]]' 

# title                                 
# AH53239 | Ensembl 87 EnsDb for Rattus Norvegicus
# AH53743 | Ensembl 88 EnsDb for Rattus Norvegicus

#hmmm... that's a little old (definitely older than the version we used for the RNA-Seq alignment)

ahEdb <- ahDb[[2]]

str(ahEdb)

## retrieve all genes
gns <- genes(ahEdb)

str(gns)

#Seeing if this improved our situation: Does it have annotation for the ENSEMBL genes that are top hits in our study but lack annotation?
gns[gns$gene_id=="ENSRNOG00000052237"]

# GRanges object with 1 range and 7 metadata columns:
#   seqnames               ranges strand |            gene_id      gene_name gene_biotype
# <Rle>            <IRanges>  <Rle> |        <character>    <character>  <character>
#   ENSRNOG00000052237        1 [94836296, 94840706]      + | ENSRNOG00000052237 AABR07071904.1      lincRNA
# seq_coord_system description         symbol entrezid
# <character> <character>    <character>   <list>
#   ENSRNOG00000052237       chromosome        NULL AABR07071904.1       NA
# -------
#   seqinfo: 162 sequences from Rnor_6.0 genome

#Perfect - this actually matches the Rnor6 annotation on the NCBI Genome Viewer

gids <- keys(ahEdb, keytype = "GENEID")
length(gids)
#[1] 32883

gids[c(1:5)]
#[1] "ENSRNOG00000000001" "ENSRNOG00000000007" "ENSRNOG00000000008" "ENSRNOG00000000009" "ENSRNOG00000000010"

columns(ahEdb)

listColumns(ahEdb)

# [1] "seq_name"              "seq_length"            "is_circular"           "gene_id"              
# [5] "entrezid"              "exon_id"               "exon_seq_start"        "exon_seq_end"         
# [9] "gene_name"             "gene_biotype"          "gene_seq_start"        "gene_seq_end"         
# [13] "seq_strand"            "seq_coord_system"      "description"           "symbol"               
# [17] "name"                  "value"                 "tx_id"                 "protein_id"           
# [21] "protein_sequence"      "protein_domain_id"     "protein_domain_source" "interpro_accession"   
# [25] "prot_dom_start"        "prot_dom_end"          "tx_biotype"            "tx_seq_start"         
# [29] "tx_seq_end"            "tx_cds_seq_start"      "tx_cds_seq_end"        "tx_support_level"     
# [33] "tx_name"               "exon_idx"              "uniprot_id"            "uniprot_db"           
# [37] "uniprot_mapping_type" 

#test run:
select(ahEdb, keys = gids[c(1:5)], keytype = "GENEID",
       columns = c("GENEID", "GENENAME", "SYMBOL", "ENTREZID", "SEQCOORDSYSTEM", "SEQNAME", "SEQSTRAND",  "GENESEQSTART", "GENESEQEND", "GENEBIOTYPE",  "EXONID", "EXONSEQSTART", "EXONSEQEND", "TXBIOTYPE",        "TXSEQSTART", "TXSEQEND", "ISCIRCULAR"))

ENSEMBL_Rnor6v88_Annotation<-select(ahEdb, keys = gids, keytype = "GENEID",
                                    columns = c("GENEID", "GENENAME", "SYMBOL", "ENTREZID", "SEQCOORDSYSTEM", "SEQNAME", "SEQSTRAND",  "GENESEQSTART", "GENESEQEND", "GENEBIOTYPE",  "EXONID", "EXONSEQSTART", "EXONSEQEND", "TXBIOTYPE",        "TXSEQSTART", "TXSEQEND", "ISCIRCULAR"))

head(ENSEMBL_Rnor6v88_Annotation)

write.csv(ENSEMBL_Rnor6v88_Annotation, "ENSEMBL_Rnor6v88_Annotation.csv")

colnames(ENSEMBL_Rnor6v88_Annotation)[1]<-"ENSEMBL"

F0_Meta_F2_4Models_withRnor6Coord<-join(F0_Meta_F2_4Models, ENSEMBL_Rnor6v88_Annotation, by="ENSEMBL", type="left")

#write.csv(F0_Meta_F2_4Models_withRnor6Coord, "F0_Meta_F2_4Models_withRnor6Coord.csv")

write.csv(F0_Meta_F2_4Models_withRnor6Coord, "F0_Meta_F2_4Models_withRnor6v88Coord_Updated.csv")


str(F0_Meta_F2_4Models_withRnor6Coord)
#'data.frame':	224543 obs. of  241 variables:
224543/13786
#[1] 16.28776
#Ouch - so most genes are represented by many many rows.
#I'm guessing that is because of the transcript level annotation.
#Let's try to trim this down.

colnames(ENSEMBL_Rnor6v88_Annotation)
# [1] "ENSEMBL"        "GENENAME"       "SYMBOL"         "ENTREZID"       "SEQCOORDSYSTEM" "SEQNAME"       
# [7] "SEQSTRAND"      "GENESEQSTART"   "GENESEQEND"     "GENEBIOTYPE"    "EXONID"         "EXONSEQSTART"  
# [13] "EXONSEQEND"     "TXBIOTYPE"      "TXSEQSTART"     "TXSEQEND"       "ISCIRCULAR"  

ENSEMBL_Rnor6v88_Annotation_Genes<-ENSEMBL_Rnor6v88_Annotation[,c(1:10)]

head(ENSEMBL_Rnor6v88_Annotation_Genes)

str(unique(ENSEMBL_Rnor6v88_Annotation_Genes))
#'data.frame':	34169 obs. of  10 variables:
#'That's better...

ENSEMBL_Rnor6v88_Annotation_Genes<-unique(ENSEMBL_Rnor6v88_Annotation_Genes)

F0_Meta_F2_4Models_withRnor6Coord_JustGenes<-join(F0_Meta_F2_4Models, ENSEMBL_Rnor6v88_Annotation_Genes, by="ENSEMBL", type="left")

str(F0_Meta_F2_4Models_withRnor6Coord_JustGenes)
#'data.frame':	14148 obs. of  234 variables:
#'Soooooo much better.
14148-13786
#[1] 362
#Just 362 Ensembl ids mapped to more than one set of annotation - workable.

#Overwriting earlier object so I can reuse my downstream code:

F0_Meta_F2_4Models_withRnor6Coord<-F0_Meta_F2_4Models_withRnor6Coord_JustGenes

#write.csv(F0_Meta_F2_4Models_withRnor6Coord, "F0_Meta_F2_4Models_withRnor6Coord_JustGenes.csv")
write.csv(F0_Meta_F2_4Models_withRnor6Coord, "F0_Meta_F2_4Models_withRnor6v88Coord_JustGenes.csv")


########################

#Let's try making a manhattan plot!

#https://www.r-graph-gallery.com/101_Manhattan_plot.html

library(qqman)

# Citation appreciated but not required:
#   Turner, (2018). qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots. Journal of Open Source Software, 3(25), 731, https://doi.org/10.21105/joss.00731.

colnames(F0_Meta_F2_4Models_withRnor6Coord)

head(F0_Meta_F2_4Models_withRnor6Coord[,c(214:222)])

# GENENAME SYMBOL ENTREZID SEQCOORDSYSTEM SEQNAME SEQSTRAND GENESEQSTART GENESEQEND    GENEBIOTYPE
# 1    Lrp11  Lrp11   292462     chromosome       1         1      1702696    1731210 protein_coding
# 2    Pcmt1  Pcmt1       NA     chromosome       1        -1      1736276    1767618 protein_coding
# 3    Nup43  Nup43   683983     chromosome       1         1      1771710    1782091 protein_coding
# 4    Lats1  Lats1   308265     chromosome       1         1      1784078    1817310 protein_coding
# 5   Katna1 Katna1   292464     chromosome       1         1      1825872    1867787 protein_coding
# 6    Ginm1  Ginm1   361448     chromosome       1        -1      1870388    1885981 protein_coding


#This is old code from when I was working with the Rnor6 annotation from the org.db package:

# #I wonder if this package can handle the - sign before the coordinates. If not, we can make a new column with abs()
# # I had to change the ChrLocChr to have the NAs be NA instead of <NA>: 
# # Error in manhattan(F0_Meta_F2_4Models_withRnor6Coord, chr = "CHRLOCCHR",  : 
# # CHRLOCCHR column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again.
#     
# is.numeric(F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR)
# #[1] FALSE
# is.character(F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR)
# #[1] TRUE
# 

table(F0_Meta_F2_4Models_withRnor6Coord$SEQNAME)
# 1             10             11             12             13             14             15 
# 1754           1153            355            472            443            479            395 

# 16             17             18             19              2             20              3 
# 404            421            374            423            925            398           1023 

# 4              5              6              7              8              9 AABR07024041.1 
# 777            899            629            906            812            577              1 

# AABR07024206.1 AABR07024291.1     KL567908.1     KL568122.1     KL568128.1     KL568199.1     KL568414.1 
# 1              2              1              1              1              1              1 

# MT              X              Y 
# 24            491              5 


F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric<-F0_Meta_F2_4Models_withRnor6Coord$SEQNAME

F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric[F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric=="X"]<-"21"

F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric[F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric=="AABR07024041.1 "]<-NA

F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric[F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric=="AABR07024206.1"]<-NA

F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric[F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric=="AABR07024291.1"]<-NA

F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric[F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric=="AABR07024041.1"]<-NA

F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric[F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric=="KL567908.1"]<-NA

F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric[F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric=="KL568122.1"]<-NA

F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric[F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric=="KL568128.1"]<-NA

F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric[F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric=="KL568199.1"]<-NA

F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric[F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric=="KL568414.1"]<-NA

F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric[F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric=="MT"]<-NA

F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric[F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric=="Y"]<-NA


table(F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric)
# 1   10   11   12   13   14   15   16   17   18   19    2   20   21    3    4    5    6    7    8    9 
# 1754 1153  355  472  443  479  395  404  421  374  423  925  398  491 1023  777  899  629  906  812  577 

head(F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric)
head(as.numeric(F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric))

tail(F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric)
tail(as.numeric(F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric))

F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric<-as.numeric(F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric)

colnames(F0_Meta_F2_4Models_withRnor6Coord)

#Code leftover from when I was working with the Rnor6 coordinates from the org.db package:

# is.numeric(F0_Meta_F2_4Models_withRnor6Coord$CHRLOC)
# #[1] TRUE
# 
# hist(F0_Meta_F2_4Models_withRnor6Coord$CHRLOC)
# 
# F0_Meta_F2_4Models_withRnor6Coord$CHRLOC_ABS<-abs(F0_Meta_F2_4Models_withRnor6Coord$CHRLOC)
# # Error in rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP, d$CHR,  : 
# #                                                             invalid 'times' value

F0_Meta_F2_4Models_withRnor6Coord$CHRLOC_ABS<-F0_Meta_F2_4Models_withRnor6Coord$GENESEQSTART

#hmm... maybe it is having issues with the Chr=NAs?

F0_Meta_F2_4Models_withRnor6Coord_noNA<-F0_Meta_F2_4Models_withRnor6Coord[is.na(F0_Meta_F2_4Models_withRnor6Coord$CHRLOCCHR_Numeric)==FALSE,]

#Error in plot.window(...) : need finite 'ylim' values

#hmmm.... there must be some weird p-values in there too.

hist(F0_Meta_F2_4Models_withRnor6Coord_noNA$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)
sum(is.na(F0_Meta_F2_4Models_withRnor6Coord_noNA$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore))
#[1] 477

head(F0_Meta_F2_4Models_withRnor6Coord_noNA[is.na(F0_Meta_F2_4Models_withRnor6Coord_noNA$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore),])
#Ah- that's probably the problem. Those would be genes found in the F0 dataset but not the F2s.

F0_Meta_F2_4Models_withRnor6Coord_noNA<-F0_Meta_F2_4Models_withRnor6Coord_noNA[is.na(F0_Meta_F2_4Models_withRnor6Coord_noNA$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore)==FALSE,]


manhattan(F0_Meta_F2_4Models_withRnor6Coord_noNA, chr="CHRLOCCHR_Numeric", bp="CHRLOC_ABS", snp="SYMBOL", p="M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore")

#Hurray - it worked!  Well that's fun...

manhattan(F0_Meta_F2_4Models_withRnor6Coord_noNA, chr="CHRLOCCHR_Numeric", bp="CHRLOC_ABS", snp="SYMBOL", p="M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore", annotatePval = 0.001)
#nifty.
#It would be nice if we could have it use ENSEMBL ID when there isn't a SYMBOL

#Trying to figure out what they are using as a symbol when there isn't a symbol:

sum(is.na(F0_Meta_F2_4Models_withRnor6Coord_noNA$SYMBOL))
#[1] 643

F0_Meta_F2_4Models_withRnor6Coord_noNA$SYMBOLorENSEMBL<-F0_Meta_F2_4Models_withRnor6Coord_noNA$SYMBOL

F0_Meta_F2_4Models_withRnor6Coord_noNA$SYMBOLorENSEMBL[is.na(F0_Meta_F2_4Models_withRnor6Coord_noNA$SYMBOLorENSEMBL)]<-F0_Meta_F2_4Models_withRnor6Coord_noNA$ENSEMBL[is.na(F0_Meta_F2_4Models_withRnor6Coord_noNA$SYMBOLorENSEMBL)]

manhattan(F0_Meta_F2_4Models_withRnor6Coord_noNA, chr="CHRLOCCHR_Numeric", bp="CHRLOC_ABS", snp="SYMBOLorENSEMBL", p="M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore", annotatePval = 0.001)
#Nice

#Also - calculating a traditional bonferroni-corrected genome-wide significance threshold:
nrow(F0_Meta_F2_4Models_withRnor6Coord_noNA)
#[1] 13633
0.05/13633

pdf("ManhattanPlot_F2Locomotor_RNASeq_Pval.pdf", width=16, height=5)
manhattan(F0_Meta_F2_4Models_withRnor6Coord_noNA, chr="CHRLOCCHR_Numeric", bp="CHRLOC_ABS", snp="SYMBOLorENSEMBL", p="M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore", annotatePval = 0.001, annotateTop=FALSE, suggestiveline=3, genomewideline=-log10(3.667571e-06))
dev.off()

pdf("ManhattanPlot_F2EPMOpenArm_RNASeq_Pval.pdf", width=16, height=5)
manhattan(F0_Meta_F2_4Models_withRnor6Coord_noNA, chr="CHRLOCCHR_Numeric", bp="CHRLOC_ABS", snp="SYMBOLorENSEMBL", p="P.value.EPM_Percent_Time_Open_Arm", annotatePval = 0.001, annotateTop=FALSE, suggestiveline=3, genomewideline=-log10(3.667571e-06))
dev.off()

pdf("ManhattanPlot_F2EPMDistance_RNASeq_Pval.pdf", width=16, height=5)
manhattan(F0_Meta_F2_4Models_withRnor6Coord_noNA, chr="CHRLOCCHR_Numeric", bp="CHRLOC_ABS", snp="SYMBOLorENSEMBL", p="P.value.EPM_DistanceTraveled", annotatePval = 0.001, annotateTop=FALSE, suggestiveline=3, genomewideline=-log10(3.667571e-06))
dev.off()

pdf("ManhattanPlot_F2STGT_GT_RNASeq_Pval.pdf", width=16, height=5)
manhattan(F0_Meta_F2_4Models_withRnor6Coord_noNA, chr="CHRLOCCHR_Numeric", bp="CHRLOC_ABS", snp="SYMBOLorENSEMBL", p="P.value.LearningClassification_AsFactorGT", annotatePval = 0.001, annotateTop=FALSE, suggestiveline=3, genomewideline=-log10(3.667571e-06))
dev.off()

pdf("ManhattanPlot_F2STGT_IN_RNASeq_Pval.pdf", width=16, height=5)
manhattan(F0_Meta_F2_4Models_withRnor6Coord_noNA, chr="CHRLOCCHR_Numeric", bp="CHRLOC_ABS", snp="SYMBOLorENSEMBL", p="P.value.LearningClassification_AsFactorIN", annotatePval = 0.001, annotateTop=FALSE, suggestiveline=3, genomewideline=-log10(3.667571e-06))
dev.off()

pdf("ManhattanPlot_F0_RNASeq_Pval.pdf", width=16, height=5)
manhattan(F0_Meta_F2_4Models_withRnor6Coord_noNA, chr="CHRLOCCHR_Numeric", bp="CHRLOC_ABS", snp="SYMBOLorENSEMBL", p="P.value.Lineage_AsFactorbLR", annotatePval = 0.0001, annotateTop=FALSE, suggestiveline=3, genomewideline=-log10(3.667571e-06))
dev.off()

F0_Meta_F2_4Models_withRnor6Coord_noNA_forMeta<-F0_Meta_F2_4Models_withRnor6Coord_noNA[is.na(F0_Meta_F2_4Models_withRnor6Coord_noNA$pval)==FALSE,]

#Also - calculating a traditional bonferroni-corrected genome-wide significance threshold:
nrow(F0_Meta_F2_4Models_withRnor6Coord_noNA_forMeta)
#[1] 11238
0.05/11238

pdf("ManhattanPlot_MetaAnalysis_Pval.pdf", width=16, height=5)
manhattan(F0_Meta_F2_4Models_withRnor6Coord_noNA_forMeta, chr="CHRLOCCHR_Numeric", bp="CHRLOC_ABS", snp="SYMBOLorENSEMBL", p="pval", annotatePval = 0.0001, annotateTop=FALSE, suggestiveline=3, genomewideline=-log10(4.44919e-06))
dev.off()

#I wonder why the F0s have so much fewer genes on chr 1 that are suggestive/sig.


#It would be nice to have some summary plots across all bHR/bLR lineage studies and behavioral traits to compare to the QTL plots. How hard would that be to do?

F0_Meta_F2_4Models_withRnor6Coord_noNA$MinHRLR_Pvalue<-apply(cbind(F0_Meta_F2_4Models_withRnor6Coord_noNA$P.value.Lineage_AsFactorbLR, F0_Meta_F2_4Models_withRnor6Coord_noNA$pval), 1, function(y) min(y,na.rm=TRUE))

pdf("ManhattanPlot_MinHRLR_Pval.pdf", width=16, height=5)
manhattan(F0_Meta_F2_4Models_withRnor6Coord_noNA, chr="CHRLOCCHR_Numeric", bp="CHRLOC_ABS", snp="SYMBOLorENSEMBL", p="MinHRLR_Pvalue", annotatePval = 0.0001, annotateTop=FALSE, suggestiveline=3, genomewideline=-log10(3.667571e-06))
dev.off()

#Can we do something similar across F2 behavioral traits?

F0_Meta_F2_4Models_withRnor6Coord_noNA$MinF2Behavior_Pvalue<-apply(cbind(F0_Meta_F2_4Models_withRnor6Coord_noNA$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore, F0_Meta_F2_4Models_withRnor6Coord_noNA$P.value.EPM_DistanceTraveled, F0_Meta_F2_4Models_withRnor6Coord_noNA$P.value.EPM_Percent_Time_Open_Arm, F0_Meta_F2_4Models_withRnor6Coord_noNA$P.value.LearningClassification_AsFactorGT, F0_Meta_F2_4Models_withRnor6Coord_noNA$P.value.LearningClassification_AsFactorIN), 1, function(y) min(y,na.rm=TRUE))

pdf("ManhattanPlot_MinF2Behavior_Pval.pdf", width=16, height=5)
manhattan(F0_Meta_F2_4Models_withRnor6Coord_noNA, chr="CHRLOCCHR_Numeric", bp="CHRLOC_ABS", snp="SYMBOLorENSEMBL", p="MinF2Behavior_Pvalue", annotatePval = 0.001, annotateTop=FALSE, suggestiveline=3, genomewideline=-log10(3.667571e-06))
dev.off()

#############################################

#Comparing the results with Apurva's QTLs:

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/F2_QTL_Output")

LocomotorQTL<-read.delim("u01_huda_akil_total_loco_score.loco.mlma", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(LocomotorQTL)
# 'data.frame':	4425349 obs. of  9 variables:
# $ Chr : int  1 1 1 1 1 1 1 1 1 1 ...
# $ SNP : chr  "chr1:36146" "chr1:38004" "chr1:49847" "chr1:49985" ...
# $ bp  : int  36146 38004 49847 49985 50070 50071 50156 57580 62485 74000 ...
# $ A1  : chr  "C" "G" "A" "T" ...
# $ A2  : chr  "A" "A" "G" "C" ...
# $ Freq: num  0.0295 0.0124 0.0373 0.0714 0.0233 ...
# $ b   : num  -0.1953 0.2471 -0.0862 -0.0375 -0.1165 ...
# $ se  : num  0.223 0.344 0.199 0.149 0.255 ...
# $ p   : num  0.381 0.472 0.666 0.801 0.648 ...

#Let's trim that down before attempting a manhattan plot so it doesn't choke:

LocomotorQTL_p10<-LocomotorQTL[LocomotorQTL$p<0.10,]

str(LocomotorQTL_p10)
# 'data.frame':	835234 obs. of  9 variables:
#   $ Chr : int  1 1 1 1 1 1 1 1 1 1 ...
# $ SNP : chr  "chr1:174701" "chr1:174865" "chr1:175664" "chr1:175665" ...
# $ bp  : int  174701 174865 175664 175665 177097 190568 191890 191991 192656 192855 ...
# $ A1  : chr  "G" "A" "T" "T" ...
# $ A2  : chr  "T" "C" "C" "C" ...
# $ Freq: num  0.0186 0.0714 0.0714 0.0714 0.0466 ...
# $ b   : num  0.465 -0.337 -0.374 -0.374 -0.38 ...
# $ se  : num  0.276 0.15 0.15 0.15 0.184 ...
# $ p   : num  0.0922 0.0244 0.0127 0.0127 0.0396 ...

#Still pretty big. Chop further?
rm(LocomotorQTL_p10)

LocomotorQTL_p05<-LocomotorQTL[LocomotorQTL$p<0.05,]
str(LocomotorQTL_p05)
# 'data.frame':	547489 obs. of  9 variables:
#   $ Chr : int  1 1 1 1 1 1 1 1 1 1 ...
# $ SNP : chr  "chr1:174865" "chr1:175664" "chr1:175665" "chr1:177097" ...
# $ bp  : int  174865 175664 175665 177097 190568 191890 191991 192656 192855 192858 ...
# $ A1  : chr  "A" "T" "T" "T" ...
# $ A2  : chr  "C" "C" "C" "A" ...
# $ Freq: num  0.0714 0.0714 0.0714 0.0466 0.045 ...
# $ b   : num  -0.337 -0.374 -0.374 -0.38 -0.604 ...
# $ se  : num  0.15 0.15 0.15 0.184 0.185 ...
# $ p   : num  0.02439 0.01268 0.01268 0.03956 0.00108 ...

#... a little more manageable. Let's try it:

# pdf("ManhattanPlot_F2Locomotor_QTL_Pval.pdf", width=16, height=5)
# manhattan(LocomotorQTL_p05, chr="Chr", bp="bp", snp="SNP", p="p", annotatePval = 0.0001, annotateTop=FALSE, suggestiveline=4, genomewideline=-log10(5e-8))
# dev.off()
# Error in rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP, d$CHR,  : 
#                                                             invalid 'times' value

sum(is.na(LocomotorQTL_p05$Chr))
#[1] 1922
#... and it is character instead of numeric

table(LocomotorQTL_p05$Chr)

# 1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20 
# 150704  13457  48316  15126   7156  47267  63489  15421  15381   8121  23215   4820   6090  10516  12733  22935  22601  38576  16666   2977 

#So the problem is just the NAs

LocomotorQTL_p05_noNA<-LocomotorQTL_p05[is.na(LocomotorQTL_p05$Chr)==FALSE,]

LocomotorQTL_p05_noNA$Chr_Numeric<-as.numeric(LocomotorQTL_p05_noNA$Chr)

#I'm going to change the width because it doesn't have an x chromosome
#I also took out the annotate argument because it is too much on this plot

#Threshold for genome-wide significance, using bonferonni correction with nrow for entire dataframe of p-value output (not sure if some of these SNPs were thrown out in the end, so I was conservative)
0.05/4425349
#[1] 1.129854e-08

pdf("ManhattanPlot_F2Locomotor_QTL_Pval.pdf", width=15, height=5)
manhattan(LocomotorQTL_p05_noNA, chr="Chr_Numeric", bp="bp", snp="SNP", p="p", suggestiveline=4, genomewideline=-log10(1.129854e-08))
dev.off()

#Oh interesting - the spacing between chromosomes in the manhattan plots isn't equivalent for these different datasets - it must not be based on actual chromosome length but on the number of SNPs sequenced in it (or transcripts from it, or filtered SNPs, in this case)

#Placing the values in MB bins:

LocomotorQTL$bp_bin<-round(LocomotorQTL$bp, digits=-6)

head(LocomotorQTL$bp_bin)
#[1] 0 0 0 0 0 0
tail(LocomotorQTL$bp_bin)
#[1] 5.6e+07 5.6e+07 5.6e+07 5.6e+07 5.6e+07 5.6e+07

LocomotorQTL_noNA<-LocomotorQTL[is.na(LocomotorQTL$Chr)==FALSE,]

#Grabbing the minimum p-value within each bin:

LocomotorQTL_noNA$Chr_bpBin<-paste(LocomotorQTL_noNA$Chr, LocomotorQTL_noNA$bp_bin, sep="_")

head(LocomotorQTL_noNA$Chr_bpBin)

LocomotorQTL_noNA_MinPperMB<-tapply(LocomotorQTL_noNA$p, LocomotorQTL_noNA$Chr_bpBin, min)

head(LocomotorQTL_noNA_MinPperMB)
# 1_0  1_1.01e+08  1_1.02e+08  1_1.03e+08  1_1.04e+08  1_1.05e+08 
# 6.31605e-04 3.57434e-01 3.20530e-02 2.82908e-11 5.30348e-09 3.99753e-07 

#Ah. That format is going to be hard to work with.

strsplit(names(LocomotorQTL_noNA_MinPperMB)[1],"_")[[1]]
#[1] "1" "0"

LocomotorQTL_noNA_MinPperMB_Coordinates<-apply(as.matrix(names(LocomotorQTL_noNA_MinPperMB), length(names(LocomotorQTL_noNA_MinPperMB)), 1), 1, function(y) strsplit(y,"_")[[1]])

str(LocomotorQTL_noNA_MinPperMB_Coordinates)   
#chr [1:2, 1:2637] "1" "0" "1" "1.01e+08" "1" "1.02e+08" "1" "1.03e+08" "1" "1.04e+08" "1" "1.05e+08" "1" ...

LocomotorQTL_noNA_MinPperMB_Coordinates[,c(1:5)]
# [,1] [,2]       [,3]       [,4]       [,5]      
# [1,] "1"  "1"        "1"        "1"        "1"       
# [2,] "0"  "1.01e+08" "1.02e+08" "1.03e+08" "1.04e+08"

LocomotorQTL_noNA_MinPperMB_Coordinates<-t(LocomotorQTL_noNA_MinPperMB_Coordinates)

LocomotorQTL_noNA_MinPperMB_DF<-data.frame("Chr"=as.numeric(LocomotorQTL_noNA_MinPperMB_Coordinates[,1]), "MB_Bin"=as.numeric(LocomotorQTL_noNA_MinPperMB_Coordinates[,2]), "MinP"=unname(LocomotorQTL_noNA_MinPperMB ))

str(LocomotorQTL_noNA_MinPperMB_DF)

# 'data.frame':	2637 obs. of  3 variables:
# $ Chr   : num  1 1 1 1 1 1 1 1 1 1 ...
# $ MB_Bin: num  0.00 1.01e+08 1.02e+08 1.03e+08 1.04e+08 1.05e+08 1.06e+08 1.07e+08 1.08e+08 1.09e+08 ...
# $ MinP  : num [1:2637(1d)] 6.32e-04 3.57e-01 3.21e-02 2.83e-11 5.30e-09 ...

write.csv(LocomotorQTL_noNA_MinPperMB_DF, "LocomotorQTL_noNA_MinPperMB_DF.csv")

#Let's get rid of some of those big memory-hogging objects:
rm(LocomotorQTL, LocomotorQTL_p05, LocomotorQTL_p05_noNA)

#####

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/F2_QTL_Output")

EPMDistanceQTL<-read.delim("u01_huda_akil_epm_distance_traveled.loco.mlma", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(EPMDistanceQTL)
# 'data.frame':	4425349 obs. of  9 variables:
# $ Chr : int  1 1 1 1 1 1 1 1 1 1 ...
# $ SNP : chr  "chr1:36146" "chr1:38004" "chr1:49847" "chr1:49985" ...
# $ bp  : int  36146 38004 49847 49985 50070 50071 50156 57580 62485 74000 ...
# $ A1  : chr  "C" "G" "A" "T" ...
# $ A2  : chr  "A" "A" "G" "C" ...
# $ Freq: num  0.0295 0.0124 0.0373 0.0714 0.0233 ...
# $ b   : num  -0.1953 0.2471 -0.0862 -0.0375 -0.1165 ...
# $ se  : num  0.223 0.344 0.199 0.149 0.255 ...
# $ p   : num  0.381 0.472 0.666 0.801 0.648 ...

#This is the same size as the locomotor QTL data.frame, so I'm going to guess it has identical coordinates.

#Let's trim that down before attempting a manhattan plot so it doesn't choke:

EPMDistanceQTL_p05<-EPMDistanceQTL[EPMDistanceQTL$p<0.05,]
str(EPMDistanceQTL_p05)
# 'data.frame':	359828 obs. of  9 variables:
#   $ Chr : int  1 1 1 1 1 1 1 1 1 1 ...
# $ SNP : chr  "chr1:400732" "chr1:400752" "chr1:1388525" "chr1:1461883" ...
# $ bp  : int  400732 400752 1388525 1461883 1461885 1479105 1481441 1483277 1537027 3013342 ...
# $ A1  : chr  "G" "T" "A" "G" ...
# $ A2  : chr  "C" "C" "G" "T" ...
# $ Freq: num  0.0326 0.0326 0.045 0.0373 0.0373 ...
# $ b   : num  0.491 0.491 0.515 0.419 0.419 ...
# $ se  : num  0.212 0.212 0.185 0.194 0.194 ...
# $ p   : num  0.02046 0.02046 0.00527 0.03096 0.03096 ...

sum(is.na(EPMDistanceQTL_p05$Chr))
#[1] 1922
#... and it is character instead of numeric

table(EPMDistanceQTL_p05$Chr)

# 1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19 
# 69857 49137 16724  5152 14459 27970 29135  7548 29812  8833 11526  1306 14578 13694  6377 10163  4529 10738 10978 
# 20 
# 15390 

#So the problem is just the NAs

EPMDistanceQTL_p05_noNA<-EPMDistanceQTL_p05[is.na(EPMDistanceQTL_p05$Chr)==FALSE,]

EPMDistanceQTL_p05_noNA$Chr_Numeric<-as.numeric(EPMDistanceQTL_p05_noNA$Chr)

#I'm going to change the width because it doesn't have an x chromosome
#I also took out the annotate argument because it is too much on this plot

#Threshold for genome-wide significance, using bonferonni correction with nrow for entire dataframe of p-value output (not sure if some of these SNPs were thrown out in the end, so I was conservative)
0.05/4425349
#[1] 1.129854e-08

pdf("ManhattanPlot_F2EPMDistance_QTL_Pval.pdf", width=15, height=5)
manhattan(EPMDistanceQTL_p05_noNA, chr="Chr_Numeric", bp="bp", snp="SNP", p="p", suggestiveline=4, genomewideline=-log10(1.129854e-08))
dev.off()

#Placing the values in MB bins:

EPMDistanceQTL$bp_bin<-round(EPMDistanceQTL$bp, digits=-6)

EPMDistanceQTL_noNA<-EPMDistanceQTL[is.na(EPMDistanceQTL$Chr)==FALSE,]

#Grabbing the minimum p-value within each bin:

EPMDistanceQTL_noNA$Chr_bpBin<-paste(EPMDistanceQTL_noNA$Chr, EPMDistanceQTL_noNA$bp_bin, sep="_")

head(EPMDistanceQTL_noNA$Chr_bpBin)

EPMDistanceQTL_noNA_MinPperMB<-tapply(EPMDistanceQTL_noNA$p, EPMDistanceQTL_noNA$Chr_bpBin, min)

head(EPMDistanceQTL_noNA_MinPperMB)
# 1_0  1_1.01e+08  1_1.02e+08  1_1.03e+08  1_1.04e+08  1_1.05e+08 
# 0.020462700 0.029647700 0.036828100 0.000498216 0.000971501 0.001807900 

#Ah. That format is going to be hard to work with.

EPMDistanceQTL_noNA_MinPperMB_Coordinates<-apply(as.matrix(names(EPMDistanceQTL_noNA_MinPperMB), length(names(EPMDistanceQTL_noNA_MinPperMB)), 1), 1, function(y) strsplit(y,"_")[[1]])

str(EPMDistanceQTL_noNA_MinPperMB_Coordinates)   
# chr [1:2, 1:2637] "1" "0" "1" "1.01e+08" "1" "1.02e+08" "1" "1.03e+08" "1" "1.04e+08" "1" "1.05e+08" "1" ...

EPMDistanceQTL_noNA_MinPperMB_Coordinates[,c(1:5)]
# [,1] [,2]       [,3]       [,4]       [,5]      
# [1,] "1"  "1"        "1"        "1"        "1"       
# [2,] "0"  "1.01e+08" "1.02e+08" "1.03e+08" "1.04e+08"

EPMDistanceQTL_noNA_MinPperMB_Coordinates<-t(EPMDistanceQTL_noNA_MinPperMB_Coordinates)

EPMDistanceQTL_noNA_MinPperMB_DF<-data.frame("Chr"=as.numeric(EPMDistanceQTL_noNA_MinPperMB_Coordinates[,1]), "MB_Bin"=as.numeric(EPMDistanceQTL_noNA_MinPperMB_Coordinates[,2]), "MinP"=unname(EPMDistanceQTL_noNA_MinPperMB ))

str(EPMDistanceQTL_noNA_MinPperMB_DF)
# 'data.frame':	2637 obs. of  3 variables:
#   $ Chr   : num  1 1 1 1 1 1 1 1 1 1 ...
# $ MB_Bin: num  0.00 1.01e+08 1.02e+08 1.03e+08 1.04e+08 1.05e+08 1.06e+08 1.07e+08 1.08e+08 1.09e+08 ...
# $ MinP  : num [1:2637(1d)] 0.020463 0.029648 0.036828 0.000498 0.000972 ...

write.csv(EPMDistanceQTL_noNA_MinPperMB_DF, "EPMDistanceQTL_noNA_MinPperMB_DF.csv")

#Let's get rid of some of those big memory-hogging objects:
rm(EPMDistanceQTL, EPMDistanceQTL_p05, EPMDistanceQTL_p05_noNA)

##############

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/F2_QTL_Output")

EPMTimeOpenQTL<-read.delim("u01_huda_akil_epm_percent_time_open_arm.loco.mlma", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(EPMTimeOpenQTL)
# 'data.frame':	4425349 obs. of  9 variables:
#   $ Chr : int  1 1 1 1 1 1 1 1 1 1 ...
# $ SNP : chr  "chr1:36146" "chr1:38004" "chr1:49847" "chr1:49985" ...
# $ bp  : int  36146 38004 49847 49985 50070 50071 50156 57580 62485 74000 ...
# $ A1  : chr  "C" "G" "A" "T" ...
# $ A2  : chr  "A" "A" "G" "C" ...
# $ Freq: num  0.0295 0.0124 0.0373 0.0714 0.0233 ...
# $ b   : num  0.2017 0.0266 -0.3643 -0.1112 -0.1474 ...
# $ se  : num  0.223 0.341 0.199 0.148 0.252 ...
# $ p   : num  0.3656 0.9378 0.0674 0.4525 0.5589 ..

#This is the same size as the locomotor QTL data.frame, so I'm going to guess it has identical coordinates.

#Let's trim that down before attempting a manhattan plot so it doesn't choke:

EPMTimeOpenQTL_p05<-EPMTimeOpenQTL[EPMTimeOpenQTL$p<0.05,]
str(EPMTimeOpenQTL_p05)
# 'data.frame':	268767 obs. of  9 variables:
# $ Chr : int  1 1 1 1 1 1 1 1 1 1 ...
# $ SNP : chr  "chr1:356244" "chr1:467068" "chr1:485124" "chr1:485670" ...
# $ bp  : int  356244 467068 485124 485670 665190 667703 670308 1427764 1464035 2152041 ...
# $ A1  : chr  "C" "G" "G" "T" ...
# $ A2  : chr  "T" "T" "C" "C" ...
# $ Freq: num  0.0714 0.0186 0.0326 0.0248 0.0575 ...
# $ b   : num  0.43 -0.669 -0.465 -0.65 -0.467 ...
# $ se  : num  0.15 0.278 0.213 0.243 0.165 ...
# $ p   : num  0.00416 0.016 0.02897 0.00744 0.00466 ...

sum(is.na(EPMTimeOpenQTL_p05$Chr))
#[1] 1922
#... and it is character instead of numeric

table(EPMTimeOpenQTL_p05$Chr)

# 1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19 
# 28807 22204 15333 17680  5966  9405 15263  8327  9755 18502  7104  1792  6779  7840  9752 42225  4636 20692  1916 
# 20 
# 12867 

#So the problem is just the NAs

EPMTimeOpenQTL_p05_noNA<-EPMTimeOpenQTL_p05[is.na(EPMTimeOpenQTL_p05$Chr)==FALSE,]

EPMTimeOpenQTL_p05_noNA$Chr_Numeric<-as.numeric(EPMTimeOpenQTL_p05_noNA$Chr)

#I'm going to change the width because it doesn't have an x chromosome
#I also took out the annotate argument because it is too much on this plot

#Threshold for genome-wide significance, using bonferonni correction with nrow for entire dataframe of p-value output (not sure if some of these SNPs were thrown out in the end, so I was conservative)
0.05/4425349
#[1] 1.129854e-08

pdf("ManhattanPlot_F2EPMTimeOpen_QTL_Pval.pdf", width=15, height=5)
manhattan(EPMTimeOpenQTL_p05_noNA, chr="Chr_Numeric", bp="bp", snp="SNP", p="p", suggestiveline=4, genomewideline=-log10(1.129854e-08))
dev.off()

#Placing the values in MB bins:

EPMTimeOpenQTL$bp_bin<-round(EPMTimeOpenQTL$bp, digits=-6)

EPMTimeOpenQTL_noNA<-EPMTimeOpenQTL[is.na(EPMTimeOpenQTL$Chr)==FALSE,]

#Grabbing the minimum p-value within each bin:

EPMTimeOpenQTL_noNA$Chr_bpBin<-paste(EPMTimeOpenQTL_noNA$Chr, EPMTimeOpenQTL_noNA$bp_bin, sep="_")

head(EPMTimeOpenQTL_noNA$Chr_bpBin)

EPMTimeOpenQTL_noNA_MinPperMB<-tapply(EPMTimeOpenQTL_noNA$p, EPMTimeOpenQTL_noNA$Chr_bpBin, min)

head(EPMTimeOpenQTL_noNA_MinPperMB)
# 1_0 1_1.01e+08 1_1.02e+08 1_1.03e+08 1_1.04e+08 1_1.05e+08 
# 0.00416352 0.47838700 0.04512720 0.01489070 0.04975640 0.01027670 

#Ah. That format is going to be hard to work with.

EPMTimeOpenQTL_noNA_MinPperMB_Coordinates<-apply(as.matrix(names(EPMTimeOpenQTL_noNA_MinPperMB), length(names(EPMTimeOpenQTL_noNA_MinPperMB)), 1), 1, function(y) strsplit(y,"_")[[1]])

str(EPMTimeOpenQTL_noNA_MinPperMB_Coordinates)   
# chr [1:2, 1:2637] "1" "0" "1" "1.01e+08" "1" "1.02e+08" "1" "1.03e+08" "1" "1.04e+08" "1" "1.05e+08" "1" ...

EPMTimeOpenQTL_noNA_MinPperMB_Coordinates[,c(1:5)]
# [,1] [,2]       [,3]       [,4]       [,5]      
# [1,] "1"  "1"        "1"        "1"        "1"       
# [2,] "0"  "1.01e+08" "1.02e+08" "1.03e+08" "1.04e+08"

EPMTimeOpenQTL_noNA_MinPperMB_Coordinates<-t(EPMTimeOpenQTL_noNA_MinPperMB_Coordinates)

EPMTimeOpenQTL_noNA_MinPperMB_DF<-data.frame("Chr"=as.numeric(EPMTimeOpenQTL_noNA_MinPperMB_Coordinates[,1]), "MB_Bin"=as.numeric(EPMTimeOpenQTL_noNA_MinPperMB_Coordinates[,2]), "MinP"=unname(EPMTimeOpenQTL_noNA_MinPperMB ))

str(EPMTimeOpenQTL_noNA_MinPperMB_DF)
# 'data.frame':	2637 obs. of  3 variables:
#   $ Chr   : num  1 1 1 1 1 1 1 1 1 1 ...
# $ MB_Bin: num  0.00 1.01e+08 1.02e+08 1.03e+08 1.04e+08 1.05e+08 1.06e+08 1.07e+08 1.08e+08 1.09e+08 ...
# $ MinP  : num [1:2637(1d)] 0.00416 0.47839 0.04513 0.01489 0.04976 ...

write.csv(EPMTimeOpenQTL_noNA_MinPperMB_DF, "EPMTimeOpenQTL_noNA_MinPperMB_DF.csv")

#Let's get rid of some of those big memory-hogging objects:
rm(EPMTimeOpenQTL, EPMTimeOpenQTL_p05, EPMTimeOpenQTL_p05_noNA)

#####

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/F2_QTL_Output")

LateralLocoQTL<-read.delim("u01_huda_akil_lateral_loco_score.loco.mlma", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(LateralLocoQTL)
# 'data.frame':	4425349 obs. of  9 variables:
#   $ Chr : int  1 1 1 1 1 1 1 1 1 1 ...
# $ SNP : chr  "chr1:36146" "chr1:38004" "chr1:49847" "chr1:49985" ...
# $ bp  : int  36146 38004 49847 49985 50070 50071 50156 57580 62485 74000 ...
# $ A1  : chr  "C" "G" "A" "T" ...
# $ A2  : chr  "A" "A" "G" "C" ...
# $ Freq: num  0.0295 0.0124 0.0373 0.0714 0.0233 ...
# $ b   : num  -0.1788 0.3314 -0.0444 -0.024 -0.1379 ...
# $ se  : num  0.228 0.351 0.204 0.152 0.26 ...
# $ p   : num  0.433 0.345 0.828 0.874 0.596 ...

#This is the same size as the locomotor QTL data.frame, so I'm going to guess it has identical coordinates.

#Let's trim that down before attempting a manhattan plot so it doesn't choke:

LateralLocoQTL_p05<-LateralLocoQTL[LateralLocoQTL$p<0.05,]
str(LateralLocoQTL_p05)
# 'data.frame':	553114 obs. of  9 variables:
# $ Chr : int  1 1 1 1 1 1 1 1 1 1 ...
# $ SNP : chr  "chr1:174865" "chr1:175664" "chr1:175665" "chr1:177097" ...
# $ bp  : int  174865 175664 175665 177097 190568 191890 191991 192656 192855 192858 ...
# $ A1  : chr  "A" "T" "T" "T" ...
# $ A2  : chr  "C" "C" "C" "A" ...
# $ Freq: num  0.0714 0.0714 0.0714 0.0466 0.045 ...
# $ b   : num  -0.369 -0.409 -0.409 -0.456 -0.729 ...
# $ se  : num  0.153 0.153 0.153 0.188 0.189 ...
# $ p   : num  0.015905 0.007669 0.007669 0.015324 0.000115 ...

sum(is.na(LateralLocoQTL_p05$Chr))
#[1] 1922
#... and it is character instead of numeric

table(LateralLocoQTL_p05$Chr)

# 1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16 
# 150709  21441  52059  16741  13922  33633  75996  16632  14061   5406  24192   1949   3558   7671   8530  28858 
# 17     18     19     20 
# 23756  33411  17587   1080 

#So the problem is just the NAs

LateralLocoQTL_p05_noNA<-LateralLocoQTL_p05[is.na(LateralLocoQTL_p05$Chr)==FALSE,]

LateralLocoQTL_p05_noNA$Chr_Numeric<-as.numeric(LateralLocoQTL_p05_noNA$Chr)

#I'm going to change the width because it doesn't have an x chromosome
#I also took out the annotate argument because it is too much on this plot

#Threshold for genome-wide significance, using bonferonni correction with nrow for entire dataframe of p-value output (not sure if some of these SNPs were thrown out in the end, so I was conservative)
0.05/4425349
#[1] 1.129854e-08

pdf("ManhattanPlot_F2LateralLoco_QTL_Pval.pdf", width=15, height=5)
manhattan(LateralLocoQTL_p05_noNA, chr="Chr_Numeric", bp="bp", snp="SNP", p="p", suggestiveline=4, genomewideline=-log10(1.129854e-08))
dev.off()

#Placing the values in MB bins:

LateralLocoQTL$bp_bin<-round(LateralLocoQTL$bp, digits=-6)

LateralLocoQTL_noNA<-LateralLocoQTL[is.na(LateralLocoQTL$Chr)==FALSE,]

#Grabbing the minimum p-value within each bin:

LateralLocoQTL_noNA$Chr_bpBin<-paste(LateralLocoQTL_noNA$Chr, LateralLocoQTL_noNA$bp_bin, sep="_")

head(LateralLocoQTL_noNA$Chr_bpBin)

LateralLocoQTL_noNA_MinPperMB<-tapply(LateralLocoQTL_noNA$p, LateralLocoQTL_noNA$Chr_bpBin, min)

head(LateralLocoQTL_noNA_MinPperMB)
# 1_0  1_1.01e+08  1_1.02e+08  1_1.03e+08  1_1.04e+08  1_1.05e+08 
# 8.46175e-05 3.04111e-01 3.16711e-02 5.93743e-11 1.92671e-09 7.36782e-08 

#Ah. That format is going to be hard to work with.

LateralLocoQTL_noNA_MinPperMB_Coordinates<-apply(as.matrix(names(LateralLocoQTL_noNA_MinPperMB), length(names(LateralLocoQTL_noNA_MinPperMB)), 1), 1, function(y) strsplit(y,"_")[[1]])

str(LateralLocoQTL_noNA_MinPperMB_Coordinates)   
#chr [1:2, 1:2637] "1" "0" "1" "1.01e+08" "1" "1.02e+08" "1" "1.03e+08" "1" "1.04e+08" "1" "1.05e+08" "1" ...

LateralLocoQTL_noNA_MinPperMB_Coordinates[,c(1:5)]
#      [,1] [,2]       [,3]       [,4]       [,5]      
#[1,] "1"  "1"        "1"        "1"        "1"       
#[2,] "0"  "1.01e+08" "1.02e+08" "1.03e+08" "1.04e+08"

LateralLocoQTL_noNA_MinPperMB_Coordinates<-t(LateralLocoQTL_noNA_MinPperMB_Coordinates)

LateralLocoQTL_noNA_MinPperMB_DF<-data.frame("Chr"=as.numeric(LateralLocoQTL_noNA_MinPperMB_Coordinates[,1]), "MB_Bin"=as.numeric(LateralLocoQTL_noNA_MinPperMB_Coordinates[,2]), "MinP"=unname(LateralLocoQTL_noNA_MinPperMB ))

str(LateralLocoQTL_noNA_MinPperMB_DF)
# 'data.frame':	2637 obs. of  3 variables:
#   $ Chr   : num  1 1 1 1 1 1 1 1 1 1 ...
# $ MB_Bin: num  0.00 1.01e+08 1.02e+08 1.03e+08 1.04e+08 1.05e+08 1.06e+08 1.07e+08 1.08e+08 1.09e+08 ...
# $ MinP  : num [1:2637(1d)] 8.46e-05 3.04e-01 3.17e-02 5.94e-11 1.93e-09 ...

write.csv(LateralLocoQTL_noNA_MinPperMB_DF, "LateralLocoQTL_noNA_MinPperMB_DF.csv")

#Let's get rid of some of those big memory-hogging objects:
rm(LateralLocoQTL, LateralLocoQTL_p05, LateralLocoQTL_p05_noNA)


#####

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/F2_QTL_Output")

RearingLocoQTL<-read.delim("u01_huda_akil_rearing_loco_score.loco.mlma", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(RearingLocoQTL)
# 'data.frame':	4425349 obs. of  9 variables:
# $ Chr : int  1 1 1 1 1 1 1 1 1 1 ...
# $ SNP : chr  "chr1:36146" "chr1:38004" "chr1:49847" "chr1:49985" ...
# $ bp  : int  36146 38004 49847 49985 50070 50071 50156 57580 62485 74000 ...
# $ A1  : chr  "C" "G" "A" "T" ...
# $ A2  : chr  "A" "A" "G" "C" ...
# $ Freq: num  0.0295 0.0124 0.0373 0.0714 0.0233 ...
# $ b   : num  -0.1709 0.1626 -0.1019 -0.0383 -0.127 ...
# $ se  : num  0.218 0.335 0.195 0.145 0.249 ...
# $ p   : num  0.433 0.628 0.601 0.792 0.61 ...

#This is the same size as the locomotor QTL data.frame, so I'm going to guess it has identical coordinates.

#Let's trim that down before attempting a manhattan plot so it doesn't choke:

RearingLocoQTL_p05<-RearingLocoQTL[RearingLocoQTL$p<0.05,]
str(RearingLocoQTL_p05)
# 'data.frame':	538859 obs. of  9 variables:
#   $ Chr : int  1 1 1 1 1 1 1 1 1 1 ...
# $ SNP : chr  "chr1:174865" "chr1:175664" "chr1:175665" "chr1:190568" ...
# $ bp  : int  174865 175664 175665 190568 191890 191991 192656 192855 192858 192910 ...
# $ A1  : chr  "A" "T" "T" "T" ...
# $ A2  : chr  "C" "C" "C" "A" ...
# $ Freq: num  0.0714 0.0714 0.0714 0.045 0.0435 ...
# $ b   : num  -0.295 -0.334 -0.334 -0.501 -0.498 ...
# $ se  : num  0.146 0.147 0.147 0.18 0.184 ...
# $ p   : num  0.04346 0.02262 0.02262 0.00546 0.00677 ...

sum(is.na(RearingLocoQTL_p05$Chr))
#[1] 1922
#... and it is character instead of numeric

table(RearingLocoQTL_p05$Chr)

# 1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16 
# 143983  12704  44078  14806   5111  51383  50964  13137  16073  19697  22209   5947   5519  14039  14724  17555 
# 17     18     19     20 
# 23497  41178  16362   3971 

#So the problem is just the NAs

RearingLocoQTL_p05_noNA<-RearingLocoQTL_p05[is.na(RearingLocoQTL_p05$Chr)==FALSE,]

RearingLocoQTL_p05_noNA$Chr_Numeric<-as.numeric(RearingLocoQTL_p05_noNA$Chr)

#I'm going to change the width because it doesn't have an x chromosome
#I also took out the annotate argument because it is too much on this plot

#Threshold for genome-wide significance, using bonferonni correction with nrow for entire dataframe of p-value output (not sure if some of these SNPs were thrown out in the end, so I was conservative)
0.05/4425349
#[1] 1.129854e-08

pdf("ManhattanPlot_F2RearingLoco_QTL_Pval.pdf", width=15, height=5)
manhattan(RearingLocoQTL_p05_noNA, chr="Chr_Numeric", bp="bp", snp="SNP", p="p", suggestiveline=4, genomewideline=-log10(1.129854e-08))
dev.off()

#Placing the values in MB bins:

RearingLocoQTL$bp_bin<-round(RearingLocoQTL$bp, digits=-6)

RearingLocoQTL_noNA<-RearingLocoQTL[is.na(RearingLocoQTL$Chr)==FALSE,]

#Grabbing the minimum p-value within each bin:

RearingLocoQTL_noNA$Chr_bpBin<-paste(RearingLocoQTL_noNA$Chr, RearingLocoQTL_noNA$bp_bin, sep="_")

head(RearingLocoQTL_noNA$Chr_bpBin)

RearingLocoQTL_noNA_MinPperMB<-tapply(RearingLocoQTL_noNA$p, RearingLocoQTL_noNA$Chr_bpBin, min)

head(RearingLocoQTL_noNA_MinPperMB)
# 1_0  1_1.01e+08  1_1.02e+08  1_1.03e+08  1_1.04e+08  1_1.05e+08 
# 3.25010e-03 3.01747e-01 3.09003e-02 1.33520e-10 7.15045e-08 1.99862e-06 

#Ah. That format is going to be hard to work with.

RearingLocoQTL_noNA_MinPperMB_Coordinates<-apply(as.matrix(names(RearingLocoQTL_noNA_MinPperMB), length(names(RearingLocoQTL_noNA_MinPperMB)), 1), 1, function(y) strsplit(y,"_")[[1]])

str(RearingLocoQTL_noNA_MinPperMB_Coordinates)   
#chr [1:2, 1:2637] "1" "0" "1" "1.01e+08" "1" "1.02e+08" "1" "1.03e+08" "1" "1.04e+08" "1" "1.05e+08" "1" ...

RearingLocoQTL_noNA_MinPperMB_Coordinates[,c(1:5)]
#      [,1] [,2]       [,3]       [,4]       [,5]      
#[1,] "1"  "1"        "1"        "1"        "1"       
#[2,] "0"  "1.01e+08" "1.02e+08" "1.03e+08" "1.04e+08"

RearingLocoQTL_noNA_MinPperMB_Coordinates<-t(RearingLocoQTL_noNA_MinPperMB_Coordinates)

RearingLocoQTL_noNA_MinPperMB_DF<-data.frame("Chr"=as.numeric(RearingLocoQTL_noNA_MinPperMB_Coordinates[,1]), "MB_Bin"=as.numeric(RearingLocoQTL_noNA_MinPperMB_Coordinates[,2]), "MinP"=unname(RearingLocoQTL_noNA_MinPperMB ))

str(RearingLocoQTL_noNA_MinPperMB_DF)
# 'data.frame':	2637 obs. of  3 variables:
#   $ Chr   : num  1 1 1 1 1 1 1 1 1 1 ...
# $ MB_Bin: num  0.00 1.01e+08 1.02e+08 1.03e+08 1.04e+08 1.05e+08 1.06e+08 1.07e+08 1.08e+08 1.09e+08 ...
# $ MinP  : num [1:2637(1d)] 3.25e-03 3.02e-01 3.09e-02 1.34e-10 7.15e-08 ...

write.csv(RearingLocoQTL_noNA_MinPperMB_DF, "RearingLocoQTL_noNA_MinPperMB_DF.csv")

#Let's get rid of some of those big memory-hogging objects:
rm(RearingLocoQTL, RearingLocoQTL_p05, RearingLocoQTL_p05_noNA)

########

#I'm not sure whether I should include EPM boli - that is one of the QTLs that Apurva identified, but the data is super messy (and non-normal) and we didn't use it for our own RNA-Seq analyses


#####

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/F2_QTL_Output")

EPMBoliQTL<-read.delim("u01_huda_akil_epm_percent_time_open_arm.loco.mlma", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(EPMBoliQTL)
# 'data.frame':	4425349 obs. of  9 variables:
#   $ Chr : int  1 1 1 1 1 1 1 1 1 1 ...
# $ SNP : chr  "chr1:36146" "chr1:38004" "chr1:49847" "chr1:49985" ...
# $ bp  : int  36146 38004 49847 49985 50070 50071 50156 57580 62485 74000 ...
# $ A1  : chr  "C" "G" "A" "T" ...
# $ A2  : chr  "A" "A" "G" "C" ...
# $ Freq: num  0.0295 0.0124 0.0373 0.0714 0.0233 ...
# $ b   : num  0.2017 0.0266 -0.3643 -0.1112 -0.1474 ...
# $ se  : num  0.223 0.341 0.199 0.148 0.252 ...
# $ p   : num  0.3656 0.9378 0.0674 0.4525 0.5589 ...

#This is the same size as the locomotor QTL data.frame, so I'm going to guess it has identical coordinates.

#Let's trim that down before attempting a manhattan plot so it doesn't choke:

EPMBoliQTL_p05<-EPMBoliQTL[EPMBoliQTL$p<0.05,]
str(EPMBoliQTL_p05)
# 'data.frame':	268767 obs. of  9 variables:
#   $ Chr : int  1 1 1 1 1 1 1 1 1 1 ...
# $ SNP : chr  "chr1:356244" "chr1:467068" "chr1:485124" "chr1:485670" ...
# $ bp  : int  356244 467068 485124 485670 665190 667703 670308 1427764 1464035 2152041 ...
# $ A1  : chr  "C" "G" "G" "T" ...
# $ A2  : chr  "T" "T" "C" "C" ...
# $ Freq: num  0.0714 0.0186 0.0326 0.0248 0.0575 ...
# $ b   : num  0.43 -0.669 -0.465 -0.65 -0.467 ...
# $ se  : num  0.15 0.278 0.213 0.243 0.165 ...
# $ p   : num  0.00416 0.016 0.02897 0.00744 0.00466 ...

sum(is.na(EPMBoliQTL_p05$Chr))
#[1] 1922
#... and it is character instead of numeric

table(EPMBoliQTL_p05$Chr)

# 1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19 
# 28807 22204 15333 17680  5966  9405 15263  8327  9755 18502  7104  1792  6779  7840  9752 42225  4636 20692  1916 
# 20 
# 12867 

#So the problem is just the NAs

EPMBoliQTL_p05_noNA<-EPMBoliQTL_p05[is.na(EPMBoliQTL_p05$Chr)==FALSE,]

EPMBoliQTL_p05_noNA$Chr_Numeric<-as.numeric(EPMBoliQTL_p05_noNA$Chr)

#I'm going to change the width because it doesn't have an x chromosome
#I also took out the annotate argument because it is too much on this plot

#Threshold for genome-wide significance, using bonferonni correction with nrow for entire dataframe of p-value output (not sure if some of these SNPs were thrown out in the end, so I was conservative)
0.05/4425349
#[1] 1.129854e-08

pdf("ManhattanPlot_F2EPMBoli_QTL_Pval.pdf", width=15, height=5)
manhattan(EPMBoliQTL_p05_noNA, chr="Chr_Numeric", bp="bp", snp="SNP", p="p", suggestiveline=4, genomewideline=-log10(1.129854e-08))
dev.off()

#Placing the values in MB bins:

EPMBoliQTL$bp_bin<-round(EPMBoliQTL$bp, digits=-6)

EPMBoliQTL_noNA<-EPMBoliQTL[is.na(EPMBoliQTL$Chr)==FALSE,]

#Grabbing the minimum p-value within each bin:

EPMBoliQTL_noNA$Chr_bpBin<-paste(EPMBoliQTL_noNA$Chr, EPMBoliQTL_noNA$bp_bin, sep="_")

head(EPMBoliQTL_noNA$Chr_bpBin)

EPMBoliQTL_noNA_MinPperMB<-tapply(EPMBoliQTL_noNA$p, EPMBoliQTL_noNA$Chr_bpBin, min)

head(EPMBoliQTL_noNA_MinPperMB)
#        1_0 1_1.01e+08 1_1.02e+08 1_1.03e+08 1_1.04e+08 1_1.05e+08 
# 0.00416352 0.47838700 0.04512720 0.01489070 0.04975640 0.01027670  

#Ah. That format is going to be hard to work with.

EPMBoliQTL_noNA_MinPperMB_Coordinates<-apply(as.matrix(names(EPMBoliQTL_noNA_MinPperMB), length(names(EPMBoliQTL_noNA_MinPperMB)), 1), 1, function(y) strsplit(y,"_")[[1]])

str(EPMBoliQTL_noNA_MinPperMB_Coordinates)   
#chr [1:2, 1:2637] "1" "0" "1" "1.01e+08" "1" "1.02e+08" "1" "1.03e+08" "1" "1.04e+08" "1" "1.05e+08" "1" ...

EPMBoliQTL_noNA_MinPperMB_Coordinates[,c(1:5)]
#      [,1] [,2]       [,3]       [,4]       [,5]      
#[1,] "1"  "1"        "1"        "1"        "1"       
#[2,] "0"  "1.01e+08" "1.02e+08" "1.03e+08" "1.04e+08"

EPMBoliQTL_noNA_MinPperMB_Coordinates<-t(EPMBoliQTL_noNA_MinPperMB_Coordinates)

EPMBoliQTL_noNA_MinPperMB_DF<-data.frame("Chr"=as.numeric(EPMBoliQTL_noNA_MinPperMB_Coordinates[,1]), "MB_Bin"=as.numeric(EPMBoliQTL_noNA_MinPperMB_Coordinates[,2]), "MinP"=unname(EPMBoliQTL_noNA_MinPperMB ))

str(EPMBoliQTL_noNA_MinPperMB_DF)
# 'data.frame':	2637 obs. of  3 variables:
# $ Chr   : num  1 1 1 1 1 1 1 1 1 1 ...
# $ MB_Bin: num  0.00 1.01e+08 1.02e+08 1.03e+08 1.04e+08 1.05e+08 1.06e+08 1.07e+08 1.08e+08 1.09e+08 ...
# $ MinP  : num [1:2637(1d)] 0.00416 0.47839 0.04513 0.01489 0.04976 ...

write.csv(EPMBoliQTL_noNA_MinPperMB_DF, "EPMBoliQTL_noNA_MinPperMB_DF.csv")

#Let's get rid of some of those big memory-hogging objects:
rm(EPMBoliQTL, EPMBoliQTL_p05, EPMBoliQTL_p05_noNA)

###########################

#Doublechecking that the coordinates are the same across these QTL-related objects:

sum(EPMBoliQTL_noNA_MinPperMB_DF$MB_Bin!=LocomotorQTL_noNA_MinPperMB_DF$MB_Bin)
#[1] 0
sum(EPMBoliQTL_noNA_MinPperMB_DF$MB_Bin==LocomotorQTL_noNA_MinPperMB_DF$MB_Bin)
#[1] 2637

F2Behavior_QTLs_MinPperMBBin_DF<-cbind.data.frame("Chr"=LocomotorQTL_noNA_MinPperMB_DF$Chr, "MB_Bin"=LocomotorQTL_noNA_MinPperMB_DF$MB_Bin, "LocomotorMinP"=LocomotorQTL_noNA_MinPperMB_DF$MinP, "RearingMinP"=RearingLocoQTL_noNA_MinPperMB_DF$MinP, "LateralMinP"=LateralLocoQTL_noNA_MinPperMB_DF$MinP, "EPMDistanceMinP"=EPMDistanceQTL_noNA_MinPperMB_DF$MinP, "EPMTimeOpenMinP"=EPMTimeOpenQTL_noNA_MinPperMB_DF$MinP, "EPMBoli"=EPMBoliQTL_noNA_MinPperMB_DF$MinP)

write.csv(F2Behavior_QTLs_MinPperMBBin_DF, "F2Behavior_QTLs_MinPperMBBin_DF.csv")

str(F2Behavior_QTLs_MinPperMBBin_DF)
# 'data.frame':	2637 obs. of  8 variables:
#   $ Chr            : num  1 1 1 1 1 1 1 1 1 1 ...
# $ MB_Bin         : num  0.00 1.01e+08 1.02e+08 1.03e+08 1.04e+08 1.05e+08 1.06e+08 1.07e+08 1.08e+08 1.09e+08 ...
# $ LocomotorMinP  : num [1:2637(1d)] 6.32e-04 3.57e-01 3.21e-02 2.83e-11 5.30e-09 ...
# $ RearingMinP    : num [1:2637(1d)] 3.25e-03 3.02e-01 3.09e-02 1.34e-10 7.15e-08 ...
# $ LateralMinP    : num [1:2637(1d)] 8.46e-05 3.04e-01 3.17e-02 5.94e-11 1.93e-09 ...
# $ EPMDistanceMinP: num [1:2637(1d)] 0.020463 0.029648 0.036828 0.000498 0.000972 ...
# $ EPMTimeOpenMinP: num [1:2637(1d)] 0.00416 0.47839 0.04513 0.01489 0.04976 ...
# $ EPMBoli        : num [1:2637(1d)] 0.00416 0.47839 0.04513 0.01489 0.04976 ...

F2Behavior_QTLs_MinPperMBBin_DF$F2Behav_MinP<-apply(F2Behavior_QTLs_MinPperMBBin_DF[,c(3:8)], 1, function(y) min(y, na.rm=TRUE))
#some errors due to bins missing data: e.g.
#In min(y, na.rm = TRUE) : no non-missing arguments to min; returning Inf

sum(F2Behavior_QTLs_MinPperMBBin_DF$F2Behav_MinP=="Inf")
#[1] 33

F2Behavior_QTLs_MinPperMBBin_DF_NoNA<-F2Behavior_QTLs_MinPperMBBin_DF[F2Behavior_QTLs_MinPperMBBin_DF$F2Behav_MinP!="Inf",]

pdf("ManhattanPlot_F2Behav_QTL_MinPval.pdf", width=15, height=5)
manhattan(F2Behavior_QTLs_MinPperMBBin_DF_NoNA, chr="Chr", bp="MB_Bin", snp="Chr", p="F2Behav_MinP", suggestiveline=4, genomewideline=-log10(1.129854e-08))
dev.off()
#LOD 4 may be too liberal of a cut-off


################################################

#Aligning across the RNA-Seq and QTL datasets:

#First question:  I used +/-1MB for the alignment with Goldman LODs.  The coding for that is a little obnoxious, so I just opened up the "F2Behavior_QTLs_MinPperMBBin_DF.csv" document and did it in Excel (calculated min for 1 MB bin across all QTLs (ignoring NAs) and then determined the min across three 1MB bins (i.e.+/-1MB each direction) while being conscious of chromosome)

F2Behavior_QTLs_MinPperMBBin_DF_PlusMinus3MB<-read.csv("F2Behavior_QTLs_MinPperMBBin_DF_PlusMinus3MB.csv")
str(F2Behavior_QTLs_MinPperMBBin_DF_PlusMinus3MB)

# 'data.frame':	2637 obs. of  11 variables:
# $ X                : int  1 2 3 4 5 6 7 8 9 10 ...
# $ Chr              : int  1 1 1 1 1 1 1 1 1 1 ...
# $ MB_Bin           : num  0.00 1.01e+08 1.02e+08 1.03e+08 1.04e+08 1.05e+08 1.06e+08 1.07e+08 1.08e+08 1.09e+08 ...
# $ LocomotorMinP    : num  6.32e-04 3.57e-01 3.21e-02 2.83e-11 5.30e-09 ...
# $ RearingMinP      : num  3.25e-03 3.02e-01 3.09e-02 1.34e-10 7.15e-08 ...
# $ LateralMinP      : num  8.46e-05 3.04e-01 3.17e-02 5.94e-11 1.93e-09 ...
# $ EPMDistanceMinP  : num  0.020463 0.029648 0.036828 0.000498 0.000972 ...
# $ EPMTimeOpenMinP  : num  0.00416 0.47839 0.04513 0.01489 0.04976 ...
# $ EPMBoli          : num  0.00416 0.47839 0.04513 0.01489 0.04976 ...
# $ MinP_QTLs_1MB_Bin: num  8.46e-05 2.96e-02 3.09e-02 2.83e-11 1.93e-09 ...
# $ MinP_QTLs_3MB_Bin: num  8.46e-05 8.46e-05 2.83e-11 2.83e-11 2.83e-11 ...

F2Behavior_QTLs_MinPperMBBin_DF_PlusMinus3MB_NoNA<-F2Behavior_QTLs_MinPperMBBin_DF_PlusMinus3MB[is.na(F2Behavior_QTLs_MinPperMBBin_DF_PlusMinus3MB$MinP_QTLs_1MB_Bin)==FALSE,]

str(F2Behavior_QTLs_MinPperMBBin_DF_PlusMinus3MB_NoNA)
#'data.frame':	2604 obs. of  11 variables:

colnames(F0_Meta_F2_4Models_withRnor6Coord_noNA)

F0_Meta_F2_4Models_withRnor6Coord_noNA_ForJoining<-F0_Meta_F2_4Models_withRnor6Coord_noNA[,c(1:3,5,10,15,20,26:32,35,57,59,77,86,198:227)]

str(F0_Meta_F2_4Models_withRnor6Coord_noNA_ForJoining)

head(F0_Meta_F2_4Models_withRnor6Coord_noNA_ForJoining)

F0_Meta_F2_4Models_withRnor6Coord_noNA_ForJoining$bp_bin<-round(F0_Meta_F2_4Models_withRnor6Coord_noNA_ForJoining$CHRLOC_ABS, digits=-6)

head(F0_Meta_F2_4Models_withRnor6Coord_noNA_ForJoining$CHRLOC_ABS)
head(F0_Meta_F2_4Models_withRnor6Coord_noNA_ForJoining$bp_bin)
tail(F0_Meta_F2_4Models_withRnor6Coord_noNA_ForJoining$bp_bin)

F0_Meta_F2_4Models_withRnor6Coord_noNA_ForJoining$Chr_bpBin<-paste(F0_Meta_F2_4Models_withRnor6Coord_noNA_ForJoining$CHRLOCCHR_Numeric, F0_Meta_F2_4Models_withRnor6Coord_noNA_ForJoining$bp_bin, sep="_")

head(F0_Meta_F2_4Models_withRnor6Coord_noNA_ForJoining$Chr_bpBin)
#[1] "1_2e+06" "1_2e+06" "1_2e+06" "1_2e+06" "1_2e+06" "1_2e+06"

F2Behavior_QTLs_MinPperMBBin_DF_PlusMinus3MB_NoNA$Chr_bpBin<-paste(F2Behavior_QTLs_MinPperMBBin_DF_PlusMinus3MB_NoNA$Chr, F2Behavior_QTLs_MinPperMBBin_DF_NoNA$MB_Bin, sep="_")

head(F2Behavior_QTLs_MinPperMBBin_DF_NoNA$Chr_bpBin)
#[1] "1_0"        "1_1.01e+08" "1_1.02e+08" "1_1.03e+08" "1_1.04e+08" "1_1.05e+08"

#I wonder if this will join properly if the scientific notation was formatted differently in the genetics data.

length(F0_Meta_F2_4Models_withRnor6Coord_noNA_ForJoining$Chr_bpBin)
#[1] 13633

sum(F0_Meta_F2_4Models_withRnor6Coord_noNA_ForJoining$Chr_bpBin%in%F2Behavior_QTLs_MinPperMBBin_DF_PlusMinus3MB_NoNA$Chr_bpBin)
#[1] 13043

#hmmm... something like 600 of the chromosomal locations in the F0/F2 data aren't found in the F2 QTL data. Missing chromosome?

table(F0_Meta_F2_4Models_withRnor6Coord_noNA_ForJoining$CHRLOCCHR_Numeric)
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21 
#1412  773  873  655  745  521  773  679  461  960  298  387  382  391  326  329  309  293  336  346  353 

table(F2Behavior_QTLs_MinPperMBBin_DF_PlusMinus3MB_NoNA$Chr)
#  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
#280 267 178 183 175 147 145 132 123 108  91  53 111 113 112  91  92  86  60  57 

#Yep - the F0/F2 expression data includes X chromosome and the QTL data does not. I forgot about that.

#Good enough for government business.

library(plyr)

F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins<-join(F0_Meta_F2_4Models_withRnor6Coord_noNA_ForJoining, F2Behavior_QTLs_MinPperMBBin_DF_PlusMinus3MB_NoNA, by="Chr_bpBin", type="left")

#Earlier analysis: We're missing the minimum pval columns for bHR/bLR and F2 BehavRNASeq:
#Newer analysis: I'm just renaming them so the old code works (and the new names are better anyway)

colnames(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins)

#Old code:
# F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$HRLR_RNASeq_minP<-apply(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins[,c(6,12)],1,function(y) min(y, na.rm=TRUE))

F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$HRLR_RNASeq_minP<-F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$MinHRLR_Pvalue

#how many NAs?
sum(is.na(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$HRLR_RNASeq_minP))
#[1] 0

# F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$F2Behavior_RNASeq_minP<-apply(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins[,c(18,22,26,32,33)],1,function(y) min(y, na.rm=TRUE))

F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$F2Behavior_RNASeq_minP<-F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$MinF2Behavior_Pvalue

write.csv(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins, "F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins.csv")

F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$SYMBOLorENSEMBL[is.na(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$SYMBOL)]


#overlap between QTLs and F0/F2 RNASeq results in 1 MB bin:
F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$SYMBOLorENSEMBL[which(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$MinP_QTLs_1MB_Bin<0.0001 & F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$HRLR_RNASeq_minP<0.0001 & F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$F2Behavior_RNASeq_minP<0.05)]
# [1] "RGD1564801"         "ENSRNOG00000052237" "Mfge8"              "C2cd3"              "Plekhb1"           
# [6] "Rps4y2"             "Sp3"                "Ttc30a1"            "Ttc30a1"            "Mal2"              
# [11] "Rpl17"              "Ist1"               "Spg7"  


#overlap between QTLs and F0/F2 RNASeq results in 1 MB bin (+/1 1MB)
F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$SYMBOLorENSEMBL[which(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$MinP_QTLs_3MB_Bin<0.0001 & F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$HRLR_RNASeq_minP<0.0001 & F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$F2Behavior_RNASeq_minP<0.05)]
# [1] "RGD1564801"         "ENSRNOG00000052237" "Mfge8"              "C2cd3"              "Plekhb1"           
# [6] "Tssc4"              "Rps4y2"             "Sp3"                "Ttc30a1"            "Ttc30a1"           
# [11] "Mal2"               "Ghdc"               "LOC302192"          "Rpl17"              "Ist1"              
# [16] "Spg7"    

#Only adds a few more - "Tssc4", "Ghdc", "LOC302192"    

#Using LOD4, but a weaker cut-off for bHR/bLR
F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$SYMBOLorENSEMBL[which(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$MinP_QTLs_1MB_Bin<0.0001 & F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$HRLR_RNASeq_minP<0.001 & F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$F2Behavior_RNASeq_minP<0.05)]
# [1] "RGD1564801"         "ENSRNOG00000052237" "Iglon5"             "Mfge8"              "Pex11a"            
# [6] "Unc45a"             "C2cd3"              "Ucp2"               "Plekhb1"            "Rps4y2"            
# [11] "Sp3"                "Ttc30a1"            "Ttc30a1"            "Ilvbl"              "Mal2"              
# [16] "ENSRNOG00000049194" "Rpl17"              "Ist1"               "Spg7"               "Vps9d1"    


#Using LOD4 in +/-1MB bin in either direction, and a weaker cut-off for bHR/bLR: 
F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$SYMBOLorENSEMBL[which(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$MinP_QTLs_3MB_Bin<0.0001 & F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$HRLR_RNASeq_minP<0.001 & F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$F2Behavior_RNASeq_minP<0.05)]
# [1] "Zfp551"             "RGD1564801"         "ENSRNOG00000052237" "Iglon5"             "Mfge8"             
# [6] "Pex11a"             "Unc45a"             "C2cd3"              "Ucp2"               "Plekhb1"           
# [11] "Hsd3b7"             "Tgfb1i1"            "Tssc4"              "Rps4y2"             "Sp3"               
# [16] "Ttc30a1"            "Ttc30a1"            "Snhg11"             "Ilvbl"              "Mal2"              
# [21] "ENSRNOG00000049194" "Ghdc"               "Sema3g"             "LOC302192"          "Rpl17"             
# [26] "Ist1"               "Spg7"               "Vps9d1"

#Huh. Well that's interesting - we would probably need to manually double check which of these is *actually* within +/-1MB.

#Outputing the genes that are near a QTL and have at least both nominal bHR/bLR and nominal F2 behavior effects:
write.csv(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins[which(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$MinP_QTLs_3MB_Bin<0.0001 & F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$HRLR_RNASeq_minP<0.05 & F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$F2Behavior_RNASeq_minP<0.05),], "F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_TopGenes.csv")

#Outputing all genes (not just near a QTL) that have a significant bHR/bLR effect (FDR<0.10) or nominal bHR/bLR effects in both datasets (p<0.05) and have at least a nominal F2 behavior effects:
write.csv(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins[which(((F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$P.value.adj.Lineage_AsFactorbLR<0.10 & is.na(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$P.value.adj.Lineage_AsFactorbLR)==FALSE)|(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$BH<0.10 & is.na(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$BH)==FALSE) |(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$P.value.Lineage_AsFactorbLR<0.05 & F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$pval<0.05 & is.na(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$P.value.adj.Lineage_AsFactorbLR)==FALSE & is.na(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$BH)==FALSE)) & (F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$F2Behavior_RNASeq_minP<0.05 & is.na(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$F2Behavior_RNASeq_minP)==FALSE)),], "F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_TopGenes_F0andF2effects.csv")


#How many genes are near these QTL regions?

sum(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$MinP_QTLs_1MB_Bin<0.001, na.rm=TRUE)
#[1] 4616

nrow(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins)
#[1] 13633

4616/13633
#[1] 0.3385902

#So something like 1/3 of the genes in our studies are within 1 MB of LD3 for *some* bHR/bLR related behavior

sum(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$MinP_QTLs_1MB_Bin<0.0001, na.rm=TRUE)
#[1] 1332
#And around 1/10 within 1 MB of LOD4

sum(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$MinP_QTLs_3MB_Bin<0.0001, na.rm=TRUE)
#[1] 2730
2730/13633
#[1] 0.2002494
#And around 20% are +/- 1MB from LOD4

#Can we use this database to make manhattan plots with similar formatting (i.e., that can be aligned vertically better?)

#We would probably need a version that doesn't have missing values.

sum(is.na(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$MinP_QTLs_1MB_Bin))
#[1] 590
sum(is.na(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$HRLR_RNASeq_minP))
#[1] 0
sum(is.na(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$F2Behavior_RNASeq_minP))
#[1] 0

F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA<-F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins[is.na(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins$MinP_QTLs_1MB_Bin)==FALSE,]

colnames(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA)

pdf("ManhattanPlot_F2Behav_QTL_MinPval_toAlign.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="SYMBOLorENSEMBL", p="MinP_QTLs_1MB_Bin", suggestiveline=4, genomewideline=-log10(1.129854e-08))
dev.off()
#Note - I used a genomewide line determined by 0.05/length of full database of QTL pvalues (since this is actually what we started with... although really we should be using something even more stringent, since we-re looking at the minimum across several behavior QTLs (but they are highly correlated, so this would be awkward to calculate))

0.05/nrow(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA)
#[1] 3.833474e-06
#Same disclaimer here. 

pdf("ManhattanPlot_HRLR_RNASeq_MinPval_toAlign.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="SYMBOLorENSEMBL", p="HRLR_RNASeq_minP", suggestiveline=3, genomewideline=-log10(3.833474e-06), annotatePval = 0.00001, annotateTop=FALSE)
dev.off()

pdf("ManhattanPlot_F2Behavior_RNASeq_MinPval_toAlign.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="SYMBOLorENSEMBL", p="F2Behavior_RNASeq_minP", suggestiveline=3, genomewideline=-log10(3.833474e-06), annotatePval = 0.001, annotateTop=FALSE)
dev.off()

#It's hard to spot convergence - what if we looked for either the sum or the maximum across all three?
colnames(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA)

F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA$Sum_Pval_QTL_HRLR_F2Behav<-apply(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA[,c(61,63:64)], 1, sum)

pdf("ManhattanPlot_F2Behavior_SumPval_QTL_F0_F2RNASeq.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="SYMBOLorENSEMBL", p="Sum_Pval_QTL_HRLR_F2Behav", suggestiveline=2, genomewideline=-log10(3.833474e-06), annotatePval = 0.01, annotateTop=FALSE)
dev.off()

#hmmm... I don't know if that is legitimately informative. The pvalues for the qtl bins are almost always less than 0.05 - they are never going to drive the result.
#What if we plotted just the F0/F2 sum and then aligned it with the QTL results to visually show the overlap instead?

F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA$Sum_Pval_RNASeq_HRLR_F2Behav<-apply(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA[,c(63:64)], 1, sum)


pdf("ManhattanPlot_F2Behavior_SumPval_F0_F2RNASeq.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="SYMBOLorENSEMBL", p="Sum_Pval_RNASeq_HRLR_F2Behav", suggestiveline=2, genomewideline=-log10(3.833474e-06), annotatePval = 0.01, annotateTop=FALSE)
dev.off()


snpsOfInterest<-as.list(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA$SYMBOLorENSEMBL[F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA$MinP_QTLs_1MB_Bin<0.0001 & F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA$HRLR_RNASeq_minP<0.05 & F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA$F2Behavior_RNASeq_minP<0.05])

pdf("ManhattanPlot_F2Behavior_SumPval_F0_F2RNASeq_LOD4highlighted.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="SYMBOLorENSEMBL", p="Sum_Pval_RNASeq_HRLR_F2Behav", suggestiveline=2, genomewideline=-log10(3.833474e-06), highlight = snpsOfInterest)
dev.off()



#A version with +/-1MB
snpsOfInterest2<-as.list(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA$SYMBOLorENSEMBL[F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA$MinP_QTLs_3MB_Bin<0.0001 & F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA$HRLR_RNASeq_minP<0.05 & F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA$F2Behavior_RNASeq_minP<0.05])

pdf("ManhattanPlot_F2Behavior_SumPval_F0_F2RNASeq_LOD4_3MBBin_highlighted.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="SYMBOLorENSEMBL", p="Sum_Pval_RNASeq_HRLR_F2Behav", suggestiveline=2, genomewideline=-log10(3.833474e-06), highlight = snpsOfInterest2)
dev.off()


#Hmm... the highlighting might make it easier to look at the HR/LR and F2 Behavior version of the plots and visually see where they align.

pdf("ManhattanPlot_HRLR_RNASeq_MinPval_toAlign_LOD4highlighted.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="SYMBOLorENSEMBL", p="HRLR_RNASeq_minP", suggestiveline=3, genomewideline=-log10(3.833474e-06), highlight = snpsOfInterest)
dev.off()

pdf("ManhattanPlot_F2Behavior_RNASeq_MinPval_toAlign_LOD4highlighted.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="SYMBOLorENSEMBL", p="F2Behavior_RNASeq_minP", suggestiveline=3, genomewideline=-log10(3.833474e-06), highlight = snpsOfInterest)
dev.off()

# snpsOfInterest<-as.list(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA$ENSEMBL[F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA$MinP_QTLs_1MB_Bin<0.0001])

pdf("ManhattanPlot_F2Behav_QTL_MinPval_toAlign_LOD4highlighted.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="SYMBOLorENSEMBL", p="MinP_QTLs_1MB_Bin", suggestiveline=4, genomewideline=-log10(1.129854e-08), highlight = snpsOfInterest)
dev.off()

#Adding some name annotation so I can pull out top overlapping genes:


pdf("ManhattanPlot_F2Behavior_SumPval_F0_F2RNASeq_LOD4highlighted_Names.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="SYMBOLorENSEMBL", p="Sum_Pval_RNASeq_HRLR_F2Behav", suggestiveline=2, genomewideline=-log10(3.833474e-06), highlight = snpsOfInterest, annotatePval = 0.01, annotateTop=FALSE)
dev.off()

pdf("ManhattanPlot_F2Behavior_SumPval_F0_F2RNASeq_LOD4_3MBBin_highlighted_Names.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="SYMBOLorENSEMBL", p="Sum_Pval_RNASeq_HRLR_F2Behav", suggestiveline=2, genomewideline=-log10(3.833474e-06), highlight = snpsOfInterest2, annotatePval = 0.01, annotateTop=FALSE)
dev.off()

pdf("ManhattanPlot_HRLR_RNASeq_MinPval_toAlign_LOD4highlighted_Names.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="SYMBOLorENSEMBL", p="HRLR_RNASeq_minP", suggestiveline=3, genomewideline=-log10(3.833474e-06), highlight = snpsOfInterest, annotatePval = 0.00001, annotateTop=FALSE)
dev.off()

pdf("ManhattanPlot_F2Behavior_RNASeq_MinPval_toAlign_LOD4highlighted_Names.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="SYMBOLorENSEMBL", p="F2Behavior_RNASeq_minP", suggestiveline=3, genomewideline=-log10(3.833474e-06), highlight = snpsOfInterest, annotatePval = 0.001, annotateTop=FALSE)
dev.off()


#Making a new version with the bHR/bLR genes that survive FDR correction highlighted:


#snpsOfInterest<-F0_Meta_F2_DEGenes_FDR10_SYMBOLorENSEMBL

snpsOfInterest<-F0_Meta_F2_DEGenes_FDR10_ENSEMBL

pdf("ManhattanPlot_F2Behav_QTL_MinPval_toAlign_FDRhighlighted.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="ENSEMBL", p="MinP_QTLs_1MB_Bin", suggestiveline=4, genomewideline=-log10(1.129854e-08), highlight = snpsOfInterest)
dev.off()

pdf("ManhattanPlot_HRLR_RNASeq_MinPval_toAlign_FDRhighlighted.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="ENSEMBL", p="HRLR_RNASeq_minP", suggestiveline=3, genomewideline=-log10(3.833474e-06), highlight = snpsOfInterest, annotateTop=FALSE)
dev.off()

pdf("ManhattanPlot_F2Behavior_RNASeq_MinPval_toAlign_FDRhighlighted.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="ENSEMBL", p="F2Behavior_RNASeq_minP", suggestiveline=3, genomewideline=-log10(3.833474e-06), highlight = snpsOfInterest, annotateTop=FALSE)
dev.off()

pdf("ManhattanPlot_HRLR_RNASeq_MinPval_toAlign_FDRhighlighted_Names.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="ENSEMBL", p="HRLR_RNASeq_minP", suggestiveline=3, genomewideline=-log10(3.833474e-06), highlight = snpsOfInterest, annotatePval = 0.00001, annotateTop=FALSE)
dev.off()

pdf("ManhattanPlot_F2Behavior_RNASeq_MinPval_toAlign_FDRhighlighted_Names.pdf", width=15, height=5)
manhattan(F0_Meta_F2_4Models_vs_F2Behavior_QTLs_MBBins_NoNA, chr="Chr", bp="MB_Bin", snp="ENSEMBL", p="F2Behavior_RNASeq_minP", suggestiveline=3, genomewideline=-log10(3.833474e-06), highlight = snpsOfInterest, annotatePval = 0.001, annotateTop=FALSE)
dev.off()
