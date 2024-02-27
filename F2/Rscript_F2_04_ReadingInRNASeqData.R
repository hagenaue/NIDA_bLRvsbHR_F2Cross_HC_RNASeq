#F2 HC RNA-Seq Dataset
#04_Reading in the RNA-Seq Data
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.


#Reading in the RNA-Seq data:

#First: What RNA-Seq data do we have?

#Files:

#Elaine_es103rn_FCgene_raw
##This one appears to be the raw counts without any normalization or filtering

#Elaine_es103rn_FCgene_4R
##This one appears to be the raw counts without any normalization or filtering, with additional annotation columns removed so that it reads into R as a data matrix


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

library(limma)

list.files()

F2_RawHitCount<-read.delim("Elaine_es103rn_FCgene_4R", sep="\t", header=TRUE, row.names=1)

str(F2_RawHitCount)
# 'data.frame':	32883 obs. of  250 variables:
#   $ SL469967: int  0 0 0 0 0 0 2 0 0 0 ...
# $ SL469968: int  0 0 0 0 0 0 4 0 0 0 ...
# $ SL469969: int  0 0 0 0 0 0 2 0 0 0 ...
# $ SL469970: int  0 0 0 0 0 0 5 0 0 0 ...
# $ SL469971: int  0 0 0 0 0 0 0 0 0 0 ...

#Rows=genes
#Columns=subjects
#values=# of reads per gene (transcript) per subject

#Double-check that the subjects are in the same order as your metadata

str(F2_PCA_wMetaData)

#Let's check visually to see if things are in the same order:
cbind(colnames(F2_RawHitCount), F2_PCA_wMetaData$row_name)


F2_PCA_wMetaData_temp<-F2_PCA_wMetaData[order(F2_PCA_wMetaData$row_name),]

cbind(colnames(F2_RawHitCount), F2_PCA_wMetaData_temp$row_name)

F2_PCA_wMetaData<-F2_PCA_wMetaData_temp

#Looks good.
# [,1]       [,2]      
# [1,] "SL469967" "sl469967"
# [2,] "SL469968" "sl469968"
# [3,] "SL469969" "sl469969"
# [4,] "SL469970" "sl469970"
# [5,] "SL469971" "sl469971"
# [6,] "SL469972" "sl469972"
# [7,] "SL469973" "sl469973"
# [8,] "SL469974" "sl469974"
# [9,] "SL469975" "sl469975"
# [10,] "SL469976" "sl469976"
# [11,] "SL469977" "sl469977"
# [12,] "SL469978" "sl469978"
# [13,] "SL469979" "sl469979"
# [14,] "SL469980" "sl469980"
# [15,] "SL469981" "sl469981"
# [16,] "SL469982" "sl469982"
# [17,] "SL469983" "sl469983"
# [18,] "SL469984" "sl469984"
# [19,] "SL469985" "sl469985"
# [20,] "SL469986" "sl469986"
# [21,] "SL469987" "sl469987"
# [22,] "SL469988" "sl469988"
# [23,] "SL469989" "sl469989"
# [24,] "SL469990" "sl469990"
# [25,] "SL469991" "sl469991"
# [26,] "SL469992" "sl469992"
# [27,] "SL469993" "sl469993"
# [28,] "SL469994" "sl469994"
# [29,] "SL469995" "sl469995"
# [30,] "SL469996" "sl469996"
# [31,] "SL469997" "sl469997"
# [32,] "SL469998" "sl469998"
# [33,] "SL469999" "sl469999"
# [34,] "SL470000" "sl470000"
# [35,] "SL470001" "sl470001"
# [36,] "SL470002" "sl470002"
# [37,] "SL470003" "sl470003"
# [38,] "SL470004" "sl470004"
# [39,] "SL470005" "sl470005"
# [40,] "SL470006" "sl470006"
# [41,] "SL470007" "sl470007"
# [42,] "SL470008" "sl470008"
# [43,] "SL470009" "sl470009"
# [44,] "SL470010" "sl470010"
# [45,] "SL470011" "sl470011"
# [46,] "SL470012" "sl470012"
# [47,] "SL470013" "sl470013"
# [48,] "SL470014" "sl470014"
# [49,] "SL470015" "sl470015"
# [50,] "SL470016" "sl470016"
# [51,] "SL470017" "sl470017"
# [52,] "SL470018" "sl470018"
# [53,] "SL470019" "sl470019"
# [54,] "SL470020" "sl470020"
# [55,] "SL470021" "sl470021"
# [56,] "SL470022" "sl470022"
# [57,] "SL470023" "sl470023"
# [58,] "SL470024" "sl470024"
# [59,] "SL470025" "sl470025"
# [60,] "SL470026" "sl470026"
# [61,] "SL470027" "sl470027"
# [62,] "SL470028" "sl470028"
# [63,] "SL470029" "sl470029"
# [64,] "SL470030" "sl470030"
# [65,] "SL470031" "sl470031"
# [66,] "SL470032" "sl470032"
# [67,] "SL470033" "sl470033"
# [68,] "SL470034" "sl470034"
# [69,] "SL470035" "sl470035"
# [70,] "SL470036" "sl470036"
# [71,] "SL470037" "sl470037"
# [72,] "SL470038" "sl470038"
# [73,] "SL470039" "sl470039"
# [74,] "SL470040" "sl470040"
# [75,] "SL470041" "sl470041"
# [76,] "SL470042" "sl470042"
# [77,] "SL470043" "sl470043"
# [78,] "SL470044" "sl470044"
# [79,] "SL470045" "sl470045"
# [80,] "SL470046" "sl470046"
# [81,] "SL470047" "sl470047"
# [82,] "SL470048" "sl470048"
# [83,] "SL470049" "sl470049"
# [84,] "SL470050" "sl470050"
# [85,] "SL470051" "sl470051"
# [86,] "SL470052" "sl470052"
# [87,] "SL470053" "sl470053"
# [88,] "SL470054" "sl470054"
# [89,] "SL470055" "sl470055"
# [90,] "SL470056" "sl470056"
# [91,] "SL470057" "sl470057"
# [92,] "SL470058" "sl470058"
# [93,] "SL470059" "sl470059"
# [94,] "SL470060" "sl470060"
# [95,] "SL470061" "sl470061"
# [96,] "SL470508" "sl470508"
# [97,] "SL470509" "sl470509"
# [98,] "SL470510" "sl470510"
# [99,] "SL470511" "sl470511"
# [100,] "SL470512" "sl470512"
# [101,] "SL470513" "sl470513"
# [102,] "SL470514" "sl470514"
# [103,] "SL470515" "sl470515"
# [104,] "SL470516" "sl470516"
# [105,] "SL470517" "sl470517"
# [106,] "SL470518" "sl470518"
# [107,] "SL470519" "sl470519"
# [108,] "SL470520" "sl470520"
# [109,] "SL470521" "sl470521"
# [110,] "SL470522" "sl470522"
# [111,] "SL470523" "sl470523"
# [112,] "SL470524" "sl470524"
# [113,] "SL470525" "sl470525"
# [114,] "SL470526" "sl470526"
# [115,] "SL470527" "sl470527"
# [116,] "SL470528" "sl470528"
# [117,] "SL470529" "sl470529"
# [118,] "SL470530" "sl470530"
# [119,] "SL470531" "sl470531"
# [120,] "SL470532" "sl470532"
# [121,] "SL470533" "sl470533"
# [122,] "SL470534" "sl470534"
# [123,] "SL470535" "sl470535"
# [124,] "SL470536" "sl470536"
# [125,] "SL470537" "sl470537"
# [126,] "SL470538" "sl470538"
# [127,] "SL470539" "sl470539"
# [128,] "SL470540" "sl470540"
# [129,] "SL470541" "sl470541"
# [130,] "SL470542" "sl470542"
# [131,] "SL470543" "sl470543"
# [132,] "SL470544" "sl470544"
# [133,] "SL470545" "sl470545"
# [134,] "SL470546" "sl470546"
# [135,] "SL470547" "sl470547"
# [136,] "SL470548" "sl470548"
# [137,] "SL470549" "sl470549"
# [138,] "SL470550" "sl470550"
# [139,] "SL470551" "sl470551"
# [140,] "SL470552" "sl470552"
# [141,] "SL470553" "sl470553"
# [142,] "SL470554" "sl470554"
# [143,] "SL470555" "sl470555"
# [144,] "SL470556" "sl470556"
# [145,] "SL470557" "sl470557"
# [146,] "SL470558" "sl470558"
# [147,] "SL470559" "sl470559"
# [148,] "SL470560" "sl470560"
# [149,] "SL470561" "sl470561"
# [150,] "SL470562" "sl470562"
# [151,] "SL470563" "sl470563"
# [152,] "SL470564" "sl470564"
# [153,] "SL470565" "sl470565"
# [154,] "SL470566" "sl470566"
# [155,] "SL470567" "sl470567"
# [156,] "SL470568" "sl470568"
# [157,] "SL470569" "sl470569"
# [158,] "SL470570" "sl470570"
# [159,] "SL470571" "sl470571"
# [160,] "SL470572" "sl470572"
# [161,] "SL470573" "sl470573"
# [162,] "SL470574" "sl470574"
# [163,] "SL470575" "sl470575"
# [164,] "SL470576" "sl470576"
# [165,] "SL470577" "sl470577"
# [166,] "SL470578" "sl470578"
# [167,] "SL470579" "sl470579"
# [168,] "SL470580" "sl470580"
# [169,] "SL470581" "sl470581"
# [170,] "SL470582" "sl470582"
# [171,] "SL470583" "sl470583"
# [172,] "SL470584" "sl470584"
# [173,] "SL470585" "sl470585"
# [174,] "SL470586" "sl470586"
# [175,] "SL470587" "sl470587"
# [176,] "SL470588" "sl470588"
# [177,] "SL470589" "sl470589"
# [178,] "SL470590" "sl470590"
# [179,] "SL470591" "sl470591"
# [180,] "SL470592" "sl470592"
# [181,] "SL470593" "sl470593"
# [182,] "SL470594" "sl470594"
# [183,] "SL470595" "sl470595"
# [184,] "SL470596" "sl470596"
# [185,] "SL470597" "sl470597"
# [186,] "SL470598" "sl470598"
# [187,] "SL470599" "sl470599"
# [188,] "SL470600" "sl470600"
# [189,] "SL470601" "sl470601"
# [190,] "SL470602" "sl470602"
# [191,] "SL470847" "sl470847"
# [192,] "SL470848" "sl470848"
# [193,] "SL470849" "sl470849"
# [194,] "SL470850" "sl470850"
# [195,] "SL470851" "sl470851"
# [196,] "SL470852" "sl470852"
# [197,] "SL470853" "sl470853"
# [198,] "SL470854" "sl470854"
# [199,] "SL470855" "sl470855"
# [200,] "SL470856" "sl470856"
# [201,] "SL470857" "sl470857"
# [202,] "SL470858" "sl470858"
# [203,] "SL470859" "sl470859"
# [204,] "SL470860" "sl470860"
# [205,] "SL470861" "sl470861"
# [206,] "SL470862" "sl470862"
# [207,] "SL470863" "sl470863"
# [208,] "SL470864" "sl470864"
# [209,] "SL470865" "sl470865"
# [210,] "SL470866" "sl470866"
# [211,] "SL470867" "sl470867"
# [212,] "SL470868" "sl470868"
# [213,] "SL470869" "sl470869"
# [214,] "SL470870" "sl470870"
# [215,] "SL470871" "sl470871"
# [216,] "SL470872" "sl470872"
# [217,] "SL470873" "sl470873"
# [218,] "SL470874" "sl470874"
# [219,] "SL470875" "sl470875"
# [220,] "SL470876" "sl470876"
# [221,] "SL470877" "sl470877"
# [222,] "SL470878" "sl470878"
# [223,] "SL470879" "sl470879"
# [224,] "SL470880" "sl470880"
# [225,] "SL470881" "sl470881"
# [226,] "SL470882" "sl470882"
# [227,] "SL470883" "sl470883"
# [228,] "SL470884" "sl470884"
# [229,] "SL470885" "sl470885"
# [230,] "SL470886" "sl470886"
# [231,] "SL470887" "sl470887"
# [232,] "SL470888" "sl470888"
# [233,] "SL470889" "sl470889"
# [234,] "SL470890" "sl470890"
# [235,] "SL470891" "sl470891"
# [236,] "SL470892" "sl470892"
# [237,] "SL470893" "sl470893"
# [238,] "SL470894" "sl470894"
# [239,] "SL470895" "sl470895"
# [240,] "SL470896" "sl470896"
# [241,] "SL470897" "sl470897"
# [242,] "SL470898" "sl470898"
# [243,] "SL470899" "sl470899"
# [244,] "SL470900" "sl470900"
# [245,] "SL470901" "sl470901"
# [246,] "SL470902" "sl470902"
# [247,] "SL470903" "sl470903"
# [248,] "SL470904" "sl470904"
# [249,] "SL470905" "sl470905"
# [250,] "SL470906" "sl470906"


################

#Working with the expression data:

F2_RawHitCount_ExpressionData<-as.matrix(F2_RawHitCount)
str(F2_RawHitCount_ExpressionData)
# int [1:32883, 1:250] 0 0 0 0 0 0 2 0 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:32883] "ENSRNOG00000046319" "ENSRNOG00000047964" "ENSRNOG00000050370" "ENSRNOG00000032365" ...
# ..$ : chr [1:250] "SL469967" "SL469968" "SL469969" "SL469970" ...


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

library(edgeR)


F2_LibrarySize<-apply(F2_RawHitCount_ExpressionData, 2, sum)

pdf("Histogram_F2_LibrarySize.pdf", height=5, width=5)
hist(F2_LibrarySize)
dev.off()
#There is one *huge* outlier mixed in here - we should probably toss that sample.

summary(F2_LibrarySize)
# Min.      1st Qu.    Median      Mean   3rd Qu.      Max. 
#17916830  20037385  21660660  22408230  23428387 109689712 

#Median Library Size = close to what was intended!
#21,660,660 
#one huge outlier.

F2_PCA_wMetaData$LibrarySize<-F2_LibrarySize

###################2/1/2022 Elaine

F2_LibrarySize<-apply(F2_RawHitCount_ExpressionData, 2, sum)

pdf("Histogram_F2_LibrarySize_NoOutlier.pdf", height=5, width=5)
hist(F2_LibrarySize)
dev.off()
#Huge outlier removed here

summary(F2_LibrarySize)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 17916830 20004950 21613379 21986830 23363734 31494272 

#Median Library Size = close to what was intended!
#21,613,379 


