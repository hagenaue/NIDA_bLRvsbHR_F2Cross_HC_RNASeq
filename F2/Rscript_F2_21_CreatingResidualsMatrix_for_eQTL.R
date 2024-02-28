#This is the code we used to subtract out the effects of nuisance variables in the F2 data prior to eQTL analysis:
#Originally part of file Rscript_F2_TMM_20220104updated_MHannotation.R that was on Elaine's server space


#Normalizing the data to subtract out the effects of nuisance variables:

library(limma)

design<-model.matrix(~Percent.Intergenic+RibosomePerc+SequencingBatch+Sex+Dissector_AsNumeric+STGT_experience, data=F2_PCA_wMetaData)
design

# (Intercept) Percent.Intergenic RibosomePerc SequencingBatch2 SequencingBatch3 SexMale Dissector_AsNumeric STGT_experienceTRUE
# 63            1           0.156375       0.0096                0                0       1                   2                   1
# 181           1           0.157720       0.0101                0                0       1                   2                   1
# 185           1           0.171269       0.0086                0                0       1                   2                   1
# 183           1           0.158797       0.0082                0                0       1                   2                   1
# 83            1           0.160969       0.0088                0                0       1                   2                   1
# 5             1           0.150889       0.0072                0                0       0                   2                   1
# 189           1           0.153309       0.0084                0                0       0                   2                   1
# 236           1           0.157996       0.0112                0                0       0                   2                   1
# 46            1           0.156587       0.0075                0                0       0                   2                   1
# 230           1           0.153681       0.0072                0                0       0                   2                   1
# 88            1           0.160669       0.0071                0                0       1                   2                   1
# 92            1           0.145419       0.0072                0                0       1                   2                   1
# 114           1           0.146102       0.0085                0                0       1                   2                   1
# 26            1           0.154343       0.0081                0                0       0                   2                   1
# 93            1           0.160282       0.0107                0                0       0                   2                   1
# 131           1           0.143910       0.0104                0                0       0                   2                   1
# 43            1           0.123024       0.0148                0                0       1                   2                   1
# 53            1           0.117843       0.0176                0                0       1                   2                   1
# 106           1           0.153867       0.0111                0                0       1                   2                   1
# 155           1           0.134842       0.0148                0                0       0                   2                   1
# 161           1           0.163016       0.0092                0                0       0                   2                   1
# 205           1           0.138626       0.0122                0                0       0                   2                   1
# 36            1           0.120871       0.0149                0                0       0                   2                   1
# 51            1           0.137124       0.0098                0                0       1                   2                   1
# 90            1           0.120061       0.0138                0                0       1                   2                   1
# 163           1           0.135500       0.0151                0                0       0                   2                   1
# 37            1           0.166493       0.0073                0                0       0                   2                   1
# 44            1           0.162847       0.0068                0                0       1                   2                   1
# 52            1           0.139808       0.0085                0                0       1                   2                   1
# 56            1           0.141776       0.0088                0                0       1                   2                   1
# 61            1           0.155381       0.0066                0                0       1                   2                   1
# 72            1           0.150220       0.0079                0                0       0                   2                   1
# 73            1           0.128873       0.0110                0                0       1                   2                   1
# 74            1           0.146298       0.0088                0                0       1                   2                   1
# 108           1           0.126816       0.0136                0                0       1                   2                   1
# 122           1           0.137929       0.0090                0                0       0                   2                   1
# 145           1           0.136469       0.0067                0                0       0                   2                   1
# 192           1           0.124155       0.0124                0                0       0                   2                   1
# 204           1           0.150867       0.0065                0                0       0                   2                   1
# 206           1           0.111698       0.0212                0                0       0                   2                   1
# 42            1           0.150344       0.0104                0                0       1                   2                   1
# 62            1           0.153273       0.0063                0                0       1                   2                   1
# 65            1           0.114503       0.0165                0                0       1                   2                   1
# 76            1           0.132645       0.0113                0                0       1                   2                   1
# 82            1           0.131465       0.0138                0                0       1                   2                   1
# 84            1           0.121308       0.0185                0                0       0                   2                   0
# 143           1           0.137532       0.0111                0                0       0                   2                   0
# 144           1           0.113010       0.0180                0                0       0                   2                   0
# 182           1           0.123315       0.0112                0                0       1                   2                   1
# 227           1           0.122716       0.0150                0                0       0                   2                   0
# 71            1           0.151292       0.0074                0                0       0                   1                   1
# 170           1           0.144730       0.0072                0                0       0                   1                   1
# 130           1           0.142805       0.0064                0                0       1                   1                   1
# 148           1           0.137009       0.0090                0                0       1                   1                   1
# 190           1           0.151835       0.0072                0                0       0                   1                   1
# 156           1           0.147129       0.0067                0                0       0                   1                   1
# 138           1           0.146912       0.0061                0                0       1                   1                   1
# 141           1           0.140912       0.0065                0                0       1                   1                   1
# 132           1           0.141332       0.0079                0                0       0                   1                   1
# 225           1           0.149765       0.0064                0                0       0                   1                   1
# 127           1           0.133808       0.0100                0                0       1                   1                   1
# 147           1           0.139199       0.0078                0                0       1                   1                   1
# 79            1           0.121247       0.0128                0                0       0                   2                   1
# 150           1           0.119169       0.0188                0                0       1                   2                   1
# 159           1           0.123190       0.0119                0                0       1                   2                   1
# 164           1           0.137472       0.0105                0                0       0                   2                   1
# 194           1           0.113718       0.0193                0                0       1                   1                   1
# 229           1           0.118761       0.0157                0                0       0                   2                   1
# 6             1           0.128534       0.0124                0                0       0                   2                   1
# 23            1           0.115266       0.0154                0                0       0                   2                   1
# 59            1           0.126156       0.0128                0                0       0                   2                   1
# 157           1           0.125274       0.0131                0                0       1                   2                   1
# 165           1           0.127183       0.0118                0                0       0                   2                   1
# 175           1           0.114503       0.0117                0                0       1                   2                   1
# 184           1           0.105300       0.0142                0                0       1                   2                   1
# 195           1           0.128233       0.0093                0                0       1                   2                   1
# 198           1           0.128915       0.0093                0                0       1                   2                   1
# 247           1           0.110550       0.0131                0                0       0                   2                   0
# 153           1           0.149044       0.0093                0                0       1                   2                   1
# 146           1           0.139454       0.0101                0                0       1                   1                   1
# 228           1           0.146336       0.0066                0                0       0                   1                   1
# 60            1           0.146411       0.0098                0                0       0                   1                   1
# 137           1           0.126908       0.0089                0                0       1                   1                   1
# 193           1           0.139840       0.0072                0                0       1                   1                   1
# 111           1           0.137615       0.0143                0                0       0                   2                   1
# 70            1           0.141413       0.0106                0                0       0                   1                   1
# 197           1           0.125328       0.0120                0                0       1                   1                   1
# 80            1           0.157354       0.0100                0                0       0                   1                   1
# 112           1           0.111215       0.0261                0                0       0                   2                   1
# 167           1           0.119616       0.0174                0                0       1                   2                   1
# 176           1           0.117566       0.0147                0                0       1                   2                   1
# 186           1           0.125427       0.0136                0                0       1                   2                   1
# 212           1           0.125719       0.0141                0                0       0                   2                   0
# 248           1           0.109089       0.0220                0                0       0                   2                   0
# 196           1           0.124757       0.0143                0                0       1                   1                   1
# 226           1           0.141558       0.0104                1                0       0                   1                   0
# 237           1           0.127319       0.0106                1                0       0                   1                   0
# 33            1           0.142121       0.0082                1                0       0                   2                   1
# 34            1           0.122231       0.0151                1                0       0                   2                   1
# 97            1           0.114593       0.0113                1                0       1                   2                   1
# 99            1           0.120654       0.0096                1                0       1                   2                   1
# 100           1           0.128887       0.0090                1                0       1                   2                   1
# 101           1           0.110921       0.0147                1                0       1                   2                   1
# 107           1           0.122754       0.0117                1                0       1                   2                   1
# 113           1           0.137847       0.0089                1                0       1                   2                   1
# 115           1           0.141659       0.0072                1                0       1                   2                   1
# 116           1           0.123551       0.0133                1                0       1                   2                   1
# 125           1           0.105398       0.0177                1                0       1                   2                   1
# 126           1           0.108537       0.0146                1                0       1                   2                   1
# 133           1           0.123358       0.0102                1                0       0                   2                   1
# 151           1           0.121652       0.0121                1                0       1                   2                   1
# 168           1           0.116026       0.0122                1                0       1                   2                   1
# 169           1           0.114466       0.0161                1                0       0                   2                   1
# 171           1           0.110549       0.0141                1                0       0                   2                   1
# 179           1           0.125089       0.0095                1                0       0                   2                   1
# 191           1           0.121172       0.0123                1                0       0                   2                   1
# 213           1           0.111479       0.0144                1                0       0                   2                   1
# 128           1           0.142683       0.0075                1                0       1                   1                   1
# 129           1           0.144721       0.0072                1                0       1                   1                   1
# 24            1           0.156471       0.0077                1                0       0                   1                   1
# 35            1           0.167130       0.0066                1                0       0                   2                   1
# 98            1           0.153066       0.0054                1                0       1                   1                   1
# 118           1           0.154093       0.0068                1                0       1                   1                   1
# 180           1           0.161632       0.0060                1                0       0                   1                   1
# 211           1           0.150445       0.0049                1                0       0                   1                   0
# [ reached getOption("max.print") -- omitted 120 rows ]
# attr(,"assign")
# [1] 0 1 2 3 3 4 5 6
# attr(,"contrasts")
# attr(,"contrasts")$SequencingBatch
# [1] "contr.treatment"
# 
# attr(,"contrasts")$Sex
# [1] "contr.treatment"
# 
# attr(,"contrasts")$STGT_experience
# [1] "contr.treatment"

v<-voom(F2_RawHitCount_dge_noLowHits_TMM, design, plot=TRUE)
v
# 
# An object of class "EList"

vfit<-lmFit(v,design)

str(vfit)
#This is the linear regression output (raw)
# 
# Formal class 'MArrayLM' [package "limma"] with 1 slot

efit<-eBayes(vfit)

ResidualsMatrix <- residuals(vfit, v)

ResidualsMatrixEbayes <- residuals(efit, v)

plot(ResidualsMatrix[1,]~ResidualsMatrixEbayes[1,])
#These are the same because the Ebayes correction is applied to the t-statistics (I think) and not the coefficients (Log2FC) for each of the predictor variables

hist(ResidualsMatrix[1,])
#The residuals are centered around 0. 
#I'm going to add back in the average gene expression for each gene

str(v)

# Formal class 'EList' [package "limma"] with 1 slot
# ..@ .Data:List of 4
# .. ..$ :'data.frame':	245 obs. of  3 variables:
#   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
# .. .. ..$ lib.size    : num [1:245] 21467043 18947762 20168506 27714579 21487895 ...
# .. .. ..$ norm.factors: num [1:245] 0.984 0.987 0.996 1.042 1.013 ...
# .. ..$ : num [1:14056, 1:245] 8.24 7.54 2.49 5.38 2.86 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:14056] "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# .. .. .. ..$ : chr [1:245] "SL469967" "SL469968" "SL469969" "SL469970" ...
# .. ..$ : num [1:14056, 1:245] 50.7 50.8 15.5 45.4 22.9 ...
# .. ..$ : num [1:245, 1:8] 1 1 1 1 1 1 1 1 1 1 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:245] "63" "181" "185" "183" ...
# .. .. .. ..$ : chr [1:8] "(Intercept)" "Percent.Intergenic" "RibosomePerc" "SequencingBatch2" ...
# .. .. ..- attr(*, "assign")= int [1:8] 0 1 2 3 3 4 5 6
# .. .. ..- attr(*, "contrasts")=List of 3
# .. .. .. ..$ SequencingBatch: chr "contr.treatment"
# .. .. .. ..$ Sex            : chr "contr.treatment"
# .. .. .. ..$ STGT_experience: chr "contr.treatment"
# ..$ names: chr [1:4] "targets" "E" "weights" "design"

#Ok, I'm going to cheat - hypothetically the values that I add back in don't matter much, as long as all values afterwards are >0




MeanLog2Expression_perGene<-apply(F2_RawHitCount_dge_noLowHits_TMM, 1, mean)

sum(is.na(MeanLog2Expression_perGene))
#[1] 0

plot(vfit$Amean~MeanLog2Expression_perGene)
#monotonic but not linear

plot(vfit$Amean~log2(MeanLog2Expression_perGene))
#Yep, vfit$Amean is basically the average gene expression for each gene

#Double-checking whether the intercept performs similarly
str(vfit)

str(vfit$coefficients)

plot(vfit$coefficients[,1]~MeanLog2Expression_perGene)
#neither monotonic or linear

ResidualsMatrix_wMean<-ResidualsMatrix+vfit$Amean

vfit$Amean[1]

# ENSRNOG00000014303 
# 7.75896 

hist(ResidualsMatrix_wMean[1,])
#Looks good - centered around 7.75

#one more to be sure:

vfit$Amean[2]
# ENSRNOG00000014330 
# 8.338438 

hist(ResidualsMatrix_wMean[2,])
#looks good - now centered on 8.something

str(ResidualsMatrix_wMean)