#F2 HC RNA-Seq Dataset
#12_Screening for a few more potential confounds or sources of noise
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.



#Digging into a few more variables that might be important sources of noise:


# GT/ST is a huge multi-day procedure - probably important to control for:

levels(F2_PCA_wMetaData$LearningClassification_AsFactor)
#[1] "ST" ""   "GT" "IN"

sum(F2_PCA_wMetaData$LearningClassification=="")

#Does ST/GT experience overlap with other batch variables?

table(F2_PCA_wMetaData$LearningClassification=="", F2_PCA_wMetaData$HPC_Dissection_Date_AsFactor_Releveled)

# 12/3/2020 12/1/2020 12/10/2020 12/11/2020 12/14/2020 12/2/2020 12/4/2020 12/7/2020 12/8/2020 12/9/2020
# FALSE        14        11         34         18          1        16        24        24        33        30
# TRUE          0         0          0         12         18         0         4         5         1         0

#Looks like the animals that didn't get ST/GT training are almost entirely in the 12/11 and 12/14 batches which look really different for all the PCs, esp. PC3
#12/7 also has a handful, so does PC4

table(F2_PCA_wMetaData$LearningClassification=="", F2_PCA_wMetaData$DOD)

# 12/10/2013 12/11/2013 12/13/2013 12/17/2013 12/18/2013 12/4/2013 12/5/2013 12/6/2013 12/9/2013
# FALSE         28         55          0         47         52         0         0         6        17
# TRUE           0          0          9         11          4         3        10         3         0

#Also partially overlaps with DOD.

table(F2_PCA_wMetaData$LearningClassification=="", F2_PCA_wMetaData$SequencingBatch)
#       1  2  3
# FALSE 88 91 26
# TRUE   7  3 30

levels(F2_PCA_wMetaData$STGT_experience)
#[1] "FALSE" "TRUE"

table(F2_PCA_wMetaData$STGT_experience)

# FALSE  TRUE 
# 40   205 
#False= no ST/GT training/experience


#######################

#Code for looking for batch effects in *behavioral variables* that might add to noise:

colnames(F2_PCA_wMetaData)

summary.lm(lm(Total_LocoScore~Sex_AsFactor+Loco_TestDate+Loco_TestDay+Loco_Batch+Loco_Rack+Loco_Box, data=F2_PCA_wMetaData))

# Call:
#   lm(formula = Total_LocoScore ~ Sex_AsFactor + Loco_TestDate + 
#        Loco_TestDay + Loco_Batch + Loco_Rack + Loco_Box, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -518.85 -208.27  -46.84  191.42  895.12 
# 
# Coefficients: (1 not defined because of singularities)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             490.535     96.556   5.080 7.64e-07 ***
#   Sex_AsFactorFemale      170.104     92.798   1.833  0.06805 .  
# Loco_TestDate9/25/2013   45.834     76.128   0.602  0.54771    
# Loco_TestDate9/26/2013 -105.341     96.467  -1.092  0.27595    
# Loco_TestDate9/27/2013  -53.780    111.879  -0.481  0.63118    
# Loco_TestDay                 NA         NA      NA       NA    
# Loco_Batch               -7.678     27.553  -0.279  0.78074    
# Loco_Rack               -70.628     37.456  -1.886  0.06057 .  
# Loco_Box                 11.799      3.799   3.106  0.00213 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 292.3 on 237 degrees of freedom
# Multiple R-squared:  0.08335,	Adjusted R-squared:  0.05628 
# F-statistic: 3.079 on 7 and 237 DF,  p-value: 0.00401
# 

summary.lm(lm(EPM_Percent_Time_Open_Arm~Sex_AsFactor+EPM_TestDate+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData))

# Call:
#   lm(formula = EPM_Percent_Time_Open_Arm ~ Sex_AsFactor + EPM_TestDate + 
#        EPM_TestDay + EPM_TestOrderByDay, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -18.329  -7.389  -2.167   6.305  29.082 
# 
# Coefficients: (1 not defined because of singularities)
# Estimate Std. Error t value Pr(>|t|)   
# (Intercept)             7.78975    2.45984   3.167  0.00175 **
#   Sex_AsFactorFemale      6.57987    2.57960   2.551  0.01139 * 
#   EPM_TestDate10/21/2013  5.88213    2.83890   2.072  0.03937 * 
#   EPM_TestDate10/22/2013  8.41663    2.99392   2.811  0.00536 **
#   EPM_TestDate10/23/2013  3.54202    3.06339   1.156  0.24877   
# EPM_TestDate10/24/2013  1.67592    3.05886   0.548  0.58429   
# EPM_TestDate10/25/2013  0.71853    2.87903   0.250  0.80314   
# EPM_TestDate10/30/2013  1.96173    3.15166   0.622  0.53426   
# EPM_TestDate10/31/2013  3.09331    2.95830   1.046  0.29681   
# EPM_TestDate11/1/2013   3.54075    2.95381   1.199  0.23186   
# EPM_TestDate11/2/2013   3.09649    3.13475   0.988  0.32428   
# EPM_TestDay                  NA         NA      NA       NA   
# EPM_TestOrderByDay     -0.04457    0.13395  -0.333  0.73965   
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 10.01 on 233 degrees of freedom
# Multiple R-squared:  0.1281,	Adjusted R-squared:  0.08689 
# F-statistic: 3.111 on 11 and 233 DF,  p-value: 0.0006259
# 

summary.lm(lm(EPM_Boli~Sex_AsFactor+EPM_TestDate+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData))

# Call:
#   lm(formula = EPM_Boli ~ Sex_AsFactor + EPM_TestDate + EPM_TestDay + 
#        EPM_TestOrderByDay, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.8067 -0.9605 -0.2686  0.1375  6.6991 
# 
# Coefficients: (1 not defined because of singularities)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             1.503903   0.372997   4.032 7.53e-05 ***
#   Sex_AsFactorFemale     -0.459747   0.398772  -1.153   0.2501    
# EPM_TestDate10/21/2013  0.117404   0.430454   0.273   0.7853    
# EPM_TestDate10/22/2013 -0.509856   0.453714  -1.124   0.2623    
# EPM_TestDate10/23/2013  0.112175   0.469915   0.239   0.8115    
# EPM_TestDate10/24/2013 -0.198133   0.469271  -0.422   0.6733    
# EPM_TestDate10/25/2013  0.041799   0.439863   0.095   0.9244    
# EPM_TestDate10/30/2013  0.685131   0.477817   1.434   0.1530    
# EPM_TestDate10/31/2013  0.005516   0.448377   0.012   0.9902    
# EPM_TestDate11/1/2013  -0.043642   0.447681  -0.097   0.9224    
# EPM_TestDate11/2/2013   0.314363   0.475075   0.662   0.5088    
# EPM_TestDay                   NA         NA      NA       NA    
# EPM_TestOrderByDay     -0.034761   0.020568  -1.690   0.0924 .  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 1.518 on 230 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.1434,	Adjusted R-squared:  0.1024 
# F-statistic: 3.501 on 11 and 230 DF,  p-value: 0.0001516

summary.lm(lm(EPM_DistanceTraveled~Sex_AsFactor+EPM_TestDate+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData))


# 
# Call:
#   lm(formula = EPM_DistanceTraveled ~ Sex_AsFactor + EPM_TestDate + 
#        EPM_TestDay + EPM_TestOrderByDay, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -927.5 -302.1  -55.4  262.5 3223.5 
# 
# Coefficients: (1 not defined because of singularities)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            1877.919    113.936  16.482  < 2e-16 ***
#   Sex_AsFactorFemale      346.220    119.483   2.898  0.00412 ** 
#   EPM_TestDate10/21/2013  185.914    131.493   1.414  0.15873    
# EPM_TestDate10/22/2013  265.319    138.674   1.913  0.05694 .  
# EPM_TestDate10/23/2013   56.740    141.891   0.400  0.68961    
# EPM_TestDate10/24/2013  -62.313    141.681  -0.440  0.66048    
# EPM_TestDate10/25/2013  -15.612    133.352  -0.117  0.90690    
# EPM_TestDate10/30/2013  125.345    145.980   0.859  0.39142    
# EPM_TestDate10/31/2013  127.921    137.024   0.934  0.35149    
# EPM_TestDate11/1/2013   241.096    136.816   1.762  0.07935 .  
# EPM_TestDate11/2/2013   251.240    145.196   1.730  0.08489 .  
# EPM_TestDay                  NA         NA      NA       NA    
# EPM_TestOrderByDay       -5.609      6.204  -0.904  0.36689    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 463.9 on 233 degrees of freedom
# Multiple R-squared:  0.1222,	Adjusted R-squared:  0.08072 
# F-statistic: 2.948 on 11 and 233 DF,  p-value: 0.001128

summary.lm(lm(EPM_Time_Immobile~Sex_AsFactor+EPM_TestDate+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData))


# Call:
#   lm(formula = EPM_Time_Immobile ~ Sex_AsFactor + EPM_TestDate + 
#        EPM_TestDay + EPM_TestOrderByDay, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -42.547 -15.602  -1.333  12.516  63.453 
# 
# Coefficients: (1 not defined because of singularities)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            119.7644     5.4177  22.106  < 2e-16 ***
# Sex_AsFactorFemale     -21.7544     5.6815  -3.829 0.000165 ***
# EPM_TestDate10/21/2013  -8.2974     6.2526  -1.327 0.185796    
# EPM_TestDate10/22/2013 -19.1869     6.5940  -2.910 0.003968 ** 
# EPM_TestDate10/23/2013 -12.7273     6.7470  -1.886 0.060492 .  
# EPM_TestDate10/24/2013  -2.3873     6.7371  -0.354 0.723390    
# EPM_TestDate10/25/2013  -5.3972     6.3410  -0.851 0.395555    
# EPM_TestDate10/30/2013   3.6969     6.9414   0.533 0.594832    
# EPM_TestDate10/31/2013   0.2716     6.5156   0.042 0.966790    
# EPM_TestDate11/1/2013   -4.3215     6.5057  -0.664 0.507177    
# EPM_TestDate11/2/2013   -1.3253     6.9042  -0.192 0.847947    
# EPM_TestDay                  NA         NA      NA       NA    
# EPM_TestOrderByDay      -0.2500     0.2950  -0.847 0.397651    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 22.06 on 233 degrees of freedom
# Multiple R-squared:  0.3079,	Adjusted R-squared:  0.2752 
# F-statistic: 9.422 on 11 and 233 DF,  p-value: 5.523e-14  

####################################Batch effects in behavioral data?

step(object=lm(Total_LocoScore~Sex_AsFactor+Loco_TestDate+Loco_TestDay+Loco_Batch+Loco_Rack+Loco_Box, data=F2_PCA_wMetaData), scope=c(lower=lm(Total_LocoScore~Sex,data=F2_PCA_wMetaData),upper=lm(Total_LocoScore~Sex_AsFactor+Loco_TestDate+Loco_TestDay+Loco_Batch+Loco_Rack+Loco_Box, data=F2_PCA_wMetaData)), scale = 0, direction = "backward", trace = 1, keep = NULL, steps = 1000, k = 2)

# Start:  AIC=2789.91
# Total_LocoScore ~ Sex_AsFactor + Loco_TestDate + Loco_TestDay + 
#   Loco_Batch + Loco_Rack + Loco_Box
# 
# 
# Step:  AIC=2789.91
# Total_LocoScore ~ Sex_AsFactor + Loco_TestDate + Loco_Batch + 
#   Loco_Rack + Loco_Box
# 
# Df Sum of Sq      RSS    AIC
# - Loco_Batch     1      6633 20250013 2788.0
# - Loco_TestDate  3    382676 20626056 2788.5
# <none>                       20243380 2789.9
# - Sex_AsFactor   1    287000 20530380 2791.4
# - Loco_Rack      1    303700 20547080 2791.6
# - Loco_Box       1    823930 21067310 2797.7
# 
# Step:  AIC=2787.99
# Total_LocoScore ~ Sex_AsFactor + Loco_TestDate + Loco_Rack + 
#   Loco_Box
# 
# Df Sum of Sq      RSS    AIC
# - Loco_TestDate  3    376395 20626409 2786.5
# <none>                       20250013 2788.0
# - Loco_Rack      1    309192 20559205 2789.7
# - Sex_AsFactor   1    310323 20560336 2789.7
# - Loco_Box       1    838519 21088533 2795.9
# 
# Step:  AIC=2786.5
# Total_LocoScore ~ Sex_AsFactor + Loco_Rack + Loco_Box
# 
# Df Sum of Sq      RSS    AIC
# <none>                      20626409 2786.5
# - Loco_Rack     1    326579 20952988 2788.3
# - Sex_AsFactor  1    327748 20954157 2788.4
# - Loco_Box      1    839068 21465477 2794.3
# 
# Call:
#   lm(formula = Total_LocoScore ~ Sex_AsFactor + Loco_Rack + Loco_Box, 
#      data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)  Sex_AsFactorFemale           Loco_Rack            Loco_Box  
# 498.80               73.19              -73.11               11.78  
# 

step(object=lm(Total_LocoScore~Sex_AsFactor+Loco_TestDay+Loco_Batch+Loco_Rack+Loco_Box, data=F2_PCA_wMetaData), scope=c(lower=lm(Total_LocoScore~Sex,data=F2_PCA_wMetaData),upper=lm(Total_LocoScore~Sex_AsFactor+Loco_TestDate+Loco_TestDay+Loco_Batch+Loco_Rack+Loco_Box, data=F2_PCA_wMetaData)), scale = 0, direction = "backward", trace = 1, keep = NULL, steps = 1000, k = 2)

# Start:  AIC=2790.49
# Total_LocoScore ~ Sex_AsFactor + Loco_TestDay + Loco_Batch + 
#   Loco_Rack + Loco_Box
# 
# Df Sum of Sq      RSS    AIC
# - Loco_TestDay  1       262 20626056 2788.5
# - Loco_Batch    1       450 20626244 2788.5
# - Sex_AsFactor  1     81674 20707468 2789.5
# <none>                      20625794 2790.5
# - Loco_Rack     1    326821 20952615 2792.3
# - Loco_Box      1    829622 21455416 2798.2
# 
# Step:  AIC=2788.5
# Total_LocoScore ~ Sex_AsFactor + Loco_Batch + Loco_Rack + Loco_Box
# 
# Df Sum of Sq      RSS    AIC
# - Loco_Batch    1       353 20626409 2786.5
# <none>                      20626056 2788.5
# - Sex_AsFactor  1    306368 20932424 2790.1
# - Loco_Rack     1    326912 20952968 2790.3
# - Loco_Box      1    835196 21461252 2796.2
# 
# Step:  AIC=2786.5
# Total_LocoScore ~ Sex_AsFactor + Loco_Rack + Loco_Box
# 
# Df Sum of Sq      RSS    AIC
# <none>                      20626409 2786.5
# - Loco_Rack     1    326579 20952988 2788.3
# - Sex_AsFactor  1    327748 20954157 2788.4
# - Loco_Box      1    839068 21465477 2794.3
# 
# Call:
#   lm(formula = Total_LocoScore ~ Sex_AsFactor + Loco_Rack + Loco_Box, 
#      data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)  Sex_AsFactorFemale           Loco_Rack            Loco_Box  
# 498.80               73.19              -73.11               11.78  
# 

summary.lm(lm(formula = Total_LocoScore ~ Sex_AsFactor + Loco_Rack + Loco_Box, data = F2_PCA_wMetaData))

# Call:
#   lm(formula = Total_LocoScore ~ Sex_AsFactor + Loco_Rack + Loco_Box, 
#      data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -496.63 -232.60  -46.41  176.80  914.64 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)         498.805     68.970   7.232 6.28e-12 ***
#   Sex_AsFactorFemale   73.193     37.403   1.957  0.05151 .  
# Loco_Rack           -73.115     37.429  -1.953  0.05193 .  
# Loco_Box             11.783      3.763   3.131  0.00196 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 292.6 on 241 degrees of freedom
# Multiple R-squared:  0.06601,	Adjusted R-squared:  0.05438 
# F-statistic: 5.678 on 3 and 241 DF,  p-value: 0.0008988




step(object=lm(EPM_Percent_Time_Open_Arm~Sex_AsFactor+EPM_TestDate+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData), scope=c(lower=lm(EPM_Percent_Time_Open_Arm~Sex, data=F2_PCA_wMetaData),upper=lm(EPM_Percent_Time_Open_Arm~Sex_AsFactor+EPM_TestDate+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData)), scale = 0, direction = "backward", trace = 1, keep = NULL, steps = 1000, k = 2)
# Start:  AIC=1140.69
# EPM_Percent_Time_Open_Arm ~ Sex_AsFactor + EPM_TestDate + EPM_TestDay + 
#   EPM_TestOrderByDay
# 
# 
# Step:  AIC=1140.69
# EPM_Percent_Time_Open_Arm ~ Sex_AsFactor + EPM_TestDate + EPM_TestOrderByDay
# 
# Df Sum of Sq   RSS    AIC
# - EPM_TestDate        9   1376.75 24746 1136.7
# - EPM_TestOrderByDay  1     11.10 23380 1138.8
# <none>                            23369 1140.7
# - Sex_AsFactor        1    652.55 24021 1145.4
# 
# Step:  AIC=1136.71
# EPM_Percent_Time_Open_Arm ~ Sex_AsFactor + EPM_TestOrderByDay
# 
# Df Sum of Sq   RSS    AIC
# - EPM_TestOrderByDay  1      9.31 24755 1134.8
# <none>                            24746 1136.7
# - Sex_AsFactor        1    694.21 25440 1141.5
# 
# Step:  AIC=1134.8
# EPM_Percent_Time_Open_Arm ~ Sex_AsFactor
# 
# Df Sum of Sq   RSS    AIC
# <none>                      24755 1134.8
# - Sex_AsFactor  1    2045.8 26801 1152.3
# 
# Call:
#   lm(formula = EPM_Percent_Time_Open_Arm ~ Sex_AsFactor, data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)  Sex_AsFactorFemale  
# 10.746               5.779  
# 

step(object=lm(EPM_Percent_Time_Open_Arm~Sex_AsFactor+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData), scope=c(lower=lm(EPM_Percent_Time_Open_Arm~Sex, data=F2_PCA_wMetaData),upper=lm(EPM_Percent_Time_Open_Arm~Sex_AsFactor+EPM_TestDate+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData)), scale = 0, direction = "backward", trace = 1, keep = NULL, steps = 1000, k = 2)


# Start:  AIC=1137.99
# EPM_Percent_Time_Open_Arm ~ Sex_AsFactor + EPM_TestDay + EPM_TestOrderByDay
# 
# Df Sum of Sq   RSS    AIC
# - EPM_TestOrderByDay  1     10.37 24683 1136.1
# - EPM_TestDay         1     73.07 24746 1136.7
# <none>                            24673 1138.0
# - Sex_AsFactor        1    705.03 25378 1142.9
# 
# Step:  AIC=1136.09
# EPM_Percent_Time_Open_Arm ~ Sex_AsFactor + EPM_TestDay
# 
# Df Sum of Sq   RSS    AIC
# - EPM_TestDay   1     72.01 24755 1134.8
# <none>                      24683 1136.1
# - Sex_AsFactor  1   2057.07 26740 1153.7
# 
# Step:  AIC=1134.8
# EPM_Percent_Time_Open_Arm ~ Sex_AsFactor
# 
# Df Sum of Sq   RSS    AIC
# <none>                      24755 1134.8
# - Sex_AsFactor  1    2045.8 26801 1152.3
# 
# Call:
#   lm(formula = EPM_Percent_Time_Open_Arm ~ Sex_AsFactor, data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)  Sex_AsFactorFemale  
# 10.746               5.779  




step(object=lm(EPM_Boli~Sex_AsFactor+EPM_TestDate+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData), scope=c(lower=lm(EPM_Boli~Sex, data=F2_PCA_wMetaData),upper=lm(EPM_Boli~Sex_AsFactor+EPM_TestDate+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData)), scale = 0, direction = "backward", trace = 1, keep = NULL, steps = 1000, k = 2)


#Step:  AIC=213.6
# EPM_Boli ~ Sex_AsFactor + EPM_TestDate + EPM_TestOrderByDay
# 
# Df Sum of Sq    RSS    AIC
# - EPM_TestDate        9   18.4688 548.23 203.90
# - Sex_AsFactor        1    3.0615 532.82 213.00
# <none>                            529.76 213.60
# - EPM_TestOrderByDay  1    6.5787 536.34 214.59
# 
# Step:  AIC=203.9
# EPM_Boli ~ Sex_AsFactor + EPM_TestOrderByDay
# 
# Df Sum of Sq    RSS    AIC
# - EPM_TestOrderByDay  1    3.5750 551.80 203.47
# <none>                            548.23 203.90
# - Sex_AsFactor        1    6.9682 555.19 204.95
# 
# Step:  AIC=203.47
# EPM_Boli ~ Sex_AsFactor
# 
# Df Sum of Sq    RSS    AIC
# <none>                      551.80 203.47
# - Sex_AsFactor  1    66.649 618.45 229.06
# 
# Call:
#   lm(formula = EPM_Boli ~ Sex_AsFactor, data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)  Sex_AsFactorFemale  
# 1.248              -1.050  



step(object=lm(EPM_DistanceTraveled~Sex_AsFactor+EPM_TestDate+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData), scope=c(lower=lm(EPM_DistanceTraveled~Sex, data=F2_PCA_wMetaData),upper=lm(EPM_DistanceTraveled~Sex_AsFactor+EPM_TestDate+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData)), scale = 0, direction = "backward", trace = 1, keep = NULL, steps = 1000, k = 2)

# Start:  AIC=3020.1
# EPM_DistanceTraveled ~ Sex_AsFactor + EPM_TestDate + EPM_TestDay + 
#   EPM_TestOrderByDay
# 
# 
# Step:  AIC=3020.1
# EPM_DistanceTraveled ~ Sex_AsFactor + EPM_TestDate + EPM_TestOrderByDay
# 
# Df Sum of Sq      RSS    AIC
# - EPM_TestDate        9   3035141 53170533 3016.5
# - EPM_TestOrderByDay  1    175875 50311266 3019.0
# <none>                            50135392 3020.1
# - Sex_AsFactor        1   1806686 51942078 3026.8
# 
# Step:  AIC=3016.5
# EPM_DistanceTraveled ~ Sex_AsFactor + EPM_TestOrderByDay
# 
# Df Sum of Sq      RSS    AIC
# - EPM_TestOrderByDay  1    177198 53347731 3015.3
# <none>                            53170533 3016.5
# - Sex_AsFactor        1   1896307 55066840 3023.1
# 
# Step:  AIC=3015.32
# EPM_DistanceTraveled ~ Sex_AsFactor
# 
# Df Sum of Sq      RSS    AIC
# <none>                      53347731 3015.3
# - Sex_AsFactor  1   3764498 57112229 3030.0
# 
# Call:
#   lm(formula = EPM_DistanceTraveled ~ Sex_AsFactor, data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)  Sex_AsFactorFemale  
# 1950.9               247.9  


step(object=lm(EPM_DistanceTraveled~Sex_AsFactor+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData), scope=c(lower=lm(EPM_DistanceTraveled~Sex, data=F2_PCA_wMetaData),upper=lm(EPM_DistanceTraveled~Sex_AsFactor+EPM_TestDate+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData)), scale = 0, direction = "backward", trace = 1, keep = NULL, steps = 1000, k = 2)

# Start:  AIC=3017.53
# EPM_DistanceTraveled ~ Sex_AsFactor + EPM_TestDay + EPM_TestOrderByDay
# 
# Df Sum of Sq      RSS    AIC
# - EPM_TestOrderByDay  1    169530 53128851 3016.3
# - EPM_TestDay         1    211211 53170533 3016.5
# <none>                            52959321 3017.5
# - Sex_AsFactor        1   1863777 54823099 3024.0
# 
# Step:  AIC=3016.31
# EPM_DistanceTraveled ~ Sex_AsFactor + EPM_TestDay
# 
# Df Sum of Sq      RSS    AIC
# - EPM_TestDay   1    218880 53347731 3015.3
# <none>                      53128851 3016.3
# - Sex_AsFactor  1   3736073 56864924 3031.0
# 
# Step:  AIC=3015.32
# EPM_DistanceTraveled ~ Sex_AsFactor
# 
# Df Sum of Sq      RSS    AIC
# <none>                      53347731 3015.3
# - Sex_AsFactor  1   3764498 57112229 3030.0
# 
# Call:
#   lm(formula = EPM_DistanceTraveled ~ Sex_AsFactor, data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)  Sex_AsFactorFemale  
# 1950.9               247.9  
# 




step(object=lm(EPM_Time_Immobile~Sex_AsFactor+EPM_TestDate+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData), scope=c(lower=lm(EPM_Time_Immobile~Sex, data=F2_PCA_wMetaData),upper=lm(EPM_Time_Immobile~Sex_AsFactor+EPM_TestDate+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData)), scale = 0, direction = "backward", trace = 1, keep = NULL, steps = 1000, k = 2)

# Start:  AIC=1527.58
# EPM_Time_Immobile ~ Sex_AsFactor + EPM_TestDate + EPM_TestDay + 
#   EPM_TestOrderByDay
# 
# 
# Step:  AIC=1527.58
# EPM_Time_Immobile ~ Sex_AsFactor + EPM_TestDate + EPM_TestOrderByDay
# 
# Df Sum of Sq    RSS    AIC
# - EPM_TestOrderByDay  1     349.4 113709 1526.3
# <none>                            113360 1527.6
# - EPM_TestDate        9    9699.6 123060 1529.7
# - Sex_AsFactor        1    7133.0 120493 1540.5
# 
# Step:  AIC=1526.33
# EPM_Time_Immobile ~ Sex_AsFactor + EPM_TestDate
# 
# Df Sum of Sq    RSS    AIC
# <none>                      113709 1526.3
# - EPM_TestDate  9      9488 123198 1528.0
# - Sex_AsFactor  1     41025 154734 1599.8
# 
# Call:
#   lm(formula = EPM_Time_Immobile ~ Sex_AsFactor + EPM_TestDate, 
#      data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)      Sex_AsFactorFemale  EPM_TestDate10/21/2013  EPM_TestDate10/22/2013  EPM_TestDate10/23/2013  EPM_TestDate10/24/2013  EPM_TestDate10/25/2013  EPM_TestDate10/30/2013  EPM_TestDate10/31/2013  
# 117.9919                -25.9319                 -9.3751                -19.0234                -13.1386                 -2.4076                 -5.8320                  2.5506                 -0.3671  
# EPM_TestDate11/1/2013   EPM_TestDate11/2/2013  
# -4.8832                 -0.9764  

summary.lm(lm(formula = EPM_Time_Immobile ~ Sex_AsFactor + EPM_TestDate, data = F2_PCA_wMetaData))

# Call:
#   lm(formula = EPM_Time_Immobile ~ Sex_AsFactor + EPM_TestDate, 
#      data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -42.590 -15.731  -1.705  12.169  61.910 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            117.9919     4.9946  23.624  < 2e-16 ***
#   Sex_AsFactorFemale     -25.9319     2.8223  -9.188  < 2e-16 ***
#   EPM_TestDate10/21/2013  -9.3751     6.1182  -1.532  0.12679    
# EPM_TestDate10/22/2013 -19.0234     6.5872  -2.888  0.00424 ** 
#   EPM_TestDate10/23/2013 -13.1386     6.7255  -1.954  0.05195 .  
# EPM_TestDate10/24/2013  -2.4076     6.7330  -0.358  0.72098    
# EPM_TestDate10/25/2013  -5.8320     6.3164  -0.923  0.35680    
# EPM_TestDate10/30/2013   2.5506     6.8042   0.375  0.70811    
# EPM_TestDate10/31/2013  -0.3671     6.4679  -0.057  0.95479    
# EPM_TestDate11/1/2013   -4.8832     6.4679  -0.755  0.45102    
# EPM_TestDate11/2/2013   -0.9764     6.8878  -0.142  0.88739    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 22.04 on 234 degrees of freedom
# Multiple R-squared:  0.3057,	Adjusted R-squared:  0.2761 
# F-statistic: 10.31 on 10 and 234 DF,  p-value: 2.303e-14


step(object=lm(EPM_Time_Immobile~Sex_AsFactor+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData), scope=c(lower=lm(EPM_Time_Immobile~Sex, data=F2_PCA_wMetaData),upper=lm(EPM_Time_Immobile~Sex_AsFactor+EPM_TestDate+EPM_TestDay+EPM_TestOrderByDay, data=F2_PCA_wMetaData)), scale = 0, direction = "backward", trace = 1, keep = NULL, steps = 1000, k = 2)

# Start:  AIC=1527.64
# EPM_Time_Immobile ~ Sex_AsFactor + EPM_TestDay + EPM_TestOrderByDay
# 
# Df Sum of Sq    RSS    AIC
# - EPM_TestOrderByDay  1     117.8 121159 1525.9
# <none>                            121041 1527.6
# - EPM_TestDay         1    2018.4 123060 1529.7
# - Sex_AsFactor        1    9386.1 130427 1543.9
# 
# Step:  AIC=1525.88
# EPM_Time_Immobile ~ Sex_AsFactor + EPM_TestDay
# 
# Df Sum of Sq    RSS    AIC
# <none>                      121159 1525.9
# - EPM_TestDay   1      2039 123198 1528.0
# - Sex_AsFactor  1     40857 162016 1595.1
# 
# Call:
#   lm(formula = EPM_Time_Immobile ~ Sex_AsFactor + EPM_TestDay, 
#      data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)  Sex_AsFactorFemale         EPM_TestDay  
# 106.802             -25.830               1.014  


summary.lm( lm(formula = EPM_Time_Immobile ~ Sex_AsFactor + EPM_TestDay, data = F2_PCA_wMetaData))

# Call:
#   lm(formula = EPM_Time_Immobile ~ Sex_AsFactor + EPM_TestDay, 
#      data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -44.673 -14.911  -1.619  13.746  64.169 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        106.8023     3.3664  31.726   <2e-16 ***
#   Sex_AsFactorFemale -25.8305     2.8594  -9.034   <2e-16 ***
#   EPM_TestDay          1.0136     0.5023   2.018   0.0447 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 22.38 on 242 degrees of freedom
# Multiple R-squared:  0.2603,	Adjusted R-squared:  0.2541 
# F-statistic: 42.57 on 2 and 242 DF,  p-value: < 2.2e-16

#################################################

#Making a version of Locoscore cleaned of batch-related effects:

Total_LocoScore_Cleaned<-residuals(lm(Total_LocoScore~Loco_Rack+Loco_Box, data=F2_PCA_wMetaData))

str(Total_LocoScore_Cleaned)

# Named num [1:245] -232 393 -110 102 140 ...
# - attr(*, "names")= chr [1:245] "63" "181" "185" "183" ...

str(unname(Total_LocoScore_Cleaned))
#num [1:245] -232 393 -110 102 140 ...


hist(unname(Total_LocoScore_Cleaned))

F2_PCA_wMetaData$Total_LocoScore_Cleaned<-unname(Total_LocoScore_Cleaned)

EPM_Time_Immobile_Cleaned1<-residuals(lm(EPM_Time_Immobile~EPM_TestDate, data=F2_PCA_wMetaData))

hist(unname(EPM_Time_Immobile_Cleaned1))

F2_PCA_wMetaData$EPM_Time_Immobile_Cleaned1<-unname(EPM_Time_Immobile_Cleaned1)

EPM_Time_Immobile_Cleaned2<-residuals(lm(EPM_Time_Immobile~EPM_TestDay, data=F2_PCA_wMetaData))

hist(unname(EPM_Time_Immobile_Cleaned2))

F2_PCA_wMetaData$EPM_Time_Immobile_Cleaned2<-unname(EPM_Time_Immobile_Cleaned2)


#################

# Looking for overlap between other batch-related variables:

str(F2_PCA_wMetaData)
# 'data.frame':	245 obs. of  65 variables:

table(F2_PCA_wMetaData$EPM_TestDate)
# 10/18/2013 10/21/2013 10/22/2013 10/23/2013 10/24/2013 10/25/2013 10/30/2013 10/31/2013  11/1/2013  11/2/2013 
# 21         34         24         22         22         29         21         26         26         20 


table(F2_PCA_wMetaData$EPM_TestDate,F2_PCA_wMetaData$HPC_Dissection_Date)
# 12/1/2020 12/10/2020 12/11/2020 12/14/2020 12/2/2020 12/3/2020 12/4/2020 12/7/2020 12/8/2020 12/9/2020
# 10/18/2013         0          6          3          3         0         0         0         1         4         4
# 10/21/2013         3          4          0          0         3         3         8         1         4         8
# 10/22/2013         2          1          3          2         2         2         4         6         0         2
# 10/23/2013         0          4          7          2         2         1         0         2         4         0
# 10/24/2013         0          5          5          1         1         0         2         4         4         0
# 10/25/2013         2          2          2          2         0         2         3         6        10         0
# 10/30/2013         2          2          6          8         0         0         1         2         0         0
# 10/31/2013         1          6          1          1         2         0         5         3         5         2
# 11/1/2013          1          1          2          0         3         3         4         2         3         7
# 11/2/2013          0          3          1          0         3         3         1         2         0         7

table(F2_PCA_wMetaData$DOD)
# 12/10/2013 12/11/2013 12/13/2013 12/17/2013 12/18/2013  12/4/2013  12/5/2013  12/6/2013  12/9/2013 
# 28         55          9         58         56          3         10          9         17

table(F2_PCA_wMetaData$DOD,F2_PCA_wMetaData$SequencingBatch)

#             1  2  3
# 12/10/2013 19  9  0
# 12/11/2013 30 25  0
# 12/13/2013  0  0  9
# 12/17/2013 16 20 22
# 12/18/2013 30 23  3
# 12/4/2013   0  0  3
# 12/5/2013   0  0 10
# 12/6/2013   0  2  7
# 12/9/2013   0 15  2

table(F2_PCA_wMetaData$STGT_experience,F2_PCA_wMetaData$SequencingBatch)

#        1  2  3
# FALSE  7  3 30
# TRUE  88 91 26

table(F2_PCA_wMetaData$Age.days.,F2_PCA_wMetaData$Sex)

#       Female Male
# 113      0    4
# 114      0    3
# 116      0    3
# 117      0    4
# 118      0    8
# 119      0   18
# 121      0   26
# 122      0   22
# 123      0   34
# 124      2    0
# 125     31    0
# 127      2    0
# 128      4    0
# 129     45    0
# 130     31    0
# 131      6    0
# 132      2    0

##############################

#At a glance, it seems like the change in RNA metrics over the different dissection dates may differ by dissector
#This makes intuitive sense, and suggests that it may be better to not treat dissector and dissection date as independent variables

HPC_Dissector_DissectionDate<-paste(F2_PCA_wMetaData$Dissector, F2_PCA_wMetaData$HPC_Dissection_Date, sep="_")

table(HPC_Dissector_DissectionDate)

# EHB_12/10/2020 EHB_12/11/2020 EHB_12/14/2020  EHB_12/4/2020  EHB_12/7/2020  EHB_12/8/2020  EHB_12/9/2020   KA_12/1/2020  KA_12/10/2020  KA_12/11/2020  KA_12/14/2020   KA_12/2/2020   KA_12/3/2020   KA_12/4/2020   KA_12/7/2020   KA_12/8/2020 
# 15              8              9             13             11             13             11             11             19             22             10             16             14             15             18             21 
# KA_12/9/2020 
# 19 



anova(lm(F2_PCA_wMetaData$mRNAperc~F2_PCA_wMetaData$Dissector * F2_PCA_wMetaData$HPC_Dissection_Date))
# 
# Analysis of Variance Table
# 
# Response: F2_PCA_wMetaData$mRNAperc
# Df    Sum Sq   Mean Sq F value    Pr(>F)    
# F2_PCA_wMetaData$Dissector                                        1 0.0115048 0.0115048 85.3985 < 2.2e-16 ***
#   F2_PCA_wMetaData$HPC_Dissection_Date                              9 0.0165383 0.0018376 13.6402 < 2.2e-16 ***
#   F2_PCA_wMetaData$Dissector:F2_PCA_wMetaData$HPC_Dissection_Date   6 0.0035263 0.0005877  4.3625 0.0003424 ***
#   Residuals                                                       228 0.0307160 0.0001347                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



#sample from Behavioral Script:    Anova(lm(EPM_Time_Immobile~Sex_AsFactor*GenPheno, data=Behavior[Behavior$RNASeq=="Elaine",], contrasts=list(GenPheno=contr.sum, Sex_AsFactor=contr.sum)), type=3)

library(car)


Anova(lm(mRNAperc~Dissector + HPC_Dissection_Date, data=F2_PCA_wMetaData, contrasts=list(Dissector=contr.sum, HPC_Dissection_Date=contr.sum)), type=3)
# Error in Anova.III.lm(mod, error, singular.ok = singular.ok, ...) : 
#   there are aliased coefficients in the model

# Anova Table (Type III tests)
# 
# Response: mRNAperc
#                       Sum Sq  Df    F value    Pr(>F)    
# (Intercept)            114.466   1 782225.277 < 2.2e-16 ***
#   Dissector             0.018   1    121.023 < 2.2e-16 ***
#   HPC_Dissection_Date   0.017   9     12.557 3.169e-16 ***
#   Residuals             0.034 234                         
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


Anova(lm(Percent.Intergenic~Dissector + HPC_Dissection_Date, data=F2_PCA_wMetaData, contrasts=list(Dissector=contr.sum, HPC_Dissection_Date=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: Percent.Intergenic
# Sum Sq  Df   F value    Pr(>F)    
# (Intercept)         3.4207   1 24621.778 < 2.2e-16 ***
#   Dissector           0.0244   1   175.653 < 2.2e-16 ***
#   HPC_Dissection_Date 0.0214   9    17.096 < 2.2e-16 ***
#   Residuals           0.0325 234                        
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


Anova(lm(RibosomePerc~Dissector + HPC_Dissection_Date, data=F2_PCA_wMetaData, contrasts=list(Dissector=contr.sum, HPC_Dissection_Date=contr.sum)), type=3)

# Anova Table (Type III tests)
# 
# Response: RibosomePerc
# Sum Sq  Df  F value    Pr(>F)    
# (Intercept)         0.0179062   1 1748.189 < 2.2e-16 ***
#   Dissector           0.0013233   1  129.193 < 2.2e-16 ***
#   HPC_Dissection_Date 0.0012802   9   13.887 < 2.2e-16 ***
#   Residuals           0.0023968 234                       
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


###################################

anova(lm(F2_PCA_wMetaData$mRNAperc~F2_PCA_wMetaData$Dissector * F2_PCA_wMetaData$HPC_Dissection_Day))

# Analysis of Variance Table
# 
# Response: F2_PCA_wMetaData$mRNAperc
# Df   Sum Sq   Mean Sq F value    Pr(>F)    
# F2_PCA_wMetaData$Dissector                                       1 0.011505 0.0115048  69.026 7.084e-15 ***
#   F2_PCA_wMetaData$HPC_Dissection_Day                              1 0.008348 0.0083477  50.084 1.597e-11 ***
#   F2_PCA_wMetaData$Dissector:F2_PCA_wMetaData$HPC_Dissection_Day   1 0.002264 0.0022644  13.586 0.0002815 ***
#   Residuals                                                      241 0.040168 0.0001667                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


anova(lm(F2_PCA_wMetaData$RNAconc~F2_PCA_wMetaData$Dissector * F2_PCA_wMetaData$HPC_Dissection_Date))
# Analysis of Variance Table
# 
# Response: F2_PCA_wMetaData$RNAconc
# Df Sum Sq Mean Sq  F value  Pr(>F)    
# F2_PCA_wMetaData$Dissector                                        1 630801  630801 234.3099 < 2e-16 ***
#   F2_PCA_wMetaData$HPC_Dissection_Date                              9 603623   67069  24.9128 < 2e-16 ***
#   F2_PCA_wMetaData$Dissector:F2_PCA_wMetaData$HPC_Dissection_Date   6  43296    7216   2.6804 0.01561 *  
#   Residuals                                                       228 613813    2692                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


anova(lm(F2_PCA_wMetaData$RNAconc~F2_PCA_wMetaData$Dissector * F2_PCA_wMetaData$HPC_Dissection_Day))
#Analysis of Variance Table
# 
# Response: F2_PCA_wMetaData$RNAconc
# Df  Sum Sq Mean Sq  F value Pr(>F)    
# F2_PCA_wMetaData$Dissector                                       1  630801  630801 121.6814 <2e-16 ***
#   F2_PCA_wMetaData$HPC_Dissection_Day                              1    8722    8722   1.6824 0.1958    
# F2_PCA_wMetaData$Dissector:F2_PCA_wMetaData$HPC_Dissection_Day   1    2659    2659   0.5129 0.4746    
# Residuals                                                      241 1249352    5184                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


anova(lm(F2_PCA_wMetaData$RibosomePerc~F2_PCA_wMetaData$Dissector * F2_PCA_wMetaData$HPC_Dissection_Date))
# Analysis of Variance Table
# 
# Response: F2_PCA_wMetaData$RibosomePerc
# Df     Sum Sq    Mean Sq  F value  Pr(>F)    
# F2_PCA_wMetaData$Dissector                                        1 0.00109339 0.00109339 109.8126 < 2e-16 ***
#   F2_PCA_wMetaData$HPC_Dissection_Date                              9 0.00128017 0.00014224  14.2857 < 2e-16 ***
#   F2_PCA_wMetaData$Dissector:F2_PCA_wMetaData$HPC_Dissection_Date   6 0.00012662 0.00002110   2.1195 0.05202 .  
# Residuals                                                       228 0.00227017 0.00000996                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

anova(lm(F2_PCA_wMetaData$RibosomePerc~F2_PCA_wMetaData$Dissector * F2_PCA_wMetaData$HPC_Dissection_Day))
# Analysis of Variance Table
# 
# Response: F2_PCA_wMetaData$RibosomePerc
# Df     Sum Sq    Mean Sq F value    Pr(>F)    
# F2_PCA_wMetaData$Dissector                                       1 0.00109339 0.00109339 84.5548 < 2.2e-16 ***
#   F2_PCA_wMetaData$HPC_Dissection_Day                              1 0.00050002 0.00050002 38.6677 2.199e-09 ***
#   F2_PCA_wMetaData$Dissector:F2_PCA_wMetaData$HPC_Dissection_Day   1 0.00006054 0.00006054  4.6814   0.03147 *  
#   Residuals                                                      241 0.00311641 0.00001293                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



pdf("F2_RibosomePercVsDissector_DissectionDate.pdf",height=5, width=10)
boxplot(F2_PCA_wMetaData$RibosomePerc~HPC_Dissector_DissectionDate)
dev.off()


pdf("F2_RibosomePercVsDissector_DissectionDay.pdf",height=5, width=10)
plot(F2_PCA_wMetaData$RibosomePerc~F2_PCA_wMetaData$HPC_Dissection_Day,col=(F2_PCA_wMetaData$Dissector_AsNumeric+1))
dev.off()


pdf("F2_mRNApercVsDissector_DissectionDate.pdf",height=5, width=10)
boxplot(F2_PCA_wMetaData$mRNAperc~HPC_Dissector_DissectionDate)
dev.off()


pdf("F2_mRNApercVsDissector_DissectionDay.pdf",height=5, width=10)
plot(F2_PCA_wMetaData$mRNAperc~F2_PCA_wMetaData$HPC_Dissection_Day,col=(F2_PCA_wMetaData$Dissector_AsNumeric+1))
dev.off()



F2_PCA_wMetaData$HPC_Dissector_DissectionDate<-HPC_Dissector_DissectionDate


tapply(F2_PCA_wMetaData$RNAconc, F2_PCA_wMetaData$HPC_Dissector_DissectionDate, mean)

tapply(F2_PCA_wMetaData$mRNAperc, F2_PCA_wMetaData$HPC_Dissector_DissectionDate, mean)

tapply(F2_PCA_wMetaData$RibosomePerc, F2_PCA_wMetaData$HPC_Dissector_DissectionDate, mean)

tapply(F2_PCA_wMetaData$RIN, F2_PCA_wMetaData$HPC_Dissector_DissectionDate, mean)

tapply(F2_PCA_wMetaData$DV200, F2_PCA_wMetaData$HPC_Dissector_DissectionDate, mean)

tapply(F2_PCA_wMetaData$PC1, F2_PCA_wMetaData$HPC_Dissector_DissectionDate, mean)

tapply(F2_PCA_wMetaData$PC2, F2_PCA_wMetaData$HPC_Dissector_DissectionDate, mean)

tapply(F2_PCA_wMetaData$PC3, F2_PCA_wMetaData$HPC_Dissector_DissectionDate, mean)

tapply(F2_PCA_wMetaData$PC4, F2_PCA_wMetaData$HPC_Dissector_DissectionDate, mean)

tapply(F2_PCA_wMetaData$PC5, F2_PCA_wMetaData$HPC_Dissector_DissectionDate, mean)

tapply(F2_PCA_wMetaData$PC6, F2_PCA_wMetaData$HPC_Dissector_DissectionDate, mean)


##################################
#Creating new collapsed variable: Dissector with DissectionDate: see Dissector_DissectionDate_PCs_trimmed.xls

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed<-as.character(F2_PCA_wMetaData$HPC_Dissector_DissectionDate)

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="EHB_12/4/2020"]<-"Batch1"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="EHB_12/7/2020"]<-"Batch2"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="EHB_12/8/2020"]<-"Batch3"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="EHB_12/9/2020"]<-"Batch3"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="EHB_12/10/2020"]<-"Batch3"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="EHB_12/11/2020"]<-"Batch4"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="EHB_12/14/2020"]<-"Batch5"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="KA_12/1/2020"]<-"Batch6"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="KA_12/2/2020"]<-"Batch7"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="KA_12/3/2020"]<-"Batch8"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="KA_12/4/2020"]<-"Batch9"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="KA_12/7/2020"]<-"Batch9"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="KA_12/8/2020"]<-"Batch9"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="KA_12/9/2020"]<-"Batch10"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="KA_12/10/2020"]<-"Batch10"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="KA_12/11/2020"]<-"Batch11"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed[F2_PCA_wMetaData$HPC_Dissector_DissectionDate=="KA_12/14/2020"]<-"Batch11"

F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed<-as.factor(F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed)
F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed<-relevel(F2_PCA_wMetaData$HPC_Dissector_DissectionDate_Collapsed, ref="Batch3")


summary.lm(lm(PC2 ~ HPC_Dissector_DissectionDate_Collapsed, data=F2_PCA_wMetaData))

# 
# Call:
#   lm(formula = PC2 ~ HPC_Dissector_DissectionDate_Collapsed, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -26.141  -8.463  -0.425   8.514  31.917 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                     6.4229     1.9356   3.318  0.00105 ** 
#   HPC_Dissector_DissectionDate_CollapsedBatch1  -16.9468     3.8711  -4.378 1.81e-05 ***
#   HPC_Dissector_DissectionDate_CollapsedBatch10  -4.0420     2.7553  -1.467  0.14371    
# HPC_Dissector_DissectionDate_CollapsedBatch11  -3.7471     2.8831  -1.300  0.19499    
# HPC_Dissector_DissectionDate_CollapsedBatch2    0.3756     4.1266   0.091  0.92755    
# HPC_Dissector_DissectionDate_CollapsedBatch4  -11.3182     4.6915  -2.412  0.01661 *  
#   HPC_Dissector_DissectionDate_CollapsedBatch5   -7.7093     4.4700  -1.725  0.08591 .  
# HPC_Dissector_DissectionDate_CollapsedBatch6  -23.0149     4.1266  -5.577 6.73e-08 ***
#   HPC_Dissector_DissectionDate_CollapsedBatch7  -17.6523     3.5886  -4.919 1.64e-06 ***
#   HPC_Dissector_DissectionDate_CollapsedBatch8    1.2950     3.7660   0.344  0.73125    
# HPC_Dissector_DissectionDate_CollapsedBatch9   -7.5283     2.5401  -2.964  0.00335 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 12.09 on 234 degrees of freedom
# Multiple R-squared:  0.229,	Adjusted R-squared:  0.196 
# F-statistic:  6.95 on 10 and 234 DF,  p-value: 1.606e-09
# 


summary.lm(lm(PC4 ~ HPC_Dissector_DissectionDate_Collapsed, data=F2_PCA_wMetaData))


# Call:
#   lm(formula = PC4 ~ HPC_Dissector_DissectionDate_Collapsed, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -38.943  -4.300   0.227   4.627  22.478 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                     2.22490    1.24842   1.782  0.07602 .  
# HPC_Dissector_DissectionDate_CollapsedBatch1   -6.98009    2.49684  -2.796  0.00561 ** 
#   HPC_Dissector_DissectionDate_CollapsedBatch10   4.26499    1.77711   2.400  0.01718 *  
#   HPC_Dissector_DissectionDate_CollapsedBatch11 -13.67942    1.85958  -7.356 3.17e-12 ***
#   HPC_Dissector_DissectionDate_CollapsedBatch2   -0.05099    2.66164  -0.019  0.98473    
# HPC_Dissector_DissectionDate_CollapsedBatch4    0.04309    3.02597   0.014  0.98865    
# HPC_Dissector_DissectionDate_CollapsedBatch5   -1.07000    2.88310  -0.371  0.71088    
# HPC_Dissector_DissectionDate_CollapsedBatch6   -7.94728    2.66164  -2.986  0.00313 ** 
#   HPC_Dissector_DissectionDate_CollapsedBatch7   -3.79218    2.31463  -1.638  0.10269    
# HPC_Dissector_DissectionDate_CollapsedBatch8   -3.15678    2.42904  -1.300  0.19502    
# HPC_Dissector_DissectionDate_CollapsedBatch9    0.43423    1.63835   0.265  0.79121    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 7.796 on 234 degrees of freedom
# Multiple R-squared:  0.3356,	Adjusted R-squared:  0.3072 
# F-statistic: 11.82 on 10 and 234 DF,  p-value: < 2.2e-16


summary.lm(lm(PC4 ~ HPC_Dissector_DissectionDate, data=F2_PCA_wMetaData))


# Call:
#   lm(formula = PC4 ~ HPC_Dissector_DissectionDate, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -28.0870  -4.1492   0.1669   3.9676  28.8474 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                  4.1376     1.8270   2.265 0.024474 *  
#   HPC_Dissector_DissectionDateEHB_12/11/2020  -1.8696     3.0979  -0.604 0.546768    
# HPC_Dissector_DissectionDateEHB_12/14/2020  -2.9827     2.9835  -1.000 0.318507    
# HPC_Dissector_DissectionDateEHB_12/4/2020   -8.8928     2.6814  -3.317 0.001061 ** 
#   HPC_Dissector_DissectionDateEHB_12/7/2020   -1.9637     2.8089  -0.699 0.485205    
# HPC_Dissector_DissectionDateEHB_12/8/2020   -6.5959     2.6814  -2.460 0.014641 *  
#   HPC_Dissector_DissectionDateEHB_12/9/2020    1.0137     2.8089   0.361 0.718512    
# HPC_Dissector_DissectionDateKA_12/1/2020    -9.8600     2.8089  -3.510 0.000539 ***
#   HPC_Dissector_DissectionDateKA_12/10/2020    1.2254     2.4441   0.501 0.616598    
# HPC_Dissector_DissectionDateKA_12/11/2020  -10.6574     2.3694  -4.498 1.09e-05 ***
#   HPC_Dissector_DissectionDateKA_12/14/2020  -26.4484     2.8888  -9.155  < 2e-16 ***
#   HPC_Dissector_DissectionDateKA_12/2/2020    -5.7049     2.5431  -2.243 0.025842 *  
#   HPC_Dissector_DissectionDateKA_12/3/2020    -5.0695     2.6296  -1.928 0.055112 .  
# HPC_Dissector_DissectionDateKA_12/4/2020    -6.7340     2.5838  -2.606 0.009758 ** 
#   HPC_Dissector_DissectionDateKA_12/7/2020    -0.5044     2.4738  -0.204 0.838604    
# HPC_Dissector_DissectionDateKA_12/8/2020     1.4406     2.3922   0.602 0.547638    
# HPC_Dissector_DissectionDateKA_12/9/2020     3.4792     2.4441   1.424 0.155948    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 7.076 on 228 degrees of freedom
# Multiple R-squared:  0.4667,	Adjusted R-squared:  0.4293 
# F-statistic: 12.47 on 16 and 228 DF,  p-value: < 2.2e-16
# 


summary.lm(lm(PC5 ~ HPC_Dissector_DissectionDate_Collapsed, data=F2_PCA_wMetaData))
# 
# Call:
#   lm(formula = PC5 ~ HPC_Dissector_DissectionDate_Collapsed, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -15.9521  -4.3533  -0.7439   3.7762  21.3852 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                    2.45088    1.04303   2.350 0.019617 *  
#   HPC_Dissector_DissectionDate_CollapsedBatch1  -5.92894    2.08606  -2.842 0.004876 ** 
#   HPC_Dissector_DissectionDate_CollapsedBatch10 -1.79416    1.48474  -1.208 0.228113    
# HPC_Dissector_DissectionDate_CollapsedBatch11 -0.48196    1.55364  -0.310 0.756677    
# HPC_Dissector_DissectionDate_CollapsedBatch2  -0.02873    2.22375  -0.013 0.989703    
# HPC_Dissector_DissectionDate_CollapsedBatch4  -1.98786    2.52814  -0.786 0.432492    
# HPC_Dissector_DissectionDate_CollapsedBatch5  -0.26149    2.40878  -0.109 0.913647    
# HPC_Dissector_DissectionDate_CollapsedBatch6  -1.96325    2.22375  -0.883 0.378219    
# HPC_Dissector_DissectionDate_CollapsedBatch7  -3.65682    1.93383  -1.891 0.059863 .  
# HPC_Dissector_DissectionDate_CollapsedBatch8  -3.81043    2.02942  -1.878 0.061680 .  
# HPC_Dissector_DissectionDate_CollapsedBatch9  -5.32897    1.36881  -3.893 0.000129 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 6.514 on 234 degrees of freedom
# Multiple R-squared:  0.103,	Adjusted R-squared:  0.0647 
# F-statistic: 2.688 on 10 and 234 DF,  p-value: 0.003915


summary.lm(lm(PC5 ~ HPC_Dissector_DissectionDate, data=F2_PCA_wMetaData))

# 
# Call:
#   lm(formula = PC5 ~ HPC_Dissector_DissectionDate, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -14.8751  -3.9914  -0.7836   3.6754  17.8240 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                  3.2189     1.6107   1.999  0.04685 *  
#   HPC_Dissector_DissectionDateEHB_12/11/2020  -2.7559     2.7310  -1.009  0.31399    
# HPC_Dissector_DissectionDateEHB_12/14/2020  -1.0295     2.6302  -0.391  0.69584    
# HPC_Dissector_DissectionDateEHB_12/4/2020   -6.6970     2.3638  -2.833  0.00502 ** 
#   HPC_Dissector_DissectionDateEHB_12/7/2020   -0.7968     2.4762  -0.322  0.74792    
# HPC_Dissector_DissectionDateEHB_12/8/2020   -3.7340     2.3638  -1.580  0.11557    
# HPC_Dissector_DissectionDateEHB_12/9/2020    1.6898     2.4762   0.682  0.49567    
# HPC_Dissector_DissectionDateKA_12/1/2020    -2.7313     2.4762  -1.103  0.27119    
# HPC_Dissector_DissectionDateKA_12/10/2020    0.9989     2.1546   0.464  0.64335    
# HPC_Dissector_DissectionDateKA_12/11/2020   -0.2857     2.0888  -0.137  0.89132    
# HPC_Dissector_DissectionDateKA_12/14/2020   -3.3715     2.5467  -1.324  0.18687    
# HPC_Dissector_DissectionDateKA_12/2/2020    -4.4249     2.2419  -1.974  0.04962 *  
#   HPC_Dissector_DissectionDateKA_12/3/2020    -4.5785     2.3181  -1.975  0.04947 *  
#   HPC_Dissector_DissectionDateKA_12/4/2020    -6.2419     2.2778  -2.740  0.00662 ** 
#   HPC_Dissector_DissectionDateKA_12/7/2020    -2.9538     2.1808  -1.354  0.17693    
# HPC_Dissector_DissectionDateKA_12/8/2020    -8.6877     2.1088  -4.120 5.31e-05 ***
#   HPC_Dissector_DissectionDateKA_12/9/2020    -6.1234     2.1546  -2.842  0.00489 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 6.238 on 228 degrees of freedom
# Multiple R-squared:  0.1984,	Adjusted R-squared:  0.1422 
# F-statistic: 3.528 on 16 and 228 DF,  p-value: 1.118e-05


table(F2_PCA_wMetaData$SequencingBatch,F2_PCA_wMetaData$HPC_Dissector_DissectionDate)

# 
#      EHB_12/10/2020 EHB_12/11/2020 EHB_12/14/2020 EHB_12/4/2020 EHB_12/7/2020 EHB_12/8/2020 EHB_12/9/2020 KA_12/1/2020 KA_12/10/2020 KA_12/11/2020 KA_12/14/2020 KA_12/2/2020 KA_12/3/2020 KA_12/4/2020 KA_12/7/2020 KA_12/8/2020
# 1              0              0              0            13             9             0             0           10             0             0             0           16           14           15           18            0
# 2             15              0              0             0             2            13            11            0            13             0             0            0            0            0            0           21
# 3              0              8              9             0             0             0             0            1             6            22            10            0            0            0            0            0
# 
# KA_12/9/2020
# 1            0
# 2           19
# 3            0

#Hmm... particular dissector-dissection dates definitely map on to particular sequencing batches.

##############################

#We discovered there are other RNA metrics included in the Sequencing Core output that might provide insights into the noise in the data:

list.files()

RNAMetricsForAllF2s<-read.csv("RNAMetrics_ForR_AllAnimals.csv",header=TRUE,stringsAsFactors = FALSE)

colnames(F2_PCA_wMetaData)

colnames(F2_lcpm_noLowHits)

str(RNAMetricsForAllF2s)

RNAMetricsForAllF2s$?..Library.ID%in%colnames(F2_lcpm_noLowHits)

RNAMetricsForFinal_F2s<-RNAMetricsForAllF2s[RNAMetricsForAllF2s$?..Library.ID%in%colnames(F2_lcpm_noLowHits),]

str(F2_PCA_wMetaData)

F2_PCA_wMetaData<-cbind.data.frame(F2_PCA_wMetaData,RNAMetricsForFinal_F2s)
