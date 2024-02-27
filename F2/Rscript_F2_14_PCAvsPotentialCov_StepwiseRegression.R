#F2 HC RNA-Seq Dataset
#14_Screening for Noise and Potential Confounds: using stepwise regression for PCA Output vs. Potential Covariates
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.


#Using stepwise regression to help us choose between potential co-variates:

library(stats)

step(object=lm(PC1~1, data=F2_PCA_wMetaData), scope=(~Sex+RibosomePerc+RNAconc+RIN+LibrarySize+SequencingBatch+Dissector_AsNumeric+HPC_Dissection_Date_Collapsed+HPC_Dissection_Day+STGT_experience+DV200+Percent.Intronic+Percent.Intergenic+Percent.Coding+Percent.UTR+MEDIAN_CV_COVERAGE), scale = 0, direction = "both", trace = 1, keep = NULL, steps = 1000, k=log(245))
# Start:  AIC=1870.49
# PC1 ~ 1
# 
# Df Sum of Sq    RSS    AIC
# + Percent.Intergenic             1    444884  50689 1317.4
# + RibosomePerc                   1    377282 118291 1525.0
# + MEDIAN_CV_COVERAGE             1    311100 184473 1633.9
# + RNAconc                        1    246685 248888 1707.3
# + Percent.Coding                 1    207724 287849 1742.9
# + Dissector_AsNumeric            1    134706 360867 1798.3
# + Percent.Intronic               1    111087 384487 1813.8
# + SequencingBatch                2     75753 419820 1840.8
# + Percent.UTR                    1     53013 442560 1848.3
# + HPC_Dissection_Day             1     44373 451201 1853.0
# + HPC_Dissection_Date_Collapsed  5     81906 413667 1853.7
# + STGT_experience                1     15764 479809 1868.1
# + RIN                            1     14240 481333 1868.8
# + LibrarySize                    1     13731 481842 1869.1
# + DV200                          1     11893 483680 1870.0
# <none>                                       495573 1870.5
# + Sex                            1      1058 494516 1875.5
# 
# Step:  AIC=1317.39
# PC1 ~ Percent.Intergenic
# 
# Df Sum of Sq    RSS    AIC
# + RibosomePerc                   1     24198  26490 1163.9
# + Percent.Coding                 1     20092  30596 1199.2
# + Percent.UTR                    1     10564  40125 1265.6
# + MEDIAN_CV_COVERAGE             1      5246  45443 1296.1
# + RNAconc                        1      4961  45727 1297.7
# + Sex                            1      2585  48104 1310.1
# + Dissector_AsNumeric            1      2575  48113 1310.1
# <none>                                        50689 1317.4
# + SequencingBatch                2      1934  48754 1318.9
# + DV200                          1       800  49889 1319.0
# + Percent.Intronic               1       563  50126 1320.2
# + STGT_experience                1       156  50533 1322.1
# + HPC_Dissection_Day             1       130  50559 1322.3
# + LibrarySize                    1       108  50581 1322.4
# + RIN                            1        14  50675 1322.8
# + HPC_Dissection_Date_Collapsed  5      1113  49576 1339.5
# - Percent.Intergenic             1    444884 495573 1870.5
# 
# Step:  AIC=1163.91
# PC1 ~ Percent.Intergenic + RibosomePerc
# 
# Df Sum of Sq    RSS    AIC
# + Percent.Coding                 1      4556  21934 1123.2
# + Percent.UTR                    1      3790  22701 1131.6
# + RIN                            1      1305  25185 1157.0
# + SequencingBatch                2      1225  25266 1163.3
# <none>                                        26490 1163.9
# + Dissector_AsNumeric            1       473  26017 1165.0
# + HPC_Dissection_Date_Collapsed  5      2625  23865 1165.8
# + MEDIAN_CV_COVERAGE             1       378  26112 1165.9
# + STGT_experience                1       259  26232 1167.0
# + RNAconc                        1       232  26259 1167.3
# + HPC_Dissection_Day             1       215  26276 1167.4
# + Sex                            1       208  26283 1167.5
# + Percent.Intronic               1       135  26356 1168.2
# + DV200                          1       114  26377 1168.4
# + LibrarySize                    1        42  26448 1169.0
# - RibosomePerc                   1     24198  50689 1317.4
# - Percent.Intergenic             1     91801 118291 1525.0
# 
# Step:  AIC=1123.17
# PC1 ~ Percent.Intergenic + RibosomePerc + Percent.Coding
# 
# Df Sum of Sq   RSS    AIC
# + SequencingBatch                2      2958 18976 1098.7
# + MEDIAN_CV_COVERAGE             1      2164 19770 1103.2
# + HPC_Dissection_Date_Collapsed  5      3780 18155 1104.3
# + RIN                            1      1225 20709 1114.6
# + Dissector_AsNumeric            1       714 21220 1120.6
# <none>                                       21934 1123.2
# + STGT_experience                1       418 21517 1124.0
# + Sex                            1       153 21781 1127.0
# + LibrarySize                    1       118 21817 1127.3
# + RNAconc                        1       112 21823 1127.4
# + Percent.UTR                    1        40 21895 1128.2
# + DV200                          1        27 21907 1128.4
# + Percent.Intronic               1        16 21918 1128.5
# + HPC_Dissection_Day             1        12 21922 1128.5
# - Percent.Coding                 1      4556 26490 1163.9
# - RibosomePerc                   1      8662 30596 1199.2
# - Percent.Intergenic             1     48861 70795 1404.8
# 
# Step:  AIC=1098.68
# PC1 ~ Percent.Intergenic + RibosomePerc + Percent.Coding + SequencingBatch
# 
# Df Sum of Sq   RSS    AIC
# + RIN                            1      1614 17362 1082.4
# + MEDIAN_CV_COVERAGE             1      1157 17819 1088.8
# + Dissector_AsNumeric            1       690 18286 1095.1
# <none>                                       18976 1098.7
# + RNAconc                        1       398 18578 1099.0
# + Sex                            1       149 18827 1102.2
# + DV200                          1        76 18901 1103.2
# + Percent.Intronic               1        51 18925 1103.5
# + Percent.UTR                    1        25 18951 1103.9
# + HPC_Dissection_Day             1        16 18960 1104.0
# + LibrarySize                    1         5 18971 1104.1
# + STGT_experience                1         1 18976 1104.2
# + HPC_Dissection_Date_Collapsed  5      1026 17951 1112.6
# - SequencingBatch                2      2958 21934 1123.2
# - Percent.Coding                 1      6289 25266 1163.3
# - RibosomePerc                   1     10689 29666 1202.7
# - Percent.Intergenic             1     48447 67423 1403.8
# 
# Step:  AIC=1082.4
# PC1 ~ Percent.Intergenic + RibosomePerc + Percent.Coding + SequencingBatch + 
#   RIN
# 
# Df Sum of Sq   RSS    AIC
# + MEDIAN_CV_COVERAGE             1       992 16370 1073.5
# + Dissector_AsNumeric            1       740 16622 1077.2
# <none>                                       17362 1082.4
# + RNAconc                        1       280 17082 1083.9
# + DV200                          1       263 17099 1084.2
# + Sex                            1       134 17228 1086.0
# + STGT_experience                1        42 17320 1087.3
# + HPC_Dissection_Day             1        25 17337 1087.5
# + Percent.UTR                    1        22 17340 1087.6
# + Percent.Intronic               1         9 17353 1087.8
# + LibrarySize                    1         6 17356 1087.8
# + HPC_Dissection_Date_Collapsed  5      1482 15881 1088.0
# - RIN                            1      1614 18976 1098.7
# - SequencingBatch                2      3347 20709 1114.6
# - Percent.Coding                 1      6142 23504 1151.1
# - RibosomePerc                   1     12114 29476 1206.6
# - Percent.Intergenic             1     42849 60211 1381.6
# 
# Step:  AIC=1073.48
# PC1 ~ Percent.Intergenic + RibosomePerc + Percent.Coding + SequencingBatch + 
#   RIN + MEDIAN_CV_COVERAGE
# 
# Df Sum of Sq   RSS    AIC
# + DV200                          1       480 15890 1071.7
# + Dissector_AsNumeric            1       479 15890 1071.7
# <none>                                       16370 1073.5
# + Sex                            1       199 16171 1076.0
# + HPC_Dissection_Date_Collapsed  5      1514 14856 1077.2
# + RNAconc                        1        43 16326 1078.3
# + STGT_experience                1        16 16354 1078.8
# + Percent.Intronic               1         7 16362 1078.9
# + Percent.UTR                    1         2 16368 1079.0
# + HPC_Dissection_Day             1         2 16368 1079.0
# + LibrarySize                    1         0 16370 1079.0
# - MEDIAN_CV_COVERAGE             1       992 17362 1082.4
# - RIN                            1      1450 17819 1088.8
# - SequencingBatch                2      2281 18651 1094.4
# - RibosomePerc                   1      6051 22420 1145.0
# - Percent.Coding                 1      7078 23448 1156.0
# - Percent.Intergenic             1     42888 59258 1383.2
# 
# Step:  AIC=1071.7
# PC1 ~ Percent.Intergenic + RibosomePerc + Percent.Coding + SequencingBatch + 
#   RIN + MEDIAN_CV_COVERAGE + DV200
# 
# Df Sum of Sq   RSS    AIC
# + Dissector_AsNumeric            1       518 15372 1069.1
# <none>                                       15890 1071.7
# + Sex                            1       249 15641 1073.3
# - DV200                          1       480 16370 1073.5
# + HPC_Dissection_Date_Collapsed  5      1533 14357 1074.3
# + STGT_experience                1        76 15814 1076.0
# + RNAconc                        1        42 15848 1076.6
# + LibrarySize                    1         5 15885 1077.1
# + Percent.UTR                    1         4 15886 1077.1
# + Percent.Intronic               1         1 15889 1077.2
# + HPC_Dissection_Day             1         0 15890 1077.2
# - SequencingBatch                2      1560 17450 1083.6
# - MEDIAN_CV_COVERAGE             1      1209 17099 1084.2
# - RIN                            1      1922 17812 1094.2
# - RibosomePerc                   1      5840 21730 1142.9
# - Percent.Coding                 1      7557 23447 1161.5
# - Percent.Intergenic             1     43360 59250 1388.6
# 
# Step:  AIC=1069.07
# PC1 ~ Percent.Intergenic + RibosomePerc + Percent.Coding + SequencingBatch + 
#   RIN + MEDIAN_CV_COVERAGE + DV200 + Dissector_AsNumeric
# 
# Df Sum of Sq   RSS    AIC
# <none>                                       15372 1069.1
# + Sex                            1       214 15158 1071.1
# - Dissector_AsNumeric            1       518 15890 1071.7
# - DV200                          1       519 15890 1071.7
# + HPC_Dissection_Day             1        63 15309 1073.6
# + STGT_experience                1        16 15355 1074.3
# + LibrarySize                    1        13 15359 1074.4
# + Percent.Intronic               1         4 15368 1074.5
# + RNAconc                        1         2 15370 1074.5
# + Percent.UTR                    1         1 15371 1074.6
# + HPC_Dissection_Date_Collapsed  5      1163 14209 1077.3
# - MEDIAN_CV_COVERAGE             1       922 16293 1077.8
# - SequencingBatch                2      1634 17006 1082.8
# - RIN                            1      2020 17392 1093.8
# - RibosomePerc                   1      5594 20966 1139.6
# - Percent.Coding                 1      7402 22774 1159.9
# - Percent.Intergenic             1     41446 56818 1383.9
# 
# Call:
#   lm(formula = PC1 ~ Percent.Intergenic + RibosomePerc + Percent.Coding + 
#        SequencingBatch + RIN + MEDIAN_CV_COVERAGE + DV200 + Dissector_AsNumeric, 
#      data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)   Percent.Intergenic         RibosomePerc       Percent.Coding     SequencingBatch2     SequencingBatch3                  RIN   MEDIAN_CV_COVERAGE                DV200  Dissector_AsNumeric  
# 11.440            -2406.502             2552.830             -764.373                6.243               -1.158              -17.658              208.900                7.307                3.846  

#k=log(nobs(PC1~1, data=F2_PCA_wMetaData))

summary.lm(lm(PC1 ~ Percent.Intergenic + RibosomePerc + Percent.Coding + SequencingBatch + RIN + MEDIAN_CV_COVERAGE + DV200 + Dissector_AsNumeric, data = F2_PCA_wMetaData))

# Call:
#   lm(formula = PC1 ~ Percent.Intergenic + RibosomePerc + Percent.Coding + 
#        SequencingBatch + RIN + MEDIAN_CV_COVERAGE + DV200 + Dissector_AsNumeric, 
#      data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -24.2908  -4.9819  -0.1093   5.4531  19.9787 
# 
# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            11.440    235.054   0.049  0.96122    
# Percent.Intergenic  -2406.502     95.603 -25.172  < 2e-16 ***
#   RibosomePerc         2552.830    276.048   9.248  < 2e-16 ***
#   Percent.Coding       -764.373     71.854 -10.638  < 2e-16 ***
#   SequencingBatch2        6.243      1.376   4.536 9.13e-06 ***
#   SequencingBatch3       -1.158      1.607  -0.720  0.47212    
# RIN                   -17.658      3.178  -5.557 7.42e-08 ***
#   MEDIAN_CV_COVERAGE    208.900     55.655   3.753  0.00022 ***
#   DV200                   7.307      2.595   2.816  0.00527 ** 
#   Dissector_AsNumeric     3.846      1.366   2.815  0.00529 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 8.088 on 235 degrees of freedom
# Multiple R-squared:  0.969,	Adjusted R-squared:  0.9678 
# F-statistic: 815.7 on 9 and 235 DF,  p-value: < 2.2e-16
#            


step(object=lm(PC2~1, data=F2_PCA_wMetaData), scope=(~Sex+RibosomePerc+RNAconc+RIN+LibrarySize+SequencingBatch+Dissector_AsNumeric+HPC_Dissection_Date_Collapsed+HPC_Dissection_Day+STGT_experience+DV200+Percent.Intronic+Percent.Intergenic+Percent.Coding+Percent.UTR+MEDIAN_CV_COVERAGE), scale = 0, direction = "both", trace = 1, keep = NULL, steps = 1000, k=log(245))
# Start:  AIC=1279.13
# PC2 ~ 1
# 
# Df Sum of Sq   RSS    AIC
# + HPC_Dissection_Date_Collapsed  5    7133.7 37210 1263.7
# + HPC_Dissection_Day             1    2785.0 41559 1268.7
# + SequencingBatch                2    2829.9 41514 1274.0
# + RNAconc                        1    1569.8 42774 1275.8
# + Percent.Intronic               1    1442.9 42901 1276.5
# + DV200                          1    1303.3 43041 1277.3
# + RibosomePerc                   1    1176.0 43168 1278.0
# + RIN                            1    1109.1 43235 1278.4
# <none>                                       44344 1279.1
# + LibrarySize                    1     691.8 43652 1280.8
# + Dissector_AsNumeric            1     352.1 43992 1282.7
# + Sex                            1     242.5 44101 1283.3
# + Percent.Coding                 1     146.9 44197 1283.8
# + Percent.Intergenic             1     114.0 44230 1284.0
# + MEDIAN_CV_COVERAGE             1      52.7 44291 1284.3
# + STGT_experience                1      42.1 44302 1284.4
# + Percent.UTR                    1       0.1 44344 1284.6
# 
# Step:  AIC=1263.66
# PC2 ~ HPC_Dissection_Date_Collapsed
# 
# Df Sum of Sq   RSS    AIC
# + RibosomePerc                   1    1200.8 36009 1261.1
# <none>                                       37210 1263.7
# + Percent.Intronic               1     738.8 36471 1264.2
# + RIN                            1     737.2 36473 1264.3
# + RNAconc                        1     483.7 36726 1266.0
# + DV200                          1     391.7 36818 1266.6
# + LibrarySize                    1     336.7 36873 1266.9
# + HPC_Dissection_Day             1     319.8 36890 1267.0
# + Sex                            1     200.1 37010 1267.8
# + MEDIAN_CV_COVERAGE             1     125.0 37085 1268.3
# + Dissector_AsNumeric            1     123.1 37087 1268.3
# + STGT_experience                1      74.5 37136 1268.7
# + Percent.UTR                    1      22.6 37188 1269.0
# + Percent.Intergenic             1       7.3 37203 1269.1
# + Percent.Coding                 1       4.4 37206 1269.1
# + SequencingBatch                2     384.8 36825 1272.1
# - HPC_Dissection_Date_Collapsed  5    7133.7 44344 1279.1
# 
# Step:  AIC=1261.13
# PC2 ~ HPC_Dissection_Date_Collapsed + RibosomePerc
# 
# Df Sum of Sq   RSS    AIC
# + Percent.Intergenic             1    1325.8 34684 1257.4
# + RIN                            1     810.6 35199 1261.0
# <none>                                       36009 1261.1
# + MEDIAN_CV_COVERAGE             1     541.7 35468 1262.9
# - RibosomePerc                   1    1200.8 37210 1263.7
# + HPC_Dissection_Day             1     356.8 35652 1264.2
# + DV200                          1     304.0 35705 1264.5
# + LibrarySize                    1     292.9 35716 1264.6
# + Percent.Intronic               1     276.2 35733 1264.7
# + Percent.Coding                 1     234.8 35774 1265.0
# + Percent.UTR                    1     125.2 35884 1265.8
# + Dissector_AsNumeric            1     117.0 35892 1265.8
# + Sex                            1     115.9 35893 1265.8
# + STGT_experience                1      49.5 35960 1266.3
# + RNAconc                        1      37.2 35972 1266.4
# + SequencingBatch                2     574.1 35435 1268.2
# - HPC_Dissection_Date_Collapsed  5    7158.6 43168 1278.0
# 
# Step:  AIC=1257.44
# PC2 ~ HPC_Dissection_Date_Collapsed + RibosomePerc + Percent.Intergenic
# 
# Df Sum of Sq   RSS    AIC
# + RIN                            1    1838.8 32845 1249.6
# + Percent.Intronic               1    1646.4 33037 1251.0
# + Percent.Coding                 1     797.7 33886 1257.2
# <none>                                       34684 1257.4
# + HPC_Dissection_Day             1     459.3 34224 1259.7
# + DV200                          1     394.0 34290 1260.1
# - Percent.Intergenic             1    1325.8 36009 1261.1
# + Percent.UTR                    1     215.3 34468 1261.4
# + LibrarySize                    1     149.3 34534 1261.9
# + RNAconc                        1      97.2 34586 1262.2
# + STGT_experience                1      65.9 34618 1262.5
# + MEDIAN_CV_COVERAGE             1      53.5 34630 1262.6
# + Sex                            1      11.5 34672 1262.9
# + Dissector_AsNumeric            1       2.2 34681 1262.9
# - HPC_Dissection_Date_Collapsed  5    5050.9 39734 1263.2
# + SequencingBatch                2     491.8 34192 1264.9
# - RibosomePerc                   1    2519.3 37203 1269.1
# 
# Step:  AIC=1249.59
# PC2 ~ HPC_Dissection_Date_Collapsed + RibosomePerc + Percent.Intergenic + 
#   RIN
# 
# Df Sum of Sq   RSS    AIC
# + Percent.Intronic               1    1025.2 31820 1247.3
# + Percent.Coding                 1     818.1 32027 1248.9
# <none>                                       32845 1249.6
# - HPC_Dissection_Date_Collapsed  5    3945.7 36790 1249.9
# + Percent.UTR                    1     320.6 32524 1252.7
# + HPC_Dissection_Day             1     308.2 32537 1252.8
# + STGT_experience                1     170.1 32675 1253.8
# + LibrarySize                    1     137.5 32707 1254.1
# + MEDIAN_CV_COVERAGE             1      87.2 32758 1254.4
# + RNAconc                        1      30.5 32814 1254.9
# + DV200                          1      13.8 32831 1255.0
# + Sex                            1       8.7 32836 1255.0
# + Dissector_AsNumeric            1       0.0 32845 1255.1
# - RIN                            1    1838.8 34684 1257.4
# + SequencingBatch                2     368.6 32476 1257.8
# - Percent.Intergenic             1    2354.0 35199 1261.0
# - RibosomePerc                   1    3622.5 36467 1269.7
# 
# Step:  AIC=1247.32
# PC2 ~ HPC_Dissection_Date_Collapsed + RibosomePerc + Percent.Intergenic + 
#   RIN + Percent.Intronic
# 
# Df Sum of Sq   RSS    AIC
# - HPC_Dissection_Date_Collapsed  5    3110.3 34930 1242.7
# + HPC_Dissection_Day             1     782.5 31037 1246.7
# <none>                                       31820 1247.3
# + Percent.Coding                 1     624.9 31195 1248.0
# + Percent.UTR                    1     609.4 31210 1248.1
# - Percent.Intronic               1    1025.2 32845 1249.6
# - RIN                            1    1217.6 33037 1251.0
# + STGT_experience                1     228.7 31591 1251.1
# + RNAconc                        1     216.1 31603 1251.2
# + LibrarySize                    1      56.7 31763 1252.4
# + Dissector_AsNumeric            1      36.6 31783 1252.5
# + Sex                            1      23.6 31796 1252.6
# + MEDIAN_CV_COVERAGE             1      14.5 31805 1252.7
# + DV200                          1      10.6 31809 1252.7
# + SequencingBatch                2     259.8 31560 1256.3
# - Percent.Intergenic             1    3326.3 35146 1266.2
# - RibosomePerc                   1    4026.2 35846 1271.0
# 
# Step:  AIC=1242.67
# PC2 ~ RibosomePerc + Percent.Intergenic + RIN + Percent.Intronic
# 
# Df Sum of Sq   RSS    AIC
# + HPC_Dissection_Day             1    1442.4 33487 1237.8
# + Percent.Coding                 1    1073.3 33857 1240.5
# + Percent.UTR                    1    1055.6 33874 1240.7
# <none>                                       34930 1242.7
# + RNAconc                        1     701.9 34228 1243.2
# + STGT_experience                1     525.5 34404 1244.5
# + Dissector_AsNumeric            1     291.5 34638 1246.1
# + LibrarySize                    1     228.6 34701 1246.6
# + HPC_Dissection_Date_Collapsed  5    3110.3 31820 1247.3
# + DV200                          1      42.1 34888 1247.9
# + MEDIAN_CV_COVERAGE             1      41.8 34888 1247.9
# - RIN                            1    1569.5 36499 1247.9
# + Sex                            1       4.8 34925 1248.1
# - Percent.Intronic               1    1860.6 36790 1249.9
# + SequencingBatch                2     503.9 34426 1250.1
# - RibosomePerc                   1    6184.9 41115 1277.1
# - Percent.Intergenic             1    6850.9 41781 1281.0
# 
# Step:  AIC=1237.84
# PC2 ~ RibosomePerc + Percent.Intergenic + RIN + Percent.Intronic + 
#   HPC_Dissection_Day
# 
# Df Sum of Sq   RSS    AIC
# <none>                                       33487 1237.8
# + Percent.Coding                 1     740.9 32747 1237.9
# + Percent.UTR                    1     723.8 32764 1238.0
# - RIN                            1    1304.3 34792 1241.7
# + RNAconc                        1     147.1 33340 1242.3
# - Percent.Intronic               1    1388.8 34876 1242.3
# - HPC_Dissection_Day             1    1442.4 34930 1242.7
# + STGT_experience                1      49.9 33438 1243.0
# + LibrarySize                    1      36.2 33451 1243.1
# + Sex                            1      13.3 33474 1243.2
# + DV200                          1       7.1 33480 1243.3
# + MEDIAN_CV_COVERAGE             1       4.0 33483 1243.3
# + Dissector_AsNumeric            1       0.2 33487 1243.3
# + SequencingBatch                2     455.7 33032 1245.5
# + HPC_Dissection_Date_Collapsed  5    2450.4 31037 1246.7
# - Percent.Intergenic             1    4789.8 38277 1265.1
# - RibosomePerc                   1    5843.4 39331 1271.7
# 
# Call:
#   lm(formula = PC2 ~ RibosomePerc + Percent.Intergenic + RIN + 
#        Percent.Intronic + HPC_Dissection_Day, data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)        RibosomePerc  Percent.Intergenic                 RIN    Percent.Intronic  HPC_Dissection_Day  
# -62.9587          -1803.5497           -477.9272             12.1489            877.9266              0.7368

summary.lm(lm(PC2 ~ RibosomePerc + Percent.Intergenic + RIN + Percent.Intronic + HPC_Dissection_Day, data = F2_PCA_wMetaData))

# Call:
#   lm(formula = PC2 ~ RibosomePerc + Percent.Intergenic + RIN + 
#        Percent.Intronic + HPC_Dissection_Day, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -29.691  -7.921  -0.185   8.501  35.078 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          -62.9587    32.9584  -1.910  0.05730 .  
# RibosomePerc       -1803.5497   279.2791  -6.458 5.88e-10 ***
#   Percent.Intergenic  -477.9272    81.7423  -5.847 1.64e-08 ***
#   RIN                   12.1489     3.9819   3.051  0.00254 ** 
#   Percent.Intronic     877.9266   278.8546   3.148  0.00185 ** 
#   HPC_Dissection_Day     0.7368     0.2296   3.208  0.00152 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 11.84 on 239 degrees of freedom
# Multiple R-squared:  0.2448,	Adjusted R-squared:  0.229 
# F-statistic:  15.5 on 5 and 239 DF,  p-value: 3.356e-13
# 

step(object=lm(PC2~1, data=F2_PCA_wMetaData), scope=(~Sex+RibosomePerc+RNAconc+RIN+LibrarySize+SequencingBatch+HPC_Dissector_DissectionDate_Collapsed+STGT_experience+DV200+Percent.Intronic+Percent.Intergenic+Percent.Coding+Percent.UTR+MEDIAN_CV_COVERAGE), scale = 0, direction = "both", trace = 1, keep = NULL, steps = 1000, k=log(245))
# Start:  AIC=1279.13
# PC2 ~ 1
# 
# Df Sum of Sq   RSS    AIC
# + HPC_Dissector_DissectionDate_Collapsed 10   10154.0 34190 1270.4
# + SequencingBatch                         2    2829.9 41514 1274.0
# + RNAconc                                 1    1569.8 42774 1275.8
# + Percent.Intronic                        1    1442.9 42901 1276.5
# + DV200                                   1    1303.3 43041 1277.3
# + RibosomePerc                            1    1176.0 43168 1278.0
# + RIN                                     1    1109.1 43235 1278.4
# <none>                                                44344 1279.1
# + LibrarySize                             1     691.8 43652 1280.8
# + Sex                                     1     242.5 44101 1283.3
# + Percent.Coding                          1     146.9 44197 1283.8
# + Percent.Intergenic                      1     114.0 44230 1284.0
# + MEDIAN_CV_COVERAGE                      1      52.7 44291 1284.3
# + STGT_experience                         1      42.1 44302 1284.4
# + Percent.UTR                             1       0.1 44344 1284.6
# 
# Step:  AIC=1270.43
# PC2 ~ HPC_Dissector_DissectionDate_Collapsed
# 
# Df Sum of Sq   RSS    AIC
# + RibosomePerc                            1    1929.8 32260 1261.7
# <none>                                                34190 1270.4
# + DV200                                   1     508.7 33681 1272.3
# + RIN                                     1     448.9 33741 1272.7
# + RNAconc                                 1     390.2 33800 1273.1
# + Percent.Intronic                        1     371.2 33819 1273.3
# + LibrarySize                             1     268.5 33921 1274.0
# + Sex                                     1     265.5 33924 1274.0
# + MEDIAN_CV_COVERAGE                      1     105.5 34084 1275.2
# + Percent.Coding                          1      67.8 34122 1275.4
# + STGT_experience                         1      64.8 34125 1275.5
# + Percent.Intergenic                      1      14.5 34175 1275.8
# + Percent.UTR                             1       0.8 34189 1275.9
# - HPC_Dissector_DissectionDate_Collapsed 10   10154.0 44344 1279.1
# + SequencingBatch                         2     126.7 34063 1280.5
# 
# Step:  AIC=1261.7
# PC2 ~ HPC_Dissector_DissectionDate_Collapsed + RibosomePerc
# 
# Df Sum of Sq   RSS    AIC
# + Percent.Intergenic                      1    1786.8 30473 1253.2
# <none>                                                32260 1261.7
# + RIN                                     1     556.1 31704 1262.9
# + MEDIAN_CV_COVERAGE                      1     417.5 31843 1264.0
# + DV200                                   1     411.2 31849 1264.0
# + Percent.Coding                          1     380.6 31879 1264.3
# + LibrarySize                             1     202.7 32057 1265.7
# + STGT_experience                         1     191.0 32069 1265.7
# + Sex                                     1     150.7 32109 1266.0
# + RNAconc                                 1     127.3 32133 1266.2
# + Percent.UTR                             1     118.4 32142 1266.3
# + Percent.Intronic                        1      24.2 32236 1267.0
# + SequencingBatch                         2     492.6 31768 1268.9
# - RibosomePerc                            1    1929.8 34190 1270.4
# - HPC_Dissector_DissectionDate_Collapsed 10   10907.8 43168 1278.0
# 
# Step:  AIC=1253.24
# PC2 ~ HPC_Dissector_DissectionDate_Collapsed + RibosomePerc + 
#   Percent.Intergenic
# 
# Df Sum of Sq   RSS    AIC
# + RIN                                     1    1566.4 28907 1245.8
# + Percent.Intronic                        1     817.0 29656 1252.1
# <none>                                                30473 1253.2
# + Percent.Coding                          1     625.4 29848 1253.7
# + DV200                                   1     499.4 29974 1254.7
# + STGT_experience                         1     258.4 30215 1256.7
# + Percent.UTR                             1     244.8 30229 1256.8
# + LibrarySize                             1      67.0 30406 1258.2
# + Sex                                     1      11.8 30461 1258.6
# + MEDIAN_CV_COVERAGE                      1       1.2 30472 1258.7
# + RNAconc                                 1       0.7 30473 1258.7
# + SequencingBatch                         2     531.2 29942 1259.9
# - Percent.Intergenic                      1    1786.8 32260 1261.7
# - HPC_Dissector_DissectionDate_Collapsed 10    9261.1 39734 1263.2
# - RibosomePerc                            1    3702.1 34175 1275.8
# 
# Step:  AIC=1245.81
# PC2 ~ HPC_Dissector_DissectionDate_Collapsed + RibosomePerc + 
#   Percent.Intergenic + RIN
# 
# Df Sum of Sq   RSS    AIC
# + Percent.Coding                          1     671.1 28236 1245.6
# <none>                                                28907 1245.8
# + Percent.Intronic                        1     443.6 28463 1247.5
# + Percent.UTR                             1     355.8 28551 1248.3
# + STGT_experience                         1     314.2 28593 1248.6
# - HPC_Dissector_DissectionDate_Collapsed 10    7883.5 36790 1249.9
# + LibrarySize                             1      85.6 28821 1250.6
# + Sex                                     1       5.1 28902 1251.3
# + MEDIAN_CV_COVERAGE                      1       0.1 28907 1251.3
# + DV200                                   1       0.0 28907 1251.3
# + RNAconc                                 1       0.0 28907 1251.3
# + SequencingBatch                         2     541.7 28365 1252.2
# - RIN                                     1    1566.4 30473 1253.2
# - Percent.Intergenic                      1    2797.2 31704 1262.9
# - RibosomePerc                            1    4761.0 33668 1277.7
# 
# Step:  AIC=1245.56
# PC2 ~ HPC_Dissector_DissectionDate_Collapsed + RibosomePerc + 
#   Percent.Intergenic + RIN + Percent.Coding
# 
# Df Sum of Sq   RSS    AIC
# <none>                                                28236 1245.6
# - HPC_Dissector_DissectionDate_Collapsed 10    7126.8 35363 1245.7
# - Percent.Coding                          1     671.1 28907 1245.8
# + STGT_experience                         1     413.0 27823 1247.5
# + Percent.UTR                             1     390.9 27845 1247.6
# + Percent.Intronic                        1     359.1 27877 1247.9
# + LibrarySize                             1      98.5 28137 1250.2
# + MEDIAN_CV_COVERAGE                      1      92.3 28144 1250.3
# + DV200                                   1      52.4 28183 1250.6
# + Sex                                     1       6.3 28230 1251.0
# + RNAconc                                 1       5.2 28231 1251.0
# - RIN                                     1    1612.0 29848 1253.7
# + SequencingBatch                         2     198.0 28038 1254.8
# - Percent.Intergenic                      1    2776.3 31012 1263.0
# - RibosomePerc                            1    5318.1 33554 1282.3
# 
# Call:
#   lm(formula = PC2 ~ HPC_Dissector_DissectionDate_Collapsed + RibosomePerc + 
#        Percent.Intergenic + RIN + Percent.Coding, data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)   HPC_Dissector_DissectionDate_CollapsedBatch1  HPC_Dissector_DissectionDate_CollapsedBatch10  HPC_Dissector_DissectionDate_CollapsedBatch11   HPC_Dissector_DissectionDate_CollapsedBatch2  
# 91.9032                                       -19.2738                                        -6.1934                                        -0.5835                                         2.1217  
# HPC_Dissector_DissectionDate_CollapsedBatch4   HPC_Dissector_DissectionDate_CollapsedBatch5   HPC_Dissector_DissectionDate_CollapsedBatch6   HPC_Dissector_DissectionDate_CollapsedBatch7   HPC_Dissector_DissectionDate_CollapsedBatch8  
# -9.6877                                        -6.3792                                       -15.4502                                       -10.9855                                         1.2397  
# HPC_Dissector_DissectionDate_CollapsedBatch9                                   RibosomePerc                             Percent.Intergenic                                            RIN                                 Percent.Coding  
# -5.7494                                     -2277.7266                                      -673.9217                                        14.3422                                      -208.9905  
# 
# 


summary.lm(lm(PC2 ~ HPC_Dissector_DissectionDate_Collapsed + RibosomePerc + Percent.Intergenic + RIN + Percent.Coding, data = F2_PCA_wMetaData))
# 
# Call:
#   lm(formula = PC2 ~ HPC_Dissector_DissectionDate_Collapsed + RibosomePerc + 
#        Percent.Intergenic + RIN + Percent.Coding, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -28.5216  -7.7660  -0.0287   7.4783  29.0208 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                      91.9032    68.3504   1.345 0.180082    
# HPC_Dissector_DissectionDate_CollapsedBatch1    -19.2738     3.6316  -5.307 2.62e-07 ***
#   HPC_Dissector_DissectionDate_CollapsedBatch10    -6.1934     3.0506  -2.030 0.043484 *  
#   HPC_Dissector_DissectionDate_CollapsedBatch11    -0.5835     3.7060  -0.157 0.875037    
# HPC_Dissector_DissectionDate_CollapsedBatch2      2.1217     3.9663   0.535 0.593215    
# HPC_Dissector_DissectionDate_CollapsedBatch4     -9.6877     4.4295  -2.187 0.029746 *  
#   HPC_Dissector_DissectionDate_CollapsedBatch5     -6.3792     4.3314  -1.473 0.142171    
# HPC_Dissector_DissectionDate_CollapsedBatch6    -15.4502     3.9373  -3.924 0.000115 ***
#   HPC_Dissector_DissectionDate_CollapsedBatch7    -10.9855     3.5587  -3.087 0.002270 ** 
#   HPC_Dissector_DissectionDate_CollapsedBatch8      1.2397     3.5305   0.351 0.725809    
# HPC_Dissector_DissectionDate_CollapsedBatch9     -5.7494     2.9359  -1.958 0.051406 .  
# RibosomePerc                                  -2277.7266   346.0680  -6.582 3.11e-10 ***
#   Percent.Intergenic                             -673.9217   141.7129  -4.756 3.49e-06 ***
#   RIN                                              14.3422     3.9579   3.624 0.000357 ***
#   Percent.Coding                                 -208.9905    89.3879  -2.338 0.020245 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 11.08 on 230 degrees of freedom
# Multiple R-squared:  0.3633,	Adjusted R-squared:  0.3245 
# F-statistic: 9.372 on 14 and 230 DF,  p-value: 2.772e-16



step(object=lm(PC3~1, data=F2_PCA_wMetaData), scope=(~Sex+RibosomePerc+RNAconc+RIN+LibrarySize+SequencingBatch+Dissector_AsNumeric+HPC_Dissection_Date_Collapsed+HPC_Dissection_Day+STGT_experience+DV200+Percent.Intronic+Percent.Intergenic+Percent.Coding+Percent.UTR+MEDIAN_CV_COVERAGE), scale = 0, direction = "both", trace = 1, keep = NULL, steps = 1000, k=log(245))
# Start:  AIC=1199.78
# PC3 ~ 1
# 
# Df Sum of Sq   RSS    AIC
# + Sex                            1   27629.6  4447  721.2
# + RIN                            1    1759.5 30317 1191.5
# + SequencingBatch                2    2366.2 29710 1192.0
# + Percent.Intronic               1    1329.6 30747 1194.9
# + Percent.Coding                 1    1156.9 30920 1196.3
# + RNAconc                        1    1016.2 31061 1197.4
# + RibosomePerc                   1     982.4 31094 1197.7
# <none>                                       32077 1199.8
# + STGT_experience                1     458.2 31619 1201.8
# + DV200                          1     280.4 31796 1203.1
# + LibrarySize                    1     157.1 31920 1204.1
# + Percent.UTR                    1     151.9 31925 1204.1
# + Percent.Intergenic             1     150.3 31926 1204.1
# + Dissector_AsNumeric            1      45.3 32031 1204.9
# + HPC_Dissection_Day             1      30.2 32046 1205.1
# + MEDIAN_CV_COVERAGE             1       0.1 32077 1205.3
# + HPC_Dissection_Date_Collapsed  5    2129.8 29947 1210.5
# 
# Step:  AIC=721.2
# PC3 ~ Sex
# 
# Df Sum of Sq   RSS     AIC
# + RIN                            1    1025.0  3422  662.51
# + SequencingBatch                2     904.9  3542  676.47
# + HPC_Dissection_Date_Collapsed  5    1007.3  3440  685.78
# + Percent.Intronic               1     513.5  3934  696.64
# + RNAconc                        1     259.2  4188  711.99
# + Percent.Coding                 1     234.6  4213  713.42
# + RibosomePerc                   1      99.7  4347  721.15
# <none>                                        4447  721.20
# + LibrarySize                    1      69.8  4377  722.83
# + Percent.Intergenic             1      59.3  4388  723.41
# + STGT_experience                1      49.7  4398  723.95
# + DV200                          1      24.5  4423  725.35
# + MEDIAN_CV_COVERAGE             1       6.3  4441  726.35
# + Percent.UTR                    1       0.8  4446  726.66
# + Dissector_AsNumeric            1       0.2  4447  726.69
# + HPC_Dissection_Day             1       0.2  4447  726.69
# - Sex                            1   27629.6 32077 1199.78
# 
# Step:  AIC=662.51
# PC3 ~ Sex + RIN
# 
# Df Sum of Sq     RSS     AIC
# + SequencingBatch                2     901.3  2520.8  598.62
# + HPC_Dissection_Date_Collapsed  5     883.9  2538.2  616.81
# + DV200                          1     641.1  2781.0  617.19
# + RNAconc                        1     401.2  3020.9  637.46
# + Percent.Intronic               1     155.7  3266.4  656.60
# + RibosomePerc                   1      96.6  3325.5  660.99
# + STGT_experience                1      89.3  3332.9  661.54
# <none>                                        3422.1  662.51
# + LibrarySize                    1      75.7  3346.4  662.53
# + Percent.Coding                 1      74.6  3347.6  662.61
# + MEDIAN_CV_COVERAGE             1      33.3  3388.8  665.61
# + Percent.UTR                    1       9.9  3412.3  667.30
# + Dissector_AsNumeric            1       7.2  3414.9  667.49
# + Percent.Intergenic             1       3.5  3418.6  667.76
# + HPC_Dissection_Day             1       2.2  3420.0  667.86
# - RIN                            1    1025.0  4447.2  721.20
# - Sex                            1   26895.1 30317.2 1191.46
# 
# Step:  AIC=598.62
# PC3 ~ Sex + RIN + SequencingBatch
# 
# Df Sum of Sq     RSS     AIC
# + Percent.Intronic               1     394.1  2126.7  562.47
# + DV200                          1     238.0  2282.8  579.82
# + Percent.UTR                    1     142.2  2378.6  589.89
# + Percent.Intergenic             1     128.2  2392.6  591.33
# + HPC_Dissection_Day             1     127.0  2393.8  591.45
# + Percent.Coding                 1      57.7  2463.1  598.45
# <none>                                        2520.8  598.62
# + STGT_experience                1      28.6  2492.1  601.32
# + RibosomePerc                   1      25.6  2495.2  601.62
# + RNAconc                        1      11.1  2509.7  603.04
# + MEDIAN_CV_COVERAGE             1       2.8  2518.0  603.85
# + Dissector_AsNumeric            1       1.8  2519.0  603.94
# + LibrarySize                    1       0.1  2520.7  604.11
# + HPC_Dissection_Date_Collapsed  5     144.9  2375.9  611.62
# - SequencingBatch                2     901.3  3422.1  662.51
# - RIN                            1    1021.5  3542.3  676.47
# - Sex                            1   25522.9 28043.7 1183.37
# 
# Step:  AIC=562.47
# PC3 ~ Sex + RIN + SequencingBatch + Percent.Intronic
# 
# Df Sum of Sq     RSS     AIC
# + DV200                          1     132.5  1994.3  552.22
# + Percent.UTR                    1      69.1  2057.6  559.88
# + MEDIAN_CV_COVERAGE             1      66.2  2060.6  560.23
# + RNAconc                        1      62.3  2064.4  560.69
# + HPC_Dissection_Day             1      54.9  2071.8  561.56
# <none>                                        2126.7  562.47
# + STGT_experience                1      22.0  2104.7  565.42
# + Percent.Coding                 1       7.8  2118.9  567.07
# + RibosomePerc                   1       6.3  2120.4  567.25
# + Dissector_AsNumeric            1       3.0  2123.8  567.63
# + Percent.Intergenic             1       1.5  2125.2  567.80
# + LibrarySize                    1       1.2  2125.5  567.84
# + HPC_Dissection_Date_Collapsed  5     119.0  2007.8  575.88
# - Percent.Intronic               1     394.1  2520.8  598.62
# - RIN                            1     543.5  2670.2  612.73
# - SequencingBatch                2    1139.7  3266.4  656.60
# - Sex                            1   24765.0 26891.8 1178.59
# 
# Step:  AIC=552.22
# PC3 ~ Sex + RIN + SequencingBatch + Percent.Intronic + DV200
# 
# Df Sum of Sq     RSS     AIC
# + Percent.UTR                    1     115.6  1878.7  543.10
# + RNAconc                        1      54.1  1940.2  550.99
# + HPC_Dissection_Day             1      46.5  1947.8  551.94
# + MEDIAN_CV_COVERAGE             1      44.8  1949.5  552.16
# <none>                                        1994.3  552.22
# + Percent.Coding                 1      13.5  1980.8  556.05
# + STGT_experience                1       4.5  1989.8  557.17
# + Percent.Intergenic             1       3.9  1990.3  557.24
# + LibrarySize                    1       3.1  1991.2  557.34
# + RibosomePerc                   1       1.2  1993.1  557.57
# + Dissector_AsNumeric            1       0.7  1993.6  557.64
# - DV200                          1     132.5  2126.7  562.47
# + HPC_Dissection_Date_Collapsed  5     101.9  1892.4  566.88
# - Percent.Intronic               1     288.5  2282.8  579.82
# - SequencingBatch                2     684.8  2679.1  613.54
# - RIN                            1     661.5  2655.8  616.90
# - Sex                            1   24361.4 26355.7 1179.16
# 
# Step:  AIC=543.1
# PC3 ~ Sex + RIN + SequencingBatch + Percent.Intronic + DV200 + 
#   Percent.UTR
# 
# Df Sum of Sq     RSS     AIC
# <none>                                        1878.7  543.10
# + HPC_Dissection_Day             1      35.8  1842.9  543.89
# + RNAconc                        1      32.1  1846.7  544.38
# + MEDIAN_CV_COVERAGE             1      20.1  1858.6  545.96
# + Percent.Coding                 1      15.3  1863.4  546.60
# + LibrarySize                    1       9.7  1869.0  547.33
# + STGT_experience                1       9.3  1869.4  547.39
# + Percent.Intergenic             1       7.4  1871.3  547.63
# + RibosomePerc                   1       2.3  1876.4  548.30
# + Dissector_AsNumeric            1       0.1  1878.6  548.58
# - Percent.UTR                    1     115.6  1994.3  552.22
# + HPC_Dissection_Date_Collapsed  5      86.8  1791.9  559.01
# - DV200                          1     178.9  2057.6  559.88
# - Percent.Intronic               1     198.5  2077.2  562.20
# - RIN                            1     711.3  2590.1  616.26
# - SequencingBatch                2     798.9  2677.6  618.91
# - Sex                            1   24437.6 26316.3 1184.30
# 
# Call:
#   lm(formula = PC3 ~ Sex + RIN + SequencingBatch + Percent.Intronic + 
#        DV200 + Percent.UTR, data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)           SexMale               RIN  SequencingBatch2  SequencingBatch3  Percent.Intronic             DV200       Percent.UTR  
# -348.144            20.287           -10.670             2.887            -2.737          -307.748             4.283            89.878  

summary.lm(lm(PC3 ~ Sex + RIN + SequencingBatch + Percent.Intronic +  DV200 + Percent.UTR, data = F2_PCA_wMetaData))
# Call:
#   lm(formula = PC3 ~ Sex + RIN + SequencingBatch + Percent.Intronic + 
#        DV200 + Percent.UTR, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -7.7939 -1.8961 -0.1209  1.6938  8.4340 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      -348.1437    84.2519  -4.132 4.99e-05 ***
#   SexMale            20.2869     0.3654  55.523  < 2e-16 ***
#   RIN               -10.6705     1.1264  -9.473  < 2e-16 ***
#   SequencingBatch2    2.8865     0.4540   6.358 1.04e-09 ***
#   SequencingBatch3   -2.7371     0.4985  -5.491 1.03e-07 ***
#   Percent.Intronic -307.7479    61.5016  -5.004 1.10e-06 ***
#   DV200               4.2832     0.9016   4.751 3.51e-06 ***
#   Percent.UTR        89.8784    23.5401   3.818 0.000172 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 2.816 on 237 degrees of freedom
# Multiple R-squared:  0.9414,	Adjusted R-squared:  0.9397 
# F-statistic: 544.2 on 7 and 237 DF,  p-value: < 2.2e-16

step(object=lm(PC4~1, data=F2_PCA_wMetaData), scope=(~Sex+RibosomePerc+RNAconc+RIN+LibrarySize+SequencingBatch+Dissector_AsNumeric+HPC_Dissection_Date_Collapsed+HPC_Dissection_Day+STGT_experience+DV200+Percent.Intronic+Percent.Intergenic+Percent.Coding+Percent.UTR+MEDIAN_CV_COVERAGE), scale = 0, direction = "both", trace = 1, keep = NULL, steps = 1000, k=log(245))
# Start:  AIC=1100.72
# PC4 ~ 1
# 
# Df Sum of Sq   RSS    AIC
# + SequencingBatch                2    4765.1 16643 1050.0
# + HPC_Dissection_Date_Collapsed  5    5410.1 15998 1056.8
# + RIN                            1    2782.5 18626 1072.1
# + RNAconc                        1    2767.0 18641 1072.3
# + Sex                            1    2377.1 19031 1077.4
# + RibosomePerc                   1    1367.2 20041 1090.0
# + MEDIAN_CV_COVERAGE             1     851.3 20557 1096.3
# + STGT_experience                1     478.9 20929 1100.7
# <none>                                       21408 1100.7
# + Percent.Intronic               1     441.3 20967 1101.1
# + LibrarySize                    1     423.8 20984 1101.3
# + Percent.UTR                    1     145.5 21263 1104.5
# + Dissector_AsNumeric            1     111.2 21297 1104.9
# + DV200                          1      92.5 21316 1105.2
# + Percent.Intergenic             1      30.5 21378 1105.9
# + HPC_Dissection_Day             1       7.7 21400 1106.1
# + Percent.Coding                 1       3.3 21405 1106.2
# 
# Step:  AIC=1050.03
# PC4 ~ SequencingBatch
# 
# Df Sum of Sq   RSS    AIC
# + Sex                            1    3177.8 13465 1003.6
# + RIN                            1    2809.0 13834 1010.2
# + Percent.UTR                    1    2436.9 14206 1016.8
# + Percent.Intronic               1    1307.5 15336 1035.5
# + MEDIAN_CV_COVERAGE             1     420.4 16223 1049.3
# <none>                                       16643 1050.0
# + RNAconc                        1     242.9 16400 1051.9
# + Percent.Intergenic             1     207.4 16436 1052.5
# + DV200                          1     184.2 16459 1052.8
# + STGT_experience                1     146.6 16496 1053.4
# + HPC_Dissection_Day             1     144.9 16498 1053.4
# + RibosomePerc                   1      20.1 16623 1055.2
# + Percent.Coding                 1       2.8 16640 1055.5
# + LibrarySize                    1       2.1 16641 1055.5
# + Dissector_AsNumeric            1       1.0 16642 1055.5
# + HPC_Dissection_Date_Collapsed  5     920.8 15722 1063.6
# - SequencingBatch                2    4765.1 21408 1100.7
# 
# Step:  AIC=1003.62
# PC4 ~ SequencingBatch + Sex
# 
# Df Sum of Sq   RSS     AIC
# + RIN                            1    3157.5 10308  943.66
# + Percent.UTR                    1    2324.9 11140  962.69
# + Percent.Intronic               1    1775.3 11690  974.49
# + Percent.Intergenic             1     344.8 13120 1002.77
# + MEDIAN_CV_COVERAGE             1     338.1 13127 1002.89
# + RNAconc                        1     336.0 13129 1002.93
# <none>                                       13465 1003.62
# + HPC_Dissection_Day             1     231.6 13234 1004.87
# + DV200                          1     119.0 13346 1006.95
# + STGT_experience                1     111.8 13354 1007.08
# + RibosomePerc                   1      84.8 13380 1007.58
# + Percent.Coding                 1      24.6 13441 1008.68
# + Dissector_AsNumeric            1       7.5 13458 1008.99
# + LibrarySize                    1       0.1 13465 1009.12
# + HPC_Dissection_Date_Collapsed  5     965.9 12499 1012.89
# - Sex                            1    3177.8 16643 1050.03
# - SequencingBatch                2    5565.8 19031 1077.38
# 
# Step:  AIC=943.66
# PC4 ~ SequencingBatch + Sex + RIN
# 
# Df Sum of Sq     RSS     AIC
# + Percent.UTR                    1    1912.5  8395.2  898.88
# + RNAconc                        1     673.1  9634.6  932.61
# + MEDIAN_CV_COVERAGE             1     651.6  9656.1  933.16
# + Percent.Intronic               1     600.1  9707.6  934.46
# + DV200                          1     536.3  9771.4  936.07
# <none>                                       10307.7  943.66
# + HPC_Dissection_Day             1     135.1 10172.6  945.93
# + RibosomePerc                   1      93.1 10214.7  946.94
# + Percent.Coding                 1      71.6 10236.1  947.45
# + Percent.Intergenic             1      39.4 10268.4  948.22
# + Dissector_AsNumeric            1      34.3 10273.4  948.34
# + STGT_experience                1      13.3 10294.4  948.84
# + LibrarySize                    1       7.1 10300.7  948.99
# + HPC_Dissection_Date_Collapsed  5     716.2  9591.6  953.52
# - RIN                            1    3157.5 13465.2 1003.62
# - Sex                            1    3526.3 13834.0 1010.24
# - SequencingBatch                2    5611.9 15919.6 1039.15
# 
# Step:  AIC=898.88
# PC4 ~ SequencingBatch + Sex + RIN + Percent.UTR
# 
# Df Sum of Sq     RSS     AIC
# + DV200                          1     892.9  7502.3  876.83
# + RNAconc                        1     463.4  7931.8  890.47
# + MEDIAN_CV_COVERAGE             1     427.9  7967.3  891.56
# + Percent.Intronic               1     272.4  8122.8  896.30
# + RibosomePerc                   1     232.8  8162.4  897.49
# <none>                                        8395.2  898.88
# + Percent.Coding                 1     139.0  8256.2  900.29
# + STGT_experience                1      78.6  8316.6  902.07
# + HPC_Dissection_Day             1      50.9  8344.4  902.89
# + Dissector_AsNumeric            1      35.8  8359.4  903.33
# + Percent.Intergenic             1      15.8  8379.4  903.92
# + LibrarySize                    1       6.9  8388.3  904.18
# + HPC_Dissection_Date_Collapsed  5     521.5  7873.7  910.67
# - Percent.UTR                    1    1912.5 10307.7  943.66
# - RIN                            1    2745.1 11140.3  962.69
# - Sex                            1    3394.8 11790.0  976.57
# - SequencingBatch                2    7501.1 15896.3 1044.29
# 
# Step:  AIC=876.83
# PC4 ~ SequencingBatch + Sex + RIN + Percent.UTR + DV200
# 
# Df Sum of Sq     RSS     AIC
# + RNAconc                        1     444.8  7057.5  867.36
# + MEDIAN_CV_COVERAGE             1     342.3  7160.0  870.89
# + RibosomePerc                   1     219.6  7282.8  875.05
# <none>                                        7502.3  876.83
# + Percent.Intronic               1      89.6  7412.7  879.38
# + Percent.Coding                 1      68.5  7433.8  880.08
# + LibrarySize                    1      28.7  7473.6  881.39
# + Dissector_AsNumeric            1      22.0  7480.3  881.61
# + HPC_Dissection_Day             1      15.1  7487.3  881.84
# + STGT_experience                1       5.8  7496.5  882.14
# + Percent.Intergenic             1       4.4  7498.0  882.19
# + HPC_Dissection_Date_Collapsed  5     443.7  7058.6  889.40
# - DV200                          1     892.9  8395.2  898.88
# - Percent.UTR                    1    2269.1  9771.4  936.07
# - RIN                            1    3629.7 11132.0  968.01
# - Sex                            1    3682.6 11184.9  969.17
# - SequencingBatch                2    5595.1 13097.4 1002.34
# 
# Step:  AIC=867.36
# PC4 ~ SequencingBatch + Sex + RIN + Percent.UTR + DV200 + RNAconc
# 
# Df Sum of Sq     RSS    AIC
# + Percent.Coding                 1     720.0  6337.5 846.49
# + Percent.Intergenic             1     571.0  6486.5 852.19
# + Percent.Intronic               1     237.5  6820.0 864.47
# <none>                                        7057.5 867.36
# + Dissector_AsNumeric            1     108.5  6949.0 869.06
# + MEDIAN_CV_COVERAGE             1      16.2  7041.3 872.29
# + LibrarySize                    1      15.3  7042.2 872.32
# + STGT_experience                1       1.9  7055.6 872.79
# + RibosomePerc                   1       1.5  7056.1 872.81
# + HPC_Dissection_Day             1       1.1  7056.4 872.82
# - RNAconc                        1     444.8  7502.3 876.83
# + HPC_Dissection_Date_Collapsed  5     371.8  6685.7 881.60
# - DV200                          1     874.3  7931.8 890.47
# - Percent.UTR                    1    2039.9  9097.4 924.06
# - SequencingBatch                2    2661.3  9718.9 934.75
# - Sex                            1    3818.5 10876.0 967.81
# - RIN                            1    3870.6 10928.1 968.98
# 
# Step:  AIC=846.49
# PC4 ~ SequencingBatch + Sex + RIN + Percent.UTR + DV200 + RNAconc + 
#   Percent.Coding
# 
# Df Sum of Sq     RSS    AIC
# + MEDIAN_CV_COVERAGE             1     312.5  6025.1 839.61
# <none>                                        6337.5 846.49
# + RibosomePerc                   1      92.6  6244.9 848.39
# + HPC_Dissection_Day             1      52.6  6284.9 849.95
# + Percent.Intergenic             1      42.7  6294.9 850.34
# + Dissector_AsNumeric            1      42.5  6295.1 850.35
# + Percent.Intronic               1       8.7  6328.9 851.66
# + STGT_experience                1       8.0  6329.5 851.68
# + LibrarySize                    1       0.0  6337.5 851.99
# + HPC_Dissection_Date_Collapsed  5     291.6  6045.9 862.46
# - DV200                          1     617.3  6954.8 863.76
# - Percent.Coding                 1     720.0  7057.5 867.36
# - RNAconc                        1    1096.3  7433.8 880.08
# - SequencingBatch                2    2321.6  8659.1 911.96
# - RIN                            1    2593.5  8931.0 925.04
# - Percent.UTR                    1    2746.7  9084.3 929.21
# - Sex                            1    4273.9 10611.4 967.28
# 
# Step:  AIC=839.61
# PC4 ~ SequencingBatch + Sex + RIN + Percent.UTR + DV200 + RNAconc + 
#   Percent.Coding + MEDIAN_CV_COVERAGE
# 
# Df Sum of Sq     RSS    AIC
# <none>                                        6025.1 839.61
# + Dissector_AsNumeric            1      62.9  5962.1 842.54
# + HPC_Dissection_Day             1      46.8  5978.3 843.20
# + Percent.Intronic               1      11.8  6013.3 844.63
# + STGT_experience                1       7.6  6017.4 844.80
# + RibosomePerc                   1       2.4  6022.7 845.01
# + LibrarySize                    1       1.4  6023.6 845.05
# + Percent.Intergenic             1       0.1  6025.0 845.10
# - RNAconc                        1     287.6  6312.6 845.53
# - MEDIAN_CV_COVERAGE             1     312.5  6337.5 846.49
# - DV200                          1     449.0  6474.1 851.72
# + HPC_Dissection_Date_Collapsed  5     247.8  5777.2 856.82
# - Percent.Coding                 1    1016.3  7041.3 872.29
# - RIN                            1    2236.7  8261.7 911.45
# - SequencingBatch                2    2608.2  8633.2 916.73
# - Percent.UTR                    1    2968.3  8993.3 932.24
# - Sex                            1    4151.5 10176.6 962.53
# 
# Call:
#   lm(formula = PC4 ~ SequencingBatch + Sex + RIN + Percent.UTR + 
#        DV200 + RNAconc + Percent.Coding + MEDIAN_CV_COVERAGE, data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)    SequencingBatch2    SequencingBatch3             SexMale                 RIN         Percent.UTR               DV200             RNAconc      Percent.Coding  MEDIAN_CV_COVERAGE  
# -706.75442             5.60915            -6.84842            -8.47777           -18.56592           505.29887             6.78128             0.02386           200.15321          -113.22619  

summary.lm(lm(PC4 ~ SequencingBatch + Sex + RIN + Percent.UTR + DV200 + RNAconc + mRNAperc + SeqID_Numeric + MEDIAN_CV_COVERAGE, data = F2_PCA_wMetaData))

# 
# Call:
#   lm(formula = PC4 ~ SequencingBatch + Sex + RIN + Percent.UTR + 
#        DV200 + RNAconc + mRNAperc + SeqID_Numeric + MEDIAN_CV_COVERAGE, 
#      data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -15.8311  -3.3550  -0.3042   2.8564  13.9218 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         2.228e+04  6.810e+03   3.272 0.001228 ** 
#   SequencingBatch2    3.112e+01  7.641e+00   4.073 6.36e-05 ***
#   SequencingBatch3    3.531e+01  1.253e+01   2.818 0.005251 ** 
#   SexMale            -8.696e+00  6.535e-01 -13.307  < 2e-16 ***
#   RIN                -1.852e+01  1.924e+00  -9.626  < 2e-16 ***
#   Percent.UTR         3.218e+02  4.200e+01   7.662 4.84e-13 ***
#   DV200               7.178e+00  1.576e+00   4.555 8.44e-06 ***
#   RNAconc             3.580e-02  7.731e-03   4.631 6.04e-06 ***
#   mRNAperc            2.374e+02  3.254e+01   7.295 4.60e-12 ***
#   SeqID_Numeric      -4.908e-02  1.453e-02  -3.379 0.000852 ***
#   MEDIAN_CV_COVERAGE -1.048e+02  3.157e+01  -3.321 0.001041 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 4.942 on 234 degrees of freedom
# Multiple R-squared:  0.7331,	Adjusted R-squared:  0.7216 
# F-statistic: 64.26 on 10 and 234 DF,  p-value: < 2.2e-16

step(object=lm(PC5~1, data=F2_PCA_wMetaData), scope=(~Sex+RibosomePerc+RNAconc+RIN+LibrarySize+SequencingBatch+Dissector_AsNumeric+HPC_Dissection_Date_Collapsed+HPC_Dissection_Day+STGT_experience+DV200+Percent.Intronic+Percent.Intergenic+Percent.Coding+Percent.UTR+MEDIAN_CV_COVERAGE), scale = 0, direction = "both", trace = 1, keep = NULL, steps = 1000, k=log(245))

# Start:  AIC=939.1
# PC5 ~ 1
# 
# Df Sum of Sq     RSS    AIC
# + Percent.UTR                    1    5333.5  5735.2 783.52
# + Percent.Coding                 1    3512.1  7556.7 851.09
# + SequencingBatch                2     713.8 10355.0 933.77
# + HPC_Dissection_Day             1     379.9 10688.9 936.05
# + MEDIAN_CV_COVERAGE             1     349.7 10719.0 936.74
# + Percent.Intergenic             1     294.3 10774.5 938.00
# <none>                                       11068.7 939.10
# + DV200                          1     204.8 10863.9 940.03
# + RNAconc                        1     204.7 10864.0 940.03
# + Dissector_AsNumeric            1     187.2 10881.6 940.43
# + RIN                            1     170.5 10898.2 940.80
# + Percent.Intronic               1      69.7 10999.0 943.06
# + Sex                            1      49.9 11018.8 943.50
# + RibosomePerc                   1      31.8 11037.0 943.90
# + LibrarySize                    1       3.1 11065.6 944.54
# + STGT_experience                1       0.6 11068.1 944.59
# + HPC_Dissection_Date_Collapsed  5     489.5 10579.2 955.53
# 
# Step:  AIC=783.52
# PC5 ~ Percent.UTR
# 
# Df Sum of Sq     RSS    AIC
# + Percent.Coding                 1    1208.3  4526.9 731.06
# + Percent.Intergenic             1    1019.6  4715.6 741.06
# + RNAconc                        1     855.3  4879.9 749.45
# + Percent.Intronic               1     810.1  4925.1 751.71
# + SequencingBatch                2     511.0  5224.2 771.66
# + RibosomePerc                   1     365.7  5369.4 772.87
# + Dissector_AsNumeric            1     337.9  5397.3 774.14
# + RIN                            1     320.2  5414.9 774.94
# + HPC_Dissection_Date_Collapsed  5     779.0  4956.2 775.26
# + MEDIAN_CV_COVERAGE             1     301.7  5433.5 775.78
# + Sex                            1     146.9  5588.3 782.66
# + STGT_experience                1     140.1  5595.1 782.96
# <none>                                        5735.2 783.52
# + HPC_Dissection_Day             1      77.3  5657.9 785.69
# + DV200                          1      42.2  5692.9 787.21
# + LibrarySize                    1      26.4  5708.8 787.89
# - Percent.UTR                    1    5333.5 11068.7 939.10
# 
# Step:  AIC=731.06
# PC5 ~ Percent.UTR + Percent.Coding
# 
# Df Sum of Sq    RSS    AIC
# + SequencingBatch                2    589.77 3937.2 707.86
# + HPC_Dissection_Date_Collapsed  5    713.70 3813.2 716.53
# + HPC_Dissection_Day             1    318.92 4208.0 718.66
# + Sex                            1    239.49 4287.4 723.24
# + RNAconc                        1    145.09 4381.8 728.58
# <none>                                       4526.9 731.06
# + MEDIAN_CV_COVERAGE             1     89.58 4437.3 731.66
# + RIN                            1     82.04 4444.9 732.08
# + STGT_experience                1     63.17 4463.8 733.11
# + Percent.Intronic               1     53.30 4473.6 733.66
# + Percent.Intergenic             1     38.45 4488.5 734.47
# + DV200                          1     19.91 4507.0 735.48
# + Dissector_AsNumeric            1     16.60 4510.3 735.66
# + RibosomePerc                   1      7.97 4519.0 736.13
# + LibrarySize                    1      0.18 4526.7 736.55
# - Percent.Coding                 1   1208.25 5735.2 783.52
# - Percent.UTR                    1   3029.74 7556.7 851.09
# 
# Step:  AIC=707.86
# PC5 ~ Percent.UTR + Percent.Coding + SequencingBatch
# 
# Df Sum of Sq    RSS    AIC
# + Sex                            1    253.61 3683.5 697.05
# + STGT_experience                1    116.26 3820.9 706.02
# <none>                                       3937.2 707.86
# + MEDIAN_CV_COVERAGE             1     70.84 3866.3 708.91
# + RibosomePerc                   1     48.34 3888.8 710.33
# + RIN                            1     42.41 3894.7 710.71
# + Percent.Intergenic             1     26.74 3910.4 711.69
# + LibrarySize                    1      3.95 3933.2 713.12
# + Dissector_AsNumeric            1      2.95 3934.2 713.18
# + DV200                          1      2.42 3934.7 713.21
# + RNAconc                        1      1.29 3935.9 713.28
# + HPC_Dissection_Day             1      0.02 3937.1 713.36
# + Percent.Intronic               1      0.01 3937.1 713.36
# + HPC_Dissection_Date_Collapsed  5    246.39 3690.8 719.53
# - SequencingBatch                2    589.77 4526.9 731.06
# - Percent.Coding                 1   1287.07 5224.2 771.66
# - Percent.UTR                    1   2613.39 6550.5 827.09
# 
# Step:  AIC=697.05
# PC5 ~ Percent.UTR + Percent.Coding + SequencingBatch + Sex
# 
# Df Sum of Sq    RSS    AIC
# + STGT_experience                1    106.56 3577.0 695.36
# + RibosomePerc                   1    104.39 3579.2 695.51
# + MEDIAN_CV_COVERAGE             1     88.53 3595.0 696.59
# <none>                                       3683.5 697.05
# + Percent.Intergenic             1     60.75 3622.8 698.48
# + RIN                            1     47.28 3636.3 699.39
# + Dissector_AsNumeric            1     10.48 3673.1 701.85
# + DV200                          1      5.66 3677.9 702.17
# + LibrarySize                    1      4.48 3679.1 702.25
# + RNAconc                        1      1.80 3681.7 702.43
# + HPC_Dissection_Day             1      0.58 3683.0 702.51
# + Percent.Intronic               1      0.23 3683.3 702.53
# + HPC_Dissection_Date_Collapsed  5    271.49 3412.1 705.80
# - Sex                            1    253.61 3937.2 707.86
# - SequencingBatch                2    603.89 4287.4 723.24
# - Percent.Coding                 1   1407.36 5090.9 770.82
# - Percent.UTR                    1   2563.74 6247.3 820.97
# 
# Step:  AIC=695.36
# PC5 ~ Percent.UTR + Percent.Coding + SequencingBatch + Sex + 
#   STGT_experience
# 
# Df Sum of Sq    RSS    AIC
# + RibosomePerc                   1    101.01 3476.0 693.84
# <none>                                       3577.0 695.36
# + MEDIAN_CV_COVERAGE             1     76.42 3500.6 695.57
# - STGT_experience                1    106.56 3683.5 697.05
# + Percent.Intergenic             1     49.10 3527.9 697.47
# + RIN                            1     30.03 3547.0 698.79
# + HPC_Dissection_Day             1      5.90 3571.1 700.45
# + LibrarySize                    1      3.54 3573.4 700.62
# + DV200                          1      1.55 3575.4 700.75
# + RNAconc                        1      0.91 3576.1 700.80
# + Percent.Intronic               1      0.91 3576.1 700.80
# + Dissector_AsNumeric            1      0.88 3576.1 700.80
# - Sex                            1    243.91 3820.9 706.02
# + HPC_Dissection_Date_Collapsed  5    216.09 3360.9 707.60
# - SequencingBatch                2    666.93 4243.9 726.24
# - Percent.Coding                 1   1409.19 4986.2 771.23
# - Percent.UTR                    1   2437.76 6014.7 817.18
# 
# Step:  AIC=693.84
# PC5 ~ Percent.UTR + Percent.Coding + SequencingBatch + Sex + 
#   STGT_experience + RibosomePerc
# 
# Df Sum of Sq    RSS    AIC
# <none>                                       3476.0 693.84
# - RibosomePerc                   1    101.01 3577.0 695.36
# - STGT_experience                1    103.18 3579.2 695.51
# + RNAconc                        1     33.64 3442.3 696.96
# + RIN                            1      9.67 3466.3 698.66
# + MEDIAN_CV_COVERAGE             1      8.51 3467.5 698.74
# + Percent.Intronic               1      7.06 3468.9 698.84
# + HPC_Dissection_Day             1      6.97 3469.0 698.85
# + Dissector_AsNumeric            1      6.41 3469.6 698.89
# + Percent.Intergenic             1      5.54 3470.4 698.95
# + DV200                          1      2.59 3473.4 699.16
# + LibrarySize                    1      0.74 3475.2 699.29
# + HPC_Dissection_Date_Collapsed  5    218.57 3257.4 705.44
# - Sex                            1    298.11 3774.1 708.50
# - SequencingBatch                2    711.23 4187.2 728.45
# - Percent.Coding                 1   1275.63 4751.6 764.93
# - Percent.UTR                    1   1650.90 5126.9 783.55
# 
# Call:
#   lm(formula = PC5 ~ Percent.UTR + Percent.Coding + SequencingBatch + 
#        Sex + STGT_experience + RibosomePerc, data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)          Percent.UTR       Percent.Coding     SequencingBatch2     SequencingBatch3              SexMale  STGT_experienceTRUE         RibosomePerc  
# -49.078              394.277             -192.017                3.783                3.222                2.282                2.117              216.108 
# 

summary.lm(lm(PC5 ~ Percent.UTR + Percent.Coding + SequencingBatch + Sex + STGT_experience + RibosomePerc, data = F2_PCA_wMetaData))

# 
# Call:
#   lm(formula = PC5 ~ Percent.UTR + Percent.Coding + SequencingBatch + 
#        Sex + STGT_experience + RibosomePerc, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -8.7216 -2.7958 -0.3894  2.7500 12.3778 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          -49.0783    19.1815  -2.559  0.01113 *  
#   Percent.UTR          394.2769    37.1625  10.610  < 2e-16 ***
#   Percent.Coding      -192.0170    20.5893  -9.326  < 2e-16 ***
#   SequencingBatch2       3.7830     0.5974   6.332 1.20e-09 ***
#   SequencingBatch3       3.2215     0.8104   3.975 9.34e-05 ***
#   SexMale                2.2817     0.5061   4.508 1.03e-05 ***
#   STGT_experienceTRUE    2.1166     0.7980   2.652  0.00853 ** 
#   RibosomePerc         216.1080    82.3497   2.624  0.00925 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 3.83 on 237 degrees of freedom
# Multiple R-squared:  0.686,	Adjusted R-squared:  0.6767 
# F-statistic: 73.96 on 7 and 237 DF,  p-value: < 2.2e-16

######PC5 with HPC_Dissector_DissectionDate_Collapsed
step(object=lm(PC5~1, data=F2_PCA_wMetaData), scope=(~Sex+RibosomePerc+RNAconc+RIN+LibrarySize+SequencingBatch+HPC_Dissector_DissectionDate_Collapsed+STGT_experience+DV200+Percent.Intronic+Percent.Intergenic+Percent.Coding+Percent.UTR+MEDIAN_CV_COVERAGE), scale = 0, direction = "both", trace = 1, keep = NULL, steps = 1000, k=log(245))
# Start:  AIC=939.1
# PC5 ~ 1
# 
# Df Sum of Sq     RSS    AIC
# + Percent.UTR                             1    5333.5  5735.2 783.52
# + Percent.Coding                          1    3512.1  7556.7 851.09
# + SequencingBatch                         2     713.8 10355.0 933.77
# + MEDIAN_CV_COVERAGE                      1     349.7 10719.0 936.74
# + Percent.Intergenic                      1     294.3 10774.5 938.00
# <none>                                                11068.7 939.10
# + DV200                                   1     204.8 10863.9 940.03
# + RNAconc                                 1     204.7 10864.0 940.03
# + RIN                                     1     170.5 10898.2 940.80
# + Percent.Intronic                        1      69.7 10999.0 943.06
# + Sex                                     1      49.9 11018.8 943.50
# + RibosomePerc                            1      31.8 11037.0 943.90
# + LibrarySize                             1       3.1 11065.6 944.54
# + STGT_experience                         1       0.6 11068.1 944.59
# + HPC_Dissector_DissectionDate_Collapsed 10    1140.4  9928.3 967.48
# 
# Step:  AIC=783.52
# PC5 ~ Percent.UTR
# 
# Df Sum of Sq     RSS    AIC
# + Percent.Coding                          1    1208.3  4526.9 731.06
# + Percent.Intergenic                      1    1019.6  4715.6 741.06
# + RNAconc                                 1     855.3  4879.9 749.45
# + Percent.Intronic                        1     810.1  4925.1 751.71
# + HPC_Dissector_DissectionDate_Collapsed 10    1418.3  4316.9 768.93
# + SequencingBatch                         2     511.0  5224.2 771.66
# + RibosomePerc                            1     365.7  5369.4 772.87
# + RIN                                     1     320.2  5414.9 774.94
# + MEDIAN_CV_COVERAGE                      1     301.7  5433.5 775.78
# + Sex                                     1     146.9  5588.3 782.66
# + STGT_experience                         1     140.1  5595.1 782.96
# <none>                                                 5735.2 783.52
# + DV200                                   1      42.2  5692.9 787.21
# + LibrarySize                             1      26.4  5708.8 787.89
# - Percent.UTR                             1    5333.5 11068.7 939.10
# 
# Step:  AIC=731.06
# PC5 ~ Percent.UTR + Percent.Coding
# 
# Df Sum of Sq    RSS    AIC
# + SequencingBatch                         2    589.77 3937.2 707.86
# + Sex                                     1    239.49 4287.4 723.24
# + RNAconc                                 1    145.09 4381.8 728.58
# <none>                                                4526.9 731.06
# + MEDIAN_CV_COVERAGE                      1     89.58 4437.3 731.66
# + RIN                                     1     82.04 4444.9 732.08
# + STGT_experience                         1     63.17 4463.8 733.11
# + Percent.Intronic                        1     53.30 4473.6 733.66
# + Percent.Intergenic                      1     38.45 4488.5 734.47
# + DV200                                   1     19.91 4507.0 735.48
# + RibosomePerc                            1      7.97 4519.0 736.13
# + LibrarySize                             1      0.18 4526.7 736.55
# + HPC_Dissector_DissectionDate_Collapsed 10    801.05 3725.9 738.36
# - Percent.Coding                          1   1208.25 5735.2 783.52
# - Percent.UTR                             1   3029.74 7556.7 851.09
# 
# Step:  AIC=707.86
# PC5 ~ Percent.UTR + Percent.Coding + SequencingBatch
# 
# Df Sum of Sq    RSS    AIC
# + Sex                                     1    253.61 3683.5 697.05
# + STGT_experience                         1    116.26 3820.9 706.02
# <none>                                                3937.2 707.86
# + MEDIAN_CV_COVERAGE                      1     70.84 3866.3 708.91
# + RibosomePerc                            1     48.34 3888.8 710.33
# + RIN                                     1     42.41 3894.7 710.71
# + Percent.Intergenic                      1     26.74 3910.4 711.69
# + LibrarySize                             1      3.95 3933.2 713.12
# + DV200                                   1      2.42 3934.7 713.21
# + RNAconc                                 1      1.29 3935.9 713.28
# + Percent.Intronic                        1      0.01 3937.1 713.36
# - SequencingBatch                         2    589.77 4526.9 731.06
# + HPC_Dissector_DissectionDate_Collapsed 10    367.02 3570.1 738.90
# - Percent.Coding                          1   1287.07 5224.2 771.66
# - Percent.UTR                             1   2613.39 6550.5 827.09
# 
# Step:  AIC=697.05
# PC5 ~ Percent.UTR + Percent.Coding + SequencingBatch + Sex
# 
# Df Sum of Sq    RSS    AIC
# + STGT_experience                         1    106.56 3577.0 695.36
# + RibosomePerc                            1    104.39 3579.2 695.51
# + MEDIAN_CV_COVERAGE                      1     88.53 3595.0 696.59
# <none>                                                3683.5 697.05
# + Percent.Intergenic                      1     60.75 3622.8 698.48
# + RIN                                     1     47.28 3636.3 699.39
# + DV200                                   1      5.66 3677.9 702.17
# + LibrarySize                             1      4.48 3679.1 702.25
# + RNAconc                                 1      1.80 3681.7 702.43
# + Percent.Intronic                        1      0.23 3683.3 702.53
# - Sex                                     1    253.61 3937.2 707.86
# - SequencingBatch                         2    603.89 4287.4 723.24
# + HPC_Dissector_DissectionDate_Collapsed 10    399.13 3284.4 723.96
# - Percent.Coding                          1   1407.36 5090.9 770.82
# - Percent.UTR                             1   2563.74 6247.3 820.97
# 
# Step:  AIC=695.36
# PC5 ~ Percent.UTR + Percent.Coding + SequencingBatch + Sex + 
#   STGT_experience
# 
# Df Sum of Sq    RSS    AIC
# + RibosomePerc                            1    101.01 3476.0 693.84
# <none>                                                3577.0 695.36
# + MEDIAN_CV_COVERAGE                      1     76.42 3500.6 695.57
# - STGT_experience                         1    106.56 3683.5 697.05
# + Percent.Intergenic                      1     49.10 3527.9 697.47
# + RIN                                     1     30.03 3547.0 698.79
# + LibrarySize                             1      3.54 3573.4 700.62
# + DV200                                   1      1.55 3575.4 700.75
# + RNAconc                                 1      0.91 3576.1 700.80
# + Percent.Intronic                        1      0.91 3576.1 700.80
# - Sex                                     1    243.91 3820.9 706.02
# - SequencingBatch                         2    666.93 4243.9 726.24
# + HPC_Dissector_DissectionDate_Collapsed 10    307.21 3269.8 728.37
# - Percent.Coding                          1   1409.19 4986.2 771.23
# - Percent.UTR                             1   2437.76 6014.7 817.18
# 
# Step:  AIC=693.84
# PC5 ~ Percent.UTR + Percent.Coding + SequencingBatch + Sex + 
#   STGT_experience + RibosomePerc
# 
# Df Sum of Sq    RSS    AIC
# <none>                                                3476.0 693.84
# - RibosomePerc                            1    101.01 3577.0 695.36
# - STGT_experience                         1    103.18 3579.2 695.51
# + RNAconc                                 1     33.64 3442.3 696.96
# + RIN                                     1      9.67 3466.3 698.66
# + MEDIAN_CV_COVERAGE                      1      8.51 3467.5 698.74
# + Percent.Intronic                        1      7.06 3468.9 698.84
# + Percent.Intergenic                      1      5.54 3470.4 698.95
# + DV200                                   1      2.59 3473.4 699.16
# + LibrarySize                             1      0.74 3475.2 699.29
# - Sex                                     1    298.11 3774.1 708.50
# + HPC_Dissector_DissectionDate_Collapsed 10    291.86 3184.1 727.37
# - SequencingBatch                         2    711.23 4187.2 728.45
# - Percent.Coding                          1   1275.63 4751.6 764.93
# - Percent.UTR                             1   1650.90 5126.9 783.55
# 
# Call:
#   lm(formula = PC5 ~ Percent.UTR + Percent.Coding + SequencingBatch + 
#        Sex + STGT_experience + RibosomePerc, data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)          Percent.UTR       Percent.Coding     SequencingBatch2     SequencingBatch3              SexMale  STGT_experienceTRUE         RibosomePerc  
# -49.078              394.277             -192.017                3.783                3.222                2.282                2.117              216.108  




summary.lm(lm(PC5 ~ Percent.UTR + Percent.Coding + SequencingBatch + Sex + STGT_experience + RibosomePerc, data = F2_PCA_wMetaData))
# 
# 
# Call:
#   lm(formula = PC5 ~ Percent.UTR + Percent.Coding + SequencingBatch + 
#        Sex + STGT_experience + RibosomePerc, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -8.7216 -2.7958 -0.3894  2.7500 12.3778 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          -49.0783    19.1815  -2.559  0.01113 *  
#   Percent.UTR          394.2769    37.1625  10.610  < 2e-16 ***
#   Percent.Coding      -192.0170    20.5893  -9.326  < 2e-16 ***
#   SequencingBatch2       3.7830     0.5974   6.332 1.20e-09 ***
#   SequencingBatch3       3.2215     0.8104   3.975 9.34e-05 ***
#   SexMale                2.2817     0.5061   4.508 1.03e-05 ***
#   STGT_experienceTRUE    2.1166     0.7980   2.652  0.00853 ** 
#   RibosomePerc         216.1080    82.3497   2.624  0.00925 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 3.83 on 237 degrees of freedom
# Multiple R-squared:  0.686,	Adjusted R-squared:  0.6767 
# F-statistic: 73.96 on 7 and 237 DF,  p-value: < 2.2e-16
# 


step(object=lm(PC6~1, data=F2_PCA_wMetaData), scope=(~Sex+RibosomePerc+RNAconc+RIN+LibrarySize+SequencingBatch+Dissector_AsNumeric+HPC_Dissection_Date_Collapsed+HPC_Dissection_Day+STGT_experience+DV200+Percent.Intronic+Percent.Intergenic+Percent.Coding+Percent.UTR+MEDIAN_CV_COVERAGE), scale = 0, direction = "both", trace = 1, keep = NULL, steps = 1000, k=log(245))
# Start:  AIC=893.12
# PC6 ~ 1
# 
# Df Sum of Sq    RSS    AIC
# + Dissector_AsNumeric            1   1150.99 8023.5 865.78
# + Percent.UTR                    1    846.65 8327.9 874.90
# + HPC_Dissection_Date_Collapsed  5   1498.89 7675.6 876.92
# + STGT_experience                1    734.72 8439.8 878.17
# + HPC_Dissection_Day             1    576.48 8598.1 882.72
# + SequencingBatch                2    616.03 8558.5 887.09
# + Percent.Coding                 1    342.51 8832.0 889.30
# <none>                                       9174.5 893.12
# + Percent.Intronic               1    123.86 9050.7 895.29
# + RIN                            1     86.87 9087.7 896.29
# + MEDIAN_CV_COVERAGE             1     58.80 9115.7 897.04
# + RibosomePerc                   1     54.38 9120.1 897.16
# + RNAconc                        1     51.22 9123.3 897.25
# + LibrarySize                    1     38.26 9136.3 897.60
# + Sex                            1     18.18 9156.3 898.13
# + Percent.Intergenic             1      9.81 9164.7 898.36
# + DV200                          1      5.96 9168.6 898.46
# 
# Step:  AIC=865.78
# PC6 ~ Dissector_AsNumeric
# 
# Df Sum of Sq    RSS    AIC
# + Percent.UTR                    1    981.38 7042.2 839.31
# + RibosomePerc                   1    723.60 7299.9 848.12
# + HPC_Dissection_Date_Collapsed  5   1291.61 6731.9 850.28
# + SequencingBatch                2    803.64 7219.9 850.92
# + STGT_experience                1    510.00 7513.5 855.19
# + Percent.Intronic               1    330.08 7693.5 860.99
# + HPC_Dissection_Day             1    303.94 7719.6 861.82
# + RNAconc                        1    232.00 7791.5 864.09
# + Percent.Intergenic             1    230.10 7793.4 864.15
# <none>                                       8023.5 865.78
# + RIN                            1    134.84 7888.7 867.13
# + MEDIAN_CV_COVERAGE             1    113.77 7909.8 867.78
# + DV200                          1     45.42 7978.1 869.89
# + Percent.Coding                 1     41.80 7981.7 870.00
# + Sex                            1     30.75 7992.8 870.34
# + LibrarySize                    1     28.70 7994.8 870.40
# - Dissector_AsNumeric            1   1150.99 9174.5 893.12
# 
# Step:  AIC=839.31
# PC6 ~ Dissector_AsNumeric + Percent.UTR
# 
# Df Sum of Sq    RSS    AIC
# + STGT_experience                1    316.22 6725.9 833.56
# + RibosomePerc                   1    302.35 6739.8 834.06
# <none>                                       7042.2 839.31
# + MEDIAN_CV_COVERAGE             1    156.29 6885.9 839.32
# + HPC_Dissection_Day             1    154.85 6887.3 839.37
# + Percent.Intronic               1    111.00 6931.2 840.92
# + RIN                            1     94.05 6948.1 841.52
# + SequencingBatch                2    247.73 6794.4 841.54
# + Percent.Intergenic             1     89.96 6952.2 841.67
# + RNAconc                        1     85.85 6956.3 841.81
# + Percent.Coding                 1     56.20 6986.0 842.85
# + HPC_Dissection_Date_Collapsed  5    651.13 6391.0 843.05
# + LibrarySize                    1     14.93 7027.2 844.30
# + Sex                            1     12.02 7030.1 844.40
# + DV200                          1      3.95 7038.2 844.68
# - Percent.UTR                    1    981.38 8023.5 865.78
# - Dissector_AsNumeric            1   1285.71 8327.9 874.90
# 
# Step:  AIC=833.56
# PC6 ~ Dissector_AsNumeric + Percent.UTR + STGT_experience
# 
# Df Sum of Sq    RSS    AIC
# + RibosomePerc                   1    149.65 6576.3 833.55
# <none>                                       6725.9 833.56
# + MEDIAN_CV_COVERAGE             1     85.89 6640.1 835.91
# + Percent.Intronic               1     74.58 6651.4 836.33
# + RIN                            1     71.81 6654.1 836.43
# + HPC_Dissection_Day             1     34.60 6691.3 837.80
# + Percent.Intergenic             1     26.80 6699.1 838.08
# + Percent.Coding                 1     16.46 6709.5 838.46
# + RNAconc                        1      8.82 6717.1 838.74
# + Sex                            1      3.97 6722.0 838.92
# + DV200                          1      2.99 6723.0 838.95
# + LibrarySize                    1      0.06 6725.9 839.06
# - STGT_experience                1    316.22 7042.2 839.31
# + SequencingBatch                2     64.29 6661.7 842.21
# + HPC_Dissection_Date_Collapsed  5    359.96 6366.0 847.59
# - Percent.UTR                    1    787.60 7513.5 855.19
# - Dissector_AsNumeric            1   1072.60 7798.6 864.31
# 
# Step:  AIC=833.55
# PC6 ~ Dissector_AsNumeric + Percent.UTR + STGT_experience + RibosomePerc
# 
# Df Sum of Sq    RSS    AIC
# <none>                                       6576.3 833.55
# - RibosomePerc                   1    149.65 6725.9 833.56
# - STGT_experience                1    163.51 6739.8 834.06
# + RIN                            1     96.33 6480.0 835.43
# + RNAconc                        1     36.44 6539.9 837.69
# + Percent.Intronic               1     35.44 6540.9 837.73
# + Percent.Intergenic             1     16.25 6560.1 838.44
# + HPC_Dissection_Day             1      9.27 6567.0 838.70
# + DV200                          1      4.84 6571.5 838.87
# + Percent.Coding                 1      4.77 6571.5 838.87
# + LibrarySize                    1      1.76 6574.5 838.98
# + MEDIAN_CV_COVERAGE             1      0.92 6575.4 839.02
# + Sex                            1      0.83 6575.5 839.02
# + SequencingBatch                2     33.22 6543.1 843.31
# - Percent.UTR                    1    531.74 7108.0 847.10
# + HPC_Dissection_Date_Collapsed  5    243.87 6332.4 851.80
# - Dissector_AsNumeric            1   1170.59 7746.9 868.18
# 
# Call:
#   lm(formula = PC6 ~ Dissector_AsNumeric + Percent.UTR + STGT_experience + 
#        RibosomePerc, data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)  Dissector_AsNumeric          Percent.UTR  STGT_experienceTRUE         RibosomePerc  
# -50.427               -5.593              174.737               -2.394              226.719  

summary.lm(lm(PC6 ~ Dissector_AsNumeric + Percent.UTR + STGT_experience + RibosomePerc, data = F2_PCA_wMetaData))

# Call:
#   lm(formula = PC6 ~ Dissector_AsNumeric + Percent.UTR + STGT_experience + 
#        RibosomePerc, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -15.5670  -3.0717   0.1145   3.5047  13.6908 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         -50.4272    13.3624  -3.774 0.000203 ***
#   Dissector_AsNumeric  -5.5931     0.8557  -6.536 3.76e-10 ***
#   Percent.UTR         174.7370    39.6663   4.405 1.59e-05 ***
#   STGT_experienceTRUE  -2.3940     0.9800  -2.443 0.015297 *  
#   RibosomePerc        226.7190    97.0154   2.337 0.020266 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 5.235 on 240 degrees of freedom
# Multiple R-squared:  0.2832,	Adjusted R-squared:  0.2713 
# F-statistic: 23.71 on 4 and 240 DF,  p-value: < 2.2e-16

step(object=lm(PC7~1, data=F2_PCA_wMetaData), scope=(~Sex+RibosomePerc+RNAconc+RIN+LibrarySize+SequencingBatch+Dissector_AsNumeric+HPC_Dissection_Date_Collapsed+HPC_Dissection_Day+STGT_experience+DV200+Percent.Intronic+Percent.Intergenic+Percent.Coding+Percent.UTR+MEDIAN_CV_COVERAGE), scale = 0, direction = "both", trace = 1, keep = NULL, steps = 1000, k=log(245))
# Start:  AIC=847.81
# PC7 ~ 1
# 
# Df Sum of Sq    RSS    AIC
# + Dissector_AsNumeric            1    532.57 7093.1 835.58
# + Percent.Coding                 1    489.65 7136.0 837.06
# + Percent.UTR                    1    387.66 7238.0 840.53
# <none>                                       7625.6 847.81
# + LibrarySize                    1    109.80 7515.8 849.76
# + RibosomePerc                   1    106.69 7518.9 849.86
# + HPC_Dissection_Day             1    101.06 7524.6 850.05
# + Percent.Intergenic             1     69.91 7555.7 851.06
# + RIN                            1     39.80 7585.8 852.03
# + RNAconc                        1     33.20 7592.4 852.25
# + STGT_experience                1     20.97 7604.7 852.64
# + Percent.Intronic               1      9.86 7615.8 853.00
# + MEDIAN_CV_COVERAGE             1      2.59 7623.0 853.23
# + DV200                          1      1.92 7623.7 853.25
# + Sex                            1      0.24 7625.4 853.31
# + SequencingBatch                2    155.34 7470.3 853.77
# + HPC_Dissection_Date_Collapsed  5    415.01 7210.6 861.61
# 
# Step:  AIC=835.58
# PC7 ~ Dissector_AsNumeric
# 
# Df Sum of Sq    RSS    AIC
# + Percent.Coding                 1   1085.37 6007.7 800.39
# + Percent.Intergenic             1    497.55 6595.5 823.26
# + Percent.UTR                    1    333.17 6759.9 829.29
# <none>                                       7093.1 835.58
# + MEDIAN_CV_COVERAGE             1    130.35 6962.7 836.53
# + LibrarySize                    1     98.36 6994.7 837.66
# + RNAconc                        1     85.86 7007.2 838.10
# + STGT_experience                1     62.16 7030.9 838.92
# + RIN                            1     61.90 7031.2 838.93
# + Percent.Intronic               1     61.45 7031.6 838.95
# + HPC_Dissection_Day             1     29.63 7063.4 840.05
# + DV200                          1     18.55 7074.5 840.44
# + RibosomePerc                   1      0.67 7092.4 841.06
# + Sex                            1      0.14 7092.9 841.07
# + SequencingBatch                2    115.21 6977.8 842.57
# - Dissector_AsNumeric            1    532.57 7625.6 847.81
# + HPC_Dissection_Date_Collapsed  5    351.21 6741.8 850.64
# 
# Step:  AIC=800.39
# PC7 ~ Dissector_AsNumeric + Percent.Coding
# 
# Df Sum of Sq    RSS    AIC
# <none>                                       6007.7 800.39
# + STGT_experience                1    133.31 5874.4 800.39
# + MEDIAN_CV_COVERAGE             1     91.23 5916.5 802.14
# + RibosomePerc                   1     89.24 5918.5 802.23
# + Percent.Intronic               1     71.45 5936.2 802.96
# + DV200                          1     61.35 5946.3 803.38
# + Percent.Intergenic             1     25.77 5981.9 804.84
# + LibrarySize                    1     25.03 5982.7 804.87
# + Sex                            1     24.29 5983.4 804.90
# + Percent.UTR                    1     22.64 5985.1 804.97
# + HPC_Dissection_Day             1     12.99 5994.7 805.36
# + RIN                            1      0.90 6006.8 805.86
# + RNAconc                        1      0.04 6007.6 805.89
# + SequencingBatch                2     92.75 5914.9 807.58
# + HPC_Dissection_Date_Collapsed  5    325.71 5682.0 814.24
# - Percent.Coding                 1   1085.37 7093.1 835.58
# - Dissector_AsNumeric            1   1128.29 7136.0 837.06
# 
# Call:
#   lm(formula = PC7 ~ Dissector_AsNumeric + Percent.Coding, data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)  Dissector_AsNumeric       Percent.Coding  
# -56.302               -4.922              134.578  

summary.lm(lm(PC7 ~ Dissector_AsNumeric + Percent.Coding, data = F2_PCA_wMetaData))
# Call:
#   lm(formula = PC7 ~ Dissector_AsNumeric + Percent.Coding, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -15.0826  -3.0939  -0.2649   3.1608  12.9036 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         -56.3022     9.3851  -5.999 7.19e-09 ***
#   Dissector_AsNumeric  -4.9223     0.7301  -6.742 1.14e-10 ***
#   Percent.Coding      134.5776    20.3531   6.612 2.40e-10 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 4.982 on 242 degrees of freedom
# Multiple R-squared:  0.2122,	Adjusted R-squared:  0.2057 
# F-statistic: 32.59 on 2 and 242 DF,  p-value: 2.939e-13

step(object=lm(PC8~1, data=F2_PCA_wMetaData), scope=(~Sex+RibosomePerc+RNAconc+RIN+LibrarySize+SequencingBatch+Dissector_AsNumeric+HPC_Dissection_Date_Collapsed+HPC_Dissection_Day+STGT_experience+DV200+Percent.Intronic+Percent.Intergenic+Percent.Coding+Percent.UTR+MEDIAN_CV_COVERAGE), scale = 0, direction = "both", trace = 1, keep = NULL, steps = 1000, k=log(245))
# Start:  AIC=810.67
# PC8 ~ 1
# 
# Df Sum of Sq    RSS    AIC
# + Percent.Intronic               1   2481.36 4071.6 699.58
# + SequencingBatch                2   2094.09 4458.9 727.35
# + HPC_Dissection_Date_Collapsed  5   2125.44 4427.5 742.12
# + HPC_Dissection_Day             1    612.17 5940.8 792.15
# + Percent.UTR                    1    353.68 6199.3 802.58
# + RNAconc                        1    259.31 6293.7 806.28
# + RIN                            1    250.98 6302.0 806.61
# <none>                                       6553.0 810.67
# + STGT_experience                1    113.90 6439.1 811.88
# + MEDIAN_CV_COVERAGE             1     52.46 6500.5 814.21
# + DV200                          1     50.82 6502.2 814.27
# + LibrarySize                    1     16.13 6536.8 815.57
# + Percent.Intergenic             1      6.93 6546.0 815.92
# + Percent.Coding                 1      4.10 6548.9 816.02
# + RibosomePerc                   1      3.37 6549.6 816.05
# + Dissector_AsNumeric            1      1.10 6551.9 816.13
# + Sex                            1      0.76 6552.2 816.15
# 
# Step:  AIC=699.58
# PC8 ~ Percent.Intronic
# 
# Df Sum of Sq    RSS    AIC
# + SequencingBatch                2   1276.67 2794.9 618.41
# + HPC_Dissection_Date_Collapsed  5   1389.90 2681.7 624.78
# + Percent.Intergenic             1    785.07 3286.5 652.61
# + Percent.Coding                 1    684.21 3387.4 660.01
# + HPC_Dissection_Day             1    662.39 3409.2 661.58
# + MEDIAN_CV_COVERAGE             1    585.35 3486.3 667.06
# + RibosomePerc                   1    313.57 3758.0 685.45
# <none>                                       4071.6 699.58
# + Dissector_AsNumeric            1     79.88 3991.7 700.23
# + Percent.UTR                    1     35.68 4035.9 702.93
# + Sex                            1     25.48 4046.1 703.55
# + STGT_experience                1     25.33 4046.3 703.56
# + LibrarySize                    1      8.97 4062.6 704.55
# + DV200                          1      5.99 4065.6 704.72
# + RIN                            1      1.99 4069.6 704.97
# + RNAconc                        1      0.02 4071.6 705.08
# - Percent.Intronic               1   2481.36 6553.0 810.67
# 
# Step:  AIC=618.41
# PC8 ~ Percent.Intronic + SequencingBatch
# 
# Df Sum of Sq    RSS    AIC
# + Percent.Intergenic             1    524.14 2270.8 573.03
# + RibosomePerc                   1    508.86 2286.1 574.67
# + MEDIAN_CV_COVERAGE             1    352.11 2442.8 590.92
# + Percent.Coding                 1    266.63 2528.3 599.35
# + Dissector_AsNumeric            1    211.10 2583.8 604.67
# + RNAconc                        1    190.88 2604.1 606.58
# <none>                                       2794.9 618.41
# + STGT_experience                1     59.36 2735.6 618.65
# + DV200                          1     58.40 2736.5 618.74
# + Sex                            1     16.48 2778.5 622.46
# + RIN                            1      4.59 2790.4 623.51
# + HPC_Dissection_Day             1      4.07 2790.9 623.56
# + Percent.UTR                    1      0.62 2794.3 623.86
# + LibrarySize                    1      0.00 2794.9 623.91
# + HPC_Dissection_Date_Collapsed  5    131.14 2663.8 634.14
# - SequencingBatch                2   1276.67 4071.6 699.58
# - Percent.Intronic               1   1663.94 4458.9 727.35
# 
# Step:  AIC=573.03
# PC8 ~ Percent.Intronic + SequencingBatch + Percent.Intergenic
# 
# Df Sum of Sq    RSS    AIC
# + RibosomePerc                   1     76.39 2194.4 570.15
# <none>                                       2270.8 573.03
# + DV200                          1     34.65 2236.2 574.77
# + STGT_experience                1     33.64 2237.2 574.88
# + Sex                            1     12.97 2257.8 577.13
# + MEDIAN_CV_COVERAGE             1      9.08 2261.7 577.55
# + Percent.Coding                 1      9.00 2261.8 577.56
# + Dissector_AsNumeric            1      8.93 2261.9 577.57
# + RNAconc                        1      7.39 2263.4 577.73
# + LibrarySize                    1      6.02 2264.8 577.88
# + RIN                            1      1.26 2269.5 578.40
# + HPC_Dissection_Day             1      0.95 2269.9 578.43
# + Percent.UTR                    1      0.20 2270.6 578.51
# + HPC_Dissection_Date_Collapsed  5     79.57 2191.2 591.80
# - Percent.Intergenic             1    524.14 2794.9 618.41
# - SequencingBatch                2   1015.74 3286.5 652.61
# - Percent.Intronic               1   2176.00 4446.8 732.18
# 
# Step:  AIC=570.15
# PC8 ~ Percent.Intronic + SequencingBatch + Percent.Intergenic + 
#   RibosomePerc
# 
# Df Sum of Sq    RSS    AIC
# <none>                                       2194.4 570.15
# + RNAconc                        1     43.92 2150.5 570.70
# + DV200                          1     38.10 2156.3 571.36
# + STGT_experience                1     34.90 2159.5 571.72
# + Sex                            1     30.13 2164.3 572.26
# - RibosomePerc                   1     76.39 2270.8 573.03
# + RIN                            1     10.52 2183.9 574.47
# - Percent.Intergenic             1     91.67 2286.1 574.67
# + Dissector_AsNumeric            1      1.91 2192.5 575.44
# + HPC_Dissection_Day             1      1.23 2193.2 575.51
# + LibrarySize                    1      1.11 2193.3 575.53
# + MEDIAN_CV_COVERAGE             1      0.86 2193.5 575.55
# + Percent.UTR                    1      0.72 2193.7 575.57
# + Percent.Coding                 1      0.59 2193.8 575.58
# + HPC_Dissection_Date_Collapsed  5     70.76 2123.7 589.62
# - SequencingBatch                2   1074.80 3269.2 656.81
# - Percent.Intronic               1   2025.23 4219.6 724.84
# 
# Call:
#   lm(formula = PC8 ~ Percent.Intronic + SequencingBatch + Percent.Intergenic + 
#        RibosomePerc, data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)    Percent.Intronic    SequencingBatch2    SequencingBatch3  Percent.Intergenic        RibosomePerc  
# 34.670           -1060.698              -5.220              -1.377              66.834            -224.003  


summary.lm(lm(PC8 ~ Percent.Intronic + SequencingBatch + Percent.Intergenic + RibosomePerc, data = F2_PCA_wMetaData))
# 
# Call:
#   lm(formula = PC8 ~ Percent.Intronic + SequencingBatch + Percent.Intergenic + 
#        RibosomePerc, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -8.6063 -1.8221 -0.1921  1.6352 11.2272 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           34.6703     3.4229  10.129  < 2e-16 ***
#   Percent.Intronic   -1060.6980    71.4192 -14.852  < 2e-16 ***
#   SequencingBatch2      -5.2198     0.4875 -10.708  < 2e-16 ***
#   SequencingBatch3      -1.3769     0.5619  -2.450  0.01499 *  
#   Percent.Intergenic    66.8341    21.1517   3.160  0.00178 ** 
#   RibosomePerc        -224.0030    77.6577  -2.884  0.00428 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 3.03 on 239 degrees of freedom
# Multiple R-squared:  0.6651,	Adjusted R-squared:  0.6581 
# F-statistic: 94.94 on 5 and 239 DF,  p-value: < 2.2e-16

step(object=lm(PC9~1, data=F2_PCA_wMetaData), scope=(~Sex+RibosomePerc+RNAconc+RIN+LibrarySize+SequencingBatch+Dissector_AsNumeric+HPC_Dissection_Date_Collapsed+HPC_Dissection_Day+STGT_experience+DV200+Percent.Intronic+Percent.Intergenic+Percent.Coding+Percent.UTR+MEDIAN_CV_COVERAGE), scale = 0, direction = "both", trace = 1, keep = NULL, steps = 1000, k=log(245))
# Start:  AIC=762.64
# PC9 ~ 1
# 
# Df Sum of Sq    RSS    AIC
# + STGT_experience                1   239.594 5146.8 757.00
# + Dissector_AsNumeric            1   167.691 5218.7 760.39
# <none>                                       5386.4 762.64
# + LibrarySize                    1    79.984 5306.4 764.48
# + Percent.Intronic               1    32.363 5354.0 766.67
# + HPC_Dissection_Day             1    31.867 5354.5 766.69
# + RibosomePerc                   1    26.816 5359.6 766.92
# + MEDIAN_CV_COVERAGE             1    12.366 5374.0 767.58
# + RNAconc                        1     6.625 5379.8 767.84
# + RIN                            1     2.543 5383.8 768.03
# + Percent.UTR                    1     2.168 5384.2 768.05
# + Percent.Coding                 1     1.895 5384.5 768.06
# + DV200                          1     1.244 5385.1 768.09
# + Sex                            1     0.369 5386.0 768.13
# + Percent.Intergenic             1     0.121 5386.3 768.14
# + SequencingBatch                2     2.886 5383.5 773.51
# + HPC_Dissection_Date_Collapsed  5   300.790 5085.6 776.07
# 
# Step:  AIC=757
# PC9 ~ STGT_experience
# 
# Df Sum of Sq    RSS    AIC
# + Dissector_AsNumeric            1    118.63 5028.2 756.78
# <none>                                       5146.8 757.00
# + Percent.Intronic               1     56.25 5090.5 759.80
# + RNAconc                        1     39.89 5106.9 760.59
# + LibrarySize                    1     29.69 5117.1 761.08
# + MEDIAN_CV_COVERAGE             1     23.41 5123.4 761.38
# + Percent.UTR                    1     14.70 5132.1 761.80
# + Percent.Intergenic             1      5.35 5141.4 762.24
# + Percent.Coding                 1      4.56 5142.2 762.28
# + Sex                            1      3.82 5143.0 762.32
# + RibosomePerc                   1      2.04 5144.8 762.40
# + HPC_Dissection_Day             1      1.17 5145.6 762.44
# + DV200                          1      0.85 5145.9 762.46
# + RIN                            1      0.29 5146.5 762.48
# - STGT_experience                1    239.59 5386.4 762.64
# + SequencingBatch                2     88.54 5058.3 763.75
# + HPC_Dissection_Date_Collapsed  5    373.28 4773.5 766.06
# 
# Step:  AIC=756.78
# PC9 ~ STGT_experience + Dissector_AsNumeric
# 
# Df Sum of Sq    RSS    AIC
# <none>                                       5028.2 756.78
# - Dissector_AsNumeric            1   118.633 5146.8 757.00
# + RibosomePerc                   1    73.580 4954.6 758.67
# - STGT_experience                1   190.536 5218.7 760.39
# + LibrarySize                    1    30.757 4997.4 760.78
# + Percent.Intronic               1    27.748 5000.4 760.93
# + Percent.Intergenic             1    15.385 5012.8 761.53
# + Percent.UTR                    1     8.419 5019.7 761.87
# + HPC_Dissection_Day             1     8.321 5019.8 761.88
# + DV200                          1     5.460 5022.7 762.02
# + Percent.Coding                 1     4.655 5023.5 762.06
# + Sex                            1     1.980 5026.2 762.19
# + RIN                            1     1.906 5026.3 762.19
# + MEDIAN_CV_COVERAGE             1     0.864 5027.3 762.24
# + RNAconc                        1     0.556 5027.6 762.26
# + SequencingBatch                2    60.626 4967.5 764.81
# + HPC_Dissection_Date_Collapsed  5   311.402 4716.8 768.63
# 
# Call:
#   lm(formula = PC9 ~ STGT_experience + Dissector_AsNumeric, data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)  STGT_experienceTRUE  Dissector_AsNumeric  
# 4.524               -2.410               -1.499  

summary.lm(lm(PC9 ~ STGT_experience + Dissector_AsNumeric, data = F2_PCA_wMetaData))
# 
# Call:
#   lm(formula = PC9 ~ STGT_experience + Dissector_AsNumeric, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -11.8812  -2.8978   0.2186   3.5831  10.1479 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           4.5241     1.1976   3.778 0.000199 ***
#   STGT_experienceTRUE  -2.4096     0.7957  -3.028 0.002726 ** 
#   Dissector_AsNumeric  -1.4986     0.6272  -2.389 0.017638 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 4.558 on 242 degrees of freedom
# Multiple R-squared:  0.06651,	Adjusted R-squared:  0.05879 
# F-statistic: 8.621 on 2 and 242 DF,  p-value: 0.0002418 

step(object=lm(PC10~1, data=F2_PCA_wMetaData), scope=(~Sex+RibosomePerc+RNAconc+RIN+LibrarySize+SequencingBatch+Dissector_AsNumeric+HPC_Dissection_Date_Collapsed+HPC_Dissection_Day+STGT_experience+DV200+Percent.Intronic+Percent.Intergenic+Percent.Coding+Percent.UTR+MEDIAN_CV_COVERAGE), scale = 0, direction = "both", trace = 1, keep = NULL, steps = 1000, k=log(245))
# Start:  AIC=726.62
# PC10 ~ 1
# 
# Df Sum of Sq    RSS    AIC
# + Percent.Intronic               1   181.189 4468.8 722.39
# <none>                                       4650.0 726.62
# + HPC_Dissection_Day             1    20.308 4629.7 731.05
# + LibrarySize                    1    17.031 4633.0 731.23
# + STGT_experience                1    14.011 4636.0 731.39
# + RibosomePerc                   1     7.759 4642.2 731.72
# + RIN                            1     7.332 4642.7 731.74
# + Dissector_AsNumeric            1     6.581 4643.4 731.78
# + RNAconc                        1     6.165 4643.8 731.80
# + Sex                            1     3.773 4646.2 731.93
# + MEDIAN_CV_COVERAGE             1     3.020 4647.0 731.97
# + Percent.Coding                 1     2.086 4647.9 732.02
# + DV200                          1     1.580 4648.4 732.04
# + Percent.UTR                    1     1.106 4648.9 732.07
# + Percent.Intergenic             1     0.016 4650.0 732.13
# + SequencingBatch                2    26.193 4623.8 736.24
# + HPC_Dissection_Date_Collapsed  5   107.327 4542.7 748.41
# 
# Step:  AIC=722.39
# PC10 ~ Percent.Intronic
# 
# Df Sum of Sq    RSS    AIC
# <none>                                       4468.8 722.39
# + Percent.Intergenic             1    68.220 4400.6 724.12
# + RIN                            1    61.237 4407.6 724.51
# + RNAconc                        1    51.447 4417.4 725.05
# + Percent.Coding                 1    36.419 4432.4 725.89
# + Dissector_AsNumeric            1    28.315 4440.5 726.33
# - Percent.Intronic               1   181.189 4650.0 726.62
# + HPC_Dissection_Day             1    22.799 4446.0 726.64
# + LibrarySize                    1    14.834 4454.0 727.08
# + MEDIAN_CV_COVERAGE             1     6.990 4461.8 727.51
# + Percent.UTR                    1     6.564 4462.2 727.53
# + RibosomePerc                   1     5.397 4463.4 727.59
# + STGT_experience                1     4.950 4463.9 727.62
# + Sex                            1     0.674 4468.1 727.85
# + DV200                          1     0.000 4468.8 727.89
# + SequencingBatch                2    48.608 4420.2 730.71
# + HPC_Dissection_Date_Collapsed  5   119.990 4348.8 743.23
# 
# Call:
#   lm(formula = PC10 ~ Percent.Intronic, data = F2_PCA_wMetaData)
# 
# Coefficients:
#   (Intercept)  Percent.Intronic  
# 9.286          -254.221  

summary.lm(lm(PC10 ~ Percent.Intronic, data = F2_PCA_wMetaData))

# Call:
#   lm(formula = PC10 ~ Percent.Intronic, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -11.9209  -3.0057  -0.0651   2.9644  16.7984 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)         9.286      2.971   3.125  0.00199 **
#   Percent.Intronic -254.221     80.991  -3.139  0.00191 **
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 4.288 on 243 degrees of freedom
# Multiple R-squared:  0.03897,	Adjusted R-squared:  0.03501 
# F-statistic: 9.853 on 1 and 243 DF,  p-value: 0.001906# 


################################################

#Why library size might not be popping up as related to the PCs in the larger models - it correlates with other variables in the dataset:

pdf("F2_LibrarySize_vs_HPC_Dissection_Day.pdf", height=6, width=4)
boxplot(LibrarySize~HPC_Dissection_Day, data=F2_PCA_wMetaData, main="F2", ylab="Library Size")
stripchart(LibrarySize~HPC_Dissection_Day, data=F2_PCA_wMetaData, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')
dev.off()


summary.lm(lm(LibrarySize~HPC_Dissection_Day, data=F2_PCA_wMetaData))

# Call:
#   lm(formula = LibrarySize ~ HPC_Dissection_Day, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -5025407 -1832707  -372555  1118761  8992082 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        20149942     392316  51.362  < 2e-16 ***
#   HPC_Dissection_Day   246898      46183   5.346 2.05e-07 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 2572000 on 247 degrees of freedom
# Multiple R-squared:  0.1037,	Adjusted R-squared:  0.1001 
# F-statistic: 28.58 on 1 and 247 DF,  p-value: 2.045e-07


summary.lm(lm(LibrarySize~HPC_Dissection_Date_AsFactor, data=F2_PCA_wMetaData))
# Call:
#   lm(formula = LibrarySize ~ HPC_Dissection_Date_AsFactor, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -5433974 -1845827  -229113  1366364  9228311 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            20656140     748530  27.596  < 2e-16 ***
#   HPC_Dissection_Date_AsFactor12/10/2020  1609822     858132   1.876 0.061880 .  
# HPC_Dissection_Date_AsFactor12/11/2020  3681762     867698   4.243 3.15e-05 ***
#   HPC_Dissection_Date_AsFactor12/14/2020  3358945     931913   3.604 0.000381 ***
#   HPC_Dissection_Date_AsFactor12/2/2020    278465     972369   0.286 0.774837    
# HPC_Dissection_Date_AsFactor12/3/2020   1087138    1000266   1.087 0.278198    
# HPC_Dissection_Date_AsFactor12/4/2020    680281     883411   0.770 0.442023    
# HPC_Dissection_Date_AsFactor12/7/2020   1100187     879104   1.251 0.211981    
# HPC_Dissection_Date_AsFactor12/8/2020    731493     861145   0.849 0.396485    
# HPC_Dissection_Date_AsFactor12/9/2020    405017     875066   0.463 0.643899    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 2483000 on 239 degrees of freedom
# Multiple R-squared:  0.192,	Adjusted R-squared:  0.1615 
# F-statistic: 6.309 on 9 and 239 DF,  p-value: 5.209e-08

#Library size is a little more closely related to dissection date than dissection day

summary.lm(lm(LibrarySize~RNAconc, data=F2_PCA_wMetaData))

# Call:
#   lm(formula = LibrarySize ~ RNAconc, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -4251139 -2122803  -358272  1303415 10234009 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 22636102     344586  65.691   <2e-16 ***
#   RNAconc        -3707       1918  -1.933   0.0544 .  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 2696000 on 247 degrees of freedom
# Multiple R-squared:  0.0149,	Adjusted R-squared:  0.01091 
# F-statistic: 3.736 on 1 and 247 DF,  p-value: 0.05439

