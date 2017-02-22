# Geothermal_ESDA
This repository contains exploratory data analysis (EDA) and exploratory spatial data analysis (ESDA) functions and scripts, which may be useful for cleaning geothermal datasets and other spatial datasets.

This repository depends on functions located in the following repositories:
calvinwhealton -> geothermal_pfa -> combining_metrics_2016 branch -> outliers
jds485 -> Geothermal_DataAnalysis_CrossSections

The methods of EDA and ESDA used in ESDA_Main.R include:

0) Identification and processing of data in the same spatial location
  i) Nugget semi-variance for data in same location
  ii) Plot of similarity for data in the same locations based on covariate (depth of measurement)

1) Local Median/Mean Deviation (LocalDevition.R)

2) Spatial Outlier Analysis
  i) Plots by depth slices

3) Q-Q plots

4) Nonparametric Local Outlier Analysis
  i) KS test on depth rank distribution in local neighborhoods
  ii) Chi-squared test for depth bins

Additional Contributions:
1) Discrete Color Function for Map Making (ColorFunctions.R)

2) Sensitivity analysis of outlier algorithm for this dataset (OutAlgoSensitivity.R)
