# Geothermal_ESDA
This repository contains exploratory spatial data analysis (ESDA) functions and scripts that were used in Smith et al. (2019) "Exploratory Spatial Data Analysis for Geothermal Datasets: An Appalachian Basin Case Study"

These functions are developed for ESDA on geothermal spatial datasets and should be applicable to other spatial datasets.

This repository depends on functions located in the following repositories:

calvinwhealton -> geothermal_pfa -> outliers

jds485 -> Geothermal_DataAnalysis_CrossSections

The methods used in ESDA_Main.R include:

0) Identification and processing of data in the same spatial location
  i) Nugget semi-variance for data in same location
  ii) Plot of similarity for data in the same locations based on covariate (depth of measurement)

1) Local Median/Mean Deviation (LocalDevition.R)

2) Spatial Outlier Analysis
  i) Plots by depth slices

3) Q-Q plots

4) Nonparametric Local Outlier Analysis
  i) KS test on depth rank distribution of outliers in local neighborhoods
  ii) Chi-squared test for depth bins

Additional Contributions:
1) Discrete Color Function for Map Making (ColorFunctions.R)

2) Sensitivity analysis of outlier algorithm for this dataset (OutAlgoSensitivity.R)
