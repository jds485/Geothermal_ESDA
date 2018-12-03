# Geothermal_ESDA
This repository contains exploratory spatial data analysis (ESDA) functions and scripts that were used in Smith et al. "Exploratory Spatial Data Analysis for Geothermal Resource Assessments: An Appalachian Basin Case Study" paper submitted to Geothermics.

These functions are developed for ESDA on geothermal spatial datasets and are applicable to other spatial datasets.

This repository depends on functions located in the following repositories:

calvinwhealton -> geothermal_pfa -> outliers

jds485 -> Geothermal_DataAnalysis_CrossSections

The methods used in ESDA_Main.R include:

0) Identification and processing of data in the same spatial location
  i) Nugget semi-variance calculation for data in same exact spatial location
  ii) Plot of similarity for data in the same locations based on a covariate (depth of measurement)

1) Local Median/Mean Deviation (LocalDevition.R)

2) Local Spatial Outlier Analysis (see geothermal_pfa repository for R script containing the functions)
  i) Plots by depth slices

3) Q-Q plots

4) Nonparametric Local Outlier Analysis
  i) KS test on depth rank distribution of outliers in local neighborhoods
  ii) Chi-squared test for depth bins

Additional Contributions:
1) Discrete Color Function for Map Making (ColorFunctions.R)

2) Sensitivity analysis of outlier algorithm for this dataset (OutAlgoSensitivity.R)

3) Jackknife confidence intervals for semi-variograms (JackknifeConfInts.R)

4) Diagnostic plots to discover operators that had systematically rogue data (bias and/or variance) as a result of their data recording behavior. (OperatorDiagnostics.R)
