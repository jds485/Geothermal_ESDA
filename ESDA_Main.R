# Description----
#Code for Smith et al. paper on Exploratory Spatial Data Analysis for Geothermal Resource Assessments: An Appalachian Basin Case Study

#Outline of Script:
# Load libraries and data
#  Note: the Loading Code and Loading Data section of this script can be skipped, and the input data loaded from the ESDA_Input_DeviatedWells.Rdata file.
# Check for and remove points with negative geothermal gradients. 
# Check for wells with the same spatial location. Only the deepest well at the same spatial location is retained.
#  Special cases of different BHT at the same depth are handled by either assigning a more likely depth to the point, or averaging the data.
# Check for local median deviation and select a minimum depth for BHTs
# Check for potentially rogue operators
# Check for local spatial outliers and analyze by depth rank.
# Check the performance of the ESDA methods using semi-variance compuations.

# Libraries ----
library(sp) # map plots
library(rgdal) #spatial data reading/writing
library(GISTools) #map making tools
library(dgof) #ks test for discrete distributions
library(Hmisc) #minor tick marks
library(readxl) #for Excel data reading
library(changepoint) #for changepoint analysis on the minimum BHT depth cutoff
library(gstat) #for variogram analysis
library(foreach) #parallel for loops
library(doParallel) #parallel package
library(lattice) #for plotting variograms

# Loading Code from Repositories ----
#From this Github repository
setwd("C:\\Users\\jsmif\\Documents\\Cornell\\Research\\Publications\\ESDA\\ESDACode\\Geothermal_ESDA")
source('LocalDeviation.R')
source('ColorFunctions.R')
source('OperatorDiagnostics.R')
#From Geothermal_DataAnalysis_CrossSections Github repository
setwd('C:\\Users\\jsmif\\Documents\\Cornell\\Research\\Publications\\DOE Grant\\ThermalConductivity\\Geothermal_DataAnalysis_CrossSections\\Geothermal_DataAnalysis_CrossSections')
source("DealingWithDataInDuplicateLocations.R")
#From geothermal_pfa Github repository
setwd('C:\\Users\\jsmif\\Documents\\Cornell\\Research\\Publications\\DOE Grant\\CombiningRiskFactorCode\\geothermal_pfa\\outliers')
source('outlier_identification.R')

# Loading Data and Map Layers ----
#  Political boundaries----
setwd('C:/Users/jsmif/Documents/Cornell/Research/Masters - Spatial Assessment/GIS/Population Density/2013_us_state_500k')
States = readOGR(dsn=getwd(), layer="us_state_WGS", stringsAsFactors=FALSE)
NY = States[which(States$STATEFP == "36"),]
PA = States[which(States$STATEFP == "42"),]
WV = States[which(States$STATEFP == "54"),]
MD = States[which(States$STATEFP == "24"),]
KY = States[which(States$STATEFP == "21"),]
VA = States[which(States$STATEFP == "51"),]
setwd('C:/Users/jsmif/Documents/Cornell/Research/Masters - Spatial Assessment/GIS/Population Density/2013_us_countys_500k')
Counties = readOGR(dsn=getwd(), layer="us_county_WGS", stringsAsFactors=FALSE) 
rm(States)

#  Geologic Regions / Spatial Interpolation Regions----
InterpRegs = readOGR(dsn="C:/Users/jsmif/Documents/Cornell/Research/Masters - Spatial Assessment/Figures/BaseCorrWells", layer="AllSectionsMerged", stringsAsFactors=FALSE) 
InterpRegs = spTransform(InterpRegs, CRS = CRS("+init=epsg:4326"))
#50 km bounded regions within the potential field edges - only applies to the WV regions
setwd("C:\\Users\\jsmif\\Documents\\Cornell\\Research\\Publications\\DOE Grant\\InterpolationDataset\\NewBoundaries\\NAD_InterpolationBounds")
VR_Bounded = readOGR(dsn=getwd(), layer = 'BoundedVREdit6') 
VR_Bounded = spTransform(VR_Bounded, CRS = CRS("+init=epsg:4326"))
CWV_Bounded = readOGR(dsn = getwd(), layer = 'BoundedCWV_Edit3')
CWV_Bounded = spTransform(CWV_Bounded, CRS = CRS("+init=epsg:4326"))
MT_Bounded = readOGR(dsn = getwd(), layer = 'BoundedMT_Edit')
MT_Bounded = spTransform(MT_Bounded, CRS = CRS("+init=epsg:4326"))

#  Wells with surface heat flow and temperatures at depth calculated----
setwd("C:\\Users\\jsmif\\Documents\\Cornell\\Research\\Publications\\ESDA\\ESDACode\\ESDA_Results")
Wells = read.csv('EDAWells_AllTempsThicksConds_BaseCorr_2018.csv', stringsAsFactors=FALSE)
coordinates(Wells) = c('LongDgr', 'LatDegr')
#WGS84 coordinates
proj4string(Wells) = CRS("+init=epsg:4326")

#  Add operator data to the well database from the original AASG spreadsheet----
Operator = read.csv('Operators.csv', stringsAsFactors = FALSE)
Wells$Operator = ''
for (i in 1:nrow(Wells)){
  Wells$Operator[i] = Operator$Operator[which(Operator$StateID == Wells$StateID[i])]
}
rm(i)

#Waco Oil and Gas operator APIs 
# these were identified as "bad" data because many logs were taken upwards, but interpreted as taken downward, 
# providing erroneously high temperatures in WV.
Wacos = read.csv('WacoOperatorAPIs.csv', stringsAsFactors = FALSE)
#WV APIs to match the Waco data
APIWV = read.csv('WV_APIs_StateIDs.csv', stringsAsFactors = FALSE)
#Get the state ID for each API and place in Waco
Wacos$StateID = ''
APIWV$APInum = as.numeric(strsplit(APIWV$API, split = '0000a', fixed = TRUE))
for (i in 1:nrow(Wacos)){
  if (length(which(APIWV$APInum == Wacos$API[i])) > 0){
    Wacos$StateID[i] = APIWV$StateID[APIWV$APInum == Wacos$API[i]]
  }
}
rm(i)

#  Horizontal and deviated wells for NY, PA, WV----
setwd("C:\\Users\\jsmif\\Documents\\Cornell\\Research\\Publications\\ESDA\\DirectionalWellInvestigations")

#Add directional well data from AASG database
AASG_DeviatedWells = read.csv('DirectionalWells_AASG.csv', stringsAsFactors = FALSE)
Wells$WellShape = ''
for (i in 1:nrow(Wells)){
  if(length(which(AASG_DeviatedWells$StateID == Wells$StateID[i])) > 0){
    Wells$WellShape[i] = AASG_DeviatedWells$WellBoreShape[which(AASG_DeviatedWells$StateID == Wells$StateID[i])]
  }
}
rm(i)

#NY - This database has some wells that are deviated by very little. May not be worth excluding these wells.
NY_DirWells = read.csv('DirectionalWellsNY.csv', stringsAsFactors = FALSE)
APINY = read.csv('NY_APIs.csv', stringsAsFactors = FALSE)
NY_DirWells$StateID = ''
for (i in 1:nrow(NY_DirWells)){
  if (length(which(APINY$API == NY_DirWells$API[i])) > 0){
    if(length(APINY$StateID[APINY$API == NY_DirWells$API[i]]) > 1){
      for (j in 1:length(APINY$StateID[APINY$API == NY_DirWells$API[i]])){
        if (j == 1){
          NY_DirWells$StateID[i] = APINY$StateID[APINY$API == NY_DirWells$API[i]][j]
        }else{
          NY_DirWells = rbind(NY_DirWells, NY_DirWells[i,])
          NY_DirWells$StateID[nrow(NY_DirWells)] = APINY$StateID[APINY$API == NY_DirWells$API[i]][j]
        }
      }
    }else{
      NY_DirWells$StateID[i] = APINY$StateID[APINY$API == NY_DirWells$API[i]]
    }
  }
}
rm(i,j)

#Only the wells that match with StateID number are in the AASG BHT dataset.
plot(Wells, pch = 16, cex = 0.3)
plot(Wells[Wells$StateID %in% NY_DirWells$StateID[NY_DirWells$StateID != ''],], pch = 16, cex = 0.3, col = 'red', add = T)
dev.off()

#Add a column for state database deviated wells
Wells$StateWellShape = ''
Wells@data[Wells$StateID %in% NY_DirWells$StateID[NY_DirWells$StateID != ''],]$StateWellShape = 'H'

#PA
PA_CDRs = read.csv('CDR_PALogsWithDirectionalLogs.csv', stringsAsFactors = FALSE)
PA_DirWells = read.csv('DirectionalWellsPA.csv', stringsAsFactors = FALSE)
APIPA = read.csv('PA_APIs.csv', stringsAsFactors = FALSE)

PA_CDRs$StateID = ''
PA_DirWells$StateID = ''
for (i in 1:nrow(PA_CDRs)){
  if (length(which(APIPA$API == PA_CDRs$API[i])) > 0){
    PA_CDRs$StateID[i] = APIPA$StateID[APIPA$API == PA_CDRs$API[i]]
  }
}
rm(i)
for (i in 1:nrow(PA_DirWells)){
  if (length(which(APIPA$API == PA_DirWells$API[i])) > 0){
    PA_DirWells$StateID[i] = APIPA$StateID[APIPA$API == PA_DirWells$API[i]]
  }
}
rm(i)

#Only the wells that match with StateID number are in the BHT dataset.
coordinates(PA_DirWells) = c('Longitude..Dec.', 'Latitude..Dec.')
proj4string(PA_DirWells) = CRS('+init=epsg:4326')

plot(Wells, pch = 16, cex = 0.3)
plot(Wells[Wells$StateID %in% PA_DirWells$StateID[PA_DirWells$StateID != ''],], pch = 16, cex = 0.3, col = 'red', add = T)
plot(Wells[Wells$StateID %in% PA_CDRs$StateID[PA_CDRs$StateID != ''],], pch = 16, cex = 0.3, col = 'red', add = T)
plot(PA_DirWells, add = T, col = 'blue')
dev.off()

#Add to column for state database deviated wells
Wells@data[Wells$StateID %in% PA_DirWells$StateID[PA_DirWells$StateID != ''],]$StateWellShape = 'H'
Wells@data[Wells$StateID %in% PA_CDRs$StateID[PA_CDRs$StateID != ''],]$StateWellShape = 'H'

#WV
WV_CompDirWells = read.csv('CompletedDirectionalWells.csv', stringsAsFactors = FALSE)
WV_PermDirWells = read.csv('PermittedDirectionalWells.csv', stringsAsFactors = FALSE)
WV_SurvsDirWells = read.csv('SurveysforDirectionalWells.csv', stringsAsFactors = FALSE)

#Get the state ID for each API and place in repsective dataset
WV_CompDirWells$StateID = ''
WV_PermDirWells$StateID = ''
WV_SurvsDirWells$StateID = ''
for (i in 1:nrow(WV_CompDirWells)){
  if (length(which(APIWV$APInum == WV_CompDirWells$API.Number[i])) > 0){
    WV_CompDirWells$StateID[i] = APIWV$StateID[APIWV$APInum == WV_CompDirWells$API.Number[i]]
  }
}
rm(i)
for (i in 1:nrow(WV_PermDirWells)){
  if (length(which(APIWV$APInum == WV_PermDirWells$API.Number[i])) > 0){
    WV_PermDirWells$StateID[i] = APIWV$StateID[APIWV$APInum == WV_PermDirWells$API.Number[i]]
  }
}
rm(i)
for (i in 1:nrow(WV_SurvsDirWells)){
  if (length(which(APIWV$APInum == WV_SurvsDirWells$API1[i])) > 0){
    WV_SurvsDirWells$StateID[i] = APIWV$StateID[APIWV$APInum == WV_SurvsDirWells$API1[i]]
  }
}
rm(i)

#Only the wells that match with StateID number are in the BHT dataset.
#WV_CompDirWells Contain all the unique directional wells for the 3 datasets. 
# All directional wells in BHT dataset were checked (16 total). Only a couple may be rogue entries.
#Check how far deviated at minimum these wells are
WV_CompDirWells$MinDev = sqrt((WV_CompDirWells$Surface.Loc.UTME - WV_CompDirWells$Btm.Hole.Loc.UTME)^2 + (WV_CompDirWells$Surface.Loc.UTMN - WV_CompDirWells$Btm.Hole.Loc.UTMN)^2)
hist(WV_CompDirWells$MinDev[WV_CompDirWells$StateID != ''])

plot(Wells, pch = 16, cex = 0.3)
plot(Wells[Wells$StateID %in% WV_CompDirWells$StateID[WV_CompDirWells$StateID != ''],], pch = 16, cex = 0.3, col = 'red', add = T)
dev.off()

#Add to column for state database deviated wells
Wells@data[Wells$StateID %in% WV_CompDirWells$StateID[WV_CompDirWells$StateID != ''],]$StateWellShape = 'H'

#Make a plot of the horizontal wells
png('DeviatedWells.png', res = 600, height = 6, width = 6, units = 'in')
par(mar = c(2,3,2,2))
plot(Wells, pch = 16, cex = 0.1)
plot(NY, add = T)
plot(PA, add = T)
plot(WV, add = T)
plot(MD, add = T)
plot(KY, add = T)
plot(VA, add = T)
north.arrow(xb = -75, yb = 37, len = 0.2, lab = 'N', col = 'black')
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
plot(Wells[grep(Wells$WellShape, pattern = 'ert'),], col = 'black', pch = 16, cex = 0.1, add = T)
plot(Wells[Wells$WellShape == '',], col = 'yellow', pch = 16, cex = 0.1, add = T)
plot(Wells[Wells$StateWellShape == 'H',], col = 'red', pch = 16, cex = 0.1, add = T)
plot(Wells[grep(Wells$WellShape, pattern = 'zon'),], col = 'green', pch = 16, cex = 0.1, add = T)
plot(Wells[grep(Wells$WellShape, pattern = 'via'),], col = 'green', pch = 16, cex = 0.1, add = T)
plot(Wells[grep(Wells$WellShape, pattern = 'Up'),], col = 'green', pch = 16, cex = 0.1, add = T)
plot(Wells[(Wells$DpthOfM - Wells$TruVrtc > 100) & Wells$TruVrtc > 0, ], col = 'blue', pch = 16, cex = 0.1, add = T)
legend('topleft', legend=c('Vertical', 'Not Specified', 'Deviated: State Data', 'Deviated: AASG Data', "BHT Depth - TVD > 20'"), col = c('black', 'yellow', 'red', 'green', 'blue'), pch = 16)
dev.off()

#Make a database that removes the deviated data of all types
Wells_NoDeviation = Wells[-unique(c(grep(Wells$WellShape, pattern = 'zon'), grep(Wells$WellShape, pattern = 'via'), grep(Wells$WellShape, pattern = 'Up'), which(Wells$StateWellShape == 'H'))),]
Wells_NoDeviation = Wells_NoDeviation[-which((Wells_NoDeviation$DpthOfM - Wells_NoDeviation$TruVrtc > 100) & Wells_NoDeviation$TruVrtc > 0), ]

#Because some NY wells are not deviated much, try keeping those in.
Wells_NoDeviationNY = Wells[-unique(c(grep(Wells$WellShape, pattern = 'zon'), grep(Wells$WellShape, pattern = 'via'), grep(Wells$WellShape, pattern = 'Up'), which((Wells$StateWellShape == 'H') & (Wells$State != 'NY')))),]
Wells_NoDeviationNY = Wells_NoDeviationNY[-which((Wells_NoDeviationNY$DpthOfM - Wells_NoDeviationNY$TruVrtc > 100) & Wells_NoDeviationNY$TruVrtc > 0), ]

#  Spicer equilibrium well temperature profiles - Saving for a different paper because these data are not public. Email authors for access.----
setwd('C:\\Users\\jsmif\\Documents\\Cornell\\Research\\Publications\\ESDA')
Spicer = read_xlsx(path = paste0(getwd(), '/EquilibriumTempProfiles.xlsx'), sheet = 'Spicers')
coordinates(Spicer) = c("Long", "Lat")
proj4string(Spicer) = CRS("+init=epsg:4326")

#Whealton MS thesis identified pseudo-equilibrium temperature profiles
Whealton = read.csv('EquilibriumTempProfiles_LocsAdded.csv')
coordinates(Whealton) = c("Long", "Lat")
proj4string(Whealton) = CRS("+init=epsg:4326")

#Whealton MS thesis identified BHTs
setwd('C:\\Users\\jsmif\\Documents\\Cornell\\Research\\Publications\\ESDA')
WhealtonBHTs = read.csv('WhealtonBHTsNYPA.csv', stringsAsFactors = FALSE)
WhealtonBHTs = WhealtonBHTs[is.na(WhealtonBHTs$BHT) == FALSE,]
WhealtonBHTs = WhealtonBHTs[is.na(WhealtonBHTs$LATITUDE) == FALSE,]
WhealtonBHTs = WhealtonBHTs[is.na(WhealtonBHTs$LONGITUDE) == FALSE,]

WhealtonBHTs$StateID = ''
for (i in 1:nrow(WhealtonBHTs)){
  if (length(which(APINY$API == WhealtonBHTs$API[i])) > 0){
    if(length(APINY$StateID[APINY$API == WhealtonBHTs$API[i]]) > 1){
      for (j in 1:length(APINY$StateID[APINY$API == WhealtonBHTs$API[i]])){
        if (j == 1){
          WhealtonBHTs$StateID[i] = APINY$StateID[APINY$API == WhealtonBHTs$API[i]][j]
        }else{
          WhealtonBHTs = rbind(WhealtonBHTs, WhealtonBHTs[i,])
          WhealtonBHTs$StateID[nrow(WhealtonBHTs)] = APINY$StateID[APINY$API == WhealtonBHTs$API[i]][j]
        }
      }
    }else{
      WhealtonBHTs$StateID[i] = APINY$StateID[APINY$API == WhealtonBHTs$API[i]]
    }
  }
}
rm(i,j)
for (i in 1:nrow(WhealtonBHTs)){
  if (length(which(APIPA$API == WhealtonBHTs$API[i])) > 0){
    if(length(APIPA$StateID[APIPA$API == WhealtonBHTs$API[i]]) > 1){
      for (j in 1:length(APIPA$StateID[APIPA$API == WhealtonBHTs$API[i]])){
        if (j == 1){
          WhealtonBHTs$StateID[i] = APIPA$StateID[APIPA$API == WhealtonBHTs$API[i]][j]
        }else{
          WhealtonBHTs = rbind(WhealtonBHTs, WhealtonBHTs[i,])
          WhealtonBHTs$StateID[nrow(WhealtonBHTs)] = APIPA$StateID[APIPA$API == WhealtonBHTs$API[i]][j]
        }
      }
    }else{
      WhealtonBHTs$StateID[i] = APIPA$StateID[APIPA$API == WhealtonBHTs$API[i]]
    }
  }
}
rm(i,j)

coordinates(WhealtonBHTs) = c('LONGITUDE', 'LATITUDE')
proj4string(WhealtonBHTs) = CRS('+init=epsg:4326')
plot(WhealtonBHTs, pch = 16, col = 'orange', cex = 0.2)
plot(Wells, pch = 16, cex = 0.2, add = T)

#Fixme: Cross check these wells for being deviated.
WhealtonBHTs$StateWellShape = ''
for (i in 1:nrow(WhealtonBHTs)){
  if (length(which(PA_CDRs$API == WhealtonBHTs$API[i])) > 0){
    WhealtonBHTs$StateWellShape[i] = 'H'
  }
}
rm(i)
for (i in 1:nrow(WhealtonBHTs)){
  if (length(which(PA_DirWells$API == WhealtonBHTs$API[i])) > 0){
    WhealtonBHTs$StateWellShape[i] = 'H'
  }
}
rm(i)

#Fixme: Check that temperatures and depths of the wells that are the same match the AASG database.

# Set working directory back to project directory----
setwd("C:\\Users\\jsmif\\Documents\\Cornell\\Research\\Publications\\ESDA\\ESDACode\\ESDA_Results")

#  Save input data to a file----
save.image("ESDA_Input_DeviatedWells.RData")

# Remove Negative Gradient Wells ----
NegsAll = Wells[Wells$Gradient <= 0,]

#Deepest well that has a negative gradient
MaxDepth_NegGrad = max(NegsAll$WellDepth)

#Location of wells that have negative Gradients.
png('WellsNegativeGradients.png', res = 300, height = 5, width = 5, units = 'in')
par(mar = c(2,3,2,2))
plot(Wells, pch=16, col='black', cex = 0.3)
plot(NY, add=TRUE)
plot(PA, add=TRUE)
plot(WV, add=TRUE)
plot(MD, add=TRUE)
plot(KY, add=TRUE)
plot(VA, add=TRUE)
plot(NegsAll, pch=16, col='red', add=TRUE, cex = 0.3)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('topleft', legend = c('Negative Gradient', 'Other'), col = c('red', 'black'), pch = 16)
dev.off()

#Map that identifies the wells listed in Table 1 that were inspected in detail
NegGradInspected = c('WV2003', 'WV3625', 'WV3621', 'WV3627', 'WV3647', 'WV2179', 'WV2181', 'WV907', 'WV66', 'NY3797', 'NY860', 'NY1001', 'NY5198', 'PA2814')
png('WellsNegativeGradients.png', res = 300, height = 5, width = 5, units = 'in')
par(mar = c(2,3,2,2))
plot(Wells, pch=16, col='black', cex = 0.3)
plot(NY, add=TRUE)
plot(PA, add=TRUE)
plot(WV, add=TRUE)
plot(MD, add=TRUE)
plot(KY, add=TRUE)
plot(VA, add=TRUE)
plot(NegsAll, pch=16, col='red', add=TRUE, cex = 0.3)
plot(Wells[Wells$StateID %in% NegGradInspected,], pch=16, col='blue', add=TRUE, cex = 0.3)
north.arrow(-75, 39, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('bottomright', legend = c('Gradient < 0', 'Gradient < 0 Source Checked', 'Gradient > 0'), col = c('red', 'blue', 'black'), pch = 16)
dev.off()

#Remove the negative gradient wells before the sorting of wells in the same spatial locations:
Wells_PosGrad = Wells[-which(Wells$Gradient <= 0),]
Wells_NoDeviation_PosGrad = Wells_NoDeviation[-which(Wells_NoDeviation$Gradient <= 0),]
Wells_NoDeviationNY_PosGrad = Wells_NoDeviationNY[-which(Wells_NoDeviationNY$Gradient <= 0),]

# Identify Wells in Same Spatial Location ----
#Note that this step is used here so that the QsDev function to calculate the
# local median surface heat flow uses only unique locations.

#Find all points that share the same location and take the deepest measurement.
Same = SameSpot(Wells_PosGrad)
SortData = SortingWells(Same$SameSpot, Same$StoreData_Frame, Wells_PosGrad, 'BHT', 'TruVrtc', 'DrllrTt', 'DpthOfM', 'WellDepth', 2)
Same_NoDev = SameSpot(Wells_NoDeviation_PosGrad)
SortData_NoDev = SortingWells(Same_NoDev$SameSpot, Same_NoDev$StoreData_Frame, Wells_NoDeviation_PosGrad, 'BHT', 'TruVrtc', 'DrllrTt', 'DpthOfM', 'WellDepth', 2)
Same_NoDevNY = SameSpot(Wells_NoDeviationNY_PosGrad)
SortData_NoDevNY = SortingWells(Same_NoDevNY$SameSpot, Same_NoDevNY$StoreData_Frame, Wells_NoDeviationNY_PosGrad, 'BHT', 'TruVrtc', 'DrllrTt', 'DpthOfM', 'WellDepth', 2)

write.csv(SortData$Sorted, "SortedUniqueSpots_AllTemps_ESDA_2018.csv")
write.csv(SortData$RerunWells, "RerunWells_AllTemps_ESDA_2018.csv")

#Number of records in full dataset
N_Wells_PosGrad = length(Wells_PosGrad)
#Number of locations with records in same spatial coordinates
N_LocsSameCoords = nrow(unique(Wells_PosGrad@coords[as.numeric(colnames(Same$StoreData_Frame)),]))
#Number of unique spatial locations after sorting
N_UniqueLocsPostSort = length(SortData$Sorted)
#Number of records in same spatial location
N_RecordsSameLoc = length(Same$StoreData_Frame)
#Number of locations with multiple BHTs at same depth that should be rerun because they did not have depth field information
N_RecordsMultipleBHTsSameDepth = length(SortData$RerunWells)

#58 points in 27 unique locations have a different BHT measurement at the same depth.
BHTsDiffSameDepth = length(unique(SortData$IndsDifferent))
BHTsDiffSameDepth_UniqueLocs = nrow(unique(Wells_PosGrad[SortData$IndsDifferent,]@coords))
#Used to see how many wells had a CensorTemp controlled output. Max of about 13 C
SortData_TestCensor = SortingWells(Same$SameSpot, Same$StoreData_Frame, Wells_PosGrad, 'BHT', 'TruVrtc', 'DrllrTt', 'DpthOfM', 'WellDepth', Inf)

#Locations dropped as a result of censoring
N_LocsDroppedTempCensor = nrow(unique(Wells_PosGrad@coords[unique(SortData$IndsCensTemp),]))

#The wells that are rerun should overwrite the last rows in SortedUniqueSopts
#Load in the wells that were rerun and add them to the SortData$Sorted
Rerun = read.csv('SortedUniqueSpots_AllTemps_ESDA_RerunAdded.csv', stringsAsFactors = FALSE)
Rerun = Rerun[c((nrow(Rerun) - (nrow(SortData$RerunWells) - 1)):nrow(Rerun)),]
Rerun$APINo = SortData$RerunWells$APINo
#Fixme: this is not generalized. Could be made cleaner.
SortData$Sorted@data[c((nrow(SortData$Sorted) - (nrow(SortData$RerunWells) - 1)):nrow(SortData$Sorted)),seq(1,ncol(Rerun)-4,1)] = Rerun[,-c(1,8,9,ncol(Rerun))]
SortData_NoDev$Sorted@data[c((nrow(SortData_NoDev$Sorted) - (nrow(SortData_NoDev$RerunWells) - 1)):nrow(SortData_NoDev$Sorted)),seq(1,ncol(Rerun)-4,1)] = Rerun[-nrow(Rerun),-c(1,8,9,ncol(Rerun))]
SortData_NoDevNY$Sorted@data[c((nrow(SortData_NoDevNY$Sorted) - (nrow(SortData_NoDevNY$RerunWells) - 1)):nrow(SortData_NoDevNY$Sorted)),seq(1,ncol(Rerun)-4,1)] = Rerun[,-c(1,8,9,ncol(Rerun))]

#Using IndsDeepSmallerBHT, check how many deeper data points have a greater temperature than shallower data points.
# Wells_PosGrad[as.numeric(colnames(Same$StoreData_Frame[which(Same$StoreData_Frame[which(as.numeric(colnames(Same$StoreData_Frame)) == unique(SortData$IndsDeepSmallerBHT)[3]),] == 1)])),]
#Check how many of these are greater than 2 degrees or so
Rows = vector('numeric')
for (i in 1:length(unique(SortData$IndsDeepSmallerBHT))){
  BHTs = Wells_PosGrad$BHT[as.numeric(colnames(Same$StoreData_Frame[which(Same$StoreData_Frame[which(as.numeric(colnames(Same$StoreData_Frame)) == unique(SortData$IndsDeepSmallerBHT)[i]),] == 1)]))]
  Depths = Wells_PosGrad$WellDepth[as.numeric(colnames(Same$StoreData_Frame[which(Same$StoreData_Frame[which(as.numeric(colnames(Same$StoreData_Frame)) == unique(SortData$IndsDeepSmallerBHT)[i]),] == 1)]))]
  if ((max(BHTs) - max(BHTs[Depths == max(Depths)])) > 2){
    Rows = c(Rows, Wells_PosGrad$RowID_[as.numeric(colnames(Same$StoreData_Frame[which(Same$StoreData_Frame[which(as.numeric(colnames(Same$StoreData_Frame)) == unique(SortData$IndsDeepSmallerBHT)[i]),] == 1)]))])
  }
}
N_DeepBHTSmaller = length(Rows)
rm(Rows,BHTs,i,Depths)

#  Line plots of the heat flow vs. depth of BHT measurement for the wells in the same spot----

#Fixme: Add equilibrium and pseudo-equilibrium well data to this plot, or make a new plot for these data

PlotSpots = function(Wells_PosGrad, #Must have a column named WellDepth
                     Same){
  # Make a copy of the database to track the wells in the same spatial location for this plot only.
  PlotSpots = Wells_PosGrad@data
  PlotSpots$LongDgr = Wells_PosGrad@coords[,1]
  PlotSpots$LatDegr = Wells_PosGrad@coords[,2]
  # The wells in the same spot will be assigned the same number in a field named SameSpot
  PlotSpots$SameSpot = NA
  
  #Track number of locations
  count = 1
  
  for (i in 1:nrow(Same$StoreData_Frame)){
    #Only take the unique spots that have more than 1 point
    if ((any(Same$StoreData_Frame[i,] == 1) & is.na(PlotSpots$SameSpot[as.numeric(colnames(Same$StoreData_Frame)[i])])) == TRUE){
      #Have not checked this spot yet. Gather all well indicies with the same spatial location.
      Indxs = as.numeric(colnames(Same$StoreData_Frame[which(Same$StoreData_Frame[i,] == 1)]))
      
      #Assign a number to these wells in the same spot
      PlotSpots$SameSpot[Indxs] = count
      count = count + 1
    }
  }
  rm(Indxs, count,i)
  
  #Sort data by the same spot location number
  PlotSpots = PlotSpots[order(PlotSpots$SameSpot),]
  
  #Retain only data in same spot as other data
  PlotSpots = PlotSpots[which(is.na(PlotSpots$SameSpot) == FALSE),]
  
  #Sort by the depth of the deepest well in the set of points
  for (i in 1:length(unique(PlotSpots$SameSpot))){
    #Determine how many wells there are in the same spot
    indxs = which(PlotSpots$SameSpot == i)
    #Sort only these wells and place the sorted data in those rows
    PlotSpots[indxs,] = PlotSpots[indxs,][order(PlotSpots$WellDepth[indxs]),]
  }
  rm(i, indxs)
  
  #Sort groups of wells by the shallowest well.
  #Obtain index of shallowest well for each location
  indxShallow = vector('numeric', length(unique(PlotSpots$SameSpot)))
  for (j in 1:length(unique(PlotSpots$SameSpot))){
    indxShallow[j] = which(PlotSpots$SameSpot == j)[1]
  }
  for (i in 1:length(unique(PlotSpots$SameSpot))){
    #Sort only these wells and place the sorted data in those rows
    NewInds = PlotSpots$SameSpot[indxShallow][order(PlotSpots$WellDepth[indxShallow])]
  }
  rm(indxShallow, i, j)
  
  #Make new database for the final plotting of points.
  #All PlotFinal records will be overwritten or deleted in the following for loop.
  PlotFinal = PlotSpots
  
  #index of data in the PlotFinal database - removing duplicate records.
  len = 1
  LocCheck = 0
  RemCheck = 0
  for (i in 1:length(unique(PlotSpots$SameSpot))){
    #Determine how many wells there are in the same spot
    indxs = which(PlotSpots$SameSpot == NewInds[i])
    
    if (anyDuplicated(PlotSpots$BHT[indxs]) != 0){
      #Find the indices that have same BHT
      res = which(PlotSpots$BHT[indxs] %in% unique(PlotSpots$BHT[indxs][duplicated(PlotSpots$BHT[indxs])]) == TRUE)
      #Of those, find indices with the same well depth
      dpth = which(PlotSpots$WellDepth[indxs][res] %in% unique(PlotSpots$WellDepth[indxs][res][duplicated(PlotSpots$WellDepth[indxs][res])]) == TRUE)
      
      #If they all have the same BHT, check if they have the same depth.
      if ((length(res) == length(indxs)) & (length(dpth) == length(indxs))){
        #Check if there are two or more sets of duplicate records (e.g. 4 total records, 2 duplicates)
        if (nrow(unique(PlotSpots[indxs, c('WellDepth', 'BHT')])) > 1){
          #Retain the unique records, and drop the remaining
          #Drop
          Drop = nrow(PlotSpots[indxs, c('WellDepth', 'BHT')]) - nrow(unique(PlotSpots[indxs, c('WellDepth', 'BHT')]))
          PlotFinal = PlotFinal[-((nrow(PlotFinal) - (Drop-1)):nrow(PlotFinal)),]
          
          #Overwrite the records in PlotFinal. These are the new data.
          indxs = indxs[which(rownames(PlotSpots[indxs, c('WellDepth', 'BHT')]) %in% rownames(unique(PlotSpots[indxs, c('WellDepth', 'BHT')])))]
          PlotFinal[len:(len + (length(indxs)-1)),] = PlotSpots[indxs,]
          len = len + length(indxs)
          
          RemCheck = RemCheck + Drop
          
        }else{
          LocCheck = LocCheck + 1
          #Do not record this in the PlotFinal database. Well has only 1 measurement. Remove the last length(indxs) rows. They are not needed.
          PlotFinal = PlotFinal[-((nrow(PlotFinal) - (length(indxs)-1)):nrow(PlotFinal)),]
          
          RemCheck = RemCheck + length(indxs)
        }
      }else{
        if (length(dpth) > 0){
          #Remove duplicate records, except 1. Keep all others.
          dups = which(duplicated(x = PlotSpots[indxs,c('BHT', 'WellDepth')]) == TRUE)
          indxs = indxs[-dups]
          # Remove the last length(dups) rows from the database. They are not needed.
          PlotFinal = PlotFinal[-((nrow(PlotFinal) - (length(dups)-1)):nrow(PlotFinal)),]
          
          RemCheck = RemCheck + length(dups)
        }
        #Overwrite the records in PlotFinal. These are the new data.
        PlotFinal[len:(len + (length(indxs)-1)),] = PlotSpots[indxs,]
        len = len + length(indxs)
      }
    }else{
      #Overwrite the records in PlotFinal. These are the new data.
      PlotFinal[len:(len + (length(indxs)-1)),] = PlotSpots[indxs,]
      len = len + length(indxs)
    }
  }
  rm(PlotSpots, i, len, indxs, dpth, res, NewInds, Drop, dups)
  
  return(PlotFinal)
}

PlotFinal = PlotSpots(Wells_PosGrad, Same)
PlotFinal_NoDev = PlotSpots(Wells_NoDeviation_PosGrad, Same_NoDev)
PlotFinal_NoDevNY = PlotSpots(Wells_NoDeviationNY_PosGrad, Same_NoDevNY)

#Colors by location, sorted by the highest Qs to lowest in shallowest measurement
PlotColPal = colorRampPalette(colors = c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))
cols = PlotColPal(length(unique(PlotFinal$SameSpot)))
cols_NoDev = PlotColPal(length(unique(PlotFinal_NoDev$SameSpot)))
cols_NoDevNY = PlotColPal(length(unique(PlotFinal_NoDevNY$SameSpot)))

#Make plot
png('SameSpotWells_QsVsDepth_colrev.png', res = 600, units = 'in', width = 7, height = 7)
par(mar = c(4.5, 5, 1.5, 1.5), xaxs='i', yaxs='i')
for (i in 1:length(unique(PlotFinal$SameSpot))){
  if (i == 1){
    plot(PlotFinal$WellDepth[which(PlotFinal$SameSpot == PlotFinal$SameSpot[length(unique(PlotFinal$SameSpot)) + 1 - i])], PlotFinal$Qs[which(PlotFinal$SameSpot == PlotFinal$SameSpot[length(unique(PlotFinal$SameSpot)) + 1 - i])], type = 'o', col = cols[length(unique(PlotFinal$SameSpot)) + 1 - i], pch = 16, xlim = c(0,3500), ylim = c(0,250), xlab = 'BHT Depth (m)', ylab = expression('Surface Heat Flow' ~ (mW/m^2)), cex.axis = 1.5, cex.lab = 1.5)
  }else{
    plot(PlotFinal$WellDepth[which(PlotFinal$SameSpot == PlotFinal$SameSpot[length(unique(PlotFinal$SameSpot)) + 1 - i])], PlotFinal$Qs[which(PlotFinal$SameSpot == PlotFinal$SameSpot[length(unique(PlotFinal$SameSpot)) + 1 - i])], type = 'o', col = cols[length(unique(PlotFinal$SameSpot)) + 1 - i], pch = 16, xlim = c(0,3500), ylim = c(0,250), axes = FALSE, xlab = '', ylab = '')
  }
  par(new = TRUE)
}
par(new = FALSE)
minor.tick(nx=5,ny=5)
legend('topright', legend = c('Shallowest BHT for Location is Shallow', '', '', 'Shallowest BHT for Location is Deep'), col = c('red', 'yellow', 'green', 'purple'), pch = 16, lty = 1)
dev.off()

png('SameSpotWells_QsVsDepth_colrev_NoDevWells.png', res = 600, units = 'in', width = 7, height = 7)
par(mar = c(4.5, 5, 1.5, 1.5), xaxs='i', yaxs='i')
for (i in 1:length(unique(PlotFinal_NoDev$SameSpot))){
  if (i == 1){
    plot(PlotFinal_NoDev$WellDepth[which(PlotFinal_NoDev$SameSpot == PlotFinal_NoDev$SameSpot[length(unique(PlotFinal_NoDev$SameSpot)) + 1 - i])], PlotFinal_NoDev$Qs[which(PlotFinal_NoDev$SameSpot == PlotFinal_NoDev$SameSpot[length(unique(PlotFinal_NoDev$SameSpot)) + 1 - i])], type = 'o', col = cols[length(unique(PlotFinal_NoDev$SameSpot)) + 1 - i], pch = 16, xlim = c(0,3500), ylim = c(0,250), xlab = 'BHT Depth (m)', ylab = expression('Surface Heat Flow' ~ (mW/m^2)), cex.axis = 1.5, cex.lab = 1.5)
  }else{
    plot(PlotFinal_NoDev$WellDepth[which(PlotFinal_NoDev$SameSpot == PlotFinal_NoDev$SameSpot[length(unique(PlotFinal_NoDev$SameSpot)) + 1 - i])], PlotFinal_NoDev$Qs[which(PlotFinal_NoDev$SameSpot == PlotFinal_NoDev$SameSpot[length(unique(PlotFinal_NoDev$SameSpot)) + 1 - i])], type = 'o', col = cols[length(unique(PlotFinal_NoDev$SameSpot)) + 1 - i], pch = 16, xlim = c(0,3500), ylim = c(0,250), axes = FALSE, xlab = '', ylab = '')
  }
  par(new = TRUE)
}
par(new = FALSE)
minor.tick(nx=5,ny=5)
legend('topright', legend = c('Shallowest BHT for Location is Shallow', '', '', 'Shallowest BHT for Location is Deep'), col = c('red', 'yellow', 'green', 'purple'), pch = 16, lty = 1)
dev.off()

png('SameSpotWells_QsVsDepth_colrev_NoDevWellsNY.png', res = 600, units = 'in', width = 7, height = 7)
par(mar = c(4.5, 5, 1.5, 1.5), xaxs='i', yaxs='i')
for (i in 1:length(unique(PlotFinal_NoDevNY$SameSpot))){
  if (i == 1){
    plot(PlotFinal_NoDevNY$WellDepth[which(PlotFinal_NoDevNY$SameSpot == PlotFinal_NoDevNY$SameSpot[length(unique(PlotFinal_NoDevNY$SameSpot)) + 1 - i])], PlotFinal_NoDevNY$Qs[which(PlotFinal_NoDevNY$SameSpot == PlotFinal_NoDevNY$SameSpot[length(unique(PlotFinal_NoDevNY$SameSpot)) + 1 - i])], type = 'o', col = cols[length(unique(PlotFinal_NoDevNY$SameSpot)) + 1 - i], pch = 16, xlim = c(0,3500), ylim = c(0,250), xlab = 'BHT Depth (m)', ylab = expression('Surface Heat Flow' ~ (mW/m^2)), cex.axis = 1.5, cex.lab = 1.5)
  }else{
    plot(PlotFinal_NoDevNY$WellDepth[which(PlotFinal_NoDevNY$SameSpot == PlotFinal_NoDevNY$SameSpot[length(unique(PlotFinal_NoDevNY$SameSpot)) + 1 - i])], PlotFinal_NoDevNY$Qs[which(PlotFinal_NoDevNY$SameSpot == PlotFinal_NoDevNY$SameSpot[length(unique(PlotFinal_NoDevNY$SameSpot)) + 1 - i])], type = 'o', col = cols[length(unique(PlotFinal_NoDevNY$SameSpot)) + 1 - i], pch = 16, xlim = c(0,3500), ylim = c(0,250), axes = FALSE, xlab = '', ylab = '')
  }
  par(new = TRUE)
}
par(new = FALSE)
minor.tick(nx=5,ny=5)
legend('topright', legend = c('Shallowest BHT for Location is Shallow', '', '', 'Shallowest BHT for Location is Deep'), col = c('red', 'yellow', 'green', 'purple'), pch = 16, lty = 1)
dev.off()

rm(PlotColPal, cols, i, cols_NoDev, cols_NoDevNY)

#Transform to spatial data
coordinates(PlotFinal) = c('LongDgr', 'LatDegr')
proj4string(PlotFinal) = CRS('+init=epsg:4326')
coordinates(PlotFinal_NoDev) = c('LongDgr', 'LatDegr')
proj4string(PlotFinal_NoDev) = CRS('+init=epsg:4326')
coordinates(PlotFinal_NoDevNY) = c('LongDgr', 'LatDegr')
proj4string(PlotFinal_NoDevNY) = CRS('+init=epsg:4326')

#   Map of which states have the most duplicate measurements----
png('Barplot_PointsSameSpatialLocationStates.png', res = 300, width = 8, height = 5, units = 'in')
layout(rbind(c(1,2)))
counts = table(PlotFinal$State)
#counts = table(Wells_PosGrad$State[as.numeric(colnames(Same$StoreData_Frame))])
barplot(counts, ylim = c(0,700), xlab = 'State', ylab = 'Frequency', main = 'Points Sharing Spatial Coordinates', cex.axis = 1.5, cex.lab = 1.5)
#map
#plot(Wells_PosGrad[as.numeric(colnames(Same$StoreData_Frame)),], pch = 16, cex = 0.3, col = 'grey')
plot(PlotFinal, pch = 16, cex = 0.3, col = 'white')
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(PlotFinal, pch = 16, cex = 0.3, col = 'grey', add = TRUE)
north.arrow(-75, 37, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
dev.off()
rm(counts)

#   Plots of BHTs in the same spot vs. unique spatial location----
#Sort by BHT 
PlotFinal = PlotFinal[rev(order(PlotFinal$BHT)),]
PlotFinal_NoDev = PlotFinal_NoDev[rev(order(PlotFinal_NoDev$BHT)),]
PlotFinal_NoDevNY = PlotFinal_NoDevNY[rev(order(PlotFinal_NoDevNY$BHT)),]

#Sort by Well Depth 
#PlotFinal = PlotFinal[rev(order(PlotFinal$WellDepth)),]

png('SameSpotWells_BHTVsLocation.png', res = 600, units = 'in', width = 7, height = 7)
par(mar = c(4.5, 5, 1.5, 1.5))
for (i in 1:length(unique(PlotFinal$SameSpot))){
  #Indices for the deepest wells at a location
  indsMaxDepth = which(PlotFinal$WellDepth[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i])] == max(PlotFinal$WellDepth[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i])]))
  if (i == 1){
    #Record used coordinates
    Coords = PlotFinal[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i]),][1,]
    #All others black
    plot(rep(i,length(PlotFinal$BHT[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i])])), PlotFinal$BHT[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i])], type = 'o', col = 'black', pch = 16, xlim = c(0,length(unique(PlotFinal$SameSpot))), ylim = c(0,140), xlab = 'Same Spot Well', ylab = expression(paste('BHT (', degree, 'C)')), cex.axis = 1.5, cex.lab = 1.5)
    par(new=TRUE)
    #Max depth wells red
    plot(rep(i, length(indsMaxDepth)), PlotFinal$BHT[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i])][indsMaxDepth], type = 'o', col = 'red', pch = 16, xlim = c(0,length(unique(PlotFinal$SameSpot))), ylim = c(0,140), xlab = '', ylab = '', axes=FALSE)
  }else if (nrow(zerodist(rbind(PlotFinal[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i]),][1,], Coords))) == 0){
    #Record used coordinates
    Coords = rbind(PlotFinal[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i]),][1,], Coords)
    #All others black
    plot(rep(i,length(PlotFinal$BHT[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i])])), PlotFinal$BHT[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i])], type = 'o', col = 'black', pch = 16, xlim = c(0,length(unique(PlotFinal$SameSpot))), ylim = c(0,140), xlab = '', ylab = '', axes = FALSE)
    par(new=TRUE)
    #Max depth wells red
    plot(rep(i, length(indsMaxDepth)), PlotFinal$BHT[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i])][indsMaxDepth], type = 'o', col = 'red', pch = 16, xlim = c(0,length(unique(PlotFinal$SameSpot))), ylim = c(0,140), axes = FALSE, xlab = '', ylab = '')
  }
  par(new = TRUE)
}
par(new = FALSE)
minor.tick(nx=5,ny=5)
dev.off()

png('SameSpotWells_BHTVsLocation_NoDeviation.png', res = 600, units = 'in', width = 7, height = 7)
par(mar = c(4.5, 5, 1.5, 1.5))
for (i in 1:length(unique(PlotFinal_NoDev$SameSpot))){
  #Indices for the deepest wells at a location
  indsMaxDepth = which(PlotFinal_NoDev$WellDepth[which(PlotFinal_NoDev$SameSpot == unique(PlotFinal_NoDev$SameSpot)[i])] == max(PlotFinal_NoDev$WellDepth[which(PlotFinal_NoDev$SameSpot == unique(PlotFinal_NoDev$SameSpot)[i])]))
  if (i == 1){
    #Record used coordinates
    Coords = PlotFinal_NoDev[which(PlotFinal_NoDev$SameSpot == unique(PlotFinal_NoDev$SameSpot)[i]),][1,]
    #All others black
    plot(rep(i,length(PlotFinal_NoDev$BHT[which(PlotFinal_NoDev$SameSpot == unique(PlotFinal_NoDev$SameSpot)[i])])), PlotFinal_NoDev$BHT[which(PlotFinal_NoDev$SameSpot == unique(PlotFinal_NoDev$SameSpot)[i])], type = 'o', col = 'black', pch = 16, xlim = c(0,length(unique(PlotFinal_NoDev$SameSpot))), ylim = c(0,140), xlab = 'Same Spot Well', ylab = expression(paste('BHT (', degree, 'C)')), cex.axis = 1.5, cex.lab = 1.5)
    par(new=TRUE)
    #Max depth wells red
    plot(rep(i, length(indsMaxDepth)), PlotFinal_NoDev$BHT[which(PlotFinal_NoDev$SameSpot == unique(PlotFinal_NoDev$SameSpot)[i])][indsMaxDepth], type = 'o', col = 'red', pch = 16, xlim = c(0,length(unique(PlotFinal_NoDev$SameSpot))), ylim = c(0,140), xlab = '', ylab = '', axes=FALSE)
  }else if (nrow(zerodist(rbind(PlotFinal_NoDev[which(PlotFinal_NoDev$SameSpot == unique(PlotFinal_NoDev$SameSpot)[i]),][1,], Coords))) == 0){
    #Record used coordinates
    Coords = rbind(PlotFinal_NoDev[which(PlotFinal_NoDev$SameSpot == unique(PlotFinal_NoDev$SameSpot)[i]),][1,], Coords)
    #All others black
    plot(rep(i,length(PlotFinal_NoDev$BHT[which(PlotFinal_NoDev$SameSpot == unique(PlotFinal_NoDev$SameSpot)[i])])), PlotFinal_NoDev$BHT[which(PlotFinal_NoDev$SameSpot == unique(PlotFinal_NoDev$SameSpot)[i])], type = 'o', col = 'black', pch = 16, xlim = c(0,length(unique(PlotFinal_NoDev$SameSpot))), ylim = c(0,140), xlab = '', ylab = '', axes = FALSE)
    par(new=TRUE)
    #Max depth wells red
    plot(rep(i, length(indsMaxDepth)), PlotFinal_NoDev$BHT[which(PlotFinal_NoDev$SameSpot == unique(PlotFinal_NoDev$SameSpot)[i])][indsMaxDepth], type = 'o', col = 'red', pch = 16, xlim = c(0,length(unique(PlotFinal_NoDev$SameSpot))), ylim = c(0,140), axes = FALSE, xlab = '', ylab = '')
  }
  par(new = TRUE)
}
par(new = FALSE)
minor.tick(nx=5,ny=5)
dev.off()
rm(indsMaxDepth, Coords, i)

#    Split by wells that have deeper BHT as smaller value----

#Figure out which wells have the deeper BHT as smaller in value
DeepBHTSmaller = function(PlotFinal){
  PlotFinal$IndsDeepSmallBHT = 0
  
  for (i in 1:length(unique(PlotFinal$SameSpot))){
    #Indices for the deepest wells at a location
    indsMaxDepth = which(PlotFinal$WellDepth[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i])] == max(PlotFinal$WellDepth[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i])]))
    if (nrow(PlotFinal[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i]),][-indsMaxDepth,]) >= 1){
      for (temp in 1:length(indsMaxDepth)){
        if (any(PlotFinal$BHT[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i])][indsMaxDepth[temp]] < PlotFinal$BHT[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i])][-indsMaxDepth])){
          PlotFinal$IndsDeepSmallBHT[which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i])] = 1
        }
      }
    }
  }
  rm(i, indsMaxDepth, temp)
  
  return(PlotFinal)
}

PlotFinal = DeepBHTSmaller(PlotFinal)
PlotFinal_DeepBHTsLarge = PlotFinal[PlotFinal$IndsDeepSmallBHT == 0,]
PlotFinal_DeepBHTsSmall = PlotFinal[PlotFinal$IndsDeepSmallBHT == 1,]

PlotFinal_NoDev = DeepBHTSmaller(PlotFinal_NoDev)
PlotFinal_NoDev_DeepBHTsLarge = PlotFinal_NoDev[PlotFinal_NoDev$IndsDeepSmallBHT == 0,]
PlotFinal_NoDev_DeepBHTsSmall = PlotFinal_NoDev[PlotFinal_NoDev$IndsDeepSmallBHT == 1,]

PlotFinal_NoDevNY = DeepBHTSmaller(PlotFinal_NoDevNY)
PlotFinal_NoDevNY_DeepBHTsLarge = PlotFinal_NoDevNY[PlotFinal_NoDevNY$IndsDeepSmallBHT == 0,]
PlotFinal_NoDevNY_DeepBHTsSmall = PlotFinal_NoDevNY[PlotFinal_NoDevNY$IndsDeepSmallBHT == 1,]

png('SameSpotWells_BHTVsLocation_SortDeepBHTSmaller.png', res = 600, units = 'in', width = 12, height = 6)
par(mar = c(4.5, 5, 1.5, 1.5))
layout(rbind(c(1,2)))

#Plot deep BHTs that are larger first
for (i in 1:length(unique(PlotFinal_DeepBHTsLarge$SameSpot))){
  #Indices for the deepest wells at a location
  indsMaxDepth = which(PlotFinal_DeepBHTsLarge$WellDepth[which(PlotFinal_DeepBHTsLarge$SameSpot == unique(PlotFinal_DeepBHTsLarge$SameSpot)[i])] == max(PlotFinal_DeepBHTsLarge$WellDepth[which(PlotFinal_DeepBHTsLarge$SameSpot == unique(PlotFinal_DeepBHTsLarge$SameSpot)[i])]))
  if (i == 1){
    #Record used coordinates
    Coords = PlotFinal_DeepBHTsLarge[which(PlotFinal_DeepBHTsLarge$SameSpot == unique(PlotFinal_DeepBHTsLarge$SameSpot)[i]),][1,]
    #All others black
    plot(rep(i,length(PlotFinal_DeepBHTsLarge$BHT[which(PlotFinal_DeepBHTsLarge$SameSpot == unique(PlotFinal_DeepBHTsLarge$SameSpot)[i])])), PlotFinal_DeepBHTsLarge$BHT[which(PlotFinal_DeepBHTsLarge$SameSpot == unique(PlotFinal_DeepBHTsLarge$SameSpot)[i])], type = 'o', col = 'black', pch = 16, xlim = c(0,length(unique(PlotFinal_DeepBHTsLarge$SameSpot))), ylim = c(0,140), xlab = 'Sorted Location ID', ylab = expression(paste('BHT (', degree, 'C)')), main = 'Deepest BHT is Largest', cex.axis = 1.5, cex.lab = 1.5)
    par(new=TRUE)
    #Max depth wells red
    plot(rep(i, length(indsMaxDepth)), PlotFinal_DeepBHTsLarge$BHT[which(PlotFinal_DeepBHTsLarge$SameSpot == unique(PlotFinal_DeepBHTsLarge$SameSpot)[i])][indsMaxDepth], type = 'o', col = 'red', pch = 16, xlim = c(0,length(unique(PlotFinal_DeepBHTsLarge$SameSpot))), ylim = c(0,140), xlab = '', ylab = '', axes=FALSE)
  }else if (nrow(zerodist(rbind(PlotFinal_DeepBHTsLarge[which(PlotFinal_DeepBHTsLarge$SameSpot == unique(PlotFinal_DeepBHTsLarge$SameSpot)[i]),][1,], Coords))) == 0){
    #Record used coordinates
    Coords = rbind(PlotFinal_DeepBHTsLarge[which(PlotFinal_DeepBHTsLarge$SameSpot == unique(PlotFinal_DeepBHTsLarge$SameSpot)[i]),][1,], Coords)
    #All others black
    plot(rep(i,length(PlotFinal_DeepBHTsLarge$BHT[which(PlotFinal_DeepBHTsLarge$SameSpot == unique(PlotFinal_DeepBHTsLarge$SameSpot)[i])])), PlotFinal_DeepBHTsLarge$BHT[which(PlotFinal_DeepBHTsLarge$SameSpot == unique(PlotFinal_DeepBHTsLarge$SameSpot)[i])], type = 'o', col = 'black', pch = 16, xlim = c(0,length(unique(PlotFinal_DeepBHTsLarge$SameSpot))), ylim = c(0,140), xlab = '', ylab = '', axes = FALSE)
    par(new=TRUE)
    #Max depth wells red
    plot(rep(i, length(indsMaxDepth)), PlotFinal_DeepBHTsLarge$BHT[which(PlotFinal_DeepBHTsLarge$SameSpot == unique(PlotFinal_DeepBHTsLarge$SameSpot)[i])][indsMaxDepth], type = 'o', col = 'red', pch = 16, xlim = c(0,length(unique(PlotFinal_DeepBHTsLarge$SameSpot))), ylim = c(0,140), axes = FALSE, xlab = '', ylab = '')
  }
  par(new = TRUE)
}
par(new = FALSE)
minor.tick(nx=5,ny=5)

#Plot deep BHTs that are smaller
for (i in 1:length(unique(PlotFinal_DeepBHTsSmall$SameSpot))){
  #Indices for the deepest wells at a location
  indsMaxDepth = which(PlotFinal_DeepBHTsSmall$WellDepth[which(PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i])] == max(PlotFinal_DeepBHTsSmall$WellDepth[which(PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i])]))
  if (i == 1){
    #Record used coordinates
    Coords = PlotFinal_DeepBHTsSmall[which(PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i]),][1,]
    #All others black
    plot(rep(i,length(PlotFinal_DeepBHTsSmall$BHT[which(PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i])])), PlotFinal_DeepBHTsSmall$BHT[which(PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i])], type = 'o', col = 'black', pch = 16, xlim = c(0,length(unique(PlotFinal_DeepBHTsSmall$SameSpot))), ylim = c(0,140), xlab = 'Sorted Location ID', ylab = expression(paste('BHT (', degree, 'C)')), main = 'Deepest BHT is Not Largest', cex.axis = 1.5, cex.lab = 1.5)
    par(new=TRUE)
    #Max depth wells red
    plot(rep(i, length(indsMaxDepth)), PlotFinal_DeepBHTsSmall$BHT[which(PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i])][indsMaxDepth], type = 'o', col = 'red', pch = 16, xlim = c(0,length(unique(PlotFinal_DeepBHTsSmall$SameSpot))), ylim = c(0,140), xlab = '', ylab = '', axes=FALSE)
  }else if (nrow(zerodist(rbind(PlotFinal_DeepBHTsSmall[which(PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i]),][1,], Coords))) == 0){
    #Record used coordinates
    Coords = rbind(PlotFinal_DeepBHTsSmall[which(PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i]),][1,], Coords)
    #All others black
    plot(rep(i,length(PlotFinal_DeepBHTsSmall$BHT[which(PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i])])), PlotFinal_DeepBHTsSmall$BHT[which(PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i])], type = 'o', col = 'black', pch = 16, xlim = c(0,length(unique(PlotFinal_DeepBHTsSmall$SameSpot))), ylim = c(0,140), xlab = '', ylab = '', axes = FALSE)
    par(new=TRUE)
    #Max depth wells red
    plot(rep(i, length(indsMaxDepth)), PlotFinal_DeepBHTsSmall$BHT[which(PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i])][indsMaxDepth], type = 'o', col = 'red', pch = 16, xlim = c(0,length(unique(PlotFinal_DeepBHTsSmall$SameSpot))), ylim = c(0,140), axes = FALSE, xlab = '', ylab = '')
  }
  par(new = TRUE)
}
par(new = FALSE)
minor.tick(nx=5,ny=5)
legend('topright', legend = c('Deepest BHTs', 'Other BHTs'), col = c('red', 'black'), pch = 16, lty = 1, cex = 1.5)
dev.off()

png('SameSpotWells_BHTVsLocation_SortDeepBHTSmaller_NoDeviation.png', res = 600, units = 'in', width = 12, height = 6)
par(mar = c(4.5, 5, 1.5, 1.5))
layout(rbind(c(1,2)))

#Plot deep BHTs that are larger first
for (i in 1:length(unique(PlotFinal_NoDev_DeepBHTsLarge$SameSpot))){
  #Indices for the deepest wells at a location
  indsMaxDepth = which(PlotFinal_NoDev_DeepBHTsLarge$WellDepth[which(PlotFinal_NoDev_DeepBHTsLarge$SameSpot == unique(PlotFinal_NoDev_DeepBHTsLarge$SameSpot)[i])] == max(PlotFinal_NoDev_DeepBHTsLarge$WellDepth[which(PlotFinal_NoDev_DeepBHTsLarge$SameSpot == unique(PlotFinal_NoDev_DeepBHTsLarge$SameSpot)[i])]))
  if (i == 1){
    #Record used coordinates
    Coords = PlotFinal_NoDev_DeepBHTsLarge[which(PlotFinal_NoDev_DeepBHTsLarge$SameSpot == unique(PlotFinal_NoDev_DeepBHTsLarge$SameSpot)[i]),][1,]
    #All others black
    plot(rep(i,length(PlotFinal_NoDev_DeepBHTsLarge$BHT[which(PlotFinal_NoDev_DeepBHTsLarge$SameSpot == unique(PlotFinal_NoDev_DeepBHTsLarge$SameSpot)[i])])), PlotFinal_NoDev_DeepBHTsLarge$BHT[which(PlotFinal_NoDev_DeepBHTsLarge$SameSpot == unique(PlotFinal_NoDev_DeepBHTsLarge$SameSpot)[i])], type = 'o', col = 'black', pch = 16, xlim = c(0,length(unique(PlotFinal_NoDev_DeepBHTsLarge$SameSpot))), ylim = c(0,140), xlab = 'Sorted Location ID', ylab = expression(paste('BHT (', degree, 'C)')), main = 'Deepest BHT is Largest', cex.axis = 1.5, cex.lab = 1.5)
    par(new=TRUE)
    #Max depth wells red
    plot(rep(i, length(indsMaxDepth)), PlotFinal_NoDev_DeepBHTsLarge$BHT[which(PlotFinal_NoDev_DeepBHTsLarge$SameSpot == unique(PlotFinal_NoDev_DeepBHTsLarge$SameSpot)[i])][indsMaxDepth], type = 'o', col = 'red', pch = 16, xlim = c(0,length(unique(PlotFinal_NoDev_DeepBHTsLarge$SameSpot))), ylim = c(0,140), xlab = '', ylab = '', axes=FALSE)
  }else if (nrow(zerodist(rbind(PlotFinal_NoDev_DeepBHTsLarge[which(PlotFinal_NoDev_DeepBHTsLarge$SameSpot == unique(PlotFinal_NoDev_DeepBHTsLarge$SameSpot)[i]),][1,], Coords))) == 0){
    #Record used coordinates
    Coords = rbind(PlotFinal_NoDev_DeepBHTsLarge[which(PlotFinal_NoDev_DeepBHTsLarge$SameSpot == unique(PlotFinal_NoDev_DeepBHTsLarge$SameSpot)[i]),][1,], Coords)
    #All others black
    plot(rep(i,length(PlotFinal_NoDev_DeepBHTsLarge$BHT[which(PlotFinal_NoDev_DeepBHTsLarge$SameSpot == unique(PlotFinal_NoDev_DeepBHTsLarge$SameSpot)[i])])), PlotFinal_NoDev_DeepBHTsLarge$BHT[which(PlotFinal_NoDev_DeepBHTsLarge$SameSpot == unique(PlotFinal_NoDev_DeepBHTsLarge$SameSpot)[i])], type = 'o', col = 'black', pch = 16, xlim = c(0,length(unique(PlotFinal_NoDev_DeepBHTsLarge$SameSpot))), ylim = c(0,140), xlab = '', ylab = '', axes = FALSE)
    par(new=TRUE)
    #Max depth wells red
    plot(rep(i, length(indsMaxDepth)), PlotFinal_NoDev_DeepBHTsLarge$BHT[which(PlotFinal_NoDev_DeepBHTsLarge$SameSpot == unique(PlotFinal_NoDev_DeepBHTsLarge$SameSpot)[i])][indsMaxDepth], type = 'o', col = 'red', pch = 16, xlim = c(0,length(unique(PlotFinal_NoDev_DeepBHTsLarge$SameSpot))), ylim = c(0,140), axes = FALSE, xlab = '', ylab = '')
  }
  par(new = TRUE)
}
par(new = FALSE)
minor.tick(nx=5,ny=5)

#Plot deep BHTs that are smaller
for (i in 1:length(unique(PlotFinal_NoDev_DeepBHTsSmall$SameSpot))){
  #Indices for the deepest wells at a location
  indsMaxDepth = which(PlotFinal_NoDev_DeepBHTsSmall$WellDepth[which(PlotFinal_NoDev_DeepBHTsSmall$SameSpot == unique(PlotFinal_NoDev_DeepBHTsSmall$SameSpot)[i])] == max(PlotFinal_NoDev_DeepBHTsSmall$WellDepth[which(PlotFinal_NoDev_DeepBHTsSmall$SameSpot == unique(PlotFinal_NoDev_DeepBHTsSmall$SameSpot)[i])]))
  if (i == 1){
    #Record used coordinates
    Coords = PlotFinal_NoDev_DeepBHTsSmall[which(PlotFinal_NoDev_DeepBHTsSmall$SameSpot == unique(PlotFinal_NoDev_DeepBHTsSmall$SameSpot)[i]),][1,]
    #All others black
    plot(rep(i,length(PlotFinal_NoDev_DeepBHTsSmall$BHT[which(PlotFinal_NoDev_DeepBHTsSmall$SameSpot == unique(PlotFinal_NoDev_DeepBHTsSmall$SameSpot)[i])])), PlotFinal_NoDev_DeepBHTsSmall$BHT[which(PlotFinal_NoDev_DeepBHTsSmall$SameSpot == unique(PlotFinal_NoDev_DeepBHTsSmall$SameSpot)[i])], type = 'o', col = 'black', pch = 16, xlim = c(0,length(unique(PlotFinal_NoDev_DeepBHTsSmall$SameSpot))), ylim = c(0,140), xlab = 'Sorted Location ID', ylab = expression(paste('BHT (', degree, 'C)')), main = 'Deepest BHT is Not Largest', cex.axis = 1.5, cex.lab = 1.5)
    par(new=TRUE)
    #Max depth wells red
    plot(rep(i, length(indsMaxDepth)), PlotFinal_NoDev_DeepBHTsSmall$BHT[which(PlotFinal_NoDev_DeepBHTsSmall$SameSpot == unique(PlotFinal_NoDev_DeepBHTsSmall$SameSpot)[i])][indsMaxDepth], type = 'o', col = 'red', pch = 16, xlim = c(0,length(unique(PlotFinal_NoDev_DeepBHTsSmall$SameSpot))), ylim = c(0,140), xlab = '', ylab = '', axes=FALSE)
  }else if (nrow(zerodist(rbind(PlotFinal_NoDev_DeepBHTsSmall[which(PlotFinal_NoDev_DeepBHTsSmall$SameSpot == unique(PlotFinal_NoDev_DeepBHTsSmall$SameSpot)[i]),][1,], Coords))) == 0){
    #Record used coordinates
    Coords = rbind(PlotFinal_NoDev_DeepBHTsSmall[which(PlotFinal_NoDev_DeepBHTsSmall$SameSpot == unique(PlotFinal_NoDev_DeepBHTsSmall$SameSpot)[i]),][1,], Coords)
    #All others black
    plot(rep(i,length(PlotFinal_NoDev_DeepBHTsSmall$BHT[which(PlotFinal_NoDev_DeepBHTsSmall$SameSpot == unique(PlotFinal_NoDev_DeepBHTsSmall$SameSpot)[i])])), PlotFinal_NoDev_DeepBHTsSmall$BHT[which(PlotFinal_NoDev_DeepBHTsSmall$SameSpot == unique(PlotFinal_NoDev_DeepBHTsSmall$SameSpot)[i])], type = 'o', col = 'black', pch = 16, xlim = c(0,length(unique(PlotFinal_NoDev_DeepBHTsSmall$SameSpot))), ylim = c(0,140), xlab = '', ylab = '', axes = FALSE)
    par(new=TRUE)
    #Max depth wells red
    plot(rep(i, length(indsMaxDepth)), PlotFinal_NoDev_DeepBHTsSmall$BHT[which(PlotFinal_NoDev_DeepBHTsSmall$SameSpot == unique(PlotFinal_NoDev_DeepBHTsSmall$SameSpot)[i])][indsMaxDepth], type = 'o', col = 'red', pch = 16, xlim = c(0,length(unique(PlotFinal_NoDev_DeepBHTsSmall$SameSpot))), ylim = c(0,140), axes = FALSE, xlab = '', ylab = '')
  }
  par(new = TRUE)
}
par(new = FALSE)
minor.tick(nx=5,ny=5)
legend('topright', legend = c('Deepest BHTs', 'Other BHTs'), col = c('red', 'black'), pch = 16, lty = 1, cex = 1.5)
dev.off()

rm(i, indsMaxDepth, Coords)

#   Check how many of the Deep BHTs that are smaller are more than 2 degrees different----
RowsDeepBHTSmallerBy2C_Counter = vector('numeric')
DiffBHTsSameDepth = vector('numeric')
for (i in 1:length(unique(PlotFinal_DeepBHTsSmall$SameSpot))){
  BHTs = PlotFinal_DeepBHTsSmall$BHT[PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i]]
  Depths = PlotFinal_DeepBHTsSmall$WellDepth[PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i]]
  DiffBHTsSameDepth = c(DiffBHTsSameDepth, (max(BHTs) - max(BHTs[Depths == max(Depths)])))
  if ((max(BHTs) - max(BHTs[Depths == max(Depths)])) > 2){
    RowsDeepBHTSmallerBy2C_Counter = c(RowsDeepBHTSmallerBy2C_Counter, PlotFinal_DeepBHTsSmall$RowID_[PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i]])
  }
}
rm(i, BHTs, RowsDeepBHTSmallerBy2C_Counter, Depths)

#Histogram of the differences between the deepest and shallower BHTs
hist(DiffBHTsSameDepth, breaks = 300)
dev.off()

#  Nugget Effect for wells in the same spatial location ----

#Fixme: Add nugget analysis for equilibrium wells. 

#Compute the Nugget Effect for points in the same spatial location.
#Make a data frame to store the locations, average nugget, number of nuggets calculated, min, max, and sd of the nugget
#Note: This is slower and provides the same result as using the PlotFinal, as done below.
# LocsNugs2 = matrix(0, ncol=9, nrow=1)
# colnames(LocsNugs2) = c('RowID_', 'POINT_X', 'POINT_Y', 'Nugget', 'Max', 'Min', 'Sd', 'PtPairs', 'NumPts')
# count=0
# #Mark the index with a 1 when it is used.
# IndsUsed = vector('numeric', length=nrow(Same$StoreData_Frame))
# for (i in 1:nrow(Same$StoreData_Frame)){
#   #Only take the unique spots that have more than 1 point
#   if (any(Same$StoreData_Frame[i,] == 1) & IndsUsed[i] != 1){
#     #Gather all well indicies with the same spatial location.
#     Indxs = as.numeric(colnames(Same$StoreData_Frame[which(Same$StoreData_Frame[i,] == 1)]))
#     #Mark that this spatial location has now been checked by marking the indices.
#     IndsUsed[which(colnames(Same$StoreData_Frame) %in% Indxs)] = 1
#     #Check to see if any of the values have the same heat flow. That means the record was a duplicate, and the nugget should not be counted for these.
#     Test = Wells_PosGrad$Qs[Indxs]
#     if (length(unique(Test)) != 1){
#       #There are unique BHTs for this well compute the nugget only for those wells that are unique records.
#       Nug = vector('numeric', length=length(unique(Test)))
#       VarioPts = matrix(0, nrow=length(Nug), ncol=length(Nug))
#       for (j in 1:length(Nug)){
#         VarioPts[j,] = ((unique(Test) - unique(Test)[j]))^2/2
#       }
#       Nug = VarioPts[lower.tri(VarioPts)]
#       #Store spatial location of point and nugget information
#       if (nrow(LocsNugs2) == 1 & Nug[1] != 0 & count == 0){
#         LocsNugs2[1,] = c(Wells_PosGrad$RowID_[Indxs[1]], Wells_PosGrad$LongDgr[Indxs[1]], Wells_PosGrad$LatDegr[Indxs[1]], mean(Nug), max(Nug), min(Nug), sd(Nug), length(Nug), length(unique(Test)))
#         count = 1
#       }
#       else if (count == 1){
#         LocsNugs2 = rbind(LocsNugs2, c(Wells_PosGrad$RowID_[Indxs[1]], Wells_PosGrad$LongDgr[Indxs[1]], Wells_PosGrad$LatDegr[Indxs[1]], mean(Nug), max(Nug), min(Nug), sd(Nug), length(Nug), length(unique(Test))))
#       }
#     }
#   }
# }
# rm(count, i, j, Indxs, Test, IndsUsed, Nug, VarioPts)
# 
# write.csv(LocsNugs2, 'NuggetLocations.csv')
LocsNugs = matrix(0, ncol=9, nrow=length(unique(PlotFinal$SameSpot)))
colnames(LocsNugs) = c('RowID_', 'POINT_X', 'POINT_Y', 'Nugget', 'Max', 'Min', 'Sd', 'PtPairs', 'NumPts')
for (i in 1:nrow(LocsNugs)){
  #Gather all well indicies with the same spatial location.
  Indxs = which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i])
  #Get the surface heat flow for these wells in the same location
  Test = PlotFinal$Qs[Indxs]
  #Compute nugget
  VarioPts = matrix(0, nrow=length(Test), ncol=length(Test))
  for (j in 1:length(Test)){
    VarioPts[j,] = ((Test - Test[j]))^2/2
  }
  Nug = VarioPts[lower.tri(VarioPts)]
  #Store spatial location of point and nugget information. Take only first index
  LocsNugs[i,] = c(PlotFinal$RowID_[Indxs[1]], PlotFinal@coords[Indxs[1],1], PlotFinal@coords[Indxs[1],2], mean(Nug), max(Nug), min(Nug), sd(Nug), length(Nug), length(Test))
}
rm(i, j, Indxs, Test, Nug, VarioPts)

write.csv(LocsNugs, 'NuggetLocations_2018.csv')

LocsNugs_Deeper1km = matrix(0, ncol=9, nrow=length(unique(PlotFinal$SameSpot)))
colnames(LocsNugs_Deeper1km) = c('RowID_', 'POINT_X', 'POINT_Y', 'Nugget', 'Max', 'Min', 'Sd', 'PtPairs', 'NumPts')
for (i in 1:nrow(LocsNugs_Deeper1km)){
  #Gather all well indicies with the same spatial location.
  Indxs = which(PlotFinal$SameSpot == unique(PlotFinal$SameSpot)[i])
  #Check if at least 2 points are deeper than 1 km
  if (length(which(PlotFinal$WellDepth[Indxs] >= 1000)) > 1){
    #Get the surface heat flow for these wells in the same location
    Test = PlotFinal$Qs[Indxs][PlotFinal$WellDepth[Indxs] >= 1000]
    #Compute nugget
    VarioPts = matrix(0, nrow=length(Test), ncol=length(Test))
    for (j in 1:length(Test)){
      VarioPts[j,] = ((Test - Test[j]))^2/2
    }
    Nug = VarioPts[lower.tri(VarioPts)]
    #Store spatial location of point and nugget information. Take only first index
    LocsNugs_Deeper1km[i,] = c(PlotFinal$RowID_[Indxs[1]], PlotFinal@coords[Indxs[1],1], PlotFinal@coords[Indxs[1],2], mean(Nug), max(Nug), min(Nug), sd(Nug), length(Nug), length(Test))
  }
}
LocsNugs_Deeper1km = LocsNugs_Deeper1km[LocsNugs_Deeper1km[,1] != 0,]
rm(i, j, Indxs, Test, Nug, VarioPts)

write.csv(LocsNugs_Deeper1km, 'NuggetLocations_Deeper1km_2018.csv')

#Add interpolation section to the nugget locations
LocsNugs = as.data.frame(LocsNugs)
coordinates(LocsNugs) = c('POINT_X', 'POINT_Y')
proj4string(LocsNugs) = CRS('+init=epsg:4326')

LocsNugs_Deeper1km = as.data.frame(LocsNugs_Deeper1km)
coordinates(LocsNugs_Deeper1km) = c('POINT_X', 'POINT_Y')
proj4string(LocsNugs_Deeper1km) = CRS('+init=epsg:4326')

LocsNugs$Reg = NA
LocsNugs$Reg[as.numeric(rownames(LocsNugs[InterpRegs[1,],]@data))] = InterpRegs$Name[1]
LocsNugs$Reg[as.numeric(rownames(LocsNugs[InterpRegs[2,],]@data))] = InterpRegs$Name[2]
LocsNugs$Reg[as.numeric(rownames(LocsNugs[InterpRegs[3,],]@data))] = InterpRegs$Name[3]
LocsNugs$Reg[as.numeric(rownames(LocsNugs[InterpRegs[5,],]@data))] = InterpRegs$Name[5]
LocsNugs$Reg[as.numeric(rownames(LocsNugs[InterpRegs[6,],]@data))] = InterpRegs$Name[6]
LocsNugs$Reg[as.numeric(rownames(LocsNugs[InterpRegs[8,],]@data))] = InterpRegs$Name[8]
LocsNugs$Reg[as.numeric(rownames(LocsNugs[InterpRegs[9,],]@data))] = InterpRegs$Name[9]
LocsNugs$Reg[as.numeric(rownames(LocsNugs[MT_Bounded,]@data))] = InterpRegs$Name[4]
LocsNugs$Reg[as.numeric(rownames(LocsNugs[CWV_Bounded,]@data))] = InterpRegs$Name[7]
LocsNugs$Reg[as.numeric(rownames(LocsNugs[VR_Bounded,]@data))] = 'VR'

LocsNugs_Deeper1km$Reg = NA
LocsNugs_Deeper1km$Reg[as.numeric(rownames(LocsNugs_Deeper1km[InterpRegs[1,],]@data))] = InterpRegs$Name[1]
LocsNugs_Deeper1km$Reg[as.numeric(rownames(LocsNugs_Deeper1km[InterpRegs[2,],]@data))] = InterpRegs$Name[2]
LocsNugs_Deeper1km$Reg[as.numeric(rownames(LocsNugs_Deeper1km[InterpRegs[3,],]@data))] = InterpRegs$Name[3]
LocsNugs_Deeper1km$Reg[as.numeric(rownames(LocsNugs_Deeper1km[InterpRegs[5,],]@data))] = InterpRegs$Name[5]
LocsNugs_Deeper1km$Reg[as.numeric(rownames(LocsNugs_Deeper1km[InterpRegs[6,],]@data))] = InterpRegs$Name[6]
LocsNugs_Deeper1km$Reg[as.numeric(rownames(LocsNugs_Deeper1km[InterpRegs[8,],]@data))] = InterpRegs$Name[8]
LocsNugs_Deeper1km$Reg[as.numeric(rownames(LocsNugs_Deeper1km[InterpRegs[9,],]@data))] = InterpRegs$Name[9]
LocsNugs_Deeper1km$Reg[as.numeric(rownames(LocsNugs_Deeper1km[MT_Bounded,]@data))] = InterpRegs$Name[4]
LocsNugs_Deeper1km$Reg[as.numeric(rownames(LocsNugs_Deeper1km[CWV_Bounded,]@data))] = InterpRegs$Name[7]
LocsNugs_Deeper1km$Reg[as.numeric(rownames(LocsNugs_Deeper1km[VR_Bounded,]@data))] = 'VR'

#One point is not in the interpolation regions. Will not be used in plots.
LocsNugs = LocsNugs[is.na(LocsNugs$Reg) == FALSE,]

#   Make boxplots for each of the interpolation sections----
png('NuggetWellsDistributions_2018.png', res=600, width=12, height=6, units='in')
par(mar=c(5,5,2,2), yaxt='n')
boxplot(Nugget ~ Reg, data = LocsNugs, at=c(4, 1, 8, 5, 6, 3, 7, 9, 2), varwidth = TRUE, col=c('green', 'red', 'purple', 'springgreen', 'skyblue', 'yellow', 'blue', 'grey', 'orange'), pch=16, cex.axis=1.5, cex.lab=1.5, ylab=expression('Sample Nugget Semivariance' ~ (mW/m^2)^2), xlab='Interpolation Region', log='y')
par(yaxt='s')
at.y <- outer(1:9, 10^(-7:7))
lab.y <- ifelse(log10(at.y) %% 1 == 0, sapply(log10(at.y), function(i) as.expression(bquote(10^ .(i)))), NA)
par(tcl = -0.2)
axis(side=2, at=at.y, labels=lab.y, cex.axis=1.5, las=1)
par(tcl = -0.5)
at.y <- 10^(-7:7)
axis(side=2, at=at.y, labels=FALSE, cex.axis=1.5, las=1)
dev.off()

png('NuggetWellsDistributions_Deeper1km_2018.png', res=600, width=12, height=6, units='in')
par(mar=c(5,5,2,2), yaxt='n')
boxplot(Nugget ~ Reg, data = LocsNugs_Deeper1km, at=c(4, 1, 8, 5, 6, 3, 7, 9, 2), varwidth = TRUE, col=c('green', 'red', 'white', 'springgreen', 'skyblue', 'yellow', 'white', 'grey', 'orange'), pch=16, cex.axis=1.5, cex.lab=1.5, ylab=expression('Sample Nugget Semivariance' ~ (mW/m^2)^2), xlab='Interpolation Region', log='y', ylim = c(0.00001,3000), xlim = c(0.5,9.5))
par(new=TRUE)
plot(x = c(7,7,7,7), y = LocsNugs_Deeper1km$Nugget[LocsNugs_Deeper1km$Reg == 'SWPA'], ylim = c(0.00001,3000), xlim = c(0.5,9.5), log='y', col = 'blue', pch = 16, axes = FALSE, ylab = '', xlab = '')
par(new=TRUE)
plot(x = c(8,8), y = LocsNugs_Deeper1km$Nugget[LocsNugs_Deeper1km$Reg == 'CWV'], ylim = c(0.00001,3000), xlim = c(0.5,9.5), log='y', col = 'purple', pch = 16, axes = FALSE, ylab = '', xlab = '')
par(yaxt='s')
at.y <- outer(1:9, 10^(-7:7))
lab.y <- ifelse(log10(at.y) %% 1 == 0, sapply(log10(at.y), function(i) as.expression(bquote(10^ .(i)))), NA)
par(tcl = -0.2)
axis(side=2, at=at.y, labels=lab.y, cex.axis=1.5, las=1)
par(tcl = -0.5)
at.y <- 10^(-7:7)
axis(side=2, at=at.y, labels=FALSE, cex.axis=1.5, las=1)
dev.off()

# With a map next to the boxplot
pin = par("pin")
dxy = apply(rbind(c(-82.64474, -74.5), c(36.75, 43.6)), 1, diff)
ratio = dxy[1]/dxy[2]
pin[1] = 3 #Modifying for margin space
png('NuggetWellsDistributions_Boxplot&Map_edit_2018.png', res=600, width=11.2, height=4.2, units='in')
layout(cbind(1,1,2))
par(mar=c(4.1,5,0.5,0), yaxt='n')
boxplot(Nugget ~ Reg, data = LocsNugs, at=c(4, 1, 8, 5, 6, 3, 7, 9, 2), varwidth = TRUE, col=c('green', 'red', 'purple', 'springgreen', 'skyblue', 'yellow', 'blue', 'grey', 'orange'), pch=16, cex.axis=1.5, cex.lab=1.5, ylab=expression('Sample Nugget Semivariance' ~ (mW/m^2)^2), xlab='Interpolation Region', log='y')
par(yaxt='s')
at.y <- outer(1:9, 10^(-5:4))
lab.y <- ifelse(log10(at.y) %% 1 == 0, sapply(log10(at.y),function(i) as.expression(bquote(10^ .(i)))), NA)
par(tcl = -0.2)
axis(side=2, at=at.y, labels=lab.y, cex.axis=1.5, las=1)
par(tcl = -0.5)
at.y <- 10^(-7:7)
axis(side=2, at=at.y, labels=FALSE, cex.axis=1.5, las=1)
#Map
par(xaxs = 'i', yaxs = 'i', mar = c(0,0,0,0))
par(pin = c(pin[1], ratio*pin[1]))
plot(InterpRegs, xlim = c(-82.64474, -74.5), ylim = c(36.75, 43.4))
plot(InterpRegs[which(InterpRegs$Name == "CT"),], col = 'red', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "WPA"),], col = 'orange', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "NWPANY"),], col = 'yellow', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "CNY"),], col = 'green', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "ENY"),], col = 'springgreen', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "ENYPA"),], col = 'skyblue', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "SWPA"),], col = 'blue', add = TRUE)
plot(CWV_Bounded, col = 'purple', add = TRUE)
plot(MT_Bounded, col = 'magenta', add = TRUE)
plot(VR_Bounded, col = 'gray', add = TRUE)
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(LocsNugs, pch = 16, cex = 0.5, add = TRUE)
north.arrow(-75, 37, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
dev.off()

png('NuggetWellsDistributions_Boxplot&Map_edit_Deeper1km_2018.png', res=600, width=11.2, height=4.2, units='in')
layout(cbind(1,1,2))
par(mar=c(4.1,5,0.5,0), yaxt='n')
boxplot(Nugget ~ Reg, data = LocsNugs_Deeper1km, at=c(4, 1, 8, 5, 6, 3, 7, 9, 2), varwidth = TRUE, col=c('green', 'red', 'white', 'springgreen', 'skyblue', 'yellow', 'white', 'grey', 'orange'), pch=16, cex.axis=1.5, cex.lab=1.5, ylab=expression('Sample Nugget Semivariance' ~ (mW/m^2)^2), xlab='Interpolation Region', log='y', ylim = c(0.00001,3000))
par(new=TRUE)
plot(x = c(7,7,7,7), y = LocsNugs_Deeper1km$Nugget[LocsNugs_Deeper1km$Reg == 'SWPA'], ylim = c(0.00001,3000), xlim = c(0.5,9.5), log='y', col = 'blue', pch = 16, axes = FALSE, ylab = '', xlab = '')
par(new=TRUE)
plot(x = c(8,8), y = LocsNugs_Deeper1km$Nugget[LocsNugs_Deeper1km$Reg == 'CWV'], ylim = c(0.00001,3000), xlim = c(0.5,9.5), log='y', col = 'purple', pch = 16, axes = FALSE, ylab = '', xlab = '')
par(yaxt='s')
at.y <- outer(1:9, 10^(-5:4))
lab.y <- ifelse(log10(at.y) %% 1 == 0, sapply(log10(at.y),function(i) as.expression(bquote(10^ .(i)))), NA)
par(tcl = -0.2)
axis(side=2, at=at.y, labels=lab.y, cex.axis=1.5, las=1)
par(tcl = -0.5)
at.y <- 10^(-7:7)
axis(side=2, at=at.y, labels=FALSE, cex.axis=1.5, las=1)
#Map
par(xaxs = 'i', yaxs = 'i', mar = c(0,0,0,0))
par(pin = c(pin[1], ratio*pin[1]))
plot(InterpRegs, xlim = c(-82.64474, -74.5), ylim = c(36.75, 43.4))
plot(InterpRegs[which(InterpRegs$Name == "CT"),], col = 'red', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "WPA"),], col = 'orange', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "NWPANY"),], col = 'yellow', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "CNY"),], col = 'green', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "ENY"),], col = 'springgreen', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "ENYPA"),], col = 'skyblue', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "SWPA"),], col = 'blue', add = TRUE)
plot(CWV_Bounded, col = 'purple', add = TRUE)
plot(MT_Bounded, col = 'magenta', add = TRUE)
plot(VR_Bounded, col = 'gray', add = TRUE)
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(LocsNugs_Deeper1km, pch = 16, cex = 0.5, add = TRUE)
north.arrow(-75, 37, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
dev.off()

rm(at.y, dxy, lab.y, pin, ratio)

# ESDA for the well database ----
#Wells must be checked for negatives and being in the same spot before this analysis.

#  EDA Plots for all heat flow data - No Map or Histogram----
sets = rbind(c(1,2,7,7), c(3,4,7,7), c(5,6,7,7))
png("HeatFlowEDA_SplitAll_2018.png", width=9, height=9, units="in", res=600)
layout(sets)
EDAPlots = function(DataAll, Var, Unit, ymin, ymax){
  plot(DataAll$WellDepth[which(DataAll$State == "MD")], DataAll@data[Var][,1][which(DataAll$State == "MD")], col = "orange", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Maryland', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(100,1300,100),labels=FALSE)
  lines(c(-1000,10000),c(0,0))
  lines(c(1000,1000),c(-1000,1700))
  lines(c(600,600),c(-1000,1700), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "KY")], DataAll@data[Var][,1][which(DataAll$State == "KY")], col = "purple", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Kentucky', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(100,1300,100),labels=FALSE)
  lines(c(-1000,10000),c(0,0))
  lines(c(1000,1000),c(-1000,1700))
  lines(c(600,600),c(-1000,1700), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "VA")], DataAll@data[Var][,1][which(DataAll$State == "VA")], col = "yellow", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab=Unit, xlab='Depth of BHT (m)', main='Virginia', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(100,1300,100),labels=FALSE)
  lines(c(-1000,10000),c(0,0))
  lines(c(1000,1000),c(-1000,1700))
  lines(c(600,600),c(-1000,1700), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "WV")], DataAll@data[Var][,1][which(DataAll$State == "WV")], col = "green", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab=Unit, xlab='Depth of BHT (m)', main='West Virginia', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(100,1300,100),labels=FALSE)
  lines(c(-1000,10000),c(0,0))
  lines(c(1000,1000),c(-1000,1700))
  lines(c(600,600),c(-1000,1700), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "PA")], DataAll@data[Var][,1][which(DataAll$State == "PA")], col = "red", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Pennsylvania', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(100,1300,100),labels=FALSE)
  lines(c(-1000,10000),c(0,0))
  lines(c(1000,1000),c(-1000,1700))
  lines(c(600,600),c(-1000,1700), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "NY")], DataAll@data[Var][,1][which(DataAll$State == "NY")], col = "blue", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='New York', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(100,1300,100),labels=FALSE)
  lines(c(-1000,10000),c(0,0))
  lines(c(1000,1000),c(-1000,1700))
  lines(c(600,600),c(-1000,1700), lty=2)
}
par(mar=c(4,5.5,3,2))
EDAPlots(Wells_PosGrad, "Qs", Unit=expression("Surface Heat Flow" ~ (mW/m^2)),-10, 1400)
EDAPlots = function(DataAll, Var, Unit, ymin, ymax){
  plot(DataAll$WellDepth[which(DataAll$State == "NY")], DataAll@data[Var][,1][which(DataAll$State == "NY")], col = "blue", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab=Unit, xlab='Depth of BHT (m)', main='All States', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "WV")], DataAll@data[Var][,1][which(DataAll$State == "WV")], col = "green", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA")], DataAll@data[Var][,1][which(DataAll$State == "PA")], col = "red", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "KY")], DataAll@data[Var][,1][which(DataAll$State == "KY")], col = "purple", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "VA")], DataAll@data[Var][,1][which(DataAll$State == "VA")], col = "yellow", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "MD")], DataAll@data[Var][,1][which(DataAll$State == "MD")], col = "orange", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  lines(c(-1000,10000),c(0,0))
}
EDAPlots(Wells_PosGrad, "Qs", Unit=expression("Surface Heat Flow" ~ (mW/m^2)),-10, 1400)
axis(side=2, at=seq(100,1300,100),labels=FALSE)
lines(c(1000,1000),c(-1000,1700))
lines(c(600,600),c(-1000,1700), lty=2)
legend('topright', legend=c("New York", "Pennsylvania", "West Virginia", "Kentucky", "Maryland", "Virginia"), pch=16, col=c("blue", "red", "green", "purple","orange","yellow"), cex=2)
dev.off()
rm(EDAPlots)

#  ESDA Procedure #1: Local Median and Local Average Deviation----
#Uses the well data that has been sorted for unique spatial locations.

#Fixme: this function is slow, even running in parallel.
cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
DataTab = QsDev(Data = SortData$Sorted@data, Var = 'Qs', xName = 'coords_x1', yName = 'coords_x2', rad = 10000, max_pts = 25)
#Add to spatial dataframe
SortData$Sorted@data = DataTab
#Rename to another variable
WellsSort = SortData$Sorted

DataTab_NoDev = QsDev(Data = SortData_NoDev$Sorted@data, Var = 'Qs', xName = 'coords_x1', yName = 'coords_x2', rad = 10000, max_pts = 25)
SortData_NoDev$Sorted@data = DataTab_NoDev
WellsSort_NoDev = SortData_NoDev$Sorted

DataTab_NoDevNY = QsDev(Data = SortData_NoDevNY$Sorted@data, Var = 'Qs', xName = 'coords_x1', yName = 'coords_x2', rad = 10000, max_pts = 25)
SortData_NoDevNY$Sorted@data = DataTab_NoDevNY
WellsSort_NoDevNY = SortData_NoDevNY$Sorted

#Local Median Deviation - wells >= 600 m depth only
DataTab_600 = QsDev(Data = SortData$Sorted@data[which(SortData$Sorted@data$WellDepth >= 600),], Var = 'Qs', xName = 'coords_x1', yName = 'coords_x2', rad = 10000, max_pts = 25)
SortData$Sorted@data[which(SortData$Sorted@data$WellDepth >= 600),] = DataTab_600[,-c((ncol(DataTab_600)-7):(ncol(DataTab_600)-4))]
WellsSort_600 = SortData$Sorted[which(SortData$Sorted@data$WellDepth >= 600),]

#Local Median Deviation - wells >= 1000 m depth only
DataTab_1000 = QsDev(Data = SortData$Sorted@data[which(SortData$Sorted@data$WellDepth >= 1000),], Var = 'Qs', xName = 'coords_x1', yName = 'coords_x2', rad = 10000, max_pts = 25)
SortData$Sorted@data[which(SortData$Sorted@data$WellDepth >= 1000),] = DataTab_1000[,-c((ncol(DataTab_1000)-7):(ncol(DataTab_1000)-4))]
WellsSort_1000 = SortData$Sorted[which(SortData$Sorted@data$WellDepth >= 1000),]

stopCluster(cl)
rm(cl, DataTab, DataTab_NoDevNY, DataTab_NoDev, DataTab_600, DataTab_1000)

#   Plot selected radius and points used in Local Median Deviation----
#the objective of this analysis is spatial coverage so that the threshold may be applied spatially
png("HeatFlowEDA_WellsUsedOrNotUsed.png", width=6, height=6, units="in", res=600)
par(mar = c(2,2.5,3,1.5), xaxs = 'i', yaxs = 'i')
plot(WellsSort, pch = 16, cex = 0.2, col = 'white', main = 'Points Tested in Local Median Deviation ESDA Procedure')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], add=TRUE, border = 'grey')
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
plot(WellsSort[WellsSort$Dist3 <= 10000,], pch = 16, cex = 0.2, add = T)
plot(WellsSort[WellsSort$Dist3 > 10000,], pch = 16, cex = 0.2, col = 'red', add = T)
legend('topleft', legend = c('Tested', 'Not Tested'), pch = 16, cex = 1.3, col = c('black', 'red'))
dev.off()

png("HeatFlowEDA_WellsUsedOrNotUsed_NoDevNY.png", width=6, height=6, units="in", res=600)
par(mar = c(2,2.5,3,1.5), xaxs = 'i', yaxs = 'i')
plot(WellsSort_NoDevNY, pch = 16, cex = 0.2, col = 'white', main = 'Points Tested in Local Median Deviation ESDA Procedure')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], add=TRUE, border = 'grey')
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
plot(WellsSort_NoDevNY[WellsSort_NoDevNY$Dist3 <= 10000,], pch = 16, cex = 0.2, add = T)
plot(WellsSort_NoDevNY[WellsSort_NoDevNY$Dist3 > 10000,], pch = 16, cex = 0.2, col = 'red', add = T)
legend('topleft', legend = c('Tested', 'Not Tested'), pch = 16, cex = 1.3, col = c('black', 'red'))
dev.off()

#Color function parameters for plotting the heat flow data map
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#   Local Median Deviation plots With map----
sets = rbind(c(1,2,7,7), c(3,4,7,7), c(5,6,8,8), c(9,9,8,8))
png("HeatFlowEDA_SplitAll_Map_Log_WellDepthSort_MedDiff_DeepWells_Box.png", width=10, height=10, units="in", res=600)
layout(sets)
EDAPlots = function(DataAll, Var, Unit, ymin, ymax){
  plot(DataAll$WellDepth[which(DataAll$State == "MD")], DataAll@data[Var][,1][which(DataAll$State == "MD")], col = "orange", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Maryland', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  par(tck = -0.03)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  par(tck = -0.015)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
  par(xpd = TRUE)
  text(x = -2000, y = 5000, expression(bold('A')), cex = 2)
  par(xpd = FALSE)
  plot(DataAll$WellDepth[which(DataAll$State == "KY")], DataAll@data[Var][,1][which(DataAll$State == "KY")], col = "purple", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Kentucky', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  par(tck = -0.03)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  par(tck = -0.015)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "VA")], DataAll@data[Var][,1][which(DataAll$State == "VA")], col = "yellow", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab=Unit, xlab='Depth of BHT (m)', main='Virginia', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  par(tck = -0.03)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  par(tck = -0.015)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "WV")], DataAll@data[Var][,1][which(DataAll$State == "WV")], col = "green", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='West Virginia', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  par(tck = -0.03)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  par(tck = -0.015)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "PA")], DataAll@data[Var][,1][which(DataAll$State == "PA")], col = "red", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Pennsylvania', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  par(tck = -0.03)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  par(tck = -0.015)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "NY")], DataAll@data[Var][,1][which(DataAll$State == "NY")], col = "blue", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='New York', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  par(tck = -0.03)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  par(tck = -0.015)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
}
par(mar=c(4,5.5,3,2), xaxs = 'i', yaxs = 'i')
EDAPlots(WellsSort, "Qs", Unit=expression("Surface Heat Flow" ~ (mW/m^2)),.5, 2000)
par(tck = NA)
# Non-Deviation from local average plot
# EDAPlots = function(DataAll, Var, Unit, ymin, ymax){
#   plot(DataAll$WellDepth[which(DataAll$State == "NY")], DataAll@data[Var][,1][which(DataAll$State == "NY")], col = "blue", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab=Unit, xlab='Depth of BHT (m)', main='All States', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes = FALSE)
#   par(new=T)
#   plot(DataAll$WellDepth[which(DataAll$State == "WV")], DataAll@data[Var][,1][which(DataAll$State == "WV")], col = "green", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, log = 'y')
#   par(new=T)
#   plot(DataAll$WellDepth[which(DataAll$State == "PA")], DataAll@data[Var][,1][which(DataAll$State == "PA")], col = "red", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, log = 'y')
#   par(new=T)
#   plot(DataAll$WellDepth[which(DataAll$State == "KY")], DataAll@data[Var][,1][which(DataAll$State == "KY")], col = "purple", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, log = 'y')
#   par(new=T)
#   plot(DataAll$WellDepth[which(DataAll$State == "VA")], DataAll@data[Var][,1][which(DataAll$State == "VA")], col = "yellow", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, log = 'y')
#   par(new=T)
#   plot(DataAll$WellDepth[which(DataAll$State == "MD")], DataAll@data[Var][,1][which(DataAll$State == "MD")], col = "orange", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, log = 'y')
# }
# EDAPlots(WellsSort, "Qs", Unit=expression("Surface Heat Flow" ~ (mW/m^2)),.5, 2000)
# axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.5)
# axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
# axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
#Local mean deviation
# EDAPlotsMean = function(DataAll, Var, Unit, ymin, ymax){
#   plot(DataAll$WellDepth[which(DataAll$State == "NY" & DataAll$RegMean > -9999)], DataAll@data[which(DataAll$State == "NY" & DataAll$RegMean > -9999), Var] - DataAll@data[which(DataAll$State == "NY" & DataAll$RegMean > -9999), 'RegMean'], col = "blue", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab=Unit, xlab='Depth of BHT (m)', main='All States', cex.main=2, cex.axis=1.5, cex.lab=1.5, axes = FALSE)
#   par(new=T)
#   plot(DataAll$WellDepth[which(DataAll$State == "WV" & DataAll$RegMean > -9999)], DataAll@data[which(DataAll$State == "WV" & DataAll$RegMean > -9999), Var] - DataAll@data[which(DataAll$State == "WV" & DataAll$RegMean > -9999), 'RegMean'], col = "green", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
#   par(new=T)
#   plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$RegMean > -9999)], DataAll@data[which(DataAll$State == "PA" & DataAll$RegMean > -9999), Var] - DataAll@data[which(DataAll$State == "PA" & DataAll$RegMean > -9999), 'RegMean'], col = "red", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
#   par(new=T)
#   plot(DataAll$WellDepth[which(DataAll$State == "KY" & DataAll$RegMean > -9999)], DataAll@data[which(DataAll$State == "KY" & DataAll$RegMean > -9999), Var] - DataAll@data[which(DataAll$State == "KY" & DataAll$RegMean > -9999), 'RegMean'], col = "purple", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
#   par(new=T)
#   plot(DataAll$WellDepth[which(DataAll$State == "VA" & DataAll$RegMean > -9999)], DataAll@data[which(DataAll$State == "VA" & DataAll$RegMean > -9999), Var] - DataAll@data[which(DataAll$State == "VA" & DataAll$RegMean > -9999), 'RegMean'], col = "yellow", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
#   par(new=T)
#   plot(DataAll$WellDepth[which(DataAll$State == "MD" & DataAll$RegMean > -9999)], DataAll@data[which(DataAll$State == "MD" & DataAll$RegMean > -9999), Var] - DataAll@data[which(DataAll$State == "MD" & DataAll$RegMean > -9999), 'RegMean'], col = "orange", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
# }
#EDAPlotsMean(WellsSort, "Qs", Unit=expression("Deviation from Local Average Surface Heat Flow" ~ (mW/m^2)),-200, 200)
#Local median deviation
EDAPlotsMed = function(DataAll, Var, Unit, ymin, ymax){
  plot(DataAll$WellDepth[which(DataAll$State == "NY" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "NY" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "NY" & DataAll$RegMed > -9999), 'RegMed'], col = "blue", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab=Unit, xlab='Depth of BHT (m)', main='Deviation from Local Median Heat Flow', cex.main=2, cex.axis=1.5, cex.lab=1.5, axes = FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "WV" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "WV" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "WV" & DataAll$RegMed > -9999), 'RegMed'], col = "green", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "PA" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "PA" & DataAll$RegMed > -9999), 'RegMed'], col = "red", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "KY" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "KY" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "KY" & DataAll$RegMed > -9999), 'RegMed'], col = "purple", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "VA" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "VA" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "VA" & DataAll$RegMed > -9999), 'RegMed'], col = "yellow", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "MD" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "MD" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "MD" & DataAll$RegMed > -9999), 'RegMed'], col = "orange", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
}
EDAPlotsMed(WellsSort, "Qs", Unit=expression(paste("Site Q"['s'], " - Local Median Q"['s'], " (mW/m"^2, ")")),-100, 600)
axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.5)
axis(side = 2, at = seq(-200,1300,50), labels=TRUE, cex.axis=1.5)
minor.tick(nx = 5, ny = 1)
lines(c(-100,10000), c(0,0))
lines(c(1000,1000),c(-1000,2000))
lines(c(600,600),c(-1000,2000), lty=2)
legend('topright', legend=c("New York", "Pennsylvania", "West Virginia", "Kentucky", "Maryland", "Virginia"), pch=16, col=c("blue", "red", "green", "purple","orange","yellow"), cex=2)
box(which = 'figure', lwd = 2)
par(xpd = TRUE)
text(x = -800, y = 625, expression(bold('C')), cex = 2)
par(xpd = FALSE)
par(mar = c(2,2.5,1.2,1), xaxs = 'i', yaxs = 'i')
plot(WellsSort, pch = 16, col = "white", cex = 0.2)
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], add=TRUE)
#plot(WellsSort[order(WellsSort$Qs),], pch = 16, col = colFun(WellsSort$Qs[order(WellsSort$Qs)]), cex = 0.5, add = TRUE)
plot(WellsSort[order(WellsSort$WellDepth, decreasing = FALSE),], pch = 16, col = colFun(WellsSort$Qs[order(WellsSort$WellDepth, decreasing = FALSE)]), cex = 0.3, add = TRUE)
#plot(WellsSort, pch = 16, col = colFun(WellsSort$Qs), cex = 0.5)
#plot(WellsSort[order(WellsSort$WellDepth),], pch = 16, col = colFun(WellsSort$Qs[order(WellsSort$WellDepth)]), cex = 0.5)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('topleft', title = expression(paste('Q'['s'], ' (mW/m'^2,')' )), legend = c('<40', '40 - 50', '50 - 60', '60 - 70', '70 - 80', '>80'), pch = 16, cex = 1.8, col = colFun(c(30, 45, 55, 65, 75, 90)))
box(which = 'figure', lwd = 2)
par(xpd = TRUE)
text(x = -83.5, y = 43.5, expression(bold('D')), cex = 2)
par(mar=c(4,5.5,3,2), xpd = FALSE)
hist(WellsSort$Qs, breaks = seq(0, 1500, 5), freq = TRUE, xlim = c(0,120), ylim = c(0,5000), xlab = expression(paste('Surface Heat Flow (mW/m'^2,')')), cex.lab = 1.5, cex.axis = 1.5, main = 'Histogram of All Data', cex.main = 2)
box(which = 'figure', lwd = 2)
par(xpd = TRUE)
text(x = -15, y = 5500, expression(bold('B')), cex = 2)
dev.off()

#Only Local Median Deviation
png("HeatFlow_LocMedDiff.png", width=8, height=8, units="in", res=600)
par(mar=c(4,5.5,3,2), xaxs = 'i', yaxs = 'i')
EDAPlotsMed = function(DataAll, Var, Unit, ymin, ymax){
  plot(DataAll$WellDepth[which(DataAll$State == "NY" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "NY" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "NY" & DataAll$RegMed > -9999), 'RegMed'], col = "blue", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab=Unit, xlab='Temperature Measurement Depth (m)', main='Local Median Deviation by Depth', cex.main=2, cex.axis=1.5, cex.lab=2, axes = FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "WV" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "WV" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "WV" & DataAll$RegMed > -9999), 'RegMed'], col = "green", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "PA" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "PA" & DataAll$RegMed > -9999), 'RegMed'], col = "red", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "KY" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "KY" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "KY" & DataAll$RegMed > -9999), 'RegMed'], col = "purple", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "VA" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "VA" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "VA" & DataAll$RegMed > -9999), 'RegMed'], col = "yellow", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "MD" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "MD" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "MD" & DataAll$RegMed > -9999), 'RegMed'], col = "orange", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
}
EDAPlotsMed(WellsSort, "Qs", Unit=expression(paste("Local Median Deviation", "(mW/m"^2, ")")),-100, 600)
axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.5)
axis(side = 2, at = seq(-200,1300,50), labels=TRUE, cex.axis=1.5)
minor.tick(nx = 5, ny = 10)
lines(c(-100,10000), c(0,0), lwd = 2)
lines(c(1000,1000),c(-1000,2000), lwd = 3)
lines(c(600,600),c(-1000,2000), lty=2, lwd = 3)
legend('topright', legend=c("New York", "Pennsylvania", "West Virginia", "Kentucky", "Maryland", "Virginia", "600 m", "1000 m"), pch=c(rep(16, 6),NA,NA), lty = c(rep(NA, 6), 2,1), lwd = 3, col=c("blue", "red", "green", "purple","orange","yellow","black", "black"), cex=2)
dev.off()

png("HeatFlow_LocMedDiff_Zoom.png", width=8, height=8, units="in", res=600)
par(mar=c(4,5.5,3,2), xaxs = 'i', yaxs = 'i')
EDAPlotsMed = function(DataAll, Var, Unit, ymin, ymax){
  plot(DataAll$WellDepth[which(DataAll$State == "NY" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "NY" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "NY" & DataAll$RegMed > -9999), 'RegMed'], col = "blue", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab=Unit, xlab='Temperature Measurement Depth (m)', main='Local Median Deviation by Depth', cex.main=2, cex.axis=1.5, cex.lab=2, axes = FALSE, cex = 0.2)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "WV" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "WV" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "WV" & DataAll$RegMed > -9999), 'RegMed'], col = "green", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, cex = 0.2)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "PA" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "PA" & DataAll$RegMed > -9999), 'RegMed'], col = "red", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, cex = 0.2)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "KY" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "KY" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "KY" & DataAll$RegMed > -9999), 'RegMed'], col = "purple", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, cex = 0.2)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "VA" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "VA" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "VA" & DataAll$RegMed > -9999), 'RegMed'], col = "yellow", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, cex = 0.2)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "MD" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "MD" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "MD" & DataAll$RegMed > -9999), 'RegMed'], col = "orange", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, cex = 0.2)
}
EDAPlotsMed(WellsSort, "Qs", Unit=expression(paste("Local Median Deviation", "(mW/m"^2, ")")),-100, 150)
axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.5)
axis(side = 2, at = seq(-200,1300,50), labels=TRUE, cex.axis=1.5)
minor.tick(nx = 5, ny = 10)
lines(c(-100,10000), c(0,0), lwd = 2)
lines(c(1000,1000),c(-1000,2000), lwd = 3)
lines(c(600,600),c(-1000,2000), lty=2, lwd = 3)
legend('topright', legend=c("New York", "Pennsylvania", "West Virginia", "Kentucky", "Maryland", "Virginia", "600 m", "1000 m"), pch=c(rep(16, 6),NA,NA), lty = c(rep(NA, 6), 2,1), lwd = 3, col=c("blue", "red", "green", "purple","orange","yellow","black", "black"), cex=1.7)
dev.off()

png("HeatFlow_LocMedDiff_Zoom_NoDevNY.png", width=8, height=8, units="in", res=600)
par(mar=c(4,5.5,3,2), xaxs = 'i', yaxs = 'i')
EDAPlotsMed(WellsSort_NoDevNY, "Qs", Unit=expression(paste("Local Median Deviation", "(mW/m"^2, ")")),-100, 150)
axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.5)
axis(side = 2, at = seq(-200,1300,50), labels=TRUE, cex.axis=1.5)
minor.tick(nx = 5, ny = 10)
lines(c(-100,10000), c(0,0), lwd = 2)
lines(c(1000,1000),c(-1000,2000), lwd = 3)
lines(c(600,600),c(-1000,2000), lty=2, lwd = 3)
legend('topright', legend=c("New York", "Pennsylvania", "West Virginia", "Kentucky", "Maryland", "Virginia", "600 m", "1000 m"), pch=c(rep(16, 6),NA,NA), lty = c(rep(NA, 6), 2,1), lwd = 3, col=c("blue", "red", "green", "purple","orange","yellow","black", "black"), cex=1.7)
dev.off()

png("HeatFlow_LocMedDiff_600.png", width=8, height=8, units="in", res=600)
par(mar=c(4,5.5,3,2), xaxs = 'i', yaxs = 'i')
EDAPlotsMed = function(DataAll, Var, Unit, ymin, ymax){
  plot(DataAll$WellDepth[which(DataAll$State == "NY" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "NY" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "NY" & DataAll$RegMed > -9999), 'RegMed'], col = "blue", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab=Unit, xlab='Temperature Measurement Depth (m)', main='Local Median Deviation by Depth', cex.main=2, cex.axis=1.5, cex.lab=2, axes = FALSE, cex = 0.2)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "WV" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "WV" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "WV" & DataAll$RegMed > -9999), 'RegMed'], col = "green", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, cex = 0.2)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "PA" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "PA" & DataAll$RegMed > -9999), 'RegMed'], col = "red", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, cex = 0.2)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "KY" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "KY" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "KY" & DataAll$RegMed > -9999), 'RegMed'], col = "purple", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, cex = 0.2)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "VA" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "VA" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "VA" & DataAll$RegMed > -9999), 'RegMed'], col = "yellow", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, cex = 0.2)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "MD" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "MD" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "MD" & DataAll$RegMed > -9999), 'RegMed'], col = "orange", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, cex = 0.2)
}
EDAPlotsMed(WellsSort_600, "Qs", Unit=expression(paste("Local Median Deviation", "(mW/m"^2, ")")),-100, 150)
axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.5)
axis(side = 2, at = seq(-200,1300,50), labels=TRUE, cex.axis=1.5)
minor.tick(nx = 5, ny = 10)
lines(c(-100,10000), c(0,0), lwd = 2)
lines(c(1000,1000),c(-1000,2000), lwd = 3)
lines(c(600,600),c(-1000,2000), lty=2, lwd = 3)
legend('topright', legend=c("New York", "Pennsylvania", "West Virginia", "Kentucky", "Maryland", "Virginia", "600 m", "1000 m"), pch=c(rep(16, 6),NA,NA), lty = c(rep(NA, 6), 2,1), lwd = 3, col=c("blue", "red", "green", "purple","orange","yellow","black", "black"), cex=1.7)
dev.off()

png("HeatFlow_LocMedDiff_1000.png", width=8, height=8, units="in", res=600)
par(mar=c(4,5.5,3,2), xaxs = 'i', yaxs = 'i')
EDAPlotsMed = function(DataAll, Var, Unit, ymin, ymax){
  plot(DataAll$WellDepth[which(DataAll$State == "NY" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "NY" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "NY" & DataAll$RegMed > -9999), 'RegMed'], col = "blue", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab=Unit, xlab='Temperature Measurement Depth (m)', main='Local Median Deviation by Depth', cex.main=2, cex.axis=1.5, cex.lab=2, axes = FALSE, cex = 0.2)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "WV" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "WV" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "WV" & DataAll$RegMed > -9999), 'RegMed'], col = "green", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, cex = 0.2)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "PA" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "PA" & DataAll$RegMed > -9999), 'RegMed'], col = "red", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, cex = 0.2)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "KY" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "KY" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "KY" & DataAll$RegMed > -9999), 'RegMed'], col = "purple", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, cex = 0.2)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "VA" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "VA" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "VA" & DataAll$RegMed > -9999), 'RegMed'], col = "yellow", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, cex = 0.2)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "MD" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "MD" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "MD" & DataAll$RegMed > -9999), 'RegMed'], col = "orange", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE, cex = 0.2)
}
EDAPlotsMed(WellsSort_1000, "Qs", Unit=expression(paste("Local Median Deviation", "(mW/m"^2, ")")),-100, 150)
axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.5)
axis(side = 2, at = seq(-200,1300,50), labels=TRUE, cex.axis=1.5)
minor.tick(nx = 5, ny = 10)
lines(c(-100,10000), c(0,0), lwd = 2)
lines(c(1000,1000),c(-1000,2000), lwd = 3)
lines(c(600,600),c(-1000,2000), lty=2, lwd = 3)
legend('topright', legend=c("New York", "Pennsylvania", "West Virginia", "Kentucky", "Maryland", "Virginia", "600 m", "1000 m"), pch=c(rep(16, 6),NA,NA), lty = c(rep(NA, 6), 2,1), lwd = 3, col=c("blue", "red", "green", "purple","orange","yellow","black", "black"), cex=1.7)
dev.off()

#   Changepoint detection for minimum BHT well depth----
#Sort the well database by the well depth
OrderedWells = WellsSort[order(WellsSort$WellDepth),]

#Compute the changepoint as a depth series. Remove NA values from points with too few neighbors to be tested.
cpt.mean((OrderedWells$Qs - OrderedWells$RegMed)[-which(is.na(OrderedWells$Qs - OrderedWells$RegMed))], test.stat = "CUSUM", penalty = 'None', method = 'AMOC')
#767 m is changepoint detected depth with AMOC
#853 m for absolute value of the difference. Same value obtained with variance using nonparametric cumulative sum of squares (CSS) test
cpt.var((OrderedWells$Qs - OrderedWells$RegMed)[-which(is.na(OrderedWells$Qs - OrderedWells$RegMed))], penalty = 'None', mu = 0, know.mean = TRUE, test.stat = 'Normal', method = "AMOC")
#973 m for variance using normal test with known mean of 0, which is likely good for this dataset.

#   Add shallow data back into dataset for PA region ----
WellsSort$LatDeg = WellsSort@coords[,2]
WellsSort$LngDegr = WellsSort@coords[,1]
WellsSort_NoDev$LatDeg = WellsSort_NoDev@coords[,2]
WellsSort_NoDev$LngDegr = WellsSort_NoDev@coords[,1]
WellsSort_NoDevNY$LatDeg = WellsSort_NoDevNY@coords[,2]
WellsSort_NoDevNY$LngDegr = WellsSort_NoDevNY@coords[,1]
  
#Northwestern PA - All constraints are to focus on only the region of interest, and excludes all other wells in these counties.
sets = rbind(c(1,2,8,8), c(3,4,7,7), c(5,6,7,7))
png('NWPennsylvaniaShallowerWells_DeepWells_2018.png', width=9, height=9, units="in", res=600)
layout(sets)
EDAPlots = function(DataAll, Var, Unit, ymin, ymax){
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "MC KEAN")], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "MC KEAN")], col = "red", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Mc Kean', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(25,175,50),labels=FALSE)
  lines(c(600,600),c(-1000,1700), lty=2)
  lines(c(750,750),c(-1000,1700), lty=1)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "FOREST" & DataAll$WellDepth < 1600)], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "FOREST" & DataAll$WellDepth < 1600)], col = "green", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Forest', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(25,175,50),labels=FALSE)
  lines(c(600,600),c(-1000,1700), lty=2)
  lines(c(750,750),c(-1000,1700), lty=1)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "ELK")], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "ELK")], col = "blue", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab=Unit, xlab='Depth of BHT (m)', main='Elk', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(25,175,50),labels=FALSE)
  lines(c(600,600),c(-1000,1700), lty=2)
  lines(c(750,750),c(-1000,1700), lty=1)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "WARREN" & DataAll$WellDepth < 1600)], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "WARREN" & DataAll$WellDepth < 1600)], col = "purple", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab=Unit, xlab='Depth of BHT (m)', main='Warren', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(25,175,50),labels=FALSE)
  lines(c(600,600),c(-1000,1700), lty=2)
  lines(c(750,750),c(-1000,1700), lty=1)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "JEFFERSON" & DataAll$LatDeg >= 41.16776)], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "JEFFERSON" & DataAll$LatDeg >= 41.16776)], col = "orange", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Jefferson', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(25,175,50),labels=FALSE)
  lines(c(600,600),c(-1000,1700), lty=2)
  lines(c(750,750),c(-1000,1700), lty=1)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "CLARION" & DataAll$LatDeg >= 41.3)], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "CLARION" & DataAll$LatDeg >= 41.3)], col = "yellow", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Clarion', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(25,175,50),labels=FALSE)
  lines(c(600,600),c(-1000,1700), lty=2)
  lines(c(750,750),c(-1000,1700), lty=1)
}
par(mar=c(4,5.5,3,2))
EDAPlots(WellsSort, "Qs", Unit=expression("Surface Heat Flow" ~ (mW/m^2)), ymin = 0, ymax = 200)
EDAPlots = function(DataAll, Var, Unit, ymin, ymax){
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "FOREST" & DataAll$WellDepth < 1000)], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "FOREST" & DataAll$WellDepth < 1000)], col = "green", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab=Unit, xlab='Depth of BHT (m)', main='All Counties', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "MC KEAN")], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "MC KEAN")], col = "red", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "ELK")], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "ELK")], col = "blue", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "WARREN" & DataAll$WellDepth < 1000)], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "WARREN" & DataAll$WellDepth < 1000)], col = "purple", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "JEFFERSON" & DataAll$LatDegr >= 41.16776)], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "JEFFERSON" & DataAll$LatDegr >= 41.16776)], col = "orange", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "CLARION" & DataAll$LatDegr >= 41.3)], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "CLARION" & DataAll$LatDegr >= 41.3)], col = "yellow", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
}
EDAPlots(WellsSort, "Qs", Unit=expression("Surface Heat Flow" ~ (mW/m^2)),0, 200)
axis(side=2, at=seq(25,175,50),labels=FALSE)
lines(c(600,600),c(-1000,1700), lty=2)
lines(c(750,750),c(-1000,1700), lty=1)
plot(PA)
plot(Counties[which(Counties$STATEFP == 42),], add=TRUE)
plot(Counties[which(Counties$STATEFP == 42 & Counties$NAME == "McKean"),], col='red', add=TRUE)
plot(Counties[which(Counties$STATEFP == 42 & Counties$NAME == "Forest"),], col='green', add=TRUE)
plot(Counties[which(Counties$STATEFP == 42 & Counties$NAME == "Elk"),], col='blue', add=TRUE)
plot(Counties[which(Counties$STATEFP == 42 & Counties$NAME == "Warren"),], col='purple', add=TRUE)
plot(Counties[which(Counties$STATEFP == 42 & Counties$NAME == "Jefferson"),], col='orange', add=TRUE)
plot(Counties[which(Counties$STATEFP == 42 & Counties$NAME == "Clarion"),], col='yellow', add=TRUE)
dev.off()

#Using local median deviation results
sets = rbind(c(1,2,8,8), c(3,4,7,7), c(5,6,7,7))
png('NWPennsylvaniaShallowerWells_DeepWells_LocMed_2018.png', width=9, height=9, units="in", res=600)
layout(sets)
EDAPlots = function(DataAll, Var, Unit, ymin, ymax){
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "MC KEAN")], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "MC KEAN")], col = "red", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Mc Kean', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(25,175,50),labels=FALSE)
  lines(c(600,600),c(-1000,1700), lty=2)
  lines(c(750,750),c(-1000,1700), lty=1)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "FOREST" & DataAll$WellDepth < 1600)], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "FOREST" & DataAll$WellDepth < 1600)], col = "green", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Forest', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(25,175,50),labels=FALSE)
  lines(c(600,600),c(-1000,1700), lty=2)
  lines(c(750,750),c(-1000,1700), lty=1)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "ELK")], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "ELK")], col = "blue", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab=Unit, xlab='Depth of BHT (m)', main='Elk', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(25,175,50),labels=FALSE)
  lines(c(600,600),c(-1000,1700), lty=2)
  lines(c(750,750),c(-1000,1700), lty=1)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "WARREN" & DataAll$WellDepth < 1000)], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "WARREN" & DataAll$WellDepth < 1000)], col = "purple", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab=Unit, xlab='Depth of BHT (m)', main='Warren', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(25,175,50),labels=FALSE)
  lines(c(600,600),c(-1000,1700), lty=2)
  lines(c(750,750),c(-1000,1700), lty=1)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "JEFFERSON" & DataAll$LatDeg >= 41.16776)], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "JEFFERSON" & DataAll$LatDeg >= 41.16776)], col = "orange", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Jefferson', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(25,175,50),labels=FALSE)
  lines(c(600,600),c(-1000,1700), lty=2)
  lines(c(750,750),c(-1000,1700), lty=1)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "CLARION" & DataAll$LatDeg >= 41.3)], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "CLARION" & DataAll$LatDeg >= 41.3)], col = "yellow", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Clarion', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  axis(side=2, at=seq(25,175,50),labels=FALSE)
  lines(c(600,600),c(-1000,1700), lty=2)
  lines(c(750,750),c(-1000,1700), lty=1)
}
par(mar=c(4,5.5,3,2))
EDAPlots(WellsSort, "Qs", Unit=expression("Surface Heat Flow" ~ (mW/m^2)), ymin = 0, ymax = 200)
EDAPlots = function(DataAll, Var, Unit, ymin, ymax){
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "FOREST" & DataAll$WellDepth < 1000)], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "FOREST" & DataAll$WellDepth < 1000)] - DataAll@data['RegMed'][,1][which(DataAll$State == "PA" & DataAll$County == "FOREST" & DataAll$WellDepth < 1000)], col = "green", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab=Unit, xlab='Depth of BHT (m)', main='All Counties', cex.main=2, cex.axis=1.5, cex.lab=1.5)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "MC KEAN")], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "MC KEAN")] - DataAll@data['RegMed'][,1][which(DataAll$State == "PA" & DataAll$County == "MC KEAN")], col = "red", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "ELK")], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "ELK")] - DataAll@data['RegMed'][,1][which(DataAll$State == "PA" & DataAll$County == "ELK")], col = "blue", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "WARREN" & DataAll$WellDepth < 1000)], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "WARREN" & DataAll$WellDepth < 1000)] - DataAll@data['RegMed'][,1][which(DataAll$State == "PA" & DataAll$County == "WARREN" & DataAll$WellDepth < 1000)], col = "purple", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "JEFFERSON" & DataAll$LatDegr >= 41.16776)], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "JEFFERSON" & DataAll$LatDegr >= 41.16776)] - DataAll@data['RegMed'][,1][which(DataAll$State == "PA" & DataAll$County == "JEFFERSON" & DataAll$LatDegr >= 41.16776)], col = "orange", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$County == "CLARION" & DataAll$LatDegr >= 41.3)], DataAll@data[Var][,1][which(DataAll$State == "PA" & DataAll$County == "CLARION" & DataAll$LatDegr >= 41.3)] - DataAll@data['RegMed'][,1][which(DataAll$State == "PA" & DataAll$County == "CLARION" & DataAll$LatDegr >= 41.3)], col = "yellow", pch=16, xlim=c(0,2500), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
}
EDAPlots(WellsSort, "Qs", Unit=expression(paste("Local Median Deviation ", "(mW/m"^2, ")")),-50, 200)
axis(side=2, at=seq(-25,175,50),labels=FALSE)
lines(c(600,600),c(-1000,1700), lty=2)
lines(c(750,750),c(-1000,1700), lty=1)
plot(PA)
plot(Counties[which(Counties$STATEFP == 42),], add=TRUE)
plot(Counties[which(Counties$STATEFP == 42 & Counties$NAME == "McKean"),], col='red', add=TRUE)
plot(Counties[which(Counties$STATEFP == 42 & Counties$NAME == "Forest"),], col='green', add=TRUE)
plot(Counties[which(Counties$STATEFP == 42 & Counties$NAME == "Elk"),], col='blue', add=TRUE)
plot(Counties[which(Counties$STATEFP == 42 & Counties$NAME == "Warren"),], col='purple', add=TRUE)
plot(Counties[which(Counties$STATEFP == 42 & Counties$NAME == "Jefferson"),], col='orange', add=TRUE)
plot(Counties[which(Counties$STATEFP == 42 & Counties$NAME == "Clarion"),], col='yellow', add=TRUE)
dev.off()

#Remove all wells with depths less than 1000 m
WellsDeep = WellsSort[-which(WellsSort$WellDepth < 1000),]
#Then add back the PA wells.
WellsDeep = rbind(WellsDeep, WellsSort[which(WellsSort$WellDepth > 750 & WellsSort$County == 'MC KEAN' & WellsSort$State == 'PA' & WellsSort$WellDepth < 1000),])
WellsDeep = rbind(WellsDeep, WellsSort[which(WellsSort$WellDepth > 750 & WellsSort$WellDepth < 1000 & WellsSort$County == 'ELK' & WellsSort$State == 'PA'),])
WellsDeep = rbind(WellsDeep, WellsSort[which(WellsSort$WellDepth > 750 & WellsSort$WellDepth < 1000 & WellsSort$County == 'WARREN' & WellsSort$State == 'PA' & WellsSort$LongDgr > -79.4),])
WellsDeep = rbind(WellsDeep, WellsSort[which(WellsSort$WellDepth > 750 & WellsSort$WellDepth < 1000 & WellsSort$County == 'FOREST' & WellsSort$State == 'PA'),])
WellsDeep = rbind(WellsDeep, WellsSort[which(WellsSort$WellDepth > 750 & WellsSort$WellDepth < 1000 & WellsSort$County == 'CLARION' & WellsSort$State == 'PA' & WellsSort$LatDegr > 41.3),])
WellsDeep = rbind(WellsDeep, WellsSort[which(WellsSort$WellDepth > 750 & WellsSort$WellDepth < 1000 & WellsSort$County == 'JEFFERSON' & WellsSort$State == 'PA' & WellsSort$LatDegr <= 41.16776),])

writeOGR(WellsDeep, dsn=getwd(), layer="WellsForOutlierTest_ESDA_2018", driver="ESRI Shapefile")

WellsDeep_NoDev = WellsSort_NoDev[-which(WellsSort_NoDev$WellDepth < 1000),]
#Then add back the PA wells.
WellsDeep_NoDev = rbind(WellsDeep_NoDev, WellsSort_NoDev[which(WellsSort_NoDev$WellDepth > 750 & WellsSort_NoDev$County == 'MC KEAN' & WellsSort_NoDev$State == 'PA' & WellsSort_NoDev$WellDepth < 1000),])
WellsDeep_NoDev = rbind(WellsDeep_NoDev, WellsSort_NoDev[which(WellsSort_NoDev$WellDepth > 750 & WellsSort_NoDev$WellDepth < 1000 & WellsSort_NoDev$County == 'ELK' & WellsSort_NoDev$State == 'PA'),])
WellsDeep_NoDev = rbind(WellsDeep_NoDev, WellsSort_NoDev[which(WellsSort_NoDev$WellDepth > 750 & WellsSort_NoDev$WellDepth < 1000 & WellsSort_NoDev$County == 'WARREN' & WellsSort_NoDev$State == 'PA' & WellsSort_NoDev$LongDgr > -79.4),])
WellsDeep_NoDev = rbind(WellsDeep_NoDev, WellsSort_NoDev[which(WellsSort_NoDev$WellDepth > 750 & WellsSort_NoDev$WellDepth < 1000 & WellsSort_NoDev$County == 'FOREST' & WellsSort_NoDev$State == 'PA'),])
WellsDeep_NoDev = rbind(WellsDeep_NoDev, WellsSort_NoDev[which(WellsSort_NoDev$WellDepth > 750 & WellsSort_NoDev$WellDepth < 1000 & WellsSort_NoDev$County == 'CLARION' & WellsSort_NoDev$State == 'PA' & WellsSort_NoDev$LatDegr > 41.3),])
WellsDeep_NoDev = rbind(WellsDeep_NoDev, WellsSort_NoDev[which(WellsSort_NoDev$WellDepth > 750 & WellsSort_NoDev$WellDepth < 1000 & WellsSort_NoDev$County == 'JEFFERSON' & WellsSort_NoDev$State == 'PA' & WellsSort_NoDev$LatDegr <= 41.16776),])

WellsDeep_NoDevNY = WellsSort_NoDevNY[-which(WellsSort_NoDevNY$WellDepth < 1000),]
#Then add back the PA wells.
WellsDeep_NoDevNY = rbind(WellsDeep_NoDevNY, WellsSort_NoDevNY[which(WellsSort_NoDevNY$WellDepth > 750 & WellsSort_NoDevNY$County == 'MC KEAN' & WellsSort_NoDevNY$State == 'PA' & WellsSort_NoDevNY$WellDepth < 1000),])
WellsDeep_NoDevNY = rbind(WellsDeep_NoDevNY, WellsSort_NoDevNY[which(WellsSort_NoDevNY$WellDepth > 750 & WellsSort_NoDevNY$WellDepth < 1000 & WellsSort_NoDevNY$County == 'ELK' & WellsSort_NoDevNY$State == 'PA'),])
WellsDeep_NoDevNY = rbind(WellsDeep_NoDevNY, WellsSort_NoDevNY[which(WellsSort_NoDevNY$WellDepth > 750 & WellsSort_NoDevNY$WellDepth < 1000 & WellsSort_NoDevNY$County == 'WARREN' & WellsSort_NoDevNY$State == 'PA' & WellsSort_NoDevNY$LongDgr > -79.4),])
WellsDeep_NoDevNY = rbind(WellsDeep_NoDevNY, WellsSort_NoDevNY[which(WellsSort_NoDevNY$WellDepth > 750 & WellsSort_NoDevNY$WellDepth < 1000 & WellsSort_NoDevNY$County == 'FOREST' & WellsSort_NoDevNY$State == 'PA'),])
WellsDeep_NoDevNY = rbind(WellsDeep_NoDevNY, WellsSort_NoDevNY[which(WellsSort_NoDevNY$WellDepth > 750 & WellsSort_NoDevNY$WellDepth < 1000 & WellsSort_NoDevNY$County == 'CLARION' & WellsSort_NoDevNY$State == 'PA' & WellsSort_NoDevNY$LatDegr > 41.3),])
WellsDeep_NoDevNY = rbind(WellsDeep_NoDevNY, WellsSort_NoDevNY[which(WellsSort_NoDevNY$WellDepth > 750 & WellsSort_NoDevNY$WellDepth < 1000 & WellsSort_NoDevNY$County == 'JEFFERSON' & WellsSort_NoDevNY$State == 'PA' & WellsSort_NoDevNY$LatDegr <= 41.16776),])


png("WellsRemoved_PAWellsAddedBack.png", width=6, height=6, units='in', res=150)
plot(WellsSort, pch=16, col='red')
plot(WellsDeep, pch=16, add=TRUE)
plot(NY, add=TRUE)
plot(PA, add=TRUE)
plot(WV, add=TRUE)
plot(MD, add=TRUE)
plot(KY, add=TRUE)
plot(VA, add=TRUE)
dev.off()

rm(sets, Pal, scaleBy, scaleRange)
#  ESDA Procedure #2: Operator Analysis----
# Some operators may have practices that make their data unusual. Use these diagnostics to detect them.
# If there are significant differences, see if operator should be removed.
# Careful that some operators may make up all of the data for a region.

#   Prepare the data ----
#Well data in WGS84 NAD83 is required for mapping
WellsDeepWGS = spTransform(WellsDeep, CRSobj = CRS('+init=epsg:4326'))
WellsDeepWGS_NoDev = spTransform(WellsDeep_NoDev, CRSobj = CRS('+init=epsg:4326'))
WellsDeepWGS_NoDevNY = spTransform(WellsDeep_NoDevNY, CRSobj = CRS('+init=epsg:4326'))

#    Get the well data for each geologic region.----
#Using dataset for points without negative gradients, and no points in same spatial location
DeepCT = spTransform(WellsDeepWGS[InterpRegs[InterpRegs$Name == 'CT',],], CRS('+init=epsg:26917'))
DeepCNY = spTransform(WellsDeepWGS[InterpRegs[InterpRegs$Name == 'CNY',],], CRS('+init=epsg:26917'))
DeepCWV = spTransform(WellsDeepWGS[CWV_Bounded,], CRS('+init=epsg:26917'))
DeepENY = spTransform(WellsDeepWGS[InterpRegs[InterpRegs$Name == 'ENY',],], CRS('+init=epsg:26917'))
DeepENYPA = spTransform(WellsDeepWGS[InterpRegs[InterpRegs$Name == 'ENYPA',],], CRS('+init=epsg:26917'))
DeepMT = spTransform(WellsDeepWGS[MT_Bounded,], CRS('+init=epsg:26917'))
DeepNWPANY = spTransform(WellsDeepWGS[InterpRegs[InterpRegs$Name == 'NWPANY',],], CRS('+init=epsg:26917'))
DeepSWPA = spTransform(WellsDeepWGS[InterpRegs[InterpRegs$Name == 'SWPA',],], CRS('+init=epsg:26917'))
DeepWPA = spTransform(WellsDeepWGS[InterpRegs[InterpRegs$Name == 'WPA',],], CRS('+init=epsg:26917'))
DeepVR = spTransform(WellsDeepWGS[VR_Bounded,], CRS('+init=epsg:26917'))
DeepFL = rbind(DeepCT, DeepCWV, DeepCNY, DeepENY, DeepENYPA, DeepMT, DeepNWPANY, DeepSWPA, DeepWPA, DeepVR)

#Seems like removing the possibly deviated data in NY doesn't make much of a difference.
#So, only do this analysis for the possibly deviated data in PA and WV
DeepCT_NoDevNY = spTransform(WellsDeepWGS_NoDevNY[InterpRegs[InterpRegs$Name == 'CT',],], CRS('+init=epsg:26917'))
DeepCNY_NoDevNY = spTransform(WellsDeepWGS_NoDevNY[InterpRegs[InterpRegs$Name == 'CNY',],], CRS('+init=epsg:26917'))
DeepCWV_NoDevNY = spTransform(WellsDeepWGS_NoDevNY[CWV_Bounded,], CRS('+init=epsg:26917'))
DeepENY_NoDevNY = spTransform(WellsDeepWGS_NoDevNY[InterpRegs[InterpRegs$Name == 'ENY',],], CRS('+init=epsg:26917'))
DeepENYPA_NoDevNY = spTransform(WellsDeepWGS_NoDevNY[InterpRegs[InterpRegs$Name == 'ENYPA',],], CRS('+init=epsg:26917'))
DeepMT_NoDevNY = spTransform(WellsDeepWGS_NoDevNY[MT_Bounded,], CRS('+init=epsg:26917'))
DeepNWPANY_NoDevNY = spTransform(WellsDeepWGS_NoDevNY[InterpRegs[InterpRegs$Name == 'NWPANY',],], CRS('+init=epsg:26917'))
DeepSWPA_NoDevNY = spTransform(WellsDeepWGS_NoDevNY[InterpRegs[InterpRegs$Name == 'SWPA',],], CRS('+init=epsg:26917'))
DeepWPA_NoDevNY = spTransform(WellsDeepWGS_NoDevNY[InterpRegs[InterpRegs$Name == 'WPA',],], CRS('+init=epsg:26917'))
DeepVR_NoDevNY = spTransform(WellsDeepWGS_NoDevNY[VR_Bounded,], CRS('+init=epsg:26917'))
DeepFL_NoDevNY = rbind(DeepCT_NoDevNY, DeepCWV_NoDevNY, DeepCNY_NoDevNY, DeepENY_NoDevNY, DeepENYPA_NoDevNY, DeepMT_NoDevNY, DeepNWPANY_NoDevNY, DeepSWPA_NoDevNY, DeepWPA_NoDevNY, DeepVR_NoDevNY)

#Project to WGS - used for mapping in the function
DeepCT_WGS = spTransform(DeepCT, CRS('+init=epsg:4326'))
DeepCT_WGS_NoDevNY = spTransform(DeepCT_NoDevNY, CRS('+init=epsg:4326'))
DeepCNY_WGS = spTransform(DeepCNY, CRS('+init=epsg:4326'))
DeepCNY_WGS_NoDevNY = spTransform(DeepCNY_NoDevNY, CRS('+init=epsg:4326'))
DeepCWV_WGS = spTransform(DeepCWV, CRS('+init=epsg:4326'))
DeepCWV_WGS_NoDevNY = spTransform(DeepCWV_NoDevNY, CRS('+init=epsg:4326'))
DeepENY_WGS = spTransform(DeepENY, CRS('+init=epsg:4326'))
DeepENY_WGS_NoDevNY = spTransform(DeepENY_NoDevNY, CRS('+init=epsg:4326'))
DeepENYPA_WGS = spTransform(DeepENYPA, CRS('+init=epsg:4326'))
DeepENYPA_WGS_NoDevNY = spTransform(DeepENYPA_NoDevNY, CRS('+init=epsg:4326'))
DeepMT_WGS = spTransform(DeepMT, CRS('+init=epsg:4326'))
DeepMT_WGS_NoDevNY = spTransform(DeepMT_NoDevNY, CRS('+init=epsg:4326'))
DeepNWPANY_WGS = spTransform(DeepNWPANY, CRS('+init=epsg:4326'))
DeepNWPANY_WGS_NoDevNY = spTransform(DeepNWPANY_NoDevNY, CRS('+init=epsg:4326'))
DeepSWPA_WGS = spTransform(DeepSWPA, CRS('+init=epsg:4326'))
DeepSWPA_WGS_NoDevNY = spTransform(DeepSWPA_NoDevNY, CRS('+init=epsg:4326'))
DeepWPA_WGS = spTransform(DeepWPA, CRS('+init=epsg:4326'))
DeepWPA_WGS_NoDevNY = spTransform(DeepWPA_NoDevNY, CRS('+init=epsg:4326'))
DeepVR_WGS = spTransform(DeepVR, CRS('+init=epsg:4326'))
DeepVR_WGS_NoDevNY = spTransform(DeepVR_NoDevNY, CRS('+init=epsg:4326'))

#    Variograms for deep data using MOM and Cressie's robust estimator----
v.DeepCT <- variogram(Qs~1, DeepCT, cutoff=60000, width=60000/50) 
v.DeepCNY <- variogram(Qs~1, DeepCNY, cutoff=60000, width=60000/15) 
v.DeepCWV <- variogram(Qs~1, DeepCWV, cutoff=60000, width=60000/50)
v.DeepENY <- variogram(Qs~1, DeepENY, cutoff=60000, width=60000/15)
v.DeepENYPA <- variogram(Qs~1, DeepENYPA, cutoff=60000, width=60000/40)
v.DeepMT <- variogram(Qs~1, DeepMT, cutoff=60000, width=60000/50)
v.DeepNWPANY <- variogram(Qs~1, DeepNWPANY, cutoff=60000, width=60000/20) 
v.DeepSWPA <- variogram(Qs~1, DeepSWPA, cutoff=60000, width=60000/50) 
v.DeepWPA <- variogram(Qs~1, DeepWPA, cutoff=60000, width=60000/50) 
v.DeepVR <- variogram(Qs~1, DeepVR, cutoff=60000, width=60000/20) 
v.DeepFL <- variogram(Qs~1, DeepFL, cutoff=60000, width=60000/200)

v.DeepCT_NoDevNY <- variogram(Qs~1, DeepCT_NoDevNY, cutoff=60000, width=60000/50) 
v.DeepCNY_NoDevNY <- variogram(Qs~1, DeepCNY_NoDevNY, cutoff=60000, width=60000/15) 
v.DeepCWV_NoDevNY <- variogram(Qs~1, DeepCWV_NoDevNY, cutoff=60000, width=60000/50)
v.DeepENY_NoDevNY <- variogram(Qs~1, DeepENY_NoDevNY, cutoff=60000, width=60000/15)
v.DeepENYPA_NoDevNY <- variogram(Qs~1, DeepENYPA_NoDevNY, cutoff=60000, width=60000/40)
v.DeepMT_NoDevNY <- variogram(Qs~1, DeepMT_NoDevNY, cutoff=60000, width=60000/50)
v.DeepNWPANY_NoDevNY <- variogram(Qs~1, DeepNWPANY_NoDevNY, cutoff=60000, width=60000/20) 
v.DeepSWPA_NoDevNY <- variogram(Qs~1, DeepSWPA_NoDevNY, cutoff=60000, width=60000/50) 
v.DeepWPA_NoDevNY <- variogram(Qs~1, DeepWPA_NoDevNY, cutoff=60000, width=60000/50) 
v.DeepVR_NoDevNY <- variogram(Qs~1, DeepVR_NoDevNY, cutoff=60000, width=60000/20) 
v.DeepFL_NoDevNY <- variogram(Qs~1, DeepFL_NoDevNY, cutoff=60000, width=60000/200)

rv.DeepCT <- variogram(Qs~1, DeepCT, cutoff=60000, width=60000/50, cressie = TRUE) 
rv.DeepCNY <- variogram(Qs~1, DeepCNY, cutoff=60000, width=60000/15, cressie = TRUE) 
rv.DeepCWV <- variogram(Qs~1, DeepCWV, cutoff=60000, width=60000/50, cressie = TRUE)
rv.DeepENY <- variogram(Qs~1, DeepENY, cutoff=60000, width=60000/15, cressie = TRUE)
rv.DeepENYPA <- variogram(Qs~1, DeepENYPA, cutoff=60000, width=60000/40, cressie = TRUE)
rv.DeepMT <- variogram(Qs~1, DeepMT, cutoff=60000, width=60000/50, cressie = TRUE)
rv.DeepNWPANY <- variogram(Qs~1, DeepNWPANY, cutoff=60000, width=60000/20, cressie = TRUE) 
rv.DeepSWPA <- variogram(Qs~1, DeepSWPA, cutoff=60000, width=60000/50, cressie = TRUE) 
rv.DeepWPA <- variogram(Qs~1, DeepWPA, cutoff=60000, width=60000/50, cressie = TRUE) 
rv.DeepVR <- variogram(Qs~1, DeepVR, cutoff=60000, width=60000/20, cressie = TRUE) 
rv.DeepFL <- variogram(Qs~1, DeepFL, cutoff=60000, width=60000/200, cressie = TRUE)

rv.DeepCT_NoDevNY <- variogram(Qs~1, DeepCT_NoDevNY, cutoff=60000, width=60000/50, cressie = TRUE) 
rv.DeepCNY_NoDevNY <- variogram(Qs~1, DeepCNY_NoDevNY, cutoff=60000, width=60000/15, cressie = TRUE) 
rv.DeepCWV_NoDevNY <- variogram(Qs~1, DeepCWV_NoDevNY, cutoff=60000, width=60000/50, cressie = TRUE)
rv.DeepENY_NoDevNY <- variogram(Qs~1, DeepENY_NoDevNY, cutoff=60000, width=60000/15, cressie = TRUE)
rv.DeepENYPA_NoDevNY <- variogram(Qs~1, DeepENYPA_NoDevNY, cutoff=60000, width=60000/40, cressie = TRUE)
rv.DeepMT_NoDevNY <- variogram(Qs~1, DeepMT_NoDevNY, cutoff=60000, width=60000/50, cressie = TRUE)
rv.DeepNWPANY_NoDevNY <- variogram(Qs~1, DeepNWPANY_NoDevNY, cutoff=60000, width=60000/20, cressie = TRUE) 
rv.DeepSWPA_NoDevNY <- variogram(Qs~1, DeepSWPA_NoDevNY, cutoff=60000, width=60000/50, cressie = TRUE) 
rv.DeepWPA_NoDevNY <- variogram(Qs~1, DeepWPA_NoDevNY, cutoff=60000, width=60000/50, cressie = TRUE) 
rv.DeepVR_NoDevNY <- variogram(Qs~1, DeepVR_NoDevNY, cutoff=60000, width=60000/20, cressie = TRUE) 
rv.DeepFL_NoDevNY <- variogram(Qs~1, DeepFL_NoDevNY, cutoff=60000, width=60000/200, cressie = TRUE)


#    Set the colors for heat flow on maps----
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)
#   Run bad operator diagnostics for each geologic region----

#    MT---- 
#Deviated wells do not look incorrect here.
cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
OpDiag_DeepMT = foreach(o = 1:length(unique(DeepMT$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepMT, o = o, LowLim = 2, MT_WGS = DeepMT_WGS, v.MT = v.DeepMT, rv.MT = rv.DeepMT, Vcut = 60000, Vbins = 50, HistSep = 10, Histylim = 500, RegName = "DeepMT", MaxLagDist = 60000, SensitivityDist = 2500, plt = TRUE)
  a
}
OpDiag_DeepMT_NoDevNY = foreach(o = 1:length(unique(DeepMT_NoDevNY$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepMT_NoDevNY, o = o, LowLim = 2, MT_WGS = DeepMT_WGS_NoDevNY, v.MT = v.DeepMT_NoDevNY, rv.MT = rv.DeepMT_NoDevNY, Vcut = 60000, Vbins = 50, HistSep = 10, Histylim = 500, RegName = "DeepMT_NoDevNY", MaxLagDist = 60000, SensitivityDist = 2500, plt = TRUE)
  a
}
stopCluster(cl)

#p-value colors
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(0,0.2)
scaleBy = 0.05
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#Sensitivity Plot for All Metrics
OpRanks_DeepMT = GetOperatorRanks(OpDiag_DeepMT)
PlotDistanceSensitivity(OpDiagRanks = OpRanks_DeepMT, res = 300, height = 10, width = 7, xlim = c(1.5,60), NumOps = 100, NumOpsStep = 20, xStep = 2.5, PlotName = 'DeepMT')
#Plot operator diagnostics for each of the distances
cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
temp = foreach (i = 1:ncol(OpRanks_DeepMT$DiffMOM), .packages = 'RColorBrewer') %dopar% {
  PlotOpsDiagnostics(OpDiag = OpDiag_DeepMT, OpDiagRanks = OpRanks_DeepMT, col = i, res = 300, NumOps = 100, NumOpsStep = 20, PlotName = paste0('DeepMT_',i), xStep = 2.5)
}
rm(temp)

OpRanks_DeepMT_NoDevNY = GetOperatorRanks(OpDiag_DeepMT_NoDevNY)
PlotDistanceSensitivity(OpDiagRanks = OpRanks_DeepMT_NoDevNY, res = 300, height = 10, width = 7, xlim = c(1.5,60), NumOps = 100, NumOpsStep = 20, xStep = 2.5, PlotName = 'DeepMT_NoDevNY')
temp = foreach (i = 1:ncol(OpRanks_DeepMT_NoDevNY$DiffMOM), .packages = 'RColorBrewer') %dopar% {
  PlotOpsDiagnostics(OpDiag = OpDiag_DeepMT_NoDevNY, OpDiagRanks = OpRanks_DeepMT_NoDevNY, col = i, res = 300, NumOps = 100, NumOpsStep = 20, PlotName = paste0('DeepMT_NoDevNY_',i), xStep = 2.5)
}
stopCluster(cl)
rm(temp)

#     Remove Waco and Equitable Production Company----
#Waco is a clear bad operator because they logged wells upwards. Equitable Production Company had wells logged by Waco. Remove both.
DeepMT_NoWaco = DeepMT[-which(DeepMT$Operator == 'Waco Oil & Gas Co., Inc.'),]
DeepMT_NoWaco = DeepMT_NoWaco[-which(DeepMT_NoWaco$Operator == 'Equitable Production Company'),]
DeepMT_WGS = spTransform(DeepMT_NoWaco, CRS('+init=epsg:4326'))
v.DeepMT_NoWaco = variogram(Qs ~ 1, DeepMT_NoWaco, cutoff = 60000, width = 60000/50)
rv.DeepMT_NoWaco = variogram(Qs ~ 1, DeepMT_NoWaco, cutoff = 60000, width = 60000/50, cressie = TRUE)

colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
OpDiag_DeepMT_NoWaco = foreach(o = 1:length(unique(DeepMT_NoWaco$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepMT_NoWaco, o = o, LowLim = 2, MT_WGS = DeepMT_WGS, v.MT = v.DeepMT_NoWaco, rv.MT = rv.DeepMT_NoWaco, Vcut = 60000, Vbins = 50, HistSep = 10, Histylim = 300, RegName = "DeepMT_NoWaco", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
#p-value colors
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(0,0.2)
scaleBy = 0.05
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#Sensitivity Plot for All Metrics
OpRanks_DeepMT_NoWaco = GetOperatorRanks(OpDiag_DeepMT_NoWaco)
PlotDistanceSensitivity(OpDiagRanks = OpRanks_DeepMT_NoWaco, res = 300, height = 10, width = 7, xlim = c(1.5,60), NumOps = 100, NumOpsStep = 20, xStep = 2.5, PlotName = 'DeepMT_NoWaco')
#Plot operator diagnostics for each of the distances
temp = foreach (i = 1:ncol(OpRanks_DeepMT_NoWaco$DiffMOM), .packages = 'RColorBrewer') %dopar% {
  PlotOpsDiagnostics(OpDiag = OpDiag_DeepMT_NoWaco, OpDiagRanks = OpRanks_DeepMT_NoWaco, col = i, res = 300, NumOps = 100, NumOpsStep = 20, PlotName = paste0('DeepMT_NoWaco',i*2500), xStep = 2.5)
}
rm(temp)
stopCluster(cl)

#Top 10 among all metrics, without p-value - These all look reasonable
TopBadOps_MT = ReturnTopN(OpDiag = OpDiag_DeepMT_NoWaco, OpRanks = OpRanks_DeepMT_NoWaco, col = 6, N = 5, ID = TRUE, WellData = DeepMT_NoWaco)
TopBadOps_MT_Names = ReturnTopN(OpDiag = OpDiag_DeepMT_NoWaco, OpRanks = OpRanks_DeepMT_NoWaco, col = 6, N = 5, ID = FALSE, WellData = DeepMT_NoWaco)

#    CWV ---- 
#Deviated wells do not look incorrect here.
#Remove the 2 high points that were checked for being incorrect.
DeepCWV_DelHigh = DeepCWV[DeepCWV$Qs < 160,]
DeepCWV_NoDevNY_DelHigh = DeepCWV_NoDevNY[DeepCWV_NoDevNY$Qs < 160,]
#Remove Waco
DeepCWV_NoWaco = DeepCWV_DelHigh[-which(DeepCWV_DelHigh$Operator == "Waco Oil & Gas Co., Inc."),]
v.DeepCWV_NoWaco = variogram(Qs~1, DeepCWV_NoWaco, cutoff = 60000, width = 60000/50)
rv.DeepCWV_NoWaco = variogram(Qs~1, DeepCWV_NoWaco, cutoff = 60000, width = 60000/50, cressie = TRUE)
DeepCWV_WGS = spTransform(DeepCWV_NoWaco, CRS('+init=epsg:4326'))

DeepCWV_NoDevNY_NoWaco = DeepCWV_NoDevNY_DelHigh[-which(DeepCWV_NoDevNY_DelHigh$Operator == "Waco Oil & Gas Co., Inc."),]
v.DeepCWV_NoDevNY_NoWaco = variogram(Qs~1, DeepCWV_NoDevNY_NoWaco, cutoff = 60000, width = 60000/50)
rv.DeepCWV_NoDevNY_NoWaco = variogram(Qs~1, DeepCWV_NoDevNY_NoWaco, cutoff = 60000, width = 60000/50, cressie = TRUE)
DeepCWV_WGS_NoDevNY = spTransform(DeepCWV_NoDevNY_NoWaco, CRS('+init=epsg:4326'))

colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
OpDiag_DeepCWV = foreach(o = 1:length(unique(DeepCWV_NoWaco$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepCWV_NoWaco, o = o, LowLim = 2, MT_WGS = DeepCWV_WGS, v.MT = v.DeepCWV_NoWaco, rv.MT = rv.DeepCWV_NoWaco, Vcut = 60000, Vbins = 50, HistSep = 10, Histylim = 1300, RegName = "DeepCWV_NoWaco", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
OpDiag_DeepCWV_NoDevNY = foreach(o = 1:length(unique(DeepCWV_NoDevNY_NoWaco$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepCWV_NoDevNY_NoWaco, o = o, LowLim = 2, MT_WGS = DeepCWV_WGS_NoDevNY, v.MT = v.DeepCWV_NoDevNY_NoWaco, rv.MT = rv.DeepCWV_NoDevNY_NoWaco, Vcut = 60000, Vbins = 50, HistSep = 10, Histylim = 1300, RegName = "DeepCWV_NoDevNY_NoWaco", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
#p-value colors
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(0,0.2)
scaleBy = 0.05
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#Sensitivity Plot for All Metrics
OpRanks_DeepCWV = GetOperatorRanks(OpDiag_DeepCWV)
PlotDistanceSensitivity(OpDiagRanks = OpRanks_DeepCWV, res = 300, height = 10, width = 7, xlim = c(1.5,60), NumOps = 220, NumOpsStep = 20, xStep = 2.5, PlotName = 'DeepCWV_NoWaco')
#Plot operator diagnostics for each of the distances
temp = foreach (i = 1:ncol(OpRanks_DeepCWV$DiffMOM), .packages = 'RColorBrewer') %dopar% {
  PlotOpsDiagnostics(OpDiag = OpDiag_DeepCWV, OpDiagRanks = OpRanks_DeepCWV, col = i, res = 300, NumOps = 220, NumOpsStep = 20, PlotName = paste0('DeepCWV_NoWaco',i*2500), xStep = 2.5)
}
rm(temp)
stopCluster(cl)

#     Split into KY and WV clusters ---- 
#KY and WV have similar variogram shapes. WV more variable, but fewer observations.
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

DeepCWV_NoWaco_KY = DeepCWV_NoWaco[-which(DeepCWV_NoWaco$LatDeg >= 38.5),]
v.DeepCWV_NoWaco_KY = variogram(Qs~1, DeepCWV_NoWaco_KY, cutoff = 60000, width = 60000/50)
rv.DeepCWV_NoWaco_KY = variogram(Qs~1, DeepCWV_NoWaco_KY, cutoff = 60000, width = 60000/50, cressie = TRUE)
DeepCWV_WGS_KY = spTransform(DeepCWV_NoWaco_KY, CRS('+init=epsg:4326'))

DeepCWV_NoWaco_WV = DeepCWV_NoWaco[-which(DeepCWV_NoWaco$LatDeg < 38.5),]
v.DeepCWV_NoWaco_WV = variogram(Qs~1, DeepCWV_NoWaco_WV, cutoff = 60000, width = 60000/50)
rv.DeepCWV_NoWaco_WV = variogram(Qs~1, DeepCWV_NoWaco_WV, cutoff = 60000, width = 60000/50, cressie = TRUE)
DeepCWV_WGS_WV = spTransform(DeepCWV_NoWaco_WV, CRS('+init=epsg:4326'))

cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
OpDiag_DeepCWV_NoWaco_KY = foreach(o = 1:length(unique(DeepCWV_NoWaco_KY$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepCWV_NoWaco_KY, o = o, LowLim = 2, MT_WGS = DeepCWV_WGS_KY, v.MT = v.DeepCWV_NoWaco_KY, rv.MT = rv.DeepCWV_NoWaco_KY, Vcut = 60000, Vbins = 50, HistSep = 10, Histylim = 1300, RegName = "DeepCWV_NoWaco_KY", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
OpDiag_DeepCWV_NoWaco_WV = foreach(o = 1:length(unique(DeepCWV_NoWaco_WV$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepCWV_NoWaco_WV, o = o, LowLim = 2, MT_WGS = DeepCWV_WGS_WV, v.MT = v.DeepCWV_NoWaco_WV, rv.MT = rv.DeepCWV_NoWaco_WV, Vcut = 60000, Vbins = 50, HistSep = 10, Histylim = 1300, RegName = "DeepCWV_NoWaco_WV", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
#p-value colors
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(0,0.2)
scaleBy = 0.05
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#Sensitivity Plot for All Metrics
OpRanks_DeepCWV_KY = GetOperatorRanks(OpDiag_DeepCWV_NoWaco_KY)
OpRanks_DeepCWV_WV = GetOperatorRanks(OpDiag_DeepCWV_NoWaco_WV)
PlotDistanceSensitivity(OpDiagRanks = OpRanks_DeepCWV_KY, res = 300, height = 10, width = 7, xlim = c(1.5,60), NumOps = 100, NumOpsStep = 20, xStep = 2.5, PlotName = 'DeepCWV_NoWaco_KY')
PlotDistanceSensitivity(OpDiagRanks = OpRanks_DeepCWV_WV, res = 300, height = 10, width = 7, xlim = c(1.5,60), NumOps = 120, NumOpsStep = 20, xStep = 2.5, PlotName = 'DeepCWV_NoWaco_WV')
#Plot operator diagnostics for each of the distances
temp = foreach (i = 1:ncol(OpRanks_DeepCWV_KY$DiffMOM), .packages = 'RColorBrewer') %dopar% {
  PlotOpsDiagnostics(OpDiag = OpDiag_DeepCWV_NoWaco_KY, OpDiagRanks = OpRanks_DeepCWV_KY, col = i, res = 300, NumOps = 100, NumOpsStep = 20, PlotName = paste0('DeepCWV_NoWaco_KY',i*2500), xStep = 2.5)
}
temp = foreach (i = 1:ncol(OpRanks_DeepCWV_WV$DiffMOM), .packages = 'RColorBrewer') %dopar% {
  PlotOpsDiagnostics(OpDiag = OpDiag_DeepCWV_NoWaco_WV, OpDiagRanks = OpRanks_DeepCWV_WV, col = i, res = 300, NumOps = 120, NumOpsStep = 20, PlotName = paste0('DeepCWV_NoWaco_WV',i*2500), xStep = 2.5)
}
rm(temp)
stopCluster(cl)

#Petroleum Resources Inc has some high and low points, but most cluster around the same value.
#Top 10 among all metrics, without p-value
TopBadOps_CWV_KY = ReturnTopN(OpDiag = OpDiag_DeepCWV_NoWaco_KY, OpRanks = OpRanks_DeepCWV_KY, col = 6, N = 5, ID = TRUE, WellData = DeepCWV_NoWaco_KY)
TopBadOps_CWV_KY_Names = ReturnTopN(OpDiag = OpDiag_DeepCWV_NoWaco_KY, OpRanks = OpRanks_DeepCWV_KY, col = 6, N = 5, ID = FALSE, WellData = DeepCWV_NoWaco_KY)
TopBadOps_CWV_WV = ReturnTopN(OpDiag = OpDiag_DeepCWV_NoWaco_WV, OpRanks = OpRanks_DeepCWV_WV, col = 6, N = 5, ID = TRUE, WellData = DeepCWV_NoWaco_WV)
TopBadOps_CWV_WV_Names = ReturnTopN(OpDiag = OpDiag_DeepCWV_NoWaco_WV, OpRanks = OpRanks_DeepCWV_WV, col = 6, N = 5, ID = FALSE, WellData = DeepCWV_NoWaco_WV)

#    SWPA----
#Deviated wells do not look incorrect here.
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
OpDiag_DeepSWPA = foreach(o = 1:length(unique(DeepSWPA$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepSWPA, o = o, LowLim = 2, MT_WGS = DeepSWPA_WGS, v.MT = v.DeepSWPA, rv.MT = rv.DeepSWPA, Vcut = 60000, Vbins = 50, HistSep = 10, Histylim = 1300, RegName = "DeepSWPA", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
OpDiag_DeepSWPA_NoDevNY = foreach(o = 1:length(unique(DeepSWPA_NoDevNY$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepSWPA_NoDevNY, o = o, LowLim = 2, MT_WGS = DeepSWPA_WGS_NoDevNY, v.MT = v.DeepSWPA_NoDevNY, rv.MT = rv.DeepSWPA_NoDevNY, Vcut = 60000, Vbins = 50, HistSep = 10, Histylim = 1300, RegName = "DeepSWPA_NoDevNY", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
#p-value colors
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(0,0.2)
scaleBy = 0.05
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#Sensitivity Plot for All Metrics
OpRanks_DeepSWPA = GetOperatorRanks(OpDiag_DeepSWPA)
PlotDistanceSensitivity(OpDiagRanks = OpRanks_DeepSWPA, res = 300, height = 10, width = 7, xlim = c(1.5,60), NumOps = 300, NumOpsStep = 50, xStep = 2.5, PlotName = 'DeepSWPA')
#Plot operator diagnostics for each of the distances
temp = foreach (i = 1:ncol(OpRanks_DeepSWPA$DiffMOM), .packages = 'RColorBrewer') %dopar% {
  PlotOpsDiagnostics(OpDiag = OpDiag_DeepSWPA, OpDiagRanks = OpRanks_DeepSWPA, col = i, res = 300, NumOps = 300, NumOpsStep = 50, PlotName = paste0('DeepSWPA',i*2500), xStep = 2.5)
}
rm(temp)
stopCluster(cl)

#     Check operators in WV for SWPA---- 
#Seems like WV may have a different variogram shape than PA in this region. Maybe splitting would be good for 2nd order stationarity.
DeepSWPA_WGS_WV = spTransform(DeepSWPA[DeepSWPA$State == 'WV',], CRS('+init=epsg:4326'))
v.DeepSWPA_WV = variogram(Qs~1, DeepSWPA[DeepSWPA$State == 'WV',], cutoff = 60000, width = 60000/50)
rv.DeepSWPA_WV = variogram(Qs~1, DeepSWPA[DeepSWPA$State == 'WV',], cutoff = 60000, width = 60000/50, cressie = TRUE)

colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
OpDiag_DeepSWPA_WV = foreach(o = 1:length(unique(DeepSWPA$Operator[DeepSWPA$State == 'WV'])), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepSWPA[DeepSWPA$State == 'WV',], o = o, LowLim = 2, MT_WGS = DeepSWPA_WGS_WV, v.MT = v.DeepSWPA_WV, rv.MT = rv.DeepSWPA_WV, Vcut = 60000, Vbins = 50, HistSep = 10, Histylim = 1300, RegName = "DeepSWPA_WV", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
#p-value colors
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(0,0.2)
scaleBy = 0.05
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#Sensitivity Plot for All Metrics
OpRanks_DeepSWPA_WV = GetOperatorRanks(OpDiag_DeepSWPA_WV)
PlotDistanceSensitivity(OpDiagRanks = OpRanks_DeepSWPA_WV, res = 300, height = 10, width = 7, xlim = c(1.5,60), NumOps = 80, NumOpsStep = 20, xStep = 2.5, PlotName = 'DeepSWPA_WV')
#Plot operator diagnostics for each of the distances
temp = foreach (i = 1:ncol(OpRanks_DeepSWPA_WV$DiffMOM), .packages = 'RColorBrewer') %dopar% {
  PlotOpsDiagnostics(OpDiag = OpDiag_DeepSWPA_WV, OpDiagRanks = OpRanks_DeepSWPA_WV, col = i, res = 300, NumOps = 80, NumOpsStep = 20, PlotName = paste0('DeepSWPA_WV',i*2500), xStep = 2.5)
}
rm(temp)
stopCluster(cl)

#Top 10 among all metrics, without p-value
TopBadOps_SWPA_WV = ReturnTopN(OpDiag = OpDiag_DeepSWPA_WV, OpRanks = OpRanks_DeepSWPA_WV, col = 12, N = 10, ID = TRUE, WellData = DeepSWPA_WGS_WV)
TopBadOps_SWPA_WV_Names = ReturnTopN(OpDiag = OpDiag_DeepSWPA_WV, OpRanks = OpRanks_DeepSWPA_WV, col = 12, N = 10, ID = FALSE, WellData = DeepSWPA_WGS_WV)

#Identified Petroleum Development Corp. as an operator that had some logs misinterpreted in AASG. Their logs are fine. 
#EMAX wells logged up in this geologic region - remove.
DeepSWPA = DeepSWPA[DeepSWPA$Operator != 'EMAX, Inc.',]
DeepSWPA_NoDevNY = DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$Operator != 'EMAX, Inc.',]
#With these adjustments, no major change in variogram happens. 
#     Check maps of this region in WV. The variogram is very high at short separation distances----
#There must be other factors affecting variogram in this region.

#Map of Heat Flow, removing extreme high values
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)
plot(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$State == 'WV' & DeepSWPA_NoDevNY$Qs < 100,][order(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$State == 'WV' & DeepSWPA_NoDevNY$Qs < 100,]$WellDepth, decreasing = FALSE),], col = colFun(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$State == 'WV' & DeepSWPA_NoDevNY$Qs < 100,]$Qs[order(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$State == 'WV' & DeepSWPA_NoDevNY$Qs < 100,]$WellDepth, decreasing = FALSE)]), cex = 1, pch = 16)
#Check only deeper than 1200 m - looks similar to > 1000 m. Depth doesn't seem to matter here.
plot(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$State == 'WV' & DeepSWPA_NoDevNY$Qs < 100,][order(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$State == 'WV' & DeepSWPA_NoDevNY$Qs < 100,]$WellDepth, decreasing = FALSE),], col = 'white', cex = 1, pch = 16)
plot(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$State == 'WV' & DeepSWPA_NoDevNY$Qs < 100 & DeepSWPA_NoDevNY$WellDepth > 1200,][order(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$State == 'WV' & DeepSWPA_NoDevNY$Qs < 100 & DeepSWPA_NoDevNY$WellDepth > 1200,]$WellDepth, decreasing = FALSE),], col = colFun(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$State == 'WV' & DeepSWPA_NoDevNY$Qs < 100 & DeepSWPA_NoDevNY$WellDepth > 1200,]$Qs[order(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$State == 'WV' & DeepSWPA_NoDevNY$Qs < 100 & DeepSWPA_NoDevNY$WellDepth > 1200,]$WellDepth, decreasing = FALSE)]), cex = 1, pch = 16, add = T)

#Map of Elevation
#Heat flow on the same spatial scale as elevation below
plot(DeepSWPA_NoDevNY[order(DeepSWPA_NoDevNY$WellDepth, decreasing = FALSE),], col = colFun(DeepSWPA_NoDevNY$Qs[order(DeepSWPA_NoDevNY$WellDepth, decreasing = FALSE)]), cex = 1, pch = 16)

colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(1000,4000)
scaleBy = 200
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#Seems to be some correlation with high elevation => high heat flow, but it's localized where the Allegheny Mtns are.
plot(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$Qs < 100,][order(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$Qs < 100,]$WellDepth, decreasing = FALSE),], 
     col = colFun(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$Qs < 100,]$ElevtnG[order(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$Qs < 100,]$WellDepth, decreasing = FALSE)]), cex = 1, pch = 16)

#WellDepth
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(1000,3000)
scaleBy = 200
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#WellDepth - Seems like some deviated wells may still be present. Not sure how to ID them. Spatial outlier analysis hopefully finds them.
plot(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$State == 'WV' & DeepSWPA_NoDevNY$Qs < 100,][order(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$State == 'WV' & DeepSWPA_NoDevNY$Qs < 100,]$WellDepth, decreasing = FALSE),], col = colFun(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$State == 'WV' & DeepSWPA_NoDevNY$Qs < 100,]$WellDepth[order(DeepSWPA_NoDevNY[DeepSWPA_NoDevNY$State == 'WV' & DeepSWPA_NoDevNY$Qs < 100,]$WellDepth, decreasing = FALSE)]), cex = 1, pch = 16)

dev.off()

#    WPA---- 
#Deviated wells do not look incorrect here.
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
OpDiag_DeepWPA = foreach(o = 1:length(unique(DeepWPA$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepWPA, o = o, LowLim = 2, MT_WGS = DeepWPA_WGS, v.MT = v.DeepWPA, rv.MT = rv.DeepWPA, Vcut = 60000, Vbins = 50, HistSep = 5, Histylim = 800, RegName = "DeepWPA", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
OpDiag_DeepWPA_NoDevNY = foreach(o = 1:length(unique(DeepWPA_NoDevNY$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepWPA_NoDevNY, o = o, LowLim = 2, MT_WGS = DeepWPA_WGS_NoDevNY, v.MT = v.DeepWPA_NoDevNY, rv.MT = rv.DeepWPA_NoDevNY, Vcut = 60000, Vbins = 50, HistSep = 5, Histylim = 800, RegName = "DeepWPA_NoDevNY", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
#p-value colors
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(0,0.2)
scaleBy = 0.05
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#Sensitivity Plot for All Metrics
OpRanks_DeepWPA = GetOperatorRanks(OpDiag_DeepWPA)
PlotDistanceSensitivity(OpDiagRanks = OpRanks_DeepWPA, res = 300, height = 10, width = 7, xlim = c(1.5,60), NumOps = 80, NumOpsStep = 20, xStep = 2.5, PlotName = 'DeepWPA')
#Plot operator diagnostics for each of the distances
temp = foreach (i = 1:ncol(OpRanks_DeepWPA$DiffMOM), .packages = 'RColorBrewer') %dopar% {
  PlotOpsDiagnostics(OpDiag = OpDiag_DeepWPA, OpDiagRanks = OpRanks_DeepWPA, col = i, res = 300, NumOps = 80, NumOpsStep = 20, PlotName = paste0('DeepWPA',i*2500), xStep = 2.5)
}
rm(temp)
stopCluster(cl)

#Top 10 among all metrics, without p-value - These all look reasonable
TopBadOps_WPA = ReturnTopN(OpDiag = OpDiag_DeepWPA, OpRanks = OpRanks_DeepWPA, col = 5, N = 5, ID = TRUE, WellData = DeepWPA)
TopBadOps_WPA_Names = ReturnTopN(OpDiag = OpDiag_DeepWPA, OpRanks = OpRanks_DeepWPA, col = 5, N = 5, ID = FALSE, WellData = DeepWPA)

#    CT---- 
#Deviated wells do not look incorrect here.
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
OpDiag_DeepCT = foreach(o = 1:length(unique(DeepCT$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepCT, o = o, LowLim = 2, MT_WGS = DeepCT_WGS, v.MT = v.DeepCT, rv.MT = rv.DeepCT, Vcut = 60000, Vbins = 50, HistSep = 10, Histylim = 800, RegName = "DeepCT", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
OpDiag_DeepCT_NoDevNY = foreach(o = 1:length(unique(DeepCT_NoDevNY$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepCT_NoDevNY, o = o, LowLim = 2, MT_WGS = DeepCT_WGS_NoDevNY, v.MT = v.DeepCT_NoDevNY, rv.MT = rv.DeepCT_NoDevNY, Vcut = 60000, Vbins = 50, HistSep = 10, Histylim = 800, RegName = "DeepCT_NoDevNY", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
#p-value colors
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(0,0.2)
scaleBy = 0.05
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#Sensitivity Plot for All Metrics
OpRanks_DeepCT = GetOperatorRanks(OpDiag_DeepCT)
PlotDistanceSensitivity(OpDiagRanks = OpRanks_DeepCT, res = 300, height = 10, width = 7, xlim = c(1.5,60), NumOps = 80, NumOpsStep = 20, xStep = 2.5, PlotName = 'DeepCT')
#Plot operator diagnostics for each of the distances
temp = foreach (i = 1:ncol(OpRanks_DeepCT$DiffMOM), .packages = 'RColorBrewer') %dopar% {
  PlotOpsDiagnostics(OpDiag = OpDiag_DeepCT, OpDiagRanks = OpRanks_DeepCT, col = i, res = 300, NumOps = 80, NumOpsStep = 20, PlotName = paste0('DeepCT',i*2500), xStep = 2.5)
}
rm(temp)
stopCluster(cl)

#Top 10 among all metrics, without p-value - These all look reasonable
TopBadOps_CT = ReturnTopN(OpDiag = OpDiag_DeepCT, OpRanks = OpRanks_DeepCT, col = 6, N = 5, ID = TRUE, WellData = DeepCT)
TopBadOps_CT_Names = ReturnTopN(OpDiag = OpDiag_DeepCT, OpRanks = OpRanks_DeepCT, col = 6, N = 5, ID = FALSE, WellData = DeepCT)

#    ENYPA----
#Deviated wells do not look incorrect here.
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
OpDiag_DeepENYPA = foreach(o = 1:length(unique(DeepENYPA$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepENYPA, o = o, LowLim = 2, MT_WGS = DeepENYPA_WGS, v.MT = v.DeepENYPA, rv.MT = rv.DeepENYPA, Vcut = 60000, Vbins = 40, HistSep = 10, Histylim = 200, RegName = "DeepENYPA", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
OpDiag_DeepENYPA_NoDevNY = foreach(o = 1:length(unique(DeepENYPA_NoDevNY$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepENYPA_NoDevNY, o = o, LowLim = 2, MT_WGS = DeepENYPA_WGS_NoDevNY, v.MT = v.DeepENYPA_NoDevNY, rv.MT = rv.DeepENYPA_NoDevNY, Vcut = 60000, Vbins = 40, HistSep = 10, Histylim = 200, RegName = "DeepENYPA_NoDevNY", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
#p-value colors
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(0,0.2)
scaleBy = 0.05
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#Sensitivity Plot for All Metrics
OpRanks_DeepENYPA = GetOperatorRanks(OpDiag_DeepENYPA)
PlotDistanceSensitivity(OpDiagRanks = OpRanks_DeepENYPA, res = 300, height = 10, width = 7, xlim = c(1.5,60), NumOps = 20, NumOpsStep = 5, xStep = 2.5, PlotName = 'DeepENYPA')
#Plot operator diagnostics for each of the distances
temp = foreach (i = 1:ncol(OpRanks_DeepENYPA$DiffMOM), .packages = 'RColorBrewer') %dopar% {
  PlotOpsDiagnostics(OpDiag = OpDiag_DeepENYPA, OpDiagRanks = OpRanks_DeepENYPA, col = i, res = 300, NumOps = 100, NumOpsStep = 20, PlotName = paste0('DeepENYPA',i*2500), xStep = 2.5)
}
rm(temp)
stopCluster(cl)

#Top 10 among all metrics, without p-value - These all look reasonable
TopBadOps_ENYPA = ReturnTopN(OpDiag = OpDiag_DeepENYPA, OpRanks = OpRanks_DeepENYPA, col = 6, N = 5, ID = TRUE, WellData = DeepENYPA)
TopBadOps_ENYPA_Names = ReturnTopN(OpDiag = OpDiag_DeepENYPA, OpRanks = OpRanks_DeepENYPA, col = 6, N = 5, ID = FALSE, WellData = DeepENYPA)

#    NWPANY---- 
#Deviated wells do not look incorrect here.
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
OpDiag_DeepNWPANY = foreach(o = 1:length(unique(DeepNWPANY$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepNWPANY, o = o, LowLim = 2, MT_WGS = DeepNWPANY_WGS, v.MT = v.DeepNWPANY, rv.MT = rv.DeepNWPANY, Vcut = 60000, Vbins = 20, HistSep = 10, Histylim = 200, RegName = "DeepNWPANY", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
OpDiag_DeepNWPANY_NoDevNY = foreach(o = 1:length(unique(DeepNWPANY_NoDevNY$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepNWPANY_NoDevNY, o = o, LowLim = 2, MT_WGS = DeepNWPANY_WGS_NoDevNY, v.MT = v.DeepNWPANY_NoDevNY, rv.MT = rv.DeepNWPANY_NoDevNY, Vcut = 60000, Vbins = 20, HistSep = 10, Histylim = 200, RegName = "DeepNWPANY_NoDevNY", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
#p-value colors
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(0,0.2)
scaleBy = 0.05
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#Sensitivity Plot for All Metrics
OpRanks_DeepNWPANY = GetOperatorRanks(OpDiag_DeepNWPANY)
PlotDistanceSensitivity(OpDiagRanks = OpRanks_DeepNWPANY, res = 300, height = 10, width = 7, xlim = c(1.5,60), NumOps = 20, NumOpsStep = 5, xStep = 2.5, PlotName = 'DeepNWPANY')
#Plot operator diagnostics for each of the distances
temp = foreach (i = 1:ncol(OpRanks_DeepNWPANY$DiffMOM), .packages = 'RColorBrewer') %dopar% {
  PlotOpsDiagnostics(OpDiag = OpDiag_DeepNWPANY, OpDiagRanks = OpRanks_DeepNWPANY, col = i, res = 300, NumOps = 20, NumOpsStep = 5, PlotName = paste0('DeepNWPANY',i*2500), xStep = 2.5)
}
rm(temp)
stopCluster(cl)

#Top 10 among all metrics, without p-value - These all look reasonable
TopBadOps_NWPANY = ReturnTopN(OpDiag = OpDiag_DeepNWPANY, OpRanks = OpRanks_DeepNWPANY, col = 6, N = 5, ID = TRUE, WellData = DeepNWPANY)
TopBadOps_NWPANY_Names = ReturnTopN(OpDiag = OpDiag_DeepNWPANY, OpRanks = OpRanks_DeepNWPANY, col = 6, N = 5, ID = FALSE, WellData = DeepNWPANY)

#    CNY - No operators listed----

#    ENY---- 
#Deviated wells do not look incorrect here.
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
OpDiag_DeepENY = foreach(o = 1:length(unique(DeepENY$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepENY, o = o, LowLim = 2, MT_WGS = DeepENY_WGS, v.MT = v.DeepENY, rv.MT = rv.DeepENY, Vcut = 60000, Vbins = 15, HistSep = 10, Histylim = 200, RegName = "DeepENY", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
OpDiag_DeepENY_NoDevNY = foreach(o = 1:length(unique(DeepENY_NoDevNY$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepENY_NoDevNY, o = o, LowLim = 2, MT_WGS = DeepENY_WGS_NoDevNY, v.MT = v.DeepENY_NoDevNY, rv.MT = rv.DeepENY_NoDevNY, Vcut = 60000, Vbins = 15, HistSep = 10, Histylim = 200, RegName = "DeepENY_NoDevNY", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
#p-value colors
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(0,0.2)
scaleBy = 0.05
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#Sensitivity Plot for All Metrics
OpRanks_DeepENY = GetOperatorRanks(OpDiag_DeepENY)
PlotDistanceSensitivity(OpDiagRanks = OpRanks_DeepENY, res = 300, height = 10, width = 7, xlim = c(1.5,60), NumOps = 20, NumOpsStep = 5, xStep = 2.5, PlotName = 'DeepENY')
#Plot operator diagnostics for each of the distances
temp = foreach (i = 1:ncol(OpRanks_DeepENY$DiffMOM), .packages = 'RColorBrewer') %dopar% {
  PlotOpsDiagnostics(OpDiag = OpDiag_DeepENY, OpDiagRanks = OpRanks_DeepENY, col = i, res = 300, NumOps = 20, NumOpsStep = 5, PlotName = paste0('DeepENY',i*2500), xStep = 2.5)
}
rm(temp)
stopCluster(cl)

#Top 10 among all metrics, without p-value - These all look reasonable
TopBadOps_ENY = ReturnTopN(OpDiag = OpDiag_DeepENY, OpRanks = OpRanks_DeepENY, col = 6, N = 2, ID = TRUE, WellData = DeepENY)
TopBadOps_ENY_Names = ReturnTopN(OpDiag = OpDiag_DeepENY, OpRanks = OpRanks_DeepENY, col = 6, N = 2, ID = FALSE, WellData = DeepENY)

#    VR---- 
#Deviated wells do not look incorrect here.
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
OpDiag_DeepVR = foreach(o = 1:length(unique(DeepVR$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepVR, o = o, LowLim = 2, MT_WGS = DeepVR_WGS, v.MT = v.DeepVR, rv.MT = rv.DeepVR, Vcut = 60000, Vbins = 20, HistSep = 10, Histylim = 200, RegName = "DeepVR", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
OpDiag_DeepVR_NoDevNY = foreach(o = 1:length(unique(DeepVR_NoDevNY$Operator)), .packages = c('gstat', 'GISTools', 'sp'), .combine = 'rbind') %dopar% {
  a = BadOperatorDiagnostics(MT = DeepVR_NoDevNY, o = o, LowLim = 2, MT_WGS = DeepVR_NoDevNY_WGS, v.MT = v.DeepVR_NoDevNY, rv.MT = rv.DeepVR_NoDevNY, Vcut = 60000, Vbins = 20, HistSep = 10, Histylim = 200, RegName = "DeepVR_NoDevNY", MaxLagDist = 60000, SensitivityDist = 2500)
  a
}
#p-value colors
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(0,0.2)
scaleBy = 0.05
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#Sensitivity Plot for All Metrics
OpRanks_DeepVR = GetOperatorRanks(OpDiag_DeepVR)
PlotDistanceSensitivity(OpDiagRanks = OpRanks_DeepVR, res = 300, height = 10, width = 7, xlim = c(1.5,60), NumOps = 40, NumOpsStep = 5, xStep = 2.5, PlotName = 'DeepVR')
#Plot operator diagnostics for each of the distances
temp = foreach (i = 1:ncol(OpRanks_DeepVR$DiffMOM), .packages = 'RColorBrewer') %dopar% {
  PlotOpsDiagnostics(OpDiag = OpDiag_DeepVR, OpDiagRanks = OpRanks_DeepVR, col = i, res = 300, NumOps = 40, NumOpsStep = 5, PlotName = paste0('DeepVR',i*2500), xStep = 2.5)
}
rm(temp)
stopCluster(cl)

#Top 10 among all metrics, without p-value - These all look reasonable
TopBadOps_VR = ReturnTopN(OpDiag = OpDiag_DeepVR, OpRanks = OpRanks_DeepVR, col = 6, N = 5, ID = TRUE, WellData = DeepVR)
TopBadOps_VR_Names = ReturnTopN(OpDiag = OpDiag_DeepVR, OpRanks = OpRanks_DeepVR, col = 6, N = 5, ID = FALSE, WellData = DeepVR)

#   Make new deep well database with bad operators removed----
WellsDeep_RmOps = WellsDeep[-which(WellsDeep$Operator == 'Waco Oil & Gas Co., Inc.'),]
WellsDeep_RmOps = WellsDeep_RmOps[-which(WellsDeep_RmOps$StateID %in% DeepMT$StateID[which(DeepMT$Operator == 'Equitable Production Company')]),]
WellsDeep_RmOps = WellsDeep_RmOps[-which(WellsDeep_RmOps$StateID %in% DeepSWPA_WGS_WV$StateID[DeepSWPA_WGS_WV$Operator == 'EMAX, Inc.']),]
WellsDeep_RmOpsWGS = spTransform(WellsDeep_RmOps, CRS('+init=epsg:4326'))

#Operator removed geologic region data
DeepCT_RmOps = spTransform(WellsDeep_RmOpsWGS[InterpRegs[InterpRegs$Name == 'CT',],], CRS('+init=epsg:26917'))
DeepCNY_RmOps = spTransform(WellsDeep_RmOpsWGS[InterpRegs[InterpRegs$Name == 'CNY',],], CRS('+init=epsg:26917'))
DeepCWV_RmOps = spTransform(WellsDeep_RmOpsWGS[CWV_Bounded,], CRS('+init=epsg:26917'))
DeepENY_RmOps = spTransform(WellsDeep_RmOpsWGS[InterpRegs[InterpRegs$Name == 'ENY',],], CRS('+init=epsg:26917'))
DeepENYPA_RmOps = spTransform(WellsDeep_RmOpsWGS[InterpRegs[InterpRegs$Name == 'ENYPA',],], CRS('+init=epsg:26917'))
DeepMT_RmOps = spTransform(WellsDeep_RmOpsWGS[MT_Bounded,], CRS('+init=epsg:26917'))
DeepNWPANY_RmOps = spTransform(WellsDeep_RmOpsWGS[InterpRegs[InterpRegs$Name == 'NWPANY',],], CRS('+init=epsg:26917'))
DeepSWPA_RmOps = spTransform(WellsDeep_RmOpsWGS[InterpRegs[InterpRegs$Name == 'SWPA',],], CRS('+init=epsg:26917'))
DeepWPA_RmOps = spTransform(WellsDeep_RmOpsWGS[InterpRegs[InterpRegs$Name == 'WPA',],], CRS('+init=epsg:26917'))
DeepVR_RmOps = spTransform(WellsDeep_RmOpsWGS[VR_Bounded,], CRS('+init=epsg:26917'))
DeepFL_RmOps = rbind(DeepCT_RmOps, DeepCWV_RmOps, DeepCNY_RmOps, DeepENY_RmOps, DeepENYPA_RmOps, DeepMT_RmOps, DeepNWPANY_RmOps, DeepSWPA_RmOps, DeepWPA_RmOps, DeepVR_RmOps)

#Resetting all deep data to before operators were removed
DeepCT = spTransform(WellsDeepWGS[InterpRegs[InterpRegs$Name == 'CT',],], CRS('+init=epsg:26917'))
DeepCNY = spTransform(WellsDeepWGS[InterpRegs[InterpRegs$Name == 'CNY',],], CRS('+init=epsg:26917'))
DeepCWV = spTransform(WellsDeepWGS[CWV_Bounded,], CRS('+init=epsg:26917'))
DeepENY = spTransform(WellsDeepWGS[InterpRegs[InterpRegs$Name == 'ENY',],], CRS('+init=epsg:26917'))
DeepENYPA = spTransform(WellsDeepWGS[InterpRegs[InterpRegs$Name == 'ENYPA',],], CRS('+init=epsg:26917'))
DeepMT = spTransform(WellsDeepWGS[MT_Bounded,], CRS('+init=epsg:26917'))
DeepNWPANY = spTransform(WellsDeepWGS[InterpRegs[InterpRegs$Name == 'NWPANY',],], CRS('+init=epsg:26917'))
DeepSWPA = spTransform(WellsDeepWGS[InterpRegs[InterpRegs$Name == 'SWPA',],], CRS('+init=epsg:26917'))
DeepWPA = spTransform(WellsDeepWGS[InterpRegs[InterpRegs$Name == 'WPA',],], CRS('+init=epsg:26917'))
DeepVR = spTransform(WellsDeepWGS[VR_Bounded,], CRS('+init=epsg:26917'))
DeepFL = rbind(DeepCT, DeepCWV, DeepCNY, DeepENY, DeepENYPA, DeepMT, DeepNWPANY, DeepSWPA, DeepWPA, DeepVR)

#  Spatial Outlier Detection ----
#This should be run after points have been reduced to unique locations and negative gradient values have been removed or corrected.
#Bad operators should also be evaluated.

#Data must have a column of UTM coordinates in m for this to work because it relies on Euclidian distances.
WellsDeep = spTransform(WellsDeep, CRS('+init=epsg:26917'))
#Add columns of UTM coordinates
WellsDeep$POINT_X = WellsDeep@coords[,1]
WellsDeep$POINT_Y = WellsDeep@coords[,2]

WellsDeep_RmOps = spTransform(WellsDeep_RmOps, CRS('+init=epsg:26917'))
#Add columns of UTM coordinates
WellsDeep_RmOps$POINT_X = WellsDeep_RmOps@coords[,1]
WellsDeep_RmOps$POINT_Y = WellsDeep_RmOps@coords[,2]

WellsDeep_NoDevNY = spTransform(WellsDeep_NoDevNY, CRS('+init=epsg:26917'))
#Add columns of UTM coordinates
WellsDeep_NoDevNY$POINT_X = WellsDeep_NoDevNY@coords[,1]
WellsDeep_NoDevNY$POINT_Y = WellsDeep_NoDevNY@coords[,2]

TestedOutliers_HeatFlow = select_out_algo(Data = WellsDeep@data, 
                                          InVarName = "Qs", OutVarName = "Qs", 
                                          X_coordName = "POINT_X", Y_coordName = "POINT_Y", 
                                          CensorName = 'WellDepth', 
                                          rad_max = 32000, 
                                          pt_min = 25, 
                                          k_loc = 3, 
                                          type = 7, 
                                          min_val = 0, max_val = 7000, rank = TRUE)

#Testing using dataset without supposedly deviated wells. NY wells all not removed.
TestedOutliers_HeatFlow_NoDevNY = select_out_algo(Data = WellsDeep_NoDevNY@data, 
                                          InVarName = "Qs", OutVarName = "Qs", 
                                          X_coordName = "POINT_X", Y_coordName = "POINT_Y", 
                                          CensorName = 'WellDepth', 
                                          rad_max = 32000, 
                                          pt_min = 25, 
                                          k_loc = 3, 
                                          type = 7, 
                                          min_val = 0, max_val = 7000, rank = TRUE)

#Testing using operator-removed dataset
TestedOutliers_HeatFlow_RmOps = select_out_algo(Data = WellsDeep_RmOps@data, 
                                          InVarName = "Qs", OutVarName = "Qs", 
                                          X_coordName = "POINT_X", Y_coordName = "POINT_Y", 
                                          CensorName = 'WellDepth', 
                                          rad_max = 32000, 
                                          pt_min = 25, 
                                          k_loc = 3, 
                                          type = 7, 
                                          min_val = 0, max_val = 7000, rank = TRUE)

#Testing with a minimum depth of 2000 m. Specifically to see if deep points in the SWPA region are still outliers among deeper data.
TestedOutliers_HeatFlow_min2k = select_out_algo(Data = WellsDeep@data, 
                                          InVarName = "Qs", OutVarName = "Qs", 
                                          X_coordName = "POINT_X", Y_coordName = "POINT_Y", 
                                          CensorName = 'WellDepth', 
                                          rad_max = 32000, 
                                          pt_min = 25, 
                                          k_loc = 3, 
                                          type = 7, 
                                          min_val = 2000, max_val = 7000, rank = TRUE)

#Testing with a maximum depth of 2000 m. Specifically to see if deep points in the SWPA region are again outliers among shallower data.
TestedOutliers_HeatFlow_max2k = select_out_algo(Data = WellsDeep@data, 
                                          InVarName = "Qs", OutVarName = "Qs", 
                                          X_coordName = "POINT_X", Y_coordName = "POINT_Y", 
                                          CensorName = 'WellDepth', 
                                          rad_max = 32000, 
                                          pt_min = 25, 
                                          k_loc = 3, 
                                          type = 7, 
                                          min_val = 0, max_val = 2000, rank = TRUE)

#Testing sensitivity to value of k
TestedOutliers_HeatFlow_k1p0 = select_out_algo(Data = WellsDeep@data, 
                                          InVarName = "Qs", OutVarName = "Qs", 
                                          X_coordName = "POINT_X", Y_coordName = "POINT_Y", 
                                          CensorName = 'WellDepth', 
                                          rad_max = 32000, 
                                          pt_min = 25, 
                                          k_loc = 1, 
                                          type = 7, 
                                          min_val = 0, max_val = 7000, rank = TRUE)

TestedOutliers_HeatFlow_k1p5 = select_out_algo(Data = WellsDeep@data, 
                                               InVarName = "Qs", OutVarName = "Qs", 
                                               X_coordName = "POINT_X", Y_coordName = "POINT_Y", 
                                               CensorName = 'WellDepth', 
                                               rad_max = 32000, 
                                               pt_min = 25, 
                                               k_loc = 1.5, 
                                               type = 7, 
                                               min_val = 0, max_val = 7000, rank = TRUE)

TestedOutliers_HeatFlow_k2p0 = select_out_algo(Data = WellsDeep@data, 
                                               InVarName = "Qs", OutVarName = "Qs", 
                                               X_coordName = "POINT_X", Y_coordName = "POINT_Y", 
                                               CensorName = 'WellDepth', 
                                               rad_max = 32000, 
                                               pt_min = 25, 
                                               k_loc = 2, 
                                               type = 7, 
                                               min_val = 0, max_val = 7000, rank = TRUE)

TestedOutliers_HeatFlow_k2p5 = select_out_algo(Data = WellsDeep@data, 
                                               InVarName = "Qs", OutVarName = "Qs", 
                                               X_coordName = "POINT_X", Y_coordName = "POINT_Y", 
                                               CensorName = 'WellDepth', 
                                               rad_max = 32000, 
                                               pt_min = 25, 
                                               k_loc = 2.5, 
                                               type = 7, 
                                               min_val = 0, max_val = 7000, rank = TRUE)

TestedOutliers_HeatFlow_k3p5 = select_out_algo(Data = WellsDeep@data, 
                                               InVarName = "Qs", OutVarName = "Qs", 
                                               X_coordName = "POINT_X", Y_coordName = "POINT_Y", 
                                               CensorName = 'WellDepth', 
                                               rad_max = 32000, 
                                               pt_min = 25, 
                                               k_loc = 3.5, 
                                               type = 7, 
                                               min_val = 0, max_val = 7000, rank = TRUE)

TestedOutliers_HeatFlow_k4p0 = select_out_algo(Data = WellsDeep@data, 
                                               InVarName = "Qs", OutVarName = "Qs", 
                                               X_coordName = "POINT_X", Y_coordName = "POINT_Y", 
                                               CensorName = 'WellDepth', 
                                               rad_max = 32000, 
                                               pt_min = 25, 
                                               k_loc = 4, 
                                               type = 7, 
                                               min_val = 0, max_val = 7000, rank = TRUE)

TestedOutliers_HeatFlow_k4p5 = select_out_algo(Data = WellsDeep@data, 
                                               InVarName = "Qs", OutVarName = "Qs", 
                                               X_coordName = "POINT_X", Y_coordName = "POINT_Y", 
                                               CensorName = 'WellDepth', 
                                               rad_max = 32000, 
                                               pt_min = 25, 
                                               k_loc = 4.5, 
                                               type = 7, 
                                               min_val = 0, max_val = 7000, rank = TRUE)

TestedOutliers_HeatFlow_k5p0 = select_out_algo(Data = WellsDeep@data, 
                                               InVarName = "Qs", OutVarName = "Qs", 
                                               X_coordName = "POINT_X", Y_coordName = "POINT_Y", 
                                               CensorName = 'WellDepth', 
                                               rad_max = 32000, 
                                               pt_min = 25, 
                                               k_loc = 5, 
                                               type = 7, 
                                               min_val = 0, max_val = 7000, rank = TRUE)

#   Convert to spatial data----
coordinates(TestedOutliers_HeatFlow$NotOutliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow$NotOutliers) = CRS('+init=epsg:26917')
coordinates(TestedOutliers_HeatFlow$Outliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow$Outliers) = CRS('+init=epsg:26917')

coordinates(TestedOutliers_HeatFlow_NoDevNY$NotOutliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_NoDevNY$NotOutliers) = CRS('+init=epsg:26917')
coordinates(TestedOutliers_HeatFlow_NoDevNY$Outliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_NoDevNY$Outliers) = CRS('+init=epsg:26917')

coordinates(TestedOutliers_HeatFlow_RmOps$NotOutliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_RmOps$NotOutliers) = CRS('+init=epsg:26917')
coordinates(TestedOutliers_HeatFlow_RmOps$Outliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_RmOps$Outliers) = CRS('+init=epsg:26917')

coordinates(TestedOutliers_HeatFlow_max2k$NotOutliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_max2k$NotOutliers) = CRS('+init=epsg:26917')
coordinates(TestedOutliers_HeatFlow_max2k$Outliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_max2k$Outliers) = CRS('+init=epsg:26917')

coordinates(TestedOutliers_HeatFlow_min2k$NotOutliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_min2k$NotOutliers) = CRS('+init=epsg:26917')
coordinates(TestedOutliers_HeatFlow_min2k$Outliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_min2k$Outliers) = CRS('+init=epsg:26917')

coordinates(TestedOutliers_HeatFlow_k1p0$NotOutliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_k1p0$NotOutliers) = CRS('+init=epsg:26917')
coordinates(TestedOutliers_HeatFlow_k1p0$Outliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_k1p0$Outliers) = CRS('+init=epsg:26917')

coordinates(TestedOutliers_HeatFlow_k1p5$NotOutliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_k1p5$NotOutliers) = CRS('+init=epsg:26917')
coordinates(TestedOutliers_HeatFlow_k1p5$Outliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_k1p5$Outliers) = CRS('+init=epsg:26917')

coordinates(TestedOutliers_HeatFlow_k2p0$NotOutliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_k2p0$NotOutliers) = CRS('+init=epsg:26917')
coordinates(TestedOutliers_HeatFlow_k2p0$Outliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_k2p0$Outliers) = CRS('+init=epsg:26917')

coordinates(TestedOutliers_HeatFlow_k2p5$NotOutliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_k2p5$NotOutliers) = CRS('+init=epsg:26917')
coordinates(TestedOutliers_HeatFlow_k2p5$Outliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_k2p5$Outliers) = CRS('+init=epsg:26917')

coordinates(TestedOutliers_HeatFlow_k3p5$NotOutliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_k3p5$NotOutliers) = CRS('+init=epsg:26917')
coordinates(TestedOutliers_HeatFlow_k3p5$Outliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_k3p5$Outliers) = CRS('+init=epsg:26917')

coordinates(TestedOutliers_HeatFlow_k4p0$NotOutliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_k4p0$NotOutliers) = CRS('+init=epsg:26917')
coordinates(TestedOutliers_HeatFlow_k4p0$Outliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_k4p0$Outliers) = CRS('+init=epsg:26917')

coordinates(TestedOutliers_HeatFlow_k4p5$NotOutliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_k4p5$NotOutliers) = CRS('+init=epsg:26917')
coordinates(TestedOutliers_HeatFlow_k4p5$Outliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_k4p5$Outliers) = CRS('+init=epsg:26917')

coordinates(TestedOutliers_HeatFlow_k5p0$NotOutliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_k5p0$NotOutliers) = CRS('+init=epsg:26917')
coordinates(TestedOutliers_HeatFlow_k5p0$Outliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_k5p0$Outliers) = CRS('+init=epsg:26917')

writeOGR(TestedOutliers_HeatFlow$NotOutliers, dsn=getwd(), layer="DeepestWells_NotOutliers_32km_Qs_CorrBase_Ranked_2018", driver = "ESRI Shapefile")
writeOGR(TestedOutliers_HeatFlow$Outliers, dsn=getwd(), layer="DeepestWells_Outliers_32km_Qs_CorrBase_Ranked_2018", driver = "ESRI Shapefile")
writeOGR(TestedOutliers_HeatFlow_RmOps$NotOutliers, dsn=getwd(), layer="DeepestWells_NotOutliers_32km_Qs_CorrBase_Ranked_2018_RmOps", driver = "ESRI Shapefile")
writeOGR(TestedOutliers_HeatFlow_RmOps$Outliers, dsn=getwd(), layer="DeepestWells_Outliers_32km_Qs_CorrBase_Ranked_2018_RmOps", driver = "ESRI Shapefile")

#   Convert to WGS coordinate system----
#Wells in NAD83 WGS84
Outs = spTransform(TestedOutliers_HeatFlow$Outliers, CRS = CRS("+init=epsg:4326"))
NoOuts = spTransform(TestedOutliers_HeatFlow$NotOutliers, CRS = CRS("+init=epsg:4326"))
AllData = rbind(NoOuts, Outs)

Outs_min2k = spTransform(TestedOutliers_HeatFlow_min2k$Outliers, CRS = CRS("+init=epsg:4326"))
NoOuts_min2k = spTransform(TestedOutliers_HeatFlow_min2k$NotOutliers, CRS = CRS("+init=epsg:4326"))

Outs_max2k = spTransform(TestedOutliers_HeatFlow_max2k$Outliers, CRS = CRS("+init=epsg:4326"))
NoOuts_max2k = spTransform(TestedOutliers_HeatFlow_max2k$NotOutliers, CRS = CRS("+init=epsg:4326"))

WellsDeepWGS = spTransform(WellsDeep, CRSobj = CRS('+init=epsg:4326'))
NotOutliersWGS = spTransform(TestedOutliers_HeatFlow$NotOutliers, CRSobj = CRS('+init=epsg:4326'))
OutliersWGS = spTransform(TestedOutliers_HeatFlow$Outliers, CRSobj = CRS('+init=epsg:4326'))

WellsDeepWGS_NoDevNY = spTransform(WellsDeep_NoDevNY, CRSobj = CRS('+init=epsg:4326'))
NotOutliersWGS_NoDevNY = spTransform(TestedOutliers_HeatFlow_NoDevNY$NotOutliers, CRSobj = CRS('+init=epsg:4326'))
OutliersWGS_NoDevNY = spTransform(TestedOutliers_HeatFlow_NoDevNY$Outliers, CRSobj = CRS('+init=epsg:4326'))

WellsDeepWGS_RmOps = spTransform(WellsDeep_RmOps, CRSobj = CRS('+init=epsg:4326'))
NotOutliersWGS_RmOps = spTransform(TestedOutliers_HeatFlow_RmOps$NotOutliers, CRSobj = CRS('+init=epsg:4326'))
OutliersWGS_RmOps = spTransform(TestedOutliers_HeatFlow_RmOps$Outliers, CRSobj = CRS('+init=epsg:4326'))

NotOutliersWGS_k1p0 = spTransform(TestedOutliers_HeatFlow_k1p0$NotOutliers, CRSobj = CRS('+init=epsg:4326'))
OutliersWGS_k1p0 = spTransform(TestedOutliers_HeatFlow_k1p0$Outliers, CRSobj = CRS('+init=epsg:4326'))
NotOutliersWGS_k1p5 = spTransform(TestedOutliers_HeatFlow_k1p5$NotOutliers, CRSobj = CRS('+init=epsg:4326'))
OutliersWGS_k1p5 = spTransform(TestedOutliers_HeatFlow_k1p5$Outliers, CRSobj = CRS('+init=epsg:4326'))
NotOutliersWGS_k2p0 = spTransform(TestedOutliers_HeatFlow_k2p0$NotOutliers, CRSobj = CRS('+init=epsg:4326'))
OutliersWGS_k2p0 = spTransform(TestedOutliers_HeatFlow_k2p0$Outliers, CRSobj = CRS('+init=epsg:4326'))
NotOutliersWGS_k2p5 = spTransform(TestedOutliers_HeatFlow_k2p5$NotOutliers, CRSobj = CRS('+init=epsg:4326'))
OutliersWGS_k2p5 = spTransform(TestedOutliers_HeatFlow_k2p5$Outliers, CRSobj = CRS('+init=epsg:4326'))
NotOutliersWGS_k3p5 = spTransform(TestedOutliers_HeatFlow_k3p5$NotOutliers, CRSobj = CRS('+init=epsg:4326'))
OutliersWGS_k3p5 = spTransform(TestedOutliers_HeatFlow_k3p5$Outliers, CRSobj = CRS('+init=epsg:4326'))
NotOutliersWGS_k4p0 = spTransform(TestedOutliers_HeatFlow_k4p0$NotOutliers, CRSobj = CRS('+init=epsg:4326'))
OutliersWGS_k4p0 = spTransform(TestedOutliers_HeatFlow_k4p0$Outliers, CRSobj = CRS('+init=epsg:4326'))
NotOutliersWGS_k4p5 = spTransform(TestedOutliers_HeatFlow_k4p5$NotOutliers, CRSobj = CRS('+init=epsg:4326'))
OutliersWGS_k4p5 = spTransform(TestedOutliers_HeatFlow_k4p5$Outliers, CRSobj = CRS('+init=epsg:4326'))
NotOutliersWGS_k5p0 = spTransform(TestedOutliers_HeatFlow_k5p0$NotOutliers, CRSobj = CRS('+init=epsg:4326'))
OutliersWGS_k5p0 = spTransform(TestedOutliers_HeatFlow_k5p0$Outliers, CRSobj = CRS('+init=epsg:4326'))

#   Clip all data into their geologic regions----
#Clip points that are not outliers to geologic regions
CT = NotOutliersWGS[InterpRegs[InterpRegs$Name == 'CT',],]
CNY = NotOutliersWGS[InterpRegs[InterpRegs$Name == 'CNY',],]
CWV = NotOutliersWGS[CWV_Bounded,]
ENY = NotOutliersWGS[InterpRegs[InterpRegs$Name == 'ENY',],]
ENYPA = NotOutliersWGS[InterpRegs[InterpRegs$Name == 'ENYPA',],]
MT = NotOutliersWGS[MT_Bounded,]
NWPANY = NotOutliersWGS[InterpRegs[InterpRegs$Name == 'NWPANY',],]
SWPA = NotOutliersWGS[InterpRegs[InterpRegs$Name == 'SWPA',],]
WPA = NotOutliersWGS[InterpRegs[InterpRegs$Name == 'WPA',],]
VR = NotOutliersWGS[VR_Bounded,]
FL = rbind(CT, CWV, CNY, ENY, ENYPA, MT, NWPANY, SWPA, WPA, VR)

CT_RmOps = NotOutliersWGS_RmOps[InterpRegs[InterpRegs$Name == 'CT',],]
CNY_RmOps = NotOutliersWGS_RmOps[InterpRegs[InterpRegs$Name == 'CNY',],]
CWV_RmOps = NotOutliersWGS_RmOps[CWV_Bounded,]
ENY_RmOps = NotOutliersWGS_RmOps[InterpRegs[InterpRegs$Name == 'ENY',],]
ENYPA_RmOps = NotOutliersWGS_RmOps[InterpRegs[InterpRegs$Name == 'ENYPA',],]
MT_RmOps = NotOutliersWGS_RmOps[MT_Bounded,]
NWPANY_RmOps = NotOutliersWGS_RmOps[InterpRegs[InterpRegs$Name == 'NWPANY',],]
SWPA_RmOps = NotOutliersWGS_RmOps[InterpRegs[InterpRegs$Name == 'SWPA',],]
WPA_RmOps = NotOutliersWGS_RmOps[InterpRegs[InterpRegs$Name == 'WPA',],]
VR_RmOps = NotOutliersWGS_RmOps[VR_Bounded,]
FL_RmOps = rbind(CT_RmOps, CWV_RmOps, CNY_RmOps, ENY_RmOps, ENYPA_RmOps, MT_RmOps, NWPANY_RmOps, SWPA_RmOps, WPA_RmOps, VR_RmOps)

CT_k1p0 = NotOutliersWGS_k1p0[InterpRegs[InterpRegs$Name == 'CT',],]
CNY_k1p0 = NotOutliersWGS_k1p0[InterpRegs[InterpRegs$Name == 'CNY',],]
CWV_k1p0 = NotOutliersWGS_k1p0[CWV_Bounded,]
ENY_k1p0 = NotOutliersWGS_k1p0[InterpRegs[InterpRegs$Name == 'ENY',],]
ENYPA_k1p0 = NotOutliersWGS_k1p0[InterpRegs[InterpRegs$Name == 'ENYPA',],]
MT_k1p0 = NotOutliersWGS_k1p0[MT_Bounded,]
NWPANY_k1p0 = NotOutliersWGS_k1p0[InterpRegs[InterpRegs$Name == 'NWPANY',],]
SWPA_k1p0 = NotOutliersWGS_k1p0[InterpRegs[InterpRegs$Name == 'SWPA',],]
WPA_k1p0 = NotOutliersWGS_k1p0[InterpRegs[InterpRegs$Name == 'WPA',],]
VR_k1p0 = NotOutliersWGS_k1p0[VR_Bounded,]
FL_k1p0 = rbind(CT_k1p0, CWV_k1p0, CNY_k1p0, ENY_k1p0, ENYPA_k1p0, MT_k1p0, NWPANY_k1p0, SWPA_k1p0, WPA_k1p0, VR_k1p0)

CT_k1p5 = NotOutliersWGS_k1p5[InterpRegs[InterpRegs$Name == 'CT',],]
CNY_k1p5 = NotOutliersWGS_k1p5[InterpRegs[InterpRegs$Name == 'CNY',],]
CWV_k1p5 = NotOutliersWGS_k1p5[CWV_Bounded,]
ENY_k1p5 = NotOutliersWGS_k1p5[InterpRegs[InterpRegs$Name == 'ENY',],]
ENYPA_k1p5 = NotOutliersWGS_k1p5[InterpRegs[InterpRegs$Name == 'ENYPA',],]
MT_k1p5 = NotOutliersWGS_k1p5[MT_Bounded,]
NWPANY_k1p5 = NotOutliersWGS_k1p5[InterpRegs[InterpRegs$Name == 'NWPANY',],]
SWPA_k1p5 = NotOutliersWGS_k1p5[InterpRegs[InterpRegs$Name == 'SWPA',],]
WPA_k1p5 = NotOutliersWGS_k1p5[InterpRegs[InterpRegs$Name == 'WPA',],]
VR_k1p5 = NotOutliersWGS_k1p5[VR_Bounded,]
FL_k1p5 = rbind(CT_k1p5, CWV_k1p5, CNY_k1p5, ENY_k1p5, ENYPA_k1p5, MT_k1p5, NWPANY_k1p5, SWPA_k1p5, WPA_k1p5, VR_k1p5)

CT_k2p0 = NotOutliersWGS_k2p0[InterpRegs[InterpRegs$Name == 'CT',],]
CNY_k2p0 = NotOutliersWGS_k2p0[InterpRegs[InterpRegs$Name == 'CNY',],]
CWV_k2p0 = NotOutliersWGS_k2p0[CWV_Bounded,]
ENY_k2p0 = NotOutliersWGS_k2p0[InterpRegs[InterpRegs$Name == 'ENY',],]
ENYPA_k2p0 = NotOutliersWGS_k2p0[InterpRegs[InterpRegs$Name == 'ENYPA',],]
MT_k2p0 = NotOutliersWGS_k2p0[MT_Bounded,]
NWPANY_k2p0 = NotOutliersWGS_k2p0[InterpRegs[InterpRegs$Name == 'NWPANY',],]
SWPA_k2p0 = NotOutliersWGS_k2p0[InterpRegs[InterpRegs$Name == 'SWPA',],]
WPA_k2p0 = NotOutliersWGS_k2p0[InterpRegs[InterpRegs$Name == 'WPA',],]
VR_k2p0 = NotOutliersWGS_k2p0[VR_Bounded,]
FL_k2p0 = rbind(CT_k2p0, CWV_k2p0, CNY_k2p0, ENY_k2p0, ENYPA_k2p0, MT_k2p0, NWPANY_k2p0, SWPA_k2p0, WPA_k2p0, VR_k2p0)

CT_k2p5 = NotOutliersWGS_k2p5[InterpRegs[InterpRegs$Name == 'CT',],]
CNY_k2p5 = NotOutliersWGS_k2p5[InterpRegs[InterpRegs$Name == 'CNY',],]
CWV_k2p5 = NotOutliersWGS_k2p5[CWV_Bounded,]
ENY_k2p5 = NotOutliersWGS_k2p5[InterpRegs[InterpRegs$Name == 'ENY',],]
ENYPA_k2p5 = NotOutliersWGS_k2p5[InterpRegs[InterpRegs$Name == 'ENYPA',],]
MT_k2p5 = NotOutliersWGS_k2p5[MT_Bounded,]
NWPANY_k2p5 = NotOutliersWGS_k2p5[InterpRegs[InterpRegs$Name == 'NWPANY',],]
SWPA_k2p5 = NotOutliersWGS_k2p5[InterpRegs[InterpRegs$Name == 'SWPA',],]
WPA_k2p5 = NotOutliersWGS_k2p5[InterpRegs[InterpRegs$Name == 'WPA',],]
VR_k2p5 = NotOutliersWGS_k2p5[VR_Bounded,]
FL_k2p5 = rbind(CT_k2p5, CWV_k2p5, CNY_k2p5, ENY_k2p5, ENYPA_k2p5, MT_k2p5, NWPANY_k2p5, SWPA_k2p5, WPA_k2p5, VR_k2p5)

CT_k3p5 = NotOutliersWGS_k3p5[InterpRegs[InterpRegs$Name == 'CT',],]
CNY_k3p5 = NotOutliersWGS_k3p5[InterpRegs[InterpRegs$Name == 'CNY',],]
CWV_k3p5 = NotOutliersWGS_k3p5[CWV_Bounded,]
ENY_k3p5 = NotOutliersWGS_k3p5[InterpRegs[InterpRegs$Name == 'ENY',],]
ENYPA_k3p5 = NotOutliersWGS_k3p5[InterpRegs[InterpRegs$Name == 'ENYPA',],]
MT_k3p5 = NotOutliersWGS_k3p5[MT_Bounded,]
NWPANY_k3p5 = NotOutliersWGS_k3p5[InterpRegs[InterpRegs$Name == 'NWPANY',],]
SWPA_k3p5 = NotOutliersWGS_k3p5[InterpRegs[InterpRegs$Name == 'SWPA',],]
WPA_k3p5 = NotOutliersWGS_k3p5[InterpRegs[InterpRegs$Name == 'WPA',],]
VR_k3p5 = NotOutliersWGS_k3p5[VR_Bounded,]
FL_k3p5 = rbind(CT_k3p5, CWV_k3p5, CNY_k3p5, ENY_k3p5, ENYPA_k3p5, MT_k3p5, NWPANY_k3p5, SWPA_k3p5, WPA_k3p5, VR_k3p5)

CT_k4p0 = NotOutliersWGS_k4p0[InterpRegs[InterpRegs$Name == 'CT',],]
CNY_k4p0 = NotOutliersWGS_k4p0[InterpRegs[InterpRegs$Name == 'CNY',],]
CWV_k4p0 = NotOutliersWGS_k4p0[CWV_Bounded,]
ENY_k4p0 = NotOutliersWGS_k4p0[InterpRegs[InterpRegs$Name == 'ENY',],]
ENYPA_k4p0 = NotOutliersWGS_k4p0[InterpRegs[InterpRegs$Name == 'ENYPA',],]
MT_k4p0 = NotOutliersWGS_k4p0[MT_Bounded,]
NWPANY_k4p0 = NotOutliersWGS_k4p0[InterpRegs[InterpRegs$Name == 'NWPANY',],]
SWPA_k4p0 = NotOutliersWGS_k4p0[InterpRegs[InterpRegs$Name == 'SWPA',],]
WPA_k4p0 = NotOutliersWGS_k4p0[InterpRegs[InterpRegs$Name == 'WPA',],]
VR_k4p0 = NotOutliersWGS_k4p0[VR_Bounded,]
FL_k4p0 = rbind(CT_k4p0, CWV_k4p0, CNY_k4p0, ENY_k4p0, ENYPA_k4p0, MT_k4p0, NWPANY_k4p0, SWPA_k4p0, WPA_k4p0, VR_k4p0)

CT_k4p5 = NotOutliersWGS_k4p5[InterpRegs[InterpRegs$Name == 'CT',],]
CNY_k4p5 = NotOutliersWGS_k4p5[InterpRegs[InterpRegs$Name == 'CNY',],]
CWV_k4p5 = NotOutliersWGS_k4p5[CWV_Bounded,]
ENY_k4p5 = NotOutliersWGS_k4p5[InterpRegs[InterpRegs$Name == 'ENY',],]
ENYPA_k4p5 = NotOutliersWGS_k4p5[InterpRegs[InterpRegs$Name == 'ENYPA',],]
MT_k4p5 = NotOutliersWGS_k4p5[MT_Bounded,]
NWPANY_k4p5 = NotOutliersWGS_k4p5[InterpRegs[InterpRegs$Name == 'NWPANY',],]
SWPA_k4p5 = NotOutliersWGS_k4p5[InterpRegs[InterpRegs$Name == 'SWPA',],]
WPA_k4p5 = NotOutliersWGS_k4p5[InterpRegs[InterpRegs$Name == 'WPA',],]
VR_k4p5 = NotOutliersWGS_k4p5[VR_Bounded,]
FL_k4p5 = rbind(CT_k4p5, CWV_k4p5, CNY_k4p5, ENY_k4p5, ENYPA_k4p5, MT_k4p5, NWPANY_k4p5, SWPA_k4p5, WPA_k4p5, VR_k4p5)

CT_k5p0 = NotOutliersWGS_k5p0[InterpRegs[InterpRegs$Name == 'CT',],]
CNY_k5p0 = NotOutliersWGS_k5p0[InterpRegs[InterpRegs$Name == 'CNY',],]
CWV_k5p0 = NotOutliersWGS_k5p0[CWV_Bounded,]
ENY_k5p0 = NotOutliersWGS_k5p0[InterpRegs[InterpRegs$Name == 'ENY',],]
ENYPA_k5p0 = NotOutliersWGS_k5p0[InterpRegs[InterpRegs$Name == 'ENYPA',],]
MT_k5p0 = NotOutliersWGS_k5p0[MT_Bounded,]
NWPANY_k5p0 = NotOutliersWGS_k5p0[InterpRegs[InterpRegs$Name == 'NWPANY',],]
SWPA_k5p0 = NotOutliersWGS_k5p0[InterpRegs[InterpRegs$Name == 'SWPA',],]
WPA_k5p0 = NotOutliersWGS_k5p0[InterpRegs[InterpRegs$Name == 'WPA',],]
VR_k5p0 = NotOutliersWGS_k5p0[VR_Bounded,]
FL_k5p0 = rbind(CT_k5p0, CWV_k5p0, CNY_k5p0, ENY_k5p0, ENYPA_k5p0, MT_k5p0, NWPANY_k5p0, SWPA_k5p0, WPA_k5p0, VR_k5p0)

#   Evaluate the impact of removing the L lowest data points and H highest data points in a region----
#Fixme: make into a function
# L and H are the number of outliers identified in the local spatial outlier algorithm for k = 3.
OutsCT = Outs[InterpRegs[InterpRegs$Name == 'CT',],]
OutsCT = DeepCT[order(DeepCT$Qs),][-c(1:sum(OutsCT$out_loc_lo),(nrow(DeepCT)-(sum(OutsCT$out_loc_hi)-1)):nrow(DeepCT)),]
t1 = plot(variogram(Qs~1, OutsCT, cutoff = 60000, width = 60000/50), ylim = c(0,20), col = 'black', pch = 16,
          ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Chautauqua, NY")
t2 = plot(variogram(Qs~1, spTransform(CT, CRS('+init=epsg:26917')), cutoff = 60000, width = 60000/50), ylim = c(0,20), col = 'blue', pch = 16,
          ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Chautauqua, NY")
plot(t1, more = T)
plot(t2)

OutsWPA = Outs[InterpRegs[InterpRegs$Name == 'WPA',],]
OutsWPA = DeepWPA[order(DeepWPA$Qs),][-c(1:sum(OutsWPA$out_loc_lo),(nrow(DeepWPA)-(sum(OutsWPA$out_loc_hi)-1)):nrow(DeepWPA)),]
t3 = plot(variogram(Qs~1, OutsWPA, cutoff = 60000, width = 60000/50), ylim = c(0,20), col = 'black', pch = 16,
          ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Western PA")
t4 = plot(variogram(Qs~1, spTransform(WPA, CRS('+init=epsg:26917')), cutoff = 60000, width = 60000/50), ylim = c(0,20), col = 'blue', pch = 16,
          ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Western PA")
plot(t3, more = T)
plot(t4)

OutsSWPA = Outs[InterpRegs[InterpRegs$Name == 'SWPA',],]
OutsSWPA = DeepSWPA[order(DeepSWPA$Qs),][-c(1:sum(OutsSWPA$out_loc_lo),(nrow(DeepSWPA)-(sum(OutsSWPA$out_loc_hi)-1)):nrow(DeepSWPA)),]
t5 = plot(variogram(Qs~1, OutsSWPA, cutoff = 60000, width = 60000/50), ylim = c(0,60), col = 'black', pch = 16,
          ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Southwestern PA")
t6 = plot(variogram(Qs~1, spTransform(SWPA, CRS('+init=epsg:26917')), cutoff = 60000, width = 60000/50), ylim = c(0,60), col = 'blue', pch = 16,
          ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Southwestern PA")
plot(t5, more = T)
plot(t6)

OutsCWV = Outs[CWV_Bounded,]
OutsCWV = DeepCWV[order(DeepCWV$Qs),][-c(1:sum(OutsCWV$out_loc_lo),(nrow(DeepCWV)-(sum(OutsCWV$out_loc_hi)-1)):nrow(DeepCWV)),]
t7 = plot(variogram(Qs~1, OutsCWV, cutoff = 60000, width = 60000/50), ylim = c(0,200), col = 'black', pch = 16,
          ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Central WV")
t8 = plot(variogram(Qs~1, spTransform(CWV, CRS('+init=epsg:26917')), cutoff = 60000, width = 60000/50), ylim = c(0,200), col = 'blue', pch = 16,
          ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Central WV")
plot(t7, more = T)
plot(t8)

OutsMT = Outs[MT_Bounded,]
OutsMT = DeepMT[order(DeepMT$Qs),][-c(1:sum(OutsMT$out_loc_lo),(nrow(DeepMT)-(sum(OutsMT$out_loc_hi)-1)):nrow(DeepMT)),]
t9 = plot(variogram(Qs~1, OutsMT, cutoff = 60000, width = 60000/50), ylim = c(0,300), col = 'black', pch = 16,
          ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Western WV")
t10 = plot(variogram(Qs~1, spTransform(MT, CRS('+init=epsg:26917')), cutoff = 60000, width = 60000/50), ylim = c(0,300), col = 'blue', pch = 16,
           ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Western WV")
plot(t9, more = T)
plot(t10)

OutsVR = Outs[VR_Bounded,]
OutsVR = DeepVR[order(DeepVR$Qs),][-c(1:sum(OutsVR$out_loc_lo),(nrow(DeepVR)-(sum(OutsVR$out_loc_hi)-1)):nrow(DeepVR)),]
t11 = plot(variogram(Qs~1, OutsVR, cutoff = 60000, width = 60000/10), ylim = c(0,300), col = 'black', pch = 16,
           ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Valley and Ridge")
t12 = plot(variogram(Qs~1, spTransform(VR, CRS('+init=epsg:26917')), cutoff = 60000, width = 60000/10), ylim = c(0,300), col = 'blue', pch = 16,
           ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Valley and Ridge")
plot(t11, more = T)
plot(t12)

OutsCNY = Outs[InterpRegs[InterpRegs$Name == 'CNY',],]
OutsCNY = DeepCNY[order(DeepCNY$Qs),][-c(1:sum(OutsCNY$out_loc_lo),(nrow(DeepCNY)-(sum(OutsCNY$out_loc_hi)-1)):nrow(DeepCNY)),]
t13 = plot(variogram(Qs~1, OutsCNY, cutoff = 60000, width = 60000/20), ylim = c(0,40), col = 'black', pch = 16,
           ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Central NY")
t14 = plot(variogram(Qs~1, spTransform(CNY, CRS('+init=epsg:26917')), cutoff = 60000, width = 60000/20), ylim = c(0,40), col = 'blue', pch = 16,
           ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Central NY")
plot(t13, more = T)
plot(t14)

OutsENY = Outs[InterpRegs[InterpRegs$Name == 'ENY',],]
OutsENY = DeepENY[order(DeepENY$Qs),][-c(1:sum(OutsENY$out_loc_lo),(nrow(DeepENY)-(sum(OutsENY$out_loc_hi)-1)):nrow(DeepENY)),]
t15 = plot(variogram(Qs~1, OutsENY, cutoff = 60000, width = 60000/20), ylim = c(0,50), col = 'black', pch = 16,
           ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Eastern NY")
t16 = plot(variogram(Qs~1, spTransform(ENY, CRS('+init=epsg:26917')), cutoff = 60000, width = 60000/20), ylim = c(0,50), col = 'blue', pch = 16,
           ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Eastern NY")
plot(t15, more = T)
plot(t16)

OutsENYPA = Outs[InterpRegs[InterpRegs$Name == 'ENYPA',],]
OutsENYPA = DeepENYPA[order(DeepENYPA$Qs),][-c(1:sum(OutsENYPA$out_loc_lo),(nrow(DeepENYPA)-(sum(OutsENYPA$out_loc_hi)-1)):nrow(DeepENYPA)),]
t17 = plot(variogram(Qs~1, OutsENYPA, cutoff = 60000, width = 60000/30), ylim = c(0,150), col = 'black', pch = 16,
           ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Eastern NY & PA")
t18 = plot(variogram(Qs~1, spTransform(ENYPA, CRS('+init=epsg:26917')), cutoff = 60000, width = 60000/30), ylim = c(0,150), col = 'blue', pch = 16,
           ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Eastern NY & PA")
plot(t17, more = T)
plot(t18)

OutsNWPANY = Outs[InterpRegs[InterpRegs$Name == 'NWPANY',],]
OutsNWPANY = DeepNWPANY[order(DeepNWPANY$Qs),][-c(1:sum(OutsNWPANY$out_loc_lo),(nrow(DeepNWPANY)-(sum(OutsNWPANY$out_loc_hi)-1)):nrow(DeepNWPANY)),]
t19 = plot(variogram(Qs~1, OutsNWPANY, cutoff = 60000, width = 60000/50), ylim = c(0,60), col = 'black', pch = 16,
           ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Northwestern PA & NY")
t20 = plot(variogram(Qs~1, spTransform(NWPANY, CRS('+init=epsg:26917')), cutoff = 60000, width = 60000/50), ylim = c(0,60), col = 'blue', pch = 16,
           ylab=expression(Semivariance ~ (mW/m^2)^2), xlab = 'Separation Distance (m)', main="Northwestern PA & NY")
plot(t19, more = T)
plot(t20)


png("Variograms_SpatialOutsVsRmExtremes.png", width=13, height=10, units="in", res=300)
plot(t1, split = c(1,1,4,3), more=T)
plot(t2, split = c(1,1,4,3), more=T)

plot(t13, split = c(4,1,4,3), more=T)
plot(t14, split = c(4,1,4,3), more=T)

plot(t15, split = c(1,2,4,3), more=T)
plot(t16, split = c(1,2,4,3), more=T)

plot(t17, split = c(2,2,4,3), more=T)
plot(t18, split = c(2,2,4,3), more=T)

plot(t19, split = c(3,1,4,3), more=T)
plot(t20, split = c(3,1,4,3), more=T)

plot(t5, split = c(3,2,4,3), more=T)
plot(t6, split = c(3,2,4,3), more=T)

plot(t3, split = c(2,1,4,3), more=T)
plot(t4, split = c(2,1,4,3), more=T)

plot(t9, split = c(1,3,4,3), more=T)
plot(t10, split = c(1,3,4,3), more=T)

plot(t7, split = c(4,2,4,3), more=T)
plot(t8, split = c(4,2,4,3), more=T)

plot(t11, split = c(2,3,4,3), more=T)
plot(t12, split = c(2,3,4,3), more=F)

dev.off()

#   Map of outliers----
#Colors
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(10,90)
scaleBy = 20
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

png('LoHiOuts_Map.png', res = 1200, units = 'in', width = 14, height = 7)
layout(cbind(1,2))
par(xaxs = 'i', yaxs = 'i', mar = c(2,3,3,1))
plot(WellsDeepWGS, col = 'white', pch = 16, cex = 0.2, main = 'Low Outliers', cex.main = 2)
plot(Counties[Counties$STATEFP %in% c(42,36,54,24,21,51),], border = 'grey', add=TRUE)
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS[NotOutliersWGS$out_loc_error == 0,], pch = 16, cex = 0.2, add = TRUE)
plot(NotOutliersWGS[NotOutliersWGS$out_loc_error == 1,], pch = 17, col = 'purple', cex = 0.4, add = TRUE)
plot(OutliersWGS[OutliersWGS$out_loc_lo == 1,], pch = 16, col = colFun(OutliersWGS$Qs[OutliersWGS$out_loc_lo == 1]), cex = 0.7, add = TRUE)
box()
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('topleft', title = expression(paste('Q'[s], ' (mW/m'^2, ')')), legend = c('< 30', '30 - 50', '50 - 70', '70 - 90', '>= 90', 'Not Tested', 'Not Outlier'), col = c(colFun(c(20,40,60,80,100)), 'purple', 'black'), pch = c(16,16,16,16,16,17,16))

par(xaxs = 'i', yaxs = 'i', mar = c(2,3,3,1))
plot(WellsDeepWGS, col = 'white', pch = 16, cex = 0.2, main = 'High Outliers', cex.main = 2)
plot(Counties[Counties$STATEFP %in% c(42,36,54,24,21,51),], border = 'grey', add=TRUE)
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS[NotOutliersWGS$out_loc_error == 0,], pch = 16, cex = 0.2, add = TRUE)
plot(NotOutliersWGS[NotOutliersWGS$out_loc_error == 1,], pch = 17, col = 'purple', cex = 0.4, add = TRUE)
plot(OutliersWGS[OutliersWGS$out_loc_hi == 1,], pch = 16, col = colFun(OutliersWGS$Qs[OutliersWGS$out_loc_hi == 1]), cex = 0.7, add = TRUE)
box()
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('topleft', title = expression(paste('Q'[s], ' (mW/m'^2, ')')), legend = c('< 30', '30 - 50', '50 - 70', '70 - 90', '>= 90', 'Not Tested', 'Not Outlier'), col = c(colFun(c(20,40,60,80,100)), 'purple', 'black'), pch = c(16,16,16,16,16,17,16))
dev.off()

png('LoHiOuts_Map_NoDevNY.png', res = 1200, units = 'in', width = 14, height = 7)
layout(cbind(1,2))
par(xaxs = 'i', yaxs = 'i', mar = c(2,3,3,1))
plot(WellsDeepWGS_NoDevNY, col = 'white', pch = 16, cex = 0.2, main = 'Low Outliers', cex.main = 2)
plot(Counties[Counties$STATEFP %in% c(42,36,54,24,21,51),], border = 'grey', add=TRUE)
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_NoDevNY[NotOutliersWGS_NoDevNY$out_loc_error == 0,], pch = 16, cex = 0.2, add = TRUE)
plot(NotOutliersWGS_NoDevNY[NotOutliersWGS_NoDevNY$out_loc_error == 1,], pch = 17, col = 'purple', cex = 0.4, add = TRUE)
plot(OutliersWGS_NoDevNY[OutliersWGS_NoDevNY$out_loc_lo == 1,], pch = 16, col = colFun(OutliersWGS_NoDevNY$Qs[OutliersWGS_NoDevNY$out_loc_lo == 1]), cex = 0.7, add = TRUE)
box()
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('topleft', title = expression(paste('Q'[s], ' (mW/m'^2, ')')), legend = c('< 30', '30 - 50', '50 - 70', '70 - 90', '>= 90', 'Not Tested', 'Not Outlier'), col = c(colFun(c(20,40,60,80,100)), 'purple', 'black'), pch = c(16,16,16,16,16,17,16))

par(xaxs = 'i', yaxs = 'i', mar = c(2,3,3,1))
plot(WellsDeepWGS_NoDevNY, col = 'white', pch = 16, cex = 0.2, main = 'High Outliers', cex.main = 2)
plot(Counties[Counties$STATEFP %in% c(42,36,54,24,21,51),], border = 'grey', add=TRUE)
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_NoDevNY[NotOutliersWGS_NoDevNY$out_loc_error == 0,], pch = 16, cex = 0.2, add = TRUE)
plot(NotOutliersWGS_NoDevNY[NotOutliersWGS_NoDevNY$out_loc_error == 1,], pch = 17, col = 'purple', cex = 0.4, add = TRUE)
plot(OutliersWGS_NoDevNY[OutliersWGS_NoDevNY$out_loc_hi == 1,], pch = 16, col = colFun(OutliersWGS_NoDevNY$Qs[OutliersWGS_NoDevNY$out_loc_hi == 1]), cex = 0.7, add = TRUE)
box()
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('topleft', title = expression(paste('Q'[s], ' (mW/m'^2, ')')), legend = c('< 30', '30 - 50', '50 - 70', '70 - 90', '>= 90', 'Not Tested', 'Not Outlier'), col = c(colFun(c(20,40,60,80,100)), 'purple', 'black'), pch = c(16,16,16,16,16,17,16))
dev.off()

png('LoHiOuts_Map_RmOps.png', res = 1200, units = 'in', width = 14, height = 7)
layout(cbind(1,2))
par(xaxs = 'i', yaxs = 'i', mar = c(2,3,3,1))
plot(WellsDeepWGS_RmOps, col = 'white', pch = 16, cex = 0.2, main = 'Low Outliers', cex.main = 2)
plot(Counties[Counties$STATEFP %in% c(42,36,54,24,21,51),], border = 'grey', add=TRUE)
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_RmOps[NotOutliersWGS_RmOps$out_loc_error == 0,], pch = 16, cex = 0.2, add = TRUE)
plot(NotOutliersWGS_RmOps[NotOutliersWGS_RmOps$out_loc_error == 1,], pch = 17, col = 'purple', cex = 0.4, add = TRUE)
plot(OutliersWGS_RmOps[OutliersWGS_RmOps$out_loc_lo == 1,], pch = 16, col = colFun(OutliersWGS_RmOps$Qs[OutliersWGS_RmOps$out_loc_lo == 1]), cex = 0.7, add = TRUE)
box()
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('topleft', title = expression(paste('Q'[s], ' (mW/m'^2, ')')), legend = c('< 30', '30 - 50', '50 - 70', '70 - 90', '>= 90', 'Not Tested', 'Not Outlier'), col = c(colFun(c(20,40,60,80,100)), 'purple', 'black'), pch = c(16,16,16,16,16,17,16))

par(xaxs = 'i', yaxs = 'i', mar = c(2,3,3,1))
plot(WellsDeepWGS_RmOps, col = 'white', pch = 16, cex = 0.2, main = 'High Outliers', cex.main = 2)
plot(Counties[Counties$STATEFP %in% c(42,36,54,24,21,51),], border = 'grey', add=TRUE)
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_RmOps[NotOutliersWGS_RmOps$out_loc_error == 0,], pch = 16, cex = 0.2, add = TRUE)
plot(NotOutliersWGS_RmOps[NotOutliersWGS_RmOps$out_loc_error == 1,], pch = 17, col = 'purple', cex = 0.4, add = TRUE)
plot(OutliersWGS_RmOps[OutliersWGS_RmOps$out_loc_hi == 1,], pch = 16, col = colFun(OutliersWGS_RmOps$Qs[OutliersWGS_RmOps$out_loc_hi == 1]), cex = 0.7, add = TRUE)
box()
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('topleft', title = expression(paste('Q'[s], ' (mW/m'^2, ')')), legend = c('< 30', '30 - 50', '50 - 70', '70 - 90', '>= 90', 'Not Tested', 'Not Outlier'), col = c(colFun(c(20,40,60,80,100)), 'purple', 'black'), pch = c(16,16,16,16,16,17,16))
dev.off()


#   Map of the depth ranks of outliers ----
scaleRange = c(1,25)
scaleBy = 5
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)
Pal = c(Pal[-3], Pal[length(Pal)])
sets = rbind(c(1,3), c(2,4))
png("OutlierRankMap_Counties.png", width=8, height=8, units="in", res=600)
layout(sets)
par(mar=c(4.5,4.5,2,1))
par(xaxs = 'i', yaxs = 'i')
hist(Outs$out_loc_drank[which(Outs$out_loc_lo == 1)]*25, breaks = seq(0,25,1), main = 'Low Outliers', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=1.5, cex.axis = 1.5)
minor.tick(nx = 5, ny = 1, tick.ratio = 0.5)
hist(Outs$out_loc_drank[which(Outs$out_loc_lo == 0)]*25, breaks = seq(0,25,1), main = 'High Outliers', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=1.5, cex.axis = 1.5)
minor.tick(nx = 5, ny = 1, tick.ratio = 0.5)
#Map Low
par(mar=c(2,2,1,1))
plot(Outs, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
#plot(NoOuts, pch = 16, add = TRUE, cex=0.5)
plot(Outs[which(Outs$out_loc_lo == 1),][order(Outs$WellDepth[which(Outs$out_loc_lo == 1)], decreasing = TRUE), ],
     pch = 16, 
     col = colFun(Outs[which(Outs$out_loc_lo == 1),][order(Outs$WellDepth[which(Outs$out_loc_lo == 1)], decreasing = TRUE), ]$out_loc_drank*25), 
     cex = 0.5, add = TRUE)
north.arrow(-83, 39.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('topleft', title = 'Depth Rank', legend = c('1 - 5', '6 - 10', '11 - 15', '16 - 20', '21 - 25'), pch = 16, cex = 1.5, col = colFun(c(1, 7, 12, 17, 22)), bty = 'n')
#Map High
plot(Outs, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
#plot(NoOuts, pch = 16, add = TRUE, cex=0.5)
plot(Outs[which(Outs$out_loc_lo == 0),][order(Outs$WellDepth[which(Outs$out_loc_lo == 0)], decreasing = TRUE), ],
     pch = 16, 
     col = colFun(Outs[which(Outs$out_loc_lo == 0),][order(Outs$WellDepth[which(Outs$out_loc_lo == 0)], decreasing = TRUE), ]$out_loc_drank*25), 
     cex = .5, add = TRUE)
north.arrow(-83, 39.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('topleft', title = 'Depth Rank', legend = c('1 - 5', '6 - 10', '11 - 15', '16 - 20', '21 - 25'), pch = 16, cex = 1.5, col = colFun(c(1, 7, 12, 17, 22)), bty = 'n')
dev.off()

png("OutlierRankMap_Counties_RmOps.png", width=8, height=8, units="in", res=600)
layout(sets)
par(mar=c(4.5,4.5,2,1))
par(xaxs = 'i', yaxs = 'i')
hist(OutliersWGS_RmOps$out_loc_drank[which(OutliersWGS_RmOps$out_loc_lo == 1)]*25, breaks = seq(0,25,1), main = 'Low Outliers', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=1.5, cex.axis = 1.5)
minor.tick(nx = 5, ny = 1, tick.ratio = 0.5)
hist(OutliersWGS_RmOps$out_loc_drank[which(OutliersWGS_RmOps$out_loc_lo == 0)]*25, breaks = seq(0,25,1), main = 'High Outliers', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=1.5, cex.axis = 1.5)
minor.tick(nx = 5, ny = 1, tick.ratio = 0.5)
#Map Low
par(mar=c(2,2,1,1))
plot(OutliersWGS_RmOps, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
#plot(NoOuts, pch = 16, add = TRUE, cex=0.5)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$out_loc_lo == 1),][order(OutliersWGS_RmOps$WellDepth[which(OutliersWGS_RmOps$out_loc_lo == 1)], decreasing = TRUE), ],
     pch = 16, 
     col = colFun(OutliersWGS_RmOps[which(OutliersWGS_RmOps$out_loc_lo == 1),][order(OutliersWGS_RmOps$WellDepth[which(OutliersWGS_RmOps$out_loc_lo == 1)], decreasing = TRUE), ]$out_loc_drank*25), 
     cex = 0.5, add = TRUE)
north.arrow(-83, 39.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('topleft', title = 'Depth Rank', legend = c('1 - 5', '6 - 10', '11 - 15', '16 - 20', '21 - 25'), pch = 16, cex = 1.5, col = colFun(c(1, 7, 12, 17, 22)), bty = 'n')
#Map High
plot(OutliersWGS_RmOps, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
#plot(NoOuts, pch = 16, add = TRUE, cex=0.5)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$out_loc_lo == 0),][order(OutliersWGS_RmOps$WellDepth[which(OutliersWGS_RmOps$out_loc_lo == 0)], decreasing = TRUE), ],
     pch = 16, 
     col = colFun(OutliersWGS_RmOps[which(OutliersWGS_RmOps$out_loc_lo == 0),][order(OutliersWGS_RmOps$WellDepth[which(OutliersWGS_RmOps$out_loc_lo == 0)], decreasing = TRUE), ]$out_loc_drank*25), 
     cex = .5, add = TRUE)
north.arrow(-83, 39.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('topleft', title = 'Depth Rank', legend = c('1 - 5', '6 - 10', '11 - 15', '16 - 20', '21 - 25'), pch = 16, cex = 1.5, col = colFun(c(1, 7, 12, 17, 22)), bty = 'n')
dev.off()

png("OutlierRankMap_Counties_NoDevNY.png", width=8, height=8, units="in", res=600)
layout(sets)
par(mar=c(4.5,4.5,2,1))
par(xaxs = 'i', yaxs = 'i')
hist(OutliersWGS_NoDevNY$out_loc_drank[which(OutliersWGS_NoDevNY$out_loc_lo == 1)]*25, breaks = seq(0,25,1), main = 'Low Outliers', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=1.5, cex.axis = 1.5)
minor.tick(nx = 5, ny = 1, tick.ratio = 0.5)
hist(OutliersWGS_NoDevNY$out_loc_drank[which(OutliersWGS_NoDevNY$out_loc_lo == 0)]*25, breaks = seq(0,25,1), main = 'High Outliers', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=1.5, cex.axis = 1.5)
minor.tick(nx = 5, ny = 1, tick.ratio = 0.5)
#Map Low
par(mar=c(2,2,1,1))
plot(OutliersWGS_NoDevNY, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
#plot(NoOuts, pch = 16, add = TRUE, cex=0.5)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$out_loc_lo == 1),][order(OutliersWGS_NoDevNY$WellDepth[which(OutliersWGS_NoDevNY$out_loc_lo == 1)], decreasing = TRUE), ],
     pch = 16, 
     col = colFun(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$out_loc_lo == 1),][order(OutliersWGS_NoDevNY$WellDepth[which(OutliersWGS_NoDevNY$out_loc_lo == 1)], decreasing = TRUE), ]$out_loc_drank*25), 
     cex = 0.5, add = TRUE)
north.arrow(-83, 39.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('topleft', title = 'Depth Rank', legend = c('1 - 5', '6 - 10', '11 - 15', '16 - 20', '21 - 25'), pch = 16, cex = 1.5, col = colFun(c(1, 7, 12, 17, 22)), bty = 'n')
#Map High
plot(OutliersWGS_NoDevNY, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
#plot(NoOuts, pch = 16, add = TRUE, cex=0.5)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$out_loc_lo == 0),][order(OutliersWGS_NoDevNY$WellDepth[which(OutliersWGS_NoDevNY$out_loc_lo == 0)], decreasing = TRUE), ],
     pch = 16, 
     col = colFun(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$out_loc_lo == 0),][order(OutliersWGS_NoDevNY$WellDepth[which(OutliersWGS_NoDevNY$out_loc_lo == 0)], decreasing = TRUE), ]$out_loc_drank*25), 
     cex = .5, add = TRUE)
north.arrow(-83, 39.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('topleft', title = 'Depth Rank', legend = c('1 - 5', '6 - 10', '11 - 15', '16 - 20', '21 - 25'), pch = 16, cex = 1.5, col = colFun(c(1, 7, 12, 17, 22)), bty = 'n')
dev.off()


#   Depth Rank Empircal CDFs and KS Test ----
png('LowOutliers_wUniform.png', res = 300, units = 'in', width = 6, height = 6)
hist(Outs$out_loc_drank[which(Outs$out_loc_lo == 1)]*25, breaks = seq(0,25,1), main = 'Low Outliers Depth Rank', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=2, cex.axis = 1.5, ylim = c(0,35))
lines(c(0,25), c(length(Outs$out_loc_drank[which(Outs$out_loc_lo == 1)])/25,length(Outs$out_loc_drank[which(Outs$out_loc_lo == 1)])/25), lwd = 2, col='blue')
dev.off()

hist(OutliersWGS_RmOps$out_loc_drank[which(OutliersWGS_RmOps$out_loc_lo == 1)]*25, breaks = seq(0,25,1), main = 'Low Outliers Depth Rank', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=2, cex.axis = 1.5, ylim = c(0,35))
lines(c(0,25), c(length(OutliersWGS_RmOps$out_loc_drank[which(OutliersWGS_RmOps$out_loc_lo == 1)])/25,length(OutliersWGS_RmOps$out_loc_drank[which(OutliersWGS_RmOps$out_loc_lo == 1)])/25), lwd = 2, col='blue')

hist(OutliersWGS_NoDevNY$out_loc_drank[which(OutliersWGS_NoDevNY$out_loc_lo == 1)]*25, breaks = seq(0,25,1), main = 'Low Outliers Depth Rank', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=2, cex.axis = 1.5, ylim = c(0,35))
lines(c(0,25), c(length(OutliersWGS_NoDevNY$out_loc_drank[which(OutliersWGS_NoDevNY$out_loc_lo == 1)])/25,length(OutliersWGS_NoDevNY$out_loc_drank[which(OutliersWGS_NoDevNY$out_loc_lo == 1)])/25), lwd = 2, col='blue')

png('HighOutliers_wUniform.png', res = 300, units = 'in', width = 6, height = 6)
hist(Outs$out_loc_drank[which(Outs$out_loc_lo == 0)]*25, breaks = seq(0,25,1), main = 'High Outliers Depth Rank', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=2, cex.axis = 1.5, ylim = c(0,35))
lines(c(0,25), c(length(Outs$out_loc_drank[which(Outs$out_loc_lo == 0)])/25,length(Outs$out_loc_drank[which(Outs$out_loc_lo == 0)])/25), lwd = 2, col='red')
dev.off()

hist(OutliersWGS_RmOps$out_loc_drank[which(OutliersWGS_RmOps$out_loc_lo == 0)]*25, breaks = seq(0,25,1), main = 'High Outliers Depth Rank', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=2, cex.axis = 1.5, ylim = c(0,35))
lines(c(0,25), c(length(OutliersWGS_RmOps$out_loc_drank[which(OutliersWGS_RmOps$out_loc_lo == 0)])/25,length(OutliersWGS_RmOps$out_loc_drank[which(OutliersWGS_RmOps$out_loc_lo == 0)])/25), lwd = 2, col='red')

hist(OutliersWGS_NoDevNY$out_loc_drank[which(OutliersWGS_NoDevNY$out_loc_lo == 0)]*25, breaks = seq(0,25,1), main = 'High Outliers Depth Rank', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=2, cex.axis = 1.5, ylim = c(0,35))
lines(c(0,25), c(length(OutliersWGS_NoDevNY$out_loc_drank[which(OutliersWGS_NoDevNY$out_loc_lo == 0)])/25,length(OutliersWGS_NoDevNY$out_loc_drank[which(OutliersWGS_NoDevNY$out_loc_lo == 0)])/25), lwd = 2, col='red')

hist(Outs$out_loc_drank*25, breaks = seq(0,25,1), main = 'All Outliers Depth Rank', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=2, cex.axis = 1.5)

#Create an empirical CDF
OutHiRank = Outs$out_loc_drank[which(Outs$out_loc_lo == 0)]*25
OutLoRank = Outs$out_loc_drank[which(Outs$out_loc_lo == 1)]*25

CumLo = cumsum(hist(OutLoRank, breaks = seq(0.5,25,.5), plot = FALSE)$counts/sum(hist(OutLoRank, breaks = seq(0.5,25,.5), plot = FALSE)$counts))
CumHi = cumsum(hist(OutHiRank, breaks = seq(0.5,25,.5), plot = FALSE)$counts/sum(hist(OutHiRank, breaks = seq(0.5,25,.5), plot = FALSE)$counts))

#Make Step Function
LoFun = stepfun(x = seq(1,25,.5), y = c(0,CumLo))
HiFun = stepfun(x = seq(1,25,.5), y = c(0,CumHi))
UnifFun = stepfun(x = seq(1,25,0.5), y = c(seq(0,1,0.5/24.5)), right = FALSE, f = 0)

#Note that the plotting characters are CLOSED circles where they are plotted, even though the (confusing) default is open circles.
png('KS_CDFs.png', res = 300, width = 5, height = 5, units = 'in')
par(xaxs = 'i', yaxs = 'i')
plot(LoFun, pch = NA, col = 'blue', lwd = 2, xlab = 'Depth Rank', ylab = 'Cumulative Frequency', xlim = c(0,25), ylim = c(0,1), cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, main = 'Empirical CDFs for K-S Test')
par(new = TRUE)
plot(HiFun, pch = NA, col = 'red', axes=FALSE, xlab='', ylab='', lwd = 2, xlim = c(0,25), ylim = c(0,1), main = '')
par(new = TRUE)
plot(UnifFun, pch = NA, col = 'black', axes=FALSE, xlab='', ylab='', lwd = 2, xlim = c(0,25), ylim = c(0,1), main = '')
legend('topleft', legend = c('Discrete Uniform Distribution', 'Low Outliers', 'High Outliers'), lwd = 2, pch = NA, col = c('black', 'blue', 'red'))
minor.tick(nx=5,ny=2)
dev.off()

#KS Test Statistics
#Needs to use the uniform as decimals because of the == comparison in ecdf
UnifFun = stepfun(x = seq(1,25,0.5)/25, y = c(seq(0,1,0.5/24.5)), right = FALSE, f = 0)
KS_StatHi = max(abs(HiFun(knots(HiFun)) - UnifFun(knots(UnifFun))))
KS_StatLo = max(abs(LoFun(knots(LoFun)) - UnifFun(knots(UnifFun))))
KSHi = ks.test(x = Outs$out_loc_drank[which(Outs$out_loc_lo == 0)], y = UnifFun, simulate.p.value = TRUE, B = 10000)
KSLo = ks.test(x = Outs$out_loc_drank[which(Outs$out_loc_lo == 1)], y = UnifFun, simulate.p.value = TRUE, B = 10000)

#Both significant at 3% level. The high distribution is different from random at 1% level.

#Kuiper or Watson test:
Kuiper_StatHi = abs(max(0, HiFun(knots(HiFun)) - UnifFun(knots(UnifFun)))) + abs(min(0, HiFun(knots(HiFun)) - UnifFun(knots(UnifFun))))
Kuiper_StatLo = abs(max(0, LoFun(knots(LoFun)) - UnifFun(knots(UnifFun)))) + abs(min(0, LoFun(knots(LoFun)) - UnifFun(knots(UnifFun))))

#   Poisson Test for Depth Bins ----
#Arrival rate in outliers per well
#lambda = nrow(Outs)/nrow(AllData) 

#Outliers in each bin class
#OutBins = hist(Outs$WellDepth, breaks = c(500,seq(1000,2000,200),3000,7000), plot = FALSE)$counts
#Total data points in each bin class
#DatBins = hist(AllData$WellDepth, breaks = c(500,seq(1000,2000,200),3000,7000), plot = FALSE)$counts

#Poisson probabilities for each bin class
#nu = lambda*DatBins
#dpois(OutBins, nu)

#   Binomial and Chi-Squared Test for Depth Bins ----
#Better than Poisson - Less Contamination in estimation of lambda

#Total number of outliers per bin based on number of samples in the bin
DatBins = hist(AllData$WellDepth[-which(AllData$out_loc_error == 1)], breaks = c(750,seq(1000,3200,200),6600), plot = FALSE)$counts
OutBins = hist(Outs$WellDepth, breaks = c(750,seq(1000,3200,200),6600), plot = FALSE)$counts
BinomP_AllData = nrow(Outs)/length(which(AllData$out_loc_error == 0))
PerfBins = round(DatBins*BinomP_AllData,0)
PerfBins[3] = 286 #round .75 down to make number of outliers match observed number.
BinomP_Bins = PerfBins/DatBins

#Chi-Squared Test difference
ChiSq = sum((OutBins - PerfBins)^2/PerfBins)
dfChiSq = length(OutBins) - 1
pVal = 1-pchisq(ChiSq, dfChiSq)
#There is preferential sorting of outliers at 1% level. Try to find out why:

#Strength of the test using Cramer's V:
CramerV = sqrt(ChiSq/length(which(AllData$out_loc_error == 0))/min(1,(length(DatBins)-1)))
BiasCorrV = sqrt(max(0, ChiSq/length(which(AllData$out_loc_error == 0)) - 1*(length(DatBins)-1)/(length(which(AllData$out_loc_error == 0)) - 1))/
                   min(2 - 1/(length(which(AllData$out_loc_error == 0)) - 1) - 1, length(DatBins) - ((length(DatBins) - 1)^2)/(length(which(AllData$out_loc_error == 0)) - 1) - 1))
  
#Upper tail test greater than 2400 m - Not Significant p = 50%
ChiSq_UpTail = sum((OutBins[c(9,10,11,12,13)] - PerfBins[c(9,10,11,12,13)])^2/PerfBins[c(9,10,11,12,13)])
dfChiSq_UpTail = length(OutBins[c(9,10,11,12,13)]) - 1
pVal_UpTail = 1-pchisq(ChiSq_UpTail, dfChiSq_UpTail)

#Lower tail test less than 1400 m - Not Significant p = 25%
ChiSq_LowTail = sum((OutBins[c(1,2,3)] - PerfBins[c(1,2,3)])^2/PerfBins[c(1,2,3)])
dfChiSq_LowTail = length(OutBins[c(1,2,3)]) - 1
pVal_LowTail = 1-pchisq(ChiSq_LowTail, dfChiSq_LowTail)

#Binomial Test on Individual Bins - 7 and 8 are significant. Maybe there's a reason for this 2000 - 2400 m anomaly
binom.test(OutBins[1], DatBins[1], BinomP_Bins[1], alternative = 'greater') #Not Sig
binom.test(OutBins[2], DatBins[2], BinomP_Bins[2], alternative = 'greater') #Not Sig
binom.test(OutBins[3], DatBins[3], BinomP_Bins[3], alternative = 'less')    #Close to significant at 10 % level
binom.test(OutBins[4], DatBins[4], BinomP_Bins[4], alternative = 'less')    #Significant at 5% level
binom.test(OutBins[5], DatBins[5], BinomP_Bins[5], alternative = 'greater') #Not Sig
binom.test(OutBins[6], DatBins[6], BinomP_Bins[6], alternative = 'greater') #Not Sig
binom.test(OutBins[7], DatBins[7], BinomP_Bins[7], alternative = 'greater') #Significant at 2% level
binom.test(OutBins[8], DatBins[8], BinomP_Bins[8], alternative = 'greater') #Significant at .1% level
binom.test(OutBins[9], DatBins[9], BinomP_Bins[9], alternative = 'greater') #Close to significant at 10 % level
binom.test(OutBins[10], DatBins[10], BinomP_Bins[10], alternative = 'greater') #Not Sig
binom.test(OutBins[11], DatBins[11], BinomP_Bins[11], alternative = 'greater') #Not Sig
binom.test(OutBins[12], DatBins[12], BinomP_Bins[12], alternative = 'greater') #Not Sig
binom.test(OutBins[13], DatBins[13], BinomP_Bins[13], alternative = 'greater') #Not Sig

#Bin 7 and 8 (with Yates correction for n = 2) test significant at .01% level
ChiSq_2k = sum(abs((OutBins[c(7,8)] - PerfBins[c(7,8)]) - 0.5)^2/PerfBins[c(7,8)])
dfChiSq_2k = length(OutBins[c(7,8)]) - 1
pVal_2k = 1-pchisq(ChiSq_2k, dfChiSq_2k)

CramerV = sqrt(ChiSq_2k/length(which(AllData$out_loc_error == 0 & AllData$WellDepth >= 2000 & AllData$WellDepth < 2400))/min(1,1))
BiasCorrV = sqrt(max(0, ChiSq_2k/length(which(AllData$out_loc_error == 0 & AllData$WellDepth >= 2000 & AllData$WellDepth < 2400)) - 1/(length(which(AllData$out_loc_error == 0 & AllData$WellDepth >= 2000 & AllData$WellDepth < 2400)) - 1))/
                   min(2 - 1/(length(which(AllData$out_loc_error == 0 & AllData$WellDepth >= 2000 & AllData$WellDepth < 2400)) - 1) - 1, 2 - 1/(length(which(AllData$out_loc_error == 0 & AllData$WellDepth >= 2000 & AllData$WellDepth < 2400)) - 1) - 1))


#Histogram of Outliers and All Data in each bin
png('AllDatOutsHist_log.png', res = 300, width = 10, height = 6, units = 'in')
par(mar = c(4.5,4.5,1,1))
barplot(height = hist(AllData$WellDepth[-which(AllData$out_loc_error == 1)], breaks = c(750,seq(1000,3200,200),6600), plot = FALSE)$counts, width = c(250,rep(200,11),3400), space = c(750/mean(c(250,rep(200,11),3400)),rep(0,12)), col = 'grey80', xlim = c(750,7000), ylim = c(1,10000), xlab='Depth (m)', ylab='Frequency', log = 'y', axes = FALSE, cex.lab=1.5)
axis(side = 1, at = c(750, seq(1000,3200,200), 6600, 7000), cex.axis=1.5)
at.y <- outer(1:9, 10^(-7:7))
lab.y <- ifelse(log10(at.y) %% 1 == 0, sapply(log10(at.y), function(i) as.expression(bquote(10^ .(i)))), NA)
par(tcl = -0.2)
axis(side=2, at=at.y, labels=lab.y, cex.axis=1.5, las=1)
at.y <- 10^(-7:7)
par(tcl = -0.5)
axis(side=2, at=at.y, labels = FALSE, cex.axis=1.5, las=1)
par(new = TRUE)
barplot(height = PerfBins, width = c(250,rep(200,11),3400), space = c(750/mean(c(250,rep(200,11),3400)),rep(0,12)), col = NA, xlim = c(750,7000), ylim = c(1,10000), xlab='', ylab='', axes = FALSE, log = 'y', border = 'green')
par(new = TRUE)
barplot(height = hist(Outs$WellDepth, breaks = c(750,seq(1000,3200,200),6600), plot = FALSE)$counts, width = c(250,rep(200,11),3400), space = c(750/mean(c(250,rep(200,11),3400)),rep(0,12)), border = 'red', xlim = c(750,7000), ylim = c(1,10000), xlab='', ylab='', axes = FALSE, log = 'y', col = NA)
par(new = TRUE)
barplot(height = hist(AllData$WellDepth[-which(AllData$out_loc_error == 1)], breaks = c(750,seq(1000,3200,200),6600), plot = FALSE)$counts, width = c(250,rep(200,11),3400), space = c(750/mean(c(250,rep(200,11),3400)),rep(0,12)), col = NA, xlim = c(750,7000), ylim = c(1,10000), xlab='', ylab='', log = 'y', axes = FALSE)
legend('topright', legend=c('All Data', 'Expected Allocation', 'Actual Allocation'), col = c('grey', 'green', 'red'), pch = 15, cex = 1.5)
dev.off()

#   Maps of the depth slices for outliers ----

#See if these are primarily high or low outliers in the significantly different depth slice
hist(Outs$out_loc_lo[which(Outs$WellDepth >= 2000 & Outs$WellDepth < 2400)])
#Primarily high outliers. See where they are located. If clustered, may have local problem
plot(NoOuts[which(NoOuts$WellDepth >= 2000 & NoOuts$WellDepth < 2400 & NoOuts$out_loc_error == 0),], pch = 16)
plot(Outs[which(Outs$WellDepth >= 2000 & Outs$WellDepth < 2400 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue')
plot(Outs[which(Outs$WellDepth >= 2000 & Outs$WellDepth < 2400 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red')

#High outliers are clustered in PA
#Check the depth rank of these high points within their local neighborhood. Seems like if they are all the deepest points, there's a relationship.
#The PA Points matter the most. Others do not appear to have a cluster. Are these good or bad data points?
hist(Outs$out_loc_drank[which(Outs$WellDepth >= 2000 & Outs$WellDepth < 3000 & Outs$out_loc_lo == 0)]*25, breaks = seq(0,25,1))
#Overwhlemingly the deep wells are the outliers. Plot spatially
plot(Outs[which(Outs$WellDepth >= 2000 & Outs$WellDepth < 3000 & Outs$out_loc_lo == 0),], col = colFun(Outs$out_loc_drank[which(Outs$WellDepth >= 2000 & Outs$WellDepth < 3000 & Outs$out_loc_lo == 0)]*25), pch = 16)

#All ranked 21 or higher in SWPA. Are they unusually deep for their region?
boxplot(Outs$out_loc_rmrank[which(Outs$out_loc_lo == 1)] ~ as.factor(Outs$out_loc_drank[which(Outs$out_loc_lo == 1)]*25), cex.axis = 1.5, cex.main = 2)
boxplot(Outs$out_loc_rmrank[which(Outs$out_loc_lo == 0)] ~ as.factor(Outs$out_loc_drank[which(Outs$out_loc_lo == 0)]*25), cex.axis = 1.5, cex.main = 2)

#Fixme: Are any ranks clustered more than would be expected under a random process?
#Could also check spatial correlation of the ranks

# Make panel plot of depth horizons with low and high outliers colored in, as well as missing datapoints
sets = rbind(c(1,2,3,4),c(5,6,7,8))
png('OutsByDepth_Panels.png', res = 300, height = 8, width = 16, units = 'in')
layout(sets)
par(xaxs = 'i', yaxs = 'i', mar = c(2,2.5,1,1))
#Data from 750 - 1000 - 4 low outliers that are not clustered. Not problematic.
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border='grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NoOuts[which(NoOuts$WellDepth >= 750 & NoOuts$WellDepth < 1000 & NoOuts$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(Outs[which(Outs$WellDepth >= 750 & Outs$WellDepth < 1000 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(Outs[which(Outs$WellDepth >= 750 & Outs$WellDepth < 1000 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NoOuts[which(NoOuts$WellDepth >= 750 & NoOuts$WellDepth < 1000 & NoOuts$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '750 - 1000 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 1000 - 1200
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NoOuts[which(NoOuts$WellDepth >= 1000 & NoOuts$WellDepth < 1200 & NoOuts$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(Outs[which(Outs$WellDepth >= 1000 & Outs$WellDepth < 1200 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(Outs[which(Outs$WellDepth >= 1000 & Outs$WellDepth < 1200 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NoOuts[which(NoOuts$WellDepth >= 1000 & NoOuts$WellDepth < 1200 & NoOuts$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '1000 - 1200 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 1200 - 1400
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NoOuts[which(NoOuts$WellDepth >= 1200 & NoOuts$WellDepth < 1400 & NoOuts$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(Outs[which(Outs$WellDepth >= 1200 & Outs$WellDepth < 1400 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(Outs[which(Outs$WellDepth >= 1200 & Outs$WellDepth < 1400 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NoOuts[which(NoOuts$WellDepth >= 1200 & NoOuts$WellDepth < 1400 & NoOuts$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '1200 - 1400 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 1400 - 1600
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NoOuts[which(NoOuts$WellDepth >= 1400 & NoOuts$WellDepth < 1600 & NoOuts$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(Outs[which(Outs$WellDepth >= 1400 & Outs$WellDepth < 1600 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(Outs[which(Outs$WellDepth >= 1400 & Outs$WellDepth < 1600 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NoOuts[which(NoOuts$WellDepth >= 1400 & NoOuts$WellDepth < 1600 & NoOuts$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '1400 - 1600 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 1600 - 2000
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NoOuts[which(NoOuts$WellDepth >= 1600 & NoOuts$WellDepth < 2000 & NoOuts$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(Outs[which(Outs$WellDepth >= 1600 & Outs$WellDepth < 2000 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(Outs[which(Outs$WellDepth >= 1600 & Outs$WellDepth < 2000 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NoOuts[which(NoOuts$WellDepth >= 1600 & NoOuts$WellDepth < 2000 & NoOuts$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '1600 - 2000 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 2000 - 2400
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NoOuts[which(NoOuts$WellDepth >= 2000 & NoOuts$WellDepth < 2400 & NoOuts$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(Outs[which(Outs$WellDepth >= 2000 & Outs$WellDepth < 2400 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(Outs[which(Outs$WellDepth >= 2000 & Outs$WellDepth < 2400 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NoOuts[which(NoOuts$WellDepth >= 2000 & NoOuts$WellDepth < 2400 & NoOuts$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '2000 - 2400 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 2400 - 3000
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NoOuts[which(NoOuts$WellDepth >= 2400 & NoOuts$WellDepth < 3000 & NoOuts$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(Outs[which(Outs$WellDepth >= 2400 & Outs$WellDepth < 3000 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(Outs[which(Outs$WellDepth >= 2400 & Outs$WellDepth < 3000 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NoOuts[which(NoOuts$WellDepth >= 2400 & NoOuts$WellDepth < 3000 & NoOuts$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '2400 - 3000 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data deeper than 3000 m are 86% in New York (Not really clustered)
nrow(AllData[which(AllData$WellDepth >= 3000 & AllData$WellDepth < 6600 & AllData$State == 'NY'),])/nrow(AllData[which(AllData$WellDepth >= 3000 & AllData$WellDepth < 6600),])
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NoOuts[which(NoOuts$WellDepth >= 3000 & NoOuts$WellDepth < 6600 & NoOuts$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(Outs[which(Outs$WellDepth >= 3000 & Outs$WellDepth < 6600 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(Outs[which(Outs$WellDepth >= 3000 & Outs$WellDepth < 6600 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NoOuts[which(NoOuts$WellDepth >= 3000 & NoOuts$WellDepth < 6600 & NoOuts$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '3000 - 6600 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')
dev.off()

#Plot outliers in the depth range from 2 - 2.4 km
plot(NoOuts_max2k[which(NoOuts_max2k$WellDepth >= 2000 & NoOuts_max2k$WellDepth < 2400 & NoOuts_max2k$out_loc_error == 0),], pch = 16)
plot(Outs_max2k[which(Outs_max2k$WellDepth >= 2000 & Outs_max2k$WellDepth < 2400 & Outs_max2k$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue')
plot(Outs_max2k[which(Outs_max2k$WellDepth >= 2000 & Outs_max2k$WellDepth < 2400 & Outs_max2k$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red')
plot(NoOuts_max2k[which(NoOuts_max2k$WellDepth >= 2000 & NoOuts_max2k$WellDepth < 2400 & NoOuts_max2k$out_loc_error == 1),], pch = 17, col='purple',add=TRUE)
plot(NY, lwd = 2, add = TRUE)
plot(PA, lwd = 2, add = TRUE)
plot(WV, lwd = 2, add = TRUE)
plot(Counties, add = TRUE)

# Test Spatial Patterns in Outliers - by Depth. Doesn't seem useful.
# scaleRange = c(750,4200)
# scaleBy = 690
# Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)
# plot(NoOuts[order(NoOuts$WellDepth),], col = colFun(NoOuts$WellDepth[order(NoOuts$WellDepth)]), pch = 16, cex = 0.5)
# plot(Outs[order(Outs$WellDepth),], col = colFun(Outs$WellDepth[order(Outs$WellDepth)]), pch = 15, cex = 0.5, add = TRUE)
# plot(NY, lwd = 2, add=TRUE)
# plot(PA, lwd = 2, add=TRUE)
# plot(WV, lwd = 2, add=TRUE)
# plot(MD, lwd = 2, add=TRUE)
# plot(KY, lwd = 2, add=TRUE)
# plot(VA, lwd = 2, add=TRUE)

png('OutsByDepth_Panels_RmOps.png', res = 300, height = 8, width = 16, units = 'in')
layout(sets)
par(xaxs = 'i', yaxs = 'i', mar = c(2,2.5,1,1))
#Data from 750 - 1000 - 4 low outliers that are not clustered. Not problematic.
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border='grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_RmOps[which(NotOutliersWGS_RmOps$WellDepth >= 750 & NotOutliersWGS_RmOps$WellDepth < 1000 & NotOutliersWGS_RmOps$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$WellDepth >= 750 & OutliersWGS_RmOps$WellDepth < 1000 & OutliersWGS_RmOps$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$WellDepth >= 750 & OutliersWGS_RmOps$WellDepth < 1000 & OutliersWGS_RmOps$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NotOutliersWGS_RmOps[which(NotOutliersWGS_RmOps$WellDepth >= 750 & NotOutliersWGS_RmOps$WellDepth < 1000 & NotOutliersWGS_RmOps$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '750 - 1000 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 1000 - 1200
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_RmOps[which(NotOutliersWGS_RmOps$WellDepth >= 1000 & NotOutliersWGS_RmOps$WellDepth < 1200 & NotOutliersWGS_RmOps$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$WellDepth >= 1000 & OutliersWGS_RmOps$WellDepth < 1200 & OutliersWGS_RmOps$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$WellDepth >= 1000 & OutliersWGS_RmOps$WellDepth < 1200 & OutliersWGS_RmOps$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NotOutliersWGS_RmOps[which(NotOutliersWGS_RmOps$WellDepth >= 1000 & NotOutliersWGS_RmOps$WellDepth < 1200 & NotOutliersWGS_RmOps$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '1000 - 1200 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 1200 - 1400
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_RmOps[which(NotOutliersWGS_RmOps$WellDepth >= 1200 & NotOutliersWGS_RmOps$WellDepth < 1400 & NotOutliersWGS_RmOps$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$WellDepth >= 1200 & OutliersWGS_RmOps$WellDepth < 1400 & OutliersWGS_RmOps$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$WellDepth >= 1200 & OutliersWGS_RmOps$WellDepth < 1400 & OutliersWGS_RmOps$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NotOutliersWGS_RmOps[which(NotOutliersWGS_RmOps$WellDepth >= 1200 & NotOutliersWGS_RmOps$WellDepth < 1400 & NotOutliersWGS_RmOps$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '1200 - 1400 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 1400 - 1600
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_RmOps[which(NotOutliersWGS_RmOps$WellDepth >= 1400 & NotOutliersWGS_RmOps$WellDepth < 1600 & NotOutliersWGS_RmOps$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$WellDepth >= 1400 & OutliersWGS_RmOps$WellDepth < 1600 & OutliersWGS_RmOps$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$WellDepth >= 1400 & OutliersWGS_RmOps$WellDepth < 1600 & OutliersWGS_RmOps$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NotOutliersWGS_RmOps[which(NotOutliersWGS_RmOps$WellDepth >= 1400 & NotOutliersWGS_RmOps$WellDepth < 1600 & NotOutliersWGS_RmOps$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '1400 - 1600 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 1600 - 2000
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_RmOps[which(NotOutliersWGS_RmOps$WellDepth >= 1600 & NotOutliersWGS_RmOps$WellDepth < 2000 & NotOutliersWGS_RmOps$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$WellDepth >= 1600 & OutliersWGS_RmOps$WellDepth < 2000 & OutliersWGS_RmOps$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$WellDepth >= 1600 & OutliersWGS_RmOps$WellDepth < 2000 & OutliersWGS_RmOps$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NotOutliersWGS_RmOps[which(NotOutliersWGS_RmOps$WellDepth >= 1600 & NotOutliersWGS_RmOps$WellDepth < 2000 & NotOutliersWGS_RmOps$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '1600 - 2000 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 2000 - 2400
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_RmOps[which(NotOutliersWGS_RmOps$WellDepth >= 2000 & NotOutliersWGS_RmOps$WellDepth < 2400 & NotOutliersWGS_RmOps$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$WellDepth >= 2000 & OutliersWGS_RmOps$WellDepth < 2400 & OutliersWGS_RmOps$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$WellDepth >= 2000 & OutliersWGS_RmOps$WellDepth < 2400 & OutliersWGS_RmOps$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NotOutliersWGS_RmOps[which(NotOutliersWGS_RmOps$WellDepth >= 2000 & NotOutliersWGS_RmOps$WellDepth < 2400 & NotOutliersWGS_RmOps$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '2000 - 2400 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 2400 - 3000
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_RmOps[which(NotOutliersWGS_RmOps$WellDepth >= 2400 & NotOutliersWGS_RmOps$WellDepth < 3000 & NotOutliersWGS_RmOps$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$WellDepth >= 2400 & OutliersWGS_RmOps$WellDepth < 3000 & OutliersWGS_RmOps$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$WellDepth >= 2400 & OutliersWGS_RmOps$WellDepth < 3000 & OutliersWGS_RmOps$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NotOutliersWGS_RmOps[which(NotOutliersWGS_RmOps$WellDepth >= 2400 & NotOutliersWGS_RmOps$WellDepth < 3000 & NotOutliersWGS_RmOps$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '2400 - 3000 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data deeper than 3000 m are 86% in New York (Not really clustered)
nrow(AllData[which(AllData$WellDepth >= 3000 & AllData$WellDepth < 6600 & AllData$State == 'NY'),])/nrow(AllData[which(AllData$WellDepth >= 3000 & AllData$WellDepth < 6600),])
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_RmOps[which(NotOutliersWGS_RmOps$WellDepth >= 3000 & NotOutliersWGS_RmOps$WellDepth < 6600 & NotOutliersWGS_RmOps$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$WellDepth >= 3000 & OutliersWGS_RmOps$WellDepth < 6600 & OutliersWGS_RmOps$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(OutliersWGS_RmOps[which(OutliersWGS_RmOps$WellDepth >= 3000 & OutliersWGS_RmOps$WellDepth < 6600 & OutliersWGS_RmOps$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NotOutliersWGS_RmOps[which(NotOutliersWGS_RmOps$WellDepth >= 3000 & NotOutliersWGS_RmOps$WellDepth < 6600 & NotOutliersWGS_RmOps$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '3000 - 6600 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')
dev.off()

png('OutsByDepth_Panels_NoDevNY.png', res = 300, height = 8, width = 16, units = 'in')
layout(sets)
par(xaxs = 'i', yaxs = 'i', mar = c(2,2.5,1,1))
#Data from 750 - 1000 - 4 low outliers that are not clustered. Not problematic.
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border='grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_NoDevNY[which(NotOutliersWGS_NoDevNY$WellDepth >= 750 & NotOutliersWGS_NoDevNY$WellDepth < 1000 & NotOutliersWGS_NoDevNY$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$WellDepth >= 750 & OutliersWGS_NoDevNY$WellDepth < 1000 & OutliersWGS_NoDevNY$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$WellDepth >= 750 & OutliersWGS_NoDevNY$WellDepth < 1000 & OutliersWGS_NoDevNY$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NotOutliersWGS_NoDevNY[which(NotOutliersWGS_NoDevNY$WellDepth >= 750 & NotOutliersWGS_NoDevNY$WellDepth < 1000 & NotOutliersWGS_NoDevNY$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '750 - 1000 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 1000 - 1200
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_NoDevNY[which(NotOutliersWGS_NoDevNY$WellDepth >= 1000 & NotOutliersWGS_NoDevNY$WellDepth < 1200 & NotOutliersWGS_NoDevNY$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$WellDepth >= 1000 & OutliersWGS_NoDevNY$WellDepth < 1200 & OutliersWGS_NoDevNY$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$WellDepth >= 1000 & OutliersWGS_NoDevNY$WellDepth < 1200 & OutliersWGS_NoDevNY$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NotOutliersWGS_NoDevNY[which(NotOutliersWGS_NoDevNY$WellDepth >= 1000 & NotOutliersWGS_NoDevNY$WellDepth < 1200 & NotOutliersWGS_NoDevNY$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '1000 - 1200 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 1200 - 1400
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_NoDevNY[which(NotOutliersWGS_NoDevNY$WellDepth >= 1200 & NotOutliersWGS_NoDevNY$WellDepth < 1400 & NotOutliersWGS_NoDevNY$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$WellDepth >= 1200 & OutliersWGS_NoDevNY$WellDepth < 1400 & OutliersWGS_NoDevNY$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$WellDepth >= 1200 & OutliersWGS_NoDevNY$WellDepth < 1400 & OutliersWGS_NoDevNY$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NotOutliersWGS_NoDevNY[which(NotOutliersWGS_NoDevNY$WellDepth >= 1200 & NotOutliersWGS_NoDevNY$WellDepth < 1400 & NotOutliersWGS_NoDevNY$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '1200 - 1400 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 1400 - 1600
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_NoDevNY[which(NotOutliersWGS_NoDevNY$WellDepth >= 1400 & NotOutliersWGS_NoDevNY$WellDepth < 1600 & NotOutliersWGS_NoDevNY$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$WellDepth >= 1400 & OutliersWGS_NoDevNY$WellDepth < 1600 & OutliersWGS_NoDevNY$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$WellDepth >= 1400 & OutliersWGS_NoDevNY$WellDepth < 1600 & OutliersWGS_NoDevNY$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NotOutliersWGS_NoDevNY[which(NotOutliersWGS_NoDevNY$WellDepth >= 1400 & NotOutliersWGS_NoDevNY$WellDepth < 1600 & NotOutliersWGS_NoDevNY$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '1400 - 1600 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 1600 - 2000
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_NoDevNY[which(NotOutliersWGS_NoDevNY$WellDepth >= 1600 & NotOutliersWGS_NoDevNY$WellDepth < 2000 & NotOutliersWGS_NoDevNY$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$WellDepth >= 1600 & OutliersWGS_NoDevNY$WellDepth < 2000 & OutliersWGS_NoDevNY$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$WellDepth >= 1600 & OutliersWGS_NoDevNY$WellDepth < 2000 & OutliersWGS_NoDevNY$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NotOutliersWGS_NoDevNY[which(NotOutliersWGS_NoDevNY$WellDepth >= 1600 & NotOutliersWGS_NoDevNY$WellDepth < 2000 & NotOutliersWGS_NoDevNY$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '1600 - 2000 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 2000 - 2400
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_NoDevNY[which(NotOutliersWGS_NoDevNY$WellDepth >= 2000 & NotOutliersWGS_NoDevNY$WellDepth < 2400 & NotOutliersWGS_NoDevNY$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$WellDepth >= 2000 & OutliersWGS_NoDevNY$WellDepth < 2400 & OutliersWGS_NoDevNY$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$WellDepth >= 2000 & OutliersWGS_NoDevNY$WellDepth < 2400 & OutliersWGS_NoDevNY$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NotOutliersWGS_NoDevNY[which(NotOutliersWGS_NoDevNY$WellDepth >= 2000 & NotOutliersWGS_NoDevNY$WellDepth < 2400 & NotOutliersWGS_NoDevNY$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '2000 - 2400 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data from 2400 - 3000
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_NoDevNY[which(NotOutliersWGS_NoDevNY$WellDepth >= 2400 & NotOutliersWGS_NoDevNY$WellDepth < 3000 & NotOutliersWGS_NoDevNY$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$WellDepth >= 2400 & OutliersWGS_NoDevNY$WellDepth < 3000 & OutliersWGS_NoDevNY$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$WellDepth >= 2400 & OutliersWGS_NoDevNY$WellDepth < 3000 & OutliersWGS_NoDevNY$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NotOutliersWGS_NoDevNY[which(NotOutliersWGS_NoDevNY$WellDepth >= 2400 & NotOutliersWGS_NoDevNY$WellDepth < 3000 & NotOutliersWGS_NoDevNY$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '2400 - 3000 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')

#Data deeper than 3000 m are 86% in New York (Not really clustered)
nrow(AllData[which(AllData$WellDepth >= 3000 & AllData$WellDepth < 6600 & AllData$State == 'NY'),])/nrow(AllData[which(AllData$WellDepth >= 3000 & AllData$WellDepth < 6600),])
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE, border = 'grey')
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NotOutliersWGS_NoDevNY[which(NotOutliersWGS_NoDevNY$WellDepth >= 3000 & NotOutliersWGS_NoDevNY$WellDepth < 6600 & NotOutliersWGS_NoDevNY$out_loc_error == 0),], pch = 16, add = TRUE, cex = 0.7)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$WellDepth >= 3000 & OutliersWGS_NoDevNY$WellDepth < 6600 & OutliersWGS_NoDevNY$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue', cex = 0.7)
plot(OutliersWGS_NoDevNY[which(OutliersWGS_NoDevNY$WellDepth >= 3000 & OutliersWGS_NoDevNY$WellDepth < 6600 & OutliersWGS_NoDevNY$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red', cex = 0.7)
plot(NotOutliersWGS_NoDevNY[which(NotOutliersWGS_NoDevNY$WellDepth >= 3000 & NotOutliersWGS_NoDevNY$WellDepth < 6600 & NotOutliersWGS_NoDevNY$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple', cex = 0.7)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
text(x = -81.2, y = 43.2, '3000 - 6600 m', cex = 2)
legend('topleft', title = ' ', legend = c('Not Outlier', 'Low Outlier', 'High Outlier', 'Not Tested'), pch = c(16,16,16,17), cex = 1.3, col = c('black','blue', 'red', 'purple'), bty = 'n')
dev.off()


#   Q-Q plot for the points that are not outliers ----
#Q-Q plots for data in each interpolation region
#Unique axes
sets = rbind(c(1,2,3), c(4,5,6),c(7,8,9))
png('QQPlotHeatFlow_NotTestedOuts.png', res=600, units='in', width=10, height=10)
layout(sets)
par(mar=c(4,5,3,2))
qqnorm(CT@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Chautauqua, NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='red')
qqline(CT@data$Qs)
qqnorm(WPA@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='orange', xlim = c(-3.5,3.5), ylim = c(10,70))
qqline(WPA@data$Qs)
temp = qqnorm(WPA@data$Qs, plot.it = FALSE)$x[which(WPA@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = WPA@data$Qs[which(WPA@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(10,70))
qqnorm(NWPANY@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Northwestern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='yellow', xlim = c(-3,3), ylim = c(5,75))
qqline(NWPANY@data$Qs)
temp = qqnorm(NWPANY@data$Qs, plot.it = FALSE)$x[which(NWPANY@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = NWPANY@data$Qs[which(NWPANY@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(5,75))
qqnorm(CNY@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='green', xlim = c(-2.75,2.75), ylim = c(35,65))
qqline(CNY@data$Qs)
temp = qqnorm(CNY@data$Qs, plot.it = FALSE)$x[which(CNY@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CNY@data$Qs[which(CNY@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-2.75,2.75), ylim = c(35,65))
qqnorm(ENY@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='springgreen', xlim = c(-3,3), ylim = c(30,75))
qqline(ENY@data$Qs)
temp = qqnorm(ENY@data$Qs, plot.it = FALSE)$x[which(ENY@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENY@data$Qs[which(ENY@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(30,75))
qqnorm(ENYPA@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='skyblue', xlim = c(-3,3), ylim = c(10,110))
qqline(ENYPA@data$Qs)
temp = qqnorm(ENYPA@data$Qs, plot.it = FALSE)$x[which(ENYPA@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENYPA@data$Qs[which(ENYPA@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(10,110))
qqnorm(SWPA@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Southwestern PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='blue', xlim = c(-4,4), ylim = c(5,110))
qqline(SWPA@data$Qs)
temp = qqnorm(SWPA@data$Qs, plot.it = FALSE)$x[which(SWPA@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = SWPA@data$Qs[which(SWPA@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,110))
qqnorm(CWV@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='purple', xlim = c(-4,4), ylim = c(5,120))
qqline(CWV@data$Qs)
temp = qqnorm(CWV@data$Qs, plot.it = FALSE)$x[which(CWV@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CWV@data$Qs[which(CWV@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,120))
qqnorm(MT@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='violet', xlim = c(-3.5,3.5), ylim = c(5,120))
qqline(MT@data$Qs)
temp = qqnorm(MT@data$Qs, plot.it = FALSE)$x[which(MT@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = MT@data$Qs[which(MT@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(5,120))
dev.off()

png('QQPlotHeatFlow_NotTestedOuts_RmOps.png', res=600, units='in', width=10, height=10)
layout(sets)
par(mar=c(4,5,3,2))
qqnorm(CT_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Chautauqua, NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='red')
qqline(CT_RmOps@data$Qs)
qqnorm(WPA_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='orange', xlim = c(-3.5,3.5), ylim = c(10,70))
qqline(WPA_RmOps@data$Qs)
temp = qqnorm(WPA_RmOps@data$Qs, plot.it = FALSE)$x[which(WPA_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = WPA_RmOps@data$Qs[which(WPA_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(10,70))
qqnorm(NWPANY_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Northwestern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='yellow', xlim = c(-3,3), ylim = c(5,75))
qqline(NWPANY_RmOps@data$Qs)
temp = qqnorm(NWPANY_RmOps@data$Qs, plot.it = FALSE)$x[which(NWPANY_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = NWPANY_RmOps@data$Qs[which(NWPANY_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(5,75))
qqnorm(CNY_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='green', xlim = c(-2.75,2.75), ylim = c(35,65))
qqline(CNY_RmOps@data$Qs)
temp = qqnorm(CNY_RmOps@data$Qs, plot.it = FALSE)$x[which(CNY_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CNY_RmOps@data$Qs[which(CNY_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-2.75,2.75), ylim = c(35,65))
qqnorm(ENY_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='springgreen', xlim = c(-3,3), ylim = c(30,75))
qqline(ENY_RmOps@data$Qs)
temp = qqnorm(ENY_RmOps@data$Qs, plot.it = FALSE)$x[which(ENY_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENY_RmOps@data$Qs[which(ENY_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(30,75))
qqnorm(ENYPA_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='skyblue', xlim = c(-3,3), ylim = c(10,110))
qqline(ENYPA_RmOps@data$Qs)
temp = qqnorm(ENYPA_RmOps@data$Qs, plot.it = FALSE)$x[which(ENYPA_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENYPA_RmOps@data$Qs[which(ENYPA_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(10,110))
qqnorm(SWPA_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Southwestern PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='blue', xlim = c(-4,4), ylim = c(5,110))
qqline(SWPA_RmOps@data$Qs)
temp = qqnorm(SWPA_RmOps@data$Qs, plot.it = FALSE)$x[which(SWPA_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = SWPA_RmOps@data$Qs[which(SWPA_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,110))
qqnorm(CWV_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='purple', xlim = c(-4,4), ylim = c(5,120))
qqline(CWV_RmOps@data$Qs)
temp = qqnorm(CWV_RmOps@data$Qs, plot.it = FALSE)$x[which(CWV_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CWV_RmOps@data$Qs[which(CWV_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,120))
qqnorm(MT_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='violet', xlim = c(-3.5,3.5), ylim = c(5,120))
qqline(MT_RmOps@data$Qs)
temp = qqnorm(MT_RmOps@data$Qs, plot.it = FALSE)$x[which(MT_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = MT_RmOps@data$Qs[which(MT_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(5,120))
dev.off()

#QQ-Student t3 - May not be great.
# png('QQPlotHeatFlow_BeforeCorr_NotOutTest_T3.png', res=600, units='in', width=10, height=10)
# layout(sets)
# par(mar=c(4,5,3,2))
# qqplot(qt(ppoints(CT@data$Qs), df = 3), y = CT@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Chautauqua, NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='red')
# qqline(CT@data$Qs)
# qqplot(qt(ppoints(WPA@data$Qs), df = 3), y = WPA@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='orange', xlim = c(-15,15), ylim = c(10,70))
# qqline(WPA@data$Qs)
# temp = qqplot(qt(ppoints(WPA@data$Qs), df = 3), y = WPA@data$Qs, plot.it = FALSE)$x[which(WPA@data$out_loc_error == 1)]
# par(new = TRUE)
# plot(x = temp, y = WPA@data$Qs[which(WPA@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-15,15), ylim = c(10,70))
# qqnorm(NWPANY@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Northwestern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='yellow', xlim = c(-3,3), ylim = c(5,75))
# qqline(NWPANY@data$Qs)
# temp = qqnorm(NWPANY@data$Qs, plot.it = FALSE)$x[which(NWPANY@data$out_loc_error == 1)]
# par(new = TRUE)
# plot(x = temp, y = NWPANY@data$Qs[which(NWPANY@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(5,75))
# qqnorm(CNY@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='green', xlim = c(-2.75,2.75), ylim = c(35,65))
# qqline(CNY@data$Qs)
# temp = qqnorm(CNY@data$Qs, plot.it = FALSE)$x[which(CNY@data$out_loc_error == 1)]
# par(new = TRUE)
# plot(x = temp, y = CNY@data$Qs[which(CNY@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-2.75,2.75), ylim = c(35,65))
# qqnorm(ENY@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='springgreen', xlim = c(-3,3), ylim = c(30,75))
# qqline(ENY@data$Qs)
# temp = qqnorm(ENY@data$Qs, plot.it = FALSE)$x[which(ENY@data$out_loc_error == 1)]
# par(new = TRUE)
# plot(x = temp, y = ENY@data$Qs[which(ENY@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(30,75))
# qqnorm(ENYPA@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='skyblue', xlim = c(-3,3), ylim = c(10,110))
# qqline(ENYPA@data$Qs)
# temp = qqnorm(ENYPA@data$Qs, plot.it = FALSE)$x[which(ENYPA@data$out_loc_error == 1)]
# par(new = TRUE)
# plot(x = temp, y = ENYPA@data$Qs[which(ENYPA@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(10,110))
# qqnorm(SWPA@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Southwestern PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='blue', xlim = c(-4,4), ylim = c(5,110))
# qqline(SWPA@data$Qs)
# temp = qqnorm(SWPA@data$Qs, plot.it = FALSE)$x[which(SWPA@data$out_loc_error == 1)]
# par(new = TRUE)
# plot(x = temp, y = SWPA@data$Qs[which(SWPA@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,110))
# qqnorm(CWV@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='purple', xlim = c(-4,4), ylim = c(5,120))
# qqline(CWV@data$Qs)
# temp = qqnorm(CWV@data$Qs, plot.it = FALSE)$x[which(CWV@data$out_loc_error == 1)]
# par(new = TRUE)
# plot(x = temp, y = CWV@data$Qs[which(CWV@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,120))
# qqnorm(MT@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='violet', xlim = c(-3.5,3.5), ylim = c(5,120))
# qqline(MT@data$Qs)
# temp = qqnorm(MT@data$Qs, plot.it = FALSE)$x[which(MT@data$out_loc_error == 1)]
# par(new = TRUE)
# plot(x = temp, y = MT@data$Qs[which(MT@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(5,120))
# dev.off()

#k = 1
sets = rbind(c(1,2,3), c(4,5,6),c(7,8,9))
png('QQPlotHeatFlow_NotTestedOuts_k1p0.png', res=600, units='in', width=10, height=10)
layout(sets)
par(mar=c(4,5,3,2))
qqnorm(CT_k1p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Chautauqua, NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='red')
qqline(CT_k1p0@data$Qs)
qqnorm(WPA_k1p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='orange', xlim = c(-3.5,3.5), ylim = c(10,70))
qqline(WPA_k1p0@data$Qs)
temp = qqnorm(WPA_k1p0@data$Qs, plot.it = FALSE)$x[which(WPA_k1p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = WPA_k1p0@data$Qs[which(WPA_k1p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(10,70))
qqnorm(NWPANY_k1p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Northwestern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='yellow', xlim = c(-3,3), ylim = c(5,75))
qqline(NWPANY_k1p0@data$Qs)
temp = qqnorm(NWPANY_k1p0@data$Qs, plot.it = FALSE)$x[which(NWPANY_k1p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = NWPANY_k1p0@data$Qs[which(NWPANY_k1p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(5,75))
qqnorm(CNY_k1p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='green', xlim = c(-2.75,2.75), ylim = c(35,65))
qqline(CNY_k1p0@data$Qs)
temp = qqnorm(CNY_k1p0@data$Qs, plot.it = FALSE)$x[which(CNY_k1p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CNY_k1p0@data$Qs[which(CNY_k1p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-2.75,2.75), ylim = c(35,65))
qqnorm(ENY_k1p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='springgreen', xlim = c(-3,3), ylim = c(30,75))
qqline(ENY_k1p0@data$Qs)
temp = qqnorm(ENY_k1p0@data$Qs, plot.it = FALSE)$x[which(ENY_k1p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENY_k1p0@data$Qs[which(ENY_k1p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(30,75))
qqnorm(ENYPA_k1p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='skyblue', xlim = c(-3,3), ylim = c(10,110))
qqline(ENYPA_k1p0@data$Qs)
temp = qqnorm(ENYPA_k1p0@data$Qs, plot.it = FALSE)$x[which(ENYPA_k1p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENYPA_k1p0@data$Qs[which(ENYPA_k1p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(10,110))
qqnorm(SWPA_k1p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Southwestern PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='blue', xlim = c(-4,4), ylim = c(5,110))
qqline(SWPA_k1p0@data$Qs)
temp = qqnorm(SWPA_k1p0@data$Qs, plot.it = FALSE)$x[which(SWPA_k1p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = SWPA_k1p0@data$Qs[which(SWPA_k1p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,110))
qqnorm(CWV_k1p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='purple', xlim = c(-4,4), ylim = c(5,120))
qqline(CWV_k1p0@data$Qs)
temp = qqnorm(CWV_k1p0@data$Qs, plot.it = FALSE)$x[which(CWV_k1p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CWV_k1p0@data$Qs[which(CWV_k1p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,120))
qqnorm(MT_k1p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='violet', xlim = c(-3.5,3.5), ylim = c(5,120))
qqline(MT_k1p0@data$Qs)
temp = qqnorm(MT_k1p0@data$Qs, plot.it = FALSE)$x[which(MT_k1p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = MT_k1p0@data$Qs[which(MT_k1p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(5,120))
dev.off()

#k = 1.5
png('QQPlotHeatFlow_NotTestedOuts_k1p5.png', res=600, units='in', width=10, height=10)
layout(sets)
par(mar=c(4,5,3,2))
qqnorm(CT_k1p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Chautauqua, NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='red')
qqline(CT_k1p5@data$Qs)
qqnorm(WPA_k1p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='orange', xlim = c(-3.5,3.5), ylim = c(10,70))
qqline(WPA_k1p5@data$Qs)
temp = qqnorm(WPA_k1p5@data$Qs, plot.it = FALSE)$x[which(WPA_k1p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = WPA_k1p5@data$Qs[which(WPA_k1p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(10,70))
qqnorm(NWPANY_k1p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Northwestern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='yellow', xlim = c(-3,3), ylim = c(5,75))
qqline(NWPANY_k1p5@data$Qs)
temp = qqnorm(NWPANY_k1p5@data$Qs, plot.it = FALSE)$x[which(NWPANY_k1p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = NWPANY_k1p5@data$Qs[which(NWPANY_k1p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(5,75))
qqnorm(CNY_k1p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='green', xlim = c(-2.75,2.75), ylim = c(35,65))
qqline(CNY_k1p5@data$Qs)
temp = qqnorm(CNY_k1p5@data$Qs, plot.it = FALSE)$x[which(CNY_k1p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CNY_k1p5@data$Qs[which(CNY_k1p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-2.75,2.75), ylim = c(35,65))
qqnorm(ENY_k1p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='springgreen', xlim = c(-3,3), ylim = c(30,75))
qqline(ENY_k1p5@data$Qs)
temp = qqnorm(ENY_k1p5@data$Qs, plot.it = FALSE)$x[which(ENY_k1p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENY_k1p5@data$Qs[which(ENY_k1p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(30,75))
qqnorm(ENYPA_k1p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='skyblue', xlim = c(-3,3), ylim = c(10,110))
qqline(ENYPA_k1p5@data$Qs)
temp = qqnorm(ENYPA_k1p5@data$Qs, plot.it = FALSE)$x[which(ENYPA_k1p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENYPA_k1p5@data$Qs[which(ENYPA_k1p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(10,110))
qqnorm(SWPA_k1p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Southwestern PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='blue', xlim = c(-4,4), ylim = c(5,110))
qqline(SWPA_k1p5@data$Qs)
temp = qqnorm(SWPA_k1p5@data$Qs, plot.it = FALSE)$x[which(SWPA_k1p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = SWPA_k1p5@data$Qs[which(SWPA_k1p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,110))
qqnorm(CWV_k1p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='purple', xlim = c(-4,4), ylim = c(5,120))
qqline(CWV_k1p5@data$Qs)
temp = qqnorm(CWV_k1p5@data$Qs, plot.it = FALSE)$x[which(CWV_k1p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CWV_k1p5@data$Qs[which(CWV_k1p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,120))
qqnorm(MT_k1p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='violet', xlim = c(-3.5,3.5), ylim = c(5,120))
qqline(MT_k1p5@data$Qs)
temp = qqnorm(MT_k1p5@data$Qs, plot.it = FALSE)$x[which(MT_k1p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = MT_k1p5@data$Qs[which(MT_k1p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(5,120))
dev.off()

png('QQPlotHeatFlow_NotTestedOuts_k2p0.png', res=600, units='in', width=10, height=10)
layout(sets)
par(mar=c(4,5,3,2))
qqnorm(CT_k2p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Chautauqua, NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='red')
qqline(CT_k2p0@data$Qs)
qqnorm(WPA_k2p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='orange', xlim = c(-3.5,3.5), ylim = c(10,70))
qqline(WPA_k2p0@data$Qs)
temp = qqnorm(WPA_k2p0@data$Qs, plot.it = FALSE)$x[which(WPA_k2p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = WPA_k2p0@data$Qs[which(WPA_k2p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(10,70))
qqnorm(NWPANY_k2p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Northwestern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='yellow', xlim = c(-3,3), ylim = c(5,75))
qqline(NWPANY_k2p0@data$Qs)
temp = qqnorm(NWPANY_k2p0@data$Qs, plot.it = FALSE)$x[which(NWPANY_k2p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = NWPANY_k2p0@data$Qs[which(NWPANY_k2p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(5,75))
qqnorm(CNY_k2p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='green', xlim = c(-2.75,2.75), ylim = c(35,65))
qqline(CNY_k2p0@data$Qs)
temp = qqnorm(CNY_k2p0@data$Qs, plot.it = FALSE)$x[which(CNY_k2p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CNY_k2p0@data$Qs[which(CNY_k2p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-2.75,2.75), ylim = c(35,65))
qqnorm(ENY_k2p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='springgreen', xlim = c(-3,3), ylim = c(30,75))
qqline(ENY_k2p0@data$Qs)
temp = qqnorm(ENY_k2p0@data$Qs, plot.it = FALSE)$x[which(ENY_k2p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENY_k2p0@data$Qs[which(ENY_k2p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(30,75))
qqnorm(ENYPA_k2p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='skyblue', xlim = c(-3,3), ylim = c(10,110))
qqline(ENYPA_k2p0@data$Qs)
temp = qqnorm(ENYPA_k2p0@data$Qs, plot.it = FALSE)$x[which(ENYPA_k2p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENYPA_k2p0@data$Qs[which(ENYPA_k2p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(10,110))
qqnorm(SWPA_k2p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Southwestern PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='blue', xlim = c(-4,4), ylim = c(5,110))
qqline(SWPA_k2p0@data$Qs)
temp = qqnorm(SWPA_k2p0@data$Qs, plot.it = FALSE)$x[which(SWPA_k2p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = SWPA_k2p0@data$Qs[which(SWPA_k2p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,110))
qqnorm(CWV_k2p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='purple', xlim = c(-4,4), ylim = c(5,120))
qqline(CWV_k2p0@data$Qs)
temp = qqnorm(CWV_k2p0@data$Qs, plot.it = FALSE)$x[which(CWV_k2p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CWV_k2p0@data$Qs[which(CWV_k2p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,120))
qqnorm(MT_k2p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='violet', xlim = c(-3.5,3.5), ylim = c(5,120))
qqline(MT_k2p0@data$Qs)
temp = qqnorm(MT_k2p0@data$Qs, plot.it = FALSE)$x[which(MT_k2p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = MT_k2p0@data$Qs[which(MT_k2p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(5,120))
dev.off()

png('QQPlotHeatFlow_NotTestedOuts_k2p5.png', res=600, units='in', width=10, height=10)
layout(sets)
par(mar=c(4,5,3,2))
qqnorm(CT_k2p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Chautauqua, NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='red')
qqline(CT_k2p5@data$Qs)
qqnorm(WPA_k2p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='orange', xlim = c(-3.5,3.5), ylim = c(10,70))
qqline(WPA_k2p5@data$Qs)
temp = qqnorm(WPA_k2p5@data$Qs, plot.it = FALSE)$x[which(WPA_k2p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = WPA_k2p5@data$Qs[which(WPA_k2p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(10,70))
qqnorm(NWPANY_k2p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Northwestern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='yellow', xlim = c(-3,3), ylim = c(5,75))
qqline(NWPANY_k2p5@data$Qs)
temp = qqnorm(NWPANY_k2p5@data$Qs, plot.it = FALSE)$x[which(NWPANY_k2p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = NWPANY_k2p5@data$Qs[which(NWPANY_k2p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(5,75))
qqnorm(CNY_k2p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='green', xlim = c(-2.75,2.75), ylim = c(35,65))
qqline(CNY_k2p5@data$Qs)
temp = qqnorm(CNY_k2p5@data$Qs, plot.it = FALSE)$x[which(CNY_k2p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CNY_k2p5@data$Qs[which(CNY_k2p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-2.75,2.75), ylim = c(35,65))
qqnorm(ENY_k2p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='springgreen', xlim = c(-3,3), ylim = c(30,75))
qqline(ENY_k2p5@data$Qs)
temp = qqnorm(ENY_k2p5@data$Qs, plot.it = FALSE)$x[which(ENY_k2p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENY_k2p5@data$Qs[which(ENY_k2p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(30,75))
qqnorm(ENYPA_k2p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='skyblue', xlim = c(-3,3), ylim = c(10,110))
qqline(ENYPA_k2p5@data$Qs)
temp = qqnorm(ENYPA_k2p5@data$Qs, plot.it = FALSE)$x[which(ENYPA_k2p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENYPA_k2p5@data$Qs[which(ENYPA_k2p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(10,110))
qqnorm(SWPA_k2p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Southwestern PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='blue', xlim = c(-4,4), ylim = c(5,110))
qqline(SWPA_k2p5@data$Qs)
temp = qqnorm(SWPA_k2p5@data$Qs, plot.it = FALSE)$x[which(SWPA_k2p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = SWPA_k2p5@data$Qs[which(SWPA_k2p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,110))
qqnorm(CWV_k2p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='purple', xlim = c(-4,4), ylim = c(5,120))
qqline(CWV_k2p5@data$Qs)
temp = qqnorm(CWV_k2p5@data$Qs, plot.it = FALSE)$x[which(CWV_k2p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CWV_k2p5@data$Qs[which(CWV_k2p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,120))
qqnorm(MT_k2p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='violet', xlim = c(-3.5,3.5), ylim = c(5,120))
qqline(MT_k2p5@data$Qs)
temp = qqnorm(MT_k2p5@data$Qs, plot.it = FALSE)$x[which(MT_k2p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = MT_k2p5@data$Qs[which(MT_k2p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(5,120))
dev.off()

png('QQPlotHeatFlow_NotTestedOuts_k3p5.png', res=600, units='in', width=10, height=10)
layout(sets)
par(mar=c(4,5,3,2))
qqnorm(CT_k3p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Chautauqua, NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='red')
qqline(CT_k3p5@data$Qs)
qqnorm(WPA_k3p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='orange', xlim = c(-3.5,3.5), ylim = c(10,70))
qqline(WPA_k3p5@data$Qs)
temp = qqnorm(WPA_k3p5@data$Qs, plot.it = FALSE)$x[which(WPA_k3p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = WPA_k3p5@data$Qs[which(WPA_k3p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(10,70))
qqnorm(NWPANY_k3p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Northwestern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='yellow', xlim = c(-3,3), ylim = c(5,75))
qqline(NWPANY_k3p5@data$Qs)
temp = qqnorm(NWPANY_k3p5@data$Qs, plot.it = FALSE)$x[which(NWPANY_k3p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = NWPANY_k3p5@data$Qs[which(NWPANY_k3p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(5,75))
qqnorm(CNY_k3p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='green', xlim = c(-2.75,2.75), ylim = c(35,65))
qqline(CNY_k3p5@data$Qs)
temp = qqnorm(CNY_k3p5@data$Qs, plot.it = FALSE)$x[which(CNY_k3p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CNY_k3p5@data$Qs[which(CNY_k3p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-2.75,2.75), ylim = c(35,65))
qqnorm(ENY_k3p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='springgreen', xlim = c(-3,3), ylim = c(30,75))
qqline(ENY_k3p5@data$Qs)
temp = qqnorm(ENY_k3p5@data$Qs, plot.it = FALSE)$x[which(ENY_k3p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENY_k3p5@data$Qs[which(ENY_k3p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(30,75))
qqnorm(ENYPA_k3p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='skyblue', xlim = c(-3,3), ylim = c(10,110))
qqline(ENYPA_k3p5@data$Qs)
temp = qqnorm(ENYPA_k3p5@data$Qs, plot.it = FALSE)$x[which(ENYPA_k3p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENYPA_k3p5@data$Qs[which(ENYPA_k3p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(10,110))
qqnorm(SWPA_k3p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Southwestern PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='blue', xlim = c(-4,4), ylim = c(5,110))
qqline(SWPA_k3p5@data$Qs)
temp = qqnorm(SWPA_k3p5@data$Qs, plot.it = FALSE)$x[which(SWPA_k3p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = SWPA_k3p5@data$Qs[which(SWPA_k3p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,110))
qqnorm(CWV_k3p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='purple', xlim = c(-4,4), ylim = c(5,120))
qqline(CWV_k3p5@data$Qs)
temp = qqnorm(CWV_k3p5@data$Qs, plot.it = FALSE)$x[which(CWV_k3p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CWV_k3p5@data$Qs[which(CWV_k3p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,120))
qqnorm(MT_k3p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='violet', xlim = c(-3.5,3.5), ylim = c(5,120))
qqline(MT_k3p5@data$Qs)
temp = qqnorm(MT_k3p5@data$Qs, plot.it = FALSE)$x[which(MT_k3p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = MT_k3p5@data$Qs[which(MT_k3p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(5,120))
dev.off()

png('QQPlotHeatFlow_NotTestedOuts_k4p0.png', res=600, units='in', width=10, height=10)
layout(sets)
par(mar=c(4,5,3,2))
qqnorm(CT_k4p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Chautauqua, NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='red')
qqline(CT_k4p0@data$Qs)
qqnorm(WPA_k4p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='orange', xlim = c(-3.5,3.5), ylim = c(10,70))
qqline(WPA_k4p0@data$Qs)
temp = qqnorm(WPA_k4p0@data$Qs, plot.it = FALSE)$x[which(WPA_k4p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = WPA_k4p0@data$Qs[which(WPA_k4p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(10,70))
qqnorm(NWPANY_k4p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Northwestern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='yellow', xlim = c(-3,3), ylim = c(5,75))
qqline(NWPANY_k4p0@data$Qs)
temp = qqnorm(NWPANY_k4p0@data$Qs, plot.it = FALSE)$x[which(NWPANY_k4p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = NWPANY_k4p0@data$Qs[which(NWPANY_k4p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(5,75))
qqnorm(CNY_k4p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='green', xlim = c(-2.75,2.75), ylim = c(35,65))
qqline(CNY_k4p0@data$Qs)
temp = qqnorm(CNY_k4p0@data$Qs, plot.it = FALSE)$x[which(CNY_k4p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CNY_k4p0@data$Qs[which(CNY_k4p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-2.75,2.75), ylim = c(35,65))
qqnorm(ENY_k4p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='springgreen', xlim = c(-3,3), ylim = c(30,75))
qqline(ENY_k4p0@data$Qs)
temp = qqnorm(ENY_k4p0@data$Qs, plot.it = FALSE)$x[which(ENY_k4p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENY_k4p0@data$Qs[which(ENY_k4p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(30,75))
qqnorm(ENYPA_k4p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='skyblue', xlim = c(-3,3), ylim = c(10,110))
qqline(ENYPA_k4p0@data$Qs)
temp = qqnorm(ENYPA_k4p0@data$Qs, plot.it = FALSE)$x[which(ENYPA_k4p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENYPA_k4p0@data$Qs[which(ENYPA_k4p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(10,110))
qqnorm(SWPA_k4p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Southwestern PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='blue', xlim = c(-4,4), ylim = c(5,110))
qqline(SWPA_k4p0@data$Qs)
temp = qqnorm(SWPA_k4p0@data$Qs, plot.it = FALSE)$x[which(SWPA_k4p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = SWPA_k4p0@data$Qs[which(SWPA_k4p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,110))
qqnorm(CWV_k4p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='purple', xlim = c(-4,4), ylim = c(5,120))
qqline(CWV_k4p0@data$Qs)
temp = qqnorm(CWV_k4p0@data$Qs, plot.it = FALSE)$x[which(CWV_k4p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CWV_k4p0@data$Qs[which(CWV_k4p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,120))
qqnorm(MT_k4p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='violet', xlim = c(-3.5,3.5), ylim = c(5,120))
qqline(MT_k4p0@data$Qs)
temp = qqnorm(MT_k4p0@data$Qs, plot.it = FALSE)$x[which(MT_k4p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = MT_k4p0@data$Qs[which(MT_k4p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(5,120))
dev.off()

png('QQPlotHeatFlow_NotTestedOuts_k4p5.png', res=600, units='in', width=10, height=10)
layout(sets)
par(mar=c(4,5,3,2))
qqnorm(CT_k4p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Chautauqua, NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='red')
qqline(CT_k4p5@data$Qs)
qqnorm(WPA_k4p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='orange', xlim = c(-3.5,3.5), ylim = c(10,70))
qqline(WPA_k4p5@data$Qs)
temp = qqnorm(WPA_k4p5@data$Qs, plot.it = FALSE)$x[which(WPA_k4p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = WPA_k4p5@data$Qs[which(WPA_k4p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(10,70))
qqnorm(NWPANY_k4p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Northwestern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='yellow', xlim = c(-3,3), ylim = c(5,75))
qqline(NWPANY_k4p5@data$Qs)
temp = qqnorm(NWPANY_k4p5@data$Qs, plot.it = FALSE)$x[which(NWPANY_k4p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = NWPANY_k4p5@data$Qs[which(NWPANY_k4p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(5,75))
qqnorm(CNY_k4p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='green', xlim = c(-2.75,2.75), ylim = c(35,65))
qqline(CNY_k4p5@data$Qs)
temp = qqnorm(CNY_k4p5@data$Qs, plot.it = FALSE)$x[which(CNY_k4p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CNY_k4p5@data$Qs[which(CNY_k4p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-2.75,2.75), ylim = c(35,65))
qqnorm(ENY_k4p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='springgreen', xlim = c(-3,3), ylim = c(30,75))
qqline(ENY_k4p5@data$Qs)
temp = qqnorm(ENY_k4p5@data$Qs, plot.it = FALSE)$x[which(ENY_k4p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENY_k4p5@data$Qs[which(ENY_k4p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(30,75))
qqnorm(ENYPA_k4p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='skyblue', xlim = c(-3,3), ylim = c(10,110))
qqline(ENYPA_k4p5@data$Qs)
temp = qqnorm(ENYPA_k4p5@data$Qs, plot.it = FALSE)$x[which(ENYPA_k4p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENYPA_k4p5@data$Qs[which(ENYPA_k4p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(10,110))
qqnorm(SWPA_k4p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Southwestern PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='blue', xlim = c(-4,4), ylim = c(5,110))
qqline(SWPA_k4p5@data$Qs)
temp = qqnorm(SWPA_k4p5@data$Qs, plot.it = FALSE)$x[which(SWPA_k4p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = SWPA_k4p5@data$Qs[which(SWPA_k4p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,110))
qqnorm(CWV_k4p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='purple', xlim = c(-4,4), ylim = c(5,120))
qqline(CWV_k4p5@data$Qs)
temp = qqnorm(CWV_k4p5@data$Qs, plot.it = FALSE)$x[which(CWV_k4p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CWV_k4p5@data$Qs[which(CWV_k4p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,120))
qqnorm(MT_k4p5@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='violet', xlim = c(-3.5,3.5), ylim = c(5,120))
qqline(MT_k4p5@data$Qs)
temp = qqnorm(MT_k4p5@data$Qs, plot.it = FALSE)$x[which(MT_k4p5@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = MT_k4p5@data$Qs[which(MT_k4p5@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(5,120))
dev.off()

png('QQPlotHeatFlow_NotTestedOuts_k5p0.png', res=600, units='in', width=10, height=10)
layout(sets)
par(mar=c(4,5,3,2))
qqnorm(CT_k5p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Chautauqua, NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='red')
qqline(CT_k5p0@data$Qs)
qqnorm(WPA_k5p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='orange', xlim = c(-3.5,3.5), ylim = c(10,70))
qqline(WPA_k5p0@data$Qs)
temp = qqnorm(WPA_k5p0@data$Qs, plot.it = FALSE)$x[which(WPA_k5p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = WPA_k5p0@data$Qs[which(WPA_k5p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(10,70))
qqnorm(NWPANY_k5p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Northwestern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='yellow', xlim = c(-3,3), ylim = c(5,75))
qqline(NWPANY_k5p0@data$Qs)
temp = qqnorm(NWPANY_k5p0@data$Qs, plot.it = FALSE)$x[which(NWPANY_k5p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = NWPANY_k5p0@data$Qs[which(NWPANY_k5p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(5,75))
qqnorm(CNY_k5p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='green', xlim = c(-2.75,2.75), ylim = c(35,65))
qqline(CNY_k5p0@data$Qs)
temp = qqnorm(CNY_k5p0@data$Qs, plot.it = FALSE)$x[which(CNY_k5p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CNY_k5p0@data$Qs[which(CNY_k5p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-2.75,2.75), ylim = c(35,65))
qqnorm(ENY_k5p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='springgreen', xlim = c(-3,3), ylim = c(30,75))
qqline(ENY_k5p0@data$Qs)
temp = qqnorm(ENY_k5p0@data$Qs, plot.it = FALSE)$x[which(ENY_k5p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENY_k5p0@data$Qs[which(ENY_k5p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(30,75))
qqnorm(ENYPA_k5p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='skyblue', xlim = c(-3,3), ylim = c(10,110))
qqline(ENYPA_k5p0@data$Qs)
temp = qqnorm(ENYPA_k5p0@data$Qs, plot.it = FALSE)$x[which(ENYPA_k5p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENYPA_k5p0@data$Qs[which(ENYPA_k5p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(10,110))
qqnorm(SWPA_k5p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Southwestern PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='blue', xlim = c(-4,4), ylim = c(5,110))
qqline(SWPA_k5p0@data$Qs)
temp = qqnorm(SWPA_k5p0@data$Qs, plot.it = FALSE)$x[which(SWPA_k5p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = SWPA_k5p0@data$Qs[which(SWPA_k5p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,110))
qqnorm(CWV_k5p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='purple', xlim = c(-4,4), ylim = c(5,120))
qqline(CWV_k5p0@data$Qs)
temp = qqnorm(CWV_k5p0@data$Qs, plot.it = FALSE)$x[which(CWV_k5p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CWV_k5p0@data$Qs[which(CWV_k5p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,120))
qqnorm(MT_k5p0@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='violet', xlim = c(-3.5,3.5), ylim = c(5,120))
qqline(MT_k5p0@data$Qs)
temp = qqnorm(MT_k5p0@data$Qs, plot.it = FALSE)$x[which(MT_k5p0@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = MT_k5p0@data$Qs[which(MT_k5p0@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(5,120))
dev.off()

#With Map for k = 3
sets = rbind(c(1,2,3,4), c(5,6,7,8),c(9,10,11,11))
png('QQPlotHeatFlow_NotTestedOuts_Map.png', res=600, units='in', width=13, height=10)
layout(sets)
par(mar=c(4,5,3,2))
qqnorm(CT@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Chautauqua, NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='red')
qqline(CT@data$Qs)
qqnorm(WPA@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='orange', xlim = c(-3.5,3.5), ylim = c(10,70))
qqline(WPA@data$Qs)
temp = qqnorm(WPA@data$Qs, plot.it = FALSE)$x[which(WPA@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = WPA@data$Qs[which(WPA@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(10,70))
qqnorm(NWPANY@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Northwestern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='yellow', xlim = c(-3,3), ylim = c(5,75))
qqline(NWPANY@data$Qs)
temp = qqnorm(NWPANY@data$Qs, plot.it = FALSE)$x[which(NWPANY@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = NWPANY@data$Qs[which(NWPANY@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(5,75))
qqnorm(CNY@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='green', xlim = c(-2.75,2.75), ylim = c(35,65))
qqline(CNY@data$Qs)
temp = qqnorm(CNY@data$Qs, plot.it = FALSE)$x[which(CNY@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CNY@data$Qs[which(CNY@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-2.75,2.75), ylim = c(35,65))
qqnorm(ENY@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='springgreen', xlim = c(-3,3), ylim = c(30,75))
qqline(ENY@data$Qs)
temp = qqnorm(ENY@data$Qs, plot.it = FALSE)$x[which(ENY@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENY@data$Qs[which(ENY@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(30,75))
qqnorm(ENYPA@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='skyblue', xlim = c(-3,3), ylim = c(10,110))
qqline(ENYPA@data$Qs)
temp = qqnorm(ENYPA@data$Qs, plot.it = FALSE)$x[which(ENYPA@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENYPA@data$Qs[which(ENYPA@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(10,110))
qqnorm(SWPA@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Southwestern PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='blue', xlim = c(-4,4), ylim = c(5,110))
qqline(SWPA@data$Qs)
temp = qqnorm(SWPA@data$Qs, plot.it = FALSE)$x[which(SWPA@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = SWPA@data$Qs[which(SWPA@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,110))
qqnorm(CWV@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='purple', xlim = c(-4,4), ylim = c(5,120))
qqline(CWV@data$Qs)
temp = qqnorm(CWV@data$Qs, plot.it = FALSE)$x[which(CWV@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CWV@data$Qs[which(CWV@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,120))
qqnorm(MT@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='violet', xlim = c(-3.5,3.5), ylim = c(5,120))
qqline(MT@data$Qs)
temp = qqnorm(MT@data$Qs, plot.it = FALSE)$x[which(MT@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = MT@data$Qs[which(MT@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(5,120))
qqnorm(VR@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Valley and Ridge', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='grey', xlim = c(-3.5,3.5), ylim = c(5,120))
qqline(VR@data$Qs)
temp = qqnorm(VR@data$Qs, plot.it = FALSE)$x[which(VR@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = VR@data$Qs[which(VR@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(5,120))

#Map
par(xaxs = 'i', yaxs = 'i', mar = c(2,10,2,10))
#par(pin = c(par('pin')[1], ratio*par('pin')[1]))
plot(InterpRegs, xlim = c(-82.64474, -74.5), ylim = c(36.75, 43.4))
plot(InterpRegs[which(InterpRegs$Name == "CT"),], col = 'red', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "WPA"),], col = 'orange', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "NWPANY"),], col = 'yellow', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "CNY"),], col = 'green', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "ENY"),], col = 'springgreen', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "ENYPA"),], col = 'skyblue', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "SWPA"),], col = 'blue', add = TRUE)
plot(CWV_Bounded, col = 'purple', add = TRUE)
plot(MT_Bounded, col = 'magenta', add = TRUE)
plot(VR_Bounded, col = 'gray', add = TRUE)
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
box()
dev.off()

png('QQPlotHeatFlow_NotTestedOuts_Map_RmOps.png', res=600, units='in', width=13, height=10)
layout(sets)
par(mar=c(4,5,3,2))
qqnorm(CT_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Chautauqua, NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='red')
qqline(CT_RmOps@data$Qs)
qqnorm(WPA_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='orange', xlim = c(-3.5,3.5), ylim = c(10,70))
qqline(WPA_RmOps@data$Qs)
temp = qqnorm(WPA_RmOps@data$Qs, plot.it = FALSE)$x[which(WPA_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = WPA_RmOps@data$Qs[which(WPA_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(10,70))
qqnorm(NWPANY_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Northwestern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='yellow', xlim = c(-3,3), ylim = c(5,75))
qqline(NWPANY_RmOps@data$Qs)
temp = qqnorm(NWPANY_RmOps@data$Qs, plot.it = FALSE)$x[which(NWPANY_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = NWPANY_RmOps@data$Qs[which(NWPANY_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(5,75))
qqnorm(CNY_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='green', xlim = c(-2.75,2.75), ylim = c(35,65))
qqline(CNY_RmOps@data$Qs)
temp = qqnorm(CNY_RmOps@data$Qs, plot.it = FALSE)$x[which(CNY_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CNY_RmOps@data$Qs[which(CNY_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-2.75,2.75), ylim = c(35,65))
qqnorm(ENY_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='springgreen', xlim = c(-3,3), ylim = c(30,75))
qqline(ENY_RmOps@data$Qs)
temp = qqnorm(ENY_RmOps@data$Qs, plot.it = FALSE)$x[which(ENY_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENY_RmOps@data$Qs[which(ENY_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(30,75))
qqnorm(ENYPA_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='skyblue', xlim = c(-3,3), ylim = c(10,110))
qqline(ENYPA_RmOps@data$Qs)
temp = qqnorm(ENYPA_RmOps@data$Qs, plot.it = FALSE)$x[which(ENYPA_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = ENYPA_RmOps@data$Qs[which(ENYPA_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(10,110))
qqnorm(SWPA_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Southwestern PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='blue', xlim = c(-4,4), ylim = c(5,110))
qqline(SWPA_RmOps@data$Qs)
temp = qqnorm(SWPA_RmOps@data$Qs, plot.it = FALSE)$x[which(SWPA_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = SWPA_RmOps@data$Qs[which(SWPA_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,110))
qqnorm(CWV_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='purple', xlim = c(-4,4), ylim = c(5,120))
qqline(CWV_RmOps@data$Qs)
temp = qqnorm(CWV_RmOps@data$Qs, plot.it = FALSE)$x[which(CWV_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = CWV_RmOps@data$Qs[which(CWV_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,120))
qqnorm(MT_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='violet', xlim = c(-3.5,3.5), ylim = c(5,120))
qqline(MT_RmOps@data$Qs)
temp = qqnorm(MT_RmOps@data$Qs, plot.it = FALSE)$x[which(MT_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = MT_RmOps@data$Qs[which(MT_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(5,120))
qqnorm(VR_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Valley and Ridge', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='grey', xlim = c(-3.5,3.5), ylim = c(5,120))
qqline(VR_RmOps@data$Qs)
temp = qqnorm(VR_RmOps@data$Qs, plot.it = FALSE)$x[which(VR_RmOps@data$out_loc_error == 1)]
par(new = TRUE)
plot(x = temp, y = VR_RmOps@data$Qs[which(VR_RmOps@data$out_loc_error == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(5,120))

#Map
par(xaxs = 'i', yaxs = 'i', mar = c(2,10,2,10))
#par(pin = c(par('pin')[1], ratio*par('pin')[1]))
plot(InterpRegs, xlim = c(-82.64474, -74.5), ylim = c(36.75, 43.4))
plot(InterpRegs[which(InterpRegs$Name == "CT"),], col = 'red', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "WPA"),], col = 'orange', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "NWPANY"),], col = 'yellow', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "CNY"),], col = 'green', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "ENY"),], col = 'springgreen', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "ENYPA"),], col = 'skyblue', add = TRUE)
plot(InterpRegs[which(InterpRegs$Name == "SWPA"),], col = 'blue', add = TRUE)
plot(CWV_Bounded, col = 'purple', add = TRUE)
plot(MT_Bounded, col = 'magenta', add = TRUE)
plot(VR_Bounded, col = 'gray', add = TRUE)
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
box()
dev.off()

rm(temp)

#    Label points as outliers and switch geologic region ownership for one point----
#Heatflow
qqnorm(CT@data$Qs)
qqline(CT@data$Qs)
qqnorm(CNY@data$Qs)
qqline(CNY@data$Qs)
qqnorm(CWV@data$Qs) #The one point above 100 mW/m^2 should probably be in MT, but it is not a spatial outlier. Place in MT and remove from CWV
qqline(CWV@data$Qs)
MT = rbind(MT, CWV[which(CWV$Qs > 110),])
CWV = CWV[which(CWV$Qs < 110),]
qqnorm(CWV@data$Qs)
qqline(CWV@data$Qs)
qqnorm(ENY@data$Qs)#Based on points that were not tested as outliers, RowID_12690 appears too high on the map based on surrounding points. It is also in the upper tail.
qqline(ENY@data$Qs)
ENY = ENY[-which(ENY$RowID_ == 12690),]
FL = FL[-which(FL$RowID_ == 12690),]
ENY_RmOps = ENY_RmOps[-which(ENY_RmOps$RowID_ == 12690),]
FL_RmOps = FL_RmOps[-which(FL_RmOps$RowID_ == 12690),]
qqnorm(ENYPA@data$Qs) #Definite outlier that did not get removed is in here. 30 mW/m^3 greater than others. Remove it from ENYPA and FL.
FL = FL[-which(FL$RowID_ == 19770),]
ENYPA = ENYPA[which(ENYPA$Qs < 100),]
FL_RmOps = FL_RmOps[-which(FL_RmOps$RowID_ == 19770),]
ENYPA_RmOps = ENYPA_RmOps[which(ENYPA_RmOps$Qs < 100),]
qqnorm(ENYPA@data$Qs)
qqline(ENYPA@data$Qs)
qqnorm(MT@data$Qs) #The two < 25 mW/m^2 points are clustered in an area with low heat flow. Looks OK
qqline(MT@data$Qs)
qqnorm(NWPANY@data$Qs) #The low heat flow of 7 was not tested for outliers. It should be removed.
FL = FL[-which(FL$RowID_ == 29908),]
NWPANY = NWPANY[which(NWPANY$Qs > 10),]
FL_RmOps = FL_RmOps[-which(FL_RmOps$RowID_ == 29908),]
NWPANY_RmOps = NWPANY_RmOps[which(NWPANY_RmOps$Qs > 10),]
qqnorm(NWPANY@data$Qs)
qqline(NWPANY@data$Qs)
qqnorm(SWPA@data$Qs)
qqline(SWPA@data$Qs)
qqnorm(WPA@data$Qs)
qqline(WPA@data$Qs)
qqnorm(VR@data$Qs)
qqline(VR@data$Qs)
qqnorm(FL@data$Qs) #Only the ENYPA, ENY, and NWPANY points are removed.
qqline(FL@data$Qs)

#Check that the number of wells in FL equals the total in all other regions
if (nrow(FL) - nrow(CT) - nrow(CNY) - nrow(CWV) - nrow(ENY) - nrow(ENYPA) - nrow(SWPA) - nrow(NWPANY) - nrow(WPA) - nrow(MT) - nrow(VR) != 0){
  print('Number of wells in FL is different than the total in other regions.')
}

FL_k1p0 = FL_k1p0[-which(FL_k1p0$RowID_ == 19770),]
ENYPA_k1p0 = ENYPA_k1p0[-which(ENYPA_k1p0$RowID_ == 19770),]
FL_k1p0 = FL_k1p0[-which(FL_k1p0$RowID_ == 29908),]
NWPANY_k1p0 = NWPANY_k1p0[-which(NWPANY_k1p0$RowID_ == 29908),]
ENY_k1p0 = ENY_k1p0[-which(ENY_k1p0$RowID_ == 12690),]
FL_k1p0 = FL_k1p0[-which(FL_k1p0$RowID_ == 12690),]

FL_k1p5 = FL_k1p5[-which(FL_k1p5$RowID_ == 19770),]
ENYPA_k1p5 = ENYPA_k1p5[-which(ENYPA_k1p5$RowID_ == 19770),]
FL_k1p5 = FL_k1p5[-which(FL_k1p5$RowID_ == 29908),]
NWPANY_k1p5 = NWPANY_k1p5[-which(NWPANY_k1p5$RowID_ == 29908),]
ENY_k1p5 = ENY_k1p5[-which(ENY_k1p5$RowID_ == 12690),]
FL_k1p5 = FL_k1p5[-which(FL_k1p5$RowID_ == 12690),]

FL_k2p0 = FL_k2p0[-which(FL_k2p0$RowID_ == 19770),]
ENYPA_k2p0 = ENYPA_k2p0[-which(ENYPA_k2p0$RowID_ == 19770),]
FL_k2p0 = FL_k2p0[-which(FL_k2p0$RowID_ == 29908),]
NWPANY_k2p0 = NWPANY_k2p0[-which(NWPANY_k2p0$RowID_ == 29908),]
ENY_k2p0 = ENY_k2p0[-which(ENY_k2p0$RowID_ == 12690),]
FL_k2p0 = FL_k2p0[-which(FL_k2p0$RowID_ == 12690),]

FL_k2p5 = FL_k2p5[-which(FL_k2p5$RowID_ == 19770),]
ENYPA_k2p5 = ENYPA_k2p5[-which(ENYPA_k2p5$RowID_ == 19770),]
FL_k2p5 = FL_k2p5[-which(FL_k2p5$RowID_ == 29908),]
NWPANY_k2p5 = NWPANY_k2p5[-which(NWPANY_k2p5$RowID_ == 29908),]
ENY_k2p5 = ENY_k2p5[-which(ENY_k2p5$RowID_ == 12690),]
FL_k2p5 = FL_k2p5[-which(FL_k2p5$RowID_ == 12690),]

FL_k3p5 = FL_k3p5[-which(FL_k3p5$RowID_ == 19770),]
ENYPA_k3p5 = ENYPA_k3p5[-which(ENYPA_k3p5$RowID_ == 19770),]
FL_k3p5 = FL_k3p5[-which(FL_k3p5$RowID_ == 29908),]
NWPANY_k3p5 = NWPANY_k3p5[-which(NWPANY_k3p5$RowID_ == 29908),]
ENY_k3p5 = ENY_k3p5[-which(ENY_k3p5$RowID_ == 12690),]
FL_k3p5 = FL_k3p5[-which(FL_k3p5$RowID_ == 12690),]

FL_k4p0 = FL_k4p0[-which(FL_k4p0$RowID_ == 19770),]
ENYPA_k4p0 = ENYPA_k4p0[-which(ENYPA_k4p0$RowID_ == 19770),]
FL_k4p0 = FL_k4p0[-which(FL_k4p0$RowID_ == 29908),]
NWPANY_k4p0 = NWPANY_k4p0[-which(NWPANY_k4p0$RowID_ == 29908),]
ENY_k4p0 = ENY_k4p0[-which(ENY_k4p0$RowID_ == 12690),]
FL_k4p0 = FL_k4p0[-which(FL_k4p0$RowID_ == 12690),]

FL_k4p5 = FL_k4p5[-which(FL_k4p5$RowID_ == 19770),]
ENYPA_k4p5 = ENYPA_k4p5[-which(ENYPA_k4p5$RowID_ == 19770),]
FL_k4p5 = FL_k4p5[-which(FL_k4p5$RowID_ == 29908),]
NWPANY_k4p5 = NWPANY_k4p5[-which(NWPANY_k4p5$RowID_ == 29908),]
ENY_k4p5 = ENY_k4p5[-which(ENY_k4p5$RowID_ == 12690),]
FL_k4p5 = FL_k4p5[-which(FL_k4p5$RowID_ == 12690),]

FL_k5p0 = FL_k5p0[-which(FL_k5p0$RowID_ == 19770),]
ENYPA_k5p0 = ENYPA_k5p0[-which(ENYPA_k5p0$RowID_ == 19770),]
FL_k5p0 = FL_k5p0[-which(FL_k5p0$RowID_ == 29908),]
NWPANY_k5p0 = NWPANY_k5p0[-which(NWPANY_k5p0$RowID_ == 29908),]
ENY_k5p0 = ENY_k5p0[-which(ENY_k5p0$RowID_ == 12690),]
FL_k5p0 = FL_k5p0[-which(FL_k5p0$RowID_ == 12690),]

#QQ Plot Figure - corrected point locations
sets = rbind(c(1,2,3), c(4,5,6),c(7,8,9))
png('QQPlotHeatFlow_corrPoints_2018.png', res=300, units='in', width=10, height=10)
layout(sets)
par(mar=c(4,5,3,2))
qqnorm(CT@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Chautauqua, NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='red')
qqline(CT@data$Qs)
qqnorm(WPA@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='orange')
qqline(WPA@data$Qs)
qqnorm(NWPANY@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Northwestern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='yellow')
qqline(NWPANY@data$Qs)
qqnorm(CNY@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='green')
qqline(CNY@data$Qs)
qqnorm(ENY@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='springgreen')
qqline(ENY@data$Qs)
qqnorm(ENYPA@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='skyblue')
qqline(ENYPA@data$Qs)
qqnorm(SWPA@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Southwestern PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='blue')
qqline(SWPA@data$Qs)
qqnorm(CWV@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='purple')
qqline(CWV@data$Qs)
qqnorm(MT@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='violet')
qqline(MT@data$Qs)
dev.off()

png('QQPlotHeatFlow_corrPoints_RmOps_2018.png', res=300, units='in', width=10, height=10)
layout(sets)
par(mar=c(4,5,3,2))
qqnorm(CT_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Chautauqua, NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='red')
qqline(CT_RmOps@data$Qs)
qqnorm(WPA_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='orange')
qqline(WPA_RmOps@data$Qs)
qqnorm(NWPANY_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Northwestern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='yellow')
qqline(NWPANY_RmOps@data$Qs)
qqnorm(CNY_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='green')
qqline(CNY_RmOps@data$Qs)
qqnorm(ENY_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='springgreen')
qqline(ENY_RmOps@data$Qs)
qqnorm(ENYPA_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='skyblue')
qqline(ENYPA_RmOps@data$Qs)
qqnorm(SWPA_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Southwestern PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='blue')
qqline(SWPA_RmOps@data$Qs)
qqnorm(CWV_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='purple')
qqline(CWV_RmOps@data$Qs)
qqnorm(MT_RmOps@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='violet')
qqline(MT_RmOps@data$Qs)
dev.off()

qqnorm(FL@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Full Region', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, ylim=c(0,120), xlim=c(-4,4))
qqline(FL@data$Qs)

#  Fixme: Post analysis of spatial correlation of outlier depth ranks----

#  Post analysis of local median deviation----
FL_NAD = spTransform(FL, CRSobj = CRS('+init=epsg:26917'))
DataTab_Post = QsDev(Data = FL_NAD@data, Var = 'Qs', xName = 'coords_x1', yName = 'coords_x2', rad = 10000, max_pts = 25)
#Add information to the original data
FL_NAD@data = DataTab_Post

#Color function parameters for plotting the heat flow data map
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#With map
FL_WGS = spTransform(FL_NAD, CRS('+init=epsg:4326'))
sets = rbind(c(1,2,7,7), c(3,4,7,7), c(5,6,8,8), c(9,9,8,8))
png("HeatFlowEDA_SplitAll_Map_Log_WellDepthSort_MedDiff_DeepWells_Box_Post.png", width=10, height=10, units="in", res=600)
layout(sets)
EDAPlots = function(DataAll, Var, Unit, ymin, ymax){
  plot(DataAll$WellDepth[which(DataAll$State == "MD")], DataAll@data[Var][,1][which(DataAll$State == "MD")], col = "orange", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Maryland', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  par(tck = -0.03)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  par(tck = -0.015)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
  par(xpd = TRUE)
  text(x = -2000, y = 5000, expression(bold('A')), cex = 2)
  par(xpd = FALSE)
  plot(DataAll$WellDepth[which(DataAll$State == "KY")], DataAll@data[Var][,1][which(DataAll$State == "KY")], col = "purple", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Kentucky', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  par(tck = -0.03)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  par(tck = -0.015)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "VA")], DataAll@data[Var][,1][which(DataAll$State == "VA")], col = "yellow", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab=Unit, xlab='Depth of BHT (m)', main='Virginia', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  par(tck = -0.03)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  par(tck = -0.015)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "WV")], DataAll@data[Var][,1][which(DataAll$State == "WV")], col = "green", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='West Virginia', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  par(tck = -0.03)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  par(tck = -0.015)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "PA")], DataAll@data[Var][,1][which(DataAll$State == "PA")], col = "red", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Pennsylvania', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  par(tck = -0.03)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  par(tck = -0.015)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "NY")], DataAll@data[Var][,1][which(DataAll$State == "NY")], col = "blue", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='New York', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  par(tck = -0.03)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  par(tck = -0.015)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
}
par(mar=c(4,5.5,3,2), xaxs = 'i', yaxs = 'i')
EDAPlots(FL_WGS, "Qs", Unit=expression("Surface Heat Flow" ~ (mW/m^2)),.5, 2000)
par(tck = NA)
EDAPlotsMed = function(DataAll, Var, Unit, ymin, ymax){
  plot(DataAll$WellDepth[which(DataAll$State == "NY" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "NY" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "NY" & DataAll$RegMed > -9999), 'RegMed'], col = "blue", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab=Unit, xlab='Depth of BHT (m)', main='Deviation from Local Median Heat Flow', cex.main=2, cex.axis=1.5, cex.lab=1.5, axes = FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "WV" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "WV" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "WV" & DataAll$RegMed > -9999), 'RegMed'], col = "green", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "PA" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "PA" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "PA" & DataAll$RegMed > -9999), 'RegMed'], col = "red", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "KY" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "KY" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "KY" & DataAll$RegMed > -9999), 'RegMed'], col = "purple", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "VA" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "VA" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "VA" & DataAll$RegMed > -9999), 'RegMed'], col = "yellow", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
  par(new=T)
  plot(DataAll$WellDepth[which(DataAll$State == "MD" & DataAll$RegMed > -9999)], DataAll@data[which(DataAll$State == "MD" & DataAll$RegMed > -9999), Var] - DataAll@data[which(DataAll$State == "MD" & DataAll$RegMed > -9999), 'RegMed'], col = "orange", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='', main='', axes=FALSE)
}
EDAPlotsMed(FL_WGS, "Qs", Unit=expression(paste("Site Q"['s'], " - Local Median Q"['s'], " (mW/m"^2, ")")),-100, 100)
axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.5)
axis(side = 2, at = seq(-200,1300,50), labels=TRUE, cex.axis=1.5)
minor.tick(nx = 5, ny = 1)
lines(c(-100,10000), c(0,0))
legend('topright', legend=c("New York", "Pennsylvania", "West Virginia", "Kentucky", "Maryland", "Virginia"), pch=16, col=c("blue", "red", "green", "purple","orange","yellow"), cex=2)
box(which = 'figure', lwd = 2)
par(xpd = TRUE)
text(x = -800, y = 625, expression(bold('C')), cex = 2)
par(xpd = FALSE)
par(mar = c(2,2.5,1.2,1), xaxs = 'i', yaxs = 'i')
plot(FL_WGS, pch = 16, col = "white", cex = 0.2)
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], add=TRUE)
plot(FL_WGS[order(FL_WGS$WellDepth, decreasing = FALSE),], pch = 16, col = colFun(FL_WGS$Qs[order(FL_WGS$WellDepth, decreasing = FALSE)]), cex = 0.3, add = TRUE)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('topleft', title = expression(paste('Q'['s'], ' (mW/m'^2,')' )), legend = c('<40', '40 - 50', '50 - 60', '60 - 70', '70 - 80', '>80'), pch = 16, cex = 1.8, col = colFun(c(30, 45, 55, 65, 75, 90)))
box(which = 'figure', lwd = 2)
par(xpd = TRUE)
text(x = -83.5, y = 43.5, expression(bold('D')), cex = 2)
par(mar=c(4,5.5,3,2), xpd = FALSE)
hist(FL_WGS$Qs, breaks = seq(0, 1500, 5), freq = TRUE, xlim = c(0,120), ylim = c(0,4000), xlab = expression(paste('Surface Heat Flow (mW/m'^2,')')), cex.lab = 1.5, cex.axis = 1.5, main = 'Histogram of All Data', cex.main = 2)
box(which = 'figure', lwd = 2)
par(xpd = TRUE)
text(x = -15, y = 5500, expression(bold('B')), cex = 2)
dev.off()

#  Post analysis of variograms in each geologic region pre and post----
#   Transform to NAD UTM17N----
CT = spTransform(CT, CRS('+init=epsg:26917'))
CNY = spTransform(CNY, CRS('+init=epsg:26917'))
CWV = spTransform(CWV, CRS('+init=epsg:26917'))
ENY = spTransform(ENY, CRS('+init=epsg:26917'))
ENYPA = spTransform(ENYPA, CRS('+init=epsg:26917'))
MT = spTransform(MT, CRS('+init=epsg:26917'))
NWPANY = spTransform(NWPANY, CRS('+init=epsg:26917'))
SWPA = spTransform(SWPA, CRS('+init=epsg:26917'))
WPA = spTransform(WPA, CRS('+init=epsg:26917'))
VR = spTransform(VR, CRS('+init=epsg:26917'))
FL = spTransform(FL, CRS('+init=epsg:26917'))

CT_RmOps = spTransform(CT_RmOps, CRS('+init=epsg:26917'))
CNY_RmOps = spTransform(CNY_RmOps, CRS('+init=epsg:26917'))
CWV_RmOps = spTransform(CWV_RmOps, CRS('+init=epsg:26917'))
ENY_RmOps = spTransform(ENY_RmOps, CRS('+init=epsg:26917'))
ENYPA_RmOps = spTransform(ENYPA_RmOps, CRS('+init=epsg:26917'))
MT_RmOps = spTransform(MT_RmOps, CRS('+init=epsg:26917'))
NWPANY_RmOps = spTransform(NWPANY_RmOps, CRS('+init=epsg:26917'))
SWPA_RmOps = spTransform(SWPA_RmOps, CRS('+init=epsg:26917'))
WPA_RmOps = spTransform(WPA_RmOps, CRS('+init=epsg:26917'))
VR_RmOps = spTransform(VR_RmOps, CRS('+init=epsg:26917'))
FL_RmOps = spTransform(FL_RmOps, CRS('+init=epsg:26917'))

CT_k1p0 = spTransform(CT_k1p0, CRS('+init=epsg:26917'))
CNY_k1p0 = spTransform(CNY_k1p0, CRS('+init=epsg:26917'))
CWV_k1p0 = spTransform(CWV_k1p0, CRS('+init=epsg:26917'))
ENY_k1p0 = spTransform(ENY_k1p0, CRS('+init=epsg:26917'))
ENYPA_k1p0 = spTransform(ENYPA_k1p0, CRS('+init=epsg:26917'))
MT_k1p0 = spTransform(MT_k1p0, CRS('+init=epsg:26917'))
NWPANY_k1p0 = spTransform(NWPANY_k1p0, CRS('+init=epsg:26917'))
SWPA_k1p0 = spTransform(SWPA_k1p0, CRS('+init=epsg:26917'))
WPA_k1p0 = spTransform(WPA_k1p0, CRS('+init=epsg:26917'))
VR_k1p0 = spTransform(VR_k1p0, CRS('+init=epsg:26917'))
FL_k1p0 = spTransform(FL_k1p0, CRS('+init=epsg:26917'))

CT_k1p5 = spTransform(CT_k1p5, CRS('+init=epsg:26917'))
CNY_k1p5 = spTransform(CNY_k1p5, CRS('+init=epsg:26917'))
CWV_k1p5 = spTransform(CWV_k1p5, CRS('+init=epsg:26917'))
ENY_k1p5 = spTransform(ENY_k1p5, CRS('+init=epsg:26917'))
ENYPA_k1p5 = spTransform(ENYPA_k1p5, CRS('+init=epsg:26917'))
MT_k1p5 = spTransform(MT_k1p5, CRS('+init=epsg:26917'))
NWPANY_k1p5 = spTransform(NWPANY_k1p5, CRS('+init=epsg:26917'))
SWPA_k1p5 = spTransform(SWPA_k1p5, CRS('+init=epsg:26917'))
WPA_k1p5 = spTransform(WPA_k1p5, CRS('+init=epsg:26917'))
VR_k1p5 = spTransform(VR_k1p5, CRS('+init=epsg:26917'))
FL_k1p5 = spTransform(FL_k1p5, CRS('+init=epsg:26917'))

CT_k2p0 = spTransform(CT_k2p0, CRS('+init=epsg:26917'))
CNY_k2p0 = spTransform(CNY_k2p0, CRS('+init=epsg:26917'))
CWV_k2p0 = spTransform(CWV_k2p0, CRS('+init=epsg:26917'))
ENY_k2p0 = spTransform(ENY_k2p0, CRS('+init=epsg:26917'))
ENYPA_k2p0 = spTransform(ENYPA_k2p0, CRS('+init=epsg:26917'))
MT_k2p0 = spTransform(MT_k2p0, CRS('+init=epsg:26917'))
NWPANY_k2p0 = spTransform(NWPANY_k2p0, CRS('+init=epsg:26917'))
SWPA_k2p0 = spTransform(SWPA_k2p0, CRS('+init=epsg:26917'))
WPA_k2p0 = spTransform(WPA_k2p0, CRS('+init=epsg:26917'))
VR_k2p0 = spTransform(VR_k2p0, CRS('+init=epsg:26917'))
FL_k2p0 = spTransform(FL_k2p0, CRS('+init=epsg:26917'))

CT_k2p5 = spTransform(CT_k2p5, CRS('+init=epsg:26917'))
CNY_k2p5 = spTransform(CNY_k2p5, CRS('+init=epsg:26917'))
CWV_k2p5 = spTransform(CWV_k2p5, CRS('+init=epsg:26917'))
ENY_k2p5 = spTransform(ENY_k2p5, CRS('+init=epsg:26917'))
ENYPA_k2p5 = spTransform(ENYPA_k2p5, CRS('+init=epsg:26917'))
MT_k2p5 = spTransform(MT_k2p5, CRS('+init=epsg:26917'))
NWPANY_k2p5 = spTransform(NWPANY_k2p5, CRS('+init=epsg:26917'))
SWPA_k2p5 = spTransform(SWPA_k2p5, CRS('+init=epsg:26917'))
WPA_k2p5 = spTransform(WPA_k2p5, CRS('+init=epsg:26917'))
VR_k2p5 = spTransform(VR_k2p5, CRS('+init=epsg:26917'))
FL_k2p5 = spTransform(FL_k2p5, CRS('+init=epsg:26917'))

CT_k3p5 = spTransform(CT_k3p5, CRS('+init=epsg:26917'))
CNY_k3p5 = spTransform(CNY_k3p5, CRS('+init=epsg:26917'))
CWV_k3p5 = spTransform(CWV_k3p5, CRS('+init=epsg:26917'))
ENY_k3p5 = spTransform(ENY_k3p5, CRS('+init=epsg:26917'))
ENYPA_k3p5 = spTransform(ENYPA_k3p5, CRS('+init=epsg:26917'))
MT_k3p5 = spTransform(MT_k3p5, CRS('+init=epsg:26917'))
NWPANY_k3p5 = spTransform(NWPANY_k3p5, CRS('+init=epsg:26917'))
SWPA_k3p5 = spTransform(SWPA_k3p5, CRS('+init=epsg:26917'))
WPA_k3p5 = spTransform(WPA_k3p5, CRS('+init=epsg:26917'))
VR_k3p5 = spTransform(VR_k3p5, CRS('+init=epsg:26917'))
FL_k3p5 = spTransform(FL_k3p5, CRS('+init=epsg:26917'))

CT_k4p0 = spTransform(CT_k4p0, CRS('+init=epsg:26917'))
CNY_k4p0 = spTransform(CNY_k4p0, CRS('+init=epsg:26917'))
CWV_k4p0 = spTransform(CWV_k4p0, CRS('+init=epsg:26917'))
ENY_k4p0 = spTransform(ENY_k4p0, CRS('+init=epsg:26917'))
ENYPA_k4p0 = spTransform(ENYPA_k4p0, CRS('+init=epsg:26917'))
MT_k4p0 = spTransform(MT_k4p0, CRS('+init=epsg:26917'))
NWPANY_k4p0 = spTransform(NWPANY_k4p0, CRS('+init=epsg:26917'))
SWPA_k4p0 = spTransform(SWPA_k4p0, CRS('+init=epsg:26917'))
WPA_k4p0 = spTransform(WPA_k4p0, CRS('+init=epsg:26917'))
VR_k4p0 = spTransform(VR_k4p0, CRS('+init=epsg:26917'))
FL_k4p0 = spTransform(FL_k4p0, CRS('+init=epsg:26917'))

CT_k4p5 = spTransform(CT_k4p5, CRS('+init=epsg:26917'))
CNY_k4p5 = spTransform(CNY_k4p5, CRS('+init=epsg:26917'))
CWV_k4p5 = spTransform(CWV_k4p5, CRS('+init=epsg:26917'))
ENY_k4p5 = spTransform(ENY_k4p5, CRS('+init=epsg:26917'))
ENYPA_k4p5 = spTransform(ENYPA_k4p5, CRS('+init=epsg:26917'))
MT_k4p5 = spTransform(MT_k4p5, CRS('+init=epsg:26917'))
NWPANY_k4p5 = spTransform(NWPANY_k4p5, CRS('+init=epsg:26917'))
SWPA_k4p5 = spTransform(SWPA_k4p5, CRS('+init=epsg:26917'))
WPA_k4p5 = spTransform(WPA_k4p5, CRS('+init=epsg:26917'))
VR_k4p5 = spTransform(VR_k4p5, CRS('+init=epsg:26917'))
FL_k4p5 = spTransform(FL_k4p5, CRS('+init=epsg:26917'))

CT_k5p0 = spTransform(CT_k5p0, CRS('+init=epsg:26917'))
CNY_k5p0 = spTransform(CNY_k5p0, CRS('+init=epsg:26917'))
CWV_k5p0 = spTransform(CWV_k5p0, CRS('+init=epsg:26917'))
ENY_k5p0 = spTransform(ENY_k5p0, CRS('+init=epsg:26917'))
ENYPA_k5p0 = spTransform(ENYPA_k5p0, CRS('+init=epsg:26917'))
MT_k5p0 = spTransform(MT_k5p0, CRS('+init=epsg:26917'))
NWPANY_k5p0 = spTransform(NWPANY_k5p0, CRS('+init=epsg:26917'))
SWPA_k5p0 = spTransform(SWPA_k5p0, CRS('+init=epsg:26917'))
WPA_k5p0 = spTransform(WPA_k5p0, CRS('+init=epsg:26917'))
VR_k5p0 = spTransform(VR_k5p0, CRS('+init=epsg:26917'))
FL_k5p0 = spTransform(FL_k5p0, CRS('+init=epsg:26917'))

#   Get dataset pre-processing into geologic regions----
#Using points without negative gradients, and no points in same spatial location
PreCT = spTransform(WellsSort[InterpRegs[InterpRegs$Name == 'CT',],], CRS('+init=epsg:26917'))
PreCNY = spTransform(WellsSort[InterpRegs[InterpRegs$Name == 'CNY',],], CRS('+init=epsg:26917'))
PreCWV = spTransform(WellsSort[CWV_Bounded,], CRS('+init=epsg:26917'))
PreENY = spTransform(WellsSort[InterpRegs[InterpRegs$Name == 'ENY',],], CRS('+init=epsg:26917'))
PreENYPA = spTransform(WellsSort[InterpRegs[InterpRegs$Name == 'ENYPA',],], CRS('+init=epsg:26917'))
PreMT = spTransform(WellsSort[MT_Bounded,], CRS('+init=epsg:26917'))
PreNWPANY = spTransform(WellsSort[InterpRegs[InterpRegs$Name == 'NWPANY',],], CRS('+init=epsg:26917'))
PreSWPA = spTransform(WellsSort[InterpRegs[InterpRegs$Name == 'SWPA',],], CRS('+init=epsg:26917'))
PreWPA = spTransform(WellsSort[InterpRegs[InterpRegs$Name == 'WPA',],], CRS('+init=epsg:26917'))
PreVR = spTransform(WellsSort[VR_Bounded,], CRS('+init=epsg:26917'))
PreFL = rbind(PreCT, PreCWV, PreCNY, PreENY, PreENYPA, PreMT, PreNWPANY, PreSWPA, PreWPA, PreVR)

#   Make variograms for each ESDA proceedure----
#Compute variograms - 1) All data, 2) Data deeper than 1 km, 3) Data that have been fully proessed. Plot all on same plot

#Possible parameters for Universal Kriging / regression kriging are:
# basement depth, surface temperature, BHT correction region, COSUNA section, Rome trough, spatial coordinates
# These are challenging to select because a trend in space can be captured by most of those parameters
#Use universal with Basement Depth for CT and WPA. Out to 30 km, not necessary to use universal for WPA.
#CT - COSUNA_ID explains nearly same variance as BasementDepth
#MT, CNY, SWPA do not need universal. NWPANY pobably not
#SWPA could use BHT correction region, but doesn't explain too much, which is acutally good.
#ENY universal with Basement Depth and COSUNA_ID. COSUNA is indicative of issues with local strat columns
#ENYPA with coordinates and Basement depth
#CWV COSUNA ID
v.CT <- variogram(Qs~1, CT, cutoff=60000, width=60000/50)
v.CNY <- variogram(Qs~1, CNY, cutoff=60000, width=60000/15)
v.CWV <- variogram(Qs~1, CWV, cutoff=60000, width=60000/50)
v.ENY <- variogram(Qs~1, ENY, cutoff=60000, width=60000/15)
v.ENYPA <- variogram(Qs~1, ENYPA, cutoff=60000, width=60000/40)
v.MT <- variogram(Qs~1, MT, cutoff=60000, width=60000/50)
v.NWPANY <- variogram(Qs~1, NWPANY, cutoff=60000, width=60000/20) 
v.SWPA <- variogram(Qs~1, SWPA, cutoff=60000, width=60000/50) 
v.WPA <- variogram(Qs~1, WPA, cutoff=60000, width=60000/50) 
v.VR <- variogram(Qs~1, VR, cutoff=60000, width=60000/20) 
v.FL <- variogram(Qs~1, FL, cutoff=60000, width=60000/200)

v.CT_RmOps <- variogram(Qs~1, CT_RmOps, cutoff=60000, width=60000/50)
v.CNY_RmOps <- variogram(Qs~1, CNY_RmOps, cutoff=60000, width=60000/15)
v.CWV_RmOps <- variogram(Qs~1, CWV_RmOps, cutoff=60000, width=60000/50)
v.ENY_RmOps <- variogram(Qs~1, ENY_RmOps, cutoff=60000, width=60000/15)
v.ENYPA_RmOps <- variogram(Qs~1, ENYPA_RmOps, cutoff=60000, width=60000/40)
v.MT_RmOps <- variogram(Qs~1, MT_RmOps, cutoff=60000, width=60000/50)
v.NWPANY_RmOps <- variogram(Qs~1, NWPANY_RmOps, cutoff=60000, width=60000/20) 
v.SWPA_RmOps <- variogram(Qs~1, SWPA_RmOps, cutoff=60000, width=60000/50) 
v.WPA_RmOps <- variogram(Qs~1, WPA_RmOps, cutoff=60000, width=60000/50) 
v.VR_RmOps <- variogram(Qs~1, VR_RmOps, cutoff=60000, width=60000/20) 
v.FL_RmOps <- variogram(Qs~1, FL_RmOps, cutoff=60000, width=60000/200)

v.CT_k1p0 <- variogram(Qs~1, CT_k1p0, cutoff=60000, width=60000/50)
v.CNY_k1p0 <- variogram(Qs~1, CNY_k1p0, cutoff=60000, width=60000/15)
v.CWV_k1p0 <- variogram(Qs~1, CWV_k1p0, cutoff=60000, width=60000/50)
v.ENY_k1p0 <- variogram(Qs~1, ENY_k1p0, cutoff=60000, width=60000/15)
v.ENYPA_k1p0 <- variogram(Qs~1, ENYPA_k1p0, cutoff=60000, width=60000/40)
v.MT_k1p0 <- variogram(Qs~1, MT_k1p0, cutoff=60000, width=60000/50)
v.NWPANY_k1p0 <- variogram(Qs~1, NWPANY_k1p0, cutoff=60000, width=60000/20) 
v.SWPA_k1p0 <- variogram(Qs~1, SWPA_k1p0, cutoff=60000, width=60000/50) 
v.WPA_k1p0 <- variogram(Qs~1, WPA_k1p0, cutoff=60000, width=60000/50) 
v.VR_k1p0 <- variogram(Qs~1, VR_k1p0, cutoff=60000, width=60000/20) 
v.FL_k1p0 <- variogram(Qs~1, FL_k1p0, cutoff=60000, width=60000/200)

v.CT_k1p5 <- variogram(Qs~1, CT_k1p5, cutoff=60000, width=60000/50)
v.CNY_k1p5 <- variogram(Qs~1, CNY_k1p5, cutoff=60000, width=60000/15)
v.CWV_k1p5 <- variogram(Qs~1, CWV_k1p5, cutoff=60000, width=60000/50)
v.ENY_k1p5 <- variogram(Qs~1, ENY_k1p5, cutoff=60000, width=60000/15)
v.ENYPA_k1p5 <- variogram(Qs~1, ENYPA_k1p5, cutoff=60000, width=60000/40)
v.MT_k1p5 <- variogram(Qs~1, MT_k1p5, cutoff=60000, width=60000/50)
v.NWPANY_k1p5 <- variogram(Qs~1, NWPANY_k1p5, cutoff=60000, width=60000/20) 
v.SWPA_k1p5 <- variogram(Qs~1, SWPA_k1p5, cutoff=60000, width=60000/50) 
v.WPA_k1p5 <- variogram(Qs~1, WPA_k1p5, cutoff=60000, width=60000/50) 
v.VR_k1p5 <- variogram(Qs~1, VR_k1p5, cutoff=60000, width=60000/20) 
v.FL_k1p5 <- variogram(Qs~1, FL_k1p5, cutoff=60000, width=60000/200)

v.CT_k2p0 <- variogram(Qs~1, CT_k2p0, cutoff=60000, width=60000/50)
v.CNY_k2p0 <- variogram(Qs~1, CNY_k2p0, cutoff=60000, width=60000/15)
v.CWV_k2p0 <- variogram(Qs~1, CWV_k2p0, cutoff=60000, width=60000/50)
v.ENY_k2p0 <- variogram(Qs~1, ENY_k2p0, cutoff=60000, width=60000/15)
v.ENYPA_k2p0 <- variogram(Qs~1, ENYPA_k2p0, cutoff=60000, width=60000/40)
v.MT_k2p0 <- variogram(Qs~1, MT_k2p0, cutoff=60000, width=60000/50)
v.NWPANY_k2p0 <- variogram(Qs~1, NWPANY_k2p0, cutoff=60000, width=60000/20) 
v.SWPA_k2p0 <- variogram(Qs~1, SWPA_k2p0, cutoff=60000, width=60000/50) 
v.WPA_k2p0 <- variogram(Qs~1, WPA_k2p0, cutoff=60000, width=60000/50) 
v.VR_k2p0 <- variogram(Qs~1, VR_k2p0, cutoff=60000, width=60000/20) 
v.FL_k2p0 <- variogram(Qs~1, FL_k2p0, cutoff=60000, width=60000/200)

v.CT_k2p5 <- variogram(Qs~1, CT_k2p5, cutoff=60000, width=60000/50)
v.CNY_k2p5 <- variogram(Qs~1, CNY_k2p5, cutoff=60000, width=60000/15)
v.CWV_k2p5 <- variogram(Qs~1, CWV_k2p5, cutoff=60000, width=60000/50)
v.ENY_k2p5 <- variogram(Qs~1, ENY_k2p5, cutoff=60000, width=60000/15)
v.ENYPA_k2p5 <- variogram(Qs~1, ENYPA_k2p5, cutoff=60000, width=60000/40)
v.MT_k2p5 <- variogram(Qs~1, MT_k2p5, cutoff=60000, width=60000/50)
v.NWPANY_k2p5 <- variogram(Qs~1, NWPANY_k2p5, cutoff=60000, width=60000/20) 
v.SWPA_k2p5 <- variogram(Qs~1, SWPA_k2p5, cutoff=60000, width=60000/50) 
v.WPA_k2p5 <- variogram(Qs~1, WPA_k2p5, cutoff=60000, width=60000/50) 
v.VR_k2p5 <- variogram(Qs~1, VR_k2p5, cutoff=60000, width=60000/20) 
v.FL_k2p5 <- variogram(Qs~1, FL_k2p5, cutoff=60000, width=60000/200)

v.CT_k3p5 <- variogram(Qs~1, CT_k3p5, cutoff=60000, width=60000/50)
v.CNY_k3p5 <- variogram(Qs~1, CNY_k3p5, cutoff=60000, width=60000/15)
v.CWV_k3p5 <- variogram(Qs~1, CWV_k3p5, cutoff=60000, width=60000/50)
v.ENY_k3p5 <- variogram(Qs~1, ENY_k3p5, cutoff=60000, width=60000/15)
v.ENYPA_k3p5 <- variogram(Qs~1, ENYPA_k3p5, cutoff=60000, width=60000/40)
v.MT_k3p5 <- variogram(Qs~1, MT_k3p5, cutoff=60000, width=60000/50)
v.NWPANY_k3p5 <- variogram(Qs~1, NWPANY_k3p5, cutoff=60000, width=60000/20) 
v.SWPA_k3p5 <- variogram(Qs~1, SWPA_k3p5, cutoff=60000, width=60000/50) 
v.WPA_k3p5 <- variogram(Qs~1, WPA_k3p5, cutoff=60000, width=60000/50) 
v.VR_k3p5 <- variogram(Qs~1, VR_k3p5, cutoff=60000, width=60000/20) 
v.FL_k3p5 <- variogram(Qs~1, FL_k3p5, cutoff=60000, width=60000/200)

v.CT_k4p0 <- variogram(Qs~1, CT_k4p0, cutoff=60000, width=60000/50)
v.CNY_k4p0 <- variogram(Qs~1, CNY_k4p0, cutoff=60000, width=60000/15)
v.CWV_k4p0 <- variogram(Qs~1, CWV_k4p0, cutoff=60000, width=60000/50)
v.ENY_k4p0 <- variogram(Qs~1, ENY_k4p0, cutoff=60000, width=60000/15)
v.ENYPA_k4p0 <- variogram(Qs~1, ENYPA_k4p0, cutoff=60000, width=60000/40)
v.MT_k4p0 <- variogram(Qs~1, MT_k4p0, cutoff=60000, width=60000/50)
v.NWPANY_k4p0 <- variogram(Qs~1, NWPANY_k4p0, cutoff=60000, width=60000/20) 
v.SWPA_k4p0 <- variogram(Qs~1, SWPA_k4p0, cutoff=60000, width=60000/50) 
v.WPA_k4p0 <- variogram(Qs~1, WPA_k4p0, cutoff=60000, width=60000/50) 
v.VR_k4p0 <- variogram(Qs~1, VR_k4p0, cutoff=60000, width=60000/20) 
v.FL_k4p0 <- variogram(Qs~1, FL_k4p0, cutoff=60000, width=60000/200)

v.CT_k4p5 <- variogram(Qs~1, CT_k4p5, cutoff=60000, width=60000/50)
v.CNY_k4p5 <- variogram(Qs~1, CNY_k4p5, cutoff=60000, width=60000/15)
v.CWV_k4p5 <- variogram(Qs~1, CWV_k4p5, cutoff=60000, width=60000/50)
v.ENY_k4p5 <- variogram(Qs~1, ENY_k4p5, cutoff=60000, width=60000/15)
v.ENYPA_k4p5 <- variogram(Qs~1, ENYPA_k4p5, cutoff=60000, width=60000/40)
v.MT_k4p5 <- variogram(Qs~1, MT_k4p5, cutoff=60000, width=60000/50)
v.NWPANY_k4p5 <- variogram(Qs~1, NWPANY_k4p5, cutoff=60000, width=60000/20) 
v.SWPA_k4p5 <- variogram(Qs~1, SWPA_k4p5, cutoff=60000, width=60000/50) 
v.WPA_k4p5 <- variogram(Qs~1, WPA_k4p5, cutoff=60000, width=60000/50) 
v.VR_k4p5 <- variogram(Qs~1, VR_k4p5, cutoff=60000, width=60000/20) 
v.FL_k4p5 <- variogram(Qs~1, FL_k4p5, cutoff=60000, width=60000/200)

v.CT_k5p0 <- variogram(Qs~1, CT_k5p0, cutoff=60000, width=60000/50)
v.CNY_k5p0 <- variogram(Qs~1, CNY_k5p0, cutoff=60000, width=60000/15)
v.CWV_k5p0 <- variogram(Qs~1, CWV_k5p0, cutoff=60000, width=60000/50)
v.ENY_k5p0 <- variogram(Qs~1, ENY_k5p0, cutoff=60000, width=60000/15)
v.ENYPA_k5p0 <- variogram(Qs~1, ENYPA_k5p0, cutoff=60000, width=60000/40)
v.MT_k5p0 <- variogram(Qs~1, MT_k5p0, cutoff=60000, width=60000/50)
v.NWPANY_k5p0 <- variogram(Qs~1, NWPANY_k5p0, cutoff=60000, width=60000/20) 
v.SWPA_k5p0 <- variogram(Qs~1, SWPA_k5p0, cutoff=60000, width=60000/50) 
v.WPA_k5p0 <- variogram(Qs~1, WPA_k5p0, cutoff=60000, width=60000/50) 
v.VR_k5p0 <- variogram(Qs~1, VR_k5p0, cutoff=60000, width=60000/20) 
v.FL_k5p0 <- variogram(Qs~1, FL_k5p0, cutoff=60000, width=60000/200)

v.PreCT <- variogram(Qs~1, PreCT, cutoff=60000, width=60000/50) 
v.PreCNY <- variogram(Qs~1, PreCNY, cutoff=60000, width=60000/15) 
v.PreCWV <- variogram(Qs~1, PreCWV, cutoff=60000, width=60000/50)
v.PreENY <- variogram(Qs~1, PreENY, cutoff=60000, width=60000/15)
v.PreENYPA <- variogram(Qs~1, PreENYPA, cutoff=60000, width=60000/40)
v.PreMT <- variogram(Qs~1, PreMT, cutoff=60000, width=60000/50)
v.PreNWPANY <- variogram(Qs~1, PreNWPANY, cutoff=60000, width=60000/20) 
v.PreSWPA <- variogram(Qs~1, PreSWPA, cutoff=60000, width=60000/50) 
v.PreWPA <- variogram(Qs~1, PreWPA, cutoff=60000, width=60000/50) 
v.PreVR <- variogram(Qs~1, PreVR, cutoff=60000, width=60000/20) 
v.PreFL <- variogram(Qs~1, PreFL, cutoff=60000, width=60000/200)

v.DeepCT <- variogram(Qs~1, DeepCT, cutoff=60000, width=60000/50) 
v.DeepCNY <- variogram(Qs~1, DeepCNY, cutoff=60000, width=60000/15) 
v.DeepCWV <- variogram(Qs~1, DeepCWV, cutoff=60000, width=60000/50)
v.DeepENY <- variogram(Qs~1, DeepENY, cutoff=60000, width=60000/15)
v.DeepENYPA <- variogram(Qs~1, DeepENYPA, cutoff=60000, width=60000/40)
v.DeepMT <- variogram(Qs~1, DeepMT, cutoff=60000, width=60000/50)
v.DeepNWPANY <- variogram(Qs~1, DeepNWPANY, cutoff=60000, width=60000/20) 
v.DeepSWPA <- variogram(Qs~1, DeepSWPA, cutoff=60000, width=60000/50) 
v.DeepWPA <- variogram(Qs~1, DeepWPA, cutoff=60000, width=60000/50) 
v.DeepVR <- variogram(Qs~1, DeepVR, cutoff=60000, width=60000/20) 
v.DeepFL <- variogram(Qs~1, DeepFL, cutoff=60000, width=60000/200)

rv.DeepCT <- variogram(Qs~1, DeepCT, cutoff=60000, width=60000/50, cressie = T) 
rv.DeepCNY <- variogram(Qs~1, DeepCNY, cutoff=60000, width=60000/15, cressie = T) 
rv.DeepCWV <- variogram(Qs~1, DeepCWV, cutoff=60000, width=60000/50, cressie = T)
rv.DeepENY <- variogram(Qs~1, DeepENY, cutoff=60000, width=60000/15, cressie = T)
rv.DeepENYPA <- variogram(Qs~1, DeepENYPA, cutoff=60000, width=60000/40, cressie = T)
rv.DeepMT <- variogram(Qs~1, DeepMT, cutoff=60000, width=60000/50, cressie = T)
rv.DeepNWPANY <- variogram(Qs~1, DeepNWPANY, cutoff=60000, width=60000/20, cressie = T) 
rv.DeepSWPA <- variogram(Qs~1, DeepSWPA, cutoff=60000, width=60000/50, cressie = T) 
rv.DeepWPA <- variogram(Qs~1, DeepWPA, cutoff=60000, width=60000/50, cressie = T) 
rv.DeepVR <- variogram(Qs~1, DeepVR, cutoff=60000, width=60000/20, cressie = T) 
rv.DeepFL <- variogram(Qs~1, DeepFL, cutoff=60000, width=60000/200, cressie = T)

v.DeepCT_RmOps <- variogram(Qs~1, DeepCT_RmOps, cutoff=60000, width=60000/50) 
v.DeepCNY_RmOps <- variogram(Qs~1, DeepCNY_RmOps, cutoff=60000, width=60000/15) 
v.DeepCWV_RmOps <- variogram(Qs~1, DeepCWV_RmOps, cutoff=60000, width=60000/50)
v.DeepENY_RmOps <- variogram(Qs~1, DeepENY_RmOps, cutoff=60000, width=60000/15)
v.DeepENYPA_RmOps <- variogram(Qs~1, DeepENYPA_RmOps, cutoff=60000, width=60000/40)
v.DeepMT_RmOps <- variogram(Qs~1, DeepMT_RmOps, cutoff=60000, width=60000/50)
v.DeepNWPANY_RmOps <- variogram(Qs~1, DeepNWPANY_RmOps, cutoff=60000, width=60000/20) 
v.DeepSWPA_RmOps <- variogram(Qs~1, DeepSWPA_RmOps, cutoff=60000, width=60000/50) 
v.DeepWPA_RmOps <- variogram(Qs~1, DeepWPA_RmOps, cutoff=60000, width=60000/50) 
v.DeepVR_RmOps <- variogram(Qs~1, DeepVR_RmOps, cutoff=60000, width=60000/20) 
v.DeepFL_RmOps <- variogram(Qs~1, DeepFL_RmOps, cutoff=60000, width=60000/200)

rv.DeepCT_RmOps <- variogram(Qs~1, DeepCT_RmOps, cutoff=60000, width=60000/50, cressie = T) 
rv.DeepCNY_RmOps <- variogram(Qs~1, DeepCNY_RmOps, cutoff=60000, width=60000/15, cressie = T) 
rv.DeepCWV_RmOps <- variogram(Qs~1, DeepCWV_RmOps, cutoff=60000, width=60000/50, cressie = T)
rv.DeepENY_RmOps <- variogram(Qs~1, DeepENY_RmOps, cutoff=60000, width=60000/15, cressie = T)
rv.DeepENYPA_RmOps <- variogram(Qs~1, DeepENYPA_RmOps, cutoff=60000, width=60000/40, cressie = T)
rv.DeepMT_RmOps <- variogram(Qs~1, DeepMT_RmOps, cutoff=60000, width=60000/50, cressie = T)
rv.DeepNWPANY_RmOps <- variogram(Qs~1, DeepNWPANY_RmOps, cutoff=60000, width=60000/20, cressie = T) 
rv.DeepSWPA_RmOps <- variogram(Qs~1, DeepSWPA_RmOps, cutoff=60000, width=60000/50, cressie = T) 
rv.DeepWPA_RmOps <- variogram(Qs~1, DeepWPA_RmOps, cutoff=60000, width=60000/50, cressie = T) 
rv.DeepVR_RmOps <- variogram(Qs~1, DeepVR_RmOps, cutoff=60000, width=60000/20, cressie = T) 
rv.DeepFL_RmOps <- variogram(Qs~1, DeepFL_RmOps, cutoff=60000, width=60000/200, cressie = T)

p1 = plot(v.CT, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Chautauqua, NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,50), col = 'red', xlim = c(0,60000), cex=0.5)
p2 = plot(v.CNY,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,1500), col = 'red', xlim = c(0,60000), cex=0.5)
p9 = plot(v.CWV, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'red', ylim = c(0,1200), xlim = c(0,60000), cex=0.5)
p3 = plot(v.ENY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'red', ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p4 = plot(v.ENYPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY & PA", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,400), cex=0.5)
p8 = plot(v.MT, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p5 = plot(v.NWPANY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,500), cex=0.5)
p6 = plot(v.SWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,850), cex=0.5)
p6z = plot(v.SWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p6zz = plot(v.SWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p10 = plot(v.VR, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p10z = plot(v.VR, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p10zz = plot(v.VR, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p7 = plot(v.WPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p7z = plot(v.WPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p7zz = plot(v.WPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p11 = plot(v.FL, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,1000), cex=0.5)

p1RmOps = plot(v.CT_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Chautauqua, NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,50), col = 'red', xlim = c(0,60000), cex=0.5)
p2RmOps = plot(v.CNY_RmOps,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,1500), col = 'red', xlim = c(0,60000), cex=0.5)
p9RmOps = plot(v.CWV_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'red', ylim = c(0,1200), xlim = c(0,60000), cex=0.5)
p3RmOps = plot(v.ENY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'red', ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p4RmOps = plot(v.ENYPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY & PA", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,400), cex=0.5)
p8RmOps = plot(v.MT_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p5RmOps = plot(v.NWPANY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,500), cex=0.5)
p6RmOps = plot(v.SWPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,850), cex=0.5)
p6zRmOps = plot(v.SWPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p6zzRmOps = plot(v.SWPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p10RmOps = plot(v.VR_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p10zRmOps = plot(v.VR_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p10zzRmOps = plot(v.VR_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p7RmOps = plot(v.WPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p7zRmOps = plot(v.WPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p7zzRmOps = plot(v.WPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p11RmOps = plot(v.FL_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,1000), cex=0.5)


kcols = c(adjustcolor('green', alpha.f = 0.1),adjustcolor('green', alpha.f = 0.2),adjustcolor('green', alpha.f = 0.3),adjustcolor('green', alpha.f = 0.4),adjustcolor('green', alpha.f = 0.6),adjustcolor('green', alpha.f = 0.7),adjustcolor('green', alpha.f = 0.8),adjustcolor('green', alpha.f = 0.9))

p1_k1p0 = plot(v.CT_k1p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Chautauqua, NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,50), col = kcols[1], xlim = c(0,60000), cex=0.5)
p2_k1p0 = plot(v.CNY_k1p0,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,1500), col = kcols[1], xlim = c(0,60000), cex=0.5)
p9_k1p0 = plot(v.CWV_k1p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[1], ylim = c(0,1200), xlim = c(0,60000), cex=0.5)
p3_k1p0 = plot(v.ENY_k1p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[1], ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p4_k1p0 = plot(v.ENYPA_k1p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY & PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[1], xlim = c(0,60000), ylim = c(0,400), cex=0.5)
p8_k1p0 = plot(v.MT_k1p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[1], xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p5_k1p0 = plot(v.NWPANY_k1p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[1], xlim = c(0,60000), ylim = c(0,500), cex=0.5)
p6_k1p0 = plot(v.SWPA_k1p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[1], xlim = c(0,60000), ylim = c(0,850), cex=0.5)
p6z_k1p0 = plot(v.SWPA_k1p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[1], xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p10_k1p0 = plot(v.VR_k1p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = kcols[1], xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p10z_k1p0 = plot(v.VR_k1p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = kcols[1], xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p7_k1p0 = plot(v.WPA_k1p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[1], xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p7z_k1p0 = plot(v.WPA_k1p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[1], xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p11_k1p0 = plot(v.FL_k1p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = kcols[1], xlim = c(0,60000), ylim = c(0,1000), cex=0.5)

p1_k1p5 = plot(v.CT_k1p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Chautauqua, NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,50), col = kcols[2], xlim = c(0,60000), cex=0.5)
p2_k1p5 = plot(v.CNY_k1p5,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,1500), col = kcols[2], xlim = c(0,60000), cex=0.5)
p9_k1p5 = plot(v.CWV_k1p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[2], ylim = c(0,1200), xlim = c(0,60000), cex=0.5)
p3_k1p5 = plot(v.ENY_k1p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[2], ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p4_k1p5 = plot(v.ENYPA_k1p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY & PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[2], xlim = c(0,60000), ylim = c(0,400), cex=0.5)
p8_k1p5 = plot(v.MT_k1p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[2], xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p5_k1p5 = plot(v.NWPANY_k1p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[2], xlim = c(0,60000), ylim = c(0,500), cex=0.5)
p6_k1p5 = plot(v.SWPA_k1p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[2], xlim = c(0,60000), ylim = c(0,850), cex=0.5)
p6z_k1p5 = plot(v.SWPA_k1p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[2], xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p10_k1p5 = plot(v.VR_k1p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = kcols[2], xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p10z_k1p5 = plot(v.VR_k1p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = kcols[2], xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p7_k1p5 = plot(v.WPA_k1p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[2], xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p7z_k1p5 = plot(v.WPA_k1p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[2], xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p11_k1p5 = plot(v.FL_k1p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = kcols[2], xlim = c(0,60000), ylim = c(0,1000), cex=0.5)

p1_k2p0 = plot(v.CT_k2p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Chautauqua, NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,50), col = kcols[3], xlim = c(0,60000), cex=0.5)
p2_k2p0 = plot(v.CNY_k2p0,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,1500), col = kcols[3], xlim = c(0,60000), cex=0.5)
p9_k2p0 = plot(v.CWV_k2p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[3], ylim = c(0,1200), xlim = c(0,60000), cex=0.5)
p3_k2p0 = plot(v.ENY_k2p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[3], ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p4_k2p0 = plot(v.ENYPA_k2p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY & PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[3], xlim = c(0,60000), ylim = c(0,400), cex=0.5)
p8_k2p0 = plot(v.MT_k2p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[3], xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p5_k2p0 = plot(v.NWPANY_k2p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[3], xlim = c(0,60000), ylim = c(0,500), cex=0.5)
p6_k2p0 = plot(v.SWPA_k2p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[3], xlim = c(0,60000), ylim = c(0,850), cex=0.5)
p6z_k2p0 = plot(v.SWPA_k2p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[3], xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p10_k2p0 = plot(v.VR_k2p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = kcols[3], xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p10z_k2p0 = plot(v.VR_k2p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = kcols[3], xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p7_k2p0 = plot(v.WPA_k2p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[3], xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p7z_k2p0 = plot(v.WPA_k2p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[3], xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p11_k2p0 = plot(v.FL_k2p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = kcols[3], xlim = c(0,60000), ylim = c(0,1000), cex=0.5)

p1_k2p5 = plot(v.CT_k2p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Chautauqua, NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,50), col = kcols[4], xlim = c(0,60000), cex=0.5)
p2_k2p5 = plot(v.CNY_k2p5,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,1500), col = kcols[4], xlim = c(0,60000), cex=0.5)
p9_k2p5 = plot(v.CWV_k2p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[4], ylim = c(0,1200), xlim = c(0,60000), cex=0.5)
p3_k2p5 = plot(v.ENY_k2p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[4], ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p4_k2p5 = plot(v.ENYPA_k2p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY & PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[4], xlim = c(0,60000), ylim = c(0,400), cex=0.5)
p8_k2p5 = plot(v.MT_k2p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[4], xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p5_k2p5 = plot(v.NWPANY_k2p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[4], xlim = c(0,60000), ylim = c(0,500), cex=0.5)
p6_k2p5 = plot(v.SWPA_k2p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[4], xlim = c(0,60000), ylim = c(0,850), cex=0.5)
p6z_k2p5 = plot(v.SWPA_k2p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[4], xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p10_k2p5 = plot(v.VR_k2p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = kcols[4], xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p10z_k2p5 = plot(v.VR_k2p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = kcols[4], xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p7_k2p5 = plot(v.WPA_k2p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[4], xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p7z_k2p5 = plot(v.WPA_k2p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[4], xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p11_k2p5 = plot(v.FL_k2p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = kcols[4], xlim = c(0,60000), ylim = c(0,1000), cex=0.5)

p1_k3p5 = plot(v.CT_k3p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Chautauqua, NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,50), col = kcols[5], xlim = c(0,60000), cex=0.5)
p2_k3p5 = plot(v.CNY_k3p5,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,1500), col = kcols[5], xlim = c(0,60000), cex=0.5)
p9_k3p5 = plot(v.CWV_k3p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[5], ylim = c(0,1200), xlim = c(0,60000), cex=0.5)
p3_k3p5 = plot(v.ENY_k3p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[5], ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p4_k3p5 = plot(v.ENYPA_k3p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY & PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[5], xlim = c(0,60000), ylim = c(0,400), cex=0.5)
p8_k3p5 = plot(v.MT_k3p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[5], xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p5_k3p5 = plot(v.NWPANY_k3p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[5], xlim = c(0,60000), ylim = c(0,500), cex=0.5)
p6_k3p5 = plot(v.SWPA_k3p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[5], xlim = c(0,60000), ylim = c(0,850), cex=0.5)
p6z_k3p5 = plot(v.SWPA_k3p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[5], xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p10_k3p5 = plot(v.VR_k3p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = kcols[5], xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p10z_k3p5 = plot(v.VR_k3p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = kcols[5], xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p7_k3p5 = plot(v.WPA_k3p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[5], xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p7z_k3p5 = plot(v.WPA_k3p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[5], xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p11_k3p5 = plot(v.FL_k3p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = kcols[5], xlim = c(0,60000), ylim = c(0,1000), cex=0.5)

p1_k4p0 = plot(v.CT_k4p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Chautauqua, NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,50), col = kcols[6], xlim = c(0,60000), cex=0.5)
p2_k4p0 = plot(v.CNY_k4p0,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,1500), col = kcols[6], xlim = c(0,60000), cex=0.5)
p9_k4p0 = plot(v.CWV_k4p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[6], ylim = c(0,1200), xlim = c(0,60000), cex=0.5)
p3_k4p0 = plot(v.ENY_k4p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[6], ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p4_k4p0 = plot(v.ENYPA_k4p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY & PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[6], xlim = c(0,60000), ylim = c(0,400), cex=0.5)
p8_k4p0 = plot(v.MT_k4p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[6], xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p5_k4p0 = plot(v.NWPANY_k4p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[6], xlim = c(0,60000), ylim = c(0,500), cex=0.5)
p6_k4p0 = plot(v.SWPA_k4p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[6], xlim = c(0,60000), ylim = c(0,850), cex=0.5)
p6z_k4p0 = plot(v.SWPA_k4p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[6], xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p10_k4p0 = plot(v.VR_k4p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = kcols[6], xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p10z_k4p0 = plot(v.VR_k4p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = kcols[6], xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p7_k4p0 = plot(v.WPA_k4p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[6], xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p7z_k4p0 = plot(v.WPA_k4p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[6], xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p11_k4p0 = plot(v.FL_k4p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = kcols[6], xlim = c(0,60000), ylim = c(0,1000), cex=0.5)

p1_k4p5 = plot(v.CT_k4p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Chautauqua, NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,50), col = kcols[7], xlim = c(0,60000), cex=0.5)
p2_k4p5 = plot(v.CNY_k4p5,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,1500), col = kcols[7], xlim = c(0,60000), cex=0.5)
p9_k4p5 = plot(v.CWV_k4p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[7], ylim = c(0,1200), xlim = c(0,60000), cex=0.5)
p3_k4p5 = plot(v.ENY_k4p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[7], ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p4_k4p5 = plot(v.ENYPA_k4p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY & PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[7], xlim = c(0,60000), ylim = c(0,400), cex=0.5)
p8_k4p5 = plot(v.MT_k4p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[7], xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p5_k4p5 = plot(v.NWPANY_k4p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[7], xlim = c(0,60000), ylim = c(0,500), cex=0.5)
p6_k4p5 = plot(v.SWPA_k4p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[7], xlim = c(0,60000), ylim = c(0,850), cex=0.5)
p6z_k4p5 = plot(v.SWPA_k4p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[7], xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p10_k4p5 = plot(v.VR_k4p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = kcols[7], xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p10z_k4p5 = plot(v.VR_k4p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = kcols[7], xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p7_k4p5 = plot(v.WPA_k4p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[7], xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p7z_k4p5 = plot(v.WPA_k4p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[7], xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p11_k4p5 = plot(v.FL_k4p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = kcols[7], xlim = c(0,60000), ylim = c(0,1000), cex=0.5)

p1_k5p0 = plot(v.CT_k5p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Chautauqua, NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,50), col = kcols[8], xlim = c(0,60000), cex=0.5)
p2_k5p0 = plot(v.CNY_k5p0,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,1500), col = kcols[8], xlim = c(0,60000), cex=0.5)
p9_k5p0 = plot(v.CWV_k5p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[8], ylim = c(0,1200), xlim = c(0,60000), cex=0.5)
p3_k5p0 = plot(v.ENY_k5p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[8], ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p4_k5p0 = plot(v.ENYPA_k5p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY & PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[8], xlim = c(0,60000), ylim = c(0,400), cex=0.5)
p8_k5p0 = plot(v.MT_k5p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[8], xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p5_k5p0 = plot(v.NWPANY_k5p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[8], xlim = c(0,60000), ylim = c(0,500), cex=0.5)
p6_k5p0 = plot(v.SWPA_k5p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[8], xlim = c(0,60000), ylim = c(0,850), cex=0.5)
p6z_k5p0 = plot(v.SWPA_k5p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[8], xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p10_k5p0 = plot(v.VR_k5p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = kcols[8], xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p10z_k5p0 = plot(v.VR_k5p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = kcols[8], xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p7_k5p0 = plot(v.WPA_k5p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[8], xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p7z_k5p0 = plot(v.WPA_k5p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = kcols[8], xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p11_k5p0 = plot(v.FL_k5p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = kcols[8], xlim = c(0,60000), ylim = c(0,1000), cex=0.5)

p1p = plot(v.PreCT, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Chautauqua, NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,50), col = 'black', xlim = c(0,60000), cex=0.5)
p2p = plot(v.PreCNY,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, col = 'black', ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p9p = plot(v.PreCWV, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'black', ylim = c(0,1200), xlim = c(0,60000), cex=0.5)
p3p = plot(v.PreENY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'black', ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p4p = plot(v.PreENYPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY & PA", xlab = 'Separation Distance (m)', pch = 16, col = 'black', xlim = c(0,60000), ylim = c(0,400), cex=0.5)
p8p = plot(v.PreMT, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = 'black', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p5p = plot(v.PreNWPANY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'black', xlim = c(0,60000), ylim = c(0,500), cex=0.5)
p6p = plot(v.PreSWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = 'black', xlim = c(0,60000), ylim = c(0,850), cex=0.5)
p10p = plot(v.PreVR, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = 'black', xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p7p = plot(v.PreWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = 'black', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p11p = plot(v.PreFL, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = 'black', xlim = c(0,60000), ylim = c(0,1000), cex=0.5)

p1d = plot(v.DeepCT, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Chautauqua, NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,50), col = 'blue', xlim = c(0,60000), cex=0.5)
p2d = plot(v.DeepCNY,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p9d = plot(v.DeepCWV, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,1200), xlim = c(0,60000), cex=0.5)
p3d = plot(v.DeepENY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p4d = plot(v.DeepENYPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY & PA", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,400), cex=0.5)
p8d = plot(v.DeepMT, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p5d = plot(v.DeepNWPANY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,500), cex=0.5)
p6d = plot(v.DeepSWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,850), cex=0.5)
p6dz = plot(v.DeepSWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p6dzz = plot(v.DeepSWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p10d = plot(v.DeepVR, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p10dz = plot(v.DeepVR, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p10dzz = plot(v.DeepVR, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p7d = plot(v.DeepWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p7dz = plot(v.DeepWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p7dzz = plot(v.DeepWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p11d = plot(v.DeepFL, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,1000), cex=0.5)

p1dRmOps = plot(v.DeepCT_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Chautauqua, NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,50), col = 'green', xlim = c(0,60000), cex=0.5)
p2dRmOps = plot(v.DeepCNY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, col = 'green', ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p9dRmOps = plot(v.DeepCWV_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'green', ylim = c(0,1200), xlim = c(0,60000), cex=0.5)
p3dRmOps = plot(v.DeepENY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'green', ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p4dRmOps = plot(v.DeepENYPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY & PA", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,400), cex=0.5)
p8dRmOps = plot(v.DeepMT_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p5dRmOps = plot(v.DeepNWPANY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,500), cex=0.5)
p6dRmOps = plot(v.DeepSWPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,850), cex=0.5)
p6dzRmOps = plot(v.DeepSWPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p6dzzRmOps = plot(v.DeepSWPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p10dRmOps = plot(v.DeepVR_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p10dzRmOps = plot(v.DeepVR_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p10dzzRmOps = plot(v.DeepVR_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p7dRmOps = plot(v.DeepWPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p7dzRmOps = plot(v.DeepWPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p7dzzRmOps = plot(v.DeepWPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p11dRmOps = plot(v.DeepFL_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,1000), cex=0.5)

p1dr = plot(rv.DeepCT, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Chautauqua, NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,50), col = 'blue', xlim = c(0,60000), cex=0.5)
p2dr = plot(rv.DeepCNY,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p9dr = plot(rv.DeepCWV, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,1200), xlim = c(0,60000), cex=0.5)
p3dr = plot(rv.DeepENY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p4dr = plot(rv.DeepENYPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY & PA", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,400), cex=0.5)
p8dr = plot(rv.DeepMT, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p5dr = plot(rv.DeepNWPANY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,500), cex=0.5)
p6dr = plot(rv.DeepSWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,850), cex=0.5)
p6dzr = plot(rv.DeepSWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p6dzzr = plot(rv.DeepSWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p10dr = plot(rv.DeepVR, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p10dzr = plot(rv.DeepVR, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p10dzzr = plot(rv.DeepVR, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p7dr = plot(rv.DeepWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p7dzr = plot(rv.DeepWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p7dzzr = plot(rv.DeepWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p11dr = plot(rv.DeepFL, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,1000), cex=0.5)

p1dRmOpsr = plot(rv.DeepCT_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Chautauqua, NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,50), col = 'green', xlim = c(0,60000), cex=0.5)
p2dRmOpsr = plot(rv.DeepCNY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, col = 'green', ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p9dRmOpsr = plot(rv.DeepCWV_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'green', ylim = c(0,1200), xlim = c(0,60000), cex=0.5)
p3dRmOpsr = plot(rv.DeepENY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'green', ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p4dRmOpsr = plot(rv.DeepENYPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY & PA", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,400), cex=0.5)
p8dRmOpsr = plot(rv.DeepMT_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p5dRmOpsr = plot(rv.DeepNWPANY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,500), cex=0.5)
p6dRmOpsr = plot(rv.DeepSWPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,850), cex=0.5)
p6dzRmOpsr = plot(rv.DeepSWPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p6dzzRmOpsr = plot(rv.DeepSWPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p10dRmOpsr = plot(rv.DeepVR_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p10dzRmOpsr = plot(rv.DeepVR_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p10dzzRmOpsr = plot(rv.DeepVR_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p7dRmOpsr = plot(rv.DeepWPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p7dzRmOpsr = plot(rv.DeepWPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p7dzzRmOpsr = plot(rv.DeepWPA_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p11dRmOpsr = plot(rv.DeepFL_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,1000), cex=0.5)


png("Variograms_UniqueYAxis_CompareESDA.png", width=13, height=10, units="in", res=300)
plot(p1p, split=c(1,1,4,3), more=T)
plot(p1d, split=c(1,1,4,3), more=T)
plot(p1, split = c(1,1,4,3), more=T)

plot(p2p, split=c(4,1,4,3), more=T)
plot(p2d, split=c(4,1,4,3), more=T)
plot(p2, split = c(4,1,4,3), more=T)

plot(p3p, split=c(1,2,4,3), more=T)
plot(p3d, split=c(1,2,4,3), more=T)
plot(p3, split = c(1,2,4,3), more=T)

plot(p4p, split=c(2,2,4,3), more=T)
plot(p4d, split=c(2,2,4,3), more=T)
plot(p4, split = c(2,2,4,3), more=T)

plot(p5p, split=c(3,1,4,3), more=T)
plot(p5d, split=c(3,1,4,3), more=T)
plot(p5, split = c(3,1,4,3), more=T)

plot(p6p, split=c(3,2,4,3), more=T)
plot(p6d, split=c(3,2,4,3), more=T)
plot(p6, split = c(3,2,4,3), more=T)

plot(p7p, split=c(2,1,4,3), more=T)
plot(p7d, split=c(2,1,4,3), more=T)
plot(p7, split = c(2,1,4,3), more=T)

plot(p8p, split=c(1,3,4,3), more=T)
plot(p8d, split=c(1,3,4,3), more=T)
plot(p8, split = c(1,3,4,3), more=T)

plot(p9p, split=c(4,2,4,3), more=T)
plot(p9d, split=c(4,2,4,3), more=T)
plot(p9, split = c(4,2,4,3), more=T)

plot(p10p, split=c(2,3,4,3), more=T)
plot(p10d, split=c(2,3,4,3), more=T)
plot(p10, split = c(2,3,4,3), more=T)

plot(p6dz, split=c(3,3,4,3), more=T)
plot(p6z, split = c(3,3,4,3), more=T)

plot(p7dz, split=c(4,3,4,3), more=T)
plot(p7z, split = c(4,3,4,3), more=F)

dev.off()

png("Variograms_UniqueYAxis_CompareESDA_RmOps.png", width=13, height=10, units="in", res=300)
plot(p1p, split=c(1,1,4,3), more=T)
plot(p1dRmOps, split=c(1,1,4,3), more=T)
plot(p1RmOps, split = c(1,1,4,3), more=T)

plot(p2p, split=c(4,1,4,3), more=T)
plot(p2dRmOps, split=c(4,1,4,3), more=T)
plot(p2RmOps, split = c(4,1,4,3), more=T)

plot(p3p, split=c(1,2,4,3), more=T)
plot(p3dRmOps, split=c(1,2,4,3), more=T)
plot(p3RmOps, split = c(1,2,4,3), more=T)

plot(p4p, split=c(2,2,4,3), more=T)
plot(p4dRmOps, split=c(2,2,4,3), more=T)
plot(p4RmOps, split = c(2,2,4,3), more=T)

plot(p5p, split=c(3,1,4,3), more=T)
plot(p5dRmOps, split=c(3,1,4,3), more=T)
plot(p5RmOps, split = c(3,1,4,3), more=T)

plot(p6p, split=c(3,2,4,3), more=T)
plot(p6dRmOps, split=c(3,2,4,3), more=T)
plot(p6RmOps, split = c(3,2,4,3), more=T)

plot(p7p, split=c(2,1,4,3), more=T)
plot(p7dRmOps, split=c(2,1,4,3), more=T)
plot(p7RmOps, split = c(2,1,4,3), more=T)

plot(p8p, split=c(1,3,4,3), more=T)
plot(p8dRmOps, split=c(1,3,4,3), more=T)
plot(p8RmOps, split = c(1,3,4,3), more=T)

plot(p9p, split=c(4,2,4,3), more=T)
plot(p9dRmOps, split=c(4,2,4,3), more=T)
plot(p9RmOps, split = c(4,2,4,3), more=T)

plot(p10p, split=c(2,3,4,3), more=T)
plot(p10dRmOps, split=c(2,3,4,3), more=T)
plot(p10RmOps, split = c(2,3,4,3), more=T)

plot(p6dzRmOps, split=c(3,3,4,3), more=T)
plot(p6zRmOps, split = c(3,3,4,3), more=T)

plot(p7dzRmOps, split=c(4,3,4,3), more=T)
plot(p7zRmOps, split = c(4,3,4,3), more=F)

dev.off()


p2 = plot(v.CNY,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,200), col = 'red', xlim = c(0,60000), cex=0.5)
p2RmOps = plot(v.CNY_RmOps,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,200), col = 'red', xlim = c(0,60000), cex=0.5)
p2_k1p0 = plot(v.CNY_k1p0,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,200), col = kcols[1], xlim = c(0,60000), cex=0.5)
p2_k1p5 = plot(v.CNY_k1p5,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,200), col = kcols[2], xlim = c(0,60000), cex=0.5)
p2_k2p0 = plot(v.CNY_k2p0,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,200), col = kcols[3], xlim = c(0,60000), cex=0.5)
p2_k2p5 = plot(v.CNY_k2p5,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,200), col = kcols[4], xlim = c(0,60000), cex=0.5)
p2_k3p5 = plot(v.CNY_k3p5,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,200), col = kcols[5], xlim = c(0,60000), cex=0.5)
p2_k4p0 = plot(v.CNY_k4p0,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,200), col = kcols[6], xlim = c(0,60000), cex=0.5)
p2_k4p5 = plot(v.CNY_k4p5,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,200), col = kcols[7], xlim = c(0,60000), cex=0.5)
p2_k5p0 = plot(v.CNY_k5p0,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,200), col = kcols[8], xlim = c(0,60000), cex=0.5)
p2d = plot(v.DeepCNY,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,200), xlim = c(0,60000), cex=0.5)
p2dRmOps = plot(v.DeepCNY_RmOps,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, col = 'green', ylim = c(0,200), xlim = c(0,60000), cex=0.5)
p2dr = plot(rv.DeepCNY,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,200), xlim = c(0,60000), cex=0.5)
p2dRmOpsr = plot(rv.DeepCNY_RmOps,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, col = 'green', ylim = c(0,200), xlim = c(0,60000), cex=0.5)

p3 = plot(v.ENY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'red', ylim = c(0,100), xlim = c(0,60000), cex=0.5)
p3RmOps = plot(v.ENY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'red', ylim = c(0,100), xlim = c(0,60000), cex=0.5)
p3_k1p0 = plot(v.ENY_k1p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[1], ylim = c(0,100), xlim = c(0,60000), cex=0.5)
p3_k1p5 = plot(v.ENY_k1p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[2], ylim = c(0,100), xlim = c(0,60000), cex=0.5)
p3_k2p0 = plot(v.ENY_k2p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[3], ylim = c(0,100), xlim = c(0,60000), cex=0.5)
p3_k2p5 = plot(v.ENY_k2p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[4], ylim = c(0,100), xlim = c(0,60000), cex=0.5)
p3_k3p5 = plot(v.ENY_k3p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[5], ylim = c(0,100), xlim = c(0,60000), cex=0.5)
p3_k4p0 = plot(v.ENY_k4p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[6], ylim = c(0,100), xlim = c(0,60000), cex=0.5)
p3_k4p5 = plot(v.ENY_k4p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[7], ylim = c(0,100), xlim = c(0,60000), cex=0.5)
p3_k5p0 = plot(v.ENY_k5p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[8], ylim = c(0,100), xlim = c(0,60000), cex=0.5)
p3d = plot(v.DeepENY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,100), xlim = c(0,60000), cex=0.5)
p3dRmOps = plot(v.DeepENY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'green', ylim = c(0,100), xlim = c(0,60000), cex=0.5)
p3dr = plot(rv.DeepENY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,100), xlim = c(0,60000), cex=0.5)
p3dRmOpsr = plot(rv.DeepENY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'green', ylim = c(0,100), xlim = c(0,60000), cex=0.5)

p5 = plot(v.NWPANY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,150), cex=0.5)
p5RmOps = plot(v.NWPANY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,150), cex=0.5)
p5_k1p0 = plot(v.NWPANY_k1p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[1], xlim = c(0,60000), ylim = c(0,150), cex=0.5)
p5_k1p5 = plot(v.NWPANY_k1p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[2], xlim = c(0,60000), ylim = c(0,150), cex=0.5)
p5_k2p0 = plot(v.NWPANY_k2p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[3], xlim = c(0,60000), ylim = c(0,150), cex=0.5)
p5_k2p5 = plot(v.NWPANY_k2p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[4], xlim = c(0,60000), ylim = c(0,150), cex=0.5)
p5_k3p5 = plot(v.NWPANY_k3p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[5], xlim = c(0,60000), ylim = c(0,150), cex=0.5)
p5_k4p0 = plot(v.NWPANY_k4p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[6], xlim = c(0,60000), ylim = c(0,150), cex=0.5)
p5_k4p5 = plot(v.NWPANY_k4p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[7], xlim = c(0,60000), ylim = c(0,150), cex=0.5)
p5_k5p0 = plot(v.NWPANY_k5p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = kcols[8], xlim = c(0,60000), ylim = c(0,150), cex=0.5)
p5d = plot(v.DeepNWPANY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,150), cex=0.5)
p5dRmOps = plot(v.DeepNWPANY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,150), cex=0.5)
p5dr = plot(rv.DeepNWPANY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,150), cex=0.5)
p5dRmOpsr = plot(rv.DeepNWPANY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'green', xlim = c(0,60000), ylim = c(0,150), cex=0.5)

p9 = plot(v.CWV, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'red', ylim = c(0,400), xlim = c(0,60000), cex=0.5)
p9RmOps = plot(v.CWV_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'red', ylim = c(0,400), xlim = c(0,60000), cex=0.5)
p9_k1p0 = plot(v.CWV_k1p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[1], ylim = c(0,400), xlim = c(0,60000), cex=0.5)
p9_k1p5 = plot(v.CWV_k1p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[2], ylim = c(0,400), xlim = c(0,60000), cex=0.5)
p9_k2p0 = plot(v.CWV_k2p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[3], ylim = c(0,400), xlim = c(0,60000), cex=0.5)
p9_k2p5 = plot(v.CWV_k2p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[4], ylim = c(0,400), xlim = c(0,60000), cex=0.5)
p9_k3p5 = plot(v.CWV_k3p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[5], ylim = c(0,400), xlim = c(0,60000), cex=0.5)
p9_k4p0 = plot(v.CWV_k4p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[6], ylim = c(0,400), xlim = c(0,60000), cex=0.5)
p9_k4p5 = plot(v.CWV_k4p5, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[7], ylim = c(0,400), xlim = c(0,60000), cex=0.5)
p9_k5p0 = plot(v.CWV_k5p0, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = kcols[8], ylim = c(0,400), xlim = c(0,60000), cex=0.5)
p9d = plot(v.DeepCWV, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,400), xlim = c(0,60000), cex=0.5)
p9dRmOps = plot(v.DeepCWV_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'green', ylim = c(0,400), xlim = c(0,60000), cex=0.5)
p9dr = plot(rv.DeepCWV, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,400), xlim = c(0,60000), cex=0.5)
p9dRmOpsr = plot(rv.DeepCWV_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'green', ylim = c(0,400), xlim = c(0,60000), cex=0.5)


png("Variograms_UniqueYAxis_CompareESDA_DeepToOutliers.png", width=13, height=10, units="in", res=300)
plot(p1d, split=c(1,1,4,3), more=T)
plot(p1, split = c(1,1,4,3), more=T)

plot(p2d, split=c(4,1,4,3), more=T)
plot(p2, split = c(4,1,4,3), more=T)

plot(p3d, split=c(1,2,4,3), more=T)
plot(p3, split = c(1,2,4,3), more=T)

plot(p4d, split=c(2,2,4,3), more=T)
plot(p4, split = c(2,2,4,3), more=T)

plot(p5d, split=c(3,1,4,3), more=T)
plot(p5, split = c(3,1,4,3), more=T)

plot(p6dzz, split=c(3,2,4,3), more=T)
plot(p6zz, split = c(3,2,4,3), more=T)

plot(p7dzz, split=c(2,1,4,3), more=T)
plot(p7zz, split = c(2,1,4,3), more=T)

plot(p8d, split=c(1,3,4,3), more=T)
plot(p8, split = c(1,3,4,3), more=T)

plot(p9d, split=c(4,2,4,3), more=T)
plot(p9, split = c(4,2,4,3), more=T)

plot(p10dzz, split=c(2,3,4,3), more=T)
plot(p10zz, split = c(2,3,4,3), more=T)

dev.off()

png("Variograms_UniqueYAxis_CompareESDA_DeepToOutliers_RmOps.png", width=13, height=10, units="in", res=300)
plot(p1dRmOps, split=c(1,1,4,3), more=T)
plot(p1RmOps, split = c(1,1,4,3), more=T)

plot(p2dRmOps, split=c(4,1,4,3), more=T)
plot(p2RmOps, split = c(4,1,4,3), more=T)

plot(p3dRmOps, split=c(1,2,4,3), more=T)
plot(p3RmOps, split = c(1,2,4,3), more=T)

plot(p4dRmOps, split=c(2,2,4,3), more=T)
plot(p4RmOps, split = c(2,2,4,3), more=T)

plot(p5dRmOps, split=c(3,1,4,3), more=T)
plot(p5RmOps, split = c(3,1,4,3), more=T)

plot(p6dzzRmOps, split=c(3,2,4,3), more=T)
plot(p6zzRmOps, split = c(3,2,4,3), more=T)

plot(p7dzzRmOps, split=c(2,1,4,3), more=T)
plot(p7zzRmOps, split = c(2,1,4,3), more=T)

plot(p8dRmOps, split=c(1,3,4,3), more=T)
plot(p8RmOps, split = c(1,3,4,3), more=T)

plot(p9dRmOps, split=c(4,2,4,3), more=T)
plot(p9RmOps, split = c(4,2,4,3), more=T)

plot(p10dzzRmOps, split=c(2,3,4,3), more=T)
plot(p10zzRmOps, split = c(2,3,4,3), more=T)

dev.off()

png("Variograms_UniqueYAxis_CompareESDA_DeepToOutliers_Robust.png", width=13, height=10, units="in", res=300)
plot(p1d, split=c(1,1,4,3), more=T)
plot(p1dr, split=c(1,1,4,3), more=T)
plot(p1, split = c(1,1,4,3), more=T)

plot(p2d, split=c(4,1,4,3), more=T)
plot(p2dr, split=c(4,1,4,3), more=T)
plot(p2, split = c(4,1,4,3), more=T)

plot(p3d, split=c(1,2,4,3), more=T)
plot(p3dr, split=c(1,2,4,3), more=T)
plot(p3, split = c(1,2,4,3), more=T)

plot(p4d, split=c(2,2,4,3), more=T)
plot(p4dr, split=c(2,2,4,3), more=T)
plot(p4, split = c(2,2,4,3), more=T)

plot(p5d, split=c(3,1,4,3), more=T)
plot(p5dr, split=c(3,1,4,3), more=T)
plot(p5, split = c(3,1,4,3), more=T)

plot(p6dzz, split=c(3,2,4,3), more=T)
plot(p6dzzr, split=c(3,2,4,3), more=T)
plot(p6zz, split = c(3,2,4,3), more=T)

plot(p7dzz, split=c(2,1,4,3), more=T)
plot(p7dzzr, split=c(2,1,4,3), more=T)
plot(p7zz, split = c(2,1,4,3), more=T)

plot(p8d, split=c(1,3,4,3), more=T)
plot(p8dr, split=c(1,3,4,3), more=T)
plot(p8, split = c(1,3,4,3), more=T)

plot(p9d, split=c(4,2,4,3), more=T)
plot(p9dr, split=c(4,2,4,3), more=T)
plot(p9, split = c(4,2,4,3), more=T)

plot(p10dzz, split=c(2,3,4,3), more=T)
plot(p10dzzr, split=c(2,3,4,3), more=T)
plot(p10zz, split = c(2,3,4,3), more=T)

dev.off()

png("Variograms_UniqueYAxis_CompareESDA_DeepToOutliers_Robust_RmOps.png", width=13, height=10, units="in", res=300)
plot(p1dRmOps, split=c(1,1,4,3), more=T)
plot(p1dRmOpsr, split=c(1,1,4,3), more=T)
plot(p1RmOps, split = c(1,1,4,3), more=T)

plot(p2dRmOps, split=c(4,1,4,3), more=T)
plot(p2dRmOpsr, split=c(4,1,4,3), more=T)
plot(p2RmOps, split = c(4,1,4,3), more=T)

plot(p3dRmOps, split=c(1,2,4,3), more=T)
plot(p3dRmOpsr, split=c(1,2,4,3), more=T)
plot(p3RmOps, split = c(1,2,4,3), more=T)

plot(p4dRmOps, split=c(2,2,4,3), more=T)
plot(p4dRmOpsr, split=c(2,2,4,3), more=T)
plot(p4RmOps, split = c(2,2,4,3), more=T)

plot(p5dRmOps, split=c(3,1,4,3), more=T)
plot(p5dRmOpsr, split=c(3,1,4,3), more=T)
plot(p5RmOps, split = c(3,1,4,3), more=T)

plot(p6dzzRmOps, split=c(3,2,4,3), more=T)
plot(p6dzzRmOpsr, split=c(3,2,4,3), more=T)
plot(p6zzRmOps, split = c(3,2,4,3), more=T)

plot(p7dzzRmOps, split=c(2,1,4,3), more=T)
plot(p7dzzRmOpsr, split=c(2,1,4,3), more=T)
plot(p7zzRmOps, split = c(2,1,4,3), more=T)

plot(p8dRmOps, split=c(1,3,4,3), more=T)
plot(p8dRmOpsr, split=c(1,3,4,3), more=T)
plot(p8RmOps, split = c(1,3,4,3), more=T)

plot(p9dRmOps, split=c(4,2,4,3), more=T)
plot(p9dRmOpsr, split=c(4,2,4,3), more=T)
plot(p9RmOps, split = c(4,2,4,3), more=T)

plot(p10dzzRmOps, split=c(2,3,4,3), more=T)
plot(p10dzzRmOpsr, split=c(2,3,4,3), more=T)
plot(p10zzRmOps, split = c(2,3,4,3), more=T)

dev.off()

png("Variograms_UniqueYAxis_CompareESDA_DeepToOutliers_kAll.png", width=13, height=10, units="in", res=300)
plot(p1d, split=c(1,1,4,3), more=T)
plot(p1_k1p0, split = c(1,1,4,3), more=T)
plot(p1_k1p5, split = c(1,1,4,3), more=T)
plot(p1_k2p0, split = c(1,1,4,3), more=T)
plot(p1_k2p5, split = c(1,1,4,3), more=T)
plot(p1_k3p5, split = c(1,1,4,3), more=T)
plot(p1_k4p0, split = c(1,1,4,3), more=T)
plot(p1_k4p5, split = c(1,1,4,3), more=T)
plot(p1_k5p0, split = c(1,1,4,3), more=T)
plot(p1, split = c(1,1,4,3), more=T)

plot(p2d, split=c(4,1,4,3), more=T)
plot(p2_k1p0, split = c(4,1,4,3), more=T)
plot(p2_k1p5, split = c(4,1,4,3), more=T)
plot(p2_k2p0, split = c(4,1,4,3), more=T)
plot(p2_k2p5, split = c(4,1,4,3), more=T)
plot(p2_k3p5, split = c(4,1,4,3), more=T)
plot(p2_k4p0, split = c(4,1,4,3), more=T)
plot(p2_k4p5, split = c(4,1,4,3), more=T)
plot(p2_k5p0, split = c(4,1,4,3), more=T)
plot(p2, split = c(4,1,4,3), more=T)

plot(p3d, split=c(1,2,4,3), more=T)
plot(p3_k1p0, split = c(1,2,4,3), more=T)
plot(p3_k1p5, split = c(1,2,4,3), more=T)
plot(p3_k2p0, split = c(1,2,4,3), more=T)
plot(p3_k2p5, split = c(1,2,4,3), more=T)
plot(p3_k3p5, split = c(1,2,4,3), more=T)
plot(p3_k4p0, split = c(1,2,4,3), more=T)
plot(p3_k4p5, split = c(1,2,4,3), more=T)
plot(p3_k5p0, split = c(1,2,4,3), more=T)
plot(p3, split = c(1,2,4,3), more=T)

plot(p4d, split=c(2,2,4,3), more=T)
plot(p4_k1p0, split = c(2,2,4,3), more=T)
plot(p4_k1p5, split = c(2,2,4,3), more=T)
plot(p4_k2p0, split = c(2,2,4,3), more=T)
plot(p4_k2p5, split = c(2,2,4,3), more=T)
plot(p4_k3p5, split = c(2,2,4,3), more=T)
plot(p4_k4p0, split = c(2,2,4,3), more=T)
plot(p4_k4p5, split = c(2,2,4,3), more=T)
plot(p4_k5p0, split = c(2,2,4,3), more=T)
plot(p4, split = c(2,2,4,3), more=T)

plot(p5d, split=c(3,1,4,3), more=T)
plot(p5_k1p0, split = c(3,1,4,3), more=T)
plot(p5_k1p5, split = c(3,1,4,3), more=T)
plot(p5_k2p0, split = c(3,1,4,3), more=T)
plot(p5_k2p5, split = c(3,1,4,3), more=T)
plot(p5_k3p5, split = c(3,1,4,3), more=T)
plot(p5_k4p0, split = c(3,1,4,3), more=T)
plot(p5_k4p5, split = c(3,1,4,3), more=T)
plot(p5_k5p0, split = c(3,1,4,3), more=T)
plot(p5, split = c(3,1,4,3), more=T)

plot(p6dzz, split=c(3,2,4,3), more=T)
plot(p6z_k1p0, split = c(3,2,4,3), more=T)
plot(p6z_k1p5, split = c(3,2,4,3), more=T)
plot(p6z_k2p0, split = c(3,2,4,3), more=T)
plot(p6z_k2p5, split = c(3,2,4,3), more=T)
plot(p6z_k3p5, split = c(3,2,4,3), more=T)
plot(p6z_k4p0, split = c(3,2,4,3), more=T)
plot(p6z_k4p5, split = c(3,2,4,3), more=T)
plot(p6z_k5p0, split = c(3,2,4,3), more=T)
plot(p6zz, split = c(3,2,4,3), more=T)

plot(p7dzz, split=c(2,1,4,3), more=T)
plot(p7z_k1p0, split = c(2,1,4,3), more=T)
plot(p7z_k1p5, split = c(2,1,4,3), more=T)
plot(p7z_k2p0, split = c(2,1,4,3), more=T)
plot(p7z_k2p5, split = c(2,1,4,3), more=T)
plot(p7z_k3p5, split = c(2,1,4,3), more=T)
plot(p7z_k4p0, split = c(2,1,4,3), more=T)
plot(p7z_k4p5, split = c(2,1,4,3), more=T)
plot(p7z_k5p0, split = c(2,1,4,3), more=T)
plot(p7zz, split = c(2,1,4,3), more=T)

plot(p8d, split=c(1,3,4,3), more=T)
plot(p8_k1p0, split = c(1,3,4,3), more=T)
plot(p8_k1p5, split = c(1,3,4,3), more=T)
plot(p8_k2p0, split = c(1,3,4,3), more=T)
plot(p8_k2p5, split = c(1,3,4,3), more=T)
plot(p8_k3p5, split = c(1,3,4,3), more=T)
plot(p8_k4p0, split = c(1,3,4,3), more=T)
plot(p8_k4p5, split = c(1,3,4,3), more=T)
plot(p8_k5p0, split = c(1,3,4,3), more=T)
plot(p8, split = c(1,3,4,3), more=T)

plot(p9d, split=c(4,2,4,3), more=T)
plot(p9_k1p0, split = c(4,2,4,3), more=T)
plot(p9_k1p5, split = c(4,2,4,3), more=T)
plot(p9_k2p0, split = c(4,2,4,3), more=T)
plot(p9_k2p5, split = c(4,2,4,3), more=T)
plot(p9_k3p5, split = c(4,2,4,3), more=T)
plot(p9_k4p0, split = c(4,2,4,3), more=T)
plot(p9_k4p5, split = c(4,2,4,3), more=T)
plot(p9_k5p0, split = c(4,2,4,3), more=T)
plot(p9, split = c(4,2,4,3), more=T)

plot(p10dzz, split=c(2,3,4,3), more=T)
plot(p10z_k1p0, split = c(2,3,4,3), more=T)
plot(p10z_k1p5, split = c(2,3,4,3), more=T)
plot(p10z_k2p0, split = c(2,3,4,3), more=T)
plot(p10z_k2p5, split = c(2,3,4,3), more=T)
plot(p10z_k3p5, split = c(2,3,4,3), more=T)
plot(p10z_k4p0, split = c(2,3,4,3), more=T)
plot(p10z_k4p5, split = c(2,3,4,3), more=T)
plot(p10z_k5p0, split = c(2,3,4,3), more=T)
plot(p10zz, split = c(2,3,4,3), more=T)

dev.off()

#   Confidence intervals for variogram lags----
#Note that if the number of point pairs in a bin is 2, then the standard deviation will be 0 because when only 2 samples are removed in that bin there is sample reuse because each point contributes the same value when it is removed. Therefore the jackknife variance is likely too small when the number of data points is small.
#Also note that point pairs is not the same as number of points.
# JackKnife = function(Dat,      #Spatial points datafame containing the points to be jackknifed. Should be in UTM coordinates.
#                       bins,    #Number of lags for the semi-variogram
#                       cut,     #Cutoff distance for the semi-variogram
#                       anis=NA, #anisotropy for semi-variogram. Currently only works for 2 directions.
#                       v.Dat    #variogram model for Dat using all of the data. Must have same bins, cut, and anis as specified.
#                       )
#   {  
#   if (is.na(anis[1]) == FALSE){
#     #Anisotropic Variogram
#     if (length(anis) < 2){
#       print('Error, need to have at least 2 angles for the anisotropy angle')
#       stop
#     }
#     #Expand the length of VarioMat matrix to accommodate the number of angles in anis
#     VarioMat = Pts = Dist = Gamma = matrix(NA, nrow=nrow(Dat)*length(anis), ncol=bins)
#     
#     #Make a vector of the number of replications of each bin from the jackknife. This number only changes if a given resampling results in no estimate for a given bin.
#     Nmat = rbind(rep(nrow(Dat), bins), rep(nrow(Dat), bins))
#     
#     #Fixme: Make a storage for the anisotropic variogram for any angle length. Currently at 2 (most common) but can be nrow(Dat)*length(anis) + i
#     for (i in 1:nrow(Dat)){
#       #Make a variogram, and save the bin information in VarioMat
#       #Use anisotropic variogram.
#       vj = variogram(Qs~1, data = Dat[-i,], cutoff=cut, width=cut/bins, alpha=anis)
#       
#       #Get the indicies in vj that are populated with information for the direction of interest.
#       Ind1 = which(vj$dir.hor == anis[1])
#       Ind2 = which(vj$dir.hor == anis[2])
#       
#       Pts[i, 1:length(Ind1)] = vj$np[Ind1]
#       Pts[(nrow(Dat) + i), 1:length(Ind2)] = vj$np[Ind2]
#       Dist[i,1:length(Ind1)] = vj$dist[Ind1]
#       Dist[(nrow(Dat) + i), 1:length(Ind2)] = vj$dist[Ind2]
#       Gamma[i, 1:length(Ind1)] = vj$gamma[Ind1]
#       Gamma[(nrow(Dat) + i), 1:length(Ind2)] = vj$gamma[Ind2]
#       
#       if (any(is.na(Pts[i,])) || any(Pts[i,] < 2)){
#         #Mark the distance as NA and reduce the number of points in Nvec. The model will have to be rerun with the new Nvec...
#         ind = which(is.na(Pts[i,]) == TRUE | Pts[i,] < 2)
#         Pts[i, ind] = NA
#         Dist[i, ind] = NA
#         Nmat[1, ind] = Nmat[1, ind] - 1
#         Gamma[i, ind] = NA
#       }
#       if (any(is.na(Pts[(nrow(Dat) + i),])) || any(Pts[(nrow(Dat) + i),] < 2)){
#         #Mark the distance as NA and reduce the number of points in Nvec. The model will have to be rerun with the new Nvec...
#         ind = which(is.na(Pts[(nrow(Dat) + i),]) == TRUE | Pts[(nrow(Dat) + i),] < 2)
#         Pts[(nrow(Dat) + i), ind] = NA
#         Dist[(nrow(Dat) + i), ind] = NA
#         Nmat[2, ind] = Nmat[2, ind] - 1
#         Gamma[(nrow(Dat) + i), ind] = NA
#       }
#       if (i == nrow(Dat)){
#         print('finished first loop')
#       }
#     }
#     rm(i)
#     
#     #Get the indices of the anisotropy for the variogram with all of the points. This will contain the maximum number of bins for each angle.
#     Ind1 = which(v.Dat$dir.hor == anis[1])
#     Ind2 = which(v.Dat$dir.hor == anis[2])
#     
#     #Fill in the matrix of the estimates. The transposition is a faster way of using sweep() to multiply a vector by row of a matrix.
#     VarioMat[1:nrow(Dat), 1:length(Ind1)] = t(Nmat[1, 1:length(Ind1)]*v.Dat$gamma[Ind1] - t(t(t(Gamma[1:nrow(Dat), 1:length(Ind1)])*(Nmat[1, 1:length(Ind1)]-1))))
#     VarioMat[(nrow(Dat)+1):(2*nrow(Dat)), 1:length(Ind2)] = t(Nmat[2, 1:length(Ind2)]*v.Dat$gamma[Ind2] - t(t(t(Gamma[(nrow(Dat)+1):(2*nrow(Dat)), 1:length(Ind2)])*(Nmat[2, 1:length(Ind2)]-1))))
#     
#     #Calculate the Jackknife mean
#     BinMean_1 = apply(VarioMat[1:nrow(Dat), ], 2, FUN=sum, na.rm=TRUE)/Nmat[1,]
#     BinMean_2 = apply(VarioMat[(nrow(Dat)+1):(nrow(Dat)*2), ], 2, FUN=sum, na.rm=TRUE)/Nmat[2,]
#     
#     BinMean = rbind(BinMean_1, BinMean_2)
#     
#     #Calculate the average bin distances for plotting purposes
#     BinMean_dist1 = apply(Dist[1:nrow(Dat), ], 2, FUN=sum, na.rm=TRUE)/Nmat[1,]
#     BinMean_dist2 = apply(Dist[(nrow(Dat)+1):(2*nrow(Dat)), ], 2, FUN=sum, na.rm=TRUE)/Nmat[2,]
#     
#     BinMean_dist = rbind(BinMean_dist1, BinMean_dist2)
#     
#     
#     #Calculate the jackknife standard error
#     VarioMatSquaredMatrix = matrix(NA, ncol=bins, nrow=nrow(Dat)*length(anis))
#     DistSd = matrix(NA, ncol=bins, nrow=nrow(Dat)*length(anis))
#     for (j in 1:(nrow(Dat)*length(anis))){
#       if (j <= nrow(Dat)){
#         VarioMatSquaredMatrix[j,] = (VarioMat[j,] - BinMean_1)^2
#         DistSd[j, ] = (Dist[j, ] - BinMean_dist1)^2
#       }
#       else{
#         VarioMatSquaredMatrix[j,] = (VarioMat[j,] - BinMean_2)^2
#         DistSd[j, ] = (Dist[j, ] - BinMean_dist2)^2
#       }
#       if (j == nrow(Dat)*length(anis)){
#         print('finished second loop')
#       }
#     }
#     rm(j)
#     
#     
#     #Calculate the jackknife variance for each bin
#     VarEst = matrix(NA, ncol=bins, nrow=length(anis))
#     VarEst[1,] = apply(VarioMatSquaredMatrix[1:nrow(Dat),], 2, FUN=sum, na.rm=TRUE)/(Nmat[1,]*(Nmat[1,]-1))
#     VarEst[2,] = apply(VarioMatSquaredMatrix[(nrow(Dat)+1):(nrow(Dat)*2),], 2, FUN=sum, na.rm=TRUE)/(Nmat[2,]*(Nmat[2,]-1))
#     
#     SdEst = sqrt(VarEst)
#     
#     #Calculate the standard deviation of the bin distances for plotting error bars on the positions
#     BinVar_dist1 = apply(DistSd[1:nrow(Dat), ], 2, FUN=sum, na.rm=TRUE)/(Nmat[1,] - 1)
#     BinVar_dist2 = apply(DistSd[(nrow(Dat)+1):(2*nrow(Dat)), ], 2, FUN=sum, na.rm=TRUE)/(Nmat[2,] - 1)
#     BinSd_dist1 = sqrt(BinVar_dist1)
#     BinSd_dist2 = sqrt(BinVar_dist2)
#     BinSd_dist = rbind(BinSd_dist1, BinSd_dist2)
#     
#     Nvec=Nmat
#   }
#   else {
#     #Make a matrix for storing the bin estimates of variance, number of points, and the distance to points.
#     VarioMat = Pts = Dist = Gamma = matrix(NA, nrow=nrow(Dat), ncol=bins)
#     
#     #Make a vector of the number of replications of each bin from the jackknife. This number only changes if a given resampling results in no estimate for a given bin.
#     Nvec = rep(nrow(Dat), bins)
#     
#     #Start jackknife
#     for (i in 1:nrow(Dat)){
#       vj = variogram(Qs~1, data = Dat[-i,], cutoff=cut, width=cut/bins)
#       
#       #Store the number of points in each lag and the average lag distance.
#       Pts[i,] = vj$np
#       Dist[i,] = vj$dist
#       Gamma[i,] = vj$gamma
#       
#       if (any(is.na(Pts[i,])) || any(Pts[i,] < 2)){
#         #Mark the distance as NA and reduce the number of points in Nvec. The model will have to be rerun with the new Nvec...
#         ind = which((is.na(Pts[i,]) == TRUE) | (Pts[i,] < 2))
#         Pts[i, ind] = NA
#         Dist[i, ind] = NA
#         Nvec[ind] = Nvec[ind] - 1
#         Gamma[i, ind] = NA
#       }
#     }
#     rm(i)
#     
#     #Calculate the estimate of the bin mean from the jackknife. The transposition is a faster way of using sweep() to multiply a vector by row of a matrix.
#     VarioMat =  t(Nvec*v.Dat$gamma - t(t(t(Gamma)*(Nvec-1))))
#     
#     #Calculate the Jackknife mean. Remove NAs generated from before.
#     BinMean = apply(VarioMat, 2, FUN=sum, na.rm=TRUE)/Nvec
#     
#     #Calculate the average bin distance for plotting purposes
#     BinMean_dist = apply(Dist, 2, FUN=sum, na.rm=TRUE)/Nvec
#     
#     #Calculate the jackknife standard error
#     VarioMatSquaredMatrix = matrix(NA, ncol=bins, nrow=nrow(Dat))
#     DistSd = matrix(NA, ncol=bins, nrow=nrow(Dat))
#     for (j in 1:nrow(Dat)){
#       VarioMatSquaredMatrix[j, ] = (VarioMat[j, ] - BinMean)^2
#       DistSd[j, ] = (Dist[j, ] - BinMean_dist)^2
#     }
#     rm(j)
#     
#     #Calculate the jackknife variance for each bin
#     VarEst = apply(VarioMatSquaredMatrix, 2, FUN=sum, na.rm=TRUE)/(Nvec*(Nvec-1))
#     
#     SdEst = sqrt(VarEst)
#     
#     #Calculate the standard deviation of the bin distances for plotting error bars on the positions
#     BinVar_dist = apply(DistSd, 2, FUN=sum, na.rm=TRUE)/(Nvec - 1)
#     BinSd_dist = sqrt(BinVar_dist)
#     
#   }
#   
#   #return a list
#   lst = list(SdEst = SdEst, BinMean = BinMean, AvgDist = BinMean_dist, SdDist = BinSd_dist, NumPts = Pts, N = Nvec)
#   return(lst)
# }

#    Parallel Jackknife function----
#Fixme: Anisotropy component still needs work, see previous commented out function
JackKnifePar = function(Dat,      #Spatial points datafame containing the points to be jackknifed. Should be in UTM coordinates.
                     bins,    #Number of lags for the semi-variogram
                     cut,     #Cutoff distance for the semi-variogram
                     anis=NA, #anisotropy for semi-variogram. Currently only works for 2 directions.
                     v.Dat    #variogram model for Dat using all of the data. Must have same bins, cut, and anis as specified.
)
{  
  if (is.na(anis[1]) == FALSE){
    #Anisotropic Variogram
    if (length(anis) < 2){
      print('Error, need to have at least 2 angles for the anisotropy angle')
      stop
    }
    #Expand the length of VarioMat matrix to accommodate the number of angles in anis
    VarioMat = Pts = Dist = Gamma = matrix(NA, nrow=nrow(Dat)*length(anis), ncol=bins)
    
    #Make a vector of the number of replications of each bin from the jackknife. This number only changes if a given resampling results in no estimate for a given bin.
    Nmat = rbind(rep(nrow(Dat), bins), rep(nrow(Dat), bins))
    
    #Fixme: Make a storage for the anisotropic variogram for any angle length. Currently at 2 (most common) but can be nrow(Dat)*length(anis) + i
    for (i in 1:nrow(Dat)){
      #Make a variogram, and save the bin information in VarioMat
      #Use anisotropic variogram.
      vj = variogram(Qs~1, data = Dat[-i,], cutoff=cut, width=cut/bins, alpha=anis)
      
      #Get the indicies in vj that are populated with information for the direction of interest.
      Ind1 = which(vj$dir.hor == anis[1])
      Ind2 = which(vj$dir.hor == anis[2])
      
      Pts[i, 1:length(Ind1)] = vj$np[Ind1]
      Pts[(nrow(Dat) + i), 1:length(Ind2)] = vj$np[Ind2]
      Dist[i,1:length(Ind1)] = vj$dist[Ind1]
      Dist[(nrow(Dat) + i), 1:length(Ind2)] = vj$dist[Ind2]
      Gamma[i, 1:length(Ind1)] = vj$gamma[Ind1]
      Gamma[(nrow(Dat) + i), 1:length(Ind2)] = vj$gamma[Ind2]
      
      if (any(is.na(Pts[i,])) || any(Pts[i,] < 2)){
        #Mark the distance as NA and reduce the number of points in Nvec. The model will have to be rerun with the new Nvec...
        ind = which(is.na(Pts[i,]) == TRUE | Pts[i,] < 2)
        Pts[i, ind] = NA
        Dist[i, ind] = NA
        Nmat[1, ind] = Nmat[1, ind] - 1
        Gamma[i, ind] = NA
      }
      if (any(is.na(Pts[(nrow(Dat) + i),])) || any(Pts[(nrow(Dat) + i),] < 2)){
        #Mark the distance as NA and reduce the number of points in Nvec. The model will have to be rerun with the new Nvec...
        ind = which(is.na(Pts[(nrow(Dat) + i),]) == TRUE | Pts[(nrow(Dat) + i),] < 2)
        Pts[(nrow(Dat) + i), ind] = NA
        Dist[(nrow(Dat) + i), ind] = NA
        Nmat[2, ind] = Nmat[2, ind] - 1
        Gamma[(nrow(Dat) + i), ind] = NA
      }
      if (i == nrow(Dat)){
        print('finished first loop')
      }
    }
    rm(i)
    
    #Get the indices of the anisotropy for the variogram with all of the points. This will contain the maximum number of bins for each angle.
    Ind1 = which(v.Dat$dir.hor == anis[1])
    Ind2 = which(v.Dat$dir.hor == anis[2])
    
    #Fill in the matrix of the estimates. The transposition is a faster way of using sweep() to multiply a vector by row of a matrix.
    VarioMat[1:nrow(Dat), 1:length(Ind1)] = t(Nmat[1, 1:length(Ind1)]*v.Dat$gamma[Ind1] - t(t(t(Gamma[1:nrow(Dat), 1:length(Ind1)])*(Nmat[1, 1:length(Ind1)]-1))))
    VarioMat[(nrow(Dat)+1):(2*nrow(Dat)), 1:length(Ind2)] = t(Nmat[2, 1:length(Ind2)]*v.Dat$gamma[Ind2] - t(t(t(Gamma[(nrow(Dat)+1):(2*nrow(Dat)), 1:length(Ind2)])*(Nmat[2, 1:length(Ind2)]-1))))
    
    #Calculate the Jackknife mean
    BinMean_1 = apply(VarioMat[1:nrow(Dat), ], 2, FUN=sum, na.rm=TRUE)/Nmat[1,]
    BinMean_2 = apply(VarioMat[(nrow(Dat)+1):(nrow(Dat)*2), ], 2, FUN=sum, na.rm=TRUE)/Nmat[2,]
    
    BinMean = rbind(BinMean_1, BinMean_2)
    
    #Calculate the average bin distances for plotting purposes
    BinMean_dist1 = apply(Dist[1:nrow(Dat), ], 2, FUN=sum, na.rm=TRUE)/Nmat[1,]
    BinMean_dist2 = apply(Dist[(nrow(Dat)+1):(2*nrow(Dat)), ], 2, FUN=sum, na.rm=TRUE)/Nmat[2,]
    
    BinMean_dist = rbind(BinMean_dist1, BinMean_dist2)
    
    
    #Calculate the jackknife standard error
    VarioMatSquaredMatrix = matrix(NA, ncol=bins, nrow=nrow(Dat)*length(anis))
    DistSd = matrix(NA, ncol=bins, nrow=nrow(Dat)*length(anis))
    for (j in 1:(nrow(Dat)*length(anis))){
      if (j <= nrow(Dat)){
        VarioMatSquaredMatrix[j,] = (VarioMat[j,] - BinMean_1)^2
        DistSd[j, ] = (Dist[j, ] - BinMean_dist1)^2
      }
      else{
        VarioMatSquaredMatrix[j,] = (VarioMat[j,] - BinMean_2)^2
        DistSd[j, ] = (Dist[j, ] - BinMean_dist2)^2
      }
      if (j == nrow(Dat)*length(anis)){
        print('finished second loop')
      }
    }
    rm(j)
    
    
    #Calculate the jackknife variance for each bin
    VarEst = matrix(NA, ncol=bins, nrow=length(anis))
    VarEst[1,] = apply(VarioMatSquaredMatrix[1:nrow(Dat),], 2, FUN=sum, na.rm=TRUE)/(Nmat[1,]*(Nmat[1,]-1))
    VarEst[2,] = apply(VarioMatSquaredMatrix[(nrow(Dat)+1):(nrow(Dat)*2),], 2, FUN=sum, na.rm=TRUE)/(Nmat[2,]*(Nmat[2,]-1))
    
    SdEst = sqrt(VarEst)
    
    #Calculate the standard deviation of the bin distances for plotting error bars on the positions
    BinVar_dist1 = apply(DistSd[1:nrow(Dat), ], 2, FUN=sum, na.rm=TRUE)/(Nmat[1,] - 1)
    BinVar_dist2 = apply(DistSd[(nrow(Dat)+1):(2*nrow(Dat)), ], 2, FUN=sum, na.rm=TRUE)/(Nmat[2,] - 1)
    BinSd_dist1 = sqrt(BinVar_dist1)
    BinSd_dist2 = sqrt(BinVar_dist2)
    BinSd_dist = rbind(BinSd_dist1, BinSd_dist2)
    
    Nvec=Nmat
  }
  else {
    #Make a vector of the number of replications of each bin from the jackknife. This number only changes if a given resampling results in no estimate for a given bin.
    Nvec = rep(nrow(Dat), bins)
    
    #Start jackknife
    Store = foreach(i = 1:nrow(Dat), .packages = "gstat", .combine = "rbind") %dopar% {
      vj = variogram(Qs~1, data = Dat[-i,], cutoff=cut, width=cut/bins)
      
      #Store the number of points in each lag and the average lag distance.
      Ptsi = vj$np
      Disti = vj$dist
      Gammai = vj$gamma
      N = Nvec
      
      if (any(is.na(Ptsi)) || any(Ptsi < 2)){
        #Mark the distance as NA and reduce the number of points in Nvec. The model will have to be rerun with the new Nvec...
        ind = which((is.na(Ptsi) == TRUE) | (Ptsi < 2))
        Ptsi[ind] = NA
        Disti[ind] = NA
        N[ind] = N[ind] - 1
        Gammai[ind] = NA
      }
      lst = data.frame(Pts = Ptsi, Dist = Disti, Gamma = Gammai, N = N)
      lst
    }
    
    #Get data from Store into matrices
    Pts = t(matrix(Store$Pts, ncol = nrow(Dat), nrow = bins))
    Dist = t(matrix(Store$Dist, ncol = nrow(Dat), nrow = bins))
    Gamma = t(matrix(Store$Gamma, ncol = nrow(Dat), nrow = bins))
    Nvec = Nvec - apply((nrow(Dat) - t(matrix(Store$N, ncol = nrow(Dat), nrow = bins))), MARGIN = 2, FUN = sum)
    
    #Calculate the estimate of the bin mean from the jackknife. The transposition is a faster way of using sweep() to multiply a vector by row of a matrix.
    VarioMat =  t(Nvec*v.Dat$gamma - t(t(t(Gamma)*(Nvec-1))))
    
    #Calculate the Jackknife mean. Remove NAs generated from before.
    BinMean = apply(VarioMat, 2, FUN=sum, na.rm=TRUE)/Nvec
    
    #Calculate the average bin distance for plotting purposes
    BinMean_dist = apply(Dist, 2, FUN=sum, na.rm=TRUE)/Nvec
    
    #Calculate the jackknife standard error
    Store2 = foreach (j = 1:nrow(Dat), .combine = "rbind") %dopar% {
      VarioMatSquaredMatrix = (VarioMat[j, ] - BinMean)^2
      DistSd = (Dist[j, ] - BinMean_dist)^2
      
      lst = data.frame(Varios = VarioMatSquaredMatrix, Dists = DistSd)
      lst
    }
    
    VarioMatSquaredMatrix = t(matrix(Store2$Varios, ncol = nrow(Dat), nrow = bins))
    DistSd = t(matrix(Store2$Dists, ncol = nrow(Dat), nrow = bins))
    
    #Calculate the jackknife variance for each bin
    VarEst = apply(VarioMatSquaredMatrix, 2, FUN=sum, na.rm=TRUE)/(Nvec*(Nvec-1))
    
    SdEst = sqrt(VarEst)
    
    #Calculate the standard deviation of the bin distances for plotting error bars on the positions
    BinVar_dist = apply(DistSd, 2, FUN=sum, na.rm=TRUE)/(Nvec - 1)
    BinSd_dist = sqrt(BinVar_dist)
  }
  
  #return a list
  lst = list(SdEst = SdEst, BinMean = BinMean, AvgDist = BinMean_dist, SdDist = BinSd_dist, NumPts = Pts, N = Nvec)
  return(lst)
}

#Detect number of computer cores:
cores = detectCores() - 1 
cl = makeCluster(cores)
registerDoParallel(cl)
#    Compute jackknife estimates of variogram----
JackPreCT  = JackKnifePar(Dat = PreCT, bins=50, cut=60000, v.Dat = v.PreCT)
JackPreCNY = JackKnifePar(Dat = PreCNY, bins=15, cut=60000, v.Dat = v.PreCNY)
JackPreCWV = JackKnifePar(Dat = PreCWV, bins=50, cut=60000, v.Dat = v.PreCWV)
JackPreENY = JackKnifePar(Dat = PreENY, bins=15, cut=60000, v.Dat = v.PreENY)
JackPreENYPA = JackKnifePar(Dat = PreENYPA, bins=40, cut=60000, v.Dat = v.PreENYPA)
JackPreWPA  = JackKnifePar(Dat = PreWPA, bins=50, cut=60000, v.Dat = v.PreWPA)
JackPreSWPA = JackKnifePar(Dat = PreSWPA, bins=50, cut=60000, v.Dat = v.PreSWPA)
JackPreNWPANY = JackKnifePar(Dat = PreNWPANY, bins=20, cut=60000, v.Dat = v.PreNWPANY)
JackPreMT = JackKnifePar(Dat = PreMT, bins=50, cut=60000, v.Dat = v.PreMT)
JackPreVR = JackKnifePar(Dat = PreVR, bins=20, cut=60000, v.Dat = v.PreVR)

JackDeepCT  = JackKnifePar(Dat = DeepCT, bins=50, cut=60000, v.Dat = v.DeepCT)
JackDeepCNY = JackKnifePar(Dat = DeepCNY, bins=15, cut=60000, v.Dat = v.DeepCNY)
JackDeepCWV = JackKnifePar(Dat = DeepCWV, bins=50, cut=60000, v.Dat = v.DeepCWV)
JackDeepENY = JackKnifePar(Dat = DeepENY, bins=15, cut=60000, v.Dat = v.DeepENY)
JackDeepENYPA = JackKnifePar(Dat = DeepENYPA, bins=40, cut=60000, v.Dat = v.DeepENYPA)
JackDeepWPA  = JackKnifePar(Dat = DeepWPA, bins=50, cut=60000, v.Dat = v.DeepWPA)
JackDeepSWPA = JackKnifePar(Dat = DeepSWPA, bins=50, cut=60000, v.Dat = v.DeepSWPA)
JackDeepNWPANY = JackKnifePar(Dat = DeepNWPANY, bins=20, cut=60000, v.Dat = v.DeepNWPANY)
JackDeepMT = JackKnifePar(Dat = DeepMT, bins=50, cut=60000, v.Dat = v.DeepMT)
JackDeepVR = JackKnifePar(Dat = DeepVR, bins=20, cut=60000, v.Dat = v.DeepVR)

JackCT  = JackKnifePar(Dat = CT, bins=50, cut=60000, v.Dat = v.CT)
JackCNY = JackKnifePar(Dat = CNY, bins=15, cut=60000, v.Dat = v.CNY)
JackCWV = JackKnifePar(Dat = CWV, bins=50, cut=60000, v.Dat = v.CWV)
JackENY = JackKnifePar(Dat = ENY, bins=15, cut=60000, v.Dat = v.ENY)
JackENYPA = JackKnifePar(Dat = ENYPA, bins=40, cut=60000, v.Dat = v.ENYPA)
JackWPA  = JackKnifePar(Dat = WPA, bins=50, cut=60000, v.Dat = v.WPA)
JackSWPA = JackKnifePar(Dat = SWPA, bins=50, cut=60000, v.Dat = v.SWPA)
JackNWPANY = JackKnifePar(Dat = NWPANY, bins=20, cut=60000, v.Dat = v.NWPANY)
JackMT = JackKnifePar(Dat = MT, bins=50, cut=60000, v.Dat = v.MT)
JackVR = JackKnifePar(Dat = VR, bins=20, cut=60000, v.Dat = v.VR)

JackCT_k1p0  = JackKnifePar(Dat = CT_k1p0, bins=50, cut=60000, v.Dat = v.CT_k1p0)
JackCNY_k1p0 = JackKnifePar(Dat = CNY_k1p0, bins=15, cut=60000, v.Dat = v.CNY_k1p0)
JackCWV_k1p0 = JackKnifePar(Dat = CWV_k1p0, bins=50, cut=60000, v.Dat = v.CWV_k1p0)
JackENY_k1p0 = JackKnifePar(Dat = ENY_k1p0, bins=15, cut=60000, v.Dat = v.ENY_k1p0)
JackENYPA_k1p0 = JackKnifePar(Dat = ENYPA_k1p0, bins=40, cut=60000, v.Dat = v.ENYPA_k1p0)
JackWPA_k1p0  = JackKnifePar(Dat = WPA_k1p0, bins=50, cut=60000, v.Dat = v.WPA_k1p0)
JackSWPA_k1p0 = JackKnifePar(Dat = SWPA_k1p0, bins=50, cut=60000, v.Dat = v.SWPA_k1p0)
JackNWPANY_k1p0 = JackKnifePar(Dat = NWPANY_k1p0, bins=20, cut=60000, v.Dat = v.NWPANY_k1p0)
JackMT_k1p0 = JackKnifePar(Dat = MT_k1p0, bins=50, cut=60000, v.Dat = v.MT_k1p0)
JackVR_k1p0 = JackKnifePar(Dat = VR_k1p0, bins=20, cut=60000, v.Dat = v.VR_k1p0)

JackCT_k1p5  = JackKnifePar(Dat = CT_k1p5, bins=50, cut=60000, v.Dat = v.CT_k1p5)
JackCNY_k1p5 = JackKnifePar(Dat = CNY_k1p5, bins=15, cut=60000, v.Dat = v.CNY_k1p5)
JackCWV_k1p5 = JackKnifePar(Dat = CWV_k1p5, bins=50, cut=60000, v.Dat = v.CWV_k1p5)
JackENY_k1p5 = JackKnifePar(Dat = ENY_k1p5, bins=15, cut=60000, v.Dat = v.ENY_k1p5)
JackENYPA_k1p5 = JackKnifePar(Dat = ENYPA_k1p5, bins=40, cut=60000, v.Dat = v.ENYPA_k1p5)
JackWPA_k1p5  = JackKnifePar(Dat = WPA_k1p5, bins=50, cut=60000, v.Dat = v.WPA_k1p5)
JackSWPA_k1p5 = JackKnifePar(Dat = SWPA_k1p5, bins=50, cut=60000, v.Dat = v.SWPA_k1p5)
JackNWPANY_k1p5 = JackKnifePar(Dat = NWPANY_k1p5, bins=20, cut=60000, v.Dat = v.NWPANY_k1p5)
JackMT_k1p5 = JackKnifePar(Dat = MT_k1p5, bins=50, cut=60000, v.Dat = v.MT_k1p5)
JackVR_k1p5 = JackKnifePar(Dat = VR_k1p5, bins=20, cut=60000, v.Dat = v.VR_k1p5)

JackCT_k2p0  = JackKnifePar(Dat = CT_k2p0, bins=50, cut=60000, v.Dat = v.CT_k2p0)
JackCNY_k2p0 = JackKnifePar(Dat = CNY_k2p0, bins=15, cut=60000, v.Dat = v.CNY_k2p0)
JackCWV_k2p0 = JackKnifePar(Dat = CWV_k2p0, bins=50, cut=60000, v.Dat = v.CWV_k2p0)
JackENY_k2p0 = JackKnifePar(Dat = ENY_k2p0, bins=15, cut=60000, v.Dat = v.ENY_k2p0)
JackENYPA_k2p0 = JackKnifePar(Dat = ENYPA_k2p0, bins=40, cut=60000, v.Dat = v.ENYPA_k2p0)
JackWPA_k2p0  = JackKnifePar(Dat = WPA_k2p0, bins=50, cut=60000, v.Dat = v.WPA_k2p0)
JackSWPA_k2p0 = JackKnifePar(Dat = SWPA_k2p0, bins=50, cut=60000, v.Dat = v.SWPA_k2p0)
JackNWPANY_k2p0 = JackKnifePar(Dat = NWPANY_k2p0, bins=20, cut=60000, v.Dat = v.NWPANY_k2p0)
JackMT_k2p0 = JackKnifePar(Dat = MT_k2p0, bins=50, cut=60000, v.Dat = v.MT_k2p0)
JackVR_k2p0 = JackKnifePar(Dat = VR_k2p0, bins=20, cut=60000, v.Dat = v.VR_k2p0)

JackCT_k2p5  = JackKnifePar(Dat = CT_k2p5, bins=50, cut=60000, v.Dat = v.CT_k2p5)
JackCNY_k2p5 = JackKnifePar(Dat = CNY_k2p5, bins=15, cut=60000, v.Dat = v.CNY_k2p5)
JackCWV_k2p5 = JackKnifePar(Dat = CWV_k2p5, bins=50, cut=60000, v.Dat = v.CWV_k2p5)
JackENY_k2p5 = JackKnifePar(Dat = ENY_k2p5, bins=15, cut=60000, v.Dat = v.ENY_k2p5)
JackENYPA_k2p5 = JackKnifePar(Dat = ENYPA_k2p5, bins=40, cut=60000, v.Dat = v.ENYPA_k2p5)
JackWPA_k2p5  = JackKnifePar(Dat = WPA_k2p5, bins=50, cut=60000, v.Dat = v.WPA_k2p5)
JackSWPA_k2p5 = JackKnifePar(Dat = SWPA_k2p5, bins=50, cut=60000, v.Dat = v.SWPA_k2p5)
JackNWPANY_k2p5 = JackKnifePar(Dat = NWPANY_k2p5, bins=20, cut=60000, v.Dat = v.NWPANY_k2p5)
JackMT_k2p5 = JackKnifePar(Dat = MT_k2p5, bins=50, cut=60000, v.Dat = v.MT_k2p5)
JackVR_k2p5 = JackKnifePar(Dat = VR_k2p5, bins=20, cut=60000, v.Dat = v.VR_k2p5)

JackCT_k3p5  = JackKnifePar(Dat = CT_k3p5, bins=50, cut=60000, v.Dat = v.CT_k3p5)
JackCNY_k3p5 = JackKnifePar(Dat = CNY_k3p5, bins=15, cut=60000, v.Dat = v.CNY_k3p5)
JackCWV_k3p5 = JackKnifePar(Dat = CWV_k3p5, bins=50, cut=60000, v.Dat = v.CWV_k3p5)
JackENY_k3p5 = JackKnifePar(Dat = ENY_k3p5, bins=15, cut=60000, v.Dat = v.ENY_k3p5)
JackENYPA_k3p5 = JackKnifePar(Dat = ENYPA_k3p5, bins=40, cut=60000, v.Dat = v.ENYPA_k3p5)
JackWPA_k3p5  = JackKnifePar(Dat = WPA_k3p5, bins=50, cut=60000, v.Dat = v.WPA_k3p5)
JackSWPA_k3p5 = JackKnifePar(Dat = SWPA_k3p5, bins=50, cut=60000, v.Dat = v.SWPA_k3p5)
JackNWPANY_k3p5 = JackKnifePar(Dat = NWPANY_k3p5, bins=20, cut=60000, v.Dat = v.NWPANY_k3p5)
JackMT_k3p5 = JackKnifePar(Dat = MT_k3p5, bins=50, cut=60000, v.Dat = v.MT_k3p5)
JackVR_k3p5 = JackKnifePar(Dat = VR_k3p5, bins=20, cut=60000, v.Dat = v.VR_k3p5)

JackCT_k4p0  = JackKnifePar(Dat = CT_k4p0, bins=50, cut=60000, v.Dat = v.CT_k4p0)
JackCNY_k4p0 = JackKnifePar(Dat = CNY_k4p0, bins=15, cut=60000, v.Dat = v.CNY_k4p0)
JackCWV_k4p0 = JackKnifePar(Dat = CWV_k4p0, bins=50, cut=60000, v.Dat = v.CWV_k4p0)
JackENY_k4p0 = JackKnifePar(Dat = ENY_k4p0, bins=15, cut=60000, v.Dat = v.ENY_k4p0)
JackENYPA_k4p0 = JackKnifePar(Dat = ENYPA_k4p0, bins=40, cut=60000, v.Dat = v.ENYPA_k4p0)
JackWPA_k4p0  = JackKnifePar(Dat = WPA_k4p0, bins=50, cut=60000, v.Dat = v.WPA_k4p0)
JackSWPA_k4p0 = JackKnifePar(Dat = SWPA_k4p0, bins=50, cut=60000, v.Dat = v.SWPA_k4p0)
JackNWPANY_k4p0 = JackKnifePar(Dat = NWPANY_k4p0, bins=20, cut=60000, v.Dat = v.NWPANY_k4p0)
JackMT_k4p0 = JackKnifePar(Dat = MT_k4p0, bins=50, cut=60000, v.Dat = v.MT_k4p0)
JackVR_k4p0 = JackKnifePar(Dat = VR_k4p0, bins=20, cut=60000, v.Dat = v.VR_k4p0)

JackCT_k4p5  = JackKnifePar(Dat = CT_k4p5, bins=50, cut=60000, v.Dat = v.CT_k4p5)
JackCNY_k4p5 = JackKnifePar(Dat = CNY_k4p5, bins=15, cut=60000, v.Dat = v.CNY_k4p5)
JackCWV_k4p5 = JackKnifePar(Dat = CWV_k4p5, bins=50, cut=60000, v.Dat = v.CWV_k4p5)
JackENY_k4p5 = JackKnifePar(Dat = ENY_k4p5, bins=15, cut=60000, v.Dat = v.ENY_k4p5)
JackENYPA_k4p5 = JackKnifePar(Dat = ENYPA_k4p5, bins=40, cut=60000, v.Dat = v.ENYPA_k4p5)
JackWPA_k4p5  = JackKnifePar(Dat = WPA_k4p5, bins=50, cut=60000, v.Dat = v.WPA_k4p5)
JackSWPA_k4p5 = JackKnifePar(Dat = SWPA_k4p5, bins=50, cut=60000, v.Dat = v.SWPA_k4p5)
JackNWPANY_k4p5 = JackKnifePar(Dat = NWPANY_k4p5, bins=20, cut=60000, v.Dat = v.NWPANY_k4p5)
JackMT_k4p5 = JackKnifePar(Dat = MT_k4p5, bins=50, cut=60000, v.Dat = v.MT_k4p5)
JackVR_k4p5 = JackKnifePar(Dat = VR_k4p5, bins=20, cut=60000, v.Dat = v.VR_k4p5)

JackCT_k5p0  = JackKnifePar(Dat = CT_k5p0, bins=50, cut=60000, v.Dat = v.CT_k5p0)
JackCNY_k5p0 = JackKnifePar(Dat = CNY_k5p0, bins=15, cut=60000, v.Dat = v.CNY_k5p0)
JackCWV_k5p0 = JackKnifePar(Dat = CWV_k5p0, bins=50, cut=60000, v.Dat = v.CWV_k5p0)
JackENY_k5p0 = JackKnifePar(Dat = ENY_k5p0, bins=15, cut=60000, v.Dat = v.ENY_k5p0)
JackENYPA_k5p0 = JackKnifePar(Dat = ENYPA_k5p0, bins=40, cut=60000, v.Dat = v.ENYPA_k5p0)
JackWPA_k5p0  = JackKnifePar(Dat = WPA_k5p0, bins=50, cut=60000, v.Dat = v.WPA_k5p0)
JackSWPA_k5p0 = JackKnifePar(Dat = SWPA_k5p0, bins=50, cut=60000, v.Dat = v.SWPA_k5p0)
JackNWPANY_k5p0 = JackKnifePar(Dat = NWPANY_k5p0, bins=20, cut=60000, v.Dat = v.NWPANY_k5p0)
JackMT_k5p0 = JackKnifePar(Dat = MT_k5p0, bins=50, cut=60000, v.Dat = v.MT_k5p0)
JackVR_k5p0 = JackKnifePar(Dat = VR_k5p0, bins=20, cut=60000, v.Dat = v.VR_k5p0)

JackDeepCT_RmOps  = JackKnifePar(Dat = DeepCT_RmOps, bins=50, cut=60000, v.Dat = v.DeepCT_RmOps)
JackDeepCNY_RmOps = JackKnifePar(Dat = DeepCNY_RmOps, bins=15, cut=60000, v.Dat = v.DeepCNY_RmOps)
JackDeepCWV_RmOps = JackKnifePar(Dat = DeepCWV_RmOps, bins=50, cut=60000, v.Dat = v.DeepCWV_RmOps)
JackDeepENY_RmOps = JackKnifePar(Dat = DeepENY_RmOps, bins=15, cut=60000, v.Dat = v.DeepENY_RmOps)
JackDeepENYPA_RmOps = JackKnifePar(Dat = DeepENYPA_RmOps, bins=40, cut=60000, v.Dat = v.DeepENYPA_RmOps)
JackDeepWPA_RmOps  = JackKnifePar(Dat = DeepWPA_RmOps, bins=50, cut=60000, v.Dat = v.DeepWPA_RmOps)
JackDeepSWPA_RmOps = JackKnifePar(Dat = DeepSWPA_RmOps, bins=50, cut=60000, v.Dat = v.DeepSWPA_RmOps)
JackDeepNWPANY_RmOps = JackKnifePar(Dat = DeepNWPANY_RmOps, bins=20, cut=60000, v.Dat = v.DeepNWPANY_RmOps)
JackDeepMT_RmOps = JackKnifePar(Dat = DeepMT_RmOps, bins=50, cut=60000, v.Dat = v.DeepMT_RmOps)
JackDeepVR_RmOps = JackKnifePar(Dat = DeepVR_RmOps, bins=20, cut=60000, v.Dat = v.DeepVR_RmOps)

JackCT_RmOps  = JackKnifePar(Dat = CT_RmOps, bins=50, cut=60000, v.Dat = v.CT_RmOps)
JackCNY_RmOps = JackKnifePar(Dat = CNY_RmOps, bins=15, cut=60000, v.Dat = v.CNY_RmOps)
JackCWV_RmOps = JackKnifePar(Dat = CWV_RmOps, bins=50, cut=60000, v.Dat = v.CWV_RmOps)
JackENY_RmOps = JackKnifePar(Dat = ENY_RmOps, bins=15, cut=60000, v.Dat = v.ENY_RmOps)
JackENYPA_RmOps = JackKnifePar(Dat = ENYPA_RmOps, bins=40, cut=60000, v.Dat = v.ENYPA_RmOps)
JackWPA_RmOps  = JackKnifePar(Dat = WPA_RmOps, bins=50, cut=60000, v.Dat = v.WPA_RmOps)
JackSWPA_RmOps = JackKnifePar(Dat = SWPA_RmOps, bins=50, cut=60000, v.Dat = v.SWPA_RmOps)
JackNWPANY_RmOps = JackKnifePar(Dat = NWPANY_RmOps, bins=20, cut=60000, v.Dat = v.NWPANY_RmOps)
JackMT_RmOps = JackKnifePar(Dat = MT_RmOps, bins=50, cut=60000, v.Dat = v.MT_RmOps)
JackVR_RmOps = JackKnifePar(Dat = VR_RmOps, bins=20, cut=60000, v.Dat = v.VR_RmOps)

stopCluster(cl)

#    Plot the variograms with confidence intervals----
png("Variograms_UniqueYAxis_CompareESDA_95ConfInt.png", width=13, height=10, units="in", res=300)
plot(p1p, split=c(1,1,4,3), more=T)
plot(p1d, split=c(1,1,4,3), more=T)
plot(p1, split = c(1,1,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackCT$AvgDist, y = (JackCT$BinMean + JackCT$SdEst*2), ylim=p1$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackCT$AvgDist, y = (JackCT$BinMean - JackCT$SdEst*2), ylim=p1$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreCT$AvgDist, y = (JackPreCT$BinMean + JackPreCT$SdEst*2), ylim=p1$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreCT$AvgDist, y = (JackPreCT$BinMean - JackPreCT$SdEst*2), ylim=p1$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepCT$AvgDist, y = (JackDeepCT$BinMean + JackDeepCT$SdEst*2), ylim=p1$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepCT$AvgDist, y = (JackDeepCT$BinMean - JackDeepCT$SdEst*2), ylim=p1$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p2p, split=c(4,1,4,3), more=T)
plot(p2d, split=c(4,1,4,3), more=T)
plot(p2, split = c(4,1,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackCNY$AvgDist, y = (JackCNY$BinMean + JackCNY$SdEst*2), ylim=p2$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackCNY$AvgDist, y = (JackCNY$BinMean - JackCNY$SdEst*2), ylim=p2$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreCNY$AvgDist, y = (JackPreCNY$BinMean + JackPreCNY$SdEst*2), ylim=p2$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreCNY$AvgDist, y = (JackPreCNY$BinMean - JackPreCNY$SdEst*2), ylim=p2$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepCNY$AvgDist, y = (JackDeepCNY$BinMean + JackDeepCNY$SdEst*2), ylim=p2$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepCNY$AvgDist, y = (JackDeepCNY$BinMean - JackDeepCNY$SdEst*2), ylim=p2$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p3p, split=c(1,2,4,3), more=T)
plot(p3d, split=c(1,2,4,3), more=T)
plot(p3, split = c(1,2,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackENY$AvgDist, y = (JackENY$BinMean + JackENY$SdEst*2), ylim=p3$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackENY$AvgDist, y = (JackENY$BinMean - JackENY$SdEst*2), ylim=p3$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreENY$AvgDist, y = (JackPreENY$BinMean + JackPreENY$SdEst*2), ylim=p3$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreENY$AvgDist, y = (JackPreENY$BinMean - JackPreENY$SdEst*2), ylim=p3$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepENY$AvgDist, y = (JackDeepENY$BinMean + JackDeepENY$SdEst*2), ylim=p3$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepENY$AvgDist, y = (JackDeepENY$BinMean - JackDeepENY$SdEst*2), ylim=p3$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p4p, split=c(2,2,4,3), more=T)
plot(p4d, split=c(2,2,4,3), more=T)
plot(p4, split = c(2,2,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackENYPA$AvgDist, y = (JackENYPA$BinMean + JackENYPA$SdEst*2), ylim=p4$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackENYPA$AvgDist, y = (JackENYPA$BinMean - JackENYPA$SdEst*2), ylim=p4$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreENYPA$AvgDist, y = (JackPreENYPA$BinMean + JackPreENYPA$SdEst*2), ylim=p4$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreENYPA$AvgDist, y = (JackPreENYPA$BinMean - JackPreENYPA$SdEst*2), ylim=p4$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepENYPA$AvgDist, y = (JackDeepENYPA$BinMean + JackDeepENYPA$SdEst*2), ylim=p4$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepENYPA$AvgDist, y = (JackDeepENYPA$BinMean - JackDeepENYPA$SdEst*2), ylim=p4$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p5p, split=c(3,1,4,3), more=T)
plot(p5d, split=c(3,1,4,3), more=T)
plot(p5, split = c(3,1,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackNWPANY$AvgDist, y = (JackNWPANY$BinMean + JackNWPANY$SdEst*2), ylim=p5$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackNWPANY$AvgDist, y = (JackNWPANY$BinMean - JackNWPANY$SdEst*2), ylim=p5$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreNWPANY$AvgDist, y = (JackPreNWPANY$BinMean + JackPreNWPANY$SdEst*2), ylim=p5$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreNWPANY$AvgDist, y = (JackPreNWPANY$BinMean - JackPreNWPANY$SdEst*2), ylim=p5$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepNWPANY$AvgDist, y = (JackDeepNWPANY$BinMean + JackDeepNWPANY$SdEst*2), ylim=p5$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepNWPANY$AvgDist, y = (JackDeepNWPANY$BinMean - JackDeepNWPANY$SdEst*2), ylim=p5$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p6p, split=c(3,2,4,3), more=T)
plot(p6d, split=c(3,2,4,3), more=T)
plot(p6, split = c(3,2,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackSWPA$AvgDist, y = (JackSWPA$BinMean + JackSWPA$SdEst*2), ylim=p6$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackSWPA$AvgDist, y = (JackSWPA$BinMean - JackSWPA$SdEst*2), ylim=p6$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreSWPA$AvgDist, y = (JackPreSWPA$BinMean + JackPreSWPA$SdEst*2), ylim=p6$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreSWPA$AvgDist, y = (JackPreSWPA$BinMean - JackPreSWPA$SdEst*2), ylim=p6$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepSWPA$AvgDist, y = (JackDeepSWPA$BinMean + JackDeepSWPA$SdEst*2), ylim=p6$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepSWPA$AvgDist, y = (JackDeepSWPA$BinMean - JackDeepSWPA$SdEst*2), ylim=p6$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p7p, split=c(2,1,4,3), more=T)
plot(p7d, split=c(2,1,4,3), more=T)
plot(p7, split = c(2,1,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackWPA$AvgDist, y = (JackWPA$BinMean + JackWPA$SdEst*2), ylim=p7$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackWPA$AvgDist, y = (JackWPA$BinMean - JackWPA$SdEst*2), ylim=p7$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreWPA$AvgDist, y = (JackPreWPA$BinMean + JackPreWPA$SdEst*2), ylim=p7$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreWPA$AvgDist, y = (JackPreWPA$BinMean - JackPreWPA$SdEst*2), ylim=p7$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepWPA$AvgDist, y = (JackDeepWPA$BinMean + JackDeepWPA$SdEst*2), ylim=p7$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepWPA$AvgDist, y = (JackDeepWPA$BinMean - JackDeepWPA$SdEst*2), ylim=p7$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p8p, split=c(1,3,4,3), more=T)
plot(p8d, split=c(1,3,4,3), more=T)
plot(p8, split = c(1,3,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackMT$AvgDist, y = (JackMT$BinMean + JackMT$SdEst*2), ylim=p8$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackMT$AvgDist, y = (JackMT$BinMean - JackMT$SdEst*2), ylim=p8$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreMT$AvgDist, y = (JackPreMT$BinMean + JackPreMT$SdEst*2), ylim=p8$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreMT$AvgDist, y = (JackPreMT$BinMean - JackPreMT$SdEst*2), ylim=p8$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepMT$AvgDist, y = (JackDeepMT$BinMean + JackDeepMT$SdEst*2), ylim=p8$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepMT$AvgDist, y = (JackDeepMT$BinMean - JackDeepMT$SdEst*2), ylim=p8$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p9p, split=c(4,2,4,3), more=T)
plot(p9d, split=c(4,2,4,3), more=T)
plot(p9, split = c(4,2,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackCWV$AvgDist, y = (JackCWV$BinMean + JackCWV$SdEst*2), ylim=p9$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackCWV$AvgDist, y = (JackCWV$BinMean - JackCWV$SdEst*2), ylim=p9$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreCWV$AvgDist, y = (JackPreCWV$BinMean + JackPreCWV$SdEst*2), ylim=p9$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreCWV$AvgDist, y = (JackPreCWV$BinMean - JackPreCWV$SdEst*2), ylim=p9$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepCWV$AvgDist, y = (JackDeepCWV$BinMean + JackDeepCWV$SdEst*2), ylim=p9$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepCWV$AvgDist, y = (JackDeepCWV$BinMean - JackDeepCWV$SdEst*2), ylim=p9$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p10p, split=c(2,3,4,3), more=T)
plot(p10d, split=c(2,3,4,3), more=T)
plot(p10, split = c(2,3,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackVR$AvgDist, y = (JackVR$BinMean + JackVR$SdEst*2), ylim=p10$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackVR$AvgDist, y = (JackVR$BinMean - JackVR$SdEst*2), ylim=p10$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreVR$AvgDist, y = (JackPreVR$BinMean + JackPreVR$SdEst*2), ylim=p10$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreVR$AvgDist, y = (JackPreVR$BinMean - JackPreVR$SdEst*2), ylim=p10$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepVR$AvgDist, y = (JackDeepVR$BinMean + JackDeepVR$SdEst*2), ylim=p10$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepVR$AvgDist, y = (JackDeepVR$BinMean - JackDeepVR$SdEst*2), ylim=p10$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

dev.off()

png("Variograms_UniqueYAxis_CompareESDA_95ConfInt_RmOps.png", width=13, height=10, units="in", res=300)
plot(p1p, split=c(1,1,4,3), more=T)
plot(p1dRmOps, split=c(1,1,4,3), more=T)
plot(p1RmOps, split = c(1,1,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackCT_RmOps$AvgDist, y = (JackCT_RmOps$BinMean + JackCT_RmOps$SdEst*2), ylim=p1RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackCT_RmOps$AvgDist, y = (JackCT_RmOps$BinMean - JackCT_RmOps$SdEst*2), ylim=p1RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreCT$AvgDist, y = (JackPreCT$BinMean + JackPreCT$SdEst*2), ylim=p1RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreCT$AvgDist, y = (JackPreCT$BinMean - JackPreCT$SdEst*2), ylim=p1RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepCT_RmOps$AvgDist, y = (JackDeepCT_RmOps$BinMean + JackDeepCT_RmOps$SdEst*2), ylim=p1RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepCT_RmOps$AvgDist, y = (JackDeepCT_RmOps$BinMean - JackDeepCT_RmOps$SdEst*2), ylim=p1RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p2p, split=c(4,1,4,3), more=T)
plot(p2dRmOps, split=c(4,1,4,3), more=T)
plot(p2RmOps, split = c(4,1,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackCNY_RmOps$AvgDist, y = (JackCNY_RmOps$BinMean + JackCNY_RmOps$SdEst*2), ylim=p2RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackCNY_RmOps$AvgDist, y = (JackCNY_RmOps$BinMean - JackCNY_RmOps$SdEst*2), ylim=p2RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreCNY$AvgDist, y = (JackPreCNY$BinMean + JackPreCNY$SdEst*2), ylim=p2RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreCNY$AvgDist, y = (JackPreCNY$BinMean - JackPreCNY$SdEst*2), ylim=p2RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepCNY_RmOps$AvgDist, y = (JackDeepCNY_RmOps$BinMean + JackDeepCNY_RmOps$SdEst*2), ylim=p2RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepCNY_RmOps$AvgDist, y = (JackDeepCNY_RmOps$BinMean - JackDeepCNY_RmOps$SdEst*2), ylim=p2RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p3p, split=c(1,2,4,3), more=T)
plot(p3dRmOps, split=c(1,2,4,3), more=T)
plot(p3RmOps, split = c(1,2,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackENY_RmOps$AvgDist, y = (JackENY_RmOps$BinMean + JackENY_RmOps$SdEst*2), ylim=p3RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackENY_RmOps$AvgDist, y = (JackENY_RmOps$BinMean - JackENY_RmOps$SdEst*2), ylim=p3RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreENY$AvgDist, y = (JackPreENY$BinMean + JackPreENY$SdEst*2), ylim=p3RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreENY$AvgDist, y = (JackPreENY$BinMean - JackPreENY$SdEst*2), ylim=p3RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepENY_RmOps$AvgDist, y = (JackDeepENY_RmOps$BinMean + JackDeepENY_RmOps$SdEst*2), ylim=p3RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepENY_RmOps$AvgDist, y = (JackDeepENY_RmOps$BinMean - JackDeepENY_RmOps$SdEst*2), ylim=p3RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p4p, split=c(2,2,4,3), more=T)
plot(p4dRmOps, split=c(2,2,4,3), more=T)
plot(p4RmOps, split = c(2,2,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackENYPA_RmOps$AvgDist, y = (JackENYPA_RmOps$BinMean + JackENYPA_RmOps$SdEst*2), ylim=p4RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackENYPA_RmOps$AvgDist, y = (JackENYPA_RmOps$BinMean - JackENYPA_RmOps$SdEst*2), ylim=p4RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreENYPA$AvgDist, y = (JackPreENYPA$BinMean + JackPreENYPA$SdEst*2), ylim=p4RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreENYPA$AvgDist, y = (JackPreENYPA$BinMean - JackPreENYPA$SdEst*2), ylim=p4RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepENYPA_RmOps$AvgDist, y = (JackDeepENYPA_RmOps$BinMean + JackDeepENYPA_RmOps$SdEst*2), ylim=p4RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepENYPA_RmOps$AvgDist, y = (JackDeepENYPA_RmOps$BinMean - JackDeepENYPA_RmOps$SdEst*2), ylim=p4RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p5p, split=c(3,1,4,3), more=T)
plot(p5dRmOps, split=c(3,1,4,3), more=T)
plot(p5RmOps, split = c(3,1,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackNWPANY_RmOps$AvgDist, y = (JackNWPANY_RmOps$BinMean + JackNWPANY_RmOps$SdEst*2), ylim=p5RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackNWPANY_RmOps$AvgDist, y = (JackNWPANY_RmOps$BinMean - JackNWPANY_RmOps$SdEst*2), ylim=p5RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreNWPANY$AvgDist, y = (JackPreNWPANY$BinMean + JackPreNWPANY$SdEst*2), ylim=p5RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreNWPANY$AvgDist, y = (JackPreNWPANY$BinMean - JackPreNWPANY$SdEst*2), ylim=p5RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepNWPANY_RmOps$AvgDist, y = (JackDeepNWPANY_RmOps$BinMean + JackDeepNWPANY_RmOps$SdEst*2), ylim=p5RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepNWPANY_RmOps$AvgDist, y = (JackDeepNWPANY_RmOps$BinMean - JackDeepNWPANY_RmOps$SdEst*2), ylim=p5RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p6p, split=c(3,2,4,3), more=T)
plot(p6dRmOps, split=c(3,2,4,3), more=T)
plot(p6RmOps, split = c(3,2,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackSWPA_RmOps$AvgDist, y = (JackSWPA_RmOps$BinMean + JackSWPA_RmOps$SdEst*2), ylim=p6RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackSWPA_RmOps$AvgDist, y = (JackSWPA_RmOps$BinMean - JackSWPA_RmOps$SdEst*2), ylim=p6RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreSWPA$AvgDist, y = (JackPreSWPA$BinMean + JackPreSWPA$SdEst*2), ylim=p6RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreSWPA$AvgDist, y = (JackPreSWPA$BinMean - JackPreSWPA$SdEst*2), ylim=p6RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepSWPA_RmOps$AvgDist, y = (JackDeepSWPA_RmOps$BinMean + JackDeepSWPA_RmOps$SdEst*2), ylim=p6RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepSWPA_RmOps$AvgDist, y = (JackDeepSWPA_RmOps$BinMean - JackDeepSWPA_RmOps$SdEst*2), ylim=p6RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p7p, split=c(2,1,4,3), more=T)
plot(p7dRmOps, split=c(2,1,4,3), more=T)
plot(p7RmOps, split = c(2,1,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackWPA_RmOps$AvgDist, y = (JackWPA_RmOps$BinMean + JackWPA_RmOps$SdEst*2), ylim=p7RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackWPA_RmOps$AvgDist, y = (JackWPA_RmOps$BinMean - JackWPA_RmOps$SdEst*2), ylim=p7RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreWPA$AvgDist, y = (JackPreWPA$BinMean + JackPreWPA$SdEst*2), ylim=p7RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreWPA$AvgDist, y = (JackPreWPA$BinMean - JackPreWPA$SdEst*2), ylim=p7RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepWPA_RmOps$AvgDist, y = (JackDeepWPA_RmOps$BinMean + JackDeepWPA_RmOps$SdEst*2), ylim=p7RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepWPA_RmOps$AvgDist, y = (JackDeepWPA_RmOps$BinMean - JackDeepWPA_RmOps$SdEst*2), ylim=p7RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p8p, split=c(1,3,4,3), more=T)
plot(p8dRmOps, split=c(1,3,4,3), more=T)
plot(p8RmOps, split = c(1,3,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackMT_RmOps$AvgDist, y = (JackMT_RmOps$BinMean + JackMT_RmOps$SdEst*2), ylim=p8RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackMT_RmOps$AvgDist, y = (JackMT_RmOps$BinMean - JackMT_RmOps$SdEst*2), ylim=p8RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreMT$AvgDist, y = (JackPreMT$BinMean + JackPreMT$SdEst*2), ylim=p8RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreMT$AvgDist, y = (JackPreMT$BinMean - JackPreMT$SdEst*2), ylim=p8RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepMT_RmOps$AvgDist, y = (JackDeepMT_RmOps$BinMean + JackDeepMT_RmOps$SdEst*2), ylim=p8RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepMT_RmOps$AvgDist, y = (JackDeepMT_RmOps$BinMean - JackDeepMT_RmOps$SdEst*2), ylim=p8RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p9p, split=c(4,2,4,3), more=T)
plot(p9dRmOps, split=c(4,2,4,3), more=T)
plot(p9RmOps, split = c(4,2,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackCWV_RmOps$AvgDist, y = (JackCWV_RmOps$BinMean + JackCWV_RmOps$SdEst*2), ylim=p9RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackCWV_RmOps$AvgDist, y = (JackCWV_RmOps$BinMean - JackCWV_RmOps$SdEst*2), ylim=p9RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreCWV$AvgDist, y = (JackPreCWV$BinMean + JackPreCWV$SdEst*2), ylim=p9RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreCWV$AvgDist, y = (JackPreCWV$BinMean - JackPreCWV$SdEst*2), ylim=p9RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepCWV_RmOps$AvgDist, y = (JackDeepCWV_RmOps$BinMean + JackDeepCWV_RmOps$SdEst*2), ylim=p9RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepCWV_RmOps$AvgDist, y = (JackDeepCWV_RmOps$BinMean - JackDeepCWV_RmOps$SdEst*2), ylim=p9RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p10p, split=c(2,3,4,3), more=T)
plot(p10dRmOps, split=c(2,3,4,3), more=T)
plot(p10RmOps, split = c(2,3,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackVR_RmOps$AvgDist, y = (JackVR_RmOps$BinMean + JackVR_RmOps$SdEst*2), ylim=p10RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackVR_RmOps$AvgDist, y = (JackVR_RmOps$BinMean - JackVR_RmOps$SdEst*2), ylim=p10RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackPreVR$AvgDist, y = (JackPreVR$BinMean + JackPreVR$SdEst*2), ylim=p10RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackPreVR$AvgDist, y = (JackPreVR$BinMean - JackPreVR$SdEst*2), ylim=p10RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='black')
lpoints(x = JackDeepVR_RmOps$AvgDist, y = (JackDeepVR_RmOps$BinMean + JackDeepVR_RmOps$SdEst*2), ylim=p10RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepVR_RmOps$AvgDist, y = (JackDeepVR_RmOps$BinMean - JackDeepVR_RmOps$SdEst*2), ylim=p10RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

dev.off()

#Only deep and post-ESDA
p2 = plot(v.CNY,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,200), col = 'red', xlim = c(0,60000), cex=0.5)
p2d = plot(v.DeepCNY,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,200), xlim = c(0,60000), cex=0.5)

p2RmOps = plot(v.CNY_RmOps,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,200), col = 'red', xlim = c(0,60000), cex=0.5)
p2dRmOps = plot(v.DeepCNY_RmOps,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,200), xlim = c(0,60000), cex=0.5)

p3 = plot(v.ENY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'red', ylim = c(0,100), xlim = c(0,60000), cex=0.5)
p3d = plot(v.DeepENY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,100), xlim = c(0,60000), cex=0.5)

p3RmOps = plot(v.ENY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'red', ylim = c(0,100), xlim = c(0,60000), cex=0.5)
p3dRmOps = plot(v.DeepENY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,100), xlim = c(0,60000), cex=0.5)

p5 = plot(v.NWPANY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,150), cex=0.5)
p5d = plot(v.DeepNWPANY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,150), cex=0.5)

p5RmOps = plot(v.NWPANY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,150), cex=0.5)
p5dRmOps = plot(v.DeepNWPANY_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,150), cex=0.5)

p9 = plot(v.CWV, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'red', ylim = c(0,400), xlim = c(0,60000), cex=0.5)
p9d = plot(v.DeepCWV, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,400), xlim = c(0,60000), cex=0.5)

p9RmOps = plot(v.CWV_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'red', ylim = c(0,400), xlim = c(0,60000), cex=0.5)
p9dRmOps = plot(v.DeepCWV_RmOps, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', ylim = c(0,400), xlim = c(0,60000), cex=0.5)

png("Variograms_UniqueYAxis_CompareESDA_DeepPost_95ConfInt.png", width=13, height=10, units="in", res=300)
plot(p1d, split=c(1,1,4,3), more=T)
plot(p1, split = c(1,1,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackCT$AvgDist, y = (JackCT$BinMean + JackCT$SdEst*2), ylim=p1$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackCT$AvgDist, y = (JackCT$BinMean - JackCT$SdEst*2), ylim=p1$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepCT$AvgDist, y = (JackDeepCT$BinMean + JackDeepCT$SdEst*2), ylim=p1$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepCT$AvgDist, y = (JackDeepCT$BinMean - JackDeepCT$SdEst*2), ylim=p1$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p2d, split=c(4,1,4,3), more=T)
plot(p2, split = c(4,1,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackCNY$AvgDist, y = (JackCNY$BinMean + JackCNY$SdEst*2), ylim=p2$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackCNY$AvgDist, y = (JackCNY$BinMean - JackCNY$SdEst*2), ylim=p2$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepCNY$AvgDist, y = (JackDeepCNY$BinMean + JackDeepCNY$SdEst*2), ylim=p2$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepCNY$AvgDist, y = (JackDeepCNY$BinMean - JackDeepCNY$SdEst*2), ylim=p2$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p3d, split=c(1,2,4,3), more=T)
plot(p3, split = c(1,2,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackENY$AvgDist, y = (JackENY$BinMean + JackENY$SdEst*2), ylim=p3$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackENY$AvgDist, y = (JackENY$BinMean - JackENY$SdEst*2), ylim=p3$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepENY$AvgDist, y = (JackDeepENY$BinMean + JackDeepENY$SdEst*2), ylim=p3$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepENY$AvgDist, y = (JackDeepENY$BinMean - JackDeepENY$SdEst*2), ylim=p3$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p4d, split=c(2,2,4,3), more=T)
plot(p4, split = c(2,2,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackENYPA$AvgDist, y = (JackENYPA$BinMean + JackENYPA$SdEst*2), ylim=p4$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackENYPA$AvgDist, y = (JackENYPA$BinMean - JackENYPA$SdEst*2), ylim=p4$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepENYPA$AvgDist, y = (JackDeepENYPA$BinMean + JackDeepENYPA$SdEst*2), ylim=p4$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepENYPA$AvgDist, y = (JackDeepENYPA$BinMean - JackDeepENYPA$SdEst*2), ylim=p4$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p5d, split=c(3,1,4,3), more=T)
plot(p5, split = c(3,1,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackNWPANY$AvgDist, y = (JackNWPANY$BinMean + JackNWPANY$SdEst*2), ylim=p5$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackNWPANY$AvgDist, y = (JackNWPANY$BinMean - JackNWPANY$SdEst*2), ylim=p5$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepNWPANY$AvgDist, y = (JackDeepNWPANY$BinMean + JackDeepNWPANY$SdEst*2), ylim=p5$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepNWPANY$AvgDist, y = (JackDeepNWPANY$BinMean - JackDeepNWPANY$SdEst*2), ylim=p5$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p6dzz, split=c(3,2,4,3), more=T)
plot(p6zz, split = c(3,2,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackSWPA$AvgDist, y = (JackSWPA$BinMean + JackSWPA$SdEst*2), ylim=p6zz$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackSWPA$AvgDist, y = (JackSWPA$BinMean - JackSWPA$SdEst*2), ylim=p6zz$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepSWPA$AvgDist, y = (JackDeepSWPA$BinMean + JackDeepSWPA$SdEst*2), ylim=p6zz$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepSWPA$AvgDist, y = (JackDeepSWPA$BinMean - JackDeepSWPA$SdEst*2), ylim=p6zz$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p7dzz, split=c(2,1,4,3), more=T)
plot(p7zz, split = c(2,1,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackWPA$AvgDist, y = (JackWPA$BinMean + JackWPA$SdEst*2), ylim=p7zz$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackWPA$AvgDist, y = (JackWPA$BinMean - JackWPA$SdEst*2), ylim=p7zz$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepWPA$AvgDist, y = (JackDeepWPA$BinMean + JackDeepWPA$SdEst*2), ylim=p7zz$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepWPA$AvgDist, y = (JackDeepWPA$BinMean - JackDeepWPA$SdEst*2), ylim=p7zz$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p8d, split=c(1,3,4,3), more=T)
plot(p8, split = c(1,3,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackMT$AvgDist, y = (JackMT$BinMean + JackMT$SdEst*2), ylim=p8$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackMT$AvgDist, y = (JackMT$BinMean - JackMT$SdEst*2), ylim=p8$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepMT$AvgDist, y = (JackDeepMT$BinMean + JackDeepMT$SdEst*2), ylim=p8$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepMT$AvgDist, y = (JackDeepMT$BinMean - JackDeepMT$SdEst*2), ylim=p8$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p9d, split=c(4,2,4,3), more=T)
plot(p9, split = c(4,2,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackCWV$AvgDist, y = (JackCWV$BinMean + JackCWV$SdEst*2), ylim=p9$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackCWV$AvgDist, y = (JackCWV$BinMean - JackCWV$SdEst*2), ylim=p9$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepCWV$AvgDist, y = (JackDeepCWV$BinMean + JackDeepCWV$SdEst*2), ylim=p9$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepCWV$AvgDist, y = (JackDeepCWV$BinMean - JackDeepCWV$SdEst*2), ylim=p9$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p10dzz, split=c(2,3,4,3), more=T)
plot(p10zz, split = c(2,3,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackVR$AvgDist, y = (JackVR$BinMean + JackVR$SdEst*2), ylim=p10zz$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackVR$AvgDist, y = (JackVR$BinMean - JackVR$SdEst*2), ylim=p10zz$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepVR$AvgDist, y = (JackDeepVR$BinMean + JackDeepVR$SdEst*2), ylim=p10zz$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepVR$AvgDist, y = (JackDeepVR$BinMean - JackDeepVR$SdEst*2), ylim=p10zz$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

dev.off()

png("Variograms_UniqueYAxis_CompareESDA_95ConfInt_DeepPost_RmOps.png", width=13, height=10, units="in", res=300)
plot(p1dRmOps, split=c(1,1,4,3), more=T)
plot(p1RmOps, split = c(1,1,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackCT_RmOps$AvgDist, y = (JackCT_RmOps$BinMean + JackCT_RmOps$SdEst*2), ylim=p1RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackCT_RmOps$AvgDist, y = (JackCT_RmOps$BinMean - JackCT_RmOps$SdEst*2), ylim=p1RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepCT_RmOps$AvgDist, y = (JackDeepCT_RmOps$BinMean + JackDeepCT_RmOps$SdEst*2), ylim=p1RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepCT_RmOps$AvgDist, y = (JackDeepCT_RmOps$BinMean - JackDeepCT_RmOps$SdEst*2), ylim=p1RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p2dRmOps, split=c(4,1,4,3), more=T)
plot(p2RmOps, split = c(4,1,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackCNY_RmOps$AvgDist, y = (JackCNY_RmOps$BinMean + JackCNY_RmOps$SdEst*2), ylim=p2RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackCNY_RmOps$AvgDist, y = (JackCNY_RmOps$BinMean - JackCNY_RmOps$SdEst*2), ylim=p2RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepCNY_RmOps$AvgDist, y = (JackDeepCNY_RmOps$BinMean + JackDeepCNY_RmOps$SdEst*2), ylim=p2RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepCNY_RmOps$AvgDist, y = (JackDeepCNY_RmOps$BinMean - JackDeepCNY_RmOps$SdEst*2), ylim=p2RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p3dRmOps, split=c(1,2,4,3), more=T)
plot(p3RmOps, split = c(1,2,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackENY_RmOps$AvgDist, y = (JackENY_RmOps$BinMean + JackENY_RmOps$SdEst*2), ylim=p3RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackENY_RmOps$AvgDist, y = (JackENY_RmOps$BinMean - JackENY_RmOps$SdEst*2), ylim=p3RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepENY_RmOps$AvgDist, y = (JackDeepENY_RmOps$BinMean + JackDeepENY_RmOps$SdEst*2), ylim=p3RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepENY_RmOps$AvgDist, y = (JackDeepENY_RmOps$BinMean - JackDeepENY_RmOps$SdEst*2), ylim=p3RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p4dRmOps, split=c(2,2,4,3), more=T)
plot(p4RmOps, split = c(2,2,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackENYPA_RmOps$AvgDist, y = (JackENYPA_RmOps$BinMean + JackENYPA_RmOps$SdEst*2), ylim=p4RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackENYPA_RmOps$AvgDist, y = (JackENYPA_RmOps$BinMean - JackENYPA_RmOps$SdEst*2), ylim=p4RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepENYPA_RmOps$AvgDist, y = (JackDeepENYPA_RmOps$BinMean + JackDeepENYPA_RmOps$SdEst*2), ylim=p4RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepENYPA_RmOps$AvgDist, y = (JackDeepENYPA_RmOps$BinMean - JackDeepENYPA_RmOps$SdEst*2), ylim=p4RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p5dRmOps, split=c(3,1,4,3), more=T)
plot(p5RmOps, split = c(3,1,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackNWPANY_RmOps$AvgDist, y = (JackNWPANY_RmOps$BinMean + JackNWPANY_RmOps$SdEst*2), ylim=p5RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackNWPANY_RmOps$AvgDist, y = (JackNWPANY_RmOps$BinMean - JackNWPANY_RmOps$SdEst*2), ylim=p5RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepNWPANY_RmOps$AvgDist, y = (JackDeepNWPANY_RmOps$BinMean + JackDeepNWPANY_RmOps$SdEst*2), ylim=p5RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepNWPANY_RmOps$AvgDist, y = (JackDeepNWPANY_RmOps$BinMean - JackDeepNWPANY_RmOps$SdEst*2), ylim=p5RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p6dRmOps, split=c(3,2,4,3), more=T)
plot(p6RmOps, split = c(3,2,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackSWPA_RmOps$AvgDist, y = (JackSWPA_RmOps$BinMean + JackSWPA_RmOps$SdEst*2), ylim=p6RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackSWPA_RmOps$AvgDist, y = (JackSWPA_RmOps$BinMean - JackSWPA_RmOps$SdEst*2), ylim=p6RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepSWPA_RmOps$AvgDist, y = (JackDeepSWPA_RmOps$BinMean + JackDeepSWPA_RmOps$SdEst*2), ylim=p6RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepSWPA_RmOps$AvgDist, y = (JackDeepSWPA_RmOps$BinMean - JackDeepSWPA_RmOps$SdEst*2), ylim=p6RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p7dRmOps, split=c(2,1,4,3), more=T)
plot(p7RmOps, split = c(2,1,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackWPA_RmOps$AvgDist, y = (JackWPA_RmOps$BinMean + JackWPA_RmOps$SdEst*2), ylim=p7RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackWPA_RmOps$AvgDist, y = (JackWPA_RmOps$BinMean - JackWPA_RmOps$SdEst*2), ylim=p7RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepWPA_RmOps$AvgDist, y = (JackDeepWPA_RmOps$BinMean + JackDeepWPA_RmOps$SdEst*2), ylim=p7RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepWPA_RmOps$AvgDist, y = (JackDeepWPA_RmOps$BinMean - JackDeepWPA_RmOps$SdEst*2), ylim=p7RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p8dRmOps, split=c(1,3,4,3), more=T)
plot(p8RmOps, split = c(1,3,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackMT_RmOps$AvgDist, y = (JackMT_RmOps$BinMean + JackMT_RmOps$SdEst*2), ylim=p8RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackMT_RmOps$AvgDist, y = (JackMT_RmOps$BinMean - JackMT_RmOps$SdEst*2), ylim=p8RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepMT_RmOps$AvgDist, y = (JackDeepMT_RmOps$BinMean + JackDeepMT_RmOps$SdEst*2), ylim=p8RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepMT_RmOps$AvgDist, y = (JackDeepMT_RmOps$BinMean - JackDeepMT_RmOps$SdEst*2), ylim=p8RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p9dRmOps, split=c(4,2,4,3), more=T)
plot(p9RmOps, split = c(4,2,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackCWV_RmOps$AvgDist, y = (JackCWV_RmOps$BinMean + JackCWV_RmOps$SdEst*2), ylim=p9RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackCWV_RmOps$AvgDist, y = (JackCWV_RmOps$BinMean - JackCWV_RmOps$SdEst*2), ylim=p9RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepCWV_RmOps$AvgDist, y = (JackDeepCWV_RmOps$BinMean + JackDeepCWV_RmOps$SdEst*2), ylim=p9RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepCWV_RmOps$AvgDist, y = (JackDeepCWV_RmOps$BinMean - JackDeepCWV_RmOps$SdEst*2), ylim=p9RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

plot(p10dRmOps, split=c(2,3,4,3), more=T)
plot(p10RmOps, split = c(2,3,4,3), more=T)
trellis.focus('panel', 1,1)
lpoints(x = JackVR_RmOps$AvgDist, y = (JackVR_RmOps$BinMean + JackVR_RmOps$SdEst*2), ylim=p10RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackVR_RmOps$AvgDist, y = (JackVR_RmOps$BinMean - JackVR_RmOps$SdEst*2), ylim=p10RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='red')
lpoints(x = JackDeepVR_RmOps$AvgDist, y = (JackDeepVR_RmOps$BinMean + JackDeepVR_RmOps$SdEst*2), ylim=p10RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
lpoints(x = JackDeepVR_RmOps$AvgDist, y = (JackDeepVR_RmOps$BinMean - JackDeepVR_RmOps$SdEst*2), ylim=p10RmOps$y.limits, pch=1, type='o', lty=2, cex=0.7, lwd=1.5, col='blue')
trellis.unfocus()

dev.off()

#    Table of variogram lag estimates at close separation distances----
#Separation distance targets
SepTargets = c(2500, 5000, 10000, 20000)
#Dataframe to store the table of values
VarLagsTable = data.frame(GeologicRegion = '', ProcessingStep = 0, SeparationDistanceLag = 0, Semivariance = 0, stringsAsFactors = FALSE)

#Get separation distance lags closest to the target values.
CTtargets = vector('numeric', length = length(SepTargets))
DeepCTtargets = vector('numeric', length = length(SepTargets))
PreCTtargets = vector('numeric', length = length(SepTargets))
for (i in 1:length(SepTargets)){
  CTtargets[i] = which(abs(v.CT$dist - SepTargets[i]) == min(abs(v.CT$dist - SepTargets[i])))
  DeepCTtargets[i] = which(abs(v.DeepCT$dist - SepTargets[i]) == min(abs(v.DeepCT$dist - SepTargets[i])))
  PreCTtargets[i] = which(abs(v.PreCT$dist - SepTargets[i]) == min(abs(v.PreCT$dist - SepTargets[i])))
}
VarLagsTable[1:12,-1] = cbind(rbind(cbind(rep(1,length(SepTargets))), cbind(rep(2,length(SepTargets))), cbind(rep(3,length(SepTargets)))), 
                        rbind(round(v.PreCT[PreCTtargets,2:3],1), round(v.DeepCT[DeepCTtargets,2:3],1), round(v.CT[CTtargets,2:3],1)))
VarLagsTable[1:12,1] = 'Chautauqua, NY'

WPAtargets = vector('numeric', length = length(SepTargets))
DeepWPAtargets = vector('numeric', length = length(SepTargets))
PreWPAtargets = vector('numeric', length = length(SepTargets))
for (i in 1:length(SepTargets)){
  WPAtargets[i] = which(abs(v.WPA$dist - SepTargets[i]) == min(abs(v.WPA$dist - SepTargets[i])))
  DeepWPAtargets[i] = which(abs(v.DeepWPA$dist - SepTargets[i]) == min(abs(v.DeepWPA$dist - SepTargets[i])))
  PreWPAtargets[i] = which(abs(v.PreWPA$dist - SepTargets[i]) == min(abs(v.PreWPA$dist - SepTargets[i])))
}
VarLagsTable[13:24,-1] = cbind(rbind(cbind(rep(1,length(SepTargets))), cbind(rep(2,length(SepTargets))), cbind(rep(3,length(SepTargets)))), 
                              rbind(round(v.PreWPA[PreWPAtargets,2:3],1), round(v.DeepWPA[DeepWPAtargets,2:3],1), round(v.WPA[WPAtargets,2:3],1)))
VarLagsTable[13:24,1] = 'Western PA'

NWPANYtargets = vector('numeric', length = length(SepTargets))
DeepNWPANYtargets = vector('numeric', length = length(SepTargets))
PreNWPANYtargets = vector('numeric', length = length(SepTargets))
for (i in 1:length(SepTargets)){
  NWPANYtargets[i] = which(abs(v.NWPANY$dist - SepTargets[i]) == min(abs(v.NWPANY$dist - SepTargets[i])))
  DeepNWPANYtargets[i] = which(abs(v.DeepNWPANY$dist - SepTargets[i]) == min(abs(v.DeepNWPANY$dist - SepTargets[i])))
  PreNWPANYtargets[i] = which(abs(v.PreNWPANY$dist - SepTargets[i]) == min(abs(v.PreNWPANY$dist - SepTargets[i])))
}
VarLagsTable[25:36,-1] = cbind(rbind(cbind(rep(1,length(SepTargets))), cbind(rep(2,length(SepTargets))), cbind(rep(3,length(SepTargets)))), 
                               rbind(round(v.PreNWPANY[PreNWPANYtargets,2:3],1), round(v.DeepNWPANY[DeepNWPANYtargets,2:3],1), round(v.NWPANY[NWPANYtargets,2:3],1)))
VarLagsTable[25:36,1] = 'Northwestern PA & NY'

CNYtargets = vector('numeric', length = length(SepTargets))
DeepCNYtargets = vector('numeric', length = length(SepTargets))
PreCNYtargets = vector('numeric', length = length(SepTargets))
for (i in 1:length(SepTargets)){
  CNYtargets[i] = which(abs(v.CNY$dist - SepTargets[i]) == min(abs(v.CNY$dist - SepTargets[i])))
  DeepCNYtargets[i] = which(abs(v.DeepCNY$dist - SepTargets[i]) == min(abs(v.DeepCNY$dist - SepTargets[i])))
  PreCNYtargets[i] = which(abs(v.PreCNY$dist - SepTargets[i]) == min(abs(v.PreCNY$dist - SepTargets[i])))
}
VarLagsTable[37:48,-1] = cbind(rbind(cbind(rep(1,length(SepTargets))), cbind(rep(2,length(SepTargets))), cbind(rep(3,length(SepTargets)))), 
                               rbind(round(v.PreCNY[PreCNYtargets,2:3],1), round(v.DeepCNY[DeepCNYtargets,2:3],1), round(v.CNY[CNYtargets,2:3],1)))
VarLagsTable[37:48,1] = 'Central NY'

ENYtargets = vector('numeric', length = length(SepTargets))
DeepENYtargets = vector('numeric', length = length(SepTargets))
PreENYtargets = vector('numeric', length = length(SepTargets))
for (i in 1:length(SepTargets)){
  ENYtargets[i] = which(abs(v.ENY$dist - SepTargets[i]) == min(abs(v.ENY$dist - SepTargets[i])))
  DeepENYtargets[i] = which(abs(v.DeepENY$dist - SepTargets[i]) == min(abs(v.DeepENY$dist - SepTargets[i])))
  PreENYtargets[i] = which(abs(v.PreENY$dist - SepTargets[i]) == min(abs(v.PreENY$dist - SepTargets[i])))
}
VarLagsTable[49:60,-1] = cbind(rbind(cbind(rep(1,length(SepTargets))), cbind(rep(2,length(SepTargets))), cbind(rep(3,length(SepTargets)))), 
                               rbind(round(v.PreENY[PreENYtargets,2:3],1), round(v.DeepENY[DeepENYtargets,2:3],1), round(v.ENY[ENYtargets,2:3],1)))
VarLagsTable[49:60,1] = 'Eastern NY'

ENYPAtargets = vector('numeric', length = length(SepTargets))
DeepENYPAtargets = vector('numeric', length = length(SepTargets))
PreENYPAtargets = vector('numeric', length = length(SepTargets))
for (i in 1:length(SepTargets)){
  ENYPAtargets[i] = which(abs(v.ENYPA$dist - SepTargets[i]) == min(abs(v.ENYPA$dist - SepTargets[i])))
  DeepENYPAtargets[i] = which(abs(v.DeepENYPA$dist - SepTargets[i]) == min(abs(v.DeepENYPA$dist - SepTargets[i])))
  PreENYPAtargets[i] = which(abs(v.PreENYPA$dist - SepTargets[i]) == min(abs(v.PreENYPA$dist - SepTargets[i])))
}
VarLagsTable[61:72,-1] = cbind(rbind(cbind(rep(1,length(SepTargets))), cbind(rep(2,length(SepTargets))), cbind(rep(3,length(SepTargets)))), 
                               rbind(round(v.PreENYPA[PreENYPAtargets,2:3],1), round(v.DeepENYPA[DeepENYPAtargets,2:3],1), round(v.ENYPA[ENYPAtargets,2:3],1)))
VarLagsTable[61:72,1] = 'Eastern NY & PA'

SWPAtargets = vector('numeric', length = length(SepTargets))
DeepSWPAtargets = vector('numeric', length = length(SepTargets))
PreSWPAtargets = vector('numeric', length = length(SepTargets))
for (i in 1:length(SepTargets)){
  SWPAtargets[i] = which(abs(v.SWPA$dist - SepTargets[i]) == min(abs(v.SWPA$dist - SepTargets[i])))
  DeepSWPAtargets[i] = which(abs(v.DeepSWPA$dist - SepTargets[i]) == min(abs(v.DeepSWPA$dist - SepTargets[i])))
  PreSWPAtargets[i] = which(abs(v.PreSWPA$dist - SepTargets[i]) == min(abs(v.PreSWPA$dist - SepTargets[i])))
}
VarLagsTable[73:84,-1] = cbind(rbind(cbind(rep(1,length(SepTargets))), cbind(rep(2,length(SepTargets))), cbind(rep(3,length(SepTargets)))), 
                               rbind(round(v.PreSWPA[PreSWPAtargets,2:3],1), round(v.DeepSWPA[DeepSWPAtargets,2:3],1), round(v.SWPA[SWPAtargets,2:3],1)))
VarLagsTable[73:84,1] = 'Southwestern PA'

CWVtargets = vector('numeric', length = length(SepTargets))
DeepCWVtargets = vector('numeric', length = length(SepTargets))
PreCWVtargets = vector('numeric', length = length(SepTargets))
for (i in 1:length(SepTargets)){
  CWVtargets[i] = which(abs(v.CWV$dist - SepTargets[i]) == min(abs(v.CWV$dist - SepTargets[i])))
  DeepCWVtargets[i] = which(abs(v.DeepCWV$dist - SepTargets[i]) == min(abs(v.DeepCWV$dist - SepTargets[i])))
  PreCWVtargets[i] = which(abs(v.PreCWV$dist - SepTargets[i]) == min(abs(v.PreCWV$dist - SepTargets[i])))
}
VarLagsTable[85:96,-1] = cbind(rbind(cbind(rep(1,length(SepTargets))), cbind(rep(2,length(SepTargets))), cbind(rep(3,length(SepTargets)))), 
                               rbind(round(v.PreCWV[PreCWVtargets,2:3],1), round(v.DeepCWV[DeepCWVtargets,2:3],1), round(v.CWV[CWVtargets,2:3],1)))
VarLagsTable[85:96,1] = 'Central WV'

MTtargets = vector('numeric', length = length(SepTargets))
DeepMTtargets = vector('numeric', length = length(SepTargets))
PreMTtargets = vector('numeric', length = length(SepTargets))
for (i in 1:length(SepTargets)){
  MTtargets[i] = which(abs(v.MT$dist - SepTargets[i]) == min(abs(v.MT$dist - SepTargets[i])))
  DeepMTtargets[i] = which(abs(v.DeepMT$dist - SepTargets[i]) == min(abs(v.DeepMT$dist - SepTargets[i])))
  PreMTtargets[i] = which(abs(v.PreMT$dist - SepTargets[i]) == min(abs(v.PreMT$dist - SepTargets[i])))
}
VarLagsTable[97:108,-1] = cbind(rbind(cbind(rep(1,length(SepTargets))), cbind(rep(2,length(SepTargets))), cbind(rep(3,length(SepTargets)))), 
                               rbind(round(v.PreMT[PreMTtargets,2:3],1), round(v.DeepMT[DeepMTtargets,2:3],1), round(v.MT[MTtargets,2:3],1)))
VarLagsTable[97:108,1] = 'Western WV'

VRtargets = vector('numeric', length = length(SepTargets))
DeepVRtargets = vector('numeric', length = length(SepTargets))
PreVRtargets = vector('numeric', length = length(SepTargets))
for (i in 1:length(SepTargets)){
  VRtargets[i] = which(abs(v.VR$dist - SepTargets[i]) == min(abs(v.VR$dist - SepTargets[i])))
  DeepVRtargets[i] = which(abs(v.DeepVR$dist - SepTargets[i]) == min(abs(v.DeepVR$dist - SepTargets[i])))
  PreVRtargets[i] = which(abs(v.PreVR$dist - SepTargets[i]) == min(abs(v.PreVR$dist - SepTargets[i])))
}
VarLagsTable[109:120,-1] = cbind(rbind(cbind(rep(1,length(SepTargets))), cbind(rep(2,length(SepTargets))), cbind(rep(3,length(SepTargets)))), 
                                rbind(round(v.PreVR[PreVRtargets,2:3],1), round(v.DeepVR[DeepVRtargets,2:3],1), round(v.VR[VRtargets,2:3],1)))
VarLagsTable[109:120,1] = 'Valley and Ridge'

#Write table
write.csv(VarLagsTable, file = 'TableVarianceLags.csv')

#  Operator data----
#Wow! This region becomes manageable with its variogram by removing these points
#Add operators and Waco drilled wells to the dataset
MT_Waco = MT
MT_Waco@data$Waco = 0
for (i in 1:nrow(MT_Waco)){
  if (length(which(MT_Waco$StateID[i] %in% Wacos$StateID)) > 0){
    MT_Waco@data$Waco[i] = 1
  }
}
MT_Waco$Waco[MT_Waco$StateID == 'WV1336'] = 1
rm(i)

#Map - the region with these wells still has data.
plot(MT_Waco, pch = 16, cex = 0.1)
plot(MT_Waco[-which(MT_Waco$Operator == 'Waco Oil & Gas Co., Inc.'),], col = 'red', add = T, pch = 16, cex = 0.2 )
plot(MT_Waco[-which(MT_Waco$Waco == 1),], col = 'purple', add = T, pch = 16, cex = 0.2 )

#Variogram cloud for Western West Virginia region
plot(variogram(Qs~1, MT_Waco[-which(MT_Waco$Operator == 'Waco Oil & Gas Co., Inc.'),], cutoff=60000, cloud = TRUE))
plot(variogram(Qs~1, MT_Waco[-which(MT_Waco$Waco == 1),], cutoff=60000, cloud = TRUE))

test = plot(variogram(Qs~1, MT_Waco[-which(MT_Waco$Operator == 'Waco Oil & Gas Co., Inc.'),], cutoff=60000, cloud = TRUE), digitize = TRUE)

test = plot(variogram(Qs~1, MT_Waco[-which(MT_Waco$Operator == 'Waco Oil & Gas Co., Inc.'),][-248,][-246,][-245,], cutoff=60000, cloud = TRUE), digitize = TRUE)
hist(c(test$head, test$tail), breaks = 100000)
mode(c(test$head, test$tail))
#251
MT_Waco[-which(MT_Waco$Operator == 'Waco Oil & Gas Co., Inc.'),][-248,][-246,][-245,][251,]

p8 = plot(variogram(Qs~1, MT_Waco, cutoff=30000, width=30000/70), plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = 'red', cex=0.5, ylim = c(0,700), xlim = c(0,60000))
plot(p8)
p8 = plot(variogram(Qs~1, MT_Waco[-which(MT_Waco$Operator == 'Waco Oil & Gas Co., Inc.'),][-248,][-246,][-245,][-251,], cutoff=30000, width=30000/70), plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = 'red', cex=0.5)
plot(p8)
p8 = plot(variogram(Qs~1, MT_Waco[-which(MT_Waco$Operator == 'Waco Oil & Gas Co., Inc.'),][-248,][-246,][-245,][-251,], cutoff=60000, width=60000/50), plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = 'red', cex=0.5, ylim = c(0,700), xlim = c(0,60000))
plot(p8)

p9 = plot(variogram(Qs~1, CWV, cutoff=60000, width=60000/500), plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = 'red', cex=0.5)
plot(p9)
p9 = plot(variogram(Qs~1, CWV[-which(CWV$Operator == 'Waco Oil & Gas Co., Inc.'),], cutoff=60000, width=60000/500), plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = 'red', cex=0.5)
plot(p9)

p8 = plot(variogram(Qs~1, MT_Waco, cutoff=60000, width=60000/50), plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = 'red', cex=0.5, ylim = c(0,700), xlim = c(0,60000))

png('VariogramMT_WacoOperatorRemoved.png', res = 300, height = 5, width = 5, units = 'in')
plot(p8p, more=T)
plot(p8d, more=T)
plot(p8, more=T)
plot(plot(variogram(Qs~1, MT_Waco[-which(MT_Waco$Waco == 1),], cutoff=60000, width = 60000/50), plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', ylim = c(0,700), xlim = c(0,60000), col = 'green', pch = 16, cex = 0.5), more=F)
dev.off()

# Fixme: Variogram point cloud for outlier analysis----
#  Uses dataset that is deeper than 1 km, and all points have unique spatial locations.
#  Analysis only for points within geologic regions
#test = plot(variogram(Qs~1, ENY, cutoff=60000, cloud = TRUE))
