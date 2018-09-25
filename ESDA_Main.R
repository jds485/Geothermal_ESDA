# Description----
#Code for Smith et al. paper on Exploratory Spatial Data Analysis Framework for Geothermal Resource Assessments: An Appalachian Basin Case Study

#Outline of Script:
# Load libraries and data
#  Note: the Loading Code and Loading Data section of this script can be skipped, and the input data loaded from the InputDataForR.Rdata file.
# Check for and remove points with negative geothermal gradients. 
# Check for wells with the same spatial location. Only the deepest well at the same spatial location is retained.
#  Special cases of different BHT at the same depth are handled by either assigning a more likely depth to the point, or averaging the data.
# Then performs a local spatial outlier detection and analysis.
# Then checks the performance of the ESDA methods using semi-variance compuations.

# Libraries ----
library(sp) # map plots
library(rgdal) #spatial data reading/writing
library(GISTools) #map making tools
library(dgof) #ks test for discrete distributions
library(Hmisc) #minor tick marks
library(readxl) #for Excel data reading
library(changepoint) #for changepoint analysis on the minimum BHT depth cutoff
library(gstat) #for variogram analysis

# Loading Code from Repositories ----
setwd("C:\\Users\\jsmif\\Documents\\Cornell\\Research\\Publications\\ESDA\\ESDACode\\Geothermal_ESDA")
source('LocalDeviation.R')
source('ColorFunctions.R')
setwd('C:\\Users\\jsmif\\Documents\\Cornell\\Research\\Publications\\DOE Grant\\ThermalConductivity\\Geothermal_DataAnalysis_CrossSections\\Geothermal_DataAnalysis_CrossSections')
source("DealingWithDataInDuplicateLocations.R")
setwd('C:\\Users\\jsmif\\Documents\\Cornell\\Research\\Publications\\DOE Grant\\CombiningRiskFactorCode\\geothermal_pfa\\outliers')
source('outlier_identification.R')
# Loading Data and Map Layers ----
#Wells with surface heat flow and temperatures at depth calculated
setwd("C:/Users/jsmif/Documents/Cornell/Research/Masters - Spatial Assessment/Figures")
#Wells = read.csv('EDAWells_AllTempsThicksConds_BaseCorr.csv', stringsAsFactors=FALSE)
Wells = read.csv('EDAWells_AllTempsThicksConds_BaseCorr_2018.csv', stringsAsFactors=FALSE)
coordinates(Wells) = c('LongDgr', 'LatDegr')
proj4string(Wells) = CRS("+init=epsg:4326")

#Spicer equilibrium well temperature profiles
setwd('C:\\Users\\jsmif\\Documents\\Cornell\\Research\\Publications\\ESDA')
Spicer = read_xlsx(path = paste0(getwd(), '/EquilibriumTempProfiles.xlsx'), sheet = 'Spicers')
coordinates(Spicer) = c("Long", "Lat")
proj4string(Spicer) = CRS("+init=epsg:4326")

#Whealton MS thesis identified pseudo-equilibrium temperature profiles
Whealton = read.csv('EquilibriumTempProfiles_LocsAdded.csv')
coordinates(Whealton) = c("Long", "Lat")
proj4string(Whealton) = CRS("+init=epsg:4326")

#Political boundaries
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

#Spatial interpolation regions
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

#Set working directory back to project directory
setwd("C:/Users/jsmif/Documents/Cornell/Research/Masters - Spatial Assessment/Figures")

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

#Remove the negative gradient wells before the sorting of wells in the same spatial locations:
Wells_PosGrad = Wells[-which(Wells$Gradient <= 0),]


# Identify Wells in Same Spatial Location ----
#Note that this step is used here so that the QsDev function to calculate the local
# median surface heat flow uses only unique locations.

#Find all points that share the same location and take the deepest measurement.
Same = SameSpot(Wells_PosGrad)
SortData = SortingWells(Same$SameSpot, Same$StoreData_Frame, Wells_PosGrad, 'BHT', 'TruVrtc', 'DrllrTt', 'DpthOfM', 'WellDepth', 2)
write.csv(SortData$Sorted, "SortedUniqueSpots_AllTemps_ESDA_2018.csv")
write.csv(SortData$RerunWells, "RerunWells_AllTemps_ESDA_2018.csv")

#Records in start of dataset
length(Wells_PosGrad)
#Number of locations with records in same spatial coordinates
nrow(unique(Wells_PosGrad@coords[as.numeric(colnames(Same$StoreData_Frame)),]))
#Unique spatial locations after sorting
length(SortData$Sorted)
#Records in same spatial location
length(Same$StoreData_Frame)
#Locations with multiple BHTs at same depth that should be rerun because they did not have depth field information
length(SortData$RerunWells)

#58 points in 27 unique locations have a different BHT measurement at the same depth.
BHTsDiffSameDepth = length(unique(SortData$IndsDifferent))
BHTsDiffSameDepth_UniqueLocs = nrow(unique(Wells_PosGrad[SortData$IndsDifferent,]@coords))
#Used to see how many wells had a CensorTemp controlled output. Max of about 13 C
SortData_TestCensor = SortingWells(Same$SameSpot, Same$StoreData_Frame, Wells_PosGrad, 'BHT', 'TruVrtc', 'DrllrTt', 'DpthOfM', 'WellDepth', Inf)

#Locations dropped as a result of censoring
nrow(unique(Wells_PosGrad@coords[unique(SortData$IndsCensTemp),]))

#The wells that are rerun should overwrite the last rows in SortedUniqueSopts
#Load in the wells that were rerun and add them to the SortData$Sorted
Rerun = read.csv('SortedUniqueSpots_AllTemps_ESDA_RerunAdded.csv', stringsAsFactors = FALSE)
Rerun = Rerun[c((nrow(Rerun) - (nrow(SortData$RerunWells) - 1)):nrow(Rerun)),]
Rerun$APINo = SortData$RerunWells$APINo
SortData$Sorted@data[c((nrow(SortData$Sorted) - (nrow(SortData$RerunWells) - 1)):nrow(SortData$Sorted)),] = Rerun[,-c(1,8,9,ncol(Rerun))]


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
rm(Rows,BHTs,i,Depths)

#  Line plot of the heat flow vs. depth of BHT measurement for the wells in the same spot----

#Fixme: Add the equilibrium and pseudo-equilibrium well data to this plot, or make a new plot for these data

# Make a copy of the database to track the wells in the same spatial location for this plot only.
PlotSpots = Wells_PosGrad@data
PlotSpots$LongDgr = Wells_PosGrad@coords[,1]
PlotSpots$LatDegr = Wells_PosGrad@coords[,2]
# The wells in the same spot will be assigned the same number in a field named SameSpot
PlotSpots$SameSpot = NA

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

#Sort by the same spot number
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

#Colors by location, sorted by the highest Qs to lowest in shallowest measurement
PlotColPal = colorRampPalette(colors = c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))
cols = PlotColPal(length(unique(PlotFinal$SameSpot)))

#Make plot
png('SameSpotWells_QsVsDepth_colrev.png', res = 600, units = 'in', width = 7, height = 7)
par(mar = c(4.5, 5, 1.5, 1.5), xaxs='i', yaxs='i')
for (i in 1:length(unique(PlotFinal$SameSpot))){
  if (i == 1){
    plot(PlotFinal$WellDepth[which(PlotFinal$SameSpot == PlotFinal$SameSpot[length(unique(PlotFinal$SameSpot)) + 1 - i])], PlotFinal$Qs[which(PlotFinal$SameSpot == PlotFinal$SameSpot[length(unique(PlotFinal$SameSpot)) + 1 - i])], type = 'o', col = cols[length(unique(PlotFinal$SameSpot)) + 1 - i], pch = 16, xlim = c(0,3500), ylim = c(0,250), xlab = 'Well Depth (m)', ylab = expression('Surface Heat Flow' ~ (mW/m^2)), cex.axis = 1.5, cex.lab = 1.5)
  }else{
    plot(PlotFinal$WellDepth[which(PlotFinal$SameSpot == PlotFinal$SameSpot[length(unique(PlotFinal$SameSpot)) + 1 - i])], PlotFinal$Qs[which(PlotFinal$SameSpot == PlotFinal$SameSpot[length(unique(PlotFinal$SameSpot)) + 1 - i])], type = 'o', col = cols[length(unique(PlotFinal$SameSpot)) + 1 - i], pch = 16, xlim = c(0,3500), ylim = c(0,250), axes = FALSE, xlab = '', ylab = '')
  }
  par(new = TRUE)
}
par(new = FALSE)
minor.tick(nx=5,ny=5)
legend('topright', legend = c('Shallowest Well is Shallow', '', '', 'Shallowest Well is Deep'), col = c('red', 'yellow', 'green', 'purple'), pch = 16, lty = 1)
dev.off()
rm(PlotColPal, cols)

#Transform to spatial data
coordinates(PlotFinal) = c('LongDgr', 'LatDegr')
proj4string(PlotFinal) = CRS('+init=epsg:4326')

#Map of which states have the most duplicate measurements
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

#Sort by BHT 
PlotFinal = PlotFinal[rev(order(PlotFinal$BHT)),]

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
rm(indsMaxDepth, Coords, i)

#Split by wells that have deeper BHT as smaller value

#Figure out which wells have the deeper BHT as smaller in value
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

PlotFinal_DeepBHTsLarge = PlotFinal[PlotFinal$IndsDeepSmallBHT == 0,]
PlotFinal_DeepBHTsSmall = PlotFinal[PlotFinal$IndsDeepSmallBHT == 1,]

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
    plot(rep(i,length(PlotFinal_DeepBHTsLarge$BHT[which(PlotFinal_DeepBHTsLarge$SameSpot == unique(PlotFinal_DeepBHTsLarge$SameSpot)[i])])), PlotFinal_DeepBHTsLarge$BHT[which(PlotFinal_DeepBHTsLarge$SameSpot == unique(PlotFinal_DeepBHTsLarge$SameSpot)[i])], type = 'o', col = 'black', pch = 16, xlim = c(0,length(unique(PlotFinal_DeepBHTsLarge$SameSpot))), ylim = c(0,140), xlab = 'Same Spot Well', ylab = expression(paste('BHT (', degree, 'C)')), main = 'Deepest BHT is Largest', cex.axis = 1.5, cex.lab = 1.5)
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
    plot(rep(i,length(PlotFinal_DeepBHTsSmall$BHT[which(PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i])])), PlotFinal_DeepBHTsSmall$BHT[which(PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i])], type = 'o', col = 'black', pch = 16, xlim = c(0,length(unique(PlotFinal_DeepBHTsSmall$SameSpot))), ylim = c(0,140), xlab = 'Same Spot Well', ylab = expression(paste('BHT (', degree, 'C)')), main = 'Deepest BHT is Not Largest', cex.axis = 1.5, cex.lab = 1.5)
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
rm(i, indsMaxDepth, Coords)

#Check how many of the Deep BHTs that are smaller are greater than 2 degrees or so
RowsDeepBHTSmallerBy2C_Counter = vector('numeric')
Diff = vector('numeric')
for (i in 1:length(unique(PlotFinal_DeepBHTsSmall$SameSpot))){
  BHTs = PlotFinal_DeepBHTsSmall$BHT[PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i]]
  Depths = PlotFinal_DeepBHTsSmall$WellDepth[PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i]]
  Diff = c(Diff, (max(BHTs) - max(BHTs[Depths == max(Depths)])))
  if ((max(BHTs) - max(BHTs[Depths == max(Depths)])) > 2){
    RowsDeepBHTSmallerBy2C_Counter = c(RowsDeepBHTSmallerBy2C_Counter, PlotFinal_DeepBHTsSmall$RowID_[PlotFinal_DeepBHTsSmall$SameSpot == unique(PlotFinal_DeepBHTsSmall$SameSpot)[i]])
  }
}
rm(i, BHTs, RowsDeepBHTSmallerBy2C_Counter, Depths)

#Histogram of the differences between the deepest and shallower BHTs
hist(Diff, breaks = 300)

#  Nugget Effect for wells in the same spatial location ----

#Fixme: Add nugget for equilibrium wells. 
#Fixme: Should this analysis use all wells, or only deep wells? Currently on all wells. Would need to move after well depth cutoff to do deep wells.

#Compute the Nugget Effect for points in the same spatial location.
#Make a data frame to store the locations, average nugget, number of nuggets calculated, min, max, and sd of the nugget
#Note: This provides the same result as using the PlotFinal, as done below.
LocsNugs2 = matrix(0, ncol=9, nrow=1)
colnames(LocsNugs2) = c('RowID_', 'POINT_X', 'POINT_Y', 'Nugget', 'Max', 'Min', 'Sd', 'PtPairs', 'NumPts')
count=0
#Mark the index with a 1 when it is used.
IndsUsed = vector('numeric', length=nrow(Same$StoreData_Frame))
for (i in 1:nrow(Same$StoreData_Frame)){
  #Only take the unique spots that have more than 1 point
  if (any(Same$StoreData_Frame[i,] == 1) & IndsUsed[i] != 1){
    #Gather all well indicies with the same spatial location.
    Indxs = as.numeric(colnames(Same$StoreData_Frame[which(Same$StoreData_Frame[i,] == 1)]))
    #Mark that this spatial location has now been checked by marking the indices.
    IndsUsed[which(colnames(Same$StoreData_Frame) %in% Indxs)] = 1
    #Check to see if any of the values have the same heat flow. That means the record was a duplicate, and the nugget should not be counted for these.
    Test = Wells_PosGrad$Qs[Indxs]
    if (length(unique(Test)) != 1){
      #There are unique BHTs for this well compute the nugget only for those wells that are unique records.
      Nug = vector('numeric', length=length(unique(Test)))
      VarioPts = matrix(0, nrow=length(Nug), ncol=length(Nug))
      for (j in 1:length(Nug)){
        VarioPts[j,] = ((unique(Test) - unique(Test)[j]))^2/2
      }
      Nug = VarioPts[lower.tri(VarioPts)]
      #Store spatial location of point and nugget information
      if (nrow(LocsNugs2) == 1 & Nug[1] != 0 & count == 0){
        LocsNugs2[1,] = c(Wells_PosGrad$RowID_[Indxs[1]], Wells_PosGrad$LongDgr[Indxs[1]], Wells_PosGrad$LatDegr[Indxs[1]], mean(Nug), max(Nug), min(Nug), sd(Nug), length(Nug), length(unique(Test)))
        count = 1
      }
      else if (count == 1){
        LocsNugs2 = rbind(LocsNugs2, c(Wells_PosGrad$RowID_[Indxs[1]], Wells_PosGrad$LongDgr[Indxs[1]], Wells_PosGrad$LatDegr[Indxs[1]], mean(Nug), max(Nug), min(Nug), sd(Nug), length(Nug), length(unique(Test))))
      }
    }
  }
}
rm(count, i, j, Indxs, Test, IndsUsed, Nug, VarioPts)

write.csv(LocsNugs2, 'NuggetLocations.csv')

#Using PlotFinal data
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

#Using PlotFinal data
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

#Make boxplots for each of the interpolation sections
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

#Plots for all data - No Map or Histogram
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

#  Local Median and Local Average Deviation----
#Uses the well data that has been sorted for unique spatial locations.
DataTab = QsDev(Data = SortData$Sorted@data, Var = 'Qs', xName = 'coords_x1', yName = 'coords_x2', rad = 10000, max_pts = 25)
#Add information to the original data
SortData$Sorted@data = DataTab
#Rename to shorter variable
WellsSort = SortData$Sorted

#These figures should be made with a dataset that has unique spatial locations, as completed above.

#Color function parameters for plotting the heat flow data map
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#With map
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
EDAPlotsMed(WellsSort, "Qs", Unit=expression(paste("Local Median Deviation", "(mW/m"^2, ")")),-100, 150)
axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.5)
axis(side = 2, at = seq(-200,1300,50), labels=TRUE, cex.axis=1.5)
minor.tick(nx = 5, ny = 10)
lines(c(-100,10000), c(0,0), lwd = 2)
lines(c(1000,1000),c(-1000,2000), lwd = 3)
lines(c(600,600),c(-1000,2000), lty=2, lwd = 3)
legend('topright', legend=c("New York", "Pennsylvania", "West Virginia", "Kentucky", "Maryland", "Virginia", "600 m", "1000 m"), pch=c(rep(16, 6),NA,NA), lty = c(rep(NA, 6), 2,1), lwd = 3, col=c("blue", "red", "green", "purple","orange","yellow","black", "black"), cex=1.7)
dev.off()

#Local Median Deviation - > 600 m well
DataTab_600 = QsDev(Data = SortData$Sorted@data[which(SortData$Sorted@data$WellDepth >= 600),], Var = 'Qs', xName = 'coords_x1', yName = 'coords_x2', rad = 10000, max_pts = 25)
#Rename to shorter variable
WellsSort_600 = WellsSort
WellsSort_600@data = DataTab_600

png("HeatFlow_LocMedDiff_600.png", width=8, height=8, units="in", res=600)
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
EDAPlotsMed(WellsSort_600, "Qs", Unit=expression(paste("Local Median Deviation", "(mW/m"^2, ")")),-100, 150)
axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.5)
axis(side = 2, at = seq(-200,1300,50), labels=TRUE, cex.axis=1.5)
minor.tick(nx = 5, ny = 10)
lines(c(-100,10000), c(0,0), lwd = 2)
lines(c(1000,1000),c(-1000,2000), lwd = 3)
lines(c(600,600),c(-1000,2000), lty=2, lwd = 3)
legend('topright', legend=c("New York", "Pennsylvania", "West Virginia", "Kentucky", "Maryland", "Virginia", "600 m", "1000 m"), pch=c(rep(16, 6),NA,NA), lty = c(rep(NA, 6), 2,1), lwd = 3, col=c("blue", "red", "green", "purple","orange","yellow","black", "black"), cex=1.7)
dev.off()

#Local Median Deviation - > 1000 m well
DataTab_1000 = QsDev(Data = SortData$Sorted@data[which(SortData$Sorted@data$WellDepth >= 1000),], Var = 'Qs', xName = 'coords_x1', yName = 'coords_x2', rad = 10000, max_pts = 25)
#Rename to shorter variable
WellsSort_1000 = WellsSort
WellsSort_1000@data = DataTab_1000

png("HeatFlow_LocMedDiff_1000.png", width=8, height=8, units="in", res=600)
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
OrderedWells$WellDepth[3855]

#   Add shallow data back into dataset for PA region ----
WellsSort$LatDeg = WellsSort@coords[,2]
WellsSort$LngDegr = WellsSort@coords[,1]
  
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
#  Spatial Outlier Detection ----
#This should be run after points have been reduced to unique locations and negative gradient values have been removed or corrected.

#Data must have a column of UTM coordinates in m for this to work because it relies on Euclidian distances.
WellsDeep = spTransform(WellsDeep, CRS('+init=epsg:26917'))
#Add columns of UTM coordinates
WellsDeep$POINT_X = WellsDeep@coords[,1]
WellsDeep$POINT_Y = WellsDeep@coords[,2]

TestedOutliers_HeatFlow = select_out_algo(Data = WellsDeep@data, 
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

#Convert to spatial data
coordinates(TestedOutliers_HeatFlow$NotOutliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow$NotOutliers) = CRS('+init=epsg:26917')
coordinates(TestedOutliers_HeatFlow$Outliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow$Outliers) = CRS('+init=epsg:26917')

coordinates(TestedOutliers_HeatFlow_max2k$NotOutliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_max2k$NotOutliers) = CRS('+init=epsg:26917')
coordinates(TestedOutliers_HeatFlow_max2k$Outliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_max2k$Outliers) = CRS('+init=epsg:26917')

coordinates(TestedOutliers_HeatFlow_min2k$NotOutliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_min2k$NotOutliers) = CRS('+init=epsg:26917')
coordinates(TestedOutliers_HeatFlow_min2k$Outliers) = c('POINT_X', 'POINT_Y')
proj4string(TestedOutliers_HeatFlow_min2k$Outliers) = CRS('+init=epsg:26917')

writeOGR(TestedOutliers_HeatFlow$NotOutliers, dsn=getwd(), layer="DeepestWells_NotOutliers_32km_Qs_CorrBase_Ranked_2018", driver = "ESRI Shapefile")
writeOGR(TestedOutliers_HeatFlow$Outliers, dsn=getwd(), layer="DeepestWells_Outliers_32km_Qs_CorrBase_Ranked_2018", driver = "ESRI Shapefile")

#Convert back to WGS coordinate system
Outs = spTransform(TestedOutliers_HeatFlow$Outliers, CRS = CRS("+init=epsg:4326"))
NoOuts = spTransform(TestedOutliers_HeatFlow$NotOutliers, CRS = CRS("+init=epsg:4326"))
AllData = rbind(NoOuts, Outs)

Outs_min2k = spTransform(TestedOutliers_HeatFlow_min2k$Outliers, CRS = CRS("+init=epsg:4326"))
NoOuts_min2k = spTransform(TestedOutliers_HeatFlow_min2k$NotOutliers, CRS = CRS("+init=epsg:4326"))

Outs_max2k = spTransform(TestedOutliers_HeatFlow_max2k$Outliers, CRS = CRS("+init=epsg:4326"))
NoOuts_max2k = spTransform(TestedOutliers_HeatFlow_max2k$NotOutliers, CRS = CRS("+init=epsg:4326"))

#   Map of outliers----
#Colors
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(10,90)
scaleBy = 20
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#Wells in NAD83 WGS84
WellsDeepWGS = spTransform(WellsDeep, CRSobj = CRS('+init=epsg:4326'))
NotOutliersWGS = spTransform(TestedOutliers_HeatFlow$NotOutliers, CRSobj = CRS('+init=epsg:4326'))
OutliersWGS = spTransform(TestedOutliers_HeatFlow$Outliers, CRSobj = CRS('+init=epsg:4326'))

png('LoHiOuts_Map.png', res = 1200, units = 'in', width = 14, height = 7)
layout(cbind(1,2))
par(xaxs = 'i', yaxs = 'i', mar = c(2,3,3,1))
plot(WellsDeepWGS, col = 'white', pch = 16, cex = 0.2, main = 'Low Outliers', cex.main = 2)
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(Counties[Counties$STATEFP %in% c(42,36,54,24,21,51),], border = 'grey', add=TRUE)
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
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(Counties[Counties$STATEFP %in% c(42,36,54,24,21,51),], border = 'grey', add=TRUE)
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
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE)
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
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE)
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

#   Depth Rank Empircal CDFs and KS Test ----
png('LowOutliers_wUniform.png', res = 300, units = 'in', width = 6, height = 6)
hist(Outs$out_loc_drank[which(Outs$out_loc_lo == 1)]*25, breaks = seq(0,25,1), main = 'Low Outliers Depth Rank', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=2, cex.axis = 1.5, ylim = c(0,35))
lines(c(0,25), c(length(Outs$out_loc_drank[which(Outs$out_loc_lo == 1)])/25,length(Outs$out_loc_drank[which(Outs$out_loc_lo == 1)])/25), lwd = 2, col='blue')
dev.off()

png('HighOutliers_wUniform.png', res = 300, units = 'in', width = 6, height = 6)
hist(Outs$out_loc_drank[which(Outs$out_loc_lo == 0)]*25, breaks = seq(0,25,1), main = 'High Outliers Depth Rank', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=2, cex.axis = 1.5, ylim = c(0,35))
lines(c(0,25), c(length(Outs$out_loc_drank[which(Outs$out_loc_lo == 0)])/25,length(Outs$out_loc_drank[which(Outs$out_loc_lo == 0)])/25), lwd = 2, col='red')
dev.off()

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

# Test Spatial Patterns in Outliers - by Depth
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


#   Q-Q plot for the points that are not outliers ----
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

#With Map
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
qqnorm(ENYPA@data$Qs) #Definite outlier that did not get removed is in here. 30 mW/m^3 greater than others. Remove it from ENYPA and FL.
FL = FL[-which(FL$RowID_ == 19770),]
ENYPA = ENYPA[which(ENYPA$Qs < 100),]
qqnorm(ENYPA@data$Qs)
qqline(ENYPA@data$Qs)
qqnorm(MT@data$Qs)
qqline(MT@data$Qs)
qqnorm(NWPANY@data$Qs) #The low heat flow of 7 was not tested for outliers. It should be removed.
FL = FL[-which(FL$RowID_ == 29908),]
NWPANY = NWPANY[which(NWPANY$Qs > 10),]
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

qqnorm(FL@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Full Region', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, ylim=c(0,120), xlim=c(-4,4))
qqline(FL@data$Qs)


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

#  Post analysis of variograms in each geologic region pre and post
#Transform to NAD UTM17N
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

#Using oringal dataset for points without negative gradients, and no points in same spatial location
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

#Using oringal dataset for points without negative gradients, and no points in same spatial location
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

#Compute variograms - All data, Data deeper than 1 km, Data that have been fully proessed. Plot all on same plot
v.CT <- variogram(Qs~1, CT, cutoff=60000, width=60000/150) 
v.CNY <- variogram(Qs~1, CNY, cutoff=60000, width=60000/10) 
v.CWV <- variogram(Qs~1, CWV, cutoff=60000, width=60000/350)
v.ENY <- variogram(Qs~1, ENY, cutoff=60000, width=60000/10)
v.ENYPA <- variogram(Qs~1, ENYPA, cutoff=40000, width=40000/120)
v.MT <- variogram(Qs~1, MT, cutoff=40000, width=40000/200)
v.NWPANY <- variogram(Qs~1, NWPANY, cutoff=60000, width=60000/20) 
v.SWPA <- variogram(Qs~1, SWPA, cutoff=60000, width=60000/200) 
v.WPA <- variogram(Qs~1, WPA, cutoff=60000, width=60000/50) 
v.VR <- variogram(Qs~1, VR, cutoff=60000, width=60000/50) 
v.FL <- variogram(Qs~1, FL, cutoff=60000, width=60000/200)

v.PreCT <- variogram(Qs~1, PreCT, cutoff=60000, width=60000/150) 
v.PreCNY <- variogram(Qs~1, PreCNY, cutoff=60000, width=60000/10) 
v.PreCWV <- variogram(Qs~1, PreCWV, cutoff=60000, width=60000/350)
v.PreENY <- variogram(Qs~1, PreENY, cutoff=60000, width=60000/10)
v.PreENYPA <- variogram(Qs~1, PreENYPA, cutoff=40000, width=40000/120)
v.PreMT <- variogram(Qs~1, PreMT, cutoff=40000, width=40000/200)
v.PreNWPANY <- variogram(Qs~1, PreNWPANY, cutoff=60000, width=60000/20) 
v.PreSWPA <- variogram(Qs~1, PreSWPA, cutoff=60000, width=60000/200) 
v.PreWPA <- variogram(Qs~1, PreWPA, cutoff=60000, width=60000/50) 
v.PreVR <- variogram(Qs~1, PreVR, cutoff=60000, width=60000/50) 
v.PreFL <- variogram(Qs~1, PreFL, cutoff=60000, width=60000/200)

v.DeepCT <- variogram(Qs~1, DeepCT, cutoff=60000, width=60000/150) 
v.DeepCNY <- variogram(Qs~1, DeepCNY, cutoff=60000, width=60000/10) 
v.DeepCWV <- variogram(Qs~1, DeepCWV, cutoff=60000, width=60000/350)
v.DeepENY <- variogram(Qs~1, DeepENY, cutoff=60000, width=60000/10)
v.DeepENYPA <- variogram(Qs~1, DeepENYPA, cutoff=40000, width=40000/120)
v.DeepMT <- variogram(Qs~1, DeepMT, cutoff=40000, width=40000/200)
v.DeepNWPANY <- variogram(Qs~1, DeepNWPANY, cutoff=60000, width=60000/20) 
v.DeepSWPA <- variogram(Qs~1, DeepSWPA, cutoff=60000, width=60000/200) 
v.DeepWPA <- variogram(Qs~1, DeepWPA, cutoff=60000, width=60000/50) 
v.DeepVR <- variogram(Qs~1, DeepVR, cutoff=60000, width=60000/50) 
v.DeepFL <- variogram(Qs~1, DeepFL, cutoff=60000, width=60000/200)

p1 = plot(v.CT, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Chautauqua, NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,50), col = 'red', xlim = c(0,60000), cex=0.5)
p2 = plot(v.CNY,plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central NY", xlab = 'Separation Distance (m)', pch = 16, ylim = c(0,1500), col = 'red', xlim = c(0,60000), cex=0.5)
p9 = plot(v.CWV, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Central WV", xlab = 'Separation Distance (m)', pch = 16, col = 'red', ylim = c(0,1200), xlim = c(0,60000), cex=0.5)
p3 = plot(v.ENY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY", xlab = 'Separation Distance (m)', pch = 16, col = 'red', ylim = c(0,1500), xlim = c(0,60000), cex=0.5)
p4 = plot(v.ENYPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Eastern NY & PA", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,400), cex=0.5)
p8 = plot(v.MT, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western WV", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p5 = plot(v.NWPANY, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Northwestern PA & NY", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,500), cex=0.5)
p6 = plot(v.SWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,850), cex=0.5)
p6z = plot(v.SWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Southwestern PA Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,60), cex=0.5)
p10 = plot(v.VR, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p10z = plot(v.VR, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p7 = plot(v.WPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p7z = plot(v.WPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p11 = plot(v.FL, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = 'red', xlim = c(0,60000), ylim = c(0,1000), cex=0.5)

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
p10d = plot(v.DeepVR, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,11000), cex=0.5)
p10dz = plot(v.DeepVR, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Valley and Ridge Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,300), cex=0.5)
p7d = plot(v.DeepWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,700), cex=0.5)
p7dz = plot(v.DeepWPA, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Western PA Zoom", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,50), cex=0.5)
p11d = plot(v.DeepFL, plot.numbers=F, ylab=expression(Semivariance ~ (mW/m^2)^2), main="Full Region", xlab = 'Separation Distance (m)', pch = 16, col = 'blue', xlim = c(0,60000), ylim = c(0,1000), cex=0.5)


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
