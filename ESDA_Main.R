#Load in the entire well database that has had the geotherms calculated. 

#This script will check for and remove wells with negative gradients. 
#Then check for wells with the same spatial location. Only the deepest well at the same spatial location will be retained.
#Then perform al local outlier analysis.

# Libraries ----
library(sp) # map plots
library(rgdal) #spatial data reading/writing
library(GISTools) #map making tools
library(dgof) #ks test for discrete distributions
library(Hmisc) #minor tick marks

# Loading Code from Repositories ----
setwd("C:\\Users\\Jared\\Documents\\Cornell\\Research\\Publications\\ESDA\\ESDACode\\Geothermal_ESDA")
source('LocalDeviation.R')
source('ColorFunctions.R')
setwd('C:\\Users\\Jared\\Documents\\Cornell\\Research\\Publications\\DOE Grant\\ThermalConductivity\\Geothermal_DataAnalysis_CrossSections\\Geothermal_DataAnalysis_CrossSections')
source("DealingWithDataInDuplicateLocations.R")
setwd('C:\\Users\\Jared\\Documents\\Cornell\\Research\\Publications\\DOE Grant\\CombiningRiskFactorCode\\geothermal_pfa\\outliers')
source('outlier_identification.R')
setwd("C:/Users/Jared/Documents/Cornell/Research/Masters - Spatial Assessment/Figures")
# Loading data and map layers ----
Data = read.csv('EDAWells_AllTempsThicksConds_BaseCorr.csv', stringsAsFactors=FALSE)
coordinates(Data) = c('LongDgr', 'LatDegr')
proj4string(Data) = CRS("+init=epsg:4326")
States = readOGR(dsn="C:/Users/Jared/Documents/Cornell/Research/Masters - Spatial Assessment/GIS/Population Density/2013_us_state_500k", layer="us_state_WGS", stringsAsFactors=FALSE)
NY = States[which(States$STATEFP == "36"),]
PA = States[which(States$STATEFP == "42"),]
WV = States[which(States$STATEFP == "54"),]
MD = States[which(States$STATEFP == "24"),]
KY = States[which(States$STATEFP == "21"),]
VA = States[which(States$STATEFP == "51"),]
Counties = readOGR(dsn="C:/Users/Jared/Documents/Cornell/Research/Masters - Spatial Assessment/GIS/Population Density/2013_us_countys_500k", layer="us_county_WGS", stringsAsFactors=FALSE) 
InterpRegs = readOGR(dsn="C:/Users/Jared/Documents/Cornell/Research/Masters - Spatial Assessment/Figures/BaseCorrWells", layer="AllSectionsMerged", stringsAsFactors=FALSE) 
InterpRegs = spTransform(InterpRegs, CRS = CRS("+init=epsg:4326"))
#50 km bounded regions within the potential field edges - only applies to the WV regions
VR_Bounded = readOGR(dsn="C:\\Users\\Jared\\Documents\\Cornell\\Research\\Publications\\DOE Grant\\InterpolationDataset\\NewBoundaries\\NAD_InterpolationBounds", layer = 'BoundedVREdit6') 
VR_Bounded = spTransform(VR_Bounded, CRS = CRS("+init=epsg:4326"))
CWV_Bounded = readOGR(dsn = 'C:\\Users\\Jared\\Documents\\Cornell\\Research\\Publications\\DOE Grant\\InterpolationDataset\\NewBoundaries\\NAD_InterpolationBounds', layer = 'BoundedCWV_Edit3')
CWV_Bounded = spTransform(CWV_Bounded, CRS = CRS("+init=epsg:4326"))
MT_Bounded = readOGR(dsn = 'C:\\Users\\Jared\\Documents\\Cornell\\Research\\Publications\\DOE Grant\\InterpolationDataset\\NewBoundaries\\NAD_InterpolationBounds', layer = 'BoundedMT_Edit')
MT_Bounded = spTransform(MT_Bounded, CRS = CRS("+init=epsg:4326"))

# Check for negative gradients ----
NegsAll = length(Data$Gradient[which(Data$Gradient < 0)])

#Deepest well that has a negative gradient
MaxDepth_NegGrad = max(Data$WellDepth[which(Data$Gradient < 0)])

#Location of wells that have negative Gradients.
plot(Data, pch=16, col='black')
plot(Data[which(Data$Gradient < 0),], pch=16, col='red', add=TRUE)
plot(NY, add=TRUE)
plot(PA, add=TRUE)
plot(WV, add=TRUE)
plot(MD, add=TRUE)
plot(KY, add=TRUE)
plot(VA, add=TRUE)

#Remove the negative gradient wells before the sorting of wells in the same spatial locations:
DataAll = Data[-which(Data$Gradient < 0),]


# Check Wells in same spatial location ----
#Note that this step is used here so that the QsDev function to calculate the local
# median surface heat flow uses only unique locations.

#Find all points that share the same location and take the deepest measurement.
Same = SameSpot(DataAll)
SortData = SortingWells(Same$SameSpot, Same$StoreData_Frame, DataAll, 'BHT', 'TruVrtc', 'DrllrTt', 'DpthOfM', 'WellDepth', 2)
write.csv(SortData$Sorted, "SortedUniqueSpots_AllTemps_ESDA.csv")
write.csv(SortData$RerunWells, "RerunWells_AllTemps_ESDA.csv")

#58 points in 27 unique locations have a different BHT measurement at the same depth.
length(unique(SortData$IndsDifferent))
nrow(unique(DataAll[SortData$IndsDifferent,]@coords))
#Used to see how many wells had a CensorTemp controlled output. Max of about 13 C
SortData_TestCensor = SortingWells(Same$SameSpot, Same$StoreData_Frame, DataAll, 'BHT', 'TruVrtc', 'DrllrTt', 'DpthOfM', 'WellDepth', 14)

#The wells that are rerun should overwrite the last rows in SortedUniqueSopts
#Load in the wells that were rerun and add them to the SortData$Sorted
Rerun = read.csv('SortedUniqueSpots_AllTemps_ESDA_RerunAdded.csv', stringsAsFactors = FALSE)
Rerun = Rerun[c((nrow(Rerun) - 4):nrow(Rerun)),]
Rerun$APINo = SortData$RerunWells$APINo
SortData$Sorted@data[c((nrow(SortData$Sorted) - 4):nrow(SortData$Sorted)),] = Rerun[,-c(1,8,9,ncol(Rerun))]

# Make a plot of the heat flow vs. depth of BHT measurement for the wells in the same spot----
# Make a copy of the database to track the wells in the same spatial location.
PlotSpots = DataAll@data
# The wells in the same spot will be assigned the same number in a field named SameSpot
PlotSpots$SameSpot = NA

count = 1

for (i in 1:nrow(Same$StoreData_Frame)){
  #Only take the unique spots that have more than 1 point
  if (any(Same$StoreData_Frame[i,] == 1) & is.na(PlotSpots$SameSpot[as.numeric(colnames(Same$StoreData_Frame)[i])]) == TRUE){
    #Have not checked this spot yet. Gather all well indicies with the same spatial location.
    Indxs = as.numeric(colnames(Same$StoreData_Frame[which(Same$StoreData_Frame[i,] == 1)]))
    
    #Assign a number to these wells in the same spot
    PlotSpots$SameSpot[Indxs] = count
    count = count + 1
  }
}

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

#Make new database for the final plotting of points
PlotFinal = PlotSpots

#index of data in the PlotFinal database
len = 1
for (i in 1:length(unique(PlotSpots$SameSpot))){
  #Determine how many wells there are in the same spot
  indxs = which(PlotSpots$SameSpot == NewInds[i])
  
  if (anyDuplicated(PlotSpots$Qs[indxs]) != 0){
    #Find the indices that have same heat flow
    res = which(PlotSpots$Qs[indxs] %in% unique(PlotSpots$Qs[indxs][duplicated(PlotSpots$Qs[indxs])]) == TRUE)
    #Of those, find indices with the same well depth
    dpth = which(PlotSpots$WellDepth[indxs][res] %in% unique(PlotSpots$WellDepth[indxs][res][duplicated(PlotSpots$WellDepth[indxs][res])]) == TRUE)
    
    #If they all have the same heat flow, check if they have the same depth.
    if ((length(res) == length(indxs)) & (length(dpth) == length(indxs))){
      #Do not record this in the PlotFinal database. Remove the last length(indxs) rows from the database.
      PlotFinal = PlotFinal[-((nrow(PlotFinal) - (length(indxs)-1)):nrow(PlotFinal)),]
    }else{
      PlotFinal[len:(len + (length(indxs)-1)),] = PlotSpots[indxs,]
      len = len + length(indxs)
    }
  }else{
    PlotFinal[len:(len + (length(indxs)-1)),] = PlotSpots[indxs,]
    len = len + length(indxs)
  }
}

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
dev.off()

# Fixme: Nugget for equilibrium wells. Need Calvin's dataset of temperature profiles. Should use only deep wells? ----
#Reliable Well API Numbers
Reliable = c("37005215920000", "37021209380000", "37021200780000", "37027206440000", "37027205360000", "37031240530000",
             "37033200710000", "37033200800000", "37033200900000", "37033201150000", "37033203260000", "37033205860000",
             "37033206090000", "37033206510000", "37033207030000", "37033212920000", "37033222810000", "37035900310000",
             "37051209960000", "37059242200000", "37059242560000", "37059242680000", "37063200530000", "37063235200000",
             "37063236380000", "37063274460000", "37063324730000", "37065203860000", "37083312520000", "37105210920000",
             "37111200630000", "37113200030000", "37121248370000", "37129224980000", "37129254590000", "37129261150000",
             "37129261160000", "37129262270000", "37129262420000", "37129263530000", "31009218090000", "31121121780000",
             "31121219000000", "37027200070000")

which(DataAll$APINo %in% Reliable)

#Compute the Nugget Effect for points in the same spatial location.
#Make a data frame to store the locations, average nugget, number of nuggets calculated, min, max, and sd of the nugget
LocsNugs = matrix(0, ncol=9, nrow=1)
colnames(LocsNugs) = c('RowID_', 'POINT_X', 'POINT_Y', 'Nugget', 'Max', 'Min', 'Sd', 'PtPairs', 'NumPts')
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
    Test = DataAll$Qs[Indxs]
    if (length(unique(Test)) != 1){
      #There are unique BHTs for this well compute the nugget only for those wells that are unique records.
      Nug = vector('numeric', length=length(unique(Test)))
      VarioPts = matrix(0, nrow=length(Nug), ncol=length(Nug))
      for (j in 1:length(Nug)){
        VarioPts[j,] = ((unique(Test) - unique(Test)[j]))^2/2
      }
      Nug = VarioPts[lower.tri(VarioPts)]
      #Store spatial location of point and nugget information
      if (nrow(LocsNugs) == 1 & Nug[1] != 0 & count == 0){
        LocsNugs[1,] = c(DataAll$RowID_[Indxs[1]], DataAll$LongDgr[Indxs[1]], DataAll$LatDegr[Indxs[1]], mean(Nug), max(Nug), min(Nug), sd(Nug), length(Nug), length(unique(Test)))
        count = 1
      }
      else if (count == 1){
        LocsNugs = rbind(LocsNugs, c(DataAll$RowID_[Indxs[1]], DataAll$LongDgr[Indxs[1]], DataAll$LatDegr[Indxs[1]], mean(Nug), max(Nug), min(Nug), sd(Nug), length(Nug), length(unique(Test))))
      }
    }
  }
}
rm(count, i, j, Indxs, Test)

write.csv(LocsNugs, 'NuggetLocations.csv')

#Make boxplots for each of the interpolation sections
NuggetWells = readOGR(dsn=getwd(), layer='Nugget_Locations_Sections')
png('NuggetWellsDistributions.png', res=600, width=12, height=6, units='in')
par(mar=c(5,5,2,2), yaxt='n')
boxplot(Nugget ~ Id, data = NuggetWells, at=c(9, 2, 7, 3, 6, 8, 1, 4, 5), varwidth = TRUE, col=c('gray', 'orange', 'blue', 'yellow', 'skyblue', 'purple', 'red', 'green', 'springgreen'), names=c('VR', 'WPA', 'SWPA', 'NWPANY', 'ENYPA', 'CWV', 'CT', 'CNY', 'ENY'), pch=16, cex.axis=1.5, cex.lab=1.5, ylab=expression('Sample Nugget Semivariance' ~ (mW/m^2)^2), xlab='Interpolation Region', log='y')
par(yaxt='s')
at.y <- outer(1:9, 10^(-7:7))
lab.y <- ifelse(log10(at.y) %% 1 == 0, sapply(log10(at.y),function(i)
  as.expression(bquote(10^ .(i)))), NA)
axis(side=2, at=at.y, labels=lab.y, cex.axis=1.5, las=1)
dev.off()

# With a map next to the boxplot
pin = par("pin")
dxy = apply(rbind(c(-82.64474, -74.5), c(36.75, 43.6)), 1, diff)
ratio = dxy[1]/dxy[2]
pin[1] = 3 #Modifying for margin space
png('NuggetWellsDistributions_Boxplot&Map_edit.png', res=600, width=11.2, height=4.2, units='in')
layout(cbind(1,1,2))
par(mar=c(4.1,5,0.5,0), yaxt='n')
boxplot(Nugget ~ Id, data = NuggetWells, at=c(9, 2, 7, 3, 6, 8, 1, 4, 5), varwidth = TRUE, col=c('gray', 'orange', 'blue', 'yellow', 'skyblue', 'purple', 'red', 'green', 'springgreen'), names=c('VR', 'WPA', 'SWPA', 'NWPANY', 'ENYPA', 'CWV', 'CT', 'CNY', 'ENY'), pch=16, cex.axis=1.5, cex.lab=1.5, ylab=expression('Average Sample Nugget Semivariance' ~ (mW/m^2)^2), xlab='Geologic Region', log='y')
par(yaxt='s')
at.y <- outer(1:9, 10^(-5:3))
lab.y <- ifelse(log10(at.y) %% 1 == 0, sapply(log10(at.y),function(i)
  as.expression(bquote(10^ .(i)))), NA)
axis(side=2, at=at.y, labels=lab.y, cex.axis=1.5, las=1)
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
plot(NuggetWells, pch = 16, cex = 0.5, add = TRUE)
north.arrow(-75, 37, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
dev.off()

# EDA for the well database ----

#Wells must be checked for negatives before this analysis.
#Plots for all data - No Map or Histogram. Can use the wells in the same spatial location.
sets = rbind(c(1,2,7,7), c(3,4,7,7), c(5,6,7,7))
layout(sets)
png("HeatFlowEDA_Thesis_SplitAll.png", width=9, height=9, units="in", res=600)
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
EDAPlots(DataAll, "Qs", Unit=expression("Surface Heat Flow" ~ (mW/m^2)),-10, 1400)
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
EDAPlots(DataAll, "Qs", Unit=expression("Surface Heat Flow" ~ (mW/m^2)),-10, 1400)
axis(side=2, at=seq(100,1300,100),labels=FALSE)
lines(c(1000,1000),c(-1000,1700))
lines(c(600,600),c(-1000,1700), lty=2)
legend('topright', legend=c("New York", "Pennsylvania", "West Virginia", "Kentucky", "Maryland", "Virginia"), pch=16, col=c("blue", "red", "green", "purple","orange","yellow"), cex=2)
dev.off()

#Color function parameters for plotting the heat flow data map
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
scaleRange = c(30,80)
scaleBy = 10
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)

#Calculate the deviation from the local median and the local average surface heat flow.
#Uses the well data that has been sorted for unique spatial locations.
DataTab = QsDev(Data = SortData$Sorted@data, Var = 'Qs', xName = 'coords_x1', yName = 'coords_x2', rad = 10000, max_pts = 25)
SortData$Sorted@data = DataTab

#Change the names of the datasets because these figures should be made with a dataset that has unique spatial locations.
DataUnsort = DataAll
DataAll = SortData$Sorted

#With map
sets = rbind(c(1,2,7,7), c(3,4,7,7), c(5,6,8,8), c(9,9,8,8))
layout(sets)
png("HeatFlowEDA_Thesis_SplitAll_Map_Log_WellDepthSort_MedDiff_DeepWells.png", width=10, height=10, units="in", res=600)
layout(sets)
EDAPlots = function(DataAll, Var, Unit, ymin, ymax){
  plot(DataAll$WellDepth[which(DataAll$State == "MD")], DataAll@data[Var][,1][which(DataAll$State == "MD")], col = "orange", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Maryland', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "KY")], DataAll@data[Var][,1][which(DataAll$State == "KY")], col = "purple", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Kentucky', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "VA")], DataAll@data[Var][,1][which(DataAll$State == "VA")], col = "yellow", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab=Unit, xlab='Depth of BHT (m)', main='Virginia', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "WV")], DataAll@data[Var][,1][which(DataAll$State == "WV")], col = "green", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='West Virginia', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "PA")], DataAll@data[Var][,1][which(DataAll$State == "PA")], col = "red", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='Pennsylvania', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
  plot(DataAll$WellDepth[which(DataAll$State == "NY")], DataAll@data[Var][,1][which(DataAll$State == "NY")], col = "blue", pch=16, xlim=c(0,7000), ylim=c(ymin,ymax), ylab='', xlab='Depth of BHT (m)', main='New York', cex.main=2, cex.axis=1.5, cex.lab=1.5, log = 'y', axes=FALSE)
  axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.2)
  axis(side = 2, at = 10^seq(0,4,1), labels=expression("10"^0, "10"^1, "10"^2, "10"^3, "10"^4), cex.axis=1.5, las = 1)
  axis(side = 2, at = c(seq(.2,.9,.1),seq(2,9,1),seq(20,90,10),seq(200,900,100), seq(2000,9000,1000)), labels = FALSE)
  lines(c(1000,1000),c(0.001,2000))
  lines(c(600,600),c(0.001,2000), lty=2)
}
par(mar=c(4,5.5,3,2), xaxs = 'i', yaxs = 'i')
EDAPlots(DataAll, "Qs", Unit=expression("Surface Heat Flow" ~ (mW/m^2)),.5, 2000)
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
# EDAPlots(DataAll, "Qs", Unit=expression("Surface Heat Flow" ~ (mW/m^2)),.5, 2000)
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
#EDAPlotsMean(DataAll, "Qs", Unit=expression("Deviation from Local Average Surface Heat Flow" ~ (mW/m^2)),-200, 200)
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
EDAPlotsMed(DataAll, "Qs", Unit=expression(paste("Site Q"['s'], " - Local Median Q"['s'], " (mW/m"^2, ")")),-100, 300)
axis(side = 1, at = seq(0,7000,1000), labels = TRUE, cex.axis = 1.5)
axis(side = 2, at = seq(-200,1300,50), labels=TRUE, cex.axis=1.5)
minor.tick(nx = 5, ny = 10)
lines(c(-100,10000), c(0,0))
lines(c(1000,1000),c(-1000,2000))
lines(c(600,600),c(-1000,2000), lty=2)
legend('topright', legend=c("New York", "Pennsylvania", "West Virginia", "Kentucky", "Maryland", "Virginia"), pch=16, col=c("blue", "red", "green", "purple","orange","yellow"), cex=2)
par(mar = c(2,2,1.2,1), xaxs = 'i', yaxs = 'i')
plot(DataAll, pch = 16, col = "white", cex = 0.5)
plot(NY, lwd = 2, add=TRUE)
plot(PA, lwd = 2, add=TRUE)
plot(WV, lwd = 2, add=TRUE)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], add=TRUE)
#plot(DataAll[order(DataAll$Qs),], pch = 16, col = colFun(DataAll$Qs[order(DataAll$Qs)]), cex = 0.5, add = TRUE)
plot(DataAll[order(DataAll$WellDepth, decreasing = FALSE),], pch = 16, col = colFun(DataAll$Qs[order(DataAll$WellDepth, decreasing = FALSE)]), cex = 0.5, add = TRUE)
#plot(DataAll, pch = 16, col = colFun(DataAll$Qs), cex = 0.5)
#plot(WellsSorted[order(WellsSorted$WellDepth),], pch = 16, col = colFun(WellsSorted$Qs[order(WellsSorted$WellDepth)]), cex = 0.5)
north.arrow(-75, 37.5, 0.1, lab = 'N', col='black', cex = 1.5)
degAxis(side = 2, seq(34, 46, 2), cex.axis = 1.5)
degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -2), cex.axis = 1.5)
degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
legend('topleft', title = expression(paste('Q'['s'], ' (mW/m'^2,')' )), legend = c('<40', '40 - 50', '50 - 60', '60 - 70', '70 - 80', '>80'), pch = 16, cex = 1.8, col = colFun(c(30, 45, 55, 65, 75, 90)))
par(mar=c(4,5.5,3,2))
hist(DataAll$Qs, breaks = c(seq(0, 130, 5), 1500), freq = FALSE, xlim = c(0,120), xlab = expression(paste('Surface Heat Flow (mW/m'^2,')')), cex.lab = 1.5, cex.axis = 1.5, main = 'Histogram of All Data', cex.main = 2)
dev.off()

# Add shallow data back into dataset for PA region ----
DataAll$LatDeg = DataAll@coords[,2]
DataAll$LngDegr = DataAll@coords[,1]
  
#Northwestern PA - All constraints are to focus on only the region of interest, and excludes all other wells in these counties.
sets = rbind(c(1,2,8,8), c(3,4,7,7), c(5,6,7,7))
layout(sets)
png('NWPennsylvaniaShallowerWells_DeepWells.png', width=9, height=9, units="in", res=600)
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
EDAPlots(DataAll, "Qs", Unit=expression("Surface Heat Flow" ~ (mW/m^2)), ymin = 0, ymax = 200)
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
EDAPlots(DataAll, "Qs", Unit=expression("Surface Heat Flow" ~ (mW/m^2)),0, 200)
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

#Remove all wells with depths less than 1000 m
WellsSorted = DataAll[-which(DataAll$WellDepth < 1000),]
#Then add back the PA wells.
WellsSorted = rbind(WellsSorted, DataAll[which(DataAll$WellDepth > 750 & DataAll$County == 'MC KEAN' & DataAll$State == 'PA' & DataAll$WellDepth < 1000),])
WellsSorted = rbind(WellsSorted, DataAll[which(DataAll$WellDepth > 750 & DataAll$WellDepth < 1000 & DataAll$County == 'ELK' & DataAll$State == 'PA'),])
WellsSorted = rbind(WellsSorted, DataAll[which(DataAll$WellDepth > 750 & DataAll$WellDepth < 1000 & DataAll$County == 'WARREN' & DataAll$State == 'PA' & DataAll$LongDgr > -79.4),])
WellsSorted = rbind(WellsSorted, DataAll[which(DataAll$WellDepth > 750 & DataAll$WellDepth < 1000 & DataAll$County == 'FOREST' & DataAll$State == 'PA'),])
WellsSorted = rbind(WellsSorted, DataAll[which(DataAll$WellDepth > 750 & DataAll$WellDepth < 1000 & DataAll$County == 'CLARION' & DataAll$State == 'PA' & DataAll$LatDegr > 41.3),])
WellsSorted = rbind(WellsSorted, DataAll[which(DataAll$WellDepth > 750 & DataAll$WellDepth < 1000 & DataAll$County == 'JEFFERSON' & DataAll$State == 'PA' & DataAll$LatDegr <= 41.16776),])

writeOGR(WellsSorted, dsn=getwd(), layer="WellsForOutlierTest_ESDA", driver="ESRI Shapefile")

dev.off()
png("WellsRemoved_PAWellsAddedBack.png", width=6, height=6, units='in', res=150)
plot(DataAll, pch=16, col='red')
plot(WellsSorted, pch=16, add=TRUE)
plot(NY, add=TRUE)
plot(PA, add=TRUE)
plot(WV, add=TRUE)
plot(MD, add=TRUE)
plot(KY, add=TRUE)
plot(VA, add=TRUE)
dev.off()

# Outlier identification ----
#Scripts for the actual outlier identification are from Whealton and Stedinger (2015)

#Load in the wells that were rerun and add them to the SortData$Sorted
Rerun = read.csv('SortedUniqueSpots_AllTemps_RerunAdded.csv', stringsAsFactors = FALSE)
Rerun = Rerun[nrow(Rerun),]
Rerun$APINo = SortData$RerunWells$APINo
SortData$Sorted@data[nrow(SortData$Sorted),] = Rerun[,-c(1,8,9,ncol(Rerun))] 

SortData$Sorted = spTransform(SortData$Sorted, CRS('+init=epsg:26917'))
#Add columns of UTM coordinates
SortData$Sorted$POINT_X = SortData$Sorted@coords[,1]
SortData$Sorted$POINT_Y = SortData$Sorted@coords[,2]

WellsSorted = spTransform(WellsSorted, CRS('+init=epsg:26917'))
#Add columns of UTM coordinates
WellsSorted$POINT_X = WellsSorted@coords[,1]
WellsSorted$POINT_Y = WellsSorted@coords[,2]

#Data must have a column of UTM coordinates in m for this to work because it relies on Euclidian distances.
#Must ensure that "test" and both coordinates are type num. 
#This should be run after duplicate points and negative gradient values have been taken care of. 

TestedOutliers_HeatFlow = select_out_algo(Data = WellsSorted@data, 
                                          InVarName = "Qs", OutVarName = "Qs", 
                                          X_coordName = "POINT_X", Y_coordName = "POINT_Y", 
                                          CensorName = 'WellDepth', 
                                          rad_max = 32000, 
                                          pt_min = 25, 
                                          k_loc = 3, 
                                          type = 7, 
                                          min_val = 0, max_val = 7000, rank = TRUE)

#Testing with a minimum depth of 2000 m. Specifically to see if deep points in the SWPA region are still outliers among deeper data.
TestedOutliers_HeatFlow_min2k = select_out_algo(Data = WellsSorted@data, 
                                          InVarName = "Qs", OutVarName = "Qs", 
                                          X_coordName = "POINT_X", Y_coordName = "POINT_Y", 
                                          CensorName = 'WellDepth', 
                                          rad_max = 32000, 
                                          pt_min = 25, 
                                          k_loc = 3, 
                                          type = 7, 
                                          min_val = 2000, max_val = 7000, rank = TRUE)

#Testing with a maximum depth of 2000 m. Specifically to see if deep points in the SWPA region are again outliers among shallower data.
TestedOutliers_HeatFlow_max2k = select_out_algo(Data = WellsSorted@data, 
                                          InVarName = "Qs", OutVarName = "Qs", 
                                          X_coordName = "POINT_X", Y_coordName = "POINT_Y", 
                                          CensorName = 'WellDepth', 
                                          rad_max = 32000, 
                                          pt_min = 25, 
                                          k_loc = 3, 
                                          type = 7, 
                                          min_val = 0, max_val = 2000, rank = TRUE)

#Convert to spatial data by adding the coords from WellsSorted.
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

writeOGR(TestedOutliers_HeatFlow$NotOutliers, dsn=getwd(), layer="DeepestWells_NotOutliers_32km_Qs_CorrBase_Ranked", driver = "ESRI Shapefile")
writeOGR(TestedOutliers_HeatFlow$Outliers, dsn=getwd(), layer="DeepestWells_Outliers_32km_Qs_CorrBase_Ranked", driver = "ESRI Shapefile")

Outs = spTransform(TestedOutliers_HeatFlow$Outliers, CRS = CRS("+init=epsg:4326"))
NoOuts = spTransform(TestedOutliers_HeatFlow$NotOutliers, CRS = CRS("+init=epsg:4326"))
AllData = rbind(NoOuts, Outs)

Outs_min2k = spTransform(TestedOutliers_HeatFlow_min2k$Outliers, CRS = CRS("+init=epsg:4326"))
NoOuts_min2k = spTransform(TestedOutliers_HeatFlow_min2k$NotOutliers, CRS = CRS("+init=epsg:4326"))

Outs_max2k = spTransform(TestedOutliers_HeatFlow_max2k$Outliers, CRS = CRS("+init=epsg:4326"))
NoOuts_max2k = spTransform(TestedOutliers_HeatFlow_max2k$NotOutliers, CRS = CRS("+init=epsg:4326"))

plot(NoOuts_max2k[which(NoOuts_max2k$WellDepth >= 2000 & NoOuts_max2k$WellDepth < 2400 & NoOuts_max2k$out_loc_error == 0),], pch = 16)
plot(Outs_max2k[which(Outs_max2k$WellDepth >= 2000 & Outs_max2k$WellDepth < 2400 & Outs_max2k$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue')
plot(Outs_max2k[which(Outs_max2k$WellDepth >= 2000 & Outs_max2k$WellDepth < 2400 & Outs_max2k$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red')
plot(NoOuts_max2k[which(NoOuts_max2k$WellDepth >= 2000 & NoOuts_max2k$WellDepth < 2400 & NoOuts_max2k$out_loc_error == 1),], pch = 17, col='purple',add=TRUE)
plot(NY, lwd = 2, add = TRUE)
plot(PA, lwd = 2, add = TRUE)
plot(WV, lwd = 2, add = TRUE)
plot(Counties, add = TRUE)

# Map of the depth ranks of outliers ----
scaleRange = c(1,25)
scaleBy = 5
colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)
Pal = c(Pal[-3], Pal[length(Pal)])
sets = rbind(c(1,3), c(2,4))
layout(sets)
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


# Depth Rank Empircal CDFs for KS Test ----
hist(Outs$out_loc_drank[which(Outs$out_loc_lo == 1)]*25, breaks = seq(0,25,1), main = 'Low Outliers Depth Rank', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=2, cex.axis = 1.5)
hist(Outs$out_loc_drank[which(Outs$out_loc_lo == 0)]*25, breaks = seq(0,25,1), main = 'High Outliers Depth Rank', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=2, cex.axis = 1.5)
hist(Outs$out_loc_drank*25, breaks = seq(0,25,1), main = 'All Outliers Depth Rank', xlab = 'Depth Rank', cex.lab = 1.5, cex.main=2, cex.axis = 1.5)

#Create an empirical CDF
CumLo = cumsum(hist(Outs$out_loc_drank[which(Outs$out_loc_lo == 1)]*25, breaks = seq(0.5,25,.5), plot = FALSE)$counts/sum(hist(Outs$out_loc_drank[which(Outs$out_loc_lo == 1)]*25, breaks = seq(0.5,25,.5), plot = FALSE)$counts))
CumHi = cumsum(hist(Outs$out_loc_drank[which(Outs$out_loc_lo == 0)]*25, breaks = seq(0.5,25,.5), plot = FALSE)$counts/sum(hist(Outs$out_loc_drank[which(Outs$out_loc_lo == 0)]*25, breaks = seq(0.5,25,.5), plot = FALSE)$counts))

#Make Step Function
LoFun = stepfun(x = seq(1,25,.5), y = c(0,CumLo))
HiFun = stepfun(x = seq(1,25,.5), y = c(0,CumHi))
UnifFun = stepfun(x = seq(1,25,.5), y = c(0,seq(1,25,.5)/25))

par(xaxs = 'i', yaxs = 'i')
plot(LoFun, pch = NA, col = 'blue', lwd = 2, xlab = 'Depth Rank', ylab = 'Cumulative Frequency', xlim = c(0,25), ylim = c(0,1), cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, main = 'Empirical CDFs for K-S Test')
par(new = TRUE)
plot(HiFun, pch = NA, col = 'red', axes=FALSE, xlab='', ylab='', lwd = 2, xlim = c(0,25), ylim = c(0,1), main = '')
par(new = TRUE)
plot(x = sort(rep(seq(1,25,.5),2))[2:length(sort(rep(seq(1,25,.5),2)))], y = sort(rep(seq(1,25,.5)/25,2))[-length(sort(rep(seq(1,25,.5)/25,2)))], col = 'black', type = 'l', axes=FALSE, xlab='', ylab='', lwd = 2, xlim = c(0,25), ylim = c(0,1), main = '')
legend('topleft', legend = c('Discrete Uniform Distribution', 'Low Outliers', 'High Outliers'), lwd = 2, pch = NA, col = c('black', 'blue', 'red'))

#KS Test
max(abs(sort(rep(CumHi,2))[-length(CumHi)*2] - sort(rep(seq(1,25,.5)/25,2))[-length(sort(rep(seq(1,25,.5)/25,2)))]))
max(abs(sort(rep(CumLo,2))[-length(CumLo)*2] - sort(rep(seq(1,25,.5)/25,2))[-length(sort(rep(seq(1,25,.5)/25,2)))]))
KSHi = ks.test(x = Outs$out_loc_drank[which(Outs$out_loc_lo == 0)]*25, y = UnifFun, simulate.p.value = TRUE, B = 10000)
KSLo = ks.test(x = Outs$out_loc_drank[which(Outs$out_loc_lo == 1)]*25, y = UnifFun, simulate.p.value = TRUE, B = 10000)

#Both significant at 1% level. The distribution is not random at 1% level.

# Poisson Test for Depth Bins ----
#Arrival rate in outliers per well
#lambda = nrow(Outs)/nrow(AllData) 

#Outliers in each bin class
#OutBins = hist(Outs$WellDepth, breaks = c(500,seq(1000,2000,200),3000,7000), plot = FALSE)$counts
#Total data points in each bin class
#DatBins = hist(AllData$WellDepth, breaks = c(500,seq(1000,2000,200),3000,7000), plot = FALSE)$counts

#Poisson probabilities for each bin class
#nu = lambda*DatBins
#dpois(OutBins, nu)

# Binomial and Chi-Squared Test for Depth Bins ----
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
ChiSq_LowHigh = sum((OutBins[c(9,10,11,12,13)] - PerfBins[c(9,10,11,12,13)])^2/PerfBins[c(9,10,11,12,13)])
dfChiSq_LowHigh = length(OutBins[c(9,10,11,12,13)]) - 1
pVal_LowHigh = 1-pchisq(ChiSq_LowHigh, dfChiSq_LowHigh)

#Lower tail test less than 1400 m - Not Significant p = 25%
ChiSq_LowHigh = sum((OutBins[c(1,2,3)] - PerfBins[c(1,2,3)])^2/PerfBins[c(1,2,3)])
dfChiSq_LowHigh = length(OutBins[c(1,2,3)]) - 1
pVal_LowHigh = 1-pchisq(ChiSq_LowHigh, dfChiSq_LowHigh)

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

#Bin 7 and 8 test significant at .01% level
ChiSq_LowHigh = sum((OutBins[c(7,8)] - PerfBins[c(7,8)])^2/PerfBins[c(7,8)])
dfChiSq_LowHigh = length(OutBins[c(7,8)]) - 1
pVal_LowHigh = 1-pchisq(ChiSq_LowHigh, dfChiSq_LowHigh)

CramerV = sqrt(ChiSq_LowHigh/length(which(AllData$out_loc_error == 0 & AllData$WellDepth >= 2000 & AllData$WellDepth < 2400))/min(1,1))
BiasCorrV = sqrt(max(0, ChiSq_LowHigh/length(which(AllData$out_loc_error == 0 & AllData$WellDepth >= 2000 & AllData$WellDepth < 2400)) - 1/(length(which(AllData$out_loc_error == 0 & AllData$WellDepth >= 2000 & AllData$WellDepth < 2400)) - 1))/
                   min(2 - 1/(length(which(AllData$out_loc_error == 0 & AllData$WellDepth >= 2000 & AllData$WellDepth < 2400)) - 1) - 1, 2 - 1/(length(which(AllData$out_loc_error == 0 & AllData$WellDepth >= 2000 & AllData$WellDepth < 2400)) - 1) - 1))


# Maps of the depth slices for outliers ----

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

#Are any ranks clustered more than would be expected under a random process?


# Make panel plot of depth horizons with low and high outliers colored in, as well as missing datapoints
sets = rbind(c(1,2,3,4),c(5,6,7,8))
png('OutsByDepth_Panels.png', res = 300, height = 8, width = 16, units = 'in')
layout(sets)
par(xaxs = 'i', yaxs = 'i', mar = c(2,2.5,1,1))
#Data from 750 - 1000 - 4 low outliers that are not clustered. Not problematic.
plot(AllData, col = 'white')
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE)
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NoOuts[which(NoOuts$WellDepth >= 750 & NoOuts$WellDepth < 1000 & NoOuts$out_loc_error == 0),], pch = 16, add = TRUE)
plot(Outs[which(Outs$WellDepth >= 750 & Outs$WellDepth < 1000 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue')
plot(Outs[which(Outs$WellDepth >= 750 & Outs$WellDepth < 1000 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red')
plot(NoOuts[which(NoOuts$WellDepth >= 750 & NoOuts$WellDepth < 1000 & NoOuts$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple')
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
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE)
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NoOuts[which(NoOuts$WellDepth >= 1000 & NoOuts$WellDepth < 1200 & NoOuts$out_loc_error == 0),], pch = 16, add = TRUE)
plot(Outs[which(Outs$WellDepth >= 1000 & Outs$WellDepth < 1200 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue')
plot(Outs[which(Outs$WellDepth >= 1000 & Outs$WellDepth < 1200 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red')
plot(NoOuts[which(NoOuts$WellDepth >= 1000 & NoOuts$WellDepth < 1200 & NoOuts$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple')
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
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE)
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NoOuts[which(NoOuts$WellDepth >= 1200 & NoOuts$WellDepth < 1400 & NoOuts$out_loc_error == 0),], pch = 16, add = TRUE)
plot(Outs[which(Outs$WellDepth >= 1200 & Outs$WellDepth < 1400 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue')
plot(Outs[which(Outs$WellDepth >= 1200 & Outs$WellDepth < 1400 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red')
plot(NoOuts[which(NoOuts$WellDepth >= 1200 & NoOuts$WellDepth < 1400 & NoOuts$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple')
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
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE)
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NoOuts[which(NoOuts$WellDepth >= 1400 & NoOuts$WellDepth < 1600 & NoOuts$out_loc_error == 0),], pch = 16, add = TRUE)
plot(Outs[which(Outs$WellDepth >= 1400 & Outs$WellDepth < 1600 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue')
plot(Outs[which(Outs$WellDepth >= 1400 & Outs$WellDepth < 1600 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red')
plot(NoOuts[which(NoOuts$WellDepth >= 1400 & NoOuts$WellDepth < 1600 & NoOuts$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple')
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
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE)
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NoOuts[which(NoOuts$WellDepth >= 1600 & NoOuts$WellDepth < 2000 & NoOuts$out_loc_error == 0),], pch = 16, add = TRUE)
plot(Outs[which(Outs$WellDepth >= 1600 & Outs$WellDepth < 2000 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue')
plot(Outs[which(Outs$WellDepth >= 1600 & Outs$WellDepth < 2000 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red')
plot(NoOuts[which(NoOuts$WellDepth >= 1600 & NoOuts$WellDepth < 2000 & NoOuts$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple')
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
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE)
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NoOuts[which(NoOuts$WellDepth >= 2000 & NoOuts$WellDepth < 2400 & NoOuts$out_loc_error == 0),], pch = 16, add = TRUE)
plot(Outs[which(Outs$WellDepth >= 2000 & Outs$WellDepth < 2400 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue')
plot(Outs[which(Outs$WellDepth >= 2000 & Outs$WellDepth < 2400 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red')
plot(NoOuts[which(NoOuts$WellDepth >= 2000 & NoOuts$WellDepth < 2400 & NoOuts$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple')
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
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE)
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NoOuts[which(NoOuts$WellDepth >= 2400 & NoOuts$WellDepth < 3000 & NoOuts$out_loc_error == 0),], pch = 16, add = TRUE)
plot(Outs[which(Outs$WellDepth >= 2400 & Outs$WellDepth < 3000 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue')
plot(Outs[which(Outs$WellDepth >= 2400 & Outs$WellDepth < 3000 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red')
plot(NoOuts[which(NoOuts$WellDepth >= 2400 & NoOuts$WellDepth < 3000 & NoOuts$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple')
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
plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], lwd = 1, add = TRUE)
plot(NY, add = TRUE, lwd = 2)
plot(PA, add = TRUE, lwd = 2)
plot(WV, add = TRUE, lwd = 2)
plot(MD, lwd = 2, add=TRUE)
plot(KY, lwd = 2, add=TRUE)
plot(VA, lwd = 2, add=TRUE)
plot(NoOuts[which(NoOuts$WellDepth >= 3000 & NoOuts$WellDepth < 6600 & NoOuts$out_loc_error == 0),], pch = 16, add = TRUE)
plot(Outs[which(Outs$WellDepth >= 3000 & Outs$WellDepth < 6600 & Outs$out_loc_lo == 1),], pch = 16, add = TRUE, col = 'blue')
plot(Outs[which(Outs$WellDepth >= 3000 & Outs$WellDepth < 6600 & Outs$out_loc_lo == 0),], pch = 16, add = TRUE, col = 'red')
plot(NoOuts[which(NoOuts$WellDepth >= 3000 & NoOuts$WellDepth < 6600 & NoOuts$out_loc_error == 1),], pch = 17, add = TRUE, col = 'purple')
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


# Q-Q plot for the points that are not outliers ----
##Well point file load, spatially referenced. Needs to be in NAD83 UTM17N coordinates
CT = readOGR(dsn=paste(getwd(), '/BaseCorrWells', sep=''), layer="WellsCT")
CNY = readOGR(dsn=paste(getwd(), '/BaseCorrWells', sep=''), layer="WellsCNY")
CWV = readOGR(dsn=paste(getwd(), '/BaseCorrWells', sep=''), layer="WellsCWV")
ENY = readOGR(dsn=paste(getwd(), '/BaseCorrWells', sep=''), layer="WellsENY")
ENYPA = readOGR(dsn=paste(getwd(), '/BaseCorrWells', sep=''), layer="WellsENYPA")
MT = readOGR(dsn=paste(getwd(), '/BaseCorrWells', sep=''), layer="WellsMT")
NWPANY = readOGR(dsn=paste(getwd(), '/BaseCorrWells', sep=''), layer="WellsNWPANY")  
SWPA = readOGR(dsn=paste(getwd(), '/BaseCorrWells', sep=''), layer="WellsSWPA")
WPA = readOGR(dsn=paste(getwd(), '/BaseCorrWells', sep=''), layer="WellsWPA")
FL = readOGR(dsn=paste(getwd(), '/BaseCorrWells',sep=''), layer="WellsAll")

#Q-Q plots for data in each interpolation region
#Unique axes
sets = rbind(c(1,2,3), c(4,5,6),c(7,8,9))
png('QQPlotHeatFlow_BeforeCorr_NotOutTest.png', res=600, units='in', width=10, height=10)
layout(sets)
par(mar=c(4,5,3,2))
qqnorm(CT@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Chautauqua, NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='red')
qqline(CT@data$Qs)
qqnorm(WPA@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='orange', xlim = c(-3.5,3.5), ylim = c(10,70))
qqline(WPA@data$Qs)
temp = qqnorm(WPA@data$Qs, plot.it = FALSE)$x[which(WPA@data$ot_lc_rr == 1)]
par(new = TRUE)
plot(x = temp, y = WPA@data$Qs[which(WPA@data$ot_lc_rr == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(10,70))
qqnorm(NWPANY@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Northwestern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='yellow', xlim = c(-3,3), ylim = c(5,75))
qqline(NWPANY@data$Qs)
temp = qqnorm(NWPANY@data$Qs, plot.it = FALSE)$x[which(NWPANY@data$ot_lc_rr == 1)]
par(new = TRUE)
plot(x = temp, y = NWPANY@data$Qs[which(NWPANY@data$ot_lc_rr == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(5,75))
qqnorm(CNY@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='green', xlim = c(-2.75,2.75), ylim = c(35,65))
qqline(CNY@data$Qs)
temp = qqnorm(CNY@data$Qs, plot.it = FALSE)$x[which(CNY@data$ot_lc_rr == 1)]
par(new = TRUE)
plot(x = temp, y = CNY@data$Qs[which(CNY@data$ot_lc_rr == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-2.75,2.75), ylim = c(35,65))
qqnorm(ENY@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='springgreen', xlim = c(-3,3), ylim = c(30,75))
qqline(ENY@data$Qs)
temp = qqnorm(ENY@data$Qs, plot.it = FALSE)$x[which(ENY@data$ot_lc_rr == 1)]
par(new = TRUE)
plot(x = temp, y = ENY@data$Qs[which(ENY@data$ot_lc_rr == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(30,75))
qqnorm(ENYPA@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Eastern NY and PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='skyblue', xlim = c(-3,3), ylim = c(10,110))
qqline(ENYPA@data$Qs)
temp = qqnorm(ENYPA@data$Qs, plot.it = FALSE)$x[which(ENYPA@data$ot_lc_rr == 1)]
par(new = TRUE)
plot(x = temp, y = ENYPA@data$Qs[which(ENYPA@data$ot_lc_rr == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3,3), ylim = c(10,110))
qqnorm(SWPA@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Southwestern PA', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='blue', xlim = c(-4,4), ylim = c(5,110))
qqline(SWPA@data$Qs)
temp = qqnorm(SWPA@data$Qs, plot.it = FALSE)$x[which(SWPA@data$ot_lc_rr == 1)]
par(new = TRUE)
plot(x = temp, y = SWPA@data$Qs[which(SWPA@data$ot_lc_rr == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,110))
qqnorm(CWV@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Central WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='purple', xlim = c(-4,4), ylim = c(5,120))
qqline(CWV@data$Qs)
temp = qqnorm(CWV@data$Qs, plot.it = FALSE)$x[which(CWV@data$ot_lc_rr == 1)]
par(new = TRUE)
plot(x = temp, y = CWV@data$Qs[which(CWV@data$ot_lc_rr == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-4,4), ylim = c(5,120))
qqnorm(MT@data$Qs, cex.axis=1.5, cex.lab=1.5, main='Western WV', ylab=expression('Sample Quantiles' ~ (mW/m^2)), cex.main=2, col='violet', xlim = c(-3.5,3.5), ylim = c(5,120))
qqline(MT@data$Qs)
temp = qqnorm(MT@data$Qs, plot.it = FALSE)$x[which(MT@data$ot_lc_rr == 1)]
par(new = TRUE)
plot(x = temp, y = MT@data$Qs[which(MT@data$ot_lc_rr == 1)], col = 'black', pch = 16, xlab = '', ylab = '', axes = FALSE, xlim = c(-3.5,3.5), ylim = c(5,120))
dev.off()
