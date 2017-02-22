#### Sensitivity Analysis for Local Points Algorithm ####
#Based on the "sensitivity analysis for local points algorithm" section of example_outlier_code.R in geothermal_pfa repository

library(rgdal)
# Load data----
DataTest = readOGR(dsn = "C:\\Users\\Jared\\Documents\\Cornell\\Research\\Masters - Spatial Assessment\\Figures", layer = "DataForOutlierTesting", stringsAsFactors = FALSE)

# Load Outlier Identification Functions----
source("C:\\Users\\Jared\\Documents\\Cornell\\Research\\Publications\\DOE Grant\\CombiningRiskFactorCode\\geothermal_pfa\\outliers\\outlier_identification.R")

# Set neighborhood parameters----
pts_sens <- c(10,25,50,100,200) # number of points for local neighborhood
rad_sens <- c(4000,8000,16000,32000,64000) # maximum size of radius

outs_iden <- matrix(0,length(pts_sens),length(rad_sens)) # matrix to hold number of outliers
outs_sparse <- matrix(0,length(pts_sens),length(rad_sens)) # number of points in sparse areas

# Calculate Outliers----
for(i in 1:length(pts_sens)){
  for(j in 1:length(rad_sens)){
    sens_data <- DataTest
    sens_data2 <- select_out_algo(Data = sens_data@data,
                                  OutVarName = "Qs",
                                  InVarName = "Qs",
                                  X_coordName = "POINT_X",
                                  Y_coordName = "POINT_Y",
                                  Threshold = 0.0,
                                  algo = 1,
                                  outcri = 1,
                                  pt_eval = pts_sens[i],
                                  rad_eval = 16000,
                                  box_size = 32000,
                                  pt_min = 25,
                                  rad_max = rad_sens[j],
                                  k_glob = 3,
                                  k_loc = 3,
                                  type = 7)
    
    outs_iden[i,j] <- nrow(sens_data2$Outliers)
    outs_sparse[i,j] <- sum(sens_data2$NotOutliers$out_loc_error)
  }
}
rm(sens_data, sens_data2, i, j, outs_iden, outs_sparse)

# Creating plot----
#File type variable. 1 = EPS, 0 = EMF
EPS = 1
setwd("C:\\Users\\Jared\\Documents\\Cornell\\Research\\Publications\\ESDA")
if (EPS == 1){
  setEPS()
  postscript(file = "outlier_sens.eps", title = "Sensitivity Outliers Local Points", width = 5, height = 5)
}else{
  #As EMF File for Word
  emf(file = "outlier_sens.emf", width = 5, height = 5)
}

#Make color ramp
Pal = colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))
cols <- Pal(max(outs_iden)+1)
par(mar =c(3,3,0,8)+0.1) 

data <- matrix(0,length(pts_sens)*length(rad_sens),4)

for(i in 1:length(rad_sens)){
  inds <- seq((i-1)*length(rad_sens)+1,i*length(rad_sens), by=1)
  
  data[c(inds),1] <- outs_iden[,i]
  data[c(inds),2] <- outs_sparse[,i]
  data[c(inds),3] <- pts_sens
  data[c(inds),4] <- rad_sens[i]/1000
}
rm(i, inds)
dataplot <- data.frame(data)
rm(data)
colnames(dataplot) <- c("iden", "sparse", "pts", "rad")

#Changing the plotting location for pts = 10 to 12.5 for equal spacing in log base 2
dataplot$pts[dataplot$pts == 10] = 12.5

#Assigning color and size of points
dataplot$cols <- cols[dataplot$iden+1]
dataplot$cex <- sqrt(nrow(DataTest) - dataplot$sparse)/30

plot(dataplot$pts
     , log(dataplot$rad, base=2)
     , log='x'
     , cex = 1.02*sqrt(nrow(DataTest))/30
     , col = "black"
     , pch = 19
     , ylab = "Maximum Neighborhood Radius (km)"
     , xlab = "Points Needed to Evaluate"
     , xaxt = "n"
     , yaxt = "n"
     , ylim = c(1.9,6.1)
     , xlim = c(10.85,230)
     , line = 2
)
points(dataplot$pts
       , log(dataplot$rad, base=2)
       , cex = 0.98*sqrt(nrow(DataTest))/30
       , col = "white"
       , pch = 19
)
points(dataplot$pts # x value
       , log(dataplot$rad, base=2)   # y value
       , cex = dataplot$cex
       , col = dataplot$cols
       , pch = 19
)

axis(1, at=c(12.5, pts_sens[-1]), labels=pts_sens, padj = -0.5)
axis(2, at=seq(2,6,1), labels=rad_sens/1000, padj = 0.5)

par(xpd = TRUE)
legend(x = 300
       , y = 6.25 # location
       , title = "% Outliers in Data"
       , legend = seq(0,8,1) # legend entries
       , pch = 19
       , col = cols[round(max(outs_iden)/(max(outs_iden)/nrow(DataTest))*seq(0,0.08,0.01),0) + 1]
       , ncol = 1
)
legend(x = 300
       , y = 3.9 # location
       , title = expression(paste(phantom("bl"), "% Data Tested", phantom("a")))
       , legend = c("20", "40", "60", "80", "100")
       , pch = 19
       , col = "gray"
       , ncol = 1
       , pt.cex = c(sqrt(nrow(DataTest)/5)/30, sqrt(nrow(DataTest)*2/5)/30, sqrt(nrow(DataTest)*3/5)/30, sqrt(nrow(DataTest)*4/5)/30, sqrt(nrow(DataTest))/30)
       , y.intersp = 1.6
       , x.intersp = 1.3
)
legend(x = 300
       , y = 3.9 # location
       , title = expression(paste(phantom("bl"), phantom("% Data Tested"), phantom("a")))
       , legend = c(expression(phantom("20")), expression(phantom("40")), expression(phantom("60")), expression(phantom("80")), expression(phantom("100")))
       , pch = 1
       , col = 'black'
       , ncol = 1
       , pt.cex = 1.02*rep(sqrt(nrow(DataTest))/30, 5)
       , bty='n'
       , y.intersp = 1.6
       , x.intersp = 1.3
       , pt.lwd = 0.4
)
dev.off()


# Global outlier test on the dataset ----
sens_data <- DataTest
GlobOut <- select_out_algo(Data = sens_data,
                             OutVarName = "Qs",
                             InVarName = "Qs",
                             X_coordName = "POINT_X",
                             Y_coordName = "POINT_Y",
                             Threshold = 0.0,
                             algo = 1,
                             outcri = 3,
                             pt_eval = 1000,
                             rad_eval = 16000,
                             box_size = 32000,
                             pt_min = 25,
                             rad_max = 200000,
                             k_glob = 3,
                             k_loc = 3,
                             type = 7)
sens_data <- DataTest
LocOut_ManyPts <- select_out_algo(Data = sens_data,
                           OutVarName = "Qs",
                           InVarName = "Qs",
                           X_coordName = "POINT_X",
                           Y_coordName = "POINT_Y",
                           Threshold = 0.0,
                           algo = 1,
                           outcri = 1,
                           pt_eval = 1000,
                           rad_eval = 16000,
                           box_size = 32000,
                           pt_min = 25,
                           rad_max = 200000,
                           k_glob = 3,
                           k_loc = 3,
                           type = 7)
rm(sens_data)