#This function is used to find the deviation of each heat flow point from its local region of points.
#The regional mean and the regional median are both calculated for each point based on the 
#neighborhood criteria. These are the radius size and the maximum number of points.

QsDev = function(Data,    # The unprojected dataset
                 Var,     # The variable to be tested for the local region
                 xName,   # Name of the x coordinate in Data in UTM coordinates
                 yName,   # Name of the y coordinate in Data in UTM coordinates
                 rad,     # Radius to search for points (m)
                 max_pts  # maximum number of points
){
  # Add a column for the local mean and median
  Data$RegMean = 0
  Data$RegMed = 0
  # Add a column for the number of points in the neighborhood
  Data$RegPts = 0
  # Add a column for the distance to the 2nd point.
  Data$Dist2 = 0
  
  # Find the distance from the point to all other points in the dataset
  for(i in 1:nrow(Data)){
    
    # initializing matrix with one column for index and the other for distance
    dist_vec <- matrix(0,nrow(Data),1) 
    
    # calculating weighted Euclidian distance for the points
    dist_vec[,1] <- sqrt((Data[,xName] - Data[i,xName])^2 + (Data[,yName] - Data[i,yName])^2)
    #Remove the point being tested from the dist_vec
    dist_vec = dist_vec[-i,1]
    
    # finding the distance corresponding to the maximum number of points
    dist_cutoff <- sort(dist_vec)[max_pts]
    # Find distance to the second point
    Data$Dist2[i] = sort(dist_vec)[2]
    
    # finding the neighbors within rad. Smaller indices are closer to the point
    if (dist_cutoff <= rad){
      Neighbors <- which(dist_vec <= dist_cutoff)
    }else{
      Neighbors <- which(dist_vec <= rad)
    }
    
    if (length(Neighbors) <= 1){
      #There are too few neighbors within rad
      Data$RegPts[i] = length(Neighbors)
      Data$RegMean[i] = -9999
      Data$RegMed[i] = -9999
    }else{
      # Adding the number of neighbors to the dataframe
      Data$RegPts[i] = length(Neighbors)
      
      #Calculate the mean and median of these points, and add it to the database
      Data$RegMean[i] = mean(Data[Neighbors, Var])
      Data$RegMed[i] = median(Data[Neighbors, Var])
    }
  }
  return(Data)
}