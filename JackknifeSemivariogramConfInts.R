#Function for jackknife computation of confidence intervals for semivariogram lags
#The serial version accepts anisotropy parameters. These features need to be worked into the parallel version still.

#Note that if the number of point pairs in a bin is 2, then the standard deviation will be 0 
# because when only 2 samples are removed in that bin there is sample reuse because each point contributes the same value when it is removed. Therefore the jackknife variance is likely too small when the number of data points is small.
#Also note that point pairs is not the same as number of points.

JackKnife = function(Dat,      #Spatial points datafame containing the points to be jackknifed. Should be in UTM coordinates.
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
    #Make a matrix for storing the bin estimates of variance, number of points, and the distance to points.
    VarioMat = Pts = Dist = Gamma = matrix(NA, nrow=nrow(Dat), ncol=bins)

    #Make a vector of the number of replications of each bin from the jackknife. This number only changes if a given resampling results in no estimate for a given bin.
    Nvec = rep(nrow(Dat), bins)

    #Start jackknife
    for (i in 1:nrow(Dat)){
      vj = variogram(Qs~1, data = Dat[-i,], cutoff=cut, width=cut/bins)

      #Store the number of points in each lag and the average lag distance.
      Pts[i,] = vj$np
      Dist[i,] = vj$dist
      Gamma[i,] = vj$gamma

      if (any(is.na(Pts[i,])) || any(Pts[i,] < 2)){
        #Mark the distance as NA and reduce the number of points in Nvec. The model will have to be rerun with the new Nvec...
        ind = which((is.na(Pts[i,]) == TRUE) | (Pts[i,] < 2))
        Pts[i, ind] = NA
        Dist[i, ind] = NA
        Nvec[ind] = Nvec[ind] - 1
        Gamma[i, ind] = NA
      }
    }
    rm(i)

    #Calculate the estimate of the bin mean from the jackknife. The transposition is a faster way of using sweep() to multiply a vector by row of a matrix.
    VarioMat =  t(Nvec*v.Dat$gamma - t(t(t(Gamma)*(Nvec-1))))

    #Calculate the Jackknife mean. Remove NAs generated from before.
    BinMean = apply(VarioMat, 2, FUN=sum, na.rm=TRUE)/Nvec

    #Calculate the average bin distance for plotting purposes
    BinMean_dist = apply(Dist, 2, FUN=sum, na.rm=TRUE)/Nvec

    #Calculate the jackknife standard error
    VarioMatSquaredMatrix = matrix(NA, ncol=bins, nrow=nrow(Dat))
    DistSd = matrix(NA, ncol=bins, nrow=nrow(Dat))
    for (j in 1:nrow(Dat)){
      VarioMatSquaredMatrix[j, ] = (VarioMat[j, ] - BinMean)^2
      DistSd[j, ] = (Dist[j, ] - BinMean_dist)^2
    }
    rm(j)

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
