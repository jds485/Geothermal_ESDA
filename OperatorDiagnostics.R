#ESDA procedure for finding possibly bad operators

BadOperatorDiagnostics = function(MT, #Spatial dataframe containing a column named "Operator" 
                                  MT_WGS, #Spatial dataframe in WGS coordinates to make map
                                  v.MT, #variogram for MT dataset using all data
                                  rv.MT, #Cressie's robust variogram for MT dataset using all data
                                  Vcut, #variogram cutoff in m
                                  Vbins, #Number of variogram bins
                                  o,  #index for the operator to be left out
                                  LowLim = 2, #Lower limit for number of wells drilled by an operator. Need at least LowLim wells to make a plot for the operator.
                                  HistSep = 10, #x-axis bar separartion on histogram
                                  Histylim = 300, #y-axis upper limit on histogram
                                  RegName, #Region name for figure
                                  res = 300, #png resolution
                                  height = 8, #png height in inches
                                  width = 8, #png width in inches
                                  MaxLagDist, #maximum lag separation distance for which to compute jackknife distance metrics
                                  SensitivityDist = NA, #If supplied, will provide results of the metrics in the specified variogram lag distance interval for MaxLagDist
                                  plt = TRUE # should diagnostic plots be made for each operator?
){
  #Unique operators
  UniOps = unique(MT$Operator)
  
  #Only plot the operator if they have more than LowLim well. It doesn't make sense otherwise.
  if(nrow(MT[MT$Operator == UniOps[o],]) >= LowLim){
    #Operator left-out variograms
    v.o = variogram(Qs~1, MT[-which(MT$Operator == UniOps[o]),], cutoff=Vcut, width = Vcut/Vbins)
    rv.o = variogram(Qs~1, MT[-which(MT$Operator == UniOps[o]),], cutoff=Vcut, width = Vcut/Vbins, cressie = TRUE)
    
    if (plt){
      png(paste0(RegName, '_OperatorESDA_', o, '.png'), res = res, height = height, width = width, units = 'in')
      #Histograms Overlain
      #layout(rbind(c(1,4), c(2,3)))
      #Histogram showing difference between current operator (red) and all data (black)
      #hist(MT$Qs, col = 'black', main = UniOps[o], xlab = expression(paste('Heat Flow (mW/m'^2,')')), ylab = 'Frequency', ylim = c(0,Histylim), xlim = c(0,round(max(MT$Qs + HistSep/2),-1)), breaks = seq(round(min(MT$Qs - HistSep/2),-1),round(max(MT$Qs + HistSep/2),-1),HistSep))
      #par(new=TRUE)
      #hist(MT$Qs[-which(MT$Operator == UniOps[o])], border = 'red', axes = FALSE, xlab = '', main = '', ylab = '', ylim = c(0,Histylim), xlim = c(0,round(max(MT$Qs + HistSep/2),-1)), breaks = seq(round(min(MT$Qs - HistSep/2),-1),round(max(MT$Qs + HistSep/2),-1),HistSep))
      #legend('topright', legend = c('All data', 'Operator Removed'), col = c('black', 'red'), pch = 15)
      
      #Histograms on separate plots
      #layout(rbind(c(1,5), c(2,5), c(3,4), c(3,4)))
      #hist(MT$Qs, col = 'black', main = UniOps[o], xlab = expression(paste('Heat Flow (mW/m'^2,')')), ylab = 'Frequency', ylim = c(0,Histylim), xlim = c(0,round(max(MT$Qs + HistSep/2),-1)), breaks = seq(round(min(MT$Qs - HistSep/2),-1),round(max(MT$Qs + HistSep/2),-1),HistSep))
      #par(new=TRUE)
      #hist(MT$Qs[-which(MT$Operator == UniOps[o])], border = 'red', axes = FALSE, xlab = '', main = '', ylab = '', ylim = c(0,Histylim), xlim = c(0,round(max(MT$Qs + HistSep/2),-1)), breaks = seq(round(min(MT$Qs - HistSep/2),-1),round(max(MT$Qs + HistSep/2),-1),HistSep))
      #legend('topright', legend = c('All data', 'Operator Removed'), col = c('black', 'red'), pch = 15)
      
      #Boxplot
      layout(rbind(c(1,4), c(2,3)))
      boxplot(horizontal = TRUE, x = MT$Qs, main = UniOps[o], xlab = expression(paste('Heat Flow (mW/m'^2,')')), ylim = c(0,round(max(MT$Qs + HistSep/2),-1)), at = 1, xlim = c(0,2))
      par(new = TRUE)
      #All red
      #plot(x = MT$Qs[which(MT$Operator == UniOps[o])], y = jitter(rep(1, length(MT$Qs[which(MT$Operator == UniOps[o])])), amount = 0.4), col = 'red', xlab = '', ylab = '', axes = FALSE, xlim = c(0,round(max(MT$Qs + HistSep/2),-1)), ylim = c(0,2))
      #Colored by value of heat flow
      plot(x = MT$Qs[which(MT$Operator == UniOps[o])], y = jitter(rep(1, length(MT$Qs[which(MT$Operator == UniOps[o])])), amount = 0.4), col = colFun(MT$Qs[which(MT$Operator == UniOps[o])]), xlab = '', ylab = '', axes = FALSE, xlim = c(0,round(max(MT$Qs + HistSep/2),-1)), ylim = c(0,2))
      legend('topright', legend = c(paste("All Wells: N wells =", nrow(MT)), paste("Operator's Wells: N wells =", nrow(MT[which(MT$Operator == UniOps[o]),]))), col = c('black', 'red'), lty = c(1,NA), pch = c(NA,1))
      
      #Variogram of region with and without operator
      par(mar = c(5,5,2,1))
      plot(v.MT$dist, v.MT$gamma, ylim = c(0,round(max(v.MT$gamma),-1)), xlim = c(0,Vcut), xlab = 'Separation Distance (m)', ylab = expression(paste('Semivariance (mW/m'^2,')'^2)), main = 'MOM Semi-variogram')
      par(new=TRUE)
      plot(v.o$dist, v.o$gamma, ylim = c(0,round(max(v.MT$gamma),-1)), xlim = c(0,Vcut), xlab = '', ylab = '', col = 'red')
      legend('bottomright', legend = c('All Wells', 'Operator Removed'), col = c('black', 'red'), pch = 1)
      
      #Robust variogram of region with and without operator
      plot(rv.MT$dist, rv.MT$gamma, ylim = c(0,round(max(v.MT$gamma),-1)), xlim = c(0,Vcut), xlab = 'Separation Distance (m)', ylab = expression(paste('Semivariance (mW/m'^2,')'^2)), main = 'Robust Semi-variogram')
      par(new=TRUE)
      plot(rv.o$dist, rv.o$gamma, ylim = c(0,round(max(v.MT$gamma),-1)), xlim = c(0,Vcut), xlab = '', ylab = '', col = 'red')
      legend('bottomright', legend = c('All Wells', 'Operator Removed'), col = c('black', 'red'), pch = 1)
      
      #Map
      plot(MT_WGS, pch = 16, cex = 0.4, col ='white')
      plot(Counties[which(Counties$STATEFP == 42 | Counties$STATEFP == 36 | Counties$STATEFP == 54 | Counties$STATEFP == 51| Counties$STATEFP == 24| Counties$STATEFP == 21),], add=TRUE, border = 'grey')
      plot(NY, lwd = 2, add=TRUE)
      plot(PA, lwd = 2, add=TRUE)
      plot(WV, lwd = 2, add=TRUE)
      plot(MD, lwd = 2, add=TRUE)
      plot(KY, lwd = 2, add=TRUE)
      plot(VA, lwd = 2, add=TRUE)
      north.arrow(-83.5, 37.8, 0.05, lab = 'N', cex.lab = 1.5, col='black', cex = 0.7)
      degAxis(side = 2, seq(34, 46, 1), cex.axis = 1.5)
      degAxis(side = 2, seq(34, 46, 1), labels = FALSE)
      degAxis(side = 4, seq(34, 46, 1), labels = FALSE)
      degAxis(side = 1, seq(-70, -86, -1), cex.axis = 1.5)
      degAxis(side = 3, seq(-70, -86, -1), labels = FALSE)
      degAxis(side = 1, seq(-70, -86, -1), labels = FALSE)
      plot(MT_WGS, pch = 16, cex = 0.4, add = T)
      #All red
      #plot(MT_WGS[which(MT$Operator == UniOps[o]),], pch = 16, cex = 0.4, add = T, col = 'red')
      #Colored by heat flow value. Deeper on top of shallower where there are overlaps in space.
      plot(MT_WGS[which(MT$Operator == UniOps[o]),][order(MT_WGS[which(MT$Operator == UniOps[o]),]$WellDepth, decreasing = FALSE),], pch = 16, cex = 0.4, add = T, col = colFun(MT_WGS$Qs[which(MT_WGS$Operator == UniOps[o])][order(MT_WGS[which(MT$Operator == UniOps[o]),]$WellDepth, decreasing = FALSE)]))
      legend('topleft', legend = c('Other Wells', "Operator's Wells"), col = c('black', 'red'), pch = 16, cex = 1)
      dev.off()
    }
    
    #Number of operator wells
    Nop = length(which(MT$Operator == UniOps[o]))
    
    #Check if a sensitivity analysis on maximum variogram lag distance is requested.
    if (!is.na(SensitivityDist)){
      #Run sensitivity analysis
      #Determine the distances to evaluate
      Dists = seq(SensitivityDist, MaxLagDist, SensitivityDist)
      
      #Compute metrics for each of the distances
      #Get the lag that should be selected for cumulative sums
      DistLags = DiffMOMW = DiffRobustW = vector('numeric', length(Dists))
      for (d in 1:length(Dists)){
        DistLags[d] = max(which(v.MT$dist <= Dists[d]))
        
        #Weighting by the separation distance lags
        #Note: not using different weights for MOM and Robust because want to compare the same number of lags for both
        weights = v.MT$np[which(v.MT$dist <= Dists[d])]*(v.MT$dist[which(v.MT$dist <= Dists[d])]^(-2))
        weights = weights/sum(weights)
        DiffMOMW[d] = sum((v.MT$gamma - v.o$gamma)[v.MT$dist <= Dists[d]]*weights)
        DiffRobustW[d] = sum((rv.MT$gamma - rv.o$gamma)[v.MT$dist <= Dists[d]]*weights)
      }
      
      #Return the sum of the difference between the variograms with and without operator over the spatial area that they operate
      #Using only the full dataset to define the lags computed because when an operator is left out is may result in changing the lag distance slightly.
      #Weighting by the difference per well
      DiffMOM = cumsum((v.MT$gamma - v.o$gamma))[DistLags]/Nop
      DiffRobust = cumsum((rv.MT$gamma - rv.o$gamma))[DistLags]/Nop
      
      #This metric is not recommended because it can make increases in semi-variances seem like decreases.
      #DiffMOM_Jack = sum((v.MT$gamma - v.o$gamma*(nrow(MT) - Nop)/nrow(MT))[v.MT$dist <= MaxLagDist])
      #DiffRobust_Jack = sum((rv.MT$gamma - rv.o$gamma*(nrow(MT) - Nop)/nrow(MT))[rv.MT$dist <= MaxLagDist])
      #This metric is better because it adjusts for number of points for each operator, but it may not be meaningful.
      DiffMOM_Jack = cumsum((v.MT$gamma - v.o$gamma))[DistLags]*(nrow(MT) - Nop)/nrow(MT)
      DiffRobust_Jack = cumsum((rv.MT$gamma - rv.o$gamma))[DistLags]*(nrow(MT) - Nop)/nrow(MT)
      
      #Return the mean of the full dataset and the operator's dataset
      DiffMeans = mean(MT$Qs) - mean(MT$Qs[which(MT$Operator == UniOps[o])])
      t_diff = t.test(x = MT$Qs[-which(MT$Operator == UniOps[o])], y = MT$Qs[which(MT$Operator == UniOps[o])], alternative = 'two.sided', mu = 0)$p.value
      w_diff = wilcox.test(x = MT$Qs[-which(MT$Operator == UniOps[o])], y = MT$Qs[which(MT$Operator == UniOps[o])], alternative = 'two.sided', mu = 0, correct = TRUE)$p.value
      
      #Return the mean and scaled mean of the operator-left-out dataset
      DiffMeanRmOp = mean(MT$Qs) - mean(MT$Qs[-which(MT$Operator == UniOps[o])])
      ScaledDiffMeanRmOp = mean(MT$Qs) - ((nrow(MT) - Nop)*mean(MT$Qs[-which(MT$Operator == UniOps[o])]))/nrow(MT)
      #2-sample t-test for difference in means. These are large samples. Reporting wilcoxon alternative anyway.
      t_RmOp = t.test(x = MT$Qs, y = MT$Qs[-which(MT$Operator == UniOps[o])], alternative = 'two.sided', mu = 0)$p.value
      w_RmOp = wilcox.test(x = MT$Qs, y = MT$Qs[-which(MT$Operator == UniOps[o])], alternative = 'two.sided', mu = 0, correct = TRUE)$p.value
      
      lst = list(Op = UniOps[o], DiffMat = data.frame(DiffMOM = DiffMOM, DiffMOMW = DiffMOMW, DiffMOM_Jack = DiffMOM_Jack, DiffRobust = DiffRobust, DiffRobustW = DiffRobustW, DiffRobust_Jack = DiffRobust_Jack), DiffMeans = DiffMeans, pt_DiffMean = t_diff, pw_DiffMean = w_diff, DiffMeanRmOp = DiffMeanRmOp, pt_DiffRmOp = t_RmOp, pw_DiffRmOp = w_RmOp, ScaledDiffMeanRmOp = ScaledDiffMeanRmOp, NumWells = Nop)
      
    }else{
      #Just report one result supplied by MaxLagDist
      #Return the sum of the difference between the variograms with and without operator over the spatial area that they operate
      #Using only the full dataset to define the lags computed because when an operator is left out is may result in changing the lag distance slightly.
      #Weighting by the difference per well
      DiffMOM = sum((v.MT$gamma - v.o$gamma)[v.MT$dist <= MaxLagDist])/Nop
      DiffRobust = sum((rv.MT$gamma - rv.o$gamma)[rv.MT$dist <= MaxLagDist])/Nop
      
      #Weighting by the separation distance lags
      weights = v.MT$np[which(v.MT$dist <= MaxLagDist)]*(v.MT$dist[which(v.MT$dist <= MaxLagDist)]^(-2))
      weights = weights/sum(weights)
      DiffMOMW = sum((v.MT$gamma - v.o$gamma)[v.MT$dist <= MaxLagDist]*weights)
      DiffRobustW = sum((rv.MT$gamma - rv.o$gamma)[rv.MT$dist <= MaxLagDist]*weights)
      
      #This metric is not recommended because it can make increases in semi-variances seem like decreases.
      #DiffMOM_Jack = sum((v.MT$gamma - v.o$gamma*(nrow(MT) - Nop)/nrow(MT))[v.MT$dist <= MaxLagDist])
      #DiffRobust_Jack = sum((rv.MT$gamma - rv.o$gamma*(nrow(MT) - Nop)/nrow(MT))[rv.MT$dist <= MaxLagDist])
      #This metric is better because it adjusts for number of points for each operator, but it may not be meaningful.
      DiffMOM_Jack = sum((v.MT$gamma - v.o$gamma)[v.MT$dist <= MaxLagDist])*(nrow(MT) - Nop)/nrow(MT)
      DiffRobust_Jack = sum((rv.MT$gamma - rv.o$gamma)[rv.MT$dist <= MaxLagDist])*(nrow(MT) - Nop)/nrow(MT)
      
      #Return the mean of the full dataset and the operator's dataset
      DiffMeans = mean(MT$Qs) - mean(MT$Qs[which(MT$Operator == UniOps[o])])
      t_diff = t.test(x = MT$Qs[-which(MT$Operator == UniOps[o])], y = MT$Qs[which(MT$Operator == UniOps[o])], alternative = 'two.sided', mu = 0)$p.value
      w_diff = wilcox.test(x = MT$Qs[-which(MT$Operator == UniOps[o])], y = MT$Qs[which(MT$Operator == UniOps[o])], alternative = 'two.sided', mu = 0, correct = TRUE)$p.value
      
      #Return the mean and scaled mean of the operator-left-out dataset
      DiffMeanRmOp = mean(MT$Qs) - mean(MT$Qs[-which(MT$Operator == UniOps[o])])
      ScaledDiffMeanRmOp = mean(MT$Qs) - ((nrow(MT) - Nop)*mean(MT$Qs[-which(MT$Operator == UniOps[o])]))/nrow(MT)
      #2-sample t-test for difference in means. These are large samples. Reporting wilcoxon alternative anyway.
      t_RmOp = t.test(x = MT$Qs, y = MT$Qs[-which(MT$Operator == UniOps[o])], alternative = 'two.sided', mu = 0)$p.value
      w_RmOp = wilcox.test(x = MT$Qs, y = MT$Qs[-which(MT$Operator == UniOps[o])], alternative = 'two.sided', mu = 0, correct = TRUE)$p.value
      
      lst = data.frame(Op = UniOps[o], OpNum = o, DiffMOM = DiffMOM, DiffMOMW = DiffMOMW, DiffMOM_Jack = DiffMOM_Jack, DiffRobust = DiffRobust, DiffRobustW = DiffRobustW, DiffRobust_Jack = DiffRobust_Jack, DiffMeans = DiffMeans, pt_DiffMean = t_diff, pw_DiffMean = w_diff, DiffMeanRmOp = DiffMeanRmOp, pt_DiffRmOp = t_RmOp, pw_DiffRmOp = w_RmOp, ScaledDiffMeanRmOp = ScaledDiffMeanRmOp, NumWells = Nop, stringsAsFactors = FALSE)
    }
    return(lst)
  }
}

GetOperatorRanks = function(OpDiag #Output from the bad operator diagnostics function
                            ){
  DiffMOM = DiffMOMW = DiffRobust = DiffRobustW = matrix(0, nrow = nrow(OpDiag), ncol = nrow(OpDiag[1,2][[1]]))
  for (i in 1:nrow(OpDiag)){
    DiffMOM[i,] = OpDiag[i,2][[1]]$DiffMOM
    DiffMOMW[i,] = OpDiag[i,2][[1]]$DiffMOMW
    DiffRobust[i,] = OpDiag[i,2][[1]]$DiffRobust
    DiffRobustW[i,] = OpDiag[i,2][[1]]$DiffRobustW
  }
  
  RankOpDiag_MOMmat = apply(X = DiffMOM, MARGIN = 2, FUN = order, decreasing = TRUE)
  RankOpDiag_MOMWmat = apply(X = DiffMOMW, MARGIN = 2, FUN = order, decreasing = TRUE)
  RankOpDiag_Robustmat = apply(X = DiffRobust, MARGIN = 2, FUN = order, decreasing = TRUE)
  RankOpDiag_RobustWmat = apply(X = DiffRobustW, MARGIN = 2, FUN = order, decreasing = TRUE)
  
  return(list(RankOpDiag_MOMmat = RankOpDiag_MOMmat, RankOpDiag_MOMWmat = RankOpDiag_MOMWmat, RankOpDiag_Robustmat = RankOpDiag_Robustmat, RankOpDiag_RobustWmat = RankOpDiag_RobustWmat, DiffMOM = DiffMOM, DiffMOMW = DiffMOMW, DiffRobust = DiffRobust, DiffRobustW =DiffRobustW))
}

PlotDistanceSensitivity = function(OpDiagRanks, #Output from the GetOperatorRanks function
                                   res = 300, #png resolution
                                   height = 10, #png height in inches
                                   width = 7, #png width in inches
                                   xlim = c(1,60), #x axis limits of plot
                                   NumOps = 100, #y axis number of operators to show
                                   NumOpsStep = 20, #y axis tick mark locations for operators
                                   xStep = 2.5, #x-axis step size
                                   PlotName, #name to be added to the plot file
                                   cols = NA #color palette function for the plots
                                   ){
  if(is.na(cols)){
    cols = colorRampPalette(brewer.pal(7, name = 'PuOr'))
  }
  
  png(paste0(PlotName, '_OpDiagnostics_DistanceSensitivity.png'), res = res, width = width, height = height, units = 'in')
  layout(cbind(c(1,2,3,4)))
  par(mar = c(2,5,3,1))
  for (i in 1:nrow(OpDiagRanks$RankOpDiag_MOMmat)){
    if (i == 1){
      plot(seq(xStep,max(xlim),xStep), which(OpDiagRanks$RankOpDiag_MOMmat == i) - (seq(0,ncol(OpDiagRanks$RankOpDiag_MOMmat)-1,1)*nrow(OpDiagRanks$RankOpDiag_MOMmat)), 
           type = 'o', lty = 1, pch = 16, 
           ylim = rev(c(0,nrow(OpDiagRanks$RankOpDiag_MOMmat))), xlim = xlim, 
           col = cols(nrow(OpDiagRanks$RankOpDiag_MOMmat))[which(OpDiagRanks$RankOpDiag_MOMmat[,1] == i)],
           ylab = 'Operator Rank', xlab = '', main = 'MOM Semi-variogram, Weighted by Number of Wells',
           cex.lab = 1.5, cex.main = 2, axes = FALSE)
      box()
      axis(side = 2, at = seq(0,NumOps,NumOpsStep), labels = TRUE, cex.axis = 1.5)
      axis(side = 4, at = seq(0,NumOps,NumOpsStep), labels = FALSE)
      axis(side = 1, at = seq(xStep,max(xlim),xStep), labels = FALSE)
    }else{
      plot(seq(xStep,max(xlim),xStep), which(OpDiagRanks$RankOpDiag_MOMmat == i) - (seq(0,ncol(OpDiagRanks$RankOpDiag_MOMmat)-1,1)*nrow(OpDiagRanks$RankOpDiag_MOMmat)), 
           type = 'o', lty = 1, pch = 16, 
           ylim = rev(c(0,nrow(OpDiagRanks$RankOpDiag_MOMmat))), xlim = xlim, 
           col = cols(nrow(OpDiagRanks$RankOpDiag_MOMmat))[which(OpDiagRanks$RankOpDiag_MOMmat[,1] == i)],
           axes = FALSE, xlab = '', ylab = '')
    }
    par(new = T)
  }
  par(new = FALSE)
  for (i in 1:nrow(OpDiagRanks$RankOpDiag_MOMmat)){
    if (i == 1){
      plot(seq(xStep,max(xlim),xStep), which(OpDiagRanks$RankOpDiag_Robustmat == i) - (seq(0,ncol(OpDiagRanks$RankOpDiag_Robustmat)-1,1)*nrow(OpDiagRanks$RankOpDiag_Robustmat)), 
           type = 'o', lty = 1, pch = 16, 
           ylim = rev(c(0,nrow(OpDiagRanks$RankOpDiag_Robustmat))), xlim = xlim, 
           col = cols(nrow(OpDiagRanks$RankOpDiag_Robustmat))[which(OpDiagRanks$RankOpDiag_Robustmat[,1] == i)],
           ylab = 'Operator Rank', xlab = '', main = 'Robust Semi-variogram, Weighted by Number of Wells',
           cex.lab = 1.5, cex.main = 2, axes = FALSE)
      box()
      axis(side = 2, at = seq(0,NumOps,NumOpsStep), labels = TRUE, cex.axis = 1.5)
      axis(side = 4, at = seq(0,NumOps,NumOpsStep), labels = FALSE)
      axis(side = 1, at = seq(xStep,max(xlim),xStep), labels = FALSE)
    }else{
      plot(seq(xStep,max(xlim),xStep), which(OpDiagRanks$RankOpDiag_Robustmat == i) - (seq(0,ncol(OpDiagRanks$RankOpDiag_Robustmat)-1,1)*nrow(OpDiagRanks$RankOpDiag_Robustmat)), 
           type = 'o', lty = 1, pch = 16, 
           ylim = rev(c(0,nrow(OpDiagRanks$RankOpDiag_Robustmat))), xlim = xlim, 
           col = cols(nrow(OpDiagRanks$RankOpDiag_Robustmat))[which(OpDiagRanks$RankOpDiag_Robustmat[,1] == i)],
           axes = FALSE, xlab = '', ylab = '')
    }
    par(new = T)
  }
  par(new = FALSE)
  for (i in 1:nrow(OpDiagRanks$RankOpDiag_MOMmat)){
    if (i == 1){
      plot(seq(xStep,max(xlim),xStep), which(OpDiagRanks$RankOpDiag_MOMWmat == i) - (seq(0,ncol(OpDiagRanks$RankOpDiag_MOMWmat)-1,1)*nrow(OpDiagRanks$RankOpDiag_MOMWmat)), 
           type = 'o', lty = 1, pch = 16, 
           ylim = rev(c(0,nrow(OpDiagRanks$RankOpDiag_MOMWmat))), xlim = xlim, 
           col = cols(nrow(OpDiagRanks$RankOpDiag_MOMWmat))[which(OpDiagRanks$RankOpDiag_MOMWmat[,1] == i)],
           ylab = 'Operator Rank', xlab = '', main = 'MOM Semi-variogram, Weighted by Separation Distance',
           cex.lab = 1.5, cex.main = 2, axes = FALSE)
      box()
      axis(side = 2, at = seq(0,NumOps,NumOpsStep), labels = TRUE, cex.axis = 1.5)
      axis(side = 4, at = seq(0,NumOps,NumOpsStep), labels = FALSE)
      axis(side = 1, at = seq(xStep,max(xlim),xStep), labels = FALSE)
    }else{
      plot(seq(xStep,max(xlim),xStep), which(OpDiagRanks$RankOpDiag_MOMWmat == i) - (seq(0,ncol(OpDiagRanks$RankOpDiag_MOMWmat)-1,1)*nrow(OpDiagRanks$RankOpDiag_MOMWmat)), 
           type = 'o', lty = 1, pch = 16, 
           ylim = rev(c(0,nrow(OpDiagRanks$RankOpDiag_MOMWmat))), xlim = xlim, 
           col = cols(nrow(OpDiagRanks$RankOpDiag_MOMWmat))[which(OpDiagRanks$RankOpDiag_MOMWmat[,1] == i)],
           axes = FALSE, xlab = '', ylab = '')
    }
    par(new = T)
  }
  par(new = FALSE)
  for (i in 1:nrow(OpDiagRanks$RankOpDiag_MOMmat)){
    if (i == 1){
      plot(seq(xStep,max(xlim),xStep), which(OpDiagRanks$RankOpDiag_RobustWmat == i) - (seq(0,ncol(OpDiagRanks$RankOpDiag_RobustWmat)-1,1)*nrow(OpDiagRanks$RankOpDiag_RobustWmat)), 
           type = 'o', lty = 1, pch = 16, 
           ylim = rev(c(0,nrow(OpDiagRanks$RankOpDiag_RobustWmat))), xlim = xlim, 
           col = cols(nrow(OpDiagRanks$RankOpDiag_RobustWmat))[which(OpDiagRanks$RankOpDiag_RobustWmat[,1] == i)],
           ylab = 'Operator Rank', xlab = 'Separation Distance (km)', main = 'Robust Semi-variogram, Weighted by Separation Distance',
           cex.lab = 1.5, cex.main = 2, axes = FALSE)
      box()
      axis(side = 2, at = seq(0,NumOps,NumOpsStep), labels = TRUE, cex.axis = 1.5)
      axis(side = 4, at = seq(0,NumOps,NumOpsStep), labels = FALSE)
      axis(side = 1, at = seq(2.5,57.5,5), labels = FALSE)
      axis(side = 1, at = seq(5,60,5), labels = seq(5,60,5), cex.axis = 1.5)
    }else{
      plot(seq(xStep,max(xlim),xStep), which(OpDiagRanks$RankOpDiag_RobustWmat == i) - (seq(0,ncol(OpDiagRanks$RankOpDiag_RobustWmat)-1,1)*nrow(OpDiagRanks$RankOpDiag_RobustWmat)), 
           type = 'o', lty = 1, pch = 16, 
           ylim = rev(c(0,nrow(OpDiagRanks$RankOpDiag_RobustWmat))), xlim = xlim, 
           col = cols(nrow(OpDiagRanks$RankOpDiag_RobustWmat))[which(OpDiagRanks$RankOpDiag_RobustWmat[,1] == i)],
           axes = FALSE, xlab = '', ylab = '')
    }
    par(new = T)
  }
  dev.off()
}

#Plots for the output
PlotOpsDiagnostics = function(OpDiag, #Output from the bad operator diagnostics function
                              OpDiagRanks, #Output from the GetOperatorRanks function
                              col, #column to plot in OpDiagRanks
                              res = 300, #png resolution
                              heightPanel = 5, #png height in inches
                              widthPanel = 10, #png width in inches
                              heightParAx = 5, #png height in inches
                              widthParAx = 10, #png width in inches
                              NumOps = 100, #y axis number of operators to show
                              NumOpsStep = 20, #y axis tick mark locations for operators
                              PlotName, #name to be added to the plot file
                              cols = NA, #color palette function for the parallel axis plots
                              xStep #separation distance step size (km) for variogram. Used to place captions on plots
){
  #Plot for variance represented by robust variogram difference, weighted by number of points the operator has
  # Colored by the bias, represented by the p-value for difference in center for operator data vs. rest of data.
  png(paste0(PlotName, '_OpDiagnostics_PanelPlot.png'), res = res, width = widthPanel, height = heightPanel, units = 'in')
  par(mar = c(5,5,5,3))
  layout(rbind(c(1,2)))
  plot(y = OpDiagRanks$DiffMOM[,col], x = OpDiagRanks$DiffMOMW[,col], 
       col = colFun(as.numeric(OpDiag[,"pw_DiffMean"])), 
       main = 'Difference in Variograms: \n All Data - Operator Removed \n 0 - 30 km Separation', ylab = 'MOM, weight = 1/number operator wells', xlab = 'MOM, weight ~ 1/(distance lag)^2')
  legend('topleft', title = 'p-value', legend = c('< 0.05', '< 0.10', '< 0.15', '< 0.20', '> 0.20'), col = colFun(c(0,0.05,0.10,0.15,0.20,0.21)), pch = 1)
  plot(y = OpDiagRanks$DiffMOM[,col], x = OpDiagRanks$DiffRobust[,col], 
       col = colFun(as.numeric(OpDiag[,"pw_DiffMean"])), 
       main = 'Difference in Variograms: \n All Data - Operator Removed \n 0 - 30 km Separation', ylab = 'MOM, weight = 1/number operator wells', xlab = 'Robust, weight = 1/number operator wells')
  legend('topleft', title = 'p-value', legend = c('< 0.05', '< 0.10', '< 0.15', '< 0.20', '> 0.20'), col = colFun(c(0,0.05,0.10,0.15,0.20,0.21)), pch = 1)
  dev.off()
  
  #plot(OpDiagRanks$DiffRobust, OpDiagRanks$DiffRobustW, col = colFun(OpDiagRanks$pw_DiffMean), main = 'Difference in MOM Semi-variance: 0 to 20 km Separation', xlab = 'Difference in MOM semi-variance: 0 - 20 km Separation', ylab = 'Difference in Robust semi-variance: 0 - 20 km Separation')
  #plot(OpDiagRanks$DiffMOMW, OpDiagRanks$DiffRobustW, col = colFun(OpDiagRanks$pw_DiffMean), main = 'Difference in MOM Semi-variance: 0 to 20 km Separation', xlab = 'Difference in MOM semi-variance: 0 - 20 km Separation', ylab = 'Difference in Robust semi-variance: 0 - 20 km Separation')
  #Univariate plots
  #plot(seq(1,length(OpDiagRanks$DiffMOM),1), abs(sort(as.numeric(OpDiagRanks$DiffMOM))), main = 'Difference in MOM Semi-variance: 0 to 20 km Separation', xlab = 'Sorted Operator ID', ylab = 'Difference')
  #plot(seq(1,length(OpDiagRanks$DiffMOM),1), abs(sort(OpDiagRanks$DiffRobust)), main = 'Difference in Robust Semi-variance: 0 to 20 km Separation', xlab = 'Sorted Operator ID', ylab = 'Difference')
  
  #Make a parallel axis plot of the rank of the operator for each metric
  #Make a matrix for the ranks
  RankMat = cbind(OpDiagRanks$RankOpDiag_MOMmat[,col], 
                  OpDiagRanks$RankOpDiag_Robustmat[,col],
                  OpDiagRanks$RankOpDiag_MOMWmat[,col],
                  OpDiagRanks$RankOpDiag_RobustWmat[,col],
                  order(as.numeric(OpDiag[,"pw_DiffMean"]), decreasing = FALSE))
  
  if(is.na(cols)){
    cols = colorRampPalette(brewer.pal(7, name = 'PuOr'))
  }
  
  png(paste0(PlotName, '_OpDiagnostics_ParallelAxis_Ranks.png'), res = res, width = widthParAx, height = heightParAx, units = 'in')
  par(mar = c(5,5,1,3))
  for (i in 1:nrow(RankMat)){
    if (i == 1){
      plot(seq(1,ncol(RankMat),1), which(RankMat == i) - (seq(0,ncol(RankMat)-1,1)*nrow(RankMat)), 
           type = 'o', lty = 1, pch = 16, 
           ylim = rev(c(0,nrow(RankMat))), xlim = c(1,ncol(RankMat)), 
           col = cols(nrow(RankMat))[which(RankMat[,1] == i)],
           ylab = 'Rank of Operator for Metric', xlab = '', cex.lab = 1.5, axes = FALSE)
      box()
      axis(side = 2, at = seq(0,NumOps,NumOpsStep), labels = TRUE, cex.axis = 1.5)
      axis(side = 1, at = seq(1,ncol(RankMat),1), line = 2, tick = FALSE, padj = -0.3,
           labels = c(paste0('MOM Variogram \n 0 - ', col*xStep, ' km \n w = N wells'),
                      paste0('Robust Variogram \n 0 - ', col*xStep, ' km \n w = N wells'),
                      paste0('MOM Variogram \n 0 - ', col*xStep, ' km \n w = spatial lag'),
                      paste0('Robust Variogram \n 0 - ', col*xStep, ' km \n w = spatial lag'),
                      'p-value \n Wilcoxon Rank Sum \n'))
    }else{
      plot(seq(1,ncol(RankMat),1), which(RankMat == i) - (seq(0,ncol(RankMat)-1,1)*nrow(RankMat)), 
           type = 'o', lty = 1, pch = 16, 
           ylim = rev(c(0,nrow(RankMat))), xlim = c(1,ncol(RankMat)), col = cols(nrow(RankMat))[which(RankMat[,1] == i)],
           axes = FALSE, xlab = '', ylab = '')
    }
    par(new = T)
  }
  dev.off()
  
  #Make a matrix for the ranks
  NormMat = cbind((OpDiagRanks$DiffMOM[,col] - min(OpDiagRanks$DiffMOM[,col])) / (max(OpDiagRanks$DiffMOM[,col]) - min(OpDiagRanks$DiffMOM[,col])),
                  (OpDiagRanks$DiffRobust[,col] - min(OpDiagRanks$DiffRobust[,col])) / (max(OpDiagRanks$DiffRobust[,col]) - min(OpDiagRanks$DiffRobust[,col])),
                  (OpDiagRanks$DiffMOMW[,col] - min(OpDiagRanks$DiffMOMW[,col])) / (max(OpDiagRanks$DiffMOMW[,col]) - min(OpDiagRanks$DiffMOMW[,col])),
                  (OpDiagRanks$DiffRobustW[,col] - min(OpDiagRanks$DiffRobustW[,col])) / (max(OpDiagRanks$DiffRobustW[,col]) - min(OpDiagRanks$DiffRobustW[,col])),
                  as.numeric(OpDiag[,"pw_DiffMean"]))
  
  png(paste0(PlotName, '_OpDiagnostics_ParallelAxis_Normalized.png'), res = res, width = widthParAx, height = heightParAx, units = 'in')
  par(mar = c(5,5,1,3))
  for (i in 1:nrow(NormMat)){
    if (i == 1){
      plot(seq(1,ncol(NormMat),1), NormMat[i,], 
           type = 'o', lty = 1, pch = 16, 
           ylim = c(0,1), xlim = c(1,ncol(NormMat)), 
           col = cols(nrow(RankMat))[which(RankMat[,1] == i)],
           ylab = 'Normalized Value of Metric', xlab = '', cex.lab = 1.5, axes = FALSE)
      box()
      axis(side = 2, at = seq(0,1,.1), labels = TRUE, cex.axis = 1.5)
      axis(side = 1, at = seq(1,ncol(RankMat),1), line = 2, tick = FALSE, padj = -0.3,
           labels = c(paste0('MOM Variogram \n 0 - ', col*xStep, ' km \n w = N wells'),
                      paste0('Robust Variogram \n 0 - ', col*xStep, ' km \n w = N wells'),
                      paste0('MOM Variogram \n 0 - ', col*xStep, ' km \n w = spatial lag'),
                      paste0('Robust Variogram \n 0 - ', col*xStep, ' km \n w = spatial lag'),
                      'p-value \n Wilcoxon Rank Sum \n'))
    }else{
      plot(seq(1,ncol(NormMat),1), NormMat[i,], 
           type = 'o', lty = 1, pch = 16, 
           ylim = c(0,1), xlim = c(1,ncol(NormMat)), col = cols(nrow(RankMat))[which(RankMat[,1] == i)],
           axes = FALSE, xlab = '', ylab = '')
    }
    par(new = T)
  }
  dev.off()
}

#How to retrieve operator number from the sorted results. Useful to check diagnostic plots for operators.
ReturnOperatorPlotNum = function(OpDiag,  #Output from the bad operator diagnostics function
                                 OpRanks, #One metric's dataframe output from the GetOperatorRanks function
                                 col, #column of OpRanks to select 
                                 rank #Rank to select. 1 is highest value in col
){
  RankOp = OpRanks[rank,col]
  as.numeric(strsplit(names(OpDiag[,1][RankOp]), split = '.', fixed = TRUE)[[1]][2])
}

ReturnTopN = function(OpDiag,  #Output from the bad operator diagnostics function
                      OpRanks, #Output from the GetOperatorRanks function
                      col, #column of OpRanks to select 
                      N, #Ranks to select. The top N are selected for each metric
                      ID = TRUE, #Should plot ID number be returned? If FALSE, returns name of operator
                      WellData #Needs a column named Operator to return the name of the operator.
){
  TopN = vector('numeric')
  for (i in 1:N){
    a = ReturnOperatorPlotNum(OpDiag, OpRanks$RankOpDiag_MOMmat, col = 12, rank = i)
    b = ReturnOperatorPlotNum(OpDiag, OpRanks$RankOpDiag_MOMWmat, col = 12, rank = i)
    c = ReturnOperatorPlotNum(OpDiag, OpRanks$RankOpDiag_Robustmat, col = 12, rank = i)
    d = ReturnOperatorPlotNum(OpDiag, OpRanks$RankOpDiag_RobustWmat, col = 12, rank = i)
    TopN = unique(c(a,b,c,d,TopN))
  }
  
  if(ID){
    #PlotID numbers
    return(sort(TopN))
  }else{
    #Operator Names
    return(unique(WellData$Operator)[sort(TopN)])
  }
}
