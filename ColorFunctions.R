#Function for custom color plotting based on the input data values.
#ColFun takes as input a vector of data that needs colors assigned based on the scheme defined by PalFun

#Variables needed to run this function include:
#Pal - Color pallete from which to get colors. Must be of length (scaleRange[2] - scaleRange[1])/scaleBy + 1
#scaleRange - The c(lower, upper) values of the data to be colored. 
             #Values greater than the upper value will be assigned the last color in Pal. 
             #Values less than the smaller value will be assigned the first color in Pal.
#scaleBy - The increment of Data for assigning colors

#Function assumes that values must be greater than the thresholds to be assigned that color.

colFun = function(Data){
  #Determine the colors for all points
  cols = Pal[as.numeric(lapply((floor((Data - scaleRange[1])/scaleBy) + 1), FUN = PalFun))]
  return(cols)
}

#Function called using lapply in the colFun. Used to select the index of colors in the pallete
PalFun = function(x){
  ifelse(x <= 0,
         #Data value smaller than the lower bound of the scaleRange. Assign the first color in Pal
         1, 
         ifelse(x > length(Pal),
                #Data value greater than the upper bound of the scaleRange. Assign the last color in Pal
                length(Pal),
                #Index will work as is.
                x
         )
  )
}

# Example Use:
#colPal = colorRampPalette(colors = rev(c('red', 'orange', 'yellow', 'green', 'blue')))
#scaleRange = c(30,80)
#scaleBy = 10
#Pal = colPal((scaleRange[2] - scaleRange[1])/scaleBy + 1)