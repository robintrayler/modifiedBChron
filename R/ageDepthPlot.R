#' Make plots of the age probability density functions vs position for several samples
#' @param ages Vector of ages
#' @param ageSds Vector of 1-sigma values for ages. Must be the same length and given in the same order as \code{ages}
#' @param positions Vector of stratigraphic positions for ages. Must be the same length and given in the same order as \code{ages}
#' @param ids Vector of sample names for each age. All samples with the same ids will be combined into a single age PDF. Must be the same length and given in the same order as \code{ages}
#' @param positionThicknesses Vector of stratigraphic uncertanties for each age. Specified as half thicknesses. Must be the same length and given in the same order as \code{ages}
#' @param distTypes Vector of distribution types to model each age as. Choices are 'G' for Gaussian, and 'U' uniform. Must be the same length and given in the same order as \code{ages}
#' @param scale scale Scaling factor for age PDFs
#' @param ... Optional arguments to be passed to plot legend
#' @details The \code{ageDepthPlot} function makes a plot of the age probability density functions vs stratigraphic position for several samples. All ages with the same name will be combined into a single PDF.
#' @export

ageDepthPlot <- function(ages,
                         ageSds,
                         positions,
                         ids,
                         positionThicknesses = rep(0,length(ages)),
                         distTypes = rep('G',length(ages)),
                         scale = 1,...){
  ##---------------------------------------------------------------------------
  ## function for creating a summed PDF
  sumPDF <- function(ages,ageSds,distType,x){
    interval <- matrix(0,nrow = length(x), ncol=length(ages))
    for (i in 1:length(ages)) {
      if(distType[i] == 'G') {
        interval[,i] <- dnorm(x,ages[i],ageSds[i])/length(ages)*mean(diff(x)) }
      else if(distType[i] == 'U') {
        interval[,i] <- dunif(x,ages[i]-ageSds[i],ages[i]+ageSds[i])/length(ages)*mean(diff(x))}
    }
    G <- apply(interval,1,sum)
    return(G)
  }

  ##---------------------------------------------------------------------------
  ## order everything
  o = order(positions)
  if(any(positions[o]!=positions)) {
    ages <- ages[o]
    ageSds <- ageSds[o]
    positions <- positions[o]
    distTypes <- distTypes[o]
    positionThicknesses <- positionThicknesses[o]
    ids = ids[o]
  }

  ##---------------------------------------------------------------------------
  ## get some useful counts
  samples <- unique(ids)
  n <- length(samples)
  currPositions <- vector()
  thicknesses <- vector()
  for(i in 1:n){
    currPositions[i] <- positions[ids == samples[i]][1]
    thicknesses[i] <- positionThicknesses[ids == samples[i]][1]
  }

  ##---------------------------------------------------------------------------
  ## create the likelihood PDFs
  prob <-matrix(0,nrow=100000,ncol=n) # create an empty matrix to store the probabilites
  ageGrid <- seq(min(ages-ageSds*10),max(ages+ageSds*10),length.out=100000) # Grid of ages to evaluate over
  for(j in 1:n){
    prob[,j] <- sumPDF(ages = ages[ids == samples[j]],
                       ageSds = ageSds[ids == samples[j]],
                       distType = distTypes[ids == samples[j]],
                       x = ageGrid)
  }

  ##---------------------------------------------------------------------------
  ## set axis limits, determine scaling and colors
  scl <- diff(range(currPositions))/n*scale # scaling factor for ploting age PDFs
  colsPal <- colorRampPalette(c('#d7191c',
                                '#fdae61',
                                '#ffffbf',
                                '#abd9e9',
                                '#2c7bb6'))
  cols <- colsPal(n)
  xlimit <- c(max(ages+ageSds*5),min(ages-ageSds*5))
  ylimit <- c(min(currPositions)-scl,max(currPositions)+scl)

  ##---------------------------------------------------------------------------
  ## open up an empty plot
  plot(NA,type = 'n',
       xlim = xlimit,
       ylim = ylimit,
       tcl = .25,
       xlab = 'Age',
       ylab = 'Position')
  grid()

  ##---------------------------------------------------------------------------
  ## add likelhood PDFs
  for(i in n:1){
    polygon(y = prob[,i]/max(prob[,i])*scl +  currPositions[i],
            x = ageGrid,
            border = NA,
            col = cols[i])
  }

  ##---------------------------------------------------------------------------
  ## add error bars for thickness uncertanties
  if(!all(thicknesses == 0)){
    for(i in n:1){
      if(thicknesses[i] > 0){
        arrows(x0 = mean(ages[ids == samples[i]]),
               x1 =  mean(ages[ids == samples[i]]),
               y0 = currPositions[i] - thicknesses[i],
               y1 = currPositions[i] + thicknesses[i],
               length = .05,,
               code = 3,
               lwd = 2,
               angle = 90)
      }
    }
  }

  ##---------------------------------------------------------------------------
  # ## add a legend
  # legend('bottomright',legend = 'Thicknesses',lwd = 2, bty = 'n')
  legend('topleft',
         legend = rev(samples),
         fill = rev(cols),
         bty = 'n',...)
}
