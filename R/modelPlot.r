#' This function generates two age model plots given the output of an ageModel run.
#'
#' @param model Output of the \code{ageModel} function
#' @param agePredictOutput Output of the agePredict function. If specified, points and error bars will be added to the age model plots for the predicted points
#' @param scale Scaling factor for age PDFs
#' @param PDF should plots be printed to the current working directory as PDF files? Defaults to FALSE
#' @param ... Optional arguments to be passed to plot legend
#' @export

modelPlot <- function(model,agePredictOutput = NA,scale = 1,PDF = F,...){
  ## This function produces plots from the output object of the modified bchron age-depth model
  ## INPUTS
  ## model = output of the modified bchron age-depth model
  ## PDF = c(T,F). If TRUE plots will be printed to a PDF in the current working directory, else plots are displated in R
  scl <- (diff(range(model$predictPositions))/ncol(model$thetas))/max(model$nAges)*scale # scaling factor for ploting age PDFs
  colsPal <- colorRampPalette(c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6'))
  cols <- colsPal(ncol(model$thetas))

  ##---------------------------------------------------------------------------
  ## if PDF is true print the plots to a PDF file instead
  if(PDF == T){
    pdf(file = 'AgeModelPlots.pdf',
        width = 10,
        height = 10)
  }
  # ##---------------------------------------------------------------------------
  # ## depth on x axis
  # ## time on y axis
  # ## set x and y limits
  # ylims <- c(model$confInt[3,1],
  #            model$confInt[1,length(model$predictPositions)]) # ylimit scaling
  # xlims <- c(min(model$predictPositions),
  #            max(model$predictPositions)+ scl)  # xlimit scaling
  # ##---------------------------------------------------------------------------
  # # open a blank plot
  # plot(NA,
  #      xlim = xlims,
  #      ylim = ylims,
  #      xlab = 'Position',
  #      ylab = 'Age',
  #      tcl = .25,
  #      type = 'n',
  #      main = 'Position - Age Plot')
  # grid() # add gridlines
  #
  # ##---------------------------------------------------------------------------
  # # add a shaded polygon of the 95% HDI
  # polygon(y = c(model$confInt[1,],
  #               rev(model$confInt[3,])),
  #         x = c(model$predictPositions,
  #               rev(model$predictPositions)),
  #         col = rgb(0,0,0,.25),
  #         border = NA)
  # ##---------------------------------------------------------------------------
  # ## add a line for the median model
  # lines(x = model$predictPositions,y = model$confInt[2,],lwd = 2)
  #
  # ##---------------------------------------------------------------------------
  # ## add the likelihood polygons
  # for(i in ncol(model$thetas):1){
  #   polygon(x = model$likelihoods[,i] / max(model$likelihoods[,i])*scl +  model$masterPositions[i],
  #           y = model$ageGrid,
  #           border = NA,
  #           col = cols[i])
  # }
  # ##---------------------------------------------------------------------------
  # ## if there is age predictions to add add them
  # if(all(!is.na(agePredictOutput))){
  #   l <- nrow(agePredictOutput$ConfInt)
  #   for(i in 1:l){
  #     level <- agePredictOutput$ConfInt[i,1]
  #     low <- agePredictOutput$ConfInt[i,2]
  #     mid <- high <- agePredictOutput$ConfInt[i,3]
  #     high <- agePredictOutput$ConfInt[i,4]
  #     arrows(level,low,level,high,length = 0,lwd = 2,col = 'steelblue')
  #     points(level,mid,pch = 19,col = 'steelblue')
  #   }
  # }
  # ##---------------------------------------------------------------------------
  # ## add some legends
  # legend('bottomright',
  #        legend = c('median','95% HDI'),
  #        fill = c('black',rgb(0,0,0,.5)),
  #        bty = 'n')
  # legend('topleft',
  #        legend = rev(model$ids),
  #        fill = rev(cols),
  #        bty = 'n',...)


  #############################################################################
  #############################################################################
  #############################################################################
  ## depth on y axis
  ## time on x axis
  ##---------------------------------------------------------------------------
  ## set x and y limits
  xlims <- c(model$confInt[3,1],
             model$confInt[1,length(model$predictPositions)])
  ylims <- c(min(model$predictPositions),
             max(model$predictPositions)+ scl)
  plot(NA,
       xlim = xlims,
       ylim = ylims,
       xlab = 'Age',
       ylab = 'Position',
       tcl = .25,
       type = 'n',
       main = 'Age - Position Plot')
  grid()
  ##---------------------------------------------------------------------------
  # add a shaded polygon of the 95% HDI
  polygon(x = c(model$confInt[1,],
                rev(model$confInt[3,])),
          y = c(model$predictPositions,
                rev(model$predictPositions)),
          col = rgb(0,0,0,.25),
          border = NA)
  ##---------------------------------------------------------------------------
  ## add a line for the median model
  lines(y = model$predictPositions,
        x = model$confInt[2,],
        lwd = 2)

  ##---------------------------------------------------------------------------
  ## add the likelihood polygons
  for(i in ncol(model$thetas):1){
    polygon(y = model$likelihoods[,i]/max(model$likelihoods[,i])*scl +  model$masterPositions[i],
            x = model$ageGrid,
            border = NA,
            col = cols[i])
  }

  ##---------------------------------------------------------------------------
  ## if there is age predictions to add add them
  if(all(!is.na(agePredictOutput))){
    l <- nrow(agePredictOutput$ConfInt)
    for(i in 1:l){
      level <- agePredictOutput$ConfInt[i, 1]
      low <- agePredictOutput$ConfInt[i, 2]
      mid <- high <- agePredictOutput$ConfInt[i, 3]
      high <- agePredictOutput$ConfInt[i, 4]
      arrows(low,
             level,
             high,
             level,
             length = 0.025,
             angle = 90,
             code = 3,
             lwd = 2,
             col = rgb(0.97, 0.46, 0.43, .5))
      points(mid,
             level,
             pch = 21,
             bg = rgb(0.97, 0.46, 0.43, .5))
    }
  }

  ##---------------------------------------------------------------------------
  ## add some legends
  legend('bottomright',
         lwd = c(2,NA),
         pch = c(NA, 15),
         legend = c('median','95% HDI'),
         col = c('black',rgb(0,0,0,.5)),
         bty = 'n')
  legend('topleft',
         legend = rev(model$ids),
         fill = rev(cols),
         bty = 'n',...)
  if(PDF == T){
    dev.off()
  }
}
