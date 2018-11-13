#' This function generates two age model plots given the output of an ageModel run.
#'
#' @param model Output of the \code{ageModel} function
#' @param agePredictOutput Output of the agePredict function. If specified, points and error bars will be added to the age model plots for the predicted points
#' @param scale Scaling factor for age PDFs
#' @param predictLabels should predicted ages for points be displayed?
#' @param type c('PDF', contour') Shound probability be displayed as likelihood input PDFs or as contours of posterior probability. Defaults to PDF
#' @param ... Optional arguments to be passed to plot (xlim, ylim, xlab, ylab, main)
#' @export
#

modelPlot <- function(model,
                      agePredictOutput = NA,
                      predictLabels = T,
                      scale = 1,
                      type = 'PDF',...){
  ##-------------------------------------------------------------------------
  ## create a scaling factor and color ramp for likelihood PDFs
  scl <- (diff(range(model$predictPositions)) / ncol(model$thetas)) / max(model$nAges) * scale
  colsPal <- colorRampPalette(c('#d7191c',
                                '#fdae61',
                                '#ffffbf',
                                '#abd9e9',
                                '#2c7bb6'))
  cols <- colsPal(ncol(model$thetas))
  ##-------------------------------------------------------------------------
  ## get a list of optional arguments
  ex <- list(...)
  ##-------------------------------------------------------------------------
  ## assign some default arguments if things arent specified
  ## default x limits
  if(is.null(ex$xlim)){ex$xlim = as.numeric(c(model$confInt[3, 1],
                                              model$confInt[1, length(model$predictPositions)]))}
  ## default y limits
  if(is.null(ex$ylim)){ex$ylim = c(min(model$predictPositions),
                                   max(model$predictPositions) + scl)}
  ## default x label
  if(is.null(ex$xlab)){ex$xlab = 'Age'}
  ## default y label
  if(is.null(ex$ylab)){ex$ylab = 'Position'}
  ## required defaults
  ex$tcl = 0.25
  ex$x = 1
  ex$y = 1
  ex$type = 'n'
  args <- ex
  ##-------------------------------------------------------------------------
  ## open up a blank plot
  do.call('plot', args)
  graphics::grid()
  ##-------------------------------------------------------------------------
  ## draw a polygon for the confidence interval
  polygon(x = c(model$confInt[1,],
                rev(model$confInt[3,])),
          y = c(model$predictPositions,
                rev(model$predictPositions)),
          col = rgb(0,0,0,.25),
          border = NA)
  ##-------------------------------------------------------------------------
  ## add a line for the median
  lines(y = model$predictPositions,
        x = model$confInt[2,],
        lwd = 2)
  ##-------------------------------------------------------------------------
  ## add polygons for each likelihood PDF
  if(type == 'PDF'){
    for(i in ncol(model$thetas):1){
      polygon(y = model$likelihoods[, i] / max(model$likelihoods[, i]) * scl +  model$masterPositions[i],
              x = model$ageGrid,
              border = NA,
              col = cols[i])
    }
  }
  ##-------------------------------------------------------------------------
  if(type == 'contour'){
    for(n in 1:ncol(model$thetas)){
      x <- MASS::kde2d(model$thetas[model$burn:model$MC, n],
                       jitter(model$positionStore[model$burn:model$MC, n],amount = 0.01),
                       n = 100)
      x$z <- x$z/max(x$z)
      #image(x,add = T,col = colorRampPalette(c(rgb(0,0,0,0),rainbow(ncol(model$thetas),alpha = .5)[n]),alpha = 1)(10))
      contour(x,
              add = T,
              nlevels = 5,
              drawlabels = F,
              col = cols[n])
    }
  }
  ##-------------------------------------------------------------------------
  if(all(!is.na(agePredictOutput))){
    l <- nrow(agePredictOutput$ConfInt)
    for(i in 1:l){
      arrows(x0 = agePredictOutput$ConfInt[i, 2],
             y0 = agePredictOutput$ConfInt[i, 1],
             x1 = agePredictOutput$ConfInt[i, 4],
             y1 = agePredictOutput$ConfInt[i, 1],
             length = 0.025,
             angle = 90,
             code = 3,
             lwd = 2,
             col = rgb(0.97, 0.46, 0.43, .5))
      points(agePredictOutput$ConfInt[i, 3],
             agePredictOutput$ConfInt[i, 1],
             pch = 21,
             bg = rgb(0.97, 0.46, 0.43, .5))
      if(predictLabels == T){
        minus <- as.numeric(round(agePredictOutput$ConfInt[i, 3] - agePredictOutput$ConfInt[i, 2], 3))
        plus <- as.numeric(round(agePredictOutput$ConfInt[i, 4] - agePredictOutput$ConfInt[i, 3], 3))
        median <- as.numeric(round(agePredictOutput$ConfInt[i, 3],3))
        text(x = agePredictOutput$ConfInt[i, 4],
             y = agePredictOutput$ConfInt[i, 1],
             paste(median, '+', plus,'/ -',minus), cex = 0.6,
             pos = 2)
      }
    }
  }


  legend('bottomright',
         lwd = c(2,NA),
         pch = c(NA, 15),
         legend = c('median','95% HDI'),
         col = c('black', rgb(0,0,0,.5)),
         bty = 'n')
  legend('topleft',
         legend = rev(model$ids),
         fill = rev(cols),
         bty = 'n',
         cex = .75,
         ncol = 4)
}
