#' This function generates two age model plots given the output of an ageModel run.
#'
#' @param model Output of the \code{ageModel} function
#' @param agePredictOutput Output of the agePredict function. If specified, points and error bars will be added to the age model plots for the predicted points
#' @param scale Scaling factor for age PDFs
#' @param predictLabels c('ages','ids','both', NA) should predicted ages and names for for points be displayed? Defaults to display both. Set to NA for no labels
#' @param type c('PDF', contour') Shound probability be displayed as likelihood input PDFs or as contours of posterior probability. Defaults to PDF
#' @param legend c('color', 'adjacent', NA) type of legend to be drawn. color draws color coded boxes for each sample, adjacent displays sample names next to each PDF. NA omits the legend.
#' @param ... Optional arguments to be passed to plot (xlim, ylim, xlab, ylab, main)
#' @import "RColorBrewer"
#' @export


modelPlot <- function(model,
                      agePredictOutput = NA,
                      predictLabels = 'both',
                      scale = 1,
                      type = 'PDF',
                      legend = 'color',
                      colorBrewerOption = 'RdYlBu',
                      ...) {
  ##-------------------------------------------------------------------------
  ## create a scaling factor and color ramp for likelihood PDFs
  scl <- (diff(range(model$predictPositions)) / ncol(model$thetas)) / max(model$nAges) * scale

  colsPal <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, colorBrewerOption))

  # colsPal <- colorRampPalette(c('#d7191c',
  #                               '#fdae61',
  #                               '#ffffbf',
  #                               '#abd9e9',
  #                               '#2c7bb6'))
  cols <- colsPal(ncol(model$thetas))
  ##-------------------------------------------------------------------------
  ## get a list of optional arguments
  ex <- list(...)
  ##-------------------------------------------------------------------------
  ## assign some default arguments if things arent specified
  ## default x limits
  if(is.null(ex$xlim)){ex$xlim = as.numeric(c(model$HDI[3, 1],
                                              model$HDI[1, length(model$predictPositions)]))}
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

  polygon(x = c(model$HDI[1, ],
                rev(model$HDI[3, ])),
          y = c(model$predictPositions,
                rev(model$predictPositions)),
          col = rgb(0,0,0,.25),
          border = NA)
  ##-------------------------------------------------------------------------
  ## add a line for the median
  lines(y = model$predictPositions,
        x = model$HDI[2, ],
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
                       jitter(model$positionStore[model$burn:model$MC, n], amount = 0.01),
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
    l <- nrow(agePredictOutput$HDI)
    for(i in 1:l){
      arrows(x0 = agePredictOutput$HDI[i, 2],
             y0 = agePredictOutput$HDI[i, 1],
             x1 = agePredictOutput$HDI[i, 4],
             y1 = agePredictOutput$HDI[i, 1],
             length = 0.025,
             angle = 90,
             code = 3,
             lwd = 2,
             col = rgb(0.97, 0.46, 0.43, .5))
      points(agePredictOutput$HDI[i, 3],
             agePredictOutput$HDI[i, 1],
             pch = 21,
             bg = rgb(0.97, 0.46, 0.43, .5))

      if(!is.na(predictLabels)){
        if(predictLabels == 'ages' | predictLabels == 'both'){
          minus <- as.numeric(round(agePredictOutput$HDI[i, 3] - agePredictOutput$HDI[i, 2], 3))
          plus <- as.numeric(round(agePredictOutput$HDI[i, 4] - agePredictOutput$HDI[i, 3], 3))
          median <- as.numeric(round(agePredictOutput$HDI[i, 3],3))
          text(x = agePredictOutput$HDI[i, 4],
               y = agePredictOutput$HDI[i, 1],
               paste(median, '+', plus,'/ -',minus), cex = 0.6,
               pos = 2)
        }
        if(predictLabels == 'ids' | predictLabels == 'both'){
          text(x = agePredictOutput$HDI[i, 2],
               y = agePredictOutput$HDI[i, 1],
               labels = agePredictOutput$HDI$ids[i], cex = 0.6,
               pos = 4)
        }
      }
    }
  }
  ##-------------------------------------------------------------------------
  legend('bottomright',
         lwd = c(2,NA),
         pch = c(NA, 15),
         legend = c('median',paste(paste(model$probability*100,'%',sep = ''),'HDI')),
         col = c('black', rgb(0,0,0,.5)),
         bty = 'n')
  ##-------------------------------------------------------------------------
  if(!is.na(legend)){
    if(legend == 'color'){
      legend('topleft',
             legend = rev(model$ids),
             fill = rev(cols),
             bty = 'n',
             cex = .75,
             ncol = 4)
    }
    if(legend == 'adjacent'){
      for(i in 1:length(model$ids)){
        x <- model$ageGrid[which(cumsum(model$likelihoods[, i]) > .01 & cumsum(model$likelihoods[, i]) < .02)[1]]
        y <- model$masterPositions[i]
        t <- model$ids[i]
        text(x = x, y = y, labels = t, pos = 4, cex = 0.6)
      }
    }
  }
  ##-------------------------------------------------------------------------
}
