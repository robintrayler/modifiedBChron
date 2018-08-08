
#' This function generates a plot of the age model confidence itervals along with a 2D kernal density estimate.
#' @details This function generates a plot of the age model confidence itervals along with a 2D kernal density estimate of the posterior distribution of each age in both time and position
#' @param model Output of the \code{ageModel} function
#' @param ... Optional arguments to be passed to plot legend
#' @export
#'

ageContourPlot <- function(model,...){
  colsPal <- colorRampPalette(c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6'))
  cols <- colsPal(ncol(model$thetas))
  scl <- (diff(range(model$predictPositions))/ncol(model$thetas))/max(model$nAges)
  xlims <- c(model$confInt[3,1],model$confInt[1,length(model$predictPositions)])
  ylims <- c(min(model$predictPositions),max(model$predictPositions))
  plot(NA,
       xlim = xlims,
       ylim = ylims,
       type='n',
       xlab = 'Age',
       ylab = 'Position',
       tcl = .25)
  grid()
  polygon(x = c(model$confInt[1,],rev(model$confInt[3,])),
          y = c(model$predictPositions,rev(model$predictPositions)),
          col = rgb(0,0,0,.25),
          border = NA)
  lines(y = model$predictPositions,
        x = model$confInt[2,],
        lwd = 2)

  for(n in 1:ncol(model$thetas)){
    x <- MASS::kde2d(model$thetas[model$burn:model$MC,n],jitter(model$positionStore[model$burn:model$MC,n],amount = 0.01),n = 100)
    x$z <- x$z/max(x$z)
    #image(x,add = T,col = colorRampPalette(c(rgb(0,0,0,0),rainbow(ncol(model$thetas),alpha = .5)[n]),alpha = 1)(10))
    contour(x,add = T,nlevels = 5,drawlabels = F,col = cols[n])
  }

  legend('bottomright',
         legend = c('median','95% C.I.'),
         fill = c('black',rgb(0,0,0,.5)),
         bty = 'n')
  legend('topleft',
         legend = rev(model$ids),
         fill = rev(cols),
         bty = 'n',...)

}


