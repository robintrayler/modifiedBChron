#' This function plots the result of a time-integrated smoothed proxy record
#' @param proxySmoothOutput the output of the proxySmooth function
#' @param model output of the age model used to predict ages and smooth the data
#' @param stratSmoothParameter optional parameter for smoothing the proxy data vs stratigraphic position. must be meet the same criteria as \code{smoothParameter} in \code{proxySmooth}.
#' Default method for spline fitting uses the same \code{smoothParameter} as proxySmooth. Default method for moving average fitting will estimate a smoothing parameter based on the median accumulation
#' rate and \code{smoothParameter} from \code{proxySmooth}.
#'
#' @export
#'
proxyPlot <- function(proxySmoothOutput,
                      model,
                      stratSmoothParameter = 'default'){
  # movingAverage Smoothing Function --------------------------------------------
  if(proxySmoothOutput$method == 'movingAverage'){
    smoothFun <- function(x, y, grid, smoothParameter){
      ymod <- vector(length = length(grid)) # preallocate
      for (i in 1:length(grid)) { # for each value in grid
        w <- dnorm(x, grid[i], smoothParameter) # weights
        ymod[i] <- sum(w * y) / sum(w) # calculate the moving weighted mean
      }
      return(ymod)
      return(as.vector(ymod))
    }
  }

  if(proxySmoothOutput$method == 'movingAverage' & stratSmoothParameter == 'default'){
    stratSmoothParameter  <- abs(mean(diff(model$predictPositions))/mean(diff(model$confInt[2, ]))*proxySmoothOutput$smoothParameter)
  } else {
    stratSmoothParameter <- proxySmoothOutput$smoothParameter
  }

  ## Spline Smoothing Function --------------------------------------------
  if(proxySmoothOutput$method == 'spline'){
    smoothFun <- function(x, y, grid, smoothParameter){
      ## INPUTS
      ## x = x values
      ## y = yvalues
      ## grid = positions to calculate statistics at. default is data positions
      ## smoothParameter = standard deviation for gaussian kernel
      fit <- smooth.spline(x,
                           y,
                           spar = smoothParameter)
      fun <- approxfun(x = fit$x,
                       y = fit$y)
      ymod <- fun(grid)
      return(as.vector(ymod))
    }
  }

  ##---------------------------------------------------------------------------
  layout(rbind(c(1,1,2,2),c(1,1,3,3)))
  modelPlot(model = model, agePredictOutput = proxySmoothOutput$agePredictOutput)
  ##---------------------------------------------------------------------------
  xlims <- range(proxySmoothOutput$agePredictOutput$ConfInt[,1])
  ylims <- range(proxySmoothOutput$proxy)

  plot(NA,
       xlim = xlims,
       ylim = ylims,
       tcl = 0.25,
       xlab = 'Position',
       ylab = 'Proxy')
  grid()
  ##---------------------------------------------------------------------------
  points(x = proxySmoothOutput$agePredictOutput$ConfInt[,1],
         y = proxySmoothOutput$proxy,
         bg = rgb(0.97, 0.46, 0.43, 1),
         pch = 21,
         cex = 1.5)
  depthGrid = seq(min(xlims),max(xlims),length = 1000)
  smooth <- smoothFun(x =  proxySmoothOutput$agePredictOutput$ConfInt[, 1],
                      y = proxySmoothOutput$proxy,
                      grid = depthGrid,
                      smoothParameter = stratSmoothParameter)
  lines(depthGrid, smooth, lwd = 2)

  ##---------------------------------------------------------------------------
  xlims <- rev(range(model$confInt[2, ]))
  ylims <- range(proxySmoothOutput$proxy)
  plot(NA,
       xlim = xlims,
       ylim = ylims,
       tcl = 0.25,
       xlab = 'Age',
       ylab = 'Proxy')
  grid()
  arrows(x0 = proxySmoothOutput$agePredictOutput$ConfInt[ ,2],
         x1 = proxySmoothOutput$agePredictOutput$ConfInt[ ,4],
         y0 = proxySmoothOutput$proxy,
         length = 0.05,
         angle = 90,
         col = rgb(0.9725490, 0.4627451, 0.4274510, .15),
         lwd = 4,
         code = 3)


  points(x = (proxySmoothOutput$agePredictOutput$ConfInt[,3]),
         y = proxySmoothOutput$proxy,
         bg = rgb(0.97, 0.46, 0.43, 1),
         pch = 21,
         cex = 1.5)



  polygon(x = c(rev(proxySmoothOutput$grid),proxySmoothOutput$grid),
          y = c(rev(proxySmoothOutput$HDI[1, ]),proxySmoothOutput$HDI[3, ]),
          col = rgb(0.97, 0.46, 0.43, .5),
          border = NA)

  lines(proxySmoothOutput$grid,
        proxySmoothOutput$HDI[2, ],
        lwd = 2)

}
