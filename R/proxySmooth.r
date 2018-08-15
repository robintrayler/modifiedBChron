#' This function uses a Monte Carlo aproach to repetedly fit a smoothing model to some proxy data incorperating the model uncertanties in age
#' @param proxy Vector proxy compostions. must be numeric
#' @param proxySD optional vector of standard deviations for \code{proxy}. If specified \code{proxy} will be varried from a normal distribution with a mean of \code{proxy} and a standard deviation of \code{proxySD}
#' @param agePredictOutput Output from \{agePredict} for the stratigraphic positions of \code{proxy}
#' @param grid vector of positions to evaluate the smoothing model at.
#' @param method Smoothing method to use. currently supported methods are 'spline', 'movingAverage' and 'polynomial'
#' @param smoothParameter parameter to use for smoothing. If method = 'spline' \code{smoothParameter} should be between 0 and 1. see \code{smooth.spline} for details. If method = 'movingAverage' \code{smoothParameter} is the width of the smoothing window in units of time used by the age-depth model. If method = 'polynomial' \code{smoothParameter} is the degree of polynomial to fit
#' @param HDI desired probability interval to return. Must be between 0 and 1
#' @export
#'
proxySmooth <- function(proxy,
                        proxySD = rep(0, length(proxy)),
                        agePredictOutput,
                        grid,
                        method = c('spline','movingAverage','polynomial'),
                        smoothParameter = 1,
                        HDI = .95){

  # movingAverage Smoothing Function --------------------------------------------
  if(method == 'movingAverage'){
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
  ## Spline Smoothing Function --------------------------------------------
  if(method == 'spline'){
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
  ## Polynomial Smoothing Function ----------------------------------------------
  if(method == 'polynomial'){
    smoothFun <- function(x, y, grid, smoothParameter){
      Df <- data.frame(x = as.vector(x), y = as.vector(y))
      fit <- lm(y ~ poly(x, degree = smoothParameter),
                data = Df)
      ymod <- predict(fit,
                      newdata = data.frame(x = grid))
      return(as.vector(ymod))
    }
  }

  ##---------------------------------------------------------------------------
  ## fit a smoothing model
  pb <- txtProgressBar(min = 0, # start up a text progress bar
                       max = nrow(agePredictOutput$raw),
                       style = 3,
                       width = 40,
                       char = '<>')
  ## preallocate some storage
  modelStore <- matrix(NA,
                       ncol = nrow(agePredictOutput$raw),
                       nrow = length(grid))
  ## apply the smoothing function for each possible set of ages and proxy parameter
  ## values
  for( i in 1:nrow(agePredictOutput$raw)){
    setTxtProgressBar(pb, i) # iterate the progress bar
    ## apply the smoothing function and store the results
    modelStore[, i] <- smoothFun(x = t(as.vector(agePredictOutput$raw[i, ])),
                                 y = rnorm(length(proxy),
                                           mean = proxy,
                                           sd = proxySD),
                                 smoothParameter = smoothParameter,
                                 grid = grid)
  }

  ## calculate the
  confInt <- apply(X = modelStore,
                   MARGIN = 1, FUN = quantile,
                   prob = c((1 - HDI)/2, 0.5, (1 + HDI)/2), na.rm = T)
  return(list(HDI = confInt, raw = modelStore))
}

