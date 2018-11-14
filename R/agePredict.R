#' This function is used to predict the age of new stratigraphic points using the output from an ageModel run.
#'
#' @param model Output of the \code{ageModel} function
#' @param newPositions Vector of new stratigraphic positions to predict the age of. New positions must be within the bounds of the orignal model run
#' @param newPositionThicknesses Vector of stratigraphic uncertanties for each position. Specified as half thicknesses. Must be the same length and given in the same order as \code{newPositions}. Default is no positional uncertanty
#' @return HDI = 95 percent Highest Density Interval for each \code{newPosition}
#' @return raw = Age predictions for each \code{newPosition} for each MCMC run.
#' @export
agePredict <- function(model,
                       newPositions,
                       newPositionThicknesses = rep(0, length(newPositions)),
                       probability = 0.95){
  ## this function is used to predict the age of new points using a modified bchron output model
  ## INPUTS
  ## model = output of the modified bchron age-depth model
  ## newPositions = new stratigraphic points to predict the age of.
  ## The new points must fall within the original stratigraphic range of the model
  ## OUTPUTS
  ## ConfInt = 95% Confidence interval for each point
  ## raw = age predictions for each point for each iteration in the MCMC chain
  pb <- utils::txtProgressBar(min = 1,
                              max = ncol(model$model),
                              style = 3,
                              width=40,
                              char = '<>') # create a progress bar
  x <- model$predictPositions
  n <- length(newPositions)
  predictStore <- matrix(nrow = ncol(model$model),
                         ncol = n)


  for(i in 1:ncol(model$model)){
    currPositions <- runif(n,
                           newPositions - newPositionThicknesses,
                           newPositions + newPositionThicknesses)
    utils::setTxtProgressBar(pb, i) # set the progress bar
    y <- model$model[, i]
    f <- approxfun(x, y)
    predictStore[i, ] <- t(f(currPositions))
  }
  HDI <- t(apply(predictStore,2,quantile,c((1 - probability) / 2, 0.5, (1 + probability) / 2)))
  HDI <- cbind(newPositions,HDI)
  colnames(HDI)[1] <- 'Position'
  predictStore <- data.frame(predictStore)
  colnames(predictStore) <- as.character(newPositions)
  return(list(HDI = HDI,
              raw = predictStore))
}

