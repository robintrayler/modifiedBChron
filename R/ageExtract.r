#' Extract the highest density intervals for the dated horizons using an ageModel output
#' @param model model Output of the \code{ageModel} function
#' @return HDI = specified probability highest density interval of the posterior distribution for each input horizon
#'
#' @export
ageExtract <- function(model, probability = 0.95){
  HDI <- t(apply(X = model$thetas,
                 MARGIN = 2,
                 FUN = quantile,
                 prob = c((1 - probability)/2, 0.5, (1 + probability)/2)))
  extract <- data.frame(ids = model$ids, HDI = as.matrix(HDI))
  colnames(extract) <- c('ids',
                         paste((1 - probability)/2*100, '%', sep = ''),
                         paste(50, '%', sep = ''),
                         paste((1 + probability)/2*100, '%', sep = ''))
  return(HDI = extract)
}
