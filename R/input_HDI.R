input_HDI <- function(ageModel,
                      prob = 0.95) {
  # number of likelihoods
  n = ncol(model$likelihoods)

  HDI_table <- data.frame(id = model$ids,
                          HDI_low  = rep(NA, n),
                          HDI_median   =  rep(NA, n),
                          HDI_high =  rep(NA, n))

  # calculate HDI for each sample
  for(i in 1:n) {
    # calculate the CDF
    cdf <- cumsum(model$likelihoods[,i])

    # low HDI
    HDI_table$HDI_low[i]  <- model$ageGrid[which.min(abs(cdf - (1 - prob)/2))] # low
    HDI_table$HDI_median[i]   <- model$ageGrid[which.min(abs(cdf - 0.5))]   # median
    HDI_table$HDI_high[i] <- model$ageGrid[which.min(abs(cdf - (1 + prob)/2))] # high
  }

  return(HDI_table)
}
