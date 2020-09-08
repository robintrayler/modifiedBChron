#' Create a single probability distribution function (PDF) from several individual ages.
#'
#' @param ages vector of ages
#' @param ageSds ageSds Vector of 1-sigma values for ages. Must be the same length and given in the same order as \code{ages}
#' @param distTypes Vector of distribution types to model each age as. Choices are 'G' for Gaussian, and 'U' uniform. Must be the same length and given in the same order as \code{ages}
#' @param ids Vector of sample names for each age. All samples with the same ids will be combined into a single age PDF. Must be the same length and given in the same order as \code{ages}
#' @export

sumPDF <- function(ages,
                   ageSds,
                   distTypes = rep('G',length(ages)),
                   ids){
  samples <- unique(ids)
  n <- length(samples)
  ageGrid <- list()
  densities <- list()
  out <- list()
  for(i in 1:n){
    ageGrid[[i]] <- seq(min(ages[ids == samples[i]] - ageSds[ids == samples[i]]*6),
                       max(ages[ids == samples[i]] + ageSds[ids == samples[i]]*6), length = 100000)
    interval <- matrix(0,nrow = 100000, ncol = length(ages[ids == samples[i]]))

    for(j in 1:length(ages[ids == samples[i]])){
      if(distTypes[j] == 'G'){
        interval[, j] <- dnorm(ageGrid[[i]],
                               mean = ages[ids == samples[i]][j],
                              sd = ageSds[ids == samples[i]][j])/
          length(ages[ids == samples[i]][j])*mean(diff(ageGrid[[i]]))
      }
      else if(distTypes[j] == 'U'){
        interval[, j] <- dunif(ageGrid[[i]],
                              min = ages[ids == samples[i]][j] - ageSds[ids == samples[i]][j],
                              max = ages[ids == samples[i]][j] + ageSds[ids == samples[i]][j])/
          length(ages[ids == samples[i]])*mean(diff(ageGrid[[i]]))
      }

    }
    densities[[i]] <- base::apply(interval, 1, sum)
    out[[i]] <- data.frame(ageGrid = ageGrid[[i]],
                     densities = densities[[i]])

  }
  names(out) <- samples
  return(out)
}

