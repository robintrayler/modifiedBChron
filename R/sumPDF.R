#' Create a single probability distribution function (PDF) from several individual ages.
#'
#' @param ages vector of ages
#' @param ageSds ageSds Vector of 1-sigma values for ages. Must be the same length and given in the same order as \code{ages}
#' @param distTypes Vector of distribution types to model each age as. Choices are 'G' for Gaussian, and 'U' uniform. Must be the same length and given in the same order as \code{ages}
#' @param ids Vector of sample names for each age. All samples with the same ids will be combined into a single age PDF. Must be the same length and given in the same order as \code{ages}
#' @export
#'
sumPDF <- function(ages,
                   ageSds,
                   distTypes = rep('G',length(ages)),
                   ids){

  samples <- unique(ids)
  n_samples <- length(samples)
  ageGrid <- list()
  densities <- list()
  out <- list()

  for(i in 1:n_samples) {
    ageGrid[[i]] <- seq(min(ages[ids == samples[i]] - ageSds[ids == samples[i]] * 10),
                        max(ages[ids == samples[i]] + ageSds[ids == samples[i]] * 10),
                        length = 100000)
    interval <- matrix(0,
                       nrow = 100000,
                       ncol = length(ages[ids == samples[i]]))
    current_n <- length(ages[ids == samples[i]])
    current_ages <- ages[ids == samples[i]]
    current_sds <- ageSds[ids == samples[i]]
    current_dists <- distTypes[ids == samples[i]]

    for(k in 1:current_n) {
      if(current_dists[k] == 'G') {
        interval[, k] <- dnorm(x = ageGrid[[i]],
                               mean = current_ages[k],
                               sd = current_sds[k])

      }
      else if(current_dists[k] == 'U') {
        interval[, k] <- dunif(x = ageGrid[[i]],
                               min = current_ages[k] - current_sds[k],
                               max = current_ages[k] + current_sds[k])
      }
    }

    densities[[i]] <- base::apply(interval, 1, sum) * mean(diff(ageGrid[[i]])) / current_n


    out[[i]] <- data.frame(ageGrid = ageGrid[[i]],
                           densities = densities[[i]])


  }

  names(out) <- samples
  return(out)
}
