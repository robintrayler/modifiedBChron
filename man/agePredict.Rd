% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/agePredict.r
\name{agePredict}
\alias{agePredict}
\title{This function is used to predict the age of new stratigraphic points using the output from an ageModel run.}
\usage{
agePredict(
  model,
  newPositions,
  ids = 1:length(newPositions),
  newPositionThicknesses = rep(0, length(newPositions)),
  probability = 0.95
)
}
\arguments{
\item{model}{Output of the \code{ageModel} function}

\item{newPositions}{Vector of new stratigraphic positions to predict the age of. New positions must be within the bounds of the orignal model run}

\item{ids}{optional vector of names for new positions. Defaults to 1, 2, 3, ...}

\item{newPositionThicknesses}{Vector of stratigraphic uncertanties for each position. Specified as half thicknesses. Must be the same length and given in the same order as \code{newPositions}. Default is no positional uncertanty}
}
\value{
HDI = 95 percent Highest Density Interval for each \code{newPosition}

raw = Age predictions for each \code{newPosition} for each MCMC run.
}
\description{
This function is used to predict the age of new stratigraphic points using the output from an ageModel run.
}
