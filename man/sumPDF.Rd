% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sumPDF.r
\name{sumPDF}
\alias{sumPDF}
\title{Create a single probability distribution function (PDF) from several individual ages.}
\usage{
sumPDF(ages, ageSds, distTypes = rep("G", length(ages)), ids)
}
\arguments{
\item{ages}{vector of ages}

\item{ageSds}{ageSds Vector of 1-sigma values for ages. Must be the same length and given in the same order as \code{ages}}

\item{distTypes}{Vector of distribution types to model each age as. Choices are 'G' for Gaussian, and 'U' uniform. Must be the same length and given in the same order as \code{ages}}

\item{ids}{Vector of sample names for each age. All samples with the same ids will be combined into a single age PDF. Must be the same length and given in the same order as \code{ages}}
}
\description{
Create a single probability distribution function (PDF) from several individual ages.
}
