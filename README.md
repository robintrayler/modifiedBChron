# modifiedBchron

This package provides a version of the [Bchron](https://github.com/andrewcparnell/Bchron/tree/master/R) R package, with features added to improve usability with deep time geochronology data (e.g., Argon-Argon, Uranium-Lead). If you use this package please cite Trayler et al.[^1] the deep time additions *and* Haslett & Parnell[^2] for the underlying `Bchron` model framework.

## Introduction

Bchron is a bayesian age-depth model implemented in R. Since it was originally designed for radiocarbon analyses It lacks some features that are desirable to those working with U-Pb or <sup>40</sup>Ar/<sup>39</sup>Ar geochronology data. In this package we have made several modifications including 

* We allow individual dates to be grouped to form complex probability density functions, reproducing the practice in both U-Pb and <sup>40</sup>Ar/<sup>39</sup>Ar geochronology of convolving many single crystal or spot analyses into a single age distribution.
* We have added an adaptive Markov Chain Monte Carlo Algorithm[^3] to ensure efficient exploration off parameter space, independent of the scale of the data (e.g., 1 Ma vs 1,000 ka, vs 1,000,000 a)
* We have removed automated outlier rejection. We instead recommend pre-screening data using established criteria based on the physical mechanisms of crystal growth and open system behavior in volcanic rocks. 

## Installation

modifiedBchron can be installed using the [`devtools`](https://github.com/r-lib/devtools) R package.  You may also need to install [Rtools](https://cran.rstudio.com/bin/windows/Rtools/) on Windows or on macOS you may need to install command line tools. 

```
# copy this code into macOS Terminal (not R!)
xcode-select --install
```

```r
# copy this code into R
# install.packages('devtools')
devtools::install_github('robintrayler/modifiedBchron')
```

Once `modifiedBchron` is installed it can be loaded as an R package by adding `library(modifiedBchron)` to the beginning of an R script. 

## Usage
### Core Functions
* `ageModel()`
* `modelPlot()`
* `agePredict()`

[^1]: Trayler, R.B., Schmitz, M.D., Cuitiño, J.I., Kohn, M.J., Bargo, M.S., Kay, R.F., Strömberg, C.A.E., and Vizcaíno, S.F., 2020, An Improved Approach To Age-Depth Modeling In Deep Time: Implications For The Santa Cruz Formation, Argentina: Geological Society of America Bulletin, v. 132, p. 233–244.

[^2]: Haslett, J., and Parnell, A.C., 2008, A Simple Monotone Process With Application To Radiocarbon-Dated Depth Chronologies: Applied Statistics, v. 57, p. 399–418, doi:doi: 10.1111/j.1467-9876.2008.00623.x.

[^3]: Haario, H., Saksman, E., and Tamminen, J., 1999, Adaptive proposal distribution for random walk Metropolis algorithm: Computational Statistics, v. 14, p. 375–396.

