#' Run a modified version of the Bchronology algorithm from the R package "Bchron" (Haslett and Parnell, 2008)
#'
#' @param ages Vector of ages
#' @param ageSds Vector of 1-sigma values for ages. Must be the same length and given in the same order as \code{ages}
#' @param positions Vector of stratigraphic positions for ages. Must be the same length and given in the same order as \code{ages}
#' @param positionThicknesses Vector of stratigraphic uncertanties for each age. Specified as half thicknesses. Must be the same length and given in the same order as \code{ages}
#' @param ids Vector of sample names for each age. All samples with the same ids will be combined into a single age PDF. Must be the same length and given in the same order as \code{ages}
#' @param distTypes c('G','U') Vector of distribution types to model each age as. Choices are 'G' for Gaussian, and 'U' uniform. Must be the same length and given in the same order as \code{ages}
#' @param MC Number of MCMC interation to run for. Defaults to 10000
#' @param burn Number of initial iterations to discard. Defaults to 2000
#' @param adapt Should the proposal standard deviation for model parameters be determined using the adaptive proposal algorithm of Haario et al., 1998. Defaults to TRUE
#' @param mhSD Metropolis-Hastings proposal standard deviation for the age parameters. If adapt = TRUE it will be determined by the model run. Otherwise a vector of standard deviations must be specified with length = \code{unique(ids)}
#' @param psiSD Metropolis-Hastings proposal standard deviation for the Compound Poisson-Gamma scale parameter. If adapt = TRUE it will be determined by the model run.
#' @param muSD Metropolis-Hastings proposal standard deviation for the Compound Poisson-Gamma mean parameter. If adapt = TRUE it will be determined by the model run.
#' @param predictPositions Vector of stratigraphic positions to evaluate the model at. Defaults to 100 evenly spaced points from the bottom to top of the section
#' @param truncateUp Truncation age for extrapolating above the top of the section. Defaults to 0
#' @param truncateDown Truncation age for extrapolation below the bottom of the section. Defaults to 1e10 (a really big number)
#' @param extrapUp Number of extrapolations to perform above the top of the section. defaults to 100
#' @param extrapDown Number of extrapolations to perform below the bottom of the section. defaults to 100
#' @return ConfInt = 95 percent confidence intervals for the model run
#' @return model = Raw model predictions for each MCMC iteration
#' @return thetas = The posterior values for each dated horizon from each MCMC run. burn is not removed.
#' @return predictPositions = Stratigraphic positions where the model was evaluated. Same as \code{predictPositions}
#' @return mu, psi = The posterior values for the Compound Poisson-Gamma mean and scale parameters. burn is not removed.
#' @return ageGrid = Grid that age PDFs were evaluated over. Useful for plotting
#' @return likelihoods = Age PDFs for each set of dated samples. Useful for plotting
#' @return nAges = Number of ages that were combined into each sample. Useful for plotting
#' @return masterPositions = Stratigraphic positions for each set of combined ages. Useful for plotting
#' @return ids = Vector of unique sample names. Useful for plotting
#' @return psiSDStore, muSDStore, mhSDStore = proposal standard deviations for each model parameter throughout the MCMC run. burn is not removed.
#' @useDynLib modifiedBChron
#' @export

ageModel <- function(ages,
                     ageSds,
                     positions,
                     positionThicknesses = rep(0, length(positions)),
                     ids,
                     distTypes = rep('G',length(ages)),
                     MC = 10000,
                     burn = 2000,
                     adapt = TRUE,
                     mhSD = rep(1,length(unique(ids))),
                     psiSD = 1,
                     muSD = 1,
                     predictPositions = seq(min(positions),max(positions),length = 100),
                     truncateUp = 0,
                     extrapUp = 100,
                     truncateDown = 1e10,
                     extrapDown = 100){

  ##############################################################################################################################
  # C function for truncated random walk
  truncatedWalk = function(old,sd,low,high) {
    if(isTRUE(all.equal(low, high, tolerance = 1e-10))) return(list(new=low,rat=1))
    new = .C('truncatedWalk',
             as.double(old),
             as.double(sd),
             as.double(low),
             as.double(high),
             as.double(0),
             PACKAGE = 'modifiedBChron')[5][[1]]
    rat =.C('truncatedRat',
            as.double(old),
            as.double(sd),
            as.double(low),
            as.double(high),
            as.double(new),
            as.double(0),
            PACKAGE = 'modifiedBChron')[6][[1]]
    #if(is.nan(rat)) rat=1 # Useful for when proposed value and new value are identical
    return(list(new=new,rat=rat))
  }

  # C functions for tweedie
  dtweediep1 = Vectorize(function(y,p,mu,phi) {
    return(.C('dtweediep1',
              as.double(y),
              as.double(p),
              as.double(mu),
              as.double(phi),
              as.double(0),
              PACKAGE = 'modifiedBChron')[5][[1]])
  })

  # Some C functions to do prediction and interpolation
  predictInterp = function(alpha,
                           lambda,
                           beta,
                           predictPositions,
                           diffPositionj,
                           currPositionsj,
                           currPositionsjp1,
                           thetaj,
                           thetajp1) {
    return(.C('predictInterp',
              as.double(alpha),
              as.double(lambda),
              as.double(beta),
              as.double(predictPositions),
              as.integer(length(predictPositions)),
              as.double(diffPositionj), as.double(currPositionsj),
              as.double(currPositionsjp1), as.double(thetaj),
              as.double(thetajp1),
              as.double(rep(0,length(predictPositions))),
              PACKAGE = 'modifiedBChron')[11][[1]])
  }
  predictExtrapUp = function(alpha,
                             lambda,
                             beta,
                             predictPositions,
                             currPositions1,
                             theta1,
                             maxExtrap,
                             extractDate) {
    return(.C('predictExtrapUp',
              as.double(alpha),
              as.double(lambda),
              as.double(beta),
              as.double(predictPositions),
              as.integer(length(predictPositions)),
              as.double(currPositions1),
              as.double(theta1),
              as.integer(maxExtrap),
              as.double(extractDate),
              as.double(rep(0,length(predictPositions))),
              PACKAGE = 'modifiedBChron')[10][[1]])
  }

  predictExtrapDown = function(alpha,
                               lambda,
                               beta,
                               predictPositions,
                               currPositions1,
                               theta1,
                               maxExtrap,
                               extractDate) {
    return(.C('predictExtrapDown',
              as.double(alpha),
              as.double(lambda),
              as.double(beta),
              as.double(predictPositions),
              as.integer(length(predictPositions)),
              as.double(currPositions1),
              as.double(theta1),
              as.integer(maxExtrap),
              as.double(extractDate),
              as.double(rep(0,length(predictPositions))),
              PACKAGE = 'modifiedBChron')[10][[1]])
  }

  ##############################################################################################################################
  ## function for creating a summed PDF of ages
  compoundProb <- function(ages,sigs,distType,x){
    interval <- matrix(0,nrow = length(x), ncol=length(ages))
    for (i in 1:length(ages)) {
      if(distType[i] == 'G') {
        interval[,i] <- dnorm(x,ages[i],sigs[i])/length(ages)*mean(diff(x)) }
      else if(distType[i] == 'U') {
        interval[,i] <- dunif(x,ages[i]-sigs[i],ages[i]+sigs[i])/length(ages)*mean(diff(x)) }
    }
    G <- apply(interval,1,sum)
    return(G)
  }



  o = order(positions)
  if(any(positions[o]!=positions)) {
    warning("positions not given in order - re-ordering")
    ages <- ages[o]
    ageSds <- ageSds[o]
    positions <- positions[o]
    distTypes <- distTypes[o]
    positionThicknesses <- positionThicknesses[o]
    ids = ids[o]
  }


  ##-------------------------------------------------------------------
  ## sort the data
  nSamples <- length(unique(ids)) # get the number of unique date layers
  masterPositions <- vector()
  nNames <- unique(ids) # get a list of each unique sample name
  for(i in 1:nSamples){
    masterPositions[i] <- positions[ids == nNames[i]][1]
  }


  nThicknesses <- vector(length = nSamples)
  for(i in 1:nSamples){
    nThicknesses[i] <- unique(positionThicknesses[positions == masterPositions[i]])
  }

  ##-------------------------------------------------------------------
  # count the number of ages at each unique stratigraphic position
  # this is used for scaling the plots later
  nAges <- rep(NA,length(nSamples))
  for(i in 1:nSamples){
    nAges[i] <- length(ids[ids == nNames[i]])
  }


  ##-------------------------------------------------------------------
  ## generate the summed probability distirbutions at each unique depth position
  prob <-matrix(0,nrow=100000,ncol=nSamples) # create an empty matrix to store the probabilites
  ageGrid <- seq(min(ages-ageSds*10),max(ages+ageSds*10),length.out=100000) # Grid of ages to evaluate over
  for(j in 1:nSamples){
    prob[,j] <- compoundProb(ages[ids == nNames[j]],
                             ageSds[ids == nNames[j]],
                             distType = distTypes[ids == nNames[j]],
                             x = ageGrid)
  }
  rm(j)


  ##-------------------------------------------------------------------
  ## prealocate some storage
  thetaStore <- matrix(NA,nrow=MC,ncol=nSamples)
  muStore <- vector(length=MC)
  psiStore <- muStore
  modelStore <- matrix(NA,ncol=MC,nrow=length(predictPositions))
  positionStore <- matrix(NA,nrow=MC,ncol=nSamples)
  predictStore <-  matrix(NA,ncol=MC,nrow=length(predictPositions))
  psiSDStore <- psiStore
  muSDStore <- psiStore
  mhSDStore <- thetaStore
  ##-------------------------------------------------------------------
  ## get some starting guesses for MH sampling and make sure they conform to superposition
  thetas <- vector(length = nSamples) # create a vector to store current thetas
  bad = TRUE
  count <- 0
  print('Finding starting ages')
  while(bad == TRUE){
    for (j in 1:nSamples){ # for each age
      thetas[j] <- sample(ageGrid, 1, prob = prob[, j],replace = TRUE) # sample the ageGrid for probabilistically for each sample
    }
    if(all(sign(diff(thetas)) == -1)) {bad = FALSE}

    if(count%%1000 == 0) {print(paste('starting ages not found after, ',count,' attempts'))}
    count <- count + 1
    if(count == 10000){stop('Unable to find acceptable starting ages. Consider checking for outliers')}
  }
  rm(j)
  print(paste('starting ages found after', count, 'attempts'))
  ##-------------------------------------------------------------------
  ## calculate some initial parameters
  ## based on Haslett and Parnell (2008)
  currPositions <- masterPositions
  mu <- abs(rnorm(1,mean=mean(diff(thetas))/mean(diff(currPositions)),sd=muSD))
  psi <- abs(rnorm(1,1,sd=psiSD))
  p = 1.2
  alpha <- (2-p)/(p-1) # haslett and parnell

  pb <- utils::txtProgressBar(min = 1, max = MC, style = 3,width=40,char = '<>') # create a progress bar

  for(n in 1:MC){
    utils::setTxtProgressBar(pb, n) # set the progress bar
    ##-----------------------------------------------------------------------------
    ## calculate some secondary model paremeters
    ## based on Haslett and Parnell (2008)
    lambda <- (mu^(2-p))/(psi*(2-p))
    beta <- 1/(psi*(p-1)*mu^(p-1))
    ##-----------------------------------------------------------------------------
    ## choose some random positions from a uniform distribution for each dated horizon.
    currPositions <- runif(nSamples, masterPositions - nThicknesses, masterPositions + nThicknesses)
    do <- order(currPositions)
    diffPositions <- diff(currPositions[do])
    thetas[do] <- sort(thetas,decreasing = T)
    positionStore[n,] <- currPositions
    ##-----------------------------------------------------------------------------
    # Interpolate between the current thetas
    for(j in 1:(nSamples - 1)){
      inRange <- predictPositions >= currPositions[do[j]] & predictPositions <= currPositions[do[j+1]]
      predval <- predictInterp(alpha = alpha,
                               lambda = lambda,
                               beta = beta,
                               predictPositions = predictPositions[inRange],
                               diffPositionj = diffPositions[j],
                               currPositionsj = currPositions[do[j]],
                               currPositionsjp1 = currPositions[do[j+1]],
                               thetaj = thetas[do[j]],
                               thetajp1 = thetas[do[j+1]])
      modelStore[inRange,n] <- predval
    }

    # Extrapolate up
    if(any(predictPositions >= currPositions[do[nSamples]])){
      inRange <- predictPositions >= currPositions[do[nSamples]]
      predval <- predictExtrapUp(alpha = alpha,
                                 lambda = lambda,
                                 beta = beta,
                                 predictPositions = predictPositions[inRange],
                                 theta1 = thetas[do[nSamples]],
                                 currPositions1 = currPositions[do[nSamples]],
                                 maxExtrap = extrapUp,
                                 extractDate = truncateUp)
      modelStore[inRange,n] <- predval
    }

    # extrapolate down
    if(any(predictPositions <= currPositions[do[1]])){
      inRange <- predictPositions <= currPositions[do[1]]
      predval <- predictExtrapDown(alpha = alpha,
                                   lambda = lambda,
                                   beta = beta,
                                   predictPositions = predictPositions[inRange],
                                   theta1 = thetas[do[1]],
                                   currPositions1 = currPositions[do[1]],
                                   maxExtrap = extrapDown,
                                   extractDate = truncateDown)
      modelStore[inRange,n] <- predval
    }



    #-----------------------------------------------------------------------------
    ## Update thetas
    for(i in 1:nSamples){
      thetaCurrent <- thetas[do[i]]
      ## choose a new theta from a truncated normal where truncations are the surrounding thetas
      thetaProposed <- truncatedWalk(old = thetaCurrent,
                                     sd = mhSD[do[i]],
                                     low = ifelse(i == nSamples,truncateUp,thetas[do[i+1]]-1e-10), # if we're at the top truncate at a really small number
                                     high = ifelse(i == 1,1e10,thetas[do[i-1]] + 1e-10)) # if we're at the bottom, truncate at a really big number

      pProposed <- log(prob[,do[i]][which.min(abs(ageGrid - thetaProposed$new))])
      pProposed <-  max(pProposed,-1000000, na.rm = T) # just incase things get really small
      pCurrent <- log(prob[,do[i]][which.min(abs(ageGrid - thetaCurrent))])
      pCurrent <-  max(pCurrent,-1000000,na.rm = T) # just incase things get really small

      priorProposed <-  ifelse(i == 1,0,
                               log(dtweediep1(thetas[do[i - 1]]-thetaProposed$new,
                                              p = p,
                                              mu = mu*(diffPositions[i-1]),
                                              phi = psi*(diffPositions[i-1])^(p-1)))) +
        ifelse(i == nSamples,0,
               log(dtweediep1(thetaProposed$new - thetas[do[i + 1]],
                              p = p,
                              mu = mu*(diffPositions[i]),
                              phi = psi*(diffPositions[i])^(p-1))))

      priorCurrent <-  ifelse(i == 1,0,
                              log(dtweediep1(thetas[do[i - 1]]-thetaCurrent,
                                             p = p,
                                             mu = mu*(diffPositions[i-1]),
                                             phi = psi*(diffPositions[i-1])^(p-1)))) +
        ifelse(i == nSamples,0,
               log(dtweediep1(thetaCurrent - thetas[do[i + 1]],
                              p = p,
                              mu = mu*(diffPositions[i]),
                              phi = psi*(diffPositions[i])^(p-1))))
      priorCurrent <- max(priorCurrent,-1000000,na.rm = T)
      priorProposed <- max(priorProposed,-1000000,na.rm = T)
      logRtheta <- priorProposed - priorCurrent + pProposed - pCurrent + log(thetaProposed$rat)
      logRtheta <- max(logRtheta,-1000000,na.rm = T)
      if(runif(1,0,1) < min(1,exp(logRtheta))){thetas[do[i]] <- thetaProposed$new}
    }
    thetaStore[n,] <- thetas
    #-----------------------------------------------------------------------------
    ## Update mu
    muCurrent <- mu
    muProposed <- truncatedWalk(old = mu, sd = muSD,low = 1e-10,high = 1e10)
    priorMu <- sum(log(dtweediep1(abs(diff(thetas[do])),
                                  p = p,
                                  mu = muProposed$new*diffPositions,
                                  phi = psi/diffPositions^(p-1)))) -
      sum(log(dtweediep1(abs(diff(thetas[do])),
                         p = p,
                         mu = mu*diffPositions,
                         phi = psi/diffPositions^(p-1)))) + log(muProposed$rat)

    mu <- ifelse(runif(1,0,1) < min(1,exp(priorMu)),muProposed$new,muCurrent)
    muStore[n] <- mu

    ##-----------------------------------------------------------------------------
    ## Update psi
    psiCurrent <- psi

    psiProposed <- truncatedWalk(old = psi,
                                 sd = psiSD,
                                 low = 1e-10,
                                 high = 1e10)

    priorPsi <- sum(log(dtweediep1(abs(diff(thetas[do])),
                                   p = p,
                                   mu = mu*diffPositions,
                                   phi = psiProposed$new/diffPositions^(p-1)))) -
      sum(log(dtweediep1(abs(diff(thetas[do])),
                         p = p,
                         mu = mu*diffPositions,
                         phi = psi/(diffPositions^(p-1))))) + log(psiProposed$rat)

    psi <- ifelse(runif(1,0,1) < min(1,exp(priorPsi)),psiProposed$new,psiCurrent)
    psiStore[n] <- psi


    ##-----------------------------------------------------------------------------
    ## Update SDs
    if(adapt == T){
      h = 200
      if(n%%h == 0){
        cd <- 2.4/sqrt(1) # Gelman et al. (1996)
        psiK <- psiStore[(n-h+1):n] - mean(psiStore[(n-h+1):n])
        psiRt <- var(psiK)
        psiSD <- sqrt((cd^2) * psiRt)
        psiSD <- ifelse(is.na(psiSD),cd*sd(psiStore,na.rm = T),psiSD)
        psiSD <- ifelse(psiSD == 0,cd*sd(psiStore,na.rm = T),psiSD)
        muK <-  muStore[(n-h+1):n] - mean(muStore[(n-h+1):n])
        muRt <- var(muK)
        muSD <- sqrt(cd^2 * muRt)
        muSD <- ifelse(is.na(muSD),cd*sd(muStore,na.rm = T),muSD)
        muSD <- ifelse(muSD == 0,cd*sd(muStore,na.rm = T),muSD)

        for(q in 1:nSamples){
          thetaK <- thetaStore[(n-h+1):n,q] - mean(thetaStore[(n-h+1):n,q])
          thetaRt <- var(thetaK)
          mhSD[q] <- ifelse(is.na(sqrt(cd^2 * thetaRt)),cd*sd(thetaStore[,q],na.rm = T),sqrt(cd^2 * thetaRt))
        }
      }
      psiSDStore[n] <- psiSD
      muSDStore[n] <- muSD
      mhSDStore[n,] <- mhSD
    }
  }
  ##-----------------------------------------------------------------------------

  return(list(confInt = apply(modelStore[,burn:MC],1,quantile,c(0.025,.5,.975)),
              model = modelStore[,burn:MC],
              thetas = thetaStore,
              positionStore = positionStore,
              psi = psiStore,
              mu = muStore,
              predictPositions = predictPositions,
              ageGrid = ageGrid,
              likelihoods = prob,
              nAges = nAges,
              masterPositions = masterPositions,
              ids = nNames,
              psiSDStore = psiSDStore,
              muSDStore = muSDStore,
              mhSDStore = mhSDStore,
              burn = burn,
              MC = MC))
}

