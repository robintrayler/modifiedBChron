#' Make plots of modeled posteriors from the output of an ageModel run
#' @param model Output of the \code{ageModel} function
#' @param prob should plots be printed to the current working directory as PDF files? Defaults to FALSE
#' @param sifF number of decimal places to display on histograms
#' @export

posteriorPlot <- function(model, prob = 0.95, sigF = 3){
  lowerP <- (1 - prob)/2
  upperP <- (1 + prob)/2
  ##----------------------------------------------------------------------------
  ## set up colors and layout
  colsPal <- colorRampPalette(c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6'))
  cols <- colsPal(ncol(model$thetas))
  layout(cbind(c(1,2,3)))

  ##----------------------------------------------------------------------------
  ## set x limits for mu plot
  xlim = c(max(0, mean(model$mu[model$burn:model$MC]) - 4*sd(model$mu[model$burn:model$MC])),
           mean(model$mu[model$burn:model$MC]) + 4*sd(model$mu[model$burn:model$MC]))
  ##----------------------------------------------------------------------------
  ## set plot mu results
  h <- hist(model$mu,
            breaks = seq(min(model$mu)-0.01,max(model$mu) + 0.01,
                         by = mean(diff(seq(min(model$mu[model$burn:model$MC]),
                                            max(model$mu[model$burn:model$MC]),
                                            length = 150)))),
            xlab = expression(mu),
            col = 'grey',
            xlim = xlim,
            border = NA,
            prob = T,
            yaxs = 'i',
            xaxs = 'i',
            main = expression(mu))
  box()
  text(x = quantile(model$mu, c(lowerP,.5,upperP)),
       y = max(h$density)*.8,
       labels = round(quantile(model$mu, c(lowerP,.5,upperP)), sigF),
       col = 'tomato', cex = 1.5,
       pos = 4)
  abline(v = quantile(model$mu,c(lowerP,.5,upperP)),
         col = 'tomato',
         lwd = 2)

  ylim = c(max(0,mean(model$mu) - 4*sd(model$mu)),
           mean(model$mu) + 8*sd(model$mu))
  plot(model$mu,
       ylim = ylim,
       type = 'l',
       xaxs = 'i',
       yaxs = 'i',
       xlab = 'MCMC iteration',
       lwd = .5,
       col = 'grey',
       ylab = expression(mu))
  polygon(x = c(0,model$burn,model$burn, 0),
          y = c(0,0,max(ylim)+0.1*max(ylim),max(ylim)+0.1*max(ylim)),
          border = NA,
          col = rgb(0,0,0,.25))
  abline(h = quantile(model$mu,.5),
         col = 'tomato',
         lwd = 2)

  ##----------------------------------------------------------------------------
  ## make mu trace plots
  ylim = c(max(0,mean(model$muSDStore) -
                 4*sd(model$muSDStore)),mean(model$muSDStore) +
             8*sd(model$muSDStore))
  ##----------------------------------------------------------------------------
  ## make mu proposal SD trace plots
  plot(model$muSDStore,
       type = 'l',
       xaxs = 'i',
       yaxs = 'i',
       lwd = 1,
       col = 'grey',
       ylab = (expression(mu~' proposal SD')),
       xlab = 'MCMC iteration',
       ylim = ylim)
  polygon(x = c(0,model$burn,model$burn, 0),
          y = c(0,0,max(ylim)+0.1*max(ylim),max(ylim)+0.1*max(ylim)),
          border = NA,
          col = rgb(0,0,0,.25))
  abline(h = quantile(model$muSDStore,.5),
         col = 'tomato',
         lwd = 2)

  ##----------------------------------------------------------------------------
  ## make psi plots
  layout(cbind(c(1,2,3)))
  xlim = c(max(0,mean(model$psi)-sd(model$psi)*4),
           mean(model$psi)+sd(model$psi)*8)
  ##----------------------------------------------------------------------------
  h <- hist(model$psi,
            breaks = length(model$psi) / 150,
            xlab = expression(psi),
            col = 'grey',
            xlim = xlim,
            border = NA,
            prob = T,
            yaxs = 'i',
            xaxs = 'i',
            main = expression(psi))
  box()
  text(x = quantile(model$psi, c(lowerP,.5,upperP)),
       y = max(h$density)*.8,
       labels = round(quantile(model$psi, c(lowerP,.5,upperP)), sigF),
       col = 'tomato', cex = 1.5,
       pos = 4)
  abline(v = quantile(model$psi,c(lowerP,.5,upperP)),
         col = 'tomato',
         lwd = 2)
  ##----------------------------------------------------------------------------
  ## make psi trace plots
  ylim = c(max(0, mean(model$psi)-sd(model$psi)*4),
           mean(model$psi)+sd(model$psi)*8)
  ##----------------------------------------------------------------------------
  plot(model$psi,
       ylim = ylim,
       type = 'l',
       xaxs = 'i',
       yaxs = 'i',
       xlab = 'MCMC iteration',
       lwd = .5,
       col = 'grey',
       ylab = expression(psi))
  polygon(x = c(0,model$burn,model$burn, 0),
          y = c(0,0,max(ylim)+0.1*max(ylim),max(ylim)+0.1*max(ylim)),
          border = NA,
          col = rgb(0,0,0,.25))
  abline(h = quantile(model$psi,.5),
         col = 'tomato',
         lwd = 2)
  ##----------------------------------------------------------------------------
  ylim = c(max(0,mean(model$psiSDStore) - 4*sd(model$psiSDStore)),
           mean(model$psiSDStore) + 8*sd(model$psiSDStore))
  ##----------------------------------------------------------------------------
  plot(model$psiSDStore,
       type = 'l',
       xaxs = 'i',
       yaxs = 'i',
       lwd = 1,
       col = 'grey',
       ylab = (expression(psi~' proposal SD')),
       xlab = 'MCMC iteration',
       ylim = ylim)
  polygon(x = c(0,model$burn,model$burn, 0),
          y = c(0,0,max(ylim)+0.1*max(ylim),max(ylim)+0.1*max(ylim)),
          border = NA,
          col = rgb(0,0,0,.25))
  ##----------------------------------------------------------------------------
  abline(h = quantile(model$psiSDStore,.5),
         col = 'tomato',
         lwd = 2)

  ##----------------------------------------------------------------------------
  ## Theta plots
  layout(c(1,2,3))
  for(i in 1:ncol(model$thetas)){
    layout(cbind(c(1,2,3)))
    xlim = c(mean(model$thetas[,i]) - 4*sd(model$thetas[,i]),
             mean(model$thetas[,i]) + 4*sd(model$thetas[,i]))
    h <-  hist(model$thetas[,i],breaks = length(model$thetas[,i])/150,
               xlab = expression(theta),
               col = cols[i],
               xlim = xlim,
               border = NA,
               prob = T,
               yaxs = 'i',
               xaxs = 'i',
               main = model$ids[i])
    box()
    abline(v = quantile(model$thetas[,i],c(lowerP,.5,upperP)),
           col = 'tomato',
           lwd = 2,
           lty = 2)
    text(x = quantile(model$thetas[, i], c(lowerP,.5,upperP)),
         y = max(h$density)*.8,
         labels = round(quantile(model$thetas[, i], c(lowerP,.5,upperP)), sigF),
         col = 'tomato', cex = 1.5,
         pos = 4)
    ##----------------------------------------------------------------------------
    ylim = c(max(0,mean(model$thetas[,i]) - 4*sd(model$thetas[,i])),
             mean(model$thetas[,i]) + 8*sd(model$thetas[,i]))
    ##----------------------------------------------------------------------------
    plot(model$thetas[,i],
         type = 'l',
         xaxs = 'i',
         yaxs = 'i',
         lwd = 1,
         col = cols[i],
         ylab = (expression(theta)),
         xlab = 'MCMC iteration',
         ylim = ylim)
    polygon(x = c(0,model$burn,model$burn, 0),
            y = c(0,0,max(ylim)+0.1*max(ylim),max(ylim)+0.1*max(ylim)),
            border = NA,
            col = rgb(0,0,0,.25))
    abline(h = quantile(model$thetas[,i],.5),
           col = 'tomato',
           lwd = 2)
    ##----------------------------------------------------------------------------
    ylim = c(max(0,mean(model$mhSDStore[,i]) - 4*sd(model$mhSDStore[,i])),
             mean(model$mhSDStore[,i]) + 8*sd(model$mhSDStore[,i]))
    ##----------------------------------------------------------------------------
    plot(model$mhSDStore[,i],
         type = 'l',
         xaxs = 'i',
         yaxs = 'i',
         lwd = 1,
         col = cols[i],
         ylab = (expression(theta ~' proposal SD')),
         xlab = 'MCMC iteration',
         ylim = ylim)
    polygon(x = c(0,model$burn,model$burn, 0),
            y = c(0,0,max(ylim)+0.1*max(ylim),max(ylim)+0.1*max(ylim)),
            border = NA,
            col = rgb(0,0,0,.25))
    abline(h = quantile(model$mhSDStore[,i],.5),
           col = 'tomato',
           lwd = 2)
  }
}
