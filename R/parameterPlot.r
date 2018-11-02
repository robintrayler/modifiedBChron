#' Make plots of model parameters from the output of an ageModel run
#' @param model Output of the \code{ageModel} function
#' @param PDF should plots be printed to the current working directory as PDF files? Defaults to FALSE
#' @export

parameterPlot <- function(model,PDF = F){
  layout(cbind(c(1,2,3)))
  xlim = c(max(0,mean(model$mu[model$burn:model$MC]) - 4*sd(model$mu[model$burn:model$MC])),mean(model$mu[model$burn:model$MC]) + 4*sd(model$mu[model$burn:model$MC]))
  colsPal <- colorRampPalette(c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6'))
  cols <- colsPal(ncol(model$thetas))

  hist(model$mu,
       breaks = seq(min(model$mu)-0.01,max(model$mu) + 0.01,
                    by = mean(diff(seq(min(model$mu[model$burn:model$MC]),max(model$mu[model$burn:model$MC]),length = 150)))),
       xlab = expression(mu),
       col = 'grey',
       xlim = xlim,
       border = NA,
       prob = T,
       yaxs = 'i',
       xaxs = 'i',
       main = expression(mu))
  box()
  abline(v = quantile(model$mu,.5),
         col = 'tomato',
         lwd = 2)

  ylim = c(max(0,mean(model$mu) - 4*sd(model$mu)),mean(model$mu) + 8*sd(model$mu))
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


  ylim = c(max(0,mean(model$muSDStore) - 4*sd(model$muSDStore)),mean(model$muSDStore) + 8*sd(model$muSDStore))
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


  ##-----------------------------------------------------------------------------------------------
  ## psi plots
  # dev.new()
  layout(cbind(c(1,2,3)))
  xlim = c(max(0,mean(model$psi)-sd(model$psi)*4),mean(model$psi)+sd(model$psi)*8)
  hist(model$psi,breaks = length(model$psi)/150,
       xlab = expression(psi),
       col = 'grey',
       xlim = xlim,
       border = NA,
       prob = T,
       yaxs = 'i',
       xaxs = 'i',
       main = expression(psi))
  box()
  abline(v = quantile(model$psi,.5),
         col = 'tomato',
         lwd = 2)
  ylim = c(max(0,mean(model$psi)-sd(model$psi)*4),mean(model$psi)+sd(model$psi)*8)
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

  ylim = c(max(0,mean(model$psiSDStore) - 4*sd(model$psiSDStore)),mean(model$psiSDStore) + 8*sd(model$psiSDStore))
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
  abline(h = quantile(model$psiSDStore,.5),
         col = 'tomato',
         lwd = 2)

  ##-----------------------------------------------------------------------------------------------
  ## Theta plots
  layout(c(1,2,3))
  for(i in 1:ncol(model$thetas)){
    layout(cbind(c(1,2,3)))
    xlim = c(mean(model$thetas[,i]) - 4*sd(model$thetas[,i]),mean(model$thetas[,i]) + 4*sd(model$thetas[,i]))
    hist(model$thetas[,i],breaks = length(model$thetas[,i])/150,
         xlab = expression(theta),
         col = cols[i],
         xlim = xlim,
         border = NA,
         prob = T,
         yaxs = 'i',
         xaxs = 'i',
         main = model$ids[i])
    box()
    abline(v = quantile(model$thetas[,i],.5),
           col = 'tomato',
           lwd = 2)

    ylim = c(max(0,mean(model$thetas[,i]) - 4*sd(model$thetas[,i])),
             mean(model$thetas[,i]) + 8*sd(model$thetas[,i]))
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

    ylim = c(max(0,mean(model$mhSDStore[,i]) - 4*sd(model$mhSDStore[,i])),
             mean(model$mhSDStore[,i]) + 8*sd(model$mhSDStore[,i]))
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
  ################################################################################
  ################################################################################
  ## If PDFs are required.
  if(PDF == T){
    pdf(file = 'parameterPlots.pdf',
        width = 10,
        height = 10)
    # dev.new()
    layout(cbind(c(1,2,3)))
    xlim = c(mean(model$mu) - 4*sd(model$mu),mean(model$mu) + 4*sd(model$mu))
    hist(model$mu,breaks = length(model$mu)/150,
         xlab = expression(mu),
         col = 'grey',
         xlim = xlim,
         border = NA,
         prob = T,
         yaxs = 'i',
         xaxs = 'i',
         main = expression(mu))
    box()
    abline(v = quantile(model$mu,.5),
           col = 'tomato',
           lwd = 2)

    ylim = c(max(0,mean(model$mu) - 4*sd(model$mu)),mean(model$mu) + 8*sd(model$mu))
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


    ylim = c(max(0,mean(model$muSDStore) - 4*sd(model$muSDStore)),mean(model$muSDStore) + 8*sd(model$muSDStore))
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


    ##-----------------------------------------------------------------------------------------------
    ## psi plots
    # dev.new()
    layout(cbind(c(1,2,3)))
    xlim = c(max(0,mean(model$psi)-sd(model$psi)*4),mean(model$psi)+sd(model$psi)*8)
    hist(model$psi,breaks = length(model$psi)/150,
         xlab = expression(psi),
         col = 'grey',
         xlim = xlim,
         border = NA,
         prob = T,
         yaxs = 'i',
         xaxs = 'i',
         main = expression(psi))
    box()
    abline(v = quantile(model$psi,.5),
           col = 'tomato',
           lwd = 2)
    ylim = c(max(0,mean(model$psi)-sd(model$psi)*4),mean(model$psi)+sd(model$psi)*8)
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

    ylim = c(max(0,mean(model$psiSDStore) - 4*sd(model$psiSDStore)),mean(model$psiSDStore) + 8*sd(model$psiSDStore))
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
    abline(h = quantile(model$psiSDStore,.5),
           col = 'tomato',
           lwd = 2)

    ##-----------------------------------------------------------------------------------------------
    ## Theta plots
    layout(c(1,2,3))
    for(i in 1:ncol(model$thetas)){
      layout(cbind(c(1,2,3)))
      xlim = c(mean(model$thetas[,i]) - 4*sd(model$thetas[,i]),mean(model$thetas[,i]) + 4*sd(model$thetas[,i]))
      hist(model$thetas[,i],breaks = length(model$thetas[,i])/150,
           xlab = expression(theta),
           col = cols[i],
           xlim = xlim,
           border = NA,
           prob = T,
           yaxs = 'i',
           xaxs = 'i',
           main = model$ids[i])
      box()
      abline(v = quantile(model$thetas[,i],.5),
             col = 'tomato',
             lwd = 2)

      ylim = c(max(0,mean(model$thetas[,i]) - 4*sd(model$thetas[,i])),
               mean(model$thetas[,i]) + 8*sd(model$thetas[,i]))
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

      ylim = c(max(0,mean(model$mhSDStore[,i]) - 4*sd(model$mhSDStore[,i])),
               mean(model$mhSDStore[,i]) + 8*sd(model$mhSDStore[,i]))
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
    dev.off()
  }

}

