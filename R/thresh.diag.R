"thresh.diag" <-
function(pks,threshs,probs=c(0.025,0.5,0.975),min.exceeds=5) {

  # prepare result matrices
  res.numex <- rep(0,times=length(threshs)) 
  res.scale.star <- matrix(NA,ncol=length(probs),nrow=length(threshs)) 
  res.shape <- matrix(NA,ncol=length(probs),nrow=length(threshs)) 
  res.mnex <- matrix(NA,ncol=length(probs),nrow=length(threshs)) 

  # loop over thresholds
  for (u in threshs) {

     # extract exceedances
     pksu <- pks[(pks[,1] > u),1]

     # don't do calculations if less than min.exceeds peaks
     if (length(pksu)<min.exceeds) {next}

     # estimate fit
     pot <- fitGPD(data.frame(peaks=pksu),lamda=1,chi=u)

     # estimate mle confidence bounds
     conf <- conf.bounds.pot(pot,probs=probs,cal="MLE")

     # extract results for shape
     ii <- which(u == threshs)
     res.shape[ii,] <- conf$conf.paras["k",]

     # safe number of exceedances
     res.numex[ii] <- length(pksu)
     
     # extract results for scale star = sigma+k*chi
     gv.scale.star <- c(0,1,u)       # partial derivs of scale.star
     var.scale.star <-  t(gv.scale.star) %*% conf$paras.cov %*% gv.scale.star
     scale.star <- pot$alpha+pot$k*u
     res.scale.star[ii,] <- scale.star+qnorm(probs)*sqrt(var.scale.star)

     # calculate results for mean exceedance 
     # here the results are not based on the standard errors of the fit,
     # they are simply determined from the approx normality of sample means.
     mnex <- mean(pksu-u)
     sigma.mnex <- sd(pksu-u)/sqrt(length(pksu))
     res.mnex[ii,] <- mnex+qnorm(probs)*sigma.mnex
  }

  # add names to result matrices
  row.names(res.scale.star) <- as.character(threshs)
  dimnames(res.scale.star)[[2]] <- as.character(probs)
  row.names(res.shape) <- as.character(threshs)
  dimnames(res.shape)[[2]] <- as.character(probs)
  row.names(res.mnex) <- as.character(threshs)
  dimnames(res.mnex)[[2]] <- as.character(probs)

  list(
    threshs = threshs,
    num.exceed = res.numex,
    q.scale.star = res.scale.star,
    q.shape = res.shape,
    q.mean.exceed = res.mnex)
}
