"conf.bounds.xval.sim" <-
function(xval,ret=c(5,10,50,100,500)/xval$lamda,
         probs=c(0.025,0.975), size=500) {
# ==================================================================   
# Calculates confidence bounds for xval object.

  kthresh <- 0.000001

  # Preparation
  est <- matrix(0,ncol=length(ret),nrow=size)
  ppp <- data.frame("chi"=rep(NA,times=size), # simul results for parameters
                    "alpha"=rep(NA,times=size),
                    "k"=rep(NA,times=size))
  len <- dim(xval$data)[1]
  estim <- xval$estim	
  lamda <- xval$lamda
		      
  # Estimate from random samples
  if (estim == "mlik") {                 # case for mlik
  for (i in seq(1,size,1)) {
       ddd <- RGEV(xval$alpha,xval$chi,xval$k,len)
       objct <- fitGEV(data.frame(GAGA=ddd),
	                      lamda=xval$lamda,
	                      ret=ret,
			      dist=xval$dist,
			      estim=xval$estim,
			      method=xval$method,
			      maxit=xval$maxit)
       est[i,] <- objct$fit[,"fitted"]
       ppp[i,"chi"] <- objct$chi
       ppp[i,"alpha"] <- objct$alpha
       ppp[i,"k"] <- objct$k
  }
  } else {                      # case for lmom
  for (i in seq(1,size,1)) {
       ddd <- RGEV(xval$alpha,xval$chi,xval$k,len)
       objct <- fitGEV(data.frame(GAGA=ddd),
	   	              lamda=xval$lamda,
	                      ret=ret,
			      dist=xval$dist,
			      estim=xval$estim)
       est[i,] <- objct$fit[,"fitted"]
       ppp[i,"chi"] <- objct$chi
       ppp[i,"alpha"] <- objct$alpha
       ppp[i,"k"] <- objct$k
  }
  }
   
  # calculate quantiles for return value simulations
  conf <- matrix(rep(0.0,length(ret)*length(probs)),
                  ncol=length(probs))
  for (i in seq(1,length(ret))) {
      conf[i,] <- quantile(est[,i],probs=probs)
  }
  row.names(conf) <- as.character(ret)
  dimnames(conf)[[2]] <- as.character(probs)

  # calculate quantiles for parameters
  conf.paras <- matrix(rep(NA,3*length(probs)),
                  ncol=length(probs))
  for (i in (1:length(names(ppp)))) {
      conf.paras[i,] <- quantile(ppp[,i],probs=probs)
  }
  row.names(conf.paras) <- names(ppp)
  dimnames(conf.paras)[[2]] <- as.character(probs)

  list(conf.retv=conf,
       conf.paras=conf.paras,
       paras.sim=ppp,
       lamda=xval$lamda,
       cal="SIM")
}

