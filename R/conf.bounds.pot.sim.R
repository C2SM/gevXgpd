"conf.bounds.pot.sim" <- 
function(pot, ret=c(5,10,50,100,500),
         probs=c(0.025,0.975), size=500) {
# ==================================================================   
# Calculates confidence bounds for pot object.

  kthresh <- 0.000001

  # Preparation
  est <- matrix(0,ncol=length(ret),nrow=size)	# simul matr for return values
  ppp <- data.frame("lamda"=rep(NA,times=size), # simul results for parameters
                    "alpha"=rep(NA,times=size),
                    "k"=rep(NA,times=size))
  
  lenex <- dim(pot$data)[1]	   # expected number of events
  estim <- pot$estim
  lamda <- pot$lamda
  nyears <- lenex/lamda

  # generate random event numbers (exceedances, Poisson-Distribution)  
  nevts <- rpois(size,lenex) 

  # Simulate random samples and get estimates of parameters and return values
  if (estim == "mlik") {                 # case for mlik
  for (i in seq(1,size,by=1)) {
       ddd <- RGPD(pot$alpha,pot$chi,pot$k,nevts[i])
       objct <- fitGPD(data.frame(GAGA=ddd),
                      chi=pot$chi,
                      lamda=lamda*nevts[i]/lenex,	   
                      ret=ret,
                      dist=pot$dist,
                      estim=pot$estim,
                      method=pot$method,
                      maxit=pot$maxit)
       est[i,] <- objct$fit[,"fitted"]
       ppp[i,"lamda"] <- objct$lamda
       ppp[i,"alpha"] <- objct$alpha
       ppp[i,"k"] <- objct$k
  }
  } else {                      # case for lmom
  for (i in seq(1,size,by=1)) {
       ddd <- RGPD(pot$alpha,pot$chi,pot$k,nevts[i])
       objct <- fitGPD(data.frame(GAGA=ddd),
                      chi=pot$chi,
                      lamda=lamda*nevts[i]/lenex,	   
                      ret=ret,
                      dist=pot$dist,
                      estim=pot$estim)
       est[i,] <- objct$fit[,"fitted"]
       ppp[i,"lamda"] <- objct$lamda
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

  # collect the results  
  list(conf.retv=conf,
       conf.paras=conf.paras,
       paras.sim=ppp,
       lamda=pot$lamda,
       cal="SIM")
}
