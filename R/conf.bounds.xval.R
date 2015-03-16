"conf.bounds.xval" <-
function(xval, ret=c(5,10,50,100,500)/xval$lamda,
               probs=c(0.025,0.975),
               size=500,cal="MLE",
               profs.out=FALSE) {
# ==================================================================   

# special case when no convergence (parameters NA, return NA)
if (is.na(xval$k)) {
  print("WARNING, conf.bounds.pot: could not calculate confidence.")
  conf <- matrix(NA,nrow=length(ret),ncol=length(probs))
  row.names(conf) <- as.character(ret)
  dimnames(conf)[[2]] <- as.character(probs)
  conf.paras <- matrix(rep(NA,3*length(probs)),
                             ncol=length(probs))
  row.names(conf.paras) <- c("nevts","alpha","k")
  dimnames(conf.paras)[[2]] <- as.character(probs)
  paras.cov <- matrix(NA,nrow=3,ncol=3)
  list(conf.retv=conf,conf.paras=conf.paras,
       paras.cov=paras.cov,lamda=xval$lamda)
} else {
  switch(cal,
    "SIM"  = conf.bounds.xval.sim(xval,ret,probs,size=size),
    "MLE"  = conf.bounds.xval.mle(xval,ret,probs),
    "LPROF"  = conf.bounds.xval.lprof(xval,ret,probs,profs.out=profs.out),
    stop("** ERROR ** (conf.bounds.xval) cal must be one of SIM, MLE, LPROF"))
}		 
}

