"conf.bounds.xval.mle" <-
function(xval, ret=c(5,10,50,100,500)/xval$lamda,
               probs=c(0.025,0.975)) {
# ==================================================================   

# GUMBEL or GEV confidence?
  dist <- xval$dist
  
# parameter estimates
  alpha <- xval$alpha
  chi <- xval$chi
  k <- xval$k
  lamda <- xval$lamda
    
# FUNCTION DEFINITIONS
   foft <- function (t) {1-1/(lamda*t)}

# check if mlik
  if((xval$estim != "mlik") & (xval$estim != "pmlik")) {
      stop("** ERROR ** estimator must be Maximum Likelihood")}
			     
# function definition for gradient vector
  gv <- function(ret) {
      yt <- (-log(foft(ret)))
      ytk <- yt^k
      if (dist == "GEV") {
	  c(1,(1-ytk)/k,-((1-ytk)/k+(ytk*log(yt)))*alpha/k)
      } else {
	  c(1,-log(yt))   
      }
  }

# covariance matrix
  xcov <- xval$cov
  
# loop over all ret and calculate confidence
  conf <- matrix(rep(NA,length(ret)*length(probs)),
                  ncol=length(probs))
  for (i in seq(1,length(ret))) {
      ggvv <- gv(ret[i])
      x0 <- XGEV(foft(ret[i]),alpha,chi,k)
      conf[i,] <- x0+qnorm(probs)*
            sqrt((t(as.matrix(ggvv)) %*% xcov %*% as.matrix(ggvv)))
  }
  row.names(conf) <- as.character(ret)
  dimnames(conf)[[2]] <- as.character(probs)

  # calculate quantiles for parameters (chi, alpha, k)
  conf.paras <- matrix(rep(NA,dim(xcov)[1]*length(probs)),
                  ncol=length(probs))
  if (dist == "GUMBEL") {
      est.paras <- c(chi,alpha)
      for (i in (1:dim(conf.paras)[1])) {
          conf.paras[i,] <- est.paras[i]+qnorm(probs)*sqrt(xcov[i,i])
      }
      row.names(conf.paras) <- c("chi","alpha")
      dimnames(conf.paras)[[2]] <- as.character(probs)
  }
  if (dist == "GEV") {
      est.paras <- c(chi,alpha,k)
      for (i in (1:dim(conf.paras)[1])) {
          conf.paras[i,] <- est.paras[i]+qnorm(probs)*sqrt(xcov[i,i])
      }
      row.names(conf.paras) <- c("chi","alpha","k")
      dimnames(conf.paras)[[2]] <- as.character(probs)
  }
  
  list(conf.retv=conf,
       conf.paras=conf.paras,
       paras.cov=xcov,
       lamda=xval$lamda,
       cal="MLE")  
}

