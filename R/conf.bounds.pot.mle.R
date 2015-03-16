"conf.bounds.pot.mle" <- 
function(pot, ret=c(5,10,50,100,500),
         probs=c(0.025,0.975)) {
# ==================================================================   

  # EXPON or GPD confidence?
  dist <- pot$dist
  
  # parameter estimates
  alpha <- pot$alpha
  chi <- pot$chi
  k <- pot$k
  lamda <- pot$lamda
  sampsize <- dim(pot$data)[1]

  # FUNCTION DEFINITIONS
  foft <- function (t) {1-1/(lamda*t)}

  # check if estimation was with mlik
  if(pot$estim != "mlik") {
      stop("** ERROR ** estimator must be Maximum Likelihood")}
			     
  # function definition for gradient vector
  gv <- function(ret) {
      yt <- (lamda*ret)
      ytk <- yt^(-k)
      if (dist == "GPD") {
          c(alpha*ytk/lamda,(1-ytk)/k,-alpha*((1-ytk)/k-log(yt)*ytk)/k)
      } else {
          c(alpha/lamda,log(yt))   
      }
  }
  
  # variance covariance matrix of parameters
  xcov <- pot$cov
  
  # loop over all ret and calculate confidence for return values
  conf <- matrix(rep(NA,length(ret)*length(probs)),
                  ncol=length(probs))
  for (i in seq(1,length(ret))) {
      ggvv <- gv(ret[i])
      x0 <- XGPD(foft(ret[i]),alpha,chi,k)
      conf[i,] <- x0+qnorm(probs)*
            sqrt((t(as.matrix(ggvv)) %*% xcov %*% as.matrix(ggvv)))
  }
  row.names(conf) <- as.character(ret)
  dimnames(conf)[[2]] <- as.character(probs)

  # calculate quantiles for parameters (lamda, alpha, k)
  conf.paras <- matrix(rep(NA,dim(xcov)[1]*length(probs)),
                  ncol=length(probs))
  if (dist == "EXPON") {
      est.paras <- c(lamda,alpha)
      conf.paras[1,] <- est.paras[1]+qnorm(probs)*sqrt(xcov[1,1])
      conf.paras[2,] <- est.paras[2]+qnorm(probs)*sqrt(xcov[2,2])
      row.names(conf.paras) <- c("lamda","alpha")
      dimnames(conf.paras)[[2]] <- as.character(probs)
  }
  if (dist == "GPD") {
      est.paras <- c(lamda,alpha,k)
      for (i in (1:dim(conf.paras)[1])) {
          conf.paras[i,] <- est.paras[i]+qnorm(probs)*sqrt(xcov[i,i])
      }
      row.names(conf.paras) <- c("lamda","alpha","k")
      dimnames(conf.paras)[[2]] <- as.character(probs)
  }
  
  list(conf.retv=conf,
       conf.paras=conf.paras,
       paras.cov=xcov,
       lamda=pot$lamda,
       cal="MLE")  
}

