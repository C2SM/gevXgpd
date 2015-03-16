"conf.bounds.xval.lprof" <-
function(xval,ret=c(5,10,50,100,500)/xval$lamda,
         probs=c(0.025,0.975),profs.out=FALSE) {

      # avoid calculation of likelihood profile for modified likelihood
      # and l-moments.
      if ((xval$estim == "lmom") | (xval$estim == "pmlik")) {
        stop("ERROR: Likelihood profile is meaningful for maximum likelihood estimation only.")
      }

      # GUMBEL or GEV confidence?
      dist <- xval$dist

      # start value for optim
      if (dist == "GEV") {
         strt <- c(xval$chi,xval$alpha,xval$k) }
      if (dist == "GUMBEL") {
         strt <- c(xval$chi,xval$alpha) }

      # extract data
      dd <- as.data.frame(xval$data)[,1]

      if (dist == "GEV") {
         likfun <- "gev.lik" }
      if (dist == "GUMBEL") {
         likfun <- "gumbel.lik" }

      # calculate profile for the three parameters
      prof.chi <- lik.prof(data=dd,
                           LikFun=likfun,index=1,
                           start=strt,probs=probs,
                           frac.incr=1/10)
      prof.alpha <- lik.prof(data=dd,
                           LikFun=likfun,index=2,
                           start=strt,probs=probs,
                           frac.incr=1/10)
      if (dist == "GEV") {
      prof.k <- lik.prof(data=dd,
                           LikFun=likfun,index=3,
                           start=strt,probs=probs,
                           frac.incr=1/10)
      }
      # collect the results
      if (dist == "GEV") {
        conf.paras <- rbind(prof.chi$conf.par,
                            prof.alpha$conf.par,
                            prof.k$conf.par)
        row.names(conf.paras) <- c("chi","alpha","k")
      }
      if (dist == "GUMBEL") {
        conf.paras <- rbind(prof.chi$conf.par,
                            prof.alpha$conf.par)
        row.names(conf.paras) <- c("chi","alpha")
      }

      # calculate profile for return values at return periods
      lamda <- xval$lamda
      conf.retv <- matrix(rep(0.0,length(ret)*length(probs)),
                          ncol=length(probs))
      for (i in seq(1,length(ret))) {
          rp <- ret[i]
          z.strt <- XGEV(F.of.T(rp,lamda=lamda),
                         alpha=xval$alpha,chi=xval$chi,k=xval$k)
          if (dist == "GEV") {
            if (xval$k < -0.25) { 
               finc <- 1/50 
            } else {
               finc <- 1/10
            }
            strt <- c(z.strt,xval$alpha,xval$k)
            prof.z <- lik.prof(data=dd,
                             LikFun=gev.lik.retv,index=1,
                             start=strt,probs=probs,
                             frac.incr=finc,ret=rp*lamda)
          }
          if (dist == "GUMBEL") {
            strt <- c(z.strt,xval$alpha)
            prof.z <- lik.prof(data=dd,
                             LikFun=gumbel.lik.retv,index=1,
                             start=strt,probs=probs,
                             frac.incr=1/10,ret=rp*lamda)
          }
          conf.retv[i,] <- prof.z$conf.par
      }
      row.names(conf.retv) <- as.character(ret)
      dimnames(conf.retv)[[2]] <- as.character(probs)

      # result data object
      conf <- list(conf.retv=conf.retv,
                   conf.paras=conf.paras,
                   lamda=xval$lamda,
                   cal="LPROF") 
      if (profs.out) {               # add profile output for parameters
         if (dist == "GEV") {
           conf <- c(conf,
                   list(
                      prof.chi = prof.chi$profile,
                      prof.alpha = prof.alpha$profile,
                      prof.k = prof.k$profile)
                   ) }
         if (dist == "GUMBEL") {
           conf <- c(conf,
                   list(
                      prof.chi = prof.chi$profile,
                      prof.alpha = prof.alpha$profile)
                   ) }
      } 
      return(conf)
}

