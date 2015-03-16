"conf.bounds.pot.lprof" <-
function(pot,ret=c(5,10,50,100,500),
         probs=c(0.025,0.975),profs.out=FALSE) {

      # avoid calculation of likelihood profile for modified likelihood
      # and l-moments.
      if ((pot$estim == "lmom") | (pot$estim == "pmlik")) {
        stop("ERROR: Likelihood profile is meaningful for maximum likelihood estimation only.")
      }

      # EXPON or GPD confidence?
      dist <- pot$dist

      #preps
      chi <- pot$chi          # get threshold
      len <- length(pot$data[,1])
      tim <- len/pot$lamda

      # start value for optim
      if (dist == "GPD") {
         strt <- c(pot$lamda,pot$alpha,pot$k) }
      if (dist == "EXPON") {
         strt <- c(pot$lamda,pot$alpha) }

      #extract data
      dd <- as.data.frame(pot$data)[,1]

      if (dist == "GPD") {
         likfun <- "gpd.lik" }
      if (dist == "EXPON") {
         likfun <- "expon.lik" }

      # calculate profile for the three parameters
      prof.lamda <- lik.prof(data=dd,
                           LikFun=likfun,index=1,
                           start=strt,probs=probs,
                           frac.incr=1/10,chi=chi,tim=tim)
      prof.alpha <- lik.prof(data=dd,
                           LikFun=likfun,index=2,
                           start=strt,probs=probs,
                           frac.incr=1/10,chi=chi,tim=tim)
      if (dist == "GPD") {
      prof.k <- lik.prof(data=dd,
                           LikFun=gpd.lik,index=3,
                           start=strt,probs=probs,
                           frac.incr=1/10,chi=chi,tim=tim)
      }
      #collect the results
      if (dist == "GPD") {
        conf.paras <- rbind(prof.lamda$conf.par,
                          prof.alpha$conf.par,
                          prof.k$conf.par)
        row.names(conf.paras) <- c("lamda","alpha","k")
      }
      if (dist == "EXPON") {
        conf.paras <- rbind(prof.lamda$conf.par,
                          prof.alpha$conf.par)
        row.names(conf.paras) <- c("lamda","alpha")
      }

      # calculate profile for return values at return periods
      lamda <- pot$lamda
      conf.retv <- matrix(rep(0.0,length(ret)*length(probs)),
                          ncol=length(probs))
      for (i in seq(1,length(ret))) {
          rp <- ret[i]
          z.strt <- XGPD(F.of.T(rp,lamda=lamda),
                         alpha=pot$alpha,chi=pot$chi,k=pot$k)
          if (dist == "GPD") {
            if (pot$k < -0.25) { 
               finc <- 1/50 
            } else {
               finc <- 1/10
            }
            strt <- c(pot$lamda,z.strt,pot$k)
            prof.z <- lik.prof(data=dd,
                           LikFun=gpd.lik.retv,index=2,
                           start=strt,probs=probs,
                           frac.incr=finc,chi=chi,tim=tim,ret=lamda*rp)
          }
          if (dist == "EXPON") {
            strt <- c(pot$lamda,z.strt)
            prof.z <- lik.prof(data=dd,
                           LikFun=expon.lik.retv,index=2,
                           start=strt,probs=probs,
                           frac.incr=1/10,chi=chi,tim=tim,ret=lamda*rp)
          }
          conf.retv[i,] <- prof.z$conf.par
      }
      row.names(conf.retv) <- as.character(ret)
      dimnames(conf.retv)[[2]] <- as.character(probs)

      #result data object
      conf <- list(conf.retv=conf.retv,
                   conf.paras=conf.paras,
                   lamda=pot$lamda,
                   cal="LPROF") 
      if (profs.out) {               #add profile output for parameters
         if (dist == "GPD") {
           conf <- c(conf,
                   list(
                      prof.lamda = prof.lamda$profile,
                      prof.alpha = prof.alpha$profile,
                      prof.k = prof.k$profile)
                   ) }
         if (dist == "EXPON") {
           conf <- c(conf,
                   list(
                      prof.lamda = prof.lamda$profile,
                      prof.alpha = prof.alpha$profile)
                   ) }
      } 
      return(conf)
}

