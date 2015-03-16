"lik.prof" <-
function(data,LikFun,index,start,probs=c(0.025,0.975),
         frac.incr=1/10, method="BFGS", maxit=10000, ...) {

     # function definition
     filt.prof <- function(prof) {
     # This function filters out multiple occurrences of the same prob column
     # in the profile data frame. It turned out that this can happen in 
     # rare cases, especially prob=0.5 can occur for several parameter values
     # near the maximum of the likelihood. Multiple values of prob will make
     # to fail the subsequent spline interpolation and invertion of the 
     # likelihood profile. This function eliminates multiple entries, only 
     # keeping that entry for which the parameter value is the closest to the
     # median of all multiple entries. 
            y <- prof[,3]
            x <- prof[,1]

            if (length(unique(y)) != length(y)) {
               yy <- unique(y)
               for (yyy in yy) {
                  ii <- which(y == yyy)
                  if(length(ii) > 1) {
                     xmed <- median(x[ii])
                     jj <- which(abs(x-xmed) == min(abs(x-xmed)))[1]
                     prof <- prof[-c(ii[ii != jj]),]
                  }
               }
            }

            prof
     }

     # preparation
     FUN <- match.fun(LikFun)

     # maximum likelihood estimate (optimize wrpt all pars)
     x <- optim(par=start, fn=FUN, data=data,
                hessian=TRUE, method=method, 
                control=list(maxit=maxit),...)
     if (x$convergence > 0) 
        stop("ERROR: iterative estimation did not converge")
     best.est <- x$par
     min.lik <- x$value
     cov.mat <- solve(x$hessian)

     # determine increment of profiling parameter
     dp <- sqrt(cov.mat[index,index])*frac.incr        

     # build function to use in optim for the profile
     # i.e. optimize constrained
     likfun <- function(free,fix,data,...) {
              ll <- length(c(free,fix))
              pars <- rep(NA,times=ll)
              pars[index] <- fix
              pars[-index] <- free
              FUN(pars=pars,data=data,...)
     }

     # critical deviance such as to reach min/max probs 
     # D = 2(l(best.estimate) - l(profile)) is chi-squared with df=1
     deviance.prob <- 1-2*min(c(probs,1-probs))
     crit.devnc <- qchisq(p=deviance.prob,df=1)

     # start profile dataset
     prof <- matrix(c(best.est[index],min.lik,0.5),ncol=3,nrow=1)
     dimnames(prof)[[2]] <- c("par","neg.loglik","prob")

     # loop for generating the profile
     for (dpp in c(-dp,dp)) {          # to the left and right of max

       lik.val <- min.lik
       strt <- c(best.est[-index])
       p.val <- best.est[index]
       devnc <- 0.0

       while(devnc < crit.devnc) {

         # increment for p.val
         p.val <- p.val+dpp

         # optimize for fixed value of parameter
         x <- optim(par=strt, fn=likfun, data=data, fix=p.val,
                    hessian=FALSE, method=method,
                    control=list(maxit=maxit),...)
         if (x$convergence > 0) 
             stop("ERROR: iterative estimation did not converge")

         # likelihood, deviance and p-value
         lik.val <- x$value
         devnc <- 2*(lik.val-min.lik)   # note that LikFun is neg log likelh.
         pval <- 0.5 + sign(dpp)*pchisq(q=devnc,df=1)/2

         # extend profile dataset
         zz <- c(p.val,lik.val,pval)
         if (dpp < 0) {prof <- rbind(zz,prof)
         } else {prof <- rbind(prof,zz)}

         # new start value for iteration
         strt <- c(x$par)
 
       }
     }
     row.names(prof) <- rep("",length=dim(prof)[1])

     # filter profile from multiple values
     prof <- filt.prof(prof)

     # interpolate at selected probs
     sf <- splinefun(x=prof[,3],y=prof[,1])
     conf.par <- matrix(sf(x=probs),nrow=1)
     dimnames(conf.par)[[2]] <- probs

     # return result
     return(
     list(
          conf.par=conf.par,
          profile=prof)
     )
}

