"gpd.lik.retv" <- 
function(pars,data,chi=0,tim=length(data),ret) {
        # computes neg log lik of gpd model
        # with a reparameterisation in terms of exceedance frequency,
        # return level, and shape parameter (in this order in pars).
        # ret is the return period in units of the average interval
        # between threshold exceedances. 
        
        kthresh <- 0.000001     # threshold to distinguish EXPON
        len <- length(data)     # sample size

        lamda <- pars[1]
        zz <- pars[2]
        kk    <- pars[3]

        # calculate alpha from the other parameters
        if (abs(kk)>kthresh) {           # gpd
           alpha <- (zz-chi)*kk/(1-(1/ret)^kk)
        } else {                         # expon
           alpha <- (zz-chi)/log(ret)
        }

        # part from gpd only
        y <- (data-chi)/alpha
        if (abs(kk)>kthresh) {           # gpd
           y <- (1 - kk * y)
           if(any(y <= 0) || (alpha <= 0)) {
               ll.gpd <- 10^9    # very large
           } else {
               ll.gpd <- len*log(alpha) + (1-1/kk)*sum(log(y))
           }
        } else {                         # expon
               ll.gpd <- len*log(alpha) + sum(y)        
        }

        # add part from peak frequency
        if (lamda <= 0) {
           ll.poisson <- 10^9 
        } else {
           ll.poisson <- -log(dpois(x=len,lambda=lamda*tim))
        }
        ll <- ll.gpd+ll.poisson

        return(ll)
}

