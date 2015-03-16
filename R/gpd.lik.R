"gpd.lik" <-
function(pars,data,chi=0,tim=length(data)) {
        # computes neg log lik of gpd model including the pdf 
        # of peak frequency
        
        kthresh <- 0.000001     #threshold to distinguish EXPON
        len <- length(data)     # sample size

        lamda  <- pars[1]
        alpha <- pars[2]
        kk    <- pars[3]

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

        #add part from peak frequency
          # expand likelihood function to account for sample size distribution
          # see theoretical details in Coles p.82. 
          # Here we deviate from the approach 
          # in Coles in the sense that lamda (numb of exceedances per 
          # time unit) is considered as the third parameter. Moreover the 
          # the variance of lamda is modelled assuming a poisson model for the
          # exceedances rather than a binomial model. This avoids knowledge of
          # the full sample size (number of total observations) which is not
          # really known anyway if the exceedances were determined after 
          # declustering. The assumption of the poisson model is just that the
          # exceedances are very rare, which is probably a good assumption for
          # for high thresholds anyway. The formuli are in my notes.  
        if (lamda <= 0) {
           ll.poisson <- 10^9 
        } else {
           ll.poisson <- -log(dpois(x=len,lambda=lamda*tim))
        }
        ll <- ll.gpd+ll.poisson

        return(ll)
}

