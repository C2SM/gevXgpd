"gev.lik" <-
function(pars,data) {
        # computes neg log lik of gev model including gumbel limit
        
        kthresh <- 0.000001     # threshold to distinguish GUMBEL
        len <- length(data)     # sample size

        chi   <- pars[1]
        alpha <- pars[2]
        kk    <- pars[3]
        y <- (data-chi)/alpha
        if (abs(kk)>kthresh) {           # gev
           y <- (1 - kk * y)
           if(any(y <= 0) || (alpha <= 0)) return(10^9) # very large
           len*log(alpha) + sum(y^(1/kk)) + sum(log(y) * (1 - 1/kk))
        } else {                         # gumbel
           len*log(alpha) + sum(y) + sum(exp(-y))        
        }
}

