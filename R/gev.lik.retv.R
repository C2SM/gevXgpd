"gev.lik.retv" <-
function(pars,data,ret) {
        # computes neg log lik of gev model including gumbel limit
        # with a reparameterisation in terms of return level, scale
        # and shape parameter (this order in pars).
        # ret is the return period in units of the block size
        
        kthresh <- 0.000001     # threshold to distinguish GUMBEL
        len <- length(data)     # sample size

        zz <- pars[1]
        alpha <- pars[2]
        kk    <- pars[3]

        # calculate chi from the return level
        if (abs(kk)>kthresh) {           # gev
           chi <- zz - alpha*(1-(-log(1-1/ret))^kk)/kk
        } else {
           chi <- zz + alpha*log(-log(1-1/ret))
        }

        y <- (data-chi)/alpha
        if (abs(kk)>kthresh) {           # gev
           y <- (1 - kk * y)
           if(any(y <= 0) || (alpha <= 0)) return(10^6) # very large
           len*log(alpha) + sum(y^(1/kk)) + sum(log(y) * (1 - 1/kk))
        } else {                         # gumbel
           len*log(alpha) + sum(y) + sum(exp(-y))        
        }
}

