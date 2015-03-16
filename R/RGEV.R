"RGEV" <-
function (alpha,chi,k,n) {
# ==================================================================
# random distributed number with Generalised Extreme Value 
# Distribution see e.g. Zwiers and Kharin 1998 
# chi  : location parameter
# alpha: scale parameter
# k    : shape parameter
# n    : number of samples

# kthresh defines a threshold for k. if abs(k)<kthresh then the
# Gumbel Distribution is taken.

      kthresh <- 0.000001

      # generate n random numbers in [0,1] (Uniform Distribution)
      rr <- runif(n)

      # deduce numbers
      XGEV(rr,alpha,chi,k)

}

