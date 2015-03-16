"RGPD" <- 
function (alpha,chi=0,k,n) {
# ==================================================================
# random distributed number with Generalised Pareto 
# Distribution see e.g. Palutikov 1999 
# chi  : location parameter
# alpha: scale parameter
# k    : shape parameter
# n    : number of samples

      # generate n random numbers in [0,1] (Uniform Distribution)
      rr <- runif(n)

      # deduce numbers
      XGPD(rr,alpha,chi,k)

}
