"XGEV" <-
function (ff,alpha,chi,k) {
# ==================================================================
# Returns extreme value for a given setting of GEV parameters and
# the cumulative frequency. 
# Invert cumulative pdf for Generalised Extreme Value Distribution
# see e.g. Zwiers and Kharin 1998.
# 
# ff   : cumulative frequency
# chi  : location parameter
# alpha: scale parameter
# k    : shape parameter

# function is listable in argument ff
# kthresh defines a threshold for k. if abs(k)<kthresh then the
# Gumbel Distribution is taken.

      kthresh <- 0.000001

      # FUNCTION DEFINITIONS
      gumb <- function (f) {-log(-log(f))}
      weib <- function (f,k) {(1.0-exp(k*log(-log(f))))/k}
      gev <- function(f,k) {
               if ((f>1) | (f<0)) {return(NA)}
               if (abs(k) < kthresh) {
                  return(alpha*gumb(f)+chi)
               } else {
                  return(alpha*weib(f,k)+chi)
               }}

      # Calculate return value (quantile) for given cumulative frequency
      sapply(ff,FUN=gev,k=k)
}

