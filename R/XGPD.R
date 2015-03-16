"XGPD" <- 
function (ff,alpha,chi=0,k) {
# ==================================================================
# Returns extreme value for a given setting of GPD parameters and
# the cumulative frequency. 
# Invert cumulative pdf for Generalised Pareto Distribution.
# Is used for random number generator.
# see e.g. Palutikov et al. 1999.
# 
# ff   : cumulative frequency
# chi  : location parameter
# alpha: scale parameter
# k    : shape parameter

# function is listable in argument ff
# kthresh defines a threshold for k. if abs(k)<kthresh then the
# Exponential Distribution is taken.

      kthresh <- 0.000001

      # FUNCTION DEFINITIONS
      expo <- function (f) {-log(1-f)}
      pare <- function (f,k) {(1-(1.0-f)^k)/k}

      # Calculate extreme value for given cumulative frequency
      if (abs(k) < kthresh) {
         alpha*expo(ff)+chi
      } else {
         alpha*pare(ff,k)+chi
      }
}
