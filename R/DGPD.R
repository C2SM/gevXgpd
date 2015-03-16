"DGPD" <- 
function (xx,alpha,chi=0,k) {
# ==================================================================
# Density of Generalized Pareto Distribution
# see e.g. Palutikov et al. 1999 
# chi  : threshold
# alpha: scale parameter
# k    : shape parameter

# function is listable in argument x
# kthresh defines a threshold for k. if abs(k)<kthresh then the
# Exponential Distribution is taken.

      kthresh <- 0.000001

      ar <- (xx-chi)/alpha

      # Exponential Distribution
      if (abs(k) <= kthresh) {1/alpha*exp(-ar)}

      # Generalized Pareto Distribution
      else if (k < -kthresh) {
         fff <- function(a) {
                  if (a >= 0) (1/alpha * (1-k*a)^(1/k-1)) else 0.0}
         apply(matrix(ar),1,FUN=fff)
      }
      else if (k > kthresh) {
         fff <- function(a) {
                  if ((a >= 0) & (a < 1/k)) {(1/alpha * (1-k*a)^(1/k-1))} 
		  else if (a < 0) {0.0}
		  else if (a >= 1/k) {0.0}}
         apply(matrix(ar),1,FUN=fff)
      }
}
