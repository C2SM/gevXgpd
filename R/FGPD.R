"FGPD" <- 
function (xx,alpha,chi=0,k) {
# ==================================================================
# cumulative pdf for Generalised Pareto Distribution
# see e.g. Palutikov et al. 1999 
# chi  : predefined threshold
# alpha: scale parameter
# k    : shape parameter

# function is listable in argument x
# kthresh defines a threshold for k. if abs(k)<kthresh then the
# Exponential Distribution is taken.

      kthresh <- 0.000001

      ar <- (xx-chi)/alpha

      # Exponential
      if (abs(k) <= kthresh) {
	 fff <- function(a) {
		 if (a > 0) (1-exp(-a)) else 0.0}
         apply(matrix(ar),1,FUN=fff)
     }
      # Gerealized Pareto Distribution
      else if (k < -kthresh) {
         fff <- function(a) {
                  if (a > 0) (1-(1-k*a)^(1/k)) else 0.0}
         apply(matrix(ar),1,FUN=fff)
      }
      else if (k > kthresh) {
         fff <- function(a) {
                  if ((a > 0) & (a < 1/k)) {1-(1-k*a)^(1/k)} 
		  else if (a <= 0) {0.0}
		  else if (a >= 1/k) {1.0}}
         apply(matrix(ar),1,FUN=fff)
      }
}