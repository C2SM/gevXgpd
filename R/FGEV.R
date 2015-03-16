"FGEV" <-
function (xx,alpha,chi,k) {
# ==================================================================
# cumulative pdf for Generalised Extreme Value Distribution
# see e.g. Zwiers and Kharin 1998 
# chi  : location parameter
# alpha: scale parameter
# k    : shape parameter
# In the notation of Coles (2001) the present (chi, alpha, k) 
# correspond to (mu,sigma,-chi)

# function is listable in argument x
# kthresh defines a threshold for k. if abs(k)<kthresh then the
# Gumbel Distribution is taken.

      kthresh <- 0.000001

      ar <- (xx-chi)/alpha

      # Gumbel Distribution
      if (abs(k) <= kthresh) {exp(-exp(-ar))}

      # Weibull Distribution
      else if (k > kthresh) {
         fff <- function(a) {
                  if (a < (1/k)) exp(-(1-k*a)^(1/k)) else 1.0}
         apply(matrix(ar),1,FUN=fff)
      }

      # Frechet Distribution
      else if (k < (-kthresh)) {
         fff <- function(a) {
                  if (a > (1/k)) exp(-(1-k*a)^(1/k)) else 0.0}
         apply(matrix(ar),1,FUN=fff)
      }

}

