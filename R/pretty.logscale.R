"pretty.logscale" <-
function(x,...) {
# ==================================================================
# determines pretty automatic labels for the return period axis of 
# the GEV plot.
# x is a vector of values. pretty values are chosen for the range
# of x. The remaining options are passed to function pretty (linear 
# scale) and can be used to define the typical number of labels.
    round.logscale <- function(x) {
         fix <- c(1,2,5)
         pp <- floor(log10(x*1.001))  # floor(log10(1000)) yields 2 !?
         dd <- abs(log10(x)-(pp+log10(fix)))
         (fix[which(dd == min(dd))[1]]*10^pp)[1]
    }
    hh <- 10^pretty(log(x[x>0],base=10),...)
    unique(sapply(hh,FUN=round.logscale))
}

