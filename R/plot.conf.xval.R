"plot.conf.xval" <-
function(xval.conf.obj, probs=c(0.025,0.975), nspline=3, ...) {
# ==================================================================
# Adds lines of confidence bound on a Gumbel diagram 
# Options are passed to lines.
# xval.conf :  confidence bounds object derived from 
#              call to conf.bounds.xval
    xval.conf <- xval.conf.obj$conf.retv
    lamda <- xval.conf.obj$lamda
    plot.lamda <- get(".Last.lamda", envir = globalenv())
    xx <- Y.of.T(as.numeric(row.names(xval.conf)),lamda=plot.lamda)
    avail.probs <- as.numeric(dimnames(xval.conf)[[2]])
    probs.ind <- which(avail.probs %in% probs)
    if (length(probs.ind)<=0) {
      print("WARNING plot.conf.xval: desired probs not available in conf.obj!")
    }
    for (k in probs.ind) {
	lines(spline(x=xx,y=xval.conf[,k],n=nspline*length(xx)),...)
    }
}

