"plot.conf.pot" <- 
function(pot.conf.obj,probs=c(0.025,0.975),nspline=3,...) {
# ==================================================================
# Adds lines of confidence bound on a Pareto diagram 
# Options are passed to lines.
# pot.conf.obj :   confidence bounds object derived from 
#                  call to conf.bounds.pot
    pot.conf <- pot.conf.obj$conf.retv
    lamda <- pot.conf.obj$lamda
    tt <- as.numeric(row.names(pot.conf))
    avail.probs <- as.numeric(dimnames(pot.conf)[[2]])
    probs.ind <- which(avail.probs %in% probs)
    if (length(probs.ind)<=0) {
      print("WARNING plot.conf.pot: desired probs not available in conf.obj!")
    }
    for (k in probs.ind) {

        # spline fit in log transformed coordinates
        yy <- pot.conf[,k]
        xx <- log(tt)
        hh <- spline(x=xx,y=yy,n=nspline*length(yy))

        # now plot that thing
	lines(x=exp(hh$x),y=hh$y,...)
    }
}
