'qqgev' <- 
function(x,fit,main="GEV Q-Q Plot",
         xlab="Theoretical Quantile",ylab="Empirical Quantile",...) {
    if (missing(fit)) { 
       stop("Argument \'fit\' is missing without defaul.") }
    # probability plotting positions
    plp <- ppoints(x)
    # empirical quantiles
    q.empir <- sort(x,decreasing=FALSE)
    # theoretical quantiles
    q.theor <- XGEV(ff=plp,alpha=fit$alpha,chi=fit$chi,k=fit$k)
    # plot
    plot(x=q.theor,y=q.empir,type="p",xlab=xlab,ylab=ylab,main=main,...)
}
