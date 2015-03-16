"plot.pot" <-
function(fit.pot,tlim=NULL,tat,tlab="T",
         ylim=NULL,ylab,col="black",col.fit=col,
         plot.data=TRUE,label.data=FALSE,label.cex=1,label.col=col,
         label.adj=c(-0.2,-0.2),cex.axis,
         main="Peaks over Threshold",add=FALSE,...) {
# ==================================================================
# Plot an extreme value diagram for a peak over threshold analysis.
# Pareto Diagram.

# fit.pot : a peak over threshold object
# tlim : the limits of the return period axis
# tat : axis ticks for return period axis
# tlab : label for return period axis
# ylab : label for y axis
# col : color for points and fit (if fit is not specified)
# fit.col : color of fit

   # preparation of data
   dat.df <- fit.pot$data
   var.nam <- names(dat.df)[1]
   sample.size  <- dim(dat.df)[1]
   lamda <- fit.pot$lamda
   
   # handling of optional arguments (graphics options)
   ylab <- if (missing(ylab)) {var.nam} else {ylab}
   pttype <- ifelse(plot.data,"p","n")
    
   # calculate plotting points via rank, empirical frequency 
   # and PARETO V for data points 
   dat.df[,"RANK"] <- rank(dat.df[,var.nam])
   dat.df[,"E.FRE"] <- dat.df[,"RANK"]/(sample.size+1)  # Weibull-Estimate
   dat.df[,"E.T"] <- T.of.F(dat.df[,"E.FRE"],lamda=lamda)
   
   # select the points really in the chosen tlim (this ensures that ylim
   # is restricted to those points really plotted) 
   if (!is.null(tlim)) {
      ii <- ((dat.df[,"E.T"] >= tlim[1]) & 
             (dat.df[,"E.T"] <= tlim[2]))
   } else { 
      ii <- rep(TRUE,times=length(dat.df[,"E.T"])) 
   }

   # prepare the plot
   if (!add) {
     plot(x=dat.df[ii,"E.T"], y=dat.df[ii,var.nam], axes=FALSE,
          xlim=tlim, ylim=ylim, xlab=tlab, ylab=ylab, 
          type="n", log="x", main=main,...)
     if (missing(cex.axis)) {cex.axis <- par("cex.axis")}
     axis(2,cex.axis=cex.axis)
     box()

     # plot the return period axis
     uuu <- par(c("usr"))
     uut <- 10^uuu[c(1,2)]
     tat  <- if (missing(tat)) {
         c(1/lamda,
           pretty.logscale(uut,n=8))
         } else {tat}
     axis(1,at=tat,labels=as.character(round(tat,2)),
          cex.axis=cex.axis)
   }

   # plot data-points
   points(x=dat.df[ii,"E.T"], y=dat.df[ii,var.nam], col=col, 
          type=pttype, ...)
   
   # plot fitted GPD distribution
   uuu <- par(c("usr"))
   xxx  <- seq(uuu[1],uuu[2],length=200)
   ttt <- 10^xxx
   lines(x=ttt,
         y=XGPD(F.of.T(ttt,lamda=lamda),
         fit.pot$alpha,fit.pot$chi,fit.pot$k),
         col=col.fit,...)

   # label the points if
   if (label.data) {
      text(x=dat.df[,"E.T"], y=dat.df[,var.nam],
           labels = as.character(dat.df[,"lab"]),
           cex = label.cex, col = label.col, adj = label.adj)
   }
     
}

