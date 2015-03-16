"plot.xval" <-
function(fit.xval,tlim=NULL,tat,tlab="T",
         ylim=NULL,ylab,col="black",col.fit=col,
         plot.data=TRUE,label.data=FALSE,label.cex=1,label.col=col,
         label.adj=c(-0.2,-0.2),cex.axis,
         main="Block Maxima",add=FALSE,...) {
# ==================================================================
# Plot an extreme value diagram for an extreme value analysis
# Options can be passed to plot (points of sample) and line (gev fit)

# fit.xval : an extreme value object
# tlim : the limits of the return period axis
# tat : axis ticks for return period axis
# tlab : label for return period axis
# ylab : label for y axis
# col : color for points and fit (if fit is not specified)
# fit.col : color of fit

   # preparation of data
   dat.df <- fit.xval$data
   var.nam <- names(dat.df)[1]
   sample.size  <- dim(dat.df)[1]
   lamda <- fit.xval$lamda
   
   # handling of optional arguments (graphics options)
   ylab <- if (missing(ylab)) {var.nam} else {ylab}
   pttype <- ifelse(plot.data,"p","n")
   xlim <- if (!is.null(tlim[1])) {
        Y.of.T(c(max(tlim[1],1.0001/lamda),tlim[2]),lamda=1)} else {NULL}

   # get the adopted lamda if there is a plot alread existing
   if (add) {                # no new plot 
      plot.lamda <- get(".Last.lamda", envir = globalenv())
   } else {
      plot.lamda <- lamda
      assign(".Last.lamda", lamda, envir = globalenv()) # save lamda of plot.
   }

   # calculate plotting points
   dat.df[,"RANK"] <- rank(dat.df[,var.nam])
   dat.df[,"E.FRE"] <- dat.df[,"RANK"]/(sample.size+1)  # Weibull-Estimate
   dat.df[,"GUMBEL.Y"] <- Y.of.T(T.of.F(t(dat.df["E.FRE"]),lamda=lamda),
                                 lamda=plot.lamda)

   # select the points really in the chosen tlim (this ensures that ylim
   # is restricted to those points really plotted) 
   if (!is.null(xlim)) {
      ii <- ((dat.df[,"GUMBEL.Y"] >= xlim[1]) & 
             (dat.df[,"GUMBEL.Y"] <= xlim[2]))
   } else { 
      ii <- rep(TRUE,times=length(dat.df[,"GUMBEL.Y"])) 
   }

   # prepare the plot
   if (!add) {
     plot(x=dat.df[ii,"GUMBEL.Y"], y=dat.df[ii,var.nam], axes=FALSE,
          xlim=xlim, ylim=ylim, xlab=tlab, ylab=ylab, 
          type="n", main=main, ...)
     if (missing(cex.axis)) {cex.axis <- par("cex.axis")}
     axis(2,cex.axis=cex.axis)
     box()

     # plot the return period axis
     uuu <- par(c("usr"))
     uut <- T.of.Y(uuu[c(1,2)],lamda=plot.lamda)
     tat  <- if (missing(tat)) {
              pretty.logscale(uut,n=8)
           } else {tat}
     axis(1,at=Y.of.T(tat,lamda=plot.lamda),
          labels=as.character(round(tat,2)),
          cex.axis=cex.axis)
   }

   # plot points
   points(x=dat.df[ii,"GUMBEL.Y"], y=dat.df[ii,var.nam], col=col, 
          type=pttype, ...)

   # plot fitted GEV distribution
   uuu <- par(c("usr"))
   xxx  <- seq(uuu[1],uuu[2],length=200)
   lines(x=xxx,
         y=XGEV(F.of.T(T.of.Y(xxx,lamda=plot.lamda),lamda=lamda),
                fit.xval$alpha,fit.xval$chi,fit.xval$k),
         col=col.fit,...)

   # label the points if
   if (label.data) {
      text(x=dat.df[,"GUMBEL.Y"], y=dat.df[,var.nam],
           labels = as.character(dat.df[,"lab"]),
           cex = label.cex, col = label.col, adj = label.adj)
   }
     
}

