"plot.thresh.diag" <-
function(diag,diagnostic="mean.exceed",
         print.exceeds=3,cex.exceeds=1,
         xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL,...) {

   # extract relevant data from diag object
   diag.dat <- switch(diagnostic,
       "mean.exceed" = diag$q.mean.exceed,
       "shape" = diag$q.shape,
       "scale.star" = diag$q.scale.star,
       stop("Selected diagnostic not available."))
   threshs <- diag$threshs

   # default x-lim, y-lim
   if(is.null(xlim)) {xlim <- c(min(threshs),max(threshs))}
   if(is.null(ylim)) {
       ylim <- c(min(diag.dat,na.rm=TRUE),max(diag.dat,na.rm=TRUE))}
   if(is.null(xlab)) {xlab <- "Threshold"}
   if(is.null(ylab)) {
        ylab <- switch(diagnostic,
                "mean.exceed" = "Mean Exceedance",
                "shape" = "Shape",
                "scale.star" = "Modified Scale",
                stop("Selected diagnostic not available."))}

   # plot estimates
   yy <- diag.dat[,"0.5"]
   plot(x=threshs,y=yy,xlim=xlim,ylim=ylim,type="b",
        xlab=xlab,ylab=ylab,...)

   # plot confidence intervals
   prbs <- dimnames(diag.dat)[[2]]
   c1 <- prbs[1]
   c2 <- prbs[length(prbs)]
   for (u in threshs) {
      ii <- which(u==threshs)
      lines(x=rep(u,2),y=diag.dat[ii,c(c1,c2)],lwd=1.5)
   }

   # print number of exceedances
   if (print.exceeds > 0) {
      ii <- seq(1,length(threshs),by=print.exceeds)
      text(x=threshs[ii],y=rep(ylim[2],times=length(ii)),
           labels=as.character(diag$num.exceed[ii]),cex=cex.exceeds,
           adj=1,srt=90)
   }
}
