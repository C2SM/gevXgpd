"fitGEV.mlik" <-
function(data.df,
         lamda=1,ret=c(5,10,50,100,500)/lamda,
         method="Nelder-Mead", maxit=10000, ...) {
# ==================================================================
#  Estimates parameters of GEV-fit to sample data.
#  The estimation is based on maximum likelihood.
#  Also estimates values with a given return period.
#  

   # PREPARATIONS
   nam <- names(data.df)[1]
   len <- dim(data.df)[1]
   ddd <- as.array(data.df[,nam])
   
   # FUNCTION DEFINITIONS
   foft <- function (t) {1-1/(lamda*t)}

   # INITIAL VALUES OF OPTIMIZATION (taken from L-moments method)
   bla.xval <- fitGEV.lmom(data.df,lamda=lamda,ret=2/lamda)
   in.alpha <- bla.xval$alpha
   in.chi <- bla.xval$chi
   in.kk <- bla.xval$k
   init <- c(in.chi,in.alpha,in.kk)

   # MINIMIZE Negative LogLik
   x <- optim(init, gev.lik, hessian = TRUE, method = method,
              control = list(maxit = maxit, ...),
              data=ddd) 
   if (x$convergence > 0) print("WARNING: iterative estimation did not converge")
    
   # Extract estimates or flag if convergence failed
   if (x$convergence > 0) {
     k <- NA
     alpha <- NA
     chi <- NA
   } else {
     k     <- x$par[3]
     alpha <- x$par[2]
     chi   <- x$par[1]
   }

   # Conduct Kolmogorov-Smirnov Goodness-of-Fit Test
   # Note: This is a misuse of the KS-Test because the distribution paremeters
   # were determined from the data itself. The result may nonetheless pinpoint to
   #Êserious misfits.
   if (x$convergence > 0) {
     ks.pval <- NA
   } else {
     ks.pval <- suppressWarnings(
             ks.test(x=data.df[,nam],y="FGEV",alpha=alpha,chi=chi,k=k)$p.value)
   } 

   # DETERMINE EXTREME VALUES FOR GIVEN RETURN PERIODS
   if (x$convergence > 0) {
     fitted <- rep(NA,times=length(ret))
   } else {
     fitted <- XGEV(foft(ret),alpha,chi,k)  
   } 

   # CONVARIENCE MATRIX OF PARAMETER ESTIMATES
   if (x$convergence > 0) {
     covmat <- matrix(NA,nrow=3,ncol=3)
   } else {
     covmat <- solve(x$hessian)
   } 

   # CONSTRUCT xval OBJECT FOR OUTPUT
   list(
       data=data.df,
       lamda=lamda,
       alpha=alpha,
       chi=chi,
       k=k,
       fit=cbind(ret,fitted),
       cov=covmat,
       dist="GEV",
       ks.pval=ks.pval,
       estim="mlik",
       method=method,
       maxit=maxit
   )
}

