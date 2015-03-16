"fitGPD.mlik" <-
function(data.df,chi=min(data.df[,1],na.rm=TRUE),
         lamda=1,ret=c(5,10,50,100,500),
         method="Nelder-Mead", maxit=10000, ...) {
# ==================================================================
#  Estimates parameters of GPD-fit to sample data.
#  The estimation is based on maximum likelihood.
#  Also estimates values with a given return period.
#  

   # PREPARATIONS
   nam <- names(data.df)[1]
   len <- dim(data.df)[1]
   ddd <- as.array(data.df[,nam])
   tim <- len/lamda
   
   # CHECK IF LARGER THAN THRESHOLD
   if (any(ddd<chi)) {
       stop("** ERROR ** data must be larger than threshold chi")}

   # FUNCTION DEFINITIONS
   foft <- function (t) {1-1/(lamda*t)}

   # INITIAL VALUES OF OPTIMIZATION 
   # (exact alpha for EXPON,
   #  small k such that distribution has no upper bound
   #  this is different from Coles)
   in.lamda <- lamda
   in.alpha <- mean(ddd-chi)
   in.kk <- -0.0
   init <- c(in.lamda,in.alpha,in.kk)

   # MINIMIZE Negative LogLik
   x <- optim(init, gpd.lik, hessian = TRUE, method = method,
              control = list(maxit = maxit, ...), 
              data=ddd, chi=chi, tim=tim) 
   if (x$convergence > 0) print("WARNING: iterative estimation did not converge")
    
   # Extract estimates or flag if convergence failed
   if (x$convergence > 0) {
     k <- NA
     alpha <- NA
   } else {
     lamda <- x$par[1]   # estimate can slightly deviate from preset value 
                       # (numerical inaccuracies from optimisation
     alpha <- x$par[2]
     k     <- x$par[3]
   }

   # Conduct Kolmogorov-Smirnov Goodness-of-Fit Test
   # Note: This is a misuse of the KS-Test because the distribution paremeters
   # were determined from the data itself. The result may nonetheless pinpoint to
   #Êserious misfits.
   if (x$convergence > 0) {
     ks.pval <- NA
   } else {
     ks.pval <- suppressWarnings(
             ks.test(x=data.df[,nam]-chi,y="FGPD",alpha=alpha,k=k)$p.value)
   } 

   # DETERMINE EXTREME VALUES FOR GIVEN RETURN PERIODS
   if (x$convergence > 0) {
     fitted <- rep(NA,times=length(ret))
   } else {
     fitted <- XGPD(foft(ret),alpha,chi,k)  
   }

   # CONVARIENCE MATRIX OF PARAMETER ESTIMATES
   if (x$convergence > 0) {
     covmat <- matrix(NA,nrow=3,ncol=3)
   } else {
     covmat <- solve(x$hessian)
   } 

   # CONSTRUCT pot OBJECT FOR OUTPUT
   list(
       data=data.df,
       lamda=lamda,
       alpha=alpha,
       chi=chi,
       k=k,
       fit=cbind(ret,fitted),
       cov=covmat,
       dist="GPD",
       ks.pval=ks.pval,
       estim="mlik",
       method=method,
       maxit=maxit
   )
}

