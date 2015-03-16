"fitEXPON.mlik" <-
function(data.df,chi=min(data.df[,1],na.rm=TRUE),
         lamda=1,ret=c(5,10,50,100,500),
         method="Nelder-Mead", maxit=10000, ...) {
# ==================================================================
#  Estimates parameters of GUMBEL-fit to sample data.
#  The estimation is based on maximum likelihood.
#  Also estimates values with a given return period.
#  
   
   # PREPARATIONS
   nam <- names(data.df)[1]
   len <- dim(data.df)[1]
   ddd <- as.array(data.df[,nam])
   
   # CHECK IF LARGER THAN THRESHOLD
   if (any(ddd<chi)) {
       stop("** ERROR ** data must be larger than threshold chi")}

   # FUNCTION DEFINITIONS
   foft <- function (t) {1-1/(lamda*t)}

   # LIKELIHOOD ESTIMATOR FOR EXPONENTIAL IS ANALYTICAL 
   alpha <- mean(ddd-chi)
   var.alpha <- alpha^2/len
   k     <- 0.0
   var.lamda <- lamda^2/len

   # Conduct Kolmogorov-Smirnov Goodness-of-Fit Test
   # Note: This is a misuse of the KS-Test because the distribution paremeters
   # were determined from the data itself. The result may nonetheless pinpoint to
   #Êserious misfits.
   ks.pval <- suppressWarnings(
           ks.test(x=data.df[,nam]-chi,y="FGPD",alpha=alpha,k=k)$p.value)

   # DETERMINE EXTREME VALUES FOR GIVEN RETURN PERIODS
   fitted <- XGPD(foft(ret),alpha,chi,k)  

   # covariance matrix
   covmat <- matrix(c(var.alpha,0,0,var.lamda),nrow=2,byrow=TRUE)
 
   # CONSTRUCT xval OBJECT FOR OUTPUT
   list(
       data=data.df,
       lamda=lamda,
       alpha=alpha,
       chi=chi,
       k=k,
       fit=cbind(ret,fitted),
       cov=covmat,
       dist="EXPON",
       ks.pval=ks.pval,
       estim="mlik",
       method=method,
       maxit=maxit
   )
}

