"fitGEV.pmlik" <-
function(data.df,lamda=1,ret=c(5,10,50,100,500)/lamda,
                        method="Nelder-Mead", maxit=10000, ...) {
# ==================================================================
#  Estimates parameters of GEV-fit to sample data.
#  The estimation is based on maximum likelihood. In contrast to
#  classical maxlik here the likelihood function is extended to
#  include a geophysical prior function for the shape parameter.
#  The idea of this function it to prevent the estimation of 
#  geophysically unrealistic values of the shape. See Martins and
#  Stedinger 2000.
#  Also estimates values with a given return period.
#  
#  ARGUMENTS:
#  data.df  : one column dataframe with data
#  OPTIONAL ARGUMENTS:
#  ret      : return periods for which extreme values will be 
#             estimated; can be a list of values > 1.
#  method   : the method to use in optim
#  maxit    : maximum number of iterations to use
#  ...      : additional arguments passed to optim
#             
#  OUTPUT: (a list of type xval with the following elements)
#  data:    original dataframe
#  alpha:   scale parameter of GEV
#  chi  :   location parameter of GEV
#  k    :   shape parameter of GEV
#  fit  :   extreme values for predefined return periods ret
#  cov  :   the covariance matrix of max likelihood estimation
#  dist :   what distribution was fitted (GEV or GUMBEL)
#  estim:   what estimator was used (mlik or lmom)
#  method:  method used in optim
#  maxit:   maximum number of iteration allowed

   # CONSTANT
   kthresh <- 0.000001
   pp <- 6  # parameters for geophysical prior
   qq <- 9  # from Martins and Stedinger p. 740
   
   # PREPARATIONS
   nam <- names(data.df)[1]
   len <- dim(data.df)[1]
   ddd <- as.array(data.df[,nam])
   
   # INITIAL VALUES OF OPTIMIZATION (taken from Coles routine)
   #in.alpha <- sqrt(6 * var(ddd))/pi
   #in.chi <- mean(ddd) - 0.57722 * in.alpha
   #in.kk <- -0.1
   #init <- c(in.chi,in.alpha,in.kk)


   # INITIAL VALUES OF OPTIMIZATION (taken from L-moments method)
   bla.xval <- fitGEV.lmom(data.df,lamda=lamda,ret=2/lamda)
   in.alpha <- bla.xval$alpha
   in.chi <- bla.xval$chi
   in.kk <- bla.xval$k
   in.kk <- min(max(-0.3,in.kk),0.1)  # prevent initial values outside +-0.5
   init <- c(in.chi,in.alpha,in.kk)

   # FUNCTION DEFINITIONS
   foft <- function (t) {1-1/(lamda*t)}

   # Geophysical prior function PPI(kk) from Martins and Stedinger p. 740
   ppi <- function(bb) {
        ppii <- function(b) {
           if (abs(b)<0.5){
              ((0.5+b)^(pp-1)*(0.5-b)^(qq-1))/beta(pp,qq)
           }
           else {0}
        }
        sapply(bb,FUN=ppii)
   }

   # LIKELIHOOD AS FUNCTION OF PARAMETERS: negative log of likelihood
   gev.lik <- function(a) {
        # computes neg log lik of gev model including gumbel limit
        chi   <- a[1]
	alpha <- a[2]
	kk    <- a[3]
	y <- (ddd-chi)/alpha
        prior <- ppi(kk)
	if (abs(kk)>kthresh) {
           # gev with prior
	   y <- (1 - kk * y)
	   if(any(y <= 0) || (alpha <= 0) || (prior == 0)) return(10^6) 
                # if unrealistic parameters, return very large value
	   len*log(alpha) + sum(y^(1/kk)) + sum(log(y) * (1 - 1/kk)) - log(prior)
        } else {    
           # gumbel (prior is also needed to warrant continuity between gumb and gev
	   if((alpha <= 0) || (prior == 0)) return(10^6) 
	   len*log(alpha) + sum(y) + sum(exp(-y)) - log(prior)
        }
   }
    
   # MINIMIZE Negative LogLik
   x <- optim(init, gev.lik, hessian = TRUE, method = method,
                    control = list(maxit = maxit, ...)) 
   if (x$convergence > 0) print("WARNING: iterative estimation did not converge")
    
   # Extract estimates  
   k     <- x$par[3]
   alpha <- x$par[2]
   chi   <- x$par[1]

   # flag if convergence failed
   if (x$convergence > 0) {
     k <- NA
     alpha <- NA
     chi <- NA
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
       estim="pmlik",
       method=method,
       maxit=maxit
   )
}

