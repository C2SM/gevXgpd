"expon.lik" <-
function(pars,data,chi=0,tim=length(data)) {
        # computes neg log lik of exponential model including the pdf 
        # of peak frequency
        gpd.lik(pars=c(pars[1],pars[2],0),data=data,chi=chi,tim=tim)
}

