"expon.lik.retv" <- 
function(pars,data,chi=0,tim=length(data),ret) {
        # computes neg log lik of exponential model
        # with a reparameterisation in terms of exceedance frequency,
        # return level (in this order in pars).
        # ret is the return period in units of the average interval
        # between threshold exceedances. 
        gpd.lik.retv(pars=c(pars[1],pars[2],0),data=data,
                     chi=chi,tim=tim,ret=ret)
}

