"gumbel.lik" <-
function(pars,data) {
        # computes neg log lik of gumbel model
        gev.lik(pars=c(pars[1],pars[2],0),data=data)
}

