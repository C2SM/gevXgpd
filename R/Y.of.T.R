"Y.of.T" <-
function (t,lamda=1) {
    "y.of.t" <- function(t,lamda) {
       if (t < lamda) {return(NA)}
       return(-log(-log(1-lamda/t)))
    }
    sapply(t,FUN=y.of.t,lamda=lamda)
}
