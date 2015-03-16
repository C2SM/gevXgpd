"all.exceeds" <-
function(dat,threshold) {
      ind <- (dat[,1] >= threshold)
      return(dat[ind,])
}

