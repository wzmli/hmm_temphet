#state seq log-likelihood

ssll <- function(mod){
  obsnum <- mod@ntimes
  vs <- viterbi(mod)$state
  transmat <- mod@trDens
  obs <- mod@response[[1]][[1]]@y
  transdf <- matrix(1,nrow=obsnum)
  transprob <- 0
  for(i in 1:(obsnum-1)){
    tempprob <- log(transmat[1,vs[i],vs[i+1]])
    transdf[i+1,] <- tempprob
    transprob <- transprob + tempprob 
    
  }
  # Removing NA states wrt obs 
  for(i in 1:(obsnum-1)){
    if(is.na(sum(obs[i],obs[i+1]))){transdf[i+1]<-0}
  }
  
  respprob <- 0
  for(i in 1:(obsnum)){
    tempe <- log(dnorm(x=obs[i],
                       mod@response[[vs[i]]][[1]]@parameters$coefficients,
                       mod@response[[vs[i]]][[1]]@parameters$sd))
    if(is.na(tempe)){tempe=0}
    respprob <- respprob + tempe
  }
  
  return(respprob + sum(transdf))
  
}

ICL <- function(mod){
  return(-2*ssll(mod) + freepars(mod)*log(nobs(mod)))
}
