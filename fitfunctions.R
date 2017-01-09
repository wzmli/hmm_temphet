library(depmixS4)

fitmix <- function(state,cat,seed=2830){
  model <- mix(LogDist~1,
               data=cat,
               prior=~1,
               nstate=state,
               family=gaussian(),
               initdata=cat)
  set.seed(seed)
  fitmodel <- fit(model,emcontrol=em.control(maxit=iter))
  return(fitmodel)
}

fitmixsin <- function(state,cat,seed=NULL){
  if(!is.null(seed))set.seed(seed)
  model <- mix(LogDist~1,
               data=cat,
               prior=~cos(2*pi*Time/24)+ sin(2*pi*Time/24),
               nstate=state,
               family=gaussian(),
               initdata=cat)
  set.seed(seed)
  fitmodel <- fit(model,emcontrol=em.control(maxit=iter))
  return(fitmodel)
}


fithomo <- function(state,cat,seed=NULL){
  if(!is.null(seed))set.seed(seed)
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~1,
                  family=gaussian())
  fitmodel <- fit(model,emcontrol=em.control(maxit=iter))
  return(fitmodel)
}


fithourly <- function(state,cat,seed=NULL){
  if(!is.null(seed))set.seed(seed)
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~factor(Time),
                  family=gaussian())
  set.seed(seed)
  fitmodel <- fit(model,emcontrol=em.control(maxit=iter))
  return(fitmodel)
}


fitquad <- function(state,cat,seed=NULL){
  if(!is.null(seed))set.seed(seed)
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~I(Time/24)+I((Time/24)^2),
                  family=gaussian())
  set.seed(seed)
  fitmodel <- fit(model,emcontrol=em.control(maxit=iter))
  return(fitmodel)
}

fitsin <- function(state,cat,seed=NULL){
  if(!is.null(seed))set.seed(seed)
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~cos(2*pi*Time/24)+ sin(2*pi*Time/24),
                  family=gaussian())
  set.seed(seed)
  fitmodel <- fit(model,emcontrol=em.control(maxit=iter))
  return(fitmodel)
}


fitblock <- function(state,cat,seed=NULL){
  if(!is.null(seed))set.seed(seed)
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~factor(Block),
                  family=gaussian())
  set.seed(seed)
  fitmodel <- fit(model,emcontrol=em.control(maxit=iter))
  return(fitmodel)
}
