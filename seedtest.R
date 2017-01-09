##fit and sim functions
library(depmixS4)
library(dplyr)
###Fit and simulate time-homogeneous HMMs
fithomo <- function(state,cat,seed=2830){
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~1,
                  family=gaussian())
  set.seed(seed)
  fitmodel <- fit(model)
  return(fitmodel)
}

#########################################################################

## Fit and Sim of Hourly-transition HMMs

fithourly <- function(state,cat,seed=2830){
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~factor(Time),
                  family=gaussian())
  set.seed(seed)
  fitmodel <- fit(model)
  return(fitmodel)
}

########################################################################

##fit and simulate time-dependent transition HMMs (sin)

fitsin <- function(state,cat,seed=2830){
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~cos(2*pi*Time/24)+ sin(2*pi*Time/24),
                  family=gaussian())
  set.seed(seed)
  fitmodel <- fit(model)
  return(fitmodel)
}
#######################################################################

##fit and simulate time-dependent transition HMMs (quadratic)

fitquad <- function(state,cat,seed=2830){
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~I(Time/24)+I((Time/24)^2),
                  family=gaussian())
  set.seed(seed)
  fitmodel <- fit(model)
  return(fitmodel)
}

###########################################

##fit and simulate time-dependent transition HMMs (block)

fitblock <- function(state,cat,seed=2830){
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~factor(Block),
                  family=gaussian())
  set.seed(seed)
  fitmodel <- fit(model)
  return(fitmodel)
}

###########################################################

##fit and simulate FMMs

fitmix <- function(state,cat,seed=2830){
  model <- mix(LogDist~1,
               data=cat,
               prior=~1,
               nstate=state,
               family=gaussian(),
               initdata=cat)
  set.seed(seed)
  fitmodel <- fit(model)
  return(fitmodel)
}

##fit and simulate time-dependent FMM (Sin)

fitmixsin <- function(state,cat,seed=2830){
  model <- mix(LogDist~1,
               data=cat,
               prior=~cos(2*pi*Time/24)+ sin(2*pi*Time/24),
               nstate=state,
               family=gaussian(),
               initdata=cat)
  set.seed(seed)
  fitmodel <- fit(model)
  return(fitmodel)
}

hh <- fithomo(3,cat,seed=2830)
ss <- fitsin(3,cat,seed=2830)
qq <- fitquad(3,cat,seed=2830)
bb <- fitblock(3,cat,seed=2830)

