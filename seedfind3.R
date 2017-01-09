##fit and sim functions
library(depmixS4)
library(dplyr)

# setup -----------------

dat <- read.csv('ArchivePantherData.csv')

block <- function(t){
  if(t %in% 7:16)return("block2")
  if(t %in% 17:20)return("block3")
  return("block1")
}

## creating suitable dataframe for depmixS4
dat <- (dat %>% rowwise() %>% transmute(
  cat = animal_id,
  Sex = Sex,
  Time = as.numeric(Time) - 1, ##as.numeric conversion
  Distance = Steplength.m.,
  LogDist = log10(Distance),
  Block = block(Time)
)
)

dat$LogDist[is.infinite(dat$LogDist)] <- NA

cats <- c(1,2,14,15)


dat <- (dat %>% 
          filter(cat %in% cats)
)

# subsetting ---------------
cat1 <- subset(dat,dat$cat == 1 )
cat2 <- subset(dat,dat$cat == 2 )
cat14 <- subset(dat,dat$cat == 14 )
cat15 <- subset(dat,dat$cat == 15 )

#######################




#Fit time-homogeneous HMMs --------------
fithomo <- function(state,cat,seed=2830){
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~1,
                  family=gaussian())
  set.seed(seed)
  fitmodel <- fit(model,verbose=TRUE,emcontrol=em.control(maxit=500))
  return(fitmodel)
}

## Fit and Sim of Hourly-transition HMMs ------

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


##fit and simulate time-dependent transition HMMs (sin) -----

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

##fit and simulate time-dependent transition HMMs (quadratic) ------

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


##fit and simulate time-dependent transition HMMs (block) ------

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


##fit and simulate FMMs -----

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

##fit and simulate time-dependent FMM (Sin) ------

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



###########################
cc <- 0
seed = 1
while(cc < 1){
  hh <- fitmix(3,cat2,seed=seed)
  if (logLik(hh) > -7200){
    cc = 1
  }
  print(seed)
  seed = seed + 1
}
hh <- fitmix(3,cat14,seed=64)
hh <- fitmix(3,cat15,seed=3)
hh <- fitmix(3,cat1,seed=43)
hh <- fitmix(3,cat2,seed=64)
