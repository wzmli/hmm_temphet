library(depmixS4)
library(dplyr)
source("mikesim.R")

#sim hmmhomo ----
simhomo <- function(state,cat,fit){
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~1,
                  family=gaussian())
  model<-setpars(model,getpars(fit))
  sim <- simhmm(model)
  df <- data.frame(obs= sim@response[[1]][[1]]@y,states=sim@states)
  return(df)
}

simhmmhomo3 <- simhomo(3,cat,fithomo3s)

#sim hmmhourly ----

simhourly <- function(state,cat,fit){
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~factor(Time),
                  family=gaussian())
  model<-setpars(model,getpars(fit))
  sim <- simhmm(model)
  df <- data.frame(obs= sim@response[[1]][[1]]@y,states=sim@states)
  return(df)}

#sim hmmsin ----

simsin <- function(state,cat,fit){
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~cos((2*pi*Time)/24)+ sin((2*pi*Time)/24),
                  family=gaussian())
  model<-setpars(model,getpars(fit))
  sim <- simhmm(model)
  df <- data.frame(obs= sim@response[[1]][[1]]@y,states=sim@states)
  return(df)}

#sim hmmquad ----

simquad <- function(state,cat,fit){
model <- depmix(LogDist~1,
                data=cat,
                nstate=state,
                transition=~I(Time/24)+I((Time/24)^2),
                family=gaussian())
model<-setpars(model,getpars(fit))
sim <- simhmm(model)
df <- data.frame(obs= sim@response[[1]][[1]]@y,states=sim@states)
return(df)}

#sim hmmblock----

simblock <- function(state,cat,fit){
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~factor(Block),
                  family=gaussian())
  model<-setpars(model,getpars(fit))
  sim <- simhmm(model)
  df <- data.frame(obs= sim@response[[1]][[1]]@y,states=sim@states)
  return(df)}

#sim fmm ----

simmix <- function(state,cat,fit){
  model <- mix(LogDist~1,
               data=cat,
               prior=~1,
               nstate=state,
               family=gaussian(),
               initdata=cat)
  model<-setpars(model,getpars(fit))
  sim <- simhmm(model)
  df <- data.frame(obs= sim@response[[1]][[1]]@y,states=sim@states)
  return(df)
}

#sim fmmsin ----

simmixsin <- function(state,cat,fit){
  model <- mix(LogDist~1,
               data=cat,
               prior=~cos((2*pi*Time)/24)+ sin((2*pi*Time)/24),
               nstate=state,
               family=gaussian(),
               initdata=cat)
  model<-setpars(model,getpars(fit))
  sim <- simulate(model)
  df <- data.frame(obs= sim@response[[1]][[1]]@y,states=sim@states)
  return(df)
}








