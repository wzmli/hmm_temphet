library(depmixS4)
library(dplyr)
library(plyr)

files <- commandArgs(trailingOnly = TRUE)
catid <- unlist(strsplit(files[1],split="[.]"))[2]
dat <- readRDS(files[1])
source(files[2])

fitmix <- function(state,cat,seed=2830){
  model <- mix(LogDist~1,
               data=cat,
               prior=~1,
               nstate=state,
               family=gaussian(),
               initdata=cat)
  set.seed(seed)
  fitmodel <- fit(model,emcontrol=em.control(maxit=500))
  return(fitmodel)
}

simmix <- function(state,cat,fit){
  model <- mix(LogDist~1,
               data=cat,
               prior=~1,
               nstate=state,
               family=gaussian(),
               initdata=cat)
  model<-setpars(model,getpars(fit))
  sim <- simulate(model)
  df <- data.frame(obs= sim@response[[1]][[1]]@y,states=sim@states)
  return(df)
}

sumdf <- function(lst){
  BIC <- ldply(lst,BIC)
  nstates <-ldply(lst,nstates)
  para <- ldply(lst,freepars)
  model <- c('FMM','FMM','FMM','FMM','FMM')
  type <- c('FMM','FMM','FMM','FMM','FMM')
  temp <- data.frame(BICS=BIC$V1,nstates=nstates$V1,parameters=para$V1,model,type)
  return(temp)
}

if(catid == 1){
  fitfmm3s <- fitmix(3,dat,seed=43)
  fitfmm4s <- fitmix(4,dat,seed=9)
  fitfmm5s <- fitmix(5,dat,seed=3)
  fitfmm6s <- fitmix(6,dat,seed=1)
  fitfmm7s <- fitmix(7,dat,seed=1)
}

if(catid == 2){
  fitfmm3s <- fitmix(3,dat,seed=562)
  fitfmm4s <- fitmix(4,dat,seed=265)
  fitfmm5s <- fitmix(5,dat,seed=11)
  fitfmm6s <- fitmix(6,dat,seed=48)
  fitfmm7s <- fitmix(7,dat,seed=12)
}

if(catid == 14){
  fitfmm3s <- fitmix(3,dat,seed=64)
  fitfmm4s <- fitmix(4,dat,seed=4)
  fitfmm5s <- fitmix(5,dat,seed=54)
  fitfmm6s <- fitmix(6,dat,seed=16)
  fitfmm7s <- fitmix(7,dat,seed=3)
}

if(catid == 15){
  fitfmm3s <- fitmix(3,dat,seed=3)
  fitfmm4s <- fitmix(4,dat,seed=1)
  fitfmm5s <- fitmix(5,dat,seed=1)
  fitfmm6s <- fitmix(6,dat,seed=1)
  fitfmm7s <- fitmix(7,dat,seed=1)
}

fitlist <- list(fitfmm3s,fitfmm4s,fitfmm5s,fitfmm6s,fitfmm7s)

catsum <- sumdf(fitlist)
simdf <- data.frame(fmm3S=simmix(3,dat,fitfmm3s)[,2]
  , fmm3obs=simmix(3,dat,fitfmm3s)[,1]
  , fmm4S=simmix(4,dat,fitfmm4s)[,2]
  , fmm4obs=simmix(4,dat,fitfmm4s)[,1]
  , fmm5S=simmix(5,dat,fitfmm5s)[,2]
  , fmm5obs=simmix(5,dat,fitfmm5s)[,1]
  , fmm6S=simmix(6,dat,fitfmm6s)[,2]
  , fmm6obs=simmix(6,dat,fitfmm6s)[,1]
  , fmm7S=simmix(7,dat,fitfmm7s)[,2]
  , fmm7obs=simmix(7,dat,fitfmm7s)[,1]
)

parslist <- list()
for(i in 1:length(fitlist)){
  parslist[[i]] <- getpars(fitlist[[i]])
}

datlist <- list(catsum,parslist,simdf)

saveRDS(datlist,file=paste("cat",catid,"fmm","RDS",sep="."))