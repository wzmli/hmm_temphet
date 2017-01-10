library(depmixS4)
library(methods)
library(dplyr)
library(plyr)

catid <- unlist(strsplit(rtargetname,split="[.]"))[2]
dat <- readRDS(input_files)

fithmm <- function(state,cat,seed=NULL){
  if(!is.null(seed))set.seed(seed)
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~1,
                  family=gaussian())
  fitmodel <- fit(model,emcontrol=em.control(maxit=500))
  return(fitmodel)
}

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

sumdf <- function(lst){
  BIC <- ldply(lst,BIC)
  nstates <-ldply(lst,nstates)
  para <- ldply(lst,freepars)
  model <- c('HMM','HMM','HMM','HMM','HMM')
  type <- c('HMM','HMM','HMM','HMM','HMM')
  temp <- data.frame(BICS=BIC$V1,nstates=nstates$V1,parameters=para$V1,model,type)
  return(temp)
}
if(catid == 1){
  fithmm3s <- fithmm(3,dat,3)
  fithmm4s <- fithmm(4,dat,2)
  fithmm5s <- fithmm(5,dat,1)
  fithmm6s <- fithmm(6,dat,1)
  fithmm7s <- fithmm(7,dat,1)
}

if(catid == 2){
  fithmm3s <- fithmm(3,dat,2830)
  fithmm4s <- fithmm(4,dat,2)
  fithmm5s <- fithmm(5,dat,1)
  fithmm6s <- fithmm(6,dat,1)
  fithmm7s <- fithmm(7,dat,1)
}

if(catid == 14){
  fithmm3s <- fithmm(3,dat,3030)
  fithmm4s <- fithmm(4,dat,2)
  fithmm5s <- fithmm(5,dat,1)
  fithmm6s <- fithmm(6,dat,1)
  fithmm7s <- fithmm(7,dat,7)
}

if(catid == 15){
  fithmm3s <- fithmm(3,dat,3030)
  fithmm4s <- fithmm(4,dat,2)
  fithmm5s <- fithmm(5,dat,1)
  fithmm6s <- fithmm(6,dat,1)
  fithmm7s <- fithmm(7,dat,1)
}
fitlist <- list(fithmm3s,fithmm4s,fithmm5s,fithmm6s,fithmm7s)

catsum <- sumdf(fitlist)
simdf <- data.frame(hmm3S=simhomo(3,dat,fithmm3s)[,2]
                    , hmm3obs=simhomo(3,dat,fithmm3s)[,1]
                    , hmm4S=simhomo(4,dat,fithmm4s)[,2]
                    , hmm4obs=simhomo(4,dat,fithmm4s)[,1]
                    , hmm5S=simhomo(5,dat,fithmm5s)[,2]
                    , hmm5obs=simhomo(5,dat,fithmm5s)[,1]
                    , hmm6S=simhomo(6,dat,fithmm6s)[,2]
                    , hmm6obs=simhomo(6,dat,fithmm6s)[,1]
                    , hmm7S=simhomo(7,dat,fithmm7s)[,2]
                    , hmm7obs=simhomo(7,dat,fithmm7s)[,1]
)
parslist <- list()
for(i in 1:length(fitlist)){
  parslist[[i]] <- getpars(fitlist[[i]])
}

datlist <- list(catsum,parslist,simdf)

saveRDS(datlist,file=paste("cat",catid,"hmm","RDS",sep="."))