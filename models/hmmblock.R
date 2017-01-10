library(depmixS4)
library(methods)
library(dplyr)
library(plyr)

files <- commandArgs(trailingOnly = TRUE)
catid <- unlist(strsplit(files[1],split="[.]"))[2]
dat <- readRDS(files[1])
source(files[2])

fithmmblock <- function(state,cat,seed=NULL){
  if(!is.null(seed))set.seed(seed)
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~factor(Block),
                  family=gaussian())
  set.seed(seed)
  fitmodel <- fit(model,emcontrol=em.control(maxit=500))
  return(fitmodel)
}

simhmmblock <- function(state,cat,fit){
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~factor(Block),
                  family=gaussian())
  model<-setpars(model,getpars(fit))
  sim <- simhmm(model)
  df <- data.frame(obs= sim@response[[1]][[1]]@y,states=sim@states)
  return(df)}

sumdf <- function(lst){
  BIC <- ldply(lst,BIC)
  nstates <-ldply(lst,nstates)
  para <- ldply(lst,freepars)
  model <- c('HMM + THblock','HMM + THblock','HMM + THblock','HMM + THblock','HMM + THblock')
  type <- c('HMM + TH','HMM + TH','HMM + TH','HMM + TH','HMM + TH')
  temp <- data.frame(BICS=BIC$V1,nstates=nstates$V1,parameters=para$V1,model,type)
  return(temp)
}

if(catid == 1){
  fithmmblock3s <- fithmmblock(3,dat,seed=2830)
  fithmmblock4s <- fithmmblock(4,dat,seed=2830)
  fithmmblock5s <- fithmmblock(5,dat,seed=2830)
  fithmmblock6s <- fithmmblock(6,dat,seed=2830)
  fithmmblock7s <- fithmmblock(7,dat,seed=3030)
}

if(catid == 2){
  fithmmblock3s <- fithmmblock(3,dat,seed=2830)
  fithmmblock4s <- fithmmblock(4,dat,seed=2830)
  fithmmblock5s <- fithmmblock(5,dat,seed=2830)
  fithmmblock6s <- fithmmblock(6,dat,seed=2830)
  fithmmblock7s <- fithmmblock(7,dat,seed=3030)
}

if(catid == 14){
  fithmmblock3s <- fithmmblock(3,dat,seed=2830)
  fithmmblock4s <- fithmmblock(4,dat,seed=2830)
  fithmmblock5s <- fithmmblock(5,dat,seed=2830)
  fithmmblock6s <- fithmmblock(6,dat,seed=2830)
  fithmmblock7s <- fithmmblock(7,dat,seed=2030)
}

if(catid == 15){
  fithmmblock3s <- fithmmblock(3,dat,seed=2830)
  fithmmblock4s <- fithmmblock(4,dat,seed=2830)
  fithmmblock5s <- fithmmblock(5,dat,seed=2830)
  fithmmblock6s <- fithmmblock(6,dat,seed=2830)
  fithmmblock7s <- fithmmblock(7,dat,seed=3030)
}

fitlist <- list(fithmmblock3s,fithmmblock4s,fithmmblock5s,fithmmblock6s,fithmmblock7s)

catsum <- sumdf(fitlist)
simdf <- data.frame(hmmblock3S=simhmmblock(3,dat,fithmmblock3s)[,2]
                    , hmmblock3obs=simhmmblock(3,dat,fithmmblock3s)[,1]
                    , hmmblock4S=simhmmblock(4,dat,fithmmblock4s)[,2]
                    , hmmblock4obs=simhmmblock(4,dat,fithmmblock4s)[,1]
                    , hmmblock5S=simhmmblock(5,dat,fithmmblock5s)[,2]
                    , hmmblock5obs=simhmmblock(5,dat,fithmmblock5s)[,1]
                    , hmmblock6S=simhmmblock(6,dat,fithmmblock6s)[,2]
                    , hmmblock6obs=simhmmblock(6,dat,fithmmblock6s)[,1]
                    , hmmblock7S=simhmmblock(7,dat,fithmmblock7s)[,2]
                    , hmmblock7obs=simhmmblock(7,dat,fithmmblock7s)[,1]
)
parslist <- list()
for(i in 1:length(fitlist)){
  parslist[[i]] <- getpars(fitlist[[i]])
}

datlist <- list(catsum,parslist,simdf)

saveRDS(datlist,file=paste("cat",catid,"hmmblock","RDS",sep="."))