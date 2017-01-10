library(depmixS4)
library(methods)
library(dplyr)
library(plyr)

catid <- unlist(strsplit(rtargetname,split="[.]"))[2]
dat <- readRDS(input_files)

fithmmsin <- function(state,cat,seed=NULL){
  if(!is.null(seed))set.seed(seed)
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~cos(2*pi*Time/24)+ sin(2*pi*Time/24),
                  family=gaussian())
  set.seed(seed)
  fitmodel <- fit(model,emcontrol=em.control(maxit=500))
  return(fitmodel)
}

simhmmsin <- function(state,cat,fit){
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~cos((2*pi*Time)/24)+ sin((2*pi*Time)/24),
                  family=gaussian())
  model<-setpars(model,getpars(fit))
  sim <- simhmm(model)
  df <- data.frame(obs= sim@response[[1]][[1]]@y,states=sim@states)
  return(df)}

sumdf <- function(lst){
  BIC <- ldply(lst,BIC)
  nstates <-ldply(lst,nstates)
  para <- ldply(lst,freepars)
  model <- c('HMM + THsin','HMM + THsin','HMM + THsin','HMM + THsin','HMM + THsin')
  type <- c('HMM + TH','HMM + TH','HMM + TH','HMM + TH','HMM + TH')
  temp <- data.frame(BICS=BIC$V1,nstates=nstates$V1,parameters=para$V1,model,type)
  return(temp)
}

if(catid == 1){
  fithmmsin3s <- fithmmsin(3,dat,seed=2030)
  fithmmsin4s <- fithmmsin(4,dat,seed=2830)
  fithmmsin5s <- fithmmsin(5,dat,seed=2030)
  fithmmsin6s <- fithmmsin(6,dat,seed=2830)
  fithmmsin7s <- fithmmsin(7,dat,seed=2830)
}

if(catid == 2){
  fithmmsin3s <- fithmmsin(3,dat,seed=2030)
  fithmmsin4s <- fithmmsin(4,dat,seed=2830)
  fithmmsin5s <- fithmmsin(5,dat,seed=2030)
  fithmmsin6s <- fithmmsin(6,dat,seed=2830)
  fithmmsin7s <- fithmmsin(7,dat,seed=2830)
}

if(catid == 14){
  fithmmsin3s <- fithmmsin(3,dat,seed=2030)
  fithmmsin4s <- fithmmsin(4,dat,seed=2830)
  fithmmsin5s <- fithmmsin(5,dat,seed=2030)
  fithmmsin6s <- fithmmsin(6,dat,seed=2830)
  fithmmsin7s <- fithmmsin(7,dat,seed=2830)
}


if(catid == 15){
  fithmmsin3s <- fithmmsin(3,dat,seed=3030)
  fithmmsin4s <- fithmmsin(4,dat,seed=2830)
  fithmmsin5s <- fithmmsin(5,dat,seed=3030)
  fithmmsin6s <- fithmmsin(6,dat,seed=2830)
  fithmmsin7s <- fithmmsin(7,dat,seed=2830)
}

fitlist <- list(fithmmsin3s,fithmmsin4s,fithmmsin5s,fithmmsin6s,fithmmsin7s)

catsum <- sumdf(fitlist)
simdf <- data.frame(hmmsin3S=simhmmsin(3,dat,fithmmsin3s)[,2]
                    , hmmsin3obs=simhmmsin(3,dat,fithmmsin3s)[,1]
                    , hmmsin4S=simhmmsin(4,dat,fithmmsin4s)[,2]
                    , hmmsin4obs=simhmmsin(4,dat,fithmmsin4s)[,1]
                    , hmmsin5S=simhmmsin(5,dat,fithmmsin5s)[,2]
                    , hmmsin5obs=simhmmsin(5,dat,fithmmsin5s)[,1]
                    , hmmsin6S=simhmmsin(6,dat,fithmmsin6s)[,2]
                    , hmmsin6obs=simhmmsin(6,dat,fithmmsin6s)[,1]
                    , hmmsin7S=simhmmsin(7,dat,fithmmsin7s)[,2]
                    , hmmsin7obs=simhmmsin(7,dat,fithmmsin7s)[,1]
)
parslist <- list()
for(i in 1:length(fitlist)){
  parslist[[i]] <- getpars(fitlist[[i]])
}

datlist <- list(catsum,parslist,simdf)

saveRDS(datlist,file=paste("cat",catid,"hmmsin","RDS",sep="."))