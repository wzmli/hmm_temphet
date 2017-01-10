library(depmixS4)
library(dplyr)
library(plyr)

catid <- unlist(strsplit(rtargetname,split="[.]"))[2]

fitmixsin <- function(state,cat,seed=NULL){
  if(!is.null(seed))set.seed(seed)
  model <- mix(LogDist~1,
               data=cat,
               prior=~cos(2*pi*Time/24)+ sin(2*pi*Time/24),
               nstate=state,
               family=gaussian(),
               initdata=cat)
  set.seed(seed)
  fitmodel <- fit(model,emcontrol=em.control(maxit=500))
  return(fitmodel)
}

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

sumdf <- function(lst){
  BIC <- ldply(lst,BIC)
  nstates <-ldply(lst,nstates)
  para <- ldply(lst,freepars)
  model <- c('FMM + THsin','FMM + THsin','FMM + THsin','FMM + THsin','FMM + THsin')
  type <- c('FMM + TH','FMM + TH','FMM + TH','FMM + TH','FMM + TH')
  temp <- data.frame(BICS=BIC$V1,nstates=nstates$V1,parameters=para$V1,model,type)
  return(temp)
}

fitfmmsin3s <- fitmixsin(3,dat,seed=3)
fitfmmsin4s <- fitmixsin(4,dat,seed=1)
fitfmmsin5s <- fitmixsin(5,dat,seed=1)
fitfmmsin6s <- fitmixsin(6,dat,seed=1)
fitfmmsin7s <- fitmixsin(7,dat,seed=1)

fitlist <- list(fitfmmsin3s,fitfmmsin4s,fitfmmsin5s,fitfmmsin6s,fitfmmsin7s)

catsum <- sumdf(fitlist)
simdf <- data.frame(fmmsin3S=simmixsin(3,dat,fitfmmsin3s)[,2]
                    , fmmsin3obs=simmixsin(3,dat,fitfmmsin3s)[,1]
                    , fmmsin4S=simmixsin(4,dat,fitfmmsin4s)[,2]
                    , fmmsin4obs=simmixsin(4,dat,fitfmmsin4s)[,1]
                    , fmmsin5S=simmixsin(5,dat,fitfmmsin5s)[,2]
                    , fmmsin5obs=simmixsin(5,dat,fitfmmsin5s)[,1]
                    , fmmsin6S=simmixsin(6,dat,fitfmmsin6s)[,2]
                    , fmmsin6obs=simmixsin(6,dat,fitfmmsin6s)[,1]
                    , fmmsin7S=simmixsin(7,dat,fitfmmsin7s)[,2]
                    , fmmsin7obs=simmixsin(7,dat,fitfmmsin7s)[,1]
)
parslist <- list()
for(i in 1:length(fitlist)){
  parslist[[i]] <- getpars(fitlist[[i]])
}

datlist <- list(catsum,parslist,simdf)

saveRDS(datlist,file=paste("cat",catid,"fmmsin","RDS",sep="."))