library(depmixS4)
library(methods)
library(dplyr)
library(plyr)

catid <- unlist(strsplit(rtargetname,split="[.]"))[2]

fithmmhourly <- function(state,cat,seed=NULL){
  if(!is.null(seed))set.seed(seed)
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~factor(Time),
                  family=gaussian())
  set.seed(seed)
  fitmodel <- fit(model,emcontrol=em.control(maxit=500))
  return(fitmodel)
}

simhmmhourly <- function(state,cat,fit){
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~factor(Time),
                  family=gaussian())
  model<-setpars(model,getpars(fit))
  sim <- simhmm(model)
  df <- data.frame(obs= sim@response[[1]][[1]]@y,states=sim@states)
  return(df)}

sumdf <- function(lst){
  BIC <- ldply(lst,BIC)
  nstates <-ldply(lst,nstates)
  para <- ldply(lst,freepars)
  model <- c('HMM + THhourly','HMM + THhourly')
  type <- c('HMM + TH','HMM + TH')
  temp <- data.frame(BICS=BIC$V1,nstates=nstates$V1,parameters=para$V1,model,type)
  return(temp)
}

fithmmhourly3s <- fithmmhourly(3,dat,seed=1)
fithmmhourly4s <- fithmmhourly(4,dat,seed=1)


fitlist <- list(fithmmhourly3s,fithmmhourly4s)

catsum <- sumdf(fitlist)
simdf <- data.frame(hmmhourly3S=simhmmhourly(3,dat,fithmmhourly3s)[,2]
                    , hmmhourly3obs=simhmmhourly(3,dat,fithmmhourly3s)[,1]
                    , hmmhourly4S=simhmmhourly(4,dat,fithmmhourly4s)[,2]
                    , hmmhourly4obs=simhmmhourly(4,dat,fithmmhourly4s)[,1]
)
parslist <- list()
for(i in 1:length(fitlist)){
  parslist[[i]] <- getpars(fitlist[[i]])
}

datlist <- list(catsum,parslist,simdf)

saveRDS(datlist,file=paste("cat",catid,"hmmhourly","RDS",sep="."))