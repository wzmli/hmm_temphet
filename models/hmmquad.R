library(depmixS4)
library(methods)
library(dplyr)
library(plyr)

files <- commandArgs(trailingOnly = TRUE)
catid <- unlist(strsplit(files[1],split="[.]"))[2]
dat <- readRDS(files[1])
source(files[2])

fithmmquad <- function(state,cat,seed=NULL){
  if(!is.null(seed))set.seed(seed)
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~I(Time/24)+I((Time/24)^2),
                  family=gaussian())
  set.seed(seed)
  fitmodel <- fit(model,emcontrol=em.control(maxit=500))
  return(fitmodel)
}
simhmmquad <- function(state,cat,fit){
  model <- depmix(LogDist~1,
                  data=cat,
                  nstate=state,
                  transition=~I(Time/24)+I((Time/24)^2),
                  family=gaussian())
  model<-setpars(model,getpars(fit))
  sim <- simhmm(model)
  df <- data.frame(obs= sim@response[[1]][[1]]@y,states=sim@states)
  return(df)}

sumdf <- function(lst){
  BIC <- ldply(lst,BIC)
  nstates <-ldply(lst,nstates)
  para <- ldply(lst,freepars)
  model <- c('HMM + THquad','HMM + THquad','HMM + THquad','HMM + THquad','HMM + THquad')
  type <- c('HMM + TH','HMM + TH','HMM + TH','HMM + TH','HMM + TH')
  temp <- data.frame(BICS=BIC$V1,nstates=nstates$V1,parameters=para$V1,model,type)
  return(temp)
}

fithmmquad3s <- fithmmquad(3,dat,seed=3030)
fithmmquad4s <- fithmmquad(4,dat,seed=3030)
fithmmquad5s <- fithmmquad(5,dat,seed=3030)
fithmmquad6s <- fithmmquad(6,dat,seed=3030)
fithmmquad7s <- fithmmquad(7,dat,seed=3030)

fitlist <- list(fithmmquad3s,fithmmquad4s,fithmmquad5s,fithmmquad6s,fithmmquad7s)

catsum <- sumdf(fitlist)
simdf <- data.frame(hmmquad3S=simhmmquad(3,dat,fithmmquad3s)[,2]
                    , hmmquad3obs=simhmmquad(3,dat,fithmmquad3s)[,1]
                    , hmmquad4S=simhmmquad(4,dat,fithmmquad4s)[,2]
                    , hmmquad4obs=simhmmquad(4,dat,fithmmquad4s)[,1]
                    , hmmquad5S=simhmmquad(5,dat,fithmmquad5s)[,2]
                    , hmmquad5obs=simhmmquad(5,dat,fithmmquad5s)[,1]
                    , hmmquad6S=simhmmquad(6,dat,fithmmquad6s)[,2]
                    , hmmquad6obs=simhmmquad(6,dat,fithmmquad6s)[,1]
                    , hmmquad7S=simhmmquad(7,dat,fithmmquad7s)[,2]
                    , hmmquad7obs=simhmmquad(7,dat,fithmmquad7s)[,1]
)
parslist <- list()
for(i in 1:length(fitlist)){
  parslist[[i]] <- getpars(fitlist[[i]])
}

datlist <- list(catsum,parslist,simdf)

saveRDS(datlist,file=paste("cat",catid,"hmmquad","RDS",sep="."))