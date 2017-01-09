library(depmixS4)
library(plyr)

t <- rep(0:23,500)
y <- rep(1,12000)

tempdat <- data.frame(y=y,t=t)

system.time(mod <- depmix(y~1
                          , data=tempdat
                          , transition=~1 #cos((2*pi*t)/24)+ sin((2*pi*t)/24)
                          , nstate=2
                          , family=gaussian())
)
getpars(mod)

seed = unlist(strsplit(rtargetname,"[.]"))[2]
set.seed(seed)

#randpars <- sample(-3:3,length(getpars(mod))-6,replace=TRUE)
transprobs <- runif(2,0.001,0.999)
randpars <- c(1-transprobs[1],transprobs[1],1-transprobs[2],transprobs[2])

newmod <- setpars(mod,c(0.5,0.5,randpars,0,1,2,1))
newmod
sim <- simhmm(newmod)
error <- rnorm(12000,0,5) ## real scale
oldsl <- sim@response[[1]][[1]]@y
newsl <- log10(abs(10^oldsl + error)) 
hist(oldsl)
hist(newsl)
df <- data.frame(obs= newsl,states=sim@states,time=t)

system.time(hmm2 <- depmix(obs~1
                           , data=df
                           , transition=~1
                           , nstate=2
                           , family=gaussian())
)

system.time(hmm2s <- fit(hmm2,verbose=FALSE))

system.time(hmm3 <- depmix(obs~1
                           , data=df
                           , transition=~1
                           , nstate=3
                           , family=gaussian())
)

system.time(hmm3s <- fit(hmm3,verbose=FALSE))

system.time(hmm4 <- depmix(obs~1
                           , data=df
                           , transition=~1
                           , nstate=4
                           , family=gaussian())
)

system.time(hmm4s <- fit(hmm4,verbose=FALSE))



# system.time(hmm5 <- depmix(obs~1
#                            , data=df
#                            , transition=~1
#                            , nstate=5
#                            , family=gaussian())
# )
# 
# system.time(hmm5s <- fit(hmm5,verbose=FALSE))

system.time(hmmsin2 <- depmix(obs~1
                              , data=df
                              , transition=~cos((2*pi*t)/24)+ sin((2*pi*t)/24)
                              , nstate=2
                              , family=gaussian())
)

system.time(hmmsin2s <- fit(hmmsin2,verbose=FALSE))

system.time(hmmsin3 <- depmix(obs~1
                              , data=df
                              , transition=~cos((2*pi*t)/24)+ sin((2*pi*t)/24)
                              , nstate=3
                              , family=gaussian())
)

system.time(hmmsin3s <- fit(hmmsin3,verbose=FALSE))

system.time(hmmsin4 <- depmix(obs~1
                              , data=df
                              , transition=~cos((2*pi*t)/24)+ sin((2*pi*t)/24)
                              , nstate=4
                              , family=gaussian())
)

system.time(hmmsin4s <- fit(hmmsin4,verbose=FALSE))



# system.time(hmmsin5 <- depmix(obs~1
#                               , data=df
#                               , transition=~1
#                               , nstate=5
#                               , family=gaussian())
# )
# 
# system.time(hmmsin5s <- fit(hmmsin5,verbose=FALSE))

sumdf <- function(lst){
  LL <- ldply(lst,logLik)
  AIC <- ldply(lst,AIC)
  BIC <- ldply(lst,BIC)
  nstates <-ldply(lst,nstates)
  para <- ldply(lst,freepars)
  model <- c('HMM','HMM','HMM','HMMsin','HMMsin','HMMsin')
  temp <- data.frame(LL=LL$V1
                     ,AIC=AIC$V1
                     ,BIC=BIC$V1
                     ,nstates=nstates$V1
                     ,parameters=para$V1
                     ,model
  )
  return(temp)
}

fitlist <- list(hmm2s,hmm3s,hmm4s,hmmsin2s,hmmsin3s,hmmsin4s)

dat <- sumdf(fitlist)
saveRDS(dat, file=paste("sim",seed,"RDS",sep="."))

# rdsave(seed)