library(depmixS4)
library(plyr)
library(dplyr)
library(circular)

t <- rep(0:23,500)
y <- rep(1,12000)

tempdat <- data.frame(y=y,t=t)

system.time(mod <- depmix(y~1
                          , data=tempdat
                          , transition=~cos((2*pi*t)/24)+ sin((2*pi*t)/24)
                          , nstate=2
                          , family=gaussian())
)
getpars(mod)

seed = unlist(strsplit(rtargetname,"[.]"))[2]
set.seed(seed)

randpars <- sample(-3:3,length(getpars(mod))-6,replace=TRUE)
newmod <- setpars(mod,c(0.5,0.5,randpars,0,1,2,1))
newmod
sim <- simhmm(newmod)

simang <- function(size,mu,ka){
  temp <- rvonmises(size,mu,ka)
  class(temp) <- "numeric"
  return(temp)
}

df <- data.frame(obs= sim@response[[1]][[1]]@y,states=sim@states,time=t)
df2 <- df %>% mutate(angle=ifelse(states==1,simang(1,0,0),simang(1,0,10)))
hist(df$obs)


xtrans <- function(oldx,r,a){
  oldx + r*cos(a)
}


ytrans <- function(oldy,r,a){
  oldy + r*sin(a)
}

x <- c(0)
y <- c(0)

for(i in 1:nrow(df2)){
  tempx <- xtrans(x[i],10^df2$obs[i],df2$angle[i])
  tempy <- ytrans(y[i],10^df2$obs[i],df2$angle[i])
  newx <- tempx + rnorm(1,0,sqrt(2)*log10(8))
  newy <- tempy + rnorm(1,0,sqrt(2)*log10(8))
  x <- c(x,newx)
  y <- c(y,newy)
}

stepdf <- c()
for(i in 2:length(x)){
  temp <- sqrt((x[i]-x[i-1])^2 + (y[i]-y[i-1])^2)
  stepdf <- c(stepdf,log10(temp))
}

library(moveHMM)
# 
# datxy <- data.frame(ID=1, E=x,N=y)
# transDistdf <- prepData(datxy,type="UTM",coordNames = c("E","N"))
newdf <- data.frame(obs=stepdf)

system.time(hmm2 <- depmix(obs~1
                           , data=newdf
                           , transition=~1
                           , nstate=2
                           , family=gaussian())
)

system.time(hmm2s <- fit(hmm2,verbose=TRUE))

system.time(hmm3 <- depmix(obs~1
                           , data=newdf
                           , transition=~1
                           , nstate=3
                           , family=gaussian())
)

system.time(hmm3s <- fit(hmm3,verbose=TRUE))

system.time(hmm4 <- depmix(obs~1
                           , data=newdf
                           , transition=~1
                           , nstate=4
                           , family=gaussian())
)

system.time(hmm4s <- fit(hmm4,verbose=TRUE))


# 
# system.time(hmm5 <- depmix(obs~1
#                            , data=newdf
#                            , transition=~1
#                            , nstate=5
#                            , family=gaussian())
# )
# 
# system.time(hmm5s <- fit(hmm5,verbose=TRUE))

system.time(hmmsin2 <- depmix(obs~1
                              , data=newdf
                              , transition=~cos((2*pi*t)/24)+ sin((2*pi*t)/24)
                              , nstate=2
                              , family=gaussian())
)

system.time(hmmsin2s <- fit(hmmsin2,verbose=TRUE))

system.time(hmmsin3 <- depmix(obs~1
                              , data=newdf
                              , transition=~cos((2*pi*t)/24)+ sin((2*pi*t)/24)
                              , nstate=3
                              , family=gaussian())
)

system.time(hmmsin3s <- fit(hmmsin3,verbose=TRUE))

system.time(hmmsin4 <- depmix(obs~1
                              , data=newdf
                              , transition=~cos((2*pi*t)/24)+ sin((2*pi*t)/24)
                              , nstate=4
                              , family=gaussian())
)

system.time(hmmsin4s <- fit(hmmsin4,verbose=TRUE))



# system.time(hmmsin5 <- depmix(obs~1
#                               , data=newdf
#                               , transition=~1
#                               , nstate=5
#                               , family=gaussian())
# )
# 
# system.time(hmmsin5s <- fit(hmmsin5,verbose=TRUE))

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
saveRDS(dat, file=paste("simxy",seed,"RDS",sep="."))
