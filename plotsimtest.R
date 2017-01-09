library(ggplot2)
library(reshape2)
library(dplyr)

dat <- readRDS("simtest.RDS")

findS <- function(Smin,S2,S3,S4,S5){
  if(Smin==S2){return(2)}
  if(Smin==S3){return(3)}
  if(Smin==S4){return(4)}
  if(Smin==S5){return(5)}
}

datABICS <- (dat 
  %>% rowwise()
  %>% mutate(minBIC = min(twoS_BIC,threeS_BIC,fourS_BIC,fiveS_BIC)
             , minAIC = min(twoS_AIC,threeS_AIC,fourS_AIC,fiveS_AIC))
  %>% mutate(minBICStates = findS(minBIC,twoS_BIC,threeS_BIC,fourS_BIC,fiveS_BIC)
             , minAICStates = findS(minAIC,twoS_AIC,threeS_AIC,fourS_AIC,fiveS_AIC))
  %>% select(c(seed,minAICStates,minBICStates))
)

print(table(datABICS$minAICStates))
print(table(datABICS$minBICStates))

AICdat <- datABICS %>% count(minAICStates) %>% transmute(States=minAICStates,n=n,type="AIC")
BICdat <- datABICS %>% count(minBICStates) %>% transmute(States=minBICStates,n=n,type="BIC")

ICdat <- rbind(AICdat,BICdat)

g1 <- (ggplot(ICdat,aes(x=States,y=n,color=type,group=type))
       + geom_point()
       + geom_line()
       + theme_bw()
       + ylab("Counts")
       + xlab("Number of States")
       + ggtitle("Simulation Results: Fitting HMMs to Two-State temphet HMM")
) 

#### fit and sim cat1sinhmm5 
# system("make catsdat.Rout")
# system("make cat1.df.Rout")
# load(".cat1.df.RData")
# source("fitfunctions.R")
# source("mikesim.R")
# source("simfunctions.R")
# fitsin5s <- fitsin(5,cat,seed=2830)
# simsin5 <- simsin(5,cat,fitsin5s)
# 
# dwelldf <- data.frame(hmmsin5obs=simsin5$obs,hmmsin5states=simsin5$states
#                       ,time = cat$Time, count=1)

load("cat1sinhmm5.RData")

for(i in 1:10283){ # hacks, not bothering with the last "group"
  j <- 1
  while(j>0){
    if(dwelldf[i,2] == dwelldf[i+j,2]){
      dwelldf[i,4] <- dwelldf[i,4]+1
      j = j+1}
    if(dwelldf[i,2] != dwelldf[i+j,2]){
      j = -1}
    }
}

dwelldf2 <- (dwelldf
  %>% transmute(States = hmmsin5states
                , time=time
                , count=count)
  %>% group_by(States,time)
  %>% summarise(averageDT=mean(count))
)

g2 <- (ggplot(dwelldf2,aes(x=time,y=averageDT,group=States,color=factor(States)))
  + geom_point()
  + geom_line()
  + theme_bw()
)

summary(fitsin5s)
