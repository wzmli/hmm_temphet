## 
library(dplyr)
library(reshape2)

catid <- unlist(strsplit(rtargetname,split="[.]"))[2]

hmm <- readRDS(input_files[1])
hmmblock <- readRDS(input_files[2])
hmmsin <- readRDS(input_files[3])
hmmquad <- readRDS(input_files[4])

bicdf <- rbind(hmm[[1]],hmmblock[[1]],hmmsin[[1]],hmmquad[[1]])

bicdf2 <- bicdf %>% mutate(deltaBIC=BICS-min(BICS))

bicp <- bicplot(bicdf2)
print(bicp)


simhmmdf <- simlnvm(dat,hmm[[2]])
simblockdf <- simlnvm(dat,hmmblock[[2]])
simsindf <- simlnvm(dat,hmmsin[[2]])
simquaddf <- simlnvm(dat,hmmquad[[2]])
vithmmdf <- vitsimlnvm(dat,hmm[[3]])

simdat <- data.frame(obs = dat$LogDist
                     , hmm = simhmmdf$steplength
                     , hmmblock = simblockdf$steplength
                     , hmmsin = simsindf$steplength
                     , hmmquad = simquaddf$steplength
                     , vithmm = vithmmdf$steplength
                     , time = dat$Time)


print(avgplot(simdat))
print(acfplot(simdat))

ll <- list(bicdf2,simdat)
saveRDS(ll,file=paste("cat",catid,"paperdf","RDS",sep="."))
