avgdf <- function(df){
  simdf <- df %>% dplyr:::select(-obs) %>% group_by(time) %>% summarise_each(funs(mean))
  obsdf <- df %>% dplyr:::select(c(obs,time)) %>% filter(!is.na(obs)) %>% group_by(time) %>% summarise_each(funs(mean,sd))
  avgdf <- cbind(simdf,obsdf)
  avgdf2 <- melt(avgdf,'time')
  return(avgdf2)
}

vmplots <- function(dat,fit){
  tempmodel <- depmix(LogDist~1,
                      data=dat,
                      nstate=nstates(fit),
                      transition=fit@transition[[1]]@formula,
                      family=gaussian())
  parindex <- length(getpars(tempmodel))
  parind <- parindex - 2*nstates(fit) 
  numstates <- nstates(fit)
  xseq <- seq(-pi,pi,0.01)
  angleden <- data.frame(xseq=xseq)
  angledat <- (angleden 
    %>% mutate(state1 = dvonmises(xseq,getpars(fit)[parind+4*(1-1)+3],getpars(fit)[parind+4*(1-1)+4])
               , state2 = dvonmises(xseq,getpars(fit)[parind+4*(2-1)+3],getpars(fit)[parind+4*(2-1)+4])
               , state3 = dvonmises(xseq,getpars(fit)[parind+4*(3-1)+3],getpars(fit)[parind+4*(3-1)+4])
               , state4 = dvonmises(xseq,getpars(fit)[parind+4*(4-1)+3],getpars(fit)[parind+4*(4-1)+4])
    )
  )
  if(numstates==5){
    angledat <- angledat %>% mutate(state5 = dvonmises(xseq,getpars(fit)[parind+4*(5-1)+3],getpars(fit)[parind+4*(5-1)+4])
    )
  }
  
  mangledat <- reshape2::melt(angledat,"xseq")
  return(return(mangledat))
}

lnplots <- function(dat,fit){
  tempmodel <- depmix(LogDist~1,
                      data=dat,
                      nstate=nstates(fit),
                      transition=fit@transition[[1]]@formula,
                      family=gaussian())
  parindex <- length(getpars(tempmodel))
  parind <- parindex - 2*nstates(fit) 
  numstates <- nstates(fit)
  xseq <- seq(-1,5,0.01)
  lndat <- data.frame(xseq=xseq)
  lndat <- (lndat 
               %>% mutate(state1 = dnorm(xseq,getpars(fit)[parind+4*(1-1)+1],getpars(fit)[parind+4*(1-1)+2])
                          , state2 = dnorm(xseq,getpars(fit)[parind+4*(2-1)+1],getpars(fit)[parind+4*(2-1)+2])
                          , state3 = dnorm(xseq,getpars(fit)[parind+4*(3-1)+1],getpars(fit)[parind+4*(3-1)+2])
                          , state4 = dnorm(xseq,getpars(fit)[parind+4*(4-1)+1],getpars(fit)[parind+4*(4-1)+2])
               )
  )
  if(numstates==5){
    lndat <- lndat %>% mutate(state5 = dnorm(xseq,getpars(fit)[parind+4*(5-1)+1],getpars(fit)[parind+4*(5-1)+2])
    )
  }
  
  mlndat <- reshape2::melt(lndat,"xseq")
  return(return(mlndat))
}
mdis <- lnplots(dat,fit)
mang <- vmplots(dat,fit)

g1 <- ggplot(mdis,aes(x=xseq,y=value,color=variable))+geom_line()
g2 <- ggplot(mang,aes(x=xseq,y=value,color=variable))+geom_line()


  