## WEIBULL VM 

library(depmixS4)
library(methods)
library(circular)
library(MASS)
library(stats)
library(dplyr)
library(plyr)

catid <- unlist(strsplit(rtargetname,split="[.]"))[2]


# VM setup ----
setClass("vonMises", contains="response")
## library(VGAM)  ## this may mess up multinomial()


setGeneric("vonMises", function(y, pstart = NULL, fixed = NULL, ...) 
  standardGeneric("vonMises"))

setMethod("vonMises", 
          signature(y = "ANY"), 
          function(y, pstart = NULL, fixed = NULL, ...) {
            y <- matrix(y, length(y))
            x <- matrix(1) 
            parameters <- list()
            npar <- 2
            if(is.null(fixed)) fixed <- as.logical(rep(0, npar))
            if(!is.null(pstart)) {
              if(length(pstart) != npar) stop("length of 'pstart' must be ", npar)
              parameters$mu <- pstart[1]
              parameters$kappa <- pstart[2]
            }
            mod <- new("vonMises", parameters = parameters, fixed = fixed,
                       x = x, y = y, npar = npar)
            mod
          }
)

setMethod("dens","vonMises",
          function(object,log=FALSE) {
            dvonmises(object@y, mu = predict(object), 
                      kappa = object@parameters$kappa, 
                      log = log)
          }
)

setMethod("getpars","response",
          function(object,which="pars",...) {
            switch(which,
                   "pars" = {
                     parameters <- numeric()
                     parameters <- unlist(object@parameters)
                     pars <- parameters
                   },
                   "fixed" = {
                     pars <- object@fixed
                   }
            )
            return(pars)
          }
)

setMethod("setpars","vonMises",
          function(object, values, which="pars", ...) {
            npar <- npar(object)
            if(length(values)!=npar) stop("length of 'values' must be",npar)
            # determine whether parameters or fixed constraints are being set
            nms <- names(object@parameters)
            switch(which,
                   "pars"= {
                     object@parameters$mu <- values[1]
                     object@parameters$kappa <- values[2]
                   },
                   "fixed" = {
                     object@fixed <- as.logical(values)
                   }
            )
            names(object@parameters) <- nms
            return(object)
          }
)

setMethod("predict","vonMises", 
          function(object) {
            ret <- object@parameters$mu
            return(ret)
          }
)


##Everything above ""should"" be fine

##*****************************
dvonmisesrad <- circular:::DvonmisesRad
setMethod("fit", "vonMises",
          function(object, w) {
            y <- object@y
            nas <- is.na(rowSums(object@y))
            start <- with(object@parameters,
                          c(mu=mu,logk=log(kappa)))
            objfun <- function(pars) {
              L <- -dvonmisesrad(as.matrix(object@y[!nas,]),pars[1],exp(pars[2]),log=TRUE)
              sum(w[!nas]*L)/sum(w[!nas])
            }
            opt <- optim(fn=objfun,par=start,method="Nelder-Mead")
            pars <- unname(c(opt$par[1],exp(opt$par[2])))
            object <- setpars(object,pars)
            object
          }
          
)

# weibull setup ----
setClass("weibull", contains="response")
library(MASS)
library(stats)
setGeneric("weibull", function(y, pstart = NULL, fixed = NULL, ...) 
  standardGeneric("weibull"))

setMethod("weibull", 
          signature(y = "ANY"), 
          function(y, pstart = NULL, fixed = NULL, ...) {
            y <- matrix(y, length(y))
            x <- matrix(1) 
            parameters <- list()
            npar <- 2
            if(is.null(fixed)) fixed <- as.logical(rep(0, npar))
            if(!is.null(pstart)) {
              if(length(pstart) != npar) stop("length of 'pstart' must be ", npar)
              parameters$shape <- pstart[1]
              parameters$scale <- pstart[2]
            }
            mod <- new("weibull", parameters = parameters, fixed = fixed,
                       x = x, y = y, npar = npar)
            mod
          }
)

setMethod("dens","weibull",
          function(object, log=FALSE) {
            dweibull(object@y, 
                     shape = object@parameters$shape, 
                     scale = object@parameters$scale, 
                     log = log)
          }
)

setMethod("getpars","response",
          function(object,which="pars",...) {
            switch(which,
                   "pars" = {
                     parameters <- numeric()
                     parameters <- unlist(object@parameters)
                     pars <- parameters
                   },
                   "fixed" = {
                     pars <- object@fixed
                   }
            )
            return(pars)
          }
)

setMethod("setpars","weibull",
          function(object, values, which="pars", ...) {
            npar <- npar(object)
            if(length(values)!=npar) stop("length of 'values' must be",npar)
            # determine whether parameters or fixed constraints are being set
            nms <- names(object@parameters)
            switch(which,
                   "pars"= {
                     object@parameters$shape <- values[1]
                     object@parameters$scale <- values[2]
                   },
                   "fixed" = {
                     object@fixed <- as.logical(values)
                   }
            )
            names(object@parameters) <- nms
            return(object)
          }
)

setMethod("fit", "weibull",
          function(object, w) {
            y <- object@y
            nas <- is.na(rowSums(y))
            start <- with(object@parameters,
                          c(logshape=log(shape),logscale=log(scale)))
            objfun <- function(pars) {
              L <- -dweibull(c(na.omit(y)),
                             exp(pars[1]),exp(pars[2]),log=TRUE)
              sum(w[!nas]*L)/sum(w[!nas])
            }
            opt <- optim(fn=objfun,par=start,method="Nelder-Mead")
            pars <- unname(exp(opt$par))
#            print(pars)
            object <- setpars(object,pars)
            object
          }
          
)

# depmix model setup 3 states ----
dist <- pmax(dat$Distance,1e-3)
rModels <- list()
rModels[[1]] <- list()
rModels[[1]][[1]] <- weibull(dist, pstart = c(0.5,0.5),data=dat)
rModels[[1]][[2]] <- vonMises(dat$Turningangle, pstart = c(0, 1),data=dat)

rModels[[2]] <- list()
rModels[[2]][[1]] <- weibull(dist, pstart = c(1.5,1.5),data=dat)
rModels[[2]][[2]] <- vonMises(dat$Turningangle, pstart = c(1, 1),data=dat)

rModels[[3]] <- list()
rModels[[3]][[1]] <- weibull(dist, pstart = c(1,4),data=dat)
rModels[[3]][[2]] <- vonMises(dat$Turningangle, pstart = c(3, 1))

trstart <- rep(1/3,9)
transition <- list()
transition[[1]] <- transInit(~ 1, nst = 3, data=dat,
                             pstart = rep(1/3, 3))
transition[[2]] <- transInit(~ 1, nst = 3,data=dat,
                             pstart = rep(1/3, 3))
transition[[3]] <- transInit(~ 1, nst = 3,data=dat,
                             pstart = rep(1/3, 3))
inMod <- transInit(~ 1, ns = 3, pstart = rep(1/3, 3),
                   data = data.frame(1))
mod3 <- makeDepmix(response = rModels, transition = transition,
                   prior=inMod,homogeneous = FALSE)

wbvmhmm3s <- fit(mod3, verbose = TRUE, emc=em.control(rand=FALSE, maxit=460))

# HMM WVM 4 states ----
rModels <- list()
rModels[[1]] <- list()
rModels[[1]][[1]] <- weibull(dist, pstart = c(0.5,1),data=dat)
rModels[[1]][[2]] <- vonMises(dat$Turningangle, pstart = c(0, 1),data=dat)

rModels[[2]] <- list()
rModels[[2]][[1]] <- weibull(dist, pstart = c(1,2),data=dat)
rModels[[2]][[2]] <- vonMises(dat$Turningangle, pstart = c(1, 1),data=dat)

rModels[[3]] <- list()
rModels[[3]][[1]] <- weibull(dist, pstart = c(3,3),data=dat)
rModels[[3]][[2]] <- vonMises(dat$Turningangle, pstart = c(2, 1))

rModels[[4]] <- list()
rModels[[4]][[1]] <- weibull(dist, pstart = c(0.6,5),data=dat)
rModels[[4]][[2]] <- vonMises(dat$Turningangle, pstart = c(3, 1))

trstart <- rep(1/4,16)
transition <- list()
transition[[1]] <- transInit(~ 1, nst = 4, data=dat,
                             pstart = rep(1/4, 4))
transition[[2]] <- transInit(~ 1, nst = 4,data=dat,
                             pstart = rep(1/4, 4))
transition[[3]] <- transInit(~ 1, nst = 4,data=dat,
                             pstart = rep(1/4, 4))
transition[[4]] <- transInit(~ 1, nst = 4,data=dat,
                             pstart = rep(1/4, 4))
inMod <- transInit(~ 1, ns = 4, pstart = rep(1/4, 4),
                   data = data.frame(1))
mod4 <- makeDepmix(response = rModels, transition = transition,
                   prior=inMod,homogeneous = FALSE)
wbvmhmm4s <- fit(mod4, verbose = TRUE, emc=em.control(rand=FALSE, maxit=460))

# HMM WVM 5 states ----

rModels <- list()
rModels[[1]] <- list()
rModels[[1]][[1]] <- weibull(dist, pstart = c(0.5,1),data=dat)
rModels[[1]][[2]] <- vonMises(dat$Turningangle, pstart = c(0, 1),data=dat)

rModels[[2]] <- list()
rModels[[2]][[1]] <- weibull(dist, pstart = c(1,1),data=dat)
rModels[[2]][[2]] <- vonMises(dat$Turningangle, pstart = c(1, 1),data=dat)

rModels[[3]] <- list()
rModels[[3]][[1]] <- weibull(dist, pstart = c(1,2),data=dat)
rModels[[3]][[2]] <- vonMises(dat$Turningangle, pstart = c(2, 1))

rModels[[4]] <- list()
rModels[[4]][[1]] <- weibull(dist, pstart = c(2,3),data=dat)
rModels[[4]][[2]] <- vonMises(dat$Turningangle, pstart = c(3, 1))

rModels[[5]] <- list()
rModels[[5]][[1]] <- weibull(dist, pstart = c(2,5),data=dat)
rModels[[5]][[2]] <- vonMises(dat$Turningangle, pstart = c(4, 1))

trstart <- rep(1/5,25)
transition <- list()
transition[[1]] <- transInit(~ 1, nst = 5, data=dat,
                             pstart = rep(1/5, 5))
transition[[2]] <- transInit(~ 1, nst = 5,data=dat,
                             pstart = rep(1/5, 5))
transition[[3]] <- transInit(~ 1, nst = 5,data=dat,
                             pstart = rep(1/5, 5))
transition[[4]] <- transInit(~ 1, nst = 5,data=dat,
                             pstart = rep(1/5, 5))
transition[[5]] <- transInit(~ 1, nst = 5,data=dat,
                             pstart = rep(1/5, 5))
inMod <- transInit(~ 1, ns = 5, pstart = rep(1/5, 5),
                   data = data.frame(1))
mod5 <- makeDepmix(response = rModels, transition = transition,
                   prior=inMod,homogeneous = FALSE)
wbvmhmm5s <- fit(mod5, verbose = TRUE, emc=em.control(rand=TRUE, maxit=460))

# HMM WVM 6 states ----

rModels <- list()
rModels[[1]] <- list()
rModels[[1]][[1]] <- weibull(dist, pstart = c(0.5,1),data=dat)
rModels[[1]][[2]] <- vonMises(dat$Turningangle, pstart = c(0, 1),data=dat)

rModels[[2]] <- list()
rModels[[2]][[1]] <- weibull(dist, pstart = c(1,1),data=dat)
rModels[[2]][[2]] <- vonMises(dat$Turningangle, pstart = c(1, 1),data=dat)

rModels[[3]] <- list()
rModels[[3]][[1]] <- weibull(dist, pstart = c(2,5),data=dat)
rModels[[3]][[2]] <- vonMises(dat$Turningangle, pstart = c(2, 1))

rModels[[4]] <- list()
rModels[[4]][[1]] <- weibull(dist, pstart = c(3,1),data=dat)
rModels[[4]][[2]] <- vonMises(dat$Turningangle, pstart = c(3, 1))

rModels[[5]] <- list()
rModels[[5]][[1]] <- weibull(dist, pstart = c(5,5),data=dat)
rModels[[5]][[2]] <- vonMises(dat$Turningangle, pstart = c(4, 1))

rModels[[6]] <- list()
rModels[[6]][[1]] <- weibull(dist, pstart = c(4,2),data=dat)
rModels[[6]][[2]] <- vonMises(dat$Turningangle, pstart = c(5, 1))

trstart <- rep(1/6,36)
transition <- list()
transition[[1]] <- transInit(~ 1, nst = 6, data=dat,
                             pstart = rep(1/6, 6))
transition[[2]] <- transInit(~ 1, nst = 6,data=dat,
                             pstart = rep(1/6, 6))
transition[[3]] <- transInit(~ 1, nst = 6,data=dat,
                             pstart = rep(1/6, 6))
transition[[4]] <- transInit(~ 1, nst = 6,data=dat,
                             pstart = rep(1/6, 6))
transition[[5]] <- transInit(~ 1, nst = 6,data=dat,
                             pstart = rep(1/6, 6))
transition[[6]] <- transInit(~ 1, nst = 6,data=dat,
                             pstart = rep(1/6, 6))
inMod <- transInit(~ 1, ns = 6, pstart = rep(1/6, 6),
                   data = data.frame(1))
mod6 <- makeDepmix(response = rModels, transition = transition,
                   prior=inMod,homogeneous = FALSE)
wbvmhmm6s <- fit(mod6, verbose = TRUE, emc=em.control(rand=FALSE, maxit=460))



# HMM WVM 7 states ----

rModels <- list()
rModels[[1]] <- list()
rModels[[1]][[1]] <- weibull(dist, pstart = c(0.5,1),data=dat)
rModels[[1]][[2]] <- vonMises(dat$Turningangle, pstart = c(0, 1),data=dat)

rModels[[2]] <- list()
rModels[[2]][[1]] <- weibull(dist, pstart = c(1,1),data=dat)
rModels[[2]][[2]] <- vonMises(dat$Turningangle, pstart = c(1, 1),data=dat)

rModels[[3]] <- list()
rModels[[3]][[1]] <- weibull(dist, pstart = c(2,5),data=dat)
rModels[[3]][[2]] <- vonMises(dat$Turningangle, pstart = c(2, 1))

rModels[[4]] <- list()
rModels[[4]][[1]] <- weibull(dist, pstart = c(3,1),data=dat)
rModels[[4]][[2]] <- vonMises(dat$Turningangle, pstart = c(3, 1))

rModels[[5]] <- list()
rModels[[5]][[1]] <- weibull(dist, pstart = c(5,5),data=dat)
rModels[[5]][[2]] <- vonMises(dat$Turningangle, pstart = c(4, 1))

rModels[[6]] <- list()
rModels[[6]][[1]] <- weibull(dist, pstart = c(4,2),data=dat)
rModels[[6]][[2]] <- vonMises(dat$Turningangle, pstart = c(5, 1))

rModels[[7]] <- list()
rModels[[7]][[1]] <- weibull(dist, pstart = c(5,2),data=dat)
rModels[[7]][[2]] <- vonMises(dat$Turningangle, pstart = c(6, 1))

trstart <- rep(1/7,49)
transition <- list()
transition[[1]] <- transInit(~ 1, nst = 7, data=dat,
                             pstart = rep(1/7, 7))
transition[[2]] <- transInit(~ 1, nst = 7,data=dat,
                             pstart = rep(1/7, 7))
transition[[3]] <- transInit(~ 1, nst = 7,data=dat,
                             pstart = rep(1/7, 7))
transition[[4]] <- transInit(~ 1, nst = 7,data=dat,
                             pstart = rep(1/7, 7))
transition[[5]] <- transInit(~ 1, nst = 7,data=dat,
                             pstart = rep(1/7, 7))
transition[[6]] <- transInit(~ 1, nst = 7,data=dat,
                             pstart = rep(1/7, 7))
transition[[7]] <- transInit(~ 1, nst = 7,data=dat,
                             pstart = rep(1/7, 7))
inMod <- transInit(~ 1, ns = 7, pstart = rep(1/7, 7),
                   data = data.frame(1))
mod7 <- makeDepmix(response = rModels, transition = transition,
                   prior=inMod,homogeneous = FALSE)
wbvmhmm7s <- fit(mod7, verbose = TRUE, emc=em.control(rand=FALSE, maxit=460))






sumdf <- function(lst){
  BIC <- ldply(lst,BIC)
  nstates <-ldply(lst,nstates)
  para <- ldply(lst,freepars)
  model <- c('HMM WeibullVM','HMM WeibullVM','HMM WeibullVM','HMM WeibullVM','HMM WeibullVM')
  type <- c('HMM','HMM','HMM','HMM','HMM')
  temp <- data.frame(BICS=BIC$V1,nstates=nstates$V1,parameters=para$V1,model,type)
  return(temp)
}

fitlist <- list(wbvmhmm3s,wbvmhmm4s,wbvmhmm5s,wbvmhmm6s,wbvmhmm7s)

catsum <- sumdf(fitlist)


saveRDS(catsum,file=paste("cat",catid,"wbvmhmm","RDS",sep="."))

