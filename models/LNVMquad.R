
library(depmixS4)
library(methods)
library(circular)
library(MASS)
library(stats)
library(dplyr)
library(plyr)

catid <- unlist(strsplit(rtargetname,split="[.]"))[2]

## library(VGAM)  ## this may mess up multinomial()

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
##****************************
## HMM LNVM 3 states ----
rModels <- list()
rModels[[1]] <- list()
rModels[[1]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(3, 1))
rModels[[1]][[2]] <- vonMises(dat$Turningangle, pstart = c(0, 1),data=dat)

rModels[[2]] <- list()
rModels[[2]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(5, 1))
rModels[[2]][[2]] <- vonMises(dat$Turningangle, pstart = c(1, 1),data=dat)

rModels[[3]] <- list()
rModels[[3]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(7, 1))
rModels[[3]][[2]] <- vonMises(dat$Turningangle, pstart = c(3, 1))

trstart <- rep(1/3,9)
transition <- list()
transition[[1]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 3, data=dat,
                             pstart = c(1/3,1/3,1/3,0,1/3,1/3,0,1/3,1/3))
transition[[2]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 3,data=dat,
                             pstart = c(1/3,1/3,1/3,0,1/3,1/3,0,1/3,1/3))
transition[[3]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 3,data=dat,
                             pstart = c(1/3,1/3,1/3,0,1/3,1/3,0,1/3,1/3))
inMod <- transInit(~ 1, ns = 3, pstart = rep(1/3, 3),
                   data = data.frame(1))
mod3 <- makeDepmix(response = rModels, transition = transition,
                   prior=inMod,homogeneous = FALSE)
lnvmhmmquad3s <- fit(mod3, verbose = TRUE, emc=em.control(rand=FALSE, maxit=460))


# HMM LNVM 4 states ----
rModels <- list()
rModels[[1]] <- list()
rModels[[1]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(3, 1))
rModels[[1]][[2]] <- vonMises(dat$Turningangle, pstart = c(0, 1),data=dat)

rModels[[2]] <- list()
rModels[[2]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(4, 1))
rModels[[2]][[2]] <- vonMises(dat$Turningangle, pstart = c(1, 1),data=dat)

rModels[[3]] <- list()
rModels[[3]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(5, 1))
rModels[[3]][[2]] <- vonMises(dat$Turningangle, pstart = c(2, 1))

rModels[[4]] <- list()
rModels[[4]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(6, 1))
rModels[[4]][[2]] <- vonMises(dat$Turningangle, pstart = c(3, 1))

trstart <- rep(1/4,16)
transition <- list()
transition[[1]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 4, data=dat,
                             pstart = c(1/4,1/4,1/4,1/4,
                                        0,1/4,1/4,1/4,
                                        0,1/4,1/4,1/4))
transition[[2]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 4, data=dat,
                             pstart = c(1/4,1/4,1/4,1/4,
                                        0,1/4,1/4,1/4,
                                        0,1/4,1/4,1/4))
transition[[3]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 4, data=dat,
                             pstart = c(1/4,1/4,1/4,1/4,
                                        0,1/4,1/4,1/4,
                                        0,1/4,1/4,1/4))
transition[[4]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 4, data=dat,
                             pstart = c(1/4,1/4,1/4,1/4,
                                        0,1/4,1/4,1/4,
                                        0,1/4,1/4,1/4))
inMod <- transInit(~ 1, ns = 4, pstart = rep(1/4, 4),
                   data = data.frame(1))
mod4 <- makeDepmix(response = rModels, transition = transition,
                   prior=inMod,homogeneous = FALSE)
lnvmhmmquad4s <- fit(mod4, verbose = TRUE, emc=em.control(rand=FALSE, maxit=460))

# HMM LNVM 5 states ----

rModels <- list()
rModels[[1]] <- list()
rModels[[1]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(3, 1))
rModels[[1]][[2]] <- vonMises(dat$Turningangle, pstart = c(0, 1),data=dat)

rModels[[2]] <- list()
rModels[[2]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(4, 1))
rModels[[2]][[2]] <- vonMises(dat$Turningangle, pstart = c(1, 1),data=dat)

rModels[[3]] <- list()
rModels[[3]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(5, 1))
rModels[[3]][[2]] <- vonMises(dat$Turningangle, pstart = c(2, 1))

rModels[[4]] <- list()
rModels[[4]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(6, 1))
rModels[[4]][[2]] <- vonMises(dat$Turningangle, pstart = c(3, 1))

rModels[[5]] <- list()
rModels[[5]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(7, 1))
rModels[[5]][[2]] <- vonMises(dat$Turningangle, pstart = c(4, 1))

trstart <- rep(1/5,25)
transition <- list()
transition[[1]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 5, data=dat,
                             pstart = c(1/5,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5))
transition[[2]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 5, data=dat,
                             pstart = c(1/5,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5))
transition[[3]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 5, data=dat,
                             pstart = c(1/5,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5))
transition[[4]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 5, data=dat,
                             pstart = c(1/5,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5))
transition[[5]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 5, data=dat,
                             pstart = c(1/5,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5))
inMod <- transInit(~ 1, ns = 5, pstart = rep(1/5, 5),
                   data = data.frame(1))
mod5 <- makeDepmix(response = rModels, transition = transition,
                   prior=inMod,homogeneous = FALSE)
lnvmhmmquad5s <- fit(mod5, verbose = TRUE, emc=em.control(rand=FALSE, maxit=460))

# HMM LNVM 6 states ----

rModels <- list()
rModels[[1]] <- list()
rModels[[1]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(3, 1))
rModels[[1]][[2]] <- vonMises(dat$Turningangle, pstart = c(0, 1),data=dat)

rModels[[2]] <- list()
rModels[[2]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(4, 1))
rModels[[2]][[2]] <- vonMises(dat$Turningangle, pstart = c(1, 1),data=dat)

rModels[[3]] <- list()
rModels[[3]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(5, 1))
rModels[[3]][[2]] <- vonMises(dat$Turningangle, pstart = c(2, 1))

rModels[[4]] <- list()
rModels[[4]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(6, 1))
rModels[[4]][[2]] <- vonMises(dat$Turningangle, pstart = c(3, 1))

rModels[[5]] <- list()
rModels[[5]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(7, 1))
rModels[[5]][[2]] <- vonMises(dat$Turningangle, pstart = c(4, 1))

rModels[[6]] <- list()
rModels[[6]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(8, 1))
rModels[[6]][[2]] <- vonMises(dat$Turningangle, pstart = c(5, 1))

trstart <- rep(1/6,36)
transition <- list()
transition[[1]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 6, data=dat,
                             pstart = c(1/6,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6))
transition[[2]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 6, data=dat,
                             pstart = c(1/6,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6))
transition[[3]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 6, data=dat,
                             pstart = c(1/6,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6))
transition[[4]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 6, data=dat,
                             pstart = c(1/6,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6))
transition[[5]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 6, data=dat,
                             pstart = c(1/6,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6))
transition[[6]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 6, data=dat,
                             pstart = c(1/6,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6))
inMod <- transInit(~ 1, ns = 6, pstart = rep(1/6, 6),
                   data = data.frame(1))
mod6 <- makeDepmix(response = rModels, transition = transition,
                   prior=inMod,homogeneous = FALSE)
lnvmhmmquad6s <- fit(mod6, verbose = TRUE, emc=em.control(rand=FALSE, maxit=460))

# HMM LNVM 7 states ----

rModels <- list()
rModels[[1]] <- list()
rModels[[1]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(3, 1))
rModels[[1]][[2]] <- vonMises(dat$Turningangle, pstart = c(0, 1),data=dat)

rModels[[2]] <- list()
rModels[[2]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(4, 1))
rModels[[2]][[2]] <- vonMises(dat$Turningangle, pstart = c(1, 1),data=dat)

rModels[[3]] <- list()
rModels[[3]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(5, 1))
rModels[[3]][[2]] <- vonMises(dat$Turningangle, pstart = c(2, 1))

rModels[[4]] <- list()
rModels[[4]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(6, 1))
rModels[[4]][[2]] <- vonMises(dat$Turningangle, pstart = c(3, 1))

rModels[[5]] <- list()
rModels[[5]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(7, 1))
rModels[[5]][[2]] <- vonMises(dat$Turningangle, pstart = c(4, 1))

rModels[[6]] <- list()
rModels[[6]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(8, 1))
rModels[[6]][[2]] <- vonMises(dat$Turningangle, pstart = c(5, 1))

rModels[[7]] <- list()
rModels[[7]][[1]] <- GLMresponse(formula = LogDist ~ 1, data = dat,
                                 family = gaussian(), pstart = c(9, 1))
rModels[[7]][[2]] <- vonMises(dat$Turningangle, pstart = c(6, 1))

trstart <- rep(1/7,49)
transition <- list()
transition[[1]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 7, data=dat,
                             pstart = c(1/7,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7))
transition[[2]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 7, data=dat,
                             pstart = c(1/7,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7))
transition[[3]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 7, data=dat,
                             pstart = c(1/7,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7))
transition[[4]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 7, data=dat,
                             pstart = c(1/7,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7))
transition[[5]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 7, data=dat,
                             pstart = c(1/7,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7))
transition[[6]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 7, data=dat,
                             pstart = c(1/7,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7))
transition[[7]] <- transInit(~ I(Time/24)+I((Time/24)^2), nst = 7, data=dat,
                             pstart = c(1/7,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7))
inMod <- transInit(~ 1, ns = 7, pstart = rep(1/7, 7),
                   data = data.frame(1))
mod7 <- makeDepmix(response = rModels, transition = transition,
                   prior=inMod,homogeneous = FALSE)
lnvmhmmquad7s <- fit(mod7, verbose = TRUE, emc=em.control(rand=FALSE, maxit=460))


sumdf <- function(lst){
  BIC <- ldply(lst,BIC)
  nstates <-ldply(lst,nstates)
  para <- ldply(lst,freepars)
  model <- c('HMM LNVMquad','HMM LNVMquad','HMM LNVMquad','HMM LNVMquad','HMM LNVMquad')
  type <- c('HMM + TH','HMM + TH','HMM + TH','HMM + TH','HMM + TH')
  temp <- data.frame(BICS=BIC$V1,nstates=nstates$V1,parameters=para$V1,model,type)
  return(temp)
}

fitlist <- list(lnvmhmmquad3s,lnvmhmmquad4s,lnvmhmmquad5s,lnvmhmmquad6s,lnvmhmmquad7s)

catsum <- sumdf(fitlist)


saveRDS(catsum,file=paste("cat",catid,"lnvmhmmquad","RDS",sep="."))




