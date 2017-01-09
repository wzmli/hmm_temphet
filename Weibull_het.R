## WEIBULL VM 

library(depmixS4)
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

# depmix model setup 3 states time het ----
dist <- pmax(cat$Distance,1e-3)
rModels <- list()
rModels[[1]] <- list()
rModels[[1]][[1]] <- weibull(dist, pstart = c(0.5,0.5),data=cat)

rModels[[2]] <- list()
rModels[[2]][[1]] <- weibull(dist, pstart = c(1.5,1.5),data=cat)

rModels[[3]] <- list()
rModels[[3]][[1]] <- weibull(dist, pstart = c(1,4),data=cat)

trstart <- rep(1/3,9)
transition <- list()
transition[[1]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 3, data=cat,
                             pstart = c(1/3,1/3,1/3,0,1/3,1/3,0,1/3,1/3))
transition[[2]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 3,data=cat,
                             pstart = c(1/3,1/3,1/3,0,1/3,1/3,0,1/3,1/3))
transition[[3]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 3,data=cat,
                             pstart = c(1/3,1/3,1/3,0,1/3,1/3,0,1/3,1/3))
inMod <- transInit(~ 1, ns = 3, pstart = rep(1/3, 3),
                   data = data.frame(1))
ws3 <- makeDepmix(response = rModels, transition = transition,
                   prior=inMod,homogeneous = FALSE)
Wsin3 <- fit(ws3, verbose = TRUE, emc=em.control(rand=FALSE, maxit=460))
summary(Wsin3)

# HMM W 4 states time het ----
rModels <- list()
rModels[[1]] <- list()
rModels[[1]][[1]] <- weibull(dist, pstart = c(0.5,1),data=cat)

rModels[[2]] <- list()
rModels[[2]][[1]] <- weibull(dist, pstart = c(1,2),data=cat)

rModels[[3]] <- list()
rModels[[3]][[1]] <- weibull(dist, pstart = c(3,3),data=cat)

rModels[[4]] <- list()
rModels[[4]][[1]] <- weibull(dist, pstart = c(0.6,5),data=cat)

trstart <- rep(1/4,16)
transition <- list()
transition[[1]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 4, data=cat,
                             pstart = c(1/4,1/4,1/4,1/4,
                                        0,1/4,1/4,1/4,
                                        0,1/4,1/4,1/4))
transition[[2]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 4, data=cat,
                             pstart = c(1/4,1/4,1/4,1/4,
                                        0,1/4,1/4,1/4,
                                        0,1/4,1/4,1/4))
transition[[3]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 4, data=cat,
                             pstart = c(1/4,1/4,1/4,1/4,
                                        0,1/4,1/4,1/4,
                                        0,1/4,1/4,1/4))
transition[[4]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 4, data=cat,
                             pstart = c(1/4,1/4,1/4,1/4,
                                        0,1/4,1/4,1/4,
                                        0,1/4,1/4,1/4))
inMod <- transInit(~ 1, ns = 4, pstart = rep(1/4, 4),
                   data = data.frame(1))
ws4 <- makeDepmix(response = rModels, transition = transition,
                   prior=inMod,homogeneous = FALSE)
Wsin4 <- fit(ws4, verbose = TRUE, emc=em.control(rand=FALSE, maxit=460))
summary(Wsin4)

# HMM W 5 states time het ----

rModels <- list()
rModels[[1]] <- list()
rModels[[1]][[1]] <- weibull(dist, pstart = c(0.5,1),data=cat)

rModels[[2]] <- list()
rModels[[2]][[1]] <- weibull(dist, pstart = c(1,1),data=cat)

rModels[[3]] <- list()
rModels[[3]][[1]] <- weibull(dist, pstart = c(1,2),data=cat)
rModels[[4]] <- list()
rModels[[4]][[1]] <- weibull(dist, pstart = c(2,3),data=cat)

rModels[[5]] <- list()
rModels[[5]][[1]] <- weibull(dist, pstart = c(2,5),data=cat)

trstart <- rep(1/5,25)
transition <- list()
transition[[1]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 5, data=cat,
                             pstart = c(1/5,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5))
transition[[2]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 5, data=cat,
                             pstart = c(1/5,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5))
transition[[3]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 5, data=cat,
                             pstart = c(1/5,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5))
transition[[4]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 5, data=cat,
                             pstart = c(1/5,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5))
transition[[5]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 5, data=cat,
                             pstart = c(1/5,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5,
                                        0,1/5,1/5,1/5,1/5))
inMod <- transInit(~ 1, ns = 5, pstart = rep(1/5, 5),
                   data = data.frame(1))
ws5 <- makeDepmix(response = rModels, transition = transition,
                   prior=inMod,homogeneous = FALSE)
Wsin5 <- fit(ws5, verbose = TRUE, emc=em.control(rand=FALSE, maxit=460))
summary(Wsin5)

# HMM W 6 states time het ----

rModels <- list()
rModels[[1]] <- list()
rModels[[1]][[1]] <- weibull(dist, pstart = c(0.5,1),data=cat)

rModels[[2]] <- list()
rModels[[2]][[1]] <- weibull(dist, pstart = c(1,1),data=cat)

rModels[[3]] <- list()
rModels[[3]][[1]] <- weibull(dist, pstart = c(2,5),data=cat)

rModels[[4]] <- list()
rModels[[4]][[1]] <- weibull(dist, pstart = c(3,1),data=cat)

rModels[[5]] <- list()
rModels[[5]][[1]] <- weibull(dist, pstart = c(5,5),data=cat)

rModels[[6]] <- list()
rModels[[6]][[1]] <- weibull(dist, pstart = c(4,2),data=cat)

trstart <- rep(1/6,36)
transition <- list()
transition[[1]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 6, data=cat,
                             pstart = c(1/6,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6))
transition[[2]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 6, data=cat,
                             pstart = c(1/6,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6))
transition[[3]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 6, data=cat,
                             pstart = c(1/6,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6))
transition[[4]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 6, data=cat,
                             pstart = c(1/6,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6))
transition[[5]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 6, data=cat,
                             pstart = c(1/6,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6))
transition[[6]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 6, data=cat,
                             pstart = c(1/6,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6,
                                        0,1/6,1/6,1/6,1/6,1/6))
inMod <- transInit(~ 1, ns = 6, pstart = rep(1/6, 6),
                   data = data.frame(1))
ws6 <- makeDepmix(response = rModels, transition = transition,
                   prior=inMod,homogeneous = FALSE)
Wsin6 <- fit(ws6, verbose = TRUE, emc=em.control(rand=FALSE, maxit=460))
summary(Wsin6)


# HMM W 7 states time het ----

rModels <- list()
rModels[[1]] <- list()
rModels[[1]][[1]] <- weibull(dist, pstart = c(0.5,1),data=cat)

rModels[[2]] <- list()
rModels[[2]][[1]] <- weibull(dist, pstart = c(1,1),data=cat)

rModels[[3]] <- list()
rModels[[3]][[1]] <- weibull(dist, pstart = c(2,5),data=cat)

rModels[[4]] <- list()
rModels[[4]][[1]] <- weibull(dist, pstart = c(3,1),data=cat)

rModels[[5]] <- list()
rModels[[5]][[1]] <- weibull(dist, pstart = c(5,5),data=cat)

rModels[[6]] <- list()
rModels[[6]][[1]] <- weibull(dist, pstart = c(4,2),data=cat)

rModels[[7]] <- list()
rModels[[7]][[1]] <- weibull(dist, pstart = c(6,2),data=cat)

trstart <- rep(1/7,49)
transition <- list()
transition[[1]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 7, data=cat,
                             pstart = c(1/7,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7))
transition[[2]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 7, data=cat,
                             pstart = c(1/7,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7))
transition[[3]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 7, data=cat,
                             pstart = c(1/7,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7))
transition[[4]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 7, data=cat,
                             pstart = c(1/7,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7))
transition[[5]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 7, data=cat,
                             pstart = c(1/7,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7))
transition[[6]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 7, data=cat,
                             pstart = c(1/7,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7))
transition[[7]] <- transInit(~ cos(2*pi*Time/24)+ sin(2*pi*Time/24), nst = 7, data=cat,
                             pstart = c(1/7,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7,
                                        0,1/7,1/7,1/7,1/7,1/7,1/7))
inMod <- transInit(~ 1, ns = 7, pstart = rep(1/7, 7),
                   data = data.frame(1))
ws7 <- makeDepmix(response = rModels, transition = transition,
                  prior=inMod,homogeneous = FALSE)
Wsin7 <- fit(ws7, verbose = TRUE, emc=em.control(rand=FALSE, maxit=280))
summary(Wsin7)




