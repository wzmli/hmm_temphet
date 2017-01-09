
library(depmixS4)
library(methods)
library(circular)
library(MASS)
library(stats)
library(dplyr)
library(plyr)

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



simlnvm <- function(dat,fit){
  tempmodel <- depmix(LogDist~1,
                  data=dat,
                  nstate=nstates(fit),
                  transition=fit@transition[[1]]@formula,
                  family=gaussian())
  parindex <- length(getpars(tempmodel))
  parind <- parindex - 2*nstates(fit) 
  logitpars <- c(getpars(fit)[1:nstates(fit)])
  pars <- exp(logitpars)/sum(exp(logitpars))
  model<-setpars(tempmodel,c(pars,getpars(fit)[(nstates(fit)+1):parindex]))
  if(fit@transition[[1]]@formula == "~1"){
    for(i in 1:(nstates(fit))){
      templogit <- fit@transition[[i]]@parameters$coefficients
      invlogit <- exp(templogit)/sum(exp(templogit))
      pars <- c(pars,invlogit)
    }
    newpars <- c(pars,tail(getpars(fit),2*nstates(fit)))
    model<-setpars(tempmodel,newpars)
  }
  sim <- simhmm(model)
  df <- data.frame(states=sim@states,obs=dat$LogDist,time=dat$Time)
  df2 <- (df %>% rowwise()
    %>% dplyr::mutate(
      steplength=rnorm(1,getpars(fit)[parind+4*(states-1)+1],getpars(fit)[parind+4*(states-1)+2])
      , angle=rvonmises(1,getpars(fit)[parind+4*(states-1)+3],getpars(fit)[parind+4*(states-1)+4])
      )
  )
  return(df2)
  }



vitsimlnvm <- function(dat,fit){
  tempmodel <- depmix(LogDist~1,
                      data=dat,
                      nstate=nstates(fit),
                      transition=fit@transition[[1]]@formula,
                      family=gaussian())
  parindex <- length(getpars(tempmodel))
  parind <- parindex - 2*nstates(fit) 
  vit <- viterbi(fit)
  df <- data.frame(states=vit$state, obs=dat$LogDist,time=dat$Time)
  df2 <- (df %>% rowwise()
          %>% dplyr::mutate(
            steplength=rnorm(1,getpars(fit)[parind+4*(states-1)+1],getpars(fit)[parind+4*(states-1)+2])
            , angle=rvonmises(1,getpars(fit)[parind+4*(states-1)+3],getpars(fit)[parind+4*(states-1)+4])
          )
  )
  return(df2)
}
