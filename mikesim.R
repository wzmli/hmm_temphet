simhmm <- function (object, nsim = 1, seed = NULL, ...) 
{
  if (!is.null(seed)) 
    set.seed(seed)
  ntim <- ntimes(object)
  nt <- sum(ntim)
  lt <- length(ntim)
  et <- cumsum(ntim)
  bt <- c(1, et[-lt] + 1)
  nr <- nresp(object)
  ns <- nstates(object)
  states <- array(, dim = c(nt, nsim))
  states[bt, ] <- simulate(object@prior, nsim = nsim, is.prior = TRUE)
  sims <- array(, dim = c(nt, ns, nsim))
  for (i in 1:ns) {
    if (depmixS4:::is.homogeneous(object)) {
      sims[, i, ] <- simulate(object@transition[[i]], nsim = nsim, 
                              times = rep(1, nt))
    }
    else {
      sims[, i, ] <- simulate(object@transition[[i]], nsim = nsim)
    }
  }
  for (case in 1:lt) {
    for (i in (bt[case] + 1):et[case]) {
      states[i, ] <- sims[cbind(i-1, states[i - 1, ], 1:nsim)]
    }
  }
  states <- as.vector(states)
  responses <- list(length = nr)
  for (i in 1:nr) {
    tmp <- matrix(, nrow = nt * nsim, ncol = NCOL(object@response[[1]][[i]]@y))
    for (j in 1:ns) {
      tmp[states == j, ] <- simulate(object@response[[j]][[i]], 
                                     nsim = nsim)[states == j, ]
    }
    responses[[i]] <- tmp
  }
  object <- as(object, "depmix.sim")
  object@states <- as.matrix(states)
  object@prior@x <- as.matrix(apply(object@prior@x, 2, rep, 
                                    nsim))
  for (j in 1:ns) {
    if (!depmixS4:::is.homogeneous(object)) 
      object@transition[[j]]@x <- as.matrix(apply(object@transition[[j]]@x, 
                                                  2, rep, nsim))
    for (i in 1:nr) {
      object@response[[j]][[i]]@y <- as.matrix(responses[[i]])
      object@response[[j]][[i]]@x <- as.matrix(apply(object@response[[j]][[i]]@x, 
                                                     2, rep, nsim))
    }
  }
  object@ntimes <- rep(object@ntimes, nsim)
  nt <- sum(object@ntimes)
  if (depmixS4:::is.homogeneous(object)) 
    trDens <- array(0, c(1, ns, ns))
  else trDens <- array(0, c(nt, ns, ns))
  dns <- array(, c(nt, nr, ns))
  for (i in 1:ns) {
    for (j in 1:nr) {
      dns[, j, i] <- dens(object@response[[i]][[j]])
    }
    trDens[, , i] <- dens(object@transition[[i]])
  }
  object@init <- dens(object@prior)
  object@trDens <- trDens
  object@dens <- dns
  return(object)
}