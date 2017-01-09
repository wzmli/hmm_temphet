load("cat11.RData")
library("depmixS4")
source("mikesim.R")  ## replace stuff from depmixS4
source("Fitting_and_Simulation_code.R")
statevec <- 2:3
fitList <- lapply(statevec,fithomo,cat=cat11)

fitModList <- list(homo=fithomo,sin=fitsin,quad=fitquad,block=fitblock,
                   mix=fitmix,block=fitblock)
fitAll <- function(dat,numstates) {
    lapply(fitModList, function(f) f(numstates,dat))
}
fitAllStates <- function(dat,statevec=3:6) {
    lapply(statevec,fitAll,dat=dat)
}
fitAllCats <- function() {
    cats
}
## instead of getting a whole bunch of objects: cat27fithomo5, cat11fitquad4, ....
## in the end we're going to get a big nested list
## cat {
##    state {
##        model { }
##    }
## }


for (cat in catVec) {
    for (state in stateVec) {
        for (model in modelVec) {
            ## fit the right thing, then save it
            saveRDS(x,file=paste0(paste("output/fit",cat,state,model,sep="_"),".rds"))
        }
    }
}

## in R, want to work with lists of things, not giant piles of objects
## i.e. use lapply(), and lists, not lots of get/assign

## outside of R, make a reasonable compromise between fitting everything all at once
## in one giant output file and fitting everything in tiny pieces with a million output
## files

## save intermediate files if they're not too big

## don't worry *too* much about optimizing everything

## don't store intermediate results on the repo unless they're
## (1) sufficiently hard to compute that people shouldn't have to regenerate them
## (2) sufficiently small that it's not horribly awkward to download them

allFiles <- list.files("output")
library("stringr")
## maybe use str_extract from stringr package instead?
##cats <- unique(gsub("fit_([0-9]+)_.*","\\1",allFiles))

ss <- strsplit(allFiles,"_")
sumData <- data.frame(cat=sapply(ss,"[[",2),  ## cats
                      state=as.numeric(lapply(ss,"[[",3)),  ## states
                      model=gsub("\\.rds","",lapply(ss,"[[",4)),
                      filename=allFiles)


allDat <- list()
for (cat in levels(sumData$cat)) {
    allDat[[cat]] <- list()
    tmpsum <- droplevels(subset(sumData,cat==cat))
    for (model in levels(tmpsum$model)) {
        allDat[[cat]][[model]] <- list()
        tmpsum2 <- droplevels(subset(tmpsum,model==model))
        for (state in unique(tmpsum2$state)) {
            tmpsum3 <- subset(tmpsum2,state==state)
            cat("loading ",tmpsum3$filename,"\n")
            ## how are these stored??
            ## might want to use saveRDS() / readRDS()
            allDat[[cat]][[model]][[state]] <- readRDS(tmpsum3$filename)

        }
    }
}
    
    
