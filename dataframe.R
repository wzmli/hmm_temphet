library(depmixS4)
library(dplyr)

# files <- commandArgs(trailingOnly = TRUE)
# catid <- unlist(strsplit(files[2],split="[.]"))[1]
rawdat <- read.csv(input_files)
catid <- unlist(strsplit(rtargetname,split="[.]"))[2]


block <- function(t){
  if(t %in% 7:16)return("block2")
  if(t %in% 17:20)return("block3")
  return("block1")
}

## creating suitable dataframe for depmixS4
rawdf <- (rawdat 
  %>% rowwise() 
  %>% transmute(cat = animal_id
    , Sex = Sex
    , Time = as.numeric(Time) - 1 ##as.numeric conversion
    , Distance = Steplength.m.
    , LogDist = log10(Distance)
    , Turningangle = Turningangle.rad.
    , Block = block(Time)
    )
)

rawdf$LogDist[is.infinite(rawdf$LogDist)] <- NA

catnum <- c(catid)

dat <- rawdf %>% filter(cat %in% catnum)


saveRDS(dat,file=paste("cat",catid,"RDS",sep="."))
