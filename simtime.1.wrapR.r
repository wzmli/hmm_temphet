# This file was generated automatically by wrapR.pl
# You probably don't want to edit it

rtargetname <- "simtime.1"
pdfname <- ".simtime.1.Rout.pdf"
csvname <- "simtime.1.Rout.csv"
rdsname <- "simtime.1.Rds"
pdf(pdfname)
# End RR preface

# Generated using wrapR file simtime.1.wrapR.r
source('mikesim.R', echo=TRUE)
source('simfunctions.R', echo=TRUE)
source('simtimetest.R', echo=TRUE)
# Wrapped output file simtime.1.wrapR.rout
# Begin RR postscript
warnings()
proc.time()

# If you see this in your terminal, the R script simtime.1.wrapR.r (or something it called) did not close properly
save(file=".simtime.1.RData", seed)

