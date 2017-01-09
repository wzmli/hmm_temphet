# This file was generated automatically by wrapR.pl
# You probably don't want to edit it

rtargetname <- "sim.1"
pdfname <- ".sim.1.Rout.pdf"
csvname <- "sim.1.Rout.csv"
rdsname <- "sim.1.Rds"
pdf(pdfname)
# End RR preface

# Generated using wrapR file sim.1.wrapR.r
source('mikesim.R', echo=TRUE)
source('simfunctions.R', echo=TRUE)
source('simtest.R', echo=TRUE)
# Wrapped output file sim.1.wrapR.rout
# Begin RR postscript
warnings()
proc.time()

# If you see this in your terminal, the R script sim.1.wrapR.r (or something it called) did not close properly
save(file=".sim.1.RData", seed)

