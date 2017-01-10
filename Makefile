###1 Take PantherDF, hack it and store it independently by cats :

R := /usr/bin/env Rscript

current:  paper3.pdf.now

#cat.%.df: dataframe.R ArchivePantherData.csv %.cat
#	$(R) $^

cat.%.df.Rout: dataframe.R ArchivePantherData.csv
	$(run-R)

cat.1.%.Rout: cat.1.df.Rout cat.1.RDS simulate.R ./models/%.R 
	$(run-R) 

cat.2.%.Rout: cat.2.df.Rout cat.2.RDS simulate.R ./models/%.R 
	$(run-R) 

cat.14.%.Rout: cat.14.df.Rout cat.14.RDS simulate.R ./models/%.R 
	$(run-R) 

cat.15.%.Rout: cat.15.df.Rout cat.15.RDS simulate.R ./models/%.R 
	$(run-R) 


%.fit.Rout: %.df.Rout fitfunctions.R mikesim.R simfunctions.R %seeds.R
	$(run-R)

%.summary.Rout: %.fit.Rout summary.R %summary.R
	$(run-R)

cat.%.lnvm.Rout: cat.%.df.Rout ./lnvm_results/cat.%.lnvmhmm.RDS ./lnvm_results/cat.%.lnvmhmmblock.RDS ./lnvm_results/cat.%.lnvmhmmsin.RDS ./lnvm_results/cat.%.lnvmhmmquad.RDS mikesim.R lnvm_setup.R plotsimfunctions.R lnvmplot.R
	$(run-R)

%.plot.Rout: %sumdat.RData %.RData plotsimfunctions.R %simplots.R
	$(run-R)

%.BIC.Rout: %.RData wrap.R
	$(run-R)

%.weibull.Rout: %.df.Rout weibull.R
	$(run-R)

%.weibullVM.Rout: %.df.Rout weibullVM.R
	$(run-R)

%.LNVM.Rout: %.RDS %.df.Rout ./models/LNVM.R
	$(run-R)

%.Weibull_het.Rout: %.df.Rout Weibull_het.R
	$(run-R)

simtest.Rout: mikesim.R simfunctions.R simtest.R 101.txt
	$(run-R)

sim.%.Rout: mikesim.R simfunctions.R simtest.R
	$(run-R)

simxy.%.Rout: mikesim.R simfunctions.R simxygps.R
	$(run-R)

simtime.%.Rout: mikesim.R simfunctions.R simtimetest.R
	$(run-R)

simsin.%.Rout: mikesim.R simfunctions.R simsin.R
	$(run-R)

paperplots.Rout: cat.1.paperdf.RDS cat.2.paperdf.RDS cat.14.paperdf.RDS cat.15.paperdf.RDS plotsimfunctions.R paperplots.R
	$(run-R)

plotsimtest.Rout: plotsimtest.R
	$(run-R)

%.tex: %.Rnw
	echo "library('knitr'); knit(\"$*.Rnw\")" | R --slave

%.pdf: %.tex paper.bib
	texi2dvi -p $*.tex

%.pdf.now: %.pdf
	evince paper3.pdf

clean:
	rm -f *.bbl *.blg *.log *.aux *.loc *~ *.txt

move_sum: 
	mv cat.*.*.RDS ./summary_stats/

rmsims:
	rm -f sim.*.Rout sim.*.Rlog sim.*.wrapR.r sim.*.wrapR.rout .sim.*.RData



#######Make stuff

ms = makestuff/

-include $(ms)/git.mk
-include $(ms)/visual.mk

-include $(ms)/wrapR.mk
-include $(ms)/oldlatex.mk

makestuff:
	git clone https://github.com/dushoff/makestuff.git
