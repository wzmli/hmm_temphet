
R := /usr/bin/env Rscript

cat.%.df.Rout: dataframe.R ArchivePantherData.csv
	$(run-R)

cat.1.%.Rout: cat.1.df.Rout simulate.R ./models/%.R 
	$(run-R) 

cat.2.%.Rout: cat.2.df.Rout simulate.R ./models/%.R 
	$(run-R) 

cat.14.%.Rout: cat.14.df.Rout simulate.R ./models/%.R 
	$(run-R) 

cat.15.%.Rout: cat.15.df.Rout simulate.R ./models/%.R 
	$(run-R) 


simtest.Rout: simulate.R simfunctions.R simtest.R 101.txt
	$(run-R)

sim.%.Rout: simulate.R simfunctions.R simtest.R
	$(run-R)

simxy.%.Rout: simulate.R simfunctions.R simxygps.R
	$(run-R)

simtime.%.Rout: simulate.R simfunctions.R simtimetest.R
	$(run-R)

simsin.%.Rout: simulate.R simfunctions.R simsin.R
	$(run-R)

plotsimtest.Rout: plotsimtest.R
	$(run-R)

clean:
	rm -f *.bbl *.blg *.log *.aux *.loc *~ *.txt

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
