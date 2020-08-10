#########################################################################
##       R PROGRAM: PSE_estimate_LLM-2.R
##
##         PROJECT: Evaluation of multiple-encounter population size
##                  estimators
##
##     DESCRIPTION: Estimate population sizes from simulated samples
##                  using loglinear models as implemented in the
##                  Rcapture package
##
##                  Mt, Mbh and Mht samples
##
##      WRITTEN BY: Steve Gutreuter, CDC/CGH/DGHT Statistics, Estimation
##                                   and Modeling Team
##                  sgutreuter@cdc.gov
##
##            D1ATE: 2020-06-15
##
#########################################################################
library(Rcapture)
basepath <- file.path(Sys.getenv("PROJ"), "PSE/PSEsim")
workpath <- file.path(basepath, "R")
datapath <- file.path(basepath, "data")
setwd(workpath)
source(file.path(workpath, "PSE_sim_functions.R"))
dataname <- file.path(datapath, "PSEsimSamples_alt.Rdata")
attach(dataname)
cgwtools::lsdata(file.path(datapath, "PSEsimSamples_alt.Rdata"))

#########################################################################
## loglinear model fitting to independent homogeneous capture
## probabilities
#########################################################################
## test <- llm.estimate(data = M0.10K.025, events = 2:5,
##                      nreps = 2)
## print(test)
## Mt data
t.start0 <- t.start <- Sys.time()
llm.Mt.10K.025 <- llm.estimate(data = Mt.10K.025, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mt.10K.050 <- llm.estimate(data = Mt.10K.050, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mt.10K.100 <- llm.estimate(data = Mt.10K.100, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mt.10K.150 <- llm.estimate(data = Mt.10K.150, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start
## Mbh data
t.start <- Sys.time()
llm.Mbh.10K.025 <- llm.estimate(data = Mbh.10K.025, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mbh.10K.050 <- llm.estimate(data = Mbh.10K.050, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mbh.10K.100 <- llm.estimate(data = Mbh.10K.100, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mbh.10K.150 <- llm.estimate(data = Mbh.10K.150, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start
## Mht data
t.start <- Sys.time()
llm.Mht.10K.025 <- llm.estimate(data = Mht.10K.025, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mht.10K.050 <- llm.estimate(data = Mht.10K.050, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mht.10K.100 <- llm.estimate(data = Mht.10K.100, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mht.10K.150 <- llm.estimate(data = Mht.10K.150, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start
Sys.time() - t.start0

rm(t.start, t.start0)
#########################################################################
## Save the results
#########################################################################
save(llm.Mt.10K.025, llm.Mt.10K.050, llm.Mt.10K.100, llm.Mt.10K.150,
     llm.Mbh.10K.025, llm.Mbh.10K.050, llm.Mbh.10K.100, llm.Mbh.10K.150,
     llm.Mht.10K.025, llm.Mht.10K.050, llm.Mht.10K.100, llm.Mht.10K.150,
     file = file.path(datapath, "LLM_estimates-2.Rdata"))
cgwtools::lsdata(file.path(datapath, "LLM_estimates-2.Rdata"))
################################   END of FILE   ########################
