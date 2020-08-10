#########################################################################
##       R PROGRAM: PSE_estimate_LLM-1.R
##
##         PROJECT: Evaluation of multiple-encounter population size
##                  estimators
##
##     DESCRIPTION: Estimate population sizes from simulated samples
##                  using loglinear models as implemented in the
##                  Rcapture package
##
##                  M0, Mh and Mb samples
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
## M0 data
tstart0 <- t.start <- Sys.time()
llm.M0.10K.025 <- llm.estimate(data = M0.10K.025, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.M0.10K.050 <- llm.estimate(data = M0.10K.050, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.M0.10K.100 <- llm.estimate(data = M0.10K.100, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.M0.10K.150 <- llm.estimate(data = M0.10K.150, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start
## Mh data
t.start <- Sys.time()
llm.Mh.10K.025 <- llm.estimate(data = Mh.10K.025, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mh.10K.050 <- llm.estimate(data = Mh.10K.050, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mh.10K.100 <- llm.estimate(data = Mh.10K.100, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mh.10K.150 <- llm.estimate(data = Mh.10K.150, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start
## Mb data
t.start <- Sys.time()
llm.Mb.10K.025 <- llm.estimate(data = Mb.10K.025, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mb.10K.050 <- llm.estimate(data = Mb.10K.050, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mb.10K.100 <- llm.estimate(data = Mb.10K.100, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mb.10K.150 <- llm.estimate(data = Mb.10K.150, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start
Sys.time() - t.start0

rm(t.start, t.start0)
#########################################################################
## Save the results
#########################################################################
save(llm.M0.10K.025, llm.M0.10K.050, llm.M0.10K.100, llm.M0.10K.150,
     llm.Mh.10K.025, llm.Mh.10K.050, llm.Mh.10K.100, llm.Mh.10K.150,
     llm.Mb.10K.025, llm.Mb.10K.050, llm.Mb.10K.100, llm.Mb.10K.150,
     file = file.path(datapath, "LLM_estimates-1.Rdata"))
cgwtools::lsdata(file.path(datapath, "LLM_estimates-1.Rdata"))
################################   END of FILE   ########################
