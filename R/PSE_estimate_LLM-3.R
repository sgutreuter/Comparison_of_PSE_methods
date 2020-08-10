#########################################################################
##       R PROGRAM: PSE_estimate_LLM-3.R
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
##            DATE: 2020-06-15
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
## Mbht data
t.start0 <- t.start <- Sys.time()
llm.Mbht.10K.025 <- llm.estimate(data = Mbht.10K.025, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mbht.10K.050 <- llm.estimate(data = Mbht.10K.050, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mbht.10K.100 <- llm.estimate(data = Mbht.10K.100, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mbht.10K.150 <- llm.estimate(data = Mbht.10K.150, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start
## tinc data
t.start <- Sys.time()
llm.tinc.10K.025 <- llm.estimate(data = tinc.10K.025, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.tinc.10K.050 <- llm.estimate(data = tinc.10K.050, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.tinc.10K.100 <- llm.estimate(data = tinc.10K.100, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.tinc.10K.150 <- llm.estimate(data = tinc.10K.150, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start
## tdec data
t.start <- Sys.time()
llm.tdec.10K.025 <- llm.estimate(data = tdec.10K.025, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.tdec.10K.050 <- llm.estimate(data = tdec.10K.050, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.tdec.10K.100 <- llm.estimate(data = tdec.10K.100, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.tdec.10K.150 <- llm.estimate(data = tdec.10K.150, events = 2:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start
Sys.time() - t.start0
rm(t.start, t.start0)
#########################################################################
## Save the results
#########################################################################
save(llm.Mbht.10K.025, llm.Mbht.10K.050, llm.Mbht.10K.100, llm.Mbht.10K.150,
     llm.tinc.10K.025, llm.tinc.10K.050, llm.tinc.10K.100, llm.tinc.10K.150,
     llm.tdec.10K.025, llm.tdec.10K.050, llm.tdec.10K.100, llm.tdec.10K.150,
     file = file.path(datapath, "LLM_estimates-3.Rdata"))
cgwtools::lsdata(file.path(datapath, "LLM_estimates-3.Rdata"))
################################   END of FILE   ########################
