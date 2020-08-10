paste0("Start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", sep = "\n"))
################################################################################
##       R PROGRAM: PSE_estimate_LLM-hetero.R
##
##         PROJECT: Evaluation of multiple-encounter population size
##                  estimators
##
##     DESCRIPTION: Estimate population sizes from simulated samples using
##                  loglinear models as implemented in the Rcapture package.
##                  This implementation is designed compare methods for handling
##                  heterogeneity.
##
##                  Mt, Mbh, Mht, and Mbth samples
##
##      WRITTEN BY: Steve Gutreuter, CDC/CGH/DGHT Statistics, Estimation
##                                   and Modeling Team
##                  sgutreuter@cdc.gov
##
##            DATE: 2020-06-12
##
################################################################################
library(Rcapture)
basepath <- file.path(Sys.getenv("PROJ"), "PSE/PSEsim")
workpath <- file.path(basepath, "R")
datapath <- file.path(basepath, "data")
setwd(workpath)
source(file.path(workpath, "PSE_sim_functions.R"))
dataname <- file.path(datapath, "PSEsimSamples_alt.Rdata")
attach(dataname, name = "currdata")
cgwtools::lsdata(file.path(datapath, "PSEsimSamples_alt.Rdata"))

################################################################################
## loglinear model fitting to independent homogeneous capture probabilities
################################################################################
## Mh data
t.start0 <- t.start <- Sys.time()
llm.Mh.10K.025.het <- llm.estimate.hetero(data = Mh.10K.025, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mh.10K.050.het <- llm.estimate.hetero(data = Mh.10K.050, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mh.10K.100.het <- llm.estimate.hetero(data = Mh.10K.100, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mh.10K.150.het <- llm.estimate.hetero(data = Mh.10K.150, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start
## Mt data
tstart0 <- t.start <- Sys.time()
llm.Mt.10K.025.het <- llm.estimate.hetero(data = Mt.10K.025, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mt.10K.050.het <- llm.estimate.hetero(data = Mt.10K.050, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mt.10K.100.het <- llm.estimate.hetero(data = Mt.10K.100, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mt.10K.150.het <- llm.estimate.hetero(data = Mt.10K.150, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start
## Mbh data
t.start <- Sys.time()
llm.Mbh.10K.025.het <- llm.estimate.hetero(data = Mbh.10K.025, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mbh.10K.050.het <- llm.estimate.hetero(data = Mbh.10K.050, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mbh.10K.100.het <- llm.estimate.hetero(data = Mbh.10K.100, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mbh.10K.150.het <- llm.estimate.hetero(data = Mbh.10K.150, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start
## Mht data
t.start <- Sys.time()
llm.Mht.10K.025.het <- llm.estimate.hetero(data = Mht.10K.025, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mht.10K.050.het <- llm.estimate.hetero(data = Mht.10K.050, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mht.10K.100.het <- llm.estimate.hetero(data = Mht.10K.100, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mht.10K.150.het <- llm.estimate.hetero(data = Mht.10K.150, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start
## Mbht data
t.start <- Sys.time()
llm.Mbht.10K.025.het <- llm.estimate.hetero(data = Mbht.10K.025, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mbht.10K.050.het <- llm.estimate.hetero(data = Mbht.10K.050, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mbht.10K.100.het <- llm.estimate.hetero(data = Mbht.10K.100, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start

llm.Mbht.10K.150.het <- llm.estimate.hetero(data = Mbht.10K.150, events = 3:5,
                               nreps = 400, progress = FALSE)
Sys.time() - t.start
Sys.time() - t.start0

rm(t.start)
################################################################################
## Save the results
################################################################################
save(llm.Mh.10K.025.het, llm.Mh.10K.050.het, llm.Mh.10K.100.het,
     llm.Mh.10K.150.het,
     llm.Mt.10K.025.het, llm.Mt.10K.050.het, llm.Mt.10K.100.het,
     llm.Mt.10K.150.het,
     llm.Mbh.10K.025.het, llm.Mbh.10K.050.het, llm.Mbh.10K.100.het,
     llm.Mbh.10K.150.het,
     llm.Mht.10K.025.het, llm.Mht.10K.050.het, llm.Mht.10K.100.het,
     llm.Mht.10K.150.het,
     llm.Mbht.10K.025.het, llm.Mbht.10K.050.het, llm.Mbht.10K.100.het,
     llm.Mbht.10K.150.het,
     file = file.path(datapath, "LLM_estimates-hetero.Rdata"))
cgwtools::lsdata(file.path(datapath, "LLM_estimates-hetero.Rdata"))
################################   END of FILE   ###############################
