#########################################################################
##       R PROGRAM: PSE_estimate_BMA-1.R
##
##         PROJECT: Evaluation of multiple-encounter population size
##                  estimators
##
##     DESCRIPTION: Estimate population sizes from simulated samples
##                  using Bayesian loglinear model averaging  as
##                  implemented in the dga package
##
##                  M0 samples
##                  priorNmax = 100000
##
##      WRITTEN BY: Steve Gutreuter, CDC/CGH/DGHT Statistics, Estimation
##                                   and Modeling Team
##                  sgutreuter@cdc.gov
##
##            DATE: 2020-05-29
##
##      DISCLAIMER: Although this program has been used by the Centers
##                  for Disease Control & Prevention (CDC), no warranty,
##                  expressed or implied, is made by the CDC or the U.S.
##                  Government as to the accuracy and functioning of the
##                  program and related program material nor shall the
##                  fact of distribution constitute any such warranty,
##                  and no responsibility is assumed by the CDC in
##                  connection therewith.
##
#########################################################################
library(HDInterval)
library(tidyr)
library(dga)
data(graphs3)
data(graphs4)
data(graphs5)
basepath <- file.path(Sys.getenv("PROJ"), "PSE/PSEsim")
workpath <- file.path(basepath, "R")
datapath <- file.path(basepath, "data")
setwd(workpath)
source(file.path(workpath, "PSE_sim_functions.R"))
dataname <- file.path(datapath, "PSEsimSamples_alt.Rdata")
attach(dataname)

cgwtools::lsdata(file.path(datapath, "PSEsimSamples_alt.Rdata"))

#########################################################################
## BMA loglinear model fitting to independent M0 data
#########################################################################
##t.start <- Sys.time()
##test <- bma.estimate(data = M0.10K.025, events = 3:5,
##                     priorNmax = 30000, nreps = 1)
##Sys.time() - t.start
##print(test)

## Runs:
t.start0 <- t.start <- Sys.time()
print(t.start)
bma.M0.10K.025 <- bma.estimate(data = M0.10K.025, events = 3:5,
                               priorNmax = 100000, nreps = 400,
                               progress = FALSE)
Sys.time() - t.start
Sys.time()

t.start <- Sys.time()
bma.M0.10K.050 <- bma.estimate(data = M0.10K.050, events = 3:5,
                               priorNmax = 100000, nreps = 400,
                               progress = FALSE)
Sys.time() - t.start
Sys.time()

t.start <- Sys.time()
bma.M0.10K.100 <- bma.estimate(data = M0.10K.100, events = 3:5,
                               priorNmax = 100000, nreps = 400,
                               progress = FALSE)
Sys.time() - t.start
Sys.time()

t.start <- Sys.time()
bma.M0.10K.150 <- bma.estimate(data = M0.10K.150, events = 3:5,
                               priorNmax = 100000, nreps = 400,
                               progress = FALSE)
Sys.time() - t.start
## Total elapsed time:
Sys.time() - t.start0

#########################################################################
## Save the results
#########################################################################
save(bma.M0.10K.025, bma.M0.10K.050, bma.M0.10K.100, bma.M0.10K.150,
     file = file.path(datapath, "BMA_M0_estimates.Rdata"))
cgwtools::lsdata(file.path(datapath, "BMA_M0_estimates.Rdata"))
################################   END of FILE   ########################
