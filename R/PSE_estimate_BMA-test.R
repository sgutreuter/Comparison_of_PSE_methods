#########################################################################
##       R PROGRAM: PSE_estimate_BMA-test.R
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
##            DATE: 2020-06-01
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
t.start <- Sys.time()
test <- bma.estimate(data = M0.10K.025, events = 3:4,
                     priorNmax = 30000, nreps = 1)
Sys.time() - t.start
print(test)


#########################################################################
## Save the results
#########################################################################
save(test, file = file.path(datapath, "BMA_test.Rdata"))
################################   END of FILE   ########################
