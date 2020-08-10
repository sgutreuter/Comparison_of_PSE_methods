#########################################################################
##       R PROGRAM: PSE_estimate_LCM-3.R
##
##         PROJECT: Evaluation of multiple-encounter population size
##                  estimators
##
##     DESCRIPTION: Estimate population sizes from simulated samples
##                  using Bayesian nonparametric latent-class models as
##                  implemented in the LCMCR package
##
##                  Mb samples
##
##      WRITTEN BY: Steve Gutreuter, CDC/CGH/DGHT Statistics, Estimation
##                                   and Modeling Team
##                  sgutreuter@cdc.gov
##
##            DATE: 2019-09-16
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
library(Rcapture)
library(LCMCR)
library(HDInterval)
basepath <- file.path(Sys.getenv("PROJ"), "PSE/PSEsim")
workpath <- file.path(basepath, "R")
datapath <- file.path(basepath, "data")
setwd(workpath)
source(file.path(workpath, "PSE_sim_functions.R"))
dataname <- file.path(datapath, "PSEsimSamples_alt.Rdata")
attach(dataname)

#########################################################################
## loglinear model fitting to independent homogeneous capture
## probabilities
#########################################################################
seed_  <- 321
## startt <- Sys.time()
## test <- lcm.estimate(data = Mb.10K.025, events = 2:5,
##                      buffer_size = 1000000, thinning = 100,
##                      burnin = 500000, samples = 50000,
##                      nreps = 2, a_alpha = 0.025,
##                      b_alpha = 0.025, seed = seed_)
## endt <- Sys.time()
## print(list(start.time = startt, end.time = endt, elapsed = endt - startt))

## NOTE: Full run time will be approximately
t.start <- Sys.time()
lcm.Mb.10K.025 <- lcm.estimate(data = Mb.10K.025, events = 2:5,
                               buffer_size = 1000000, thinning = 100,
                               burnin = 500000, samples = 50000,
                               nreps = 400, a_alpha = 0.025,
                               b_alpha = 0.025, seed = seed_)
Sys.time() - t.start
Sys.time()

t.start <- Sys.time()
lcm.Mb.10K.050 <- lcm.estimate(data = Mb.10K.050, events = 2:5,
                               buffer_size = 1000000, thinning = 100,
                               burnin = 500000, samples = 50000,
                               nreps = 400, a_alpha = 0.025,
                               b_alpha = 0.025, seed = seed_)
Sys.time() - t.start
Sys.time()

t.start <- Sys.time()
lcm.Mb.10K.100 <- lcm.estimate(data = Mb.10K.100, events = 2:5,
                               buffer_size = 1000000, thinning = 100,
                               burnin = 500000, samples = 50000,
                               nreps = 400, a_alpha = 0.025,
                               b_alpha = 0.025, seed = seed_)
Sys.time() - t.start
Sys.time()

t.start <- Sys.time()
lcm.Mb.10K.150 <- lcm.estimate(data = Mb.10K.150, events = 2:5,
                               buffer_size = 1000000, thinning = 100,
                               burnin = 500000, samples = 50000,
                               nreps = 400, a_alpha = 0.025,
                               b_alpha = 0.025, seed = seed_)
Sys.time() - t.start
Sys.time()


#########################################################################
## Save the results
#########################################################################
save(lcm.Mb.10K.025, lcm.Mb.10K.050, lcm.Mb.10K.100, lcm.Mb.10K.150,
     file = file.path(datapath, "LCM_Mb_estimates.Rdata"))
################################   END of FILE   ########################
