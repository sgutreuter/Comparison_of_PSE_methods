#################################################################################
##     R CODE FILE: summarize_estimate_failure_frequencies.R
##
##         PROJECT: Evaluation of multiple-encounter population size
##                  estimators
##
##     DESCRIPTION: Tabulate frequencies of failed estimates
##
##           INPUT: PSE_estimates.rds
##
##          OUTPUT:
##
##
##      WRITTEN BY: Steve Gutreuter                  E-mail:  sgutreuter@cdc.gov
##                  Statistics, Estimation and Modeling Team
##                  Division of Global HIV & TB
##                  Center for Global Health
##                  Centers for Disease Control & Prevention
##
##            DATE: 2020-07-13
##
#################################################################################

#################################################################################
## Change basepath to the location of the files on your computer
#################################################################################
basepath <- file.path(Sys.getenv("PROJ"), "PSE/PSEsim")

#################################################################################
### Define additional file paths.  These should not be changed.
#################################################################################
workpath <- file.path(basepath, "R")
datapath <- file.path(basepath, "data")
outpath <- file.path(basepath, "output")
setwd(workpath)
library(tidyverse)
library(data.table)

dat <- readRDS(file.path(datapath, "PSE_estimates.rds"))
#################################################################################
## Assign all log-linear model variants to "Loglinear"
#################################################################################
smry <- dat %>%
    dplyr::select(Events, Rep, Model, Nest, SE, gen.modl, pmean, Ifail) %>%
    dplyr::filter(!is.na(Nest) & !is.na(pmean)) %>%
    dplyr::mutate(
               Model = ifelse(!(Model == "BMA" | Model == "LCM"), "EcM", Model),
               I_inf = ifelse(Nest > 1E9, 1, 0),
               I_neg = ifelse(Nest < 0, 1, 0),
               I_dgn = ifelse(SE < 0.1, 1, 0)) %>%
    group_by(Model, pmean, Events) %>%
    dplyr::summarize(
               Ntotal = n(),
               Ninf = sum(I_inf),
               Pinf = mean(I_inf),
               Nneg = sum(I_neg),
               Pneg = mean(I_neg),
               Ndgn = sum(I_dgn),
               Pdgn = mean(I_dgn))
data.frame(smry)

#################################  END of FILE  #################################
