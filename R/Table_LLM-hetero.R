paste0("Start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", sep = "\n"))
################################################################################
##     R CODE FILE: Table_LLM-heter0.R
##
##         PROJECT: Evaluation of multiple-encounter population size
##                  estimators
##
##     DESCRIPTION: Create a summary table of the RMSEs of the various
##                  heterogeneity corrections for models Mh and Mth
##
##           INPUT: Hetero_estimate_summaries_by_pmean.rds
##
##          OUTPUT: LaTeX table in R transcript.
##
##
##      WRITTEN BY: Steve Gutreuter                  E-mail:  sgutreuter@cdc.gov
##                  Statistics, Estimation and Modeling Team
##                  Division of Global HIV & TB
##                  Center for Global Health
##                  Centers for Disease Control & Prevention
##
##            DATE: 2020-07-14
##
################################################################################

################################################################################
## Change basepath to the location of the files on your computer
################################################################################
basepath <- file.path(Sys.getenv("PROJ"), "PSE/PSEsim")

################################################################################
### Define additional file paths.  These should not be changed.
################################################################################
workpath <- file.path(basepath, "R")
datapath <- file.path(basepath, "data")
outpath <- file.path(basepath, "output")
graphs <- file.path(basepath, "output/graphs")
setwd(workpath)


library(tidyverse)
library(stringr)
library(reshape2)
library(memisc)  ## For toLatex

################################################################################
## Create a table of estimate counts and RMSE values
################################################################################
## Read the log-linear model estimates
dat <- readRDS(file.path(outpath, "Hetero_estimate_summaries_by_pmean.rds"))
dat <- ungroup(dat)
dat <- dat %>%
    dplyr::select(Model, pmean, Events, Hetero, n, RMSE)

#################################################################################
## Reshape the resulting dataframe and round digits.  The latter requires
## conversion to character.
#################################################################################
mltn <- melt(dat, id.vars = c("Model", "pmean", "Events", "Hetero"))
rcst <- dcast(mltn, Model + pmean + Events ~ Hetero + variable)
tb3 <- rcst %>%
    mutate(Darroch = round(ifelse(Darroch_RMSE > 1E9, Inf, Darroch_RMSE),
                          digits = 0),
           Gamma = round(ifelse(Gamma3.5_RMSE > 1E9, Inf, Gamma3.5_RMSE),
                         digits = 0),
           Poisson = round(ifelse(Poisson2_RMSE > 1E9, Inf, Poisson2_RMSE),
                           digits = 0),
           Darroch = as.character(Darroch),
           Gamma = as.character(Gamma),
           Poisson = as.character(Poisson),
           n_Poisson = as.character(Poisson2_n),
           n_Darroch = as.character(Darroch_n),
           n_Gamma = as.character(Gamma3.5_n)) %>%
    dplyr::select(Model, pmean, Events, n_Poisson, Poisson,
                  n_Darroch, Darroch, n_Gamma, Gamma)

#################################################################################
## Generate rough LaTeX table using toLatex for copy/paste to *.tex
#################################################################################
(table3 <- toLatex(tb3))

#################################  END of FILE  ################################
