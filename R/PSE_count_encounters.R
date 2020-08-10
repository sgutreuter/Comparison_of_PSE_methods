#################################################################################
##       R PROGRAM: PSE_count_encounters.R
##
##         PROJECT: Evaluation of multiple-encounter population size estimators
##
##     DESCRIPTION: Count the average numbers of individuals encountered and
##                  their population fraction during each of 2-, 3-, 4- and 5-
##                  source sampling encounters for each combination of true
##                  population size, data-generating model and model parameters.
##
##           INPUT: PSEsimSamples_alt.Rdata
##
##          OUTPUT: PSE_sample_counts.rds
##                  PSE_sample_counts.csv
##
##      WRITTEN BY: Steve Gutreuter, CDC/CGH/DGHT Statistics, Estimation
##                                   and Modeling Team
##                  sgutreuter@cdc.gov
##
##            DATE: 2019-09-16
##
#################################################################################
library(tidyverse)
library(data.table)
basepath <- file.path(Sys.getenv("PROJ"), "PSE/PSEsim")
workpath <- file.path(basepath, "R")
datapath <- file.path(basepath, "data")
outpath <- file.path(basepath, "output")
setwd(workpath)
source(file.path(workpath, "PSE_sim_functions.R"))

#################################################################################
## Load the data object containing the simulated samples and provide a vector of
## PSEsimSample object names
################################################################################load(file.path(datapath, "PSEsimSamples_alt.Rdata"))
datnames <- c(  "M0.10K.025",   "M0.10K.050",   "M0.10K.100",   "M0.10K.150",
                "Mh.10K.025",   "Mh.10K.050",   "Mh.10K.100",   "Mh.10K.150",
                "Mb.10K.025",   "Mb.10K.050",   "Mb.10K.100",   "Mb.10K.150",
                "Mt.10K.025",   "Mt.10K.050",   "Mt.10K.100",   "Mt.10K.150",
               "Mbh.10K.025",  "Mbh.10K.050",  "Mbh.10K.100",  "Mbh.10K.150",
               "Mht.10K.025",  "Mht.10K.050",  "Mht.10K.100",  "Mht.10K.150",
              "Mbht.10K.025", "Mbht.10K.050", "Mbht.10K.100", "Mbht.10K.150",
              "tinc.10K.025", "tinc.10K.050", "tinc.10K.100", "tinc.10K.150",
              "tdec.10K.025", "tdec.10K.050", "tdec.10K.100", "tdec.10K.150")

wtf <- countEncounters(datnames)
saveRDS(wtf, file = file.path(outpath, "PSE_sample_counts.rds"))
fwrite(wtf, file = file.path(outpath, "PSE_sample_counts.csv"))

####################################   END of FILE   ############################
