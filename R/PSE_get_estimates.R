#########################################################################
##       R PROGRAM: PSE_get_estimates.R
##
##         PROJECT: Evaluation of multiple-encounter population size
##                  estimators
##
##     DESCRIPTION: Extract and consolidate PSE estimates into a single
##                  data frame and save as both rds and csv.
##
##           INPUT: .Rdata files containing raw estimates
##
##          OUTPUT: PSE_estimates.rds
##                  PSE_estimates.csv
##
##      WRITTEN BY: Steve Gutreuter, CDC/CGH/DGHT Statistics, Estimation
##                                   and Modeling Team
##                  sgutreuter@cdc.gov
##
##            DATE: 2020-06-11
##
#########################################################################
library(tidyverse)
library(data.table)
basepath <- file.path(Sys.getenv("PROJ"), "PSE/PSEsim")
workpath <- file.path(basepath, "R")
datapath <- file.path(basepath, "data")
outpath <- file.path(basepath, "output")
setwd(workpath)
source(file.path(workpath, "PSE_sim_functions.R"))

est.objs <- c("BMA_M0_estimates.Rdata", "BMA_Mb_estimates.Rdata",
              "BMA_Mbh_estimates.Rdata", "BMA_Mbht_estimates.Rdata",
              "BMA_Mh_estimates.Rdata", "BMA_Mht_estimates.Rdata",
              "BMA_Mt_estimates.Rdata",
              "LCM_M0_estimates-1.Rdata", "LCM_M0_estimates-2.Rdata",
              "LCM_M0_estimates-3.Rdata", "LCM_M0_estimates-4.Rdata",
              "LCM_Mb_estimates.Rdata",
              "LCM_Mbh_estimates.Rdata", "LCM_Mbht_estimates.Rdata",
              "LCM_Mh_estimates.Rdata", "LCM_Mht_estimates.Rdata",
              "LCM_Mt_estimates.Rdata",
              "LLM_estimates-1.Rdata", "LLM_estimates-2.Rdata",
              "LLM_estimates-3.Rdata")

combineRdata <- function(Rdatanames, datapath){
    combined <- data.frame(NULL)
    counts <- data.frame(NULL)
    for(i in seq_along(Rdatanames)){
        efname <- file.path(datapath, Rdatanames[i])
        attach(efname, name = "currdata" )
        dfnames <- cgwtools::lsdata(efname)
        for(j in seq_along(dfnames)){
            dfj <- get(dfnames[j])
            if(substring(dfnames[j], 1, 3) == "llm"){
                dfj$ucl.truncated <- NULL
                dfj$HPDIcov <- dfj$HPDCIwidth <- dfj$HPDlcl <- dfj$HPDucl <-
                    rep(NA, dim(dfj)[1])
                dfj <- dfj[c("Events", "Rep","Ntrue", "Model", "Nest", "SE",
                             "Icoverage", "CIwidth", "lcl",  "ucl", "HPDIcov",
                             "HPDCIwidth", "HPDlcl", "HPDucl", "gen.modl", "gen.parms")]
            }
            combined <- rbind(combined, dfj)
            x <- data.frame(with(dfj, table(Events)))
            counts <- rbind(counts,
                            cbind(rep(dfnames[j], nrow(x)), x))
        }
        detach(currdata)
    }
    names(counts) <- c("Estimates frame", "Events", "Estimates")
    invisible(list(combined = combined, counts = counts))
}

#########################################################################
## Combine all estimates
#########################################################################
comb.list <- combineRdata(est.objs, datapath = datapath)
comb <- comb.list$combined
comb$Bias <- NULL

#########################################################################
## Recoding
#########################################################################
## Code pmean = mean probability of detectection per observation event
Idx1 <- comb$gen.parms == "p = 0.025" |
    comb$gen.parms == "betaparms=(1.32448,51.6548)" |
    comb$gen.parms == "p0=0.025, frac=0.5" |
    comb$gen.parms == "betaparms=(1.32448,51.6548), frac=0.5"
Idx2 <- comb$gen.parms == "p = 0.05" |
    comb$gen.parms == "betaparms=(1.26488,24.0327)" |
    comb$gen.parms == "p0=0.05, frac=0.5" |
    comb$gen.parms == "betaparms=(1.26488,24.0327), frac=0.5"
Idx3 <- comb$gen.parms == "p = 0.1" |
    comb$gen.parms == "betaparms=(1.14567,10.3111)" |
    comb$gen.parms == "p0=0.1, frac=0.5" |
    comb$gen.parms == "betaparms=(1.14567,10.3111), frac=0.5"
Idx4 <- comb$gen.parms == "p = 0.15" |
    comb$gen.parms == "betaparms=(1.02647,5.81667)" |
    comb$gen.parms == "p0=0.15, frac=0.5" |
    comb$gen.parms == "betaparms=(1.02647,5.81667), frac=0.5"
comb$pmean[Idx1] <- "0.025"
comb$pmean[Idx2] <- "0.050"
comb$pmean[Idx3] <- "0.100"
comb$pmean[Idx4] <- "0.150"
## Drop "Poisson 2" from Model values
comb$Model <- str_split(comb$Model, pattern = " ", simplify = TRUE)[, 1]
## Correct Model name Mth to Mht
comb$Model[comb$Model == "Mth"] <- "Mht"            ## Fix Model name Mht
comb$Model[comb$Model == "LCMCR"] <- "LCM"
comb$ModelType <- comb$Model
comb$ModelType[!(comb$Model == "BMA" | comb$Model == "LCM")] <- "LLM"
## Code an indicator of estimation failure (Nest <= 0)
comb$Ifail <- ifelse((comb$Nest <= 0 | is.na(comb$Nest)), 1, 0)

#########################################################################
## Frequencies by data-generating model
#########################################################################
table(comb$gen.modl)

#########################################################################
## Frequencies by data-generating model and generating parameters
#########################################################################
table(comb$gen.modl, comb$gen.parms)

#########################################################################
## Save the data
#########################################################################
comb$Model <- as.character(comb$Model)
comb$gen.modl <- as.character(comb$gen.modl)
comb$gen.parms <- as.character(comb$gen.parms)
saveRDS(comb, file = file.path(datapath, "PSE_estimates.rds"))
fwrite(comb, file = file.path(datapath, "PSE_estimates.csv"))
fwrite(comb.list$counts, file = file.path(outpath, "Estimate_counts.csv"))
################################   END of FILE   ########################
