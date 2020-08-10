paste0("Start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", sep = "\n"))
#################################################################################
##     R CODE FILE: Table_2.R
##
##         PROJECT: Evaluation of multiple-encounter population size
##                  estimators
##
##     DESCRIPTION: Create Table 2 containing a comparative summary of the
##                  performances of the estimators
##
##
##           INPUT: Estimate_summaries_by_generator.rds
##
##          OUTPUT: LaTeX-ready table content
##
##
##      WRITTEN BY: Steve Gutreuter                  E-mail:  sgutreuter@cdc.gov
##                  Statistics, Estimation and Modeling Team
##                  Division of Global HIV & TB
##                  Center for Global Health
##                  Centers for Disease Control & Prevention
##
##            DATE: 2020-06-11
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
graphs <- file.path(basepath, "output/graphs")
setwd(workpath)

library(tidyverse)
library(stringr)
library(reshape2)
library(memisc)  ## For toLatex

#################################################################################
## Extract estimates---excluding those from data-generating models M0 and the two
## linear trend models---and arrange, filter, select and transform variables.
#################################################################################
dat <- readRDS(file.path(outpath, "Estimate_summaries_sansM0trend_by_pmean.rds"))
dat <- ungroup(dat)
dat$Model <- as.character(dat$Model)
dat$Model[dat$Model == "Loglinear"]  <- "a_EcM"
dat$Model[dat$Model == "BMA"]  <- "b_BMA"
dat$Model[dat$Model == "LCM"]  <- "c_LCM"
dat <- dat %>%
    arrange(Events, Model, pmean) %>%
    dplyr::select(Model, pmean, Events, RMSE, Bias, frac_gt_2True,
                  frac_lt_.5True, Coverage, CoverageHPD, n) %>%
    dplyr::mutate(RMSE = ifelse(RMSE > 1E9, Inf, RMSE),
                  Bias = ifelse(Bias > 1E9, Inf, Bias))

#################################################################################
## Reshape the resulting dataframe and round digits.  The latter requires
## conversion to character.
#################################################################################
mltn <- melt(dat, id.vars = c("Events", "Model", "pmean"))
rcst <- dcast(mltn, variable + pmean ~ Events + Model, mean)
tb2top <- cbind(rcst[1:8, 1:2], round(rcst[1:8, 3:13], digits = 0))
tb2top <- tb2top %>%
    mutate_all(as.character)
tb2mid <- cbind(rcst[9:24, 1:2], round(rcst[9:24, 3:13], digits = 3))
tb2mid <- tb2mid%>%
    mutate_all(as.character)
tb2bot <- cbind(rcst[25:28, 1:2], round(rcst[25:28, 3:13], digits = 0))
tb2bot <- tb2bot%>%
    mutate_all(as.character)
tb2 <- bind_rows(tb2top, tb2mid, tb2bot)

#################################################################################
## Generate rough LaTeX table using toLatex for copy/paste to *.tex
#################################################################################
(table2 <- toLatex(tb2))

#################################  END of FILE  #################################
