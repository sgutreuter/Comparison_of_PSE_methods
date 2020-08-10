paste0("Start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", sep = "\n"))
################################################################################
##     R CODE FILE: computeSummaries.R
##
##         PROJECT: Evaluation of multiple-encounter population size
##                  estimators
##
##     DESCRIPTION: Compute summary measures of the performance of the estimators
##                  for each partition of the simulation.
##
##            NOTE: Estimates from the trend data-generating models are
##                  excluded because the values of pmean from those are not
##                  compatible with the other models.
##
##           INPUT: PSE_estimates.rds
##
##          OUTPUT: Estimate_summaries_by_generator.rds
##                  Estimate_summaries_by_generator.csv
##                  Estimate_summaries_sansM0trend_by_generator.rds
##                  Estimate_summaries_sansM0trend_by_generator.csv
##                  Estimate_summaries_by_pmean.rds
##                  Estimate_summaries_by_pmean.csv
##                  Estimate_summaries_sansM0trend_by_pmean.rds
##                  Estimate_summaries_sansM0trend_by_pmean.csv
##                  Estimate_summaries.rds
##                  Estimate_summaries.csv
##                  Estimate_summaries_sansM0trend.rds
##                  Estimate_summaries_sansM0trend.csv
##
##      WRITTEN BY: Steve Gutreuter                  E-mail:  sgutreuter@cdc.gov
##                  Statistics, Estimation and Modeling Team
##                  Division of Global HIV & TB
##                  Center for Global Health
##                  Centers for Disease Control & Prevention
##
##            DATE: 2020-06-11
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
setwd(workpath)
source(file.path(workpath, "PSE_sim_functions.R"))
library(tidyverse)
library(data.table)
library(e1071)

df <- readRDS(file.path(datapath, "PSE_estimates.rds"))
################################################################################
## Assign all log-linear model variants to "Loglinear"
################################################################################
df$Model <- as.character(df$Model)
df$Model[!(df$Model == "BMA" | df$Model == "LCM")] <- "Loglinear"
## Set global indicator for na.rm
narm <- TRUE

################################################################################
## Flag failed variance estimation in loglinear models for exclusion from
## summaries
################################################################################
df <- df %>%
    mutate(ISEfail = ifelse(Model == "Loglinear" & SE < 1, 1, 0))
with(df, table(Model, ISEfail))
################################################################################
## Compute summaries by Events, pmean, Model, gen.modl and gen.parms excluding
## estimates from the trend data-generating models which do not have compatible
## values of pmean.
################################################################################
sumry1 <- df %>%
    filter(!(gen.modl == "trend" | ISEfail == 1)) %>%
    group_by(Model, pmean, gen.modl, gen.parms, Events) %>%
    mutate(I_inf = ifelse(Nest > 1E9, 1, 0),
           I_neg = ifelse(Nest < 0, 1, 0)) %>%
    summarise(n_ = n(),
              RMSE = sqrt(mean((Nest - Ntrue)^2, na.rm = narm)),
              MAE = mean(abs(Nest - Ntrue), na.rm = narm),
              bias = mean(Nest, na.rm = narm) - mean(Ntrue),
              SE = sd(Nest, na.rm = narm),
              frac_gt_2True = mean((Nest > (2*Ntrue)), na.rm = narm),
              frac_lt_.5True = mean((Nest < (0.5*Ntrue)), na.rm = narm),
              Q90ratio = quantile(Nest, prob = 0.90, na.rm = narm)/mean(Ntrue),
              Coverage = mean(Icoverage, na.rm = narm),
              CoverageHPD = mean(HPDIcov, na.rm = narm),
              MedMOE = median(0.5 * CIwidth, na.rm = narm),
              MedMOEhpd = median(0.5 * (HPDucl - HPDlcl), na.rm = narm),
              skewness = skewness(Nest, na.rm = narm),
              Max = max(Nest, na.rm = narm),
              Ninf = sum(I_inf),
              Nneg = sum(I_neg),
              Ntrue = mean(Ntrue, na.rm = narm)
              )
names(sumry1) <- c("Model", "pmean", "gen.modl", "gen.parms", "Events",
                   "n", "RMSE", "MAE", "Bias", "SE", "frac_gt_2True",
                   "frac_lt_.5True", "Q90ratio", "Coverage", "CoverageHPD",
                   "MedMOE", "MedMOEhpd", "skewness", "Max", "n_inf", "n_neg",
                   "Ntrue")
sumry1$gen.modl <- as.character(sumry1$gen.modl)
sumry1$gen.parms <- as.character(sumry1$gen.parms)
sumry1$Model <- factor(sumry1$Model)
saveRDS(sumry1, file = file.path(outpath,
                                 "Estimate_summaries_by_generator.rds"))
fwrite(sumry1, file = file.path(outpath, "Estimate_summaries_by_generator.csv"))

################################################################################
## Compute summaries by Events, pmean, Model, gen.modl and gen.parms, excluding
## estimates from the MO and trend data-generating models
################################################################################
sumry2 <- df %>%
    filter(!(gen.modl == "M0" | gen.modl == "trend" | ISEfail == 1)) %>%
    group_by(Model, pmean, gen.modl, gen.parms, Events) %>%
    mutate(I_inf = ifelse(Nest > 1E9, 1, 0),
           I_neg = ifelse(Nest < 0, 1, 0)) %>%
    summarise(n_ = n(),
              RMSE = sqrt(mean((Nest - Ntrue)^2, na.rm = narm)),
              MAE = mean(abs(Nest - Ntrue), na.rm = narm),
              bias = mean(Nest, na.rm = narm) - mean(Ntrue),
              SE = sd(Nest, na.rm = narm),
              frac_gt_2True = mean((Nest > (2*Ntrue)), na.rm = narm),
              frac_lt_.5True = mean((Nest < (0.5*Ntrue)), na.rm = narm),
              Q90ratio = quantile(Nest, prob = 0.90, na.rm = narm)/mean(Ntrue),
              Coverage = mean(Icoverage, na.rm = narm),
              CoverageHPD = mean(HPDIcov, na.rm = narm),
              MedMOE = median(0.5 * CIwidth, na.rm = narm),
              MedMOEhpd = median(0.5 * (HPDucl - HPDlcl), na.rm = narm),
              skewness = skewness(Nest, na.rm = narm),
              Max = max(Nest, na.rm = narm),
              Ninf = sum(I_inf),
              Nneg = sum(I_neg),
              Ntrue = mean(Ntrue, na.rm = narm)
              )
names(sumry2) <- c("Model", "pmean", "gen.modl", "gen.parms", "Events",
                   "n", "RMSE", "MAE", "Bias", "SE", "frac_gt_2True",
                   "frac_lt_.5True", "Q90ratio", "Coverage", "CoverageHPD",
                   "MedMOE", "MedMOEhpd", "skewness", "Max", "n_inf", "n_neg",
                   "Ntrue")
sumry2$gen.modl <- as.character(sumry2$gen.modl)
sumry2$gen.parms <- as.character(sumry2$gen.parms)
sumry2$Model <- factor(sumry2$Model)
saveRDS(sumry2,
        file = file.path(outpath,
                         "Estimate_summaries_sansM0trend_by_generator.rds"))
fwrite(sumry2,
       file = file.path(outpath,
                        "Estimate_summaries_sansM0trend_by_generator.csv"))

################################################################################
## Compute summaries by Events, pmean and Model
################################################################################
sumry3 <- df %>%
    filter(!(gen.modl == "trend" | ISEfail == 1)) %>%
    group_by(Model, pmean, Events) %>%
    mutate(I_inf = ifelse(Nest > 1E9, 1, 0),
           I_neg = ifelse(Nest < 0, 1, 0)) %>%
    summarise(n_ = n(),
              RMSE = sqrt(mean((Nest - Ntrue)^2, na.rm = narm)),
              MAE = mean(abs(Nest - Ntrue), na.rm = narm),
              bias = mean(Nest, na.rm = narm) - mean(Ntrue),
              SE = sd(Nest, na.rm = narm),
              frac_gt_2True = mean((Nest > (2*Ntrue)), na.rm = narm),
              frac_lt_.5True = mean((Nest < (0.5*Ntrue)), na.rm = narm),
              Q90ratio = quantile(Nest, prob = 0.90, na.rm = narm)/mean(Ntrue),
              Coverage = mean(Icoverage, na.rm = narm),
              CoverageHPD = mean(HPDIcov, na.rm = narm),
              MedMOE = median(0.5 * CIwidth, na.rm = narm),
              MedMOEhpd = median(0.5 * (HPDucl - HPDlcl), na.rm = narm),
              skewness = skewness(Nest, na.rm = narm),
              Max = max(Nest, na.rm = narm),
              Ninf = sum(I_inf),
              Nneg = sum(I_neg),
              Ntrue = mean(Ntrue, na.rm = narm)
              )
names(sumry3) <- c("Model", "pmean", "Events", "n", "RMSE", "MAE", "Bias", "SE",
                   "frac_gt_2True", "frac_lt_.5True", "Q90ratio", "Coverage",
                   "CoverageHPD", "MedMOE", "MedMOEhpd", "skewness", "Max",
                   "n_inf", "n_neg", "Ntrue")
sumry3$Model <- factor(sumry3$Model)
saveRDS(sumry3, file = file.path(outpath, "Estimate_summaries_by_pmean.rds"))
fwrite(sumry3, file = file.path(outpath, "Estimate_summaries_by_pmean.csv"))

################################################################################
## Compute summaries by Events, pmean and Model, excluding estimates from the
## MO and trend data-generating models
################################################################################
sumry4 <- df %>%
    filter(!(gen.modl == "M0" | gen.modl == "trend" | ISEfail == 1)) %>%
    group_by(Model, pmean, Events) %>%
    mutate(I_inf = ifelse(Nest > 1E9, 1, 0),
           I_neg = ifelse(Nest < 0, 1, 0)) %>%
    summarise(n_ = n(),
              RMSE = sqrt(mean((Nest - Ntrue)^2, na.rm = narm)),
              MAE = mean(abs(Nest - Ntrue), na.rm = narm),
              bias = mean(Nest, na.rm = narm) - mean(Ntrue),
              SE = sd(Nest, na.rm = narm),
              frac_gt_2True = mean((Nest > (2*Ntrue)), na.rm = narm),
              frac_lt_.5True = mean((Nest < (0.5*Ntrue)), na.rm = narm),
              Q90ratio = quantile(Nest, prob = 0.90, na.rm = narm)/mean(Ntrue),
              Coverage = mean(Icoverage, na.rm = narm),
              CoverageHPD = mean(HPDIcov, na.rm = narm),
              MedMOE = median(0.5 * CIwidth, na.rm = narm),
              MedMOEhpd = median(0.5 * (HPDucl - HPDlcl), na.rm = narm),
              skewness = skewness(Nest, na.rm = narm),
              Max = max(Nest, na.rm = narm),
              Ninf = sum(I_inf),
              Nneg = sum(I_neg),
              Ntrue = mean(Ntrue, na.rm = narm)
              )
names(sumry4) <- c("Model", "pmean", "Events", "n", "RMSE", "MAE", "Bias", "SE",
                   "frac_gt_2True", "frac_lt_.5True", "Q90ratio", "Coverage",
                   "CoverageHPD", "MedMOE", "MedMOEhpd", "skewness", "Max",
                   "n_inf", "n_neg", "Ntrue")
sumry4$Model <- factor(sumry4$Model)
saveRDS(sumry4,
        file = file.path(outpath,
                         "Estimate_summaries_sansM0trend_by_pmean.rds"))
fwrite(sumry4,
       file = file.path(outpath, "Estimate_summaries_sansM0trendby_pmean.csv"))


################################################################################
## Compute summaries by Events and Model
################################################################################
sumry5 <- df %>%
    filter(!(gen.modl == "trend" | ISEfail == 1)) %>%
    group_by(Model, Events) %>%
    mutate(I_inf = ifelse(Nest > 1E9, 1, 0),
           I_neg = ifelse(Nest < 0, 1, 0)) %>%
    summarise(n_ = n(),
              RMSE = sqrt(mean((Nest - Ntrue)^2, na.rm = narm)),
              MAE = mean(abs(Nest - Ntrue), na.rm = narm),
              bias = mean(Nest, na.rm = narm) - mean(Ntrue),
              SE = sd(Nest, na.rm = narm),
              frac_gt_2True = mean((Nest > (2*Ntrue)), na.rm = narm),
              frac_lt_.5True = mean((Nest < (0.5*Ntrue)), na.rm = narm),
              Q90ratio = quantile(Nest, prob = 0.90, na.rm = narm)/mean(Ntrue),
              Coverage = mean(Icoverage, na.rm = narm),
              CoverageHPD = mean(HPDIcov, na.rm = narm),
              MedMOE = median(0.5 * CIwidth, na.rm = narm),
              MedMOEhpd = median(0.5 * (HPDucl - HPDlcl), na.rm = narm),
              skewness = skewness(Nest, na.rm = narm),
              Max = max(Nest, na.rm = narm),
              Ninf = sum(I_inf),
              Nneg = sum(I_neg),
              Ntrue = mean(Ntrue, na.rm = narm)
              )
names(sumry5) <- c("Model", "Events", "n", "RMSE", "MAE", "Bias", "SE",
                   "frac_gt_2True", "frac_lt_.5True", "Q90ratio", "Coverage",
                   "CoverageHPD", "MedMOE", "MedMOEhpd", "skewness", "Max",
                   "n_inf", "n_neg", "Ntrue")
sumry5$Model <- factor(sumry5$Model)
saveRDS(sumry5, file = file.path(outpath, "Estimate_summaries.rds"))
fwrite(sumry5, file = file.path(outpath, "Estimate_summaries.csv"))

################################################################################
## Compute summaries by Events and Model, excluding estimates from the
## MO and trend data-generating models
################################################################################
sumry6 <- df %>%
    filter(!(gen.modl == "M0" | gen.modl == "trend" | ISEfail == 1)) %>%
    group_by(Model, Events) %>%
    mutate(I_inf = ifelse(Nest > 1E9, 1, 0),
           I_neg = ifelse(Nest < 0, 1, 0)) %>%
    summarise(n_ = n(),
              RMSE = sqrt(mean((Nest - Ntrue)^2, na.rm = narm)),
              MAE = mean(abs(Nest - Ntrue), na.rm = narm),
              bias = mean(Nest, na.rm = narm) - mean(Ntrue),
              SE = sd(Nest, na.rm = narm),
              frac_gt_2True = mean((Nest > (2*Ntrue)), na.rm = narm),
              frac_lt_.5True = mean((Nest < (0.5*Ntrue)), na.rm = narm),
              Q90ratio = quantile(Nest, prob = 0.90, na.rm = narm)/mean(Ntrue),
              Coverage = mean(Icoverage, na.rm = narm),
              CoverageHPD = mean(HPDIcov, na.rm = narm),
              MedMOE = median(0.5 * CIwidth, na.rm = narm),
              MedMOEhpd = median(0.5 * (HPDucl - HPDlcl), na.rm = narm),
              skewness = skewness(Nest, na.rm = narm),
              Max = max(Nest, na.rm = narm),
              Ninf = sum(I_inf),
              Nneg = sum(I_neg),
              Ntrue = mean(Ntrue, na.rm = narm)
              )
names(sumry6) <- c("Model", "Events", "n", "RMSE", "MAE", "Bias", "SE",
                   "frac_gt_2True", "frac_lt_.5True", "Q90ratio", "Coverage",
                   "CoverageHPD", "MedMOE", "MedMOEhpd", "skewness", "Max",
                   "n_inf", "n_neg", "Ntrue")
sumry6$Model <- factor(sumry6$Model)
saveRDS(sumry6, file = file.path(outpath, "Estimate_summaries_sansM0trend.rds"))
fwrite(sumry6, file = file.path(outpath, "Estimate_summaries_sansM0trend.csv"))

#################################  END of FILE  ################################
