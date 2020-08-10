#################################################################################
##     R CODE FILE: Table_3.R
##
##         PROJECT: Evaluation of multiple-encounter population size
##                  estimators
##
##     DESCRIPTION: Create a summary table of the RMSEs of the various
##                  heterogeneity corrections for models Mh and Mth
##
##           INPUT: Hetero_estimate_summaries_by_pmean.rds
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
##            DATE: 2019-08-27
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
## Create a table of percentages of  gen.modl correctly identified by the best-
## fitting log-linear model
#################################################################################
## Read the log-linear model estimates
dat <- readRDS(file.path(datapath, "PSE_estimates.rds"))
## Select only LLM results
dat <- dat[dat$ModelType == "LLM", ]
## Remove trend and Mbht because Rcapture does not estimate those
dat <- dat[!(dat$gen.modl == "Mbht" | dat$gen.modl == 'trend'), ]
dat <- dat[dat$Events > 2L, ]
## Code an indicator for identification of the correct data-generating model
dat$ImodID <- ifelse(dat$Model == dat$gen.modl, 1, 0)
## Compute proportions of data types which were correctly identified by the best
## fitting model.
tb0 <- dat %>%
    group_by(pmean, gen.modl, Events) %>%
    summarise(pCorrect = mean(ImodID, na.rm = TRUE))
## Reshape the resulting dataframe
tb0 <- data.frame(tb0)
tb0 <- melt(tb0, id.vars = c("pmean", "gen.modl", "Events"))
tb0$value <- round(100*tb0$value, 1)
tb0 <- dcast(tb0, pmean + Events ~ gen.modl)
## Generate rough LaTeX table using xtable for copy/paste to
(table3 <- toLatex(tb0))

#################################################################################
## Generate Appendix tables detailing frequencies of each type of
## (mis)identification.
#################################################################################
## Upper rows of pairs: frequencies
## Cross-tab dat$Model and dat$gen.modl by pmean and Events (one table for each
## value of pmean)
tb2 <- with(dat, table(pmean, gen.modl, Events, Model))
tb2f <- ftable(tb2)
tb2a <- dcast(as.data.frame(tb2f),
              as.formula(paste(paste(names(attr(tb2f, "row.vars")), collapse="+"),
                               "~", paste(names(attr(tb2f, "col.vars"))))))
tb2a <- mutate(tb2a, M0 = as.character(M0),  Mb = as.character(Mb),
               Mbh = as.character(Mbh), Mh = as.character(Mh),
               Mht = as.character(Mht), Mt = as.character(Mt))
tb2a$part  <- 1
## Lower rows of pairs: percentages
tb2p <- prop.table(tb2, margin = c(1:3))
tb2pf <- ftable(100*tb2p, row.vars = c("pmean", "gen.modl", "Events"))
tb2b <- dcast(as.data.frame(tb2pf),
              as.formula(paste(paste(names(attr(tb2pf, "row.vars")), collapse="+"),
                               "~", paste(names(attr(tb2pf, "col.vars"))))))[, 1:3]
tb2b <- cbind(tb2b,
              matrix(paste('(', sprintf(tb2pf, fmt = '%#.1f'), ')', sep = ''),
                     nrow = 72))
names(tb2b)[4:9] <- c("M0", "Mb", "Mbh", "Mh", "Mht", "Mt")
tb2b <- mutate(tb2b, M0 = as.character(M0),  Mb = as.character(Mb),
               Mbh = as.character(Mbh), Mh = as.character(Mh),
               Mht = as.character(Mht), Mt = as.character(Mt))
tb2b$part <- 2
## Interleave the rows of tb2a and tb2b producing tb2c
tb2c <- tb2a %>%
    mutate(row_number = row_number()) %>%
    bind_rows(tb2b %>% mutate(row_number = row_number())) %>%
    arrange(row_number, part)
tb2c$part <- NULL
tb2c$row_number <- NULL
View(tb2c)
names(tb2c)
## Write to LaTeX
(AppendixTable1 <- toLatex(tb2c))


#################################  END of FILE  #################################
