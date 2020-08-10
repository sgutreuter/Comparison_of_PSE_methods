#################################################################################
##     R CODE FILE: plotEstimates.R
##
##         PROJECT: Evaluation of multiple-encounter population size
##                  estimators
##
##     DESCRIPTION: Contruct violin plots showing the distributions of the
##                  PSE estimates.
##
##            NOTE: Estimates from M0 and the trend data-generating models are
##                  omitted
##
##           INPUT: Estimates_summaries.rds
##
##          OUTPUT: PSE_violin_plots.png
##                  PSE_violin_plots.pdf
##
##      WRITTEN BY: Steve Gutreuter                  E-mail:  sgutreuter@cdc.gov
##                  Statistics, Estimation and Modeling Team
##                  Division of Global HIV & TB
##                  Center for Global Health
##                  Centers for Disease Control & Prevention
##
##            DATE: 2020-06-12
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
library(ggplot2)
library(ggthemes)
library(scales)


#################################################################################
## Extract the estimates
#################################################################################
df1 <- readRDS(file.path(datapath, "PSE_estimates.rds"))
df1 <- df1 %>%
    filter(!is.na(pmean) & Nest < 1E9 & Nest > 0 & !gen.modl == "trend") %>%
    mutate(ModelType = factor(ModelType, ordered = TRUE,
                              levels = c("LLM", "BMA", "LCM")),
           Events = factor(Events, ordered = TRUE, levels = c(2, 3, 4, 5)),
           pmean = factor(pmean, ordered = TRUE,
                          levels = c("0.025", "0.050", "0.100", "0.150")))

#################################################################################
## Construct violin plots of the estimates
#################################################################################
vp <- ggplot(df1, aes(Events, Nest)) +
    ylab(expression(hat(italic(N)))) + xlab("Encounter events") +
    facet_grid(pmean ~ ModelType) +
    geom_hline(yintercept = 10000, size = 0.5, col = "gray") +
    geom_violin() +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    theme_base()
png(file.path(graphs, "PSE_violin_plots.png"), width = 5, height = 5,
    units = "in", res = 600)
plot(vp)
dev.off()
pdf(file.path(graphs, "PSE_violin_plots.pdf"), width = 5, height = 5)
plot(vp)
dev.off()

#################################  END of FILE  #################################
