#################################################################################
##     R CODE FILE: PSE_LCM_chains.R
##
##         PROJECT:
##
## INVESTIGATOR(S):
##
##     DESCRIPTION: Generate a sample of MCMC chains and examine evidence for
##                  convergence.
##
##    REFERENCE(S):
##
##           INPUT:
##
##          OUTPUT:
##
##      WRITTEN BY: Steve Gutreuter                  E-mail:  sgutreuter@cdc.gov
##                  Statistics, Estimation and Modeling Team
##                  Division of Global HIV & TB
##                  Center for Global Health
##                  Centers for Disease Control & Prevention
##
##            DATE: 201_-__-__
##
##      CHANGE LOG: Date        Change
#################################################################################

library(LCMCR)
library(tidyverse)
library(lattice)

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
source(file.path(workpath, "PSE_sim_functions.R"))
dataname <- file.path(datapath, "PSEsimSamples_alt.Rdata")
attach(dataname)


i <- 3               ## Use 3 lists
repno <- 1           ## Select replicate number (5 is fair)
##nreps <- 1
samples <- 500000
burnin <- 0

## Mbht -- No thinning
thinning <- 1
dat <- Mbht.10K.100$sample %>%
    dplyr::filter(.data$repl == repno)
dati <- sum_histories(dat, 1:i)
dati[ , seq(1, i)] <- lapply(dati[ , seq(1, i)], factor)
class(dati) <- "data.frame"
datf <- dati[dati$repl == repno, 1:i]
is.BadList(datf)
datf$Freq <- dati[dati$repl == repno, "Freq"]

chains <- data.frame(Repl = integer(), Iter = integer(), Nest = numeric())
set.seed(97531)        ## RNG control
for(j in 1:3){
    smplr <- lcmCR(datf, tabular = TRUE, K = 10,
                   a_alpha = 0.025, b_alpha = 0.025,
                   buffer_size = 1000000,
                   thinning = thinning)
    post <- lcmCR_PostSampl(smplr,
                            burnin = burnin,
                            samples = samples,
                            thinning = thinning,
                            output = FALSE)
    chains <- rbind(chains, data.frame(Repl = as.integer(rep(j, length(post))),
                                       Iter = as.integer(1:length(post)),
                                       Nest = post))

}
write.csv(chains, file.path(datapath, "Mbht_chains.csv"), row.names = FALSE)

#################################################################################
## Trace plots
#################################################################################
chains <- read.csv(file.path(datapath, "Mbht_chains.csv"))

oscipen <- options(scipen=10)

## One chain
Mbht_chainplot_1 <- xyplot(Nest ~ Iter,
                         data = chains[chains$Repl == 1, ],
                         type = "l",
                         ylab = "Population size",
                         xlab = "MCMC iteration",
                         col = c("blue"), alpha = 0.5,
                         par.settings = list(superpose.line = list(lwd=1)))
png(file.path(outpath, "graphs/Mbht_LCM_chainplot-1.png"),
    width = 5, height = 3, units = "in", res = 600)
plot(Mbht_chainplot_1)
dev.off()

## Three chains
Mbht_chainplot_3 <- xyplot(Nest ~ Iter, group = Repl,
                         data = chains,
                         type = "l",
                         ylab = "Population size",
                         xlab = "MCMC iteration",
                         col = c("blue", "red", "green"), alpha = 0.5,
                         par.settings = list(superpose.line = list(lwd=0.5)))
png(file.path(outpath, "graphs/Mbht_LCM_chainplot-3.png"),
    width = 5, height = 3, units = "in", res = 600)
plot(Mbht_chainplot_3)
dev.off()

## Three chains, expanded (last 100000 iterations)
Mbht_chainplot_3x <- xyplot(Nest ~ Iter, group = Repl,
                         data = chains[chains$Iter > 480000, ],
                         type = "l",
                         ylab = "Population size",
                         xlab = "MCMC iteration",
                         col = c("blue", "red", "green"), alpha = 0.5,
                         par.settings = list(superpose.line = list(lwd=0.5)))
png(file.path(outpath, "graphs/Mbht_LCM_chainplot-3x.png"),
    width = 5, height = 3, units = "in", res = 600)
plot(Mbht_chainplot_3x)
dev.off()

## Three chains, expanded (last 100000 iterations)
Mbht_chainplot_1x <- xyplot(Nest ~ Iter, group = Repl,
                         data = chains[chains$Iter > 480000 & chains$Repl == 2, ],
                         type = "l",
                         ylab = "Population size",
                         xlab = "MCMC iteration",
                         col = c("red"), alpha = 1,
                         par.settings = list(superpose.line = list(lwd=0.5)))
png(file.path(outpath, "graphs/Mbht_LCM_chainplot-1x.png"),
    width = 5, height = 3, units = "in", res = 600)
plot(Mbht_chainplot_1x)
dev.off()

options(oscipen)
#################################  END of FILE  #################################
