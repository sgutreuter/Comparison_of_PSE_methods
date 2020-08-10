#################################################################################
##     R CODE FILE: plotBetaDist.R
##
##         PROJECT: Evaluation of multiple-encounter population size
##                  estimators
##
##     DESCRIPTION: Plot and summarize the beta distributions used for data
##                  generation.
##
##           INPUT: None
##
##          OUTPUT: Encounter_probability_beta_distributions-alt.pdf
##                  Encounter_probability_beta_distributions-alt.png
##
##      WRITTEN BY: Steve Gutreuter                  E-mail:  sgutreuter@cdc.gov
##                  Statistics, Estimation and Modeling Team
##                  Division of Global HIV & TB
##                  Center for Global Health
##                  Centers for Disease Control & Prevention
##
##            DATE: 2019-09-16
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
graphs <- file.path(outpath, "graphs")
setwd(workpath)
source(file.path(workpath, "PSE_sim_functions.R"))
library(lattice)

mean  <- 0.025
##var <- 0.002
CV <- 0.85
var <- (CV*mean)^2
cat("SD = ", var^0.5, "\n")
(parms1 <- BetaParms(mean, var))
BetaMoments(parms1[[1]], parms1[[2]])
a <- parms1[[1]]
b <- parms1[[2]]
x <- (1:1000)/1000
betadist1 <- data.frame(x = x, y = dbeta(x, a, b),
                       parms = paste("a = ", a, ", b = ", b, sep = ""))
mean  <-  0.05
var <- (CV*mean)^2
cat("SD = ", var^0.5, "\n")
(parms2 <- BetaParms(mean, var))
BetaMoments(parms2[[1]], parms2[[2]])
a <- parms2[[1]]
b <- parms2[[2]]
betadist2 <- data.frame(x = x, y = dbeta(x, a, b),
                       parms = paste("a = ", a, ", b = ", b, sep = ""))
mean  <-  0.10
var <- (CV*mean)^2
cat("SD = ", var^0.5, "\n")
(parms3 <- BetaParms(mean, var))
BetaMoments(parms3[[1]], parms3[[2]])
a <- parms3[[1]]
b <- parms3[[2]]
betadist3 <- data.frame(x = x, y = dbeta(x, a, b),
                       parms = paste("a = ", a, ", b = ", b, sep = ""))
mean <- 0.15
var <- (CV*mean)^2
cat("SD = ", var^0.5, "\n")
(parms4 <- BetaParms(mean, var))
BetaMoments(parms4[[1]], parms4[[2]])
a <- parms4[[1]]
b <- parms4[[2]]
betadist4 <- data.frame(x = x, y = dbeta(x, a, b),
                       parms = paste("a = ", a, ", b = ", b, sep = ""))
betadist <- rbind(betadist1, betadist2, betadist3, betadist4)
betadist$parms <- ordered(betadist$parms)

## Plot the distributions
betaplot <- xyplot(y ~ x, group = parms, data = betadist,
                   type = c("l", "l"), col = rep("black", 2), lty = c(1,5,2,4),
                   ylab = "Density", xlab = expression(italic(p)),
                   key = list(lines = list(lty = c(1,5,2,4)),
                              text = list(c(expression(paste("E("*italic(p)*") = 0.025; ",
                                                             "SD = 0.021",
                                                             sep = "")),
                                            expression(paste("E("*italic(p)*") = 0.050; ",
                                                             "SD = 0.042",
                                                             sep = "")),
                                            expression(paste("E("*italic(p)*") = 0.100; ",
                                                             "SD = 0.085",
                                                             sep = "")),
                                            expression(paste("E("*italic(p)*") = 0.150; ",
                                                             "SD = 0.127",
                                                             sep = "")))),
                              corner = c(0.95, 0.95)))
pdf(file.path(graphs, "Encounter_probability_beta_distributions-alt.pdf"),
    width = 5, height = 5)
plot(betaplot)
dev.off()

png(file.path(graphs, "Encounter_probability_beta_distributions-alt.png"),
    width = 5, height = 5, units = "in", res = 600)
plot(betaplot)
dev.off()


#################################  END of FILE  #################################
