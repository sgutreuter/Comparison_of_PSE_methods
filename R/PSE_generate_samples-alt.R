#################################################################################
##       R PROGRAM: PSE_generate_samples-alt.R
##
##         PROJECT: Evaluation of multiple-encounter population size
##                  estimators
##
## INVESTIGATOR(S):
##
##     DESCRIPTION: Compare estimates from Rcapture and LCMCR using
##                  simulated data.
##
##                  Part 1.  Generate samples
##
##   INPUT FILE(S): None
##
##  OUTPUT FILE(S): PSEsimSamples_alt.Rdata
##
##      WRITTEN BY: Steve Gutreuter, CDC/CGH/DGHT Statistics, Estimation
##                                   and Modeling Team
##                  sgutreuter@cdc.gov
##
##            NOTE: Effects of sample size are revealed by varying the encounter
##                  probabilities.  Model M0 is simulated with encounter
##                  probabilties of 0.025, 0.05, 0.10 and 0.15.  Inhomogeneities
##                  in the encounter probabilities are generated from the Beta
##                  distributions having means of 0.025, 0.50, 0.10 and 0.15,
##                  and standard deviation equal to 0.85 x mean (See
##                  plotBetaDist.R).  The data generators for trends in encounter
##                  probabilities use length-5 vectors of encounter probabilities
##                  centered at (event 3) 0.025, 0.05, 0.10 and 0.15.
##
##            DATE: 2019-09-16
##
#################################################################################

basepath <- file.path(Sys.getenv("PROJ"), "PSE/PSEsim")
workpath <- file.path(basepath, "R")
datapath <- file.path(basepath, "data")
setwd(workpath)
source(file.path(workpath, "PSE_sim_functions.R"))

seed0 <- 321
#################################################################################
## Generate random encounter histories from homogeneous (M0) populations
#################################################################################
M0.10K.025 <- M0_gen(10000, reps = 400, p = 0.025, K = 5, seed = seed0)
M0.10K.050 <- M0_gen(10000, reps = 400, p = 0.050, K = 5, seed = seed0)
M0.10K.100 <- M0_gen(10000, reps = 400, p = 0.100, K = 5, seed = seed0)
M0.10K.150 <- M0_gen(10000, reps = 400, p = 0.150, K = 5, seed = seed0)

#################################################################################
## Generate random encounter histories from Mh populations
#################################################################################
Mh.10K.025 <- Mh_gen(10000, reps = 400, K = 5, betaparms = c(1.32448, 51.6548),
                     seed = seed0)
Mh.10K.050 <- Mh_gen(10000, reps = 400, K = 5, betaparms = c(1.26488, 24.0327),
                       seed = seed0)
Mh.10K.100 <- Mh_gen(10000, reps = 400, K = 5, betaparms = c(1.14567, 10.3111),
                       seed = seed0)
Mh.10K.150 <- Mh_gen(10000, reps = 400, K = 5, betaparms = c(1.02647, 5.81667),
                       seed = seed0)

#################################################################################
## Generate random encounter histories from Mb populations where P(encounter)
## drops by a factor of 50% after first contact.
#################################################################################
Mb.10K.025  <- Mb_gen(10000, reps = 400, K = 5, p  = 0.025, frac = 0.5,
                        seed = seed0)
Mb.10K.050  <- Mb_gen(10000, reps = 400, K = 5, p  = 0.050, frac = 0.5,
                        seed = seed0)
Mb.10K.100  <- Mb_gen(10000, reps = 400, K = 5, p  = 0.100, frac = 0.5,
                        seed = seed0)
Mb.10K.150  <- Mb_gen(10000, reps = 400, K = 5, p  = 0.150, frac = 0.5,
                        seed = seed0)

#################################################################################
## Generate random encounter histories from Mbh populations where P(encounter)
## drops by 50% after first contact.
#################################################################################
Mbh.10K.025 <- Mbh_gen(10000, reps = 400, K = 5,
                       betaparms = c(1.32448, 51.6548), frac = 0.5,
                       seed = seed0)
Mbh.10K.050 <- Mbh_gen(10000, reps = 400, K = 5,
                       betaparms = c(1.26488, 24.0327),frac = 0.5,
                       seed = seed0)
Mbh.10K.100 <- Mbh_gen(10000, reps = 400, K = 5,
                       betaparms = c(1.14567, 10.3111), frac = 0.5,
                       seed = seed0)
Mbh.10K.150 <- Mbh_gen(10000, reps = 400, K = 5,
                       betaparms = c(1.02647, 5.81667), frac = 0.5,
                       seed = seed0)

#################################################################################
## Generate random encounter histories from Mbht populations where P(encounter)
## drops by 50% after first contact.
#################################################################################
Mbht.10K.025 <- Mbht_gen(10000, reps = 400, K = 5,
                       betaparms = c(1.32448, 51.6548), frac = 0.5,
                       seed = seed0)
Mbht.10K.050 <- Mbht_gen(10000, reps = 400, K = 5,
                       betaparms = c(1.26488, 24.0327),frac = 0.5,
                       seed = seed0)
Mbht.10K.100 <- Mbht_gen(10000, reps = 400, K = 5,
                       betaparms = c(1.14567, 10.3111), frac = 0.5,
                       seed = seed0)
Mbht.10K.150 <- Mbht_gen(10000, reps = 400, K = 5,
                       betaparms = c(1.02647, 5.81667), frac = 0.5,
                       seed = seed0)

#################################################################################
## Generate random encounter histories from Mt populations
#################################################################################
Mt.10K.025  <- Mt_gen(10000, reps = 400, K = 5, betaparms = c(1.32448, 51.6548),
                      seed = seed0)
Mt.10K.050  <- Mt_gen(10000, reps = 400, K = 5, betaparms = c(1.26488, 24.0327),
                      seed = seed0)
Mt.10K.100  <- Mt_gen(10000, reps = 400, K = 5, betaparms = c(1.14567, 10.3111),
                      seed = seed0)
Mt.10K.150  <- Mt_gen(10000, reps = 400, K = 5, betaparms = c(1.02647, 5.81667),
                   seed = seed0)

#################################################################################
## Generate random encounter histories from Mht populations
#################################################################################
Mht.10K.025  <- Mht_gen(10000, reps = 400, K = 5, betaparms = c(1.32448, 51.6548),
                        seed = seed0)
Mht.10K.050  <- Mht_gen(10000, reps = 400, K = 5, betaparms = c(1.26488, 24.0327),
                        seed = seed0)
Mht.10K.100  <- Mht_gen(10000, reps = 400, K = 5, betaparms = c(1.14567, 10.3111),
                        seed = seed0)
Mht.10K.150  <- Mht_gen(10000, reps = 400, K = 5, betaparms = c(1.02647, 5.81667),
                     seed = seed0)

#################################################################################
## Generate 5-encounter histories for increasing tend in p over encounters
#################################################################################
(probI <-  matrix(c(0.025, 0.050, 0.100, 0.150), nrow = 4) %*%
    matrix(c(0.25, 0.50, 1.0, 1.25, 1.5), nrow = 1))
tinc.10K.025 <- trend_gen(10000, reps = 400, pvec = as.vector(probI[1, ]),
                          seed = seed0)
tinc.10K.050 <- trend_gen(10000, reps = 400, pvec = as.vector(probI[2, ]),
                          seed = seed0)
tinc.10K.100 <- trend_gen(10000, reps = 400, pvec = as.vector(probI[3, ]),
                          seed = seed0)
tinc.10K.150 <- trend_gen(10000, reps = 400, pvec = as.vector(probI[4, ]),
                        seed = seed0)


#################################################################################
## Generate 5-encounter histories for decreasing tend in p over encounters
#################################################################################
(probD <- probI[, c(5:1)])
tdec.10K.025 <- trend_gen(10000, reps = 400, pvec = as.vector(probD[1, ]),
                          seed = seed0)
tdec.10K.050 <- trend_gen(10000, reps = 400, pvec = as.vector(probD[2, ]),
                          seed = seed0)
tdec.10K.100 <- trend_gen(10000, reps = 400, pvec = as.vector(probD[3, ]),
                          seed = seed0)
tdec.10K.150 <- trend_gen(10000, reps = 400, pvec = as.vector(probD[4, ]),
                          seed = seed0)

#################################################################################
## Save the generated samples
#################################################################################
save(  M0.10K.025,    M0.10K.050,    M0.10K.100,    M0.10K.150,
       Mh.10K.025,    Mh.10K.050,    Mh.10K.100,    Mh.10K.150,
       Mb.10K.025,    Mb.10K.050,    Mb.10K.100,    Mb.10K.150,
       Mt.10K.025,    Mt.10K.050,    Mt.10K.100,    Mt.10K.150,
      Mbh.10K.025,   Mbh.10K.050,   Mbh.10K.100,   Mbh.10K.150,
      Mht.10K.025,   Mht.10K.050,   Mht.10K.100,   Mht.10K.150,
     Mbht.10K.025,  Mbht.10K.050,  Mbht.10K.100,  Mbht.10K.150,
     tinc.10K.025,  tinc.10K.050,  tinc.10K.100,  tinc.10K.150,
     tdec.10K.025,  tdec.10K.050,  tdec.10K.100,  tdec.10K.150,
     file = file.path(datapath, "PSEsimSamples_alt.Rdata"))

################################   END of FILE   ################################
