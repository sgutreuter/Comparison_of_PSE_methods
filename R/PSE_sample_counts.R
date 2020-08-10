################################################################################
##     R CODE FILE: PSE_sample_counts.R
##
##     DESCRIPTION: Compute numbers of unique individuals encountered up through
##                  list k as a function of encounter probability
##
##    REFERENCE(S):
##
##      WRITTEN BY: Steve Gutreuter                  E-mail:  sgutreuter@cdc.gov
##                  Statistics, Estimation and Modeling Team
##                  Division of Global HIV & TB
##                  Center for Global Health
##                  Centers for Disease Control & Prevention
##
##            DATE: 202_-__-__
################################################################################

library(tidyverse)
library(memisc)         ## For toLatex

################################################################################
## Change basepath to the location of the files on your computer
################################################################################
basepath <- file.path(Sys.getenv("PROJ"), "PSE/PSEsim")

################################################################################
### Define additional file paths.  These should not be changed.
################################################################################
datapath <- file.path(basepath, "data")
workpath <- file.path(basepath, "R")
output <- file.path(basepath, "output")
setwd(workpath)


## PSE_sample_counts
#' Compute sample counts in multiple-list (capture-recapture) sampling.
#'
#' Given encounter probability \emph{p}, number of lists (encounter events)
#' \emph{K} and population size \emph{N}, compute the expected number and
#' percentage of unique population members encountered on each list,
#' the cumulative number and cumulative percentage of unique individuals
#' encountered, and the number and percentage of the population detected
#' through list \emph{k}, \emph{k} = 1, ..., \emph{K}.
#'
#' @param N Integer-valued population size
#' @param p Real-valued encounter probability in (0,1)
#' @param K Integer-valued number of lists (encounters)
#'
#' @return A tibble or dataframe containing columns:
#' \describe{
#' \item{\code{p}}{List-wise encounter probability.}
#' \item{\code{k}}{Number of lists.}
#' \item{\code{N_unique}}{Number encountered for the first time.}
#' \item{\code{Cum_unique}}{Cumulative number encountered at least once.}
#' \item{\code{Pct_cum_unique}}{Cumulative percentage encountered at least once.}
#' \item{\code{N_undetected}}{Number remaining undetected.}
#' \item{\code{Pct_undetected}}{Percent remaining undetected.}
#' }
#'
#' @author Steve Gutreuter, \email{sgutreuter@@cdc.gov}
#'
PSE_sample_counts <- function(N, p, K){
    stopifnot(p > 0 & p < 1)
    counts <- data.frame(p = rep(p, K),
                         k = 1:K)
    counts  <- counts %>%
        mutate(N_unique = N * p * (1 - p)^(k - 1),
               N_undetected = N * (1 - p)^k,
               Cum_unique = cumsum(N_unique),
               Pct_unique = 100 * N_unique / N,
               Pct_cum_unique = 100 * Cum_unique / N,
               Pct_undetected = 100 * N_undetected / N) %>%
        dplyr::select(p, k, N_unique, Pct_unique, Cum_unique, Pct_cum_unique,
               N_undetected, Pct_undetected)
    counts
}

b1 <- PSE_sample_counts(10000, 0.025, 5)
b2 <- PSE_sample_counts(10000, 0.050, 5)
b3 <- PSE_sample_counts(10000, 0.100, 5)
b4 <- PSE_sample_counts(10000, 0.150, 5)
res <- rbind(b1, b2, b3, b4)
tbl <- res %>%
    dplyr::select(p, k, Pct_unique, Pct_cum_unique) %>%
    mutate(p = as.character(round(p, 3)),
           k = as.character(k),
           Pct_unique = as.character(round(Pct_unique, 2)),
           Pct_cum_unique = as.character(round(Pct_cum_unique, 2)))
(tbl <- toLatex(tbl))


#################################  END of FILE  ################################
