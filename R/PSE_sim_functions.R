#########################################################################
##       R PROGRAM: PSE_sim_functions.R
##
##         PROJECT: Evaluation of multiple-encounter population size
##                  estimators
##
##     DESCRIPTION: Sample-generator, estimation and helper functions
##
##      WRITTEN BY: Steve Gutreuter, CDC/CGH/DGHT Statistics, Estimation
##                                   and Modeling Team
##                  sgutreuter@cdc.gov
##
#########################################################################

## Function M0_gen
#' Generate encounter histories from model M0
#'
#' Generate K-source event encounter histories from homogeneous encounter
#' probabilities (Model M0), where all individuals share a common and
#' known encounter probability over all encounters.
#'
#' @param N Integer-valued true population size
#' @param reps  Number of replicate samples (integer)
#' @param p  Per-venue encounter probability (numeric)
#' @param K  Number of observation events (integer)
#' @param aggregate  (Logical) return aggregated event histories if TRUE
#' @param seed RNG seed (defaults to clock time)
#'
#' @return A list containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{sample}}{A dataframe containing the maximum sample.}
#' \item{\code{type}}{Envent-history type; "aggregated" or "individual".
#' \item{\code{Ntrue}}{True population size}
#' \item{\code{gen.modl}}{Name of the data-generating model.}
#' \item{\code{gen.parms}}{Parameters for the data-generating model.}
#' \item{\code{maxenc}}{The number of individual observed in the maximum sample.}
#' }
#' @export
M0_gen <- function(N = NULL, reps = NULL,  p = NULL, K = 2,
                   aggregate = TRUE, seed = NULL){
    if(!(K > 1)){stop("K must be > 1\n")}
    if(aggregate){type <- "aggregated"
                 } else {
                  type <- "individual"
                 }
    Ntrue <- N
    N <- rep(N, reps)
    e <- matrix(rep(numeric(), K), ncol = K)
    y <- data.frame(e, repl = numeric(), Freq = numeric())
    repl <- 0
    set.seed(seed)
        for(j in seq_along(N)){
            repl <- repl + 1
            e <- matrix(rbinom(K*N[j], size = 1, prob = p), ncol = K)
            this <- data.frame(e = e,
                               repl = rep(repl, N[j]),
                               Freq = rep(1, N[j]))
            if(aggregate) this <- sum_histories(this, eventcols = 1:K)
            y <- rbind(y, this)
        }
    y <- y[!rowSums(y[, (1:K)]) == 0, ]
    class(y) <- c("EventHistory", "data.frame")
    result <- list(call = match.call(), sample = y, type = type,
                   Ntrue = Ntrue, gen.modl = "M0",
                   gen.parms = paste("p =", p))
    class(result) <- c("PSEsimSample")
    invisible(result)
}

## Function Mh_gen
#' Generate encounter histories from model Mh
#'
#' Generate K-source event encounter histories from heterogeneous encounter
#' probabilities (Model Mh).  Among-individual heterogeneity is generated using
#' a beta distribution.
#'
#' @param N Integer-valued true population size
#' @param reps  Number of replicate samples (integer)
#' @param K  Number of observation events (integer)
#' @param aggregate  (Logical) return aggregated event histories if TRUE
#' @param betaparns Parameters of the beta distribution for hetergeneity
#' @param seed RNG seed (defaults to clock time)
#'
#' @return A list containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{sample}}{A dataframe containing the maximum sample.}
#' \item{\code{type}}{Envent-history type; "aggregated" or "individual".
#' \item{\code{Ntrue}}{True population size}
#' \item{\code{gen.modl}}{Name of the data-generating model.}
#' \item{\code{gen.parms}}{Parameters for the data-generating model.}
#' \item{\code{maxenc}}{The number of individual observed in the maximum sample.}
#' }
#' @export
Mh_gen <- function(N = NULL, reps = NULL,  K = 2,  betaparms = c(1,1),
                   aggregate = TRUE, seed = NULL){
    if(!(length(betaparms) == 2)){
        stop("betaparms must be a length 2 vector\n")}
    if(!(K > 1)){stop("K must be > 1\n")}
    if(aggregate){type <- "aggregated"
                 } else {
                  type <- "individual"
                 }
    Ntrue <- N
    N <- rep(N, reps)
    e <- matrix(rep(numeric(), K), ncol = K)
    y <- data.frame(e, repl = numeric(), Freq = numeric())
    repl <- 0
    set.seed(seed)
    for(i in seq_along(N)){
        repl <- repl + 1
        p <- rbeta(N[i], betaparms[1], betaparms[2])
        e <- matrix(rbinom(K*N[i], size = 1, prob = p),
                    byrow = FALSE, ncol = K)
        this <- data.frame(e = e,
                           repl = rep(repl, N[i]),
                           Freq = rep(1, N[i]))
        if(aggregate) this <- sum_histories(this, eventcols = 1:K)
        y <- rbind(y, this)
    }
    y <- y[!rowSums(y[, (1:K)]) == 0, ]
    class(y) <- c("EventHistory", "data.frame")
        result <- list(call = match.call(), sample = y, type = type,
                       Ntrue = Ntrue, gen.modl = "Mh",
                       gen.parms = paste("betaparms=(", betaparms[1],
                                         ",", betaparms[2], ")", sep =""))
    class(result) <- c("PSEsimSample")
    invisible(result)
}

## Function Mb_gen
#' Generate encounter histories from model Mb
#'
#' Generate K-source event encounter histories from the "behavioral" model
#' (Model Mb), where all individuals share a common initial encounter
#' probability, which is reduced by \code{frac} after the first encounter.
#'
#' @param N Integer-valued true population size
#' @param reps  Number of replicate samples (integer)
#' @param K  Number of observation events (integer)
#' @param p Probability of first encounter
#' @param frac Fraction reduction in p following the first encounter
#' @param aggregate  (Logical) return aggregated event histories if TRUE
#' @param seed RNG seed (defaults to clock time)
#'
#' @return A list containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{sample}}{A dataframe containing the maximum sample.}
#' \item{\code{type}}{Envent-history type; "aggregated" or "individual".
#' \item{\code{Ntrue}}{True population size}
#' \item{\code{gen.modl}}{Name of the data-generating model.}
#' \item{\code{gen.parms}}{Parameters for the data-generating model.}
#' \item{\code{maxenc}}{The number of individual observed in the maximum sample.}
#' }
#' @export
Mb_gen <- function(N = NULL, reps = NULL, K = 3, p = NULL,
                   frac = NULL, aggregate = TRUE, seed = NULL){
    if(!(p > 0 & p < 1)) stop("p not in (0,1)\n")
    if(!(frac > 0 & frac < 1)){stop("frac must be in (1,1)\n")}
    if(!(K > 2)){stop("K must be > 2; for K = 2 use Mt_gen\n")}
    if(aggregate){type <- "aggregated"
                 } else {
                  type <- "individual"
                 }
    Ntrue <- N
    N <- rep(N, reps)
    e0 <- matrix(rep(numeric(), K), ncol = K)
    y <- data.frame(e0, repl = numeric(), Freq = numeric())
    repl <- 0
    set.seed(seed)
    p2 <- p * frac
    for(i in seq_along(N)){
        e <- matrix(rep(NA, (N[i] * K)), ncol = K)
        repl <- repl + 1
        for(j in 1:N[i]){
            ienc <- TRUE
            for(k in 1:K){
                if(ienc) {
                    e[j,k] <- rbinom(1, size = 1, prob = p)
                    if(e[j,k] == 1) ienc <- FALSE
                } else {
                    e[j,k] <- rbinom(1, size = 1, prob = p2)
                }
            }
        }
        this <- data.frame(e = e,
                           repl = rep(repl, N[i]),
                           Freq = rep(1, N[i]))
        if(aggregate) this <- sum_histories(this, eventcols = 1:K)
        y <- rbind(y, this)
    }
    y <- y[!rowSums(y[, (1:K)]) == 0, ]
    class(y) <- c("EventHistory", "data.frame")
    result <- list(call = match.call(), sample = y, type = type,
                   Ntrue = Ntrue, gen.modl = "Mb",
                   gen.parms = paste("p0=", p, ", frac=", frac, sep=""))
    class(result) <- c("PSEsimSample")
    invisible(result)
}

## Function Mbh_gen
#' Generate encounter histories from model Mbh
#'
#' Generate K-source event encounter histories from heterogeneous encounter
#' probabilities with a behavioral effect (Model Mbh).  Among-individual
#' heterogeneity is generated using a beta distribution. After the first
#' encounter, an individuals probability of encounter is reduced by \code{frac}.
#'
#' @param N Integer-valued true population size
#' @param reps  Number of replicate samples (integer)
#' @param K  Number of observation events (integer)
#' @param betaparns Parameters of the beta distribution for hetergeneity
#' @param frac Fraction reduction in p following the first encounter
#' @param aggregate  (Logical) return aggregated event histories if TRUE
#' @param seed RNG seed (defaults to clock time)
#'
#' @return A list containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{sample}}{A dataframe containing the maximum sample.}
#' \item{\code{type}}{Envent-history type; "aggregated" or "individual".
#' \item{\code{Ntrue}}{True population size}
#' \item{\code{gen.modl}}{Name of the data-generating model.}
#' \item{\code{gen.parms}}{Parameters for the data-generating model.}
#' \item{\code{maxenc}}{The number of individual observed in the maximum sample.}
#' }
#' @export
Mbh_gen <- function(N = NULL, reps = NULL, K = 3, betaparms = c(1,1),
                    frac = NULL, aggregate = TRUE, seed = NULL){
    if(!(length(betaparms) == 2)){
        stop("betaparms must be a length 2 vector\n")
    }
    if(!(frac > 0 & frac < 1)){stop("frac must be in (1,1)\n")}
    if(!(K > 2)){stop("K must be > 2; for K = 2 use Mt_gen\n")}
    if(aggregate){type <- "aggregated"
                 } else {
                  type <- "individual"
                 }
    Ntrue <- N
    N <- rep(N, reps)
    e0 <- matrix(rep(numeric(), K), ncol = K)
    y <- data.frame(e0, repl = numeric(), Freq = numeric())
    repl <- 0
    set.seed(seed)
    for(i in seq_along(N)){
        e <- matrix(rep(NA, (N[i] * K)), ncol = K)
        repl <- repl + 1
        for(j in 1:N[i]){
            ienc <- TRUE
            p0 <- rbeta(1, betaparms[1], betaparms[2])
            p <- p0 * frac
            for(k in 1:K){
                if(ienc) {
                    e[j,k] <- rbinom(1, size = 1, prob = p0)
                    if(e[j,k] == 1) ienc <- FALSE
                } else {
                    e[j,k] <- rbinom(1, size = 1, prob = p)
                }
            }
        }
        this <- data.frame(e = e,
                           repl = rep(repl, N[i]),
                           Freq = rep(1, N[i]))
        if(aggregate) this <- sum_histories(this, eventcols = 1:K)
        y <- rbind(y, this)
    }
    y <- y[!rowSums(y[, (1:K)]) == 0, ]
    class(y) <- c("EventHistory", "data.frame")
    result <- list(call = match.call(), sample = y, type = type,
                   Ntrue = Ntrue, gen.modl = "Mbh",
                   gen.parms = paste("betaparms=(", betaparms[1],
                                     ",", betaparms[2],
                                     "), frac=", frac, sep =""))
    class(result) <- c("PSEsimSample")
    invisible(result)
}

## Function Mbht_gen
#' Generate encounter histories from model Mbht
#'
#' Generate K-source event encounter histories from heterogeneous encounter
#' probabilities with both behavioral and temporal effects (Model Mbht).
#' Initial encounter probabilities are independently beta-distributed over
#' individuals and encounter venues.  After the first encounter, an
#' individuals probabilities of encounter is reduced by \code{frac} for
#' any subsequent encounters.
#'
#' @param N Integer-valued true population size
#' @param reps  Number of replicate samples (integer)
#' @param K  Number of observation events (integer)
#' @param betaparms Parameters of the beta distribution for hetergeneity
#' @param frac Fraction reduction in p following the first encounter
#' @param aggregate  (Logical) return aggregated event histories if TRUE
#' @param seed RNG seed (defaults to clock time)
#'
#' @return A list containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{sample}}{A dataframe containing the maximum sample.}
#' \item{\code{type}}{Envent-history type; "aggregated" or "individual".
#' \item{\code{Ntrue}}{True population size}
#' \item{\code{gen.modl}}{Name of the data-generating model.}
#' \item{\code{gen.parms}}{Parameters for the data-generating model.}
#' \item{\code{maxenc}}{The number of individual observed in the maximum sample.}
#' }
#' @export
Mbht_gen <- function(N = NULL, reps = NULL, K = 3, betaparms = c(1,1),
                    frac = NULL, aggregate = TRUE, seed = NULL){
    if(!(length(betaparms) == 2)){
        stop("betaparms must be a length 2 vector\n")
    }
    if(!(frac > 0 & frac < 1)){stop("frac must be in (1,1)\n")}
    if(!(K > 2)){stop("K must be > 2; for K = 2 use Mt_gen\n")}
    if(aggregate){type <- "aggregated"
                 } else {
                  type <- "individual"
                 }
    Ntrue <- N
    N <- rep(N, reps)
    e0 <- matrix(rep(numeric(), K), ncol = K)
    y <- data.frame(e0, repl = numeric(), Freq = numeric())
    repl <- 0
    set.seed(seed)
    for(i in seq_along(N)){
        e <- matrix(rep(NA, (N[i] * K)), ncol = K)
        repl <- repl + 1
        for(j in 1:N[i]){
            ienc <- TRUE
            for(k in 1:K){
                p0 <- rbeta(1, betaparms[1], betaparms[2])
                p <- p0 * frac
                if(ienc) {
                    e[j,k] <- rbinom(1, size = 1, prob = p0)
                    if(e[j,k] == 1) ienc <- FALSE
                } else {
                    e[j,k] <- rbinom(1, size = 1, prob = p)
                }
            }
        }
        this <- data.frame(e = e,
                           repl = rep(repl, N[i]),
                           Freq = rep(1, N[i]))
        if(aggregate) this <- sum_histories(this, eventcols = 1:K)
        y <- rbind(y, this)
    }
    y <- y[!rowSums(y[, (1:K)]) == 0, ]
    class(y) <- c("EventHistory", "data.frame")
    result <- list(call = match.call(), sample = y, type = type,
                   Ntrue = Ntrue, gen.modl = "Mbht",
                   gen.parms = paste("betaparms=(", betaparms[1],
                                     ",", betaparms[2],
                                     "), frac=", frac, sep =""))
    class(result) <- c("PSEsimSample")
    invisible(result)
}

## Function Mt_gen
#########################################################################
#' Generate encounter histories from model Mt
#'
#' Generate K-source event encounter histories from heterogeneous encounter
#' probabilities (Model Mt).  All individuals share a common temporal sequence
#' of encounter probabilities which follow a beta distribution.
#'
#' @param N Integer-valued true population size
#' @param reps  Number of replicate samples (integer)
#' @param K  Number of observation events (integer)
#' @param betaparns Parameters of the beta distribution for temporal hetergeneity
#' @param aggregate  (Logical) return aggregated event histories if TRUE
#' @param seed RNG seed (defaults to clock time)
#'
#' @return A list containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{sample}}{A dataframe containing the maximum sample.}
#' \item{\code{type}}{Envent-history type; "aggregated" or "individual".
#' \item{\code{Ntrue}}{True population size}
#' \item{\code{gen.modl}}{Name of the data-generating model.}
#' \item{\code{gen.parms}}{Parameters for the data-generating model.}
#' \item{\code{maxenc}}{The number of individual observed in the maximum sample.}
#' }
#' @export
Mt_gen <- function(N = NULL, reps = NULL, K = 2, betaparms = c(1,1),
                       aggregate = TRUE, seed = NULL){
    if(!(length(betaparms) == 2)){
        stop("betaparms must be a length 2 vector\n")
    }
    if(!(K > 1)){stop("K must be > 1\n")}
    if(aggregate){type <- "aggregated"
                 } else {
                  type <- "individual"
                 }
    Ntrue <- N
    N <- rep(N, reps)
    e <- matrix(rep(numeric(), K), ncol = K)
    y <- data.frame(e, repl = numeric(), N = numeric())
    repl <- 0
    set.seed(seed)
    for(i in seq_along(N)){
        pt <- rbeta(K, betaparms[1], betaparms[2])
        repl <- repl + 1
        e <- matrix(rbinom(K*N[i], size = 1, prob = pt),
                    byrow = TRUE, ncol = K)
        pmat <- matrix(rep(pt, N[i]), byrow = TRUE, ncol = K)
        this <- data.frame(e = e,
                           repl = rep(repl, N[i]),
                           Freq = rep(1, N[i]))
        if(aggregate) this <- sum_histories(this, eventcols = 1:K)
        y <- rbind(y, this)
    }
    y <- y[!rowSums(y[, (1:K)]) == 0, ]
    class(y) <- c("EventHistory", "data.frame")
        result <- list(call = match.call(), sample = y, type = type,
                       Ntrue = Ntrue, gen.modl = "Mt",
                       gen.parms = paste("betaparms=(", betaparms[1],
                                         ",", betaparms[2], ")", sep =""))
    class(result) <- c("PSEsimSample")
    invisible(result)
}

## Function Mht_gen
#' Generate encounter histories from model Mht
#'
#' Generate K-source event encounter histories from independently
#' beta-distributed heterogeneous encounter probabilities among individuals
#' and times (Model Mht).
#'
#' @param N Integer-valued true population size
#' @param reps  Number of replicate samples (integer)
#' @param K  Number of observation events (integer)
#' @param betaparns Parameters of the beta distribution for hetergeneity
#' @param aggregate  (Logical) return aggregated event histories if TRUE
#' @param seed RNG seed (defaults to clock time)
#'
#' @return A list containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{sample}}{A dataframe containing the maximum sample.}
#' \item{\code{type}}{Envent-history type; "aggregated" or "individual".
#' \item{\code{Ntrue}}{True population size}
#' \item{\code{gen.modl}}{Name of the data-generating model.}
#' \item{\code{gen.parms}}{Parameters for the data-generating model.}
#' \item{\code{maxenc}}{The number of individual observed in the maximum sample.}
#' }
#' @export
Mht_gen <- function(N = NULL, reps = NULL, K = 2, betaparms = c(1,1),
                    aggregate = TRUE, seed = NULL){
    if(!(length(betaparms) == 2)){
        stop("betaparms must be a length 2 vector\n")
    }
    if(!(K > 1)){stop("K must be > 1\n")}
    if(aggregate){type <- "aggregated"
                 } else {
                  type <- "individual"
                 }
    Ntrue <- N
    N <- rep(N, reps)
    e <- matrix(rep(numeric(), K), ncol = K)
    y <- data.frame(e, repl = numeric(), N = numeric())
    repl <- 0
    set.seed(seed)
    for(i in seq_along(N)){
        repl <- repl + 1
        p <- rbeta(K*N[i], betaparms[1], betaparms[2])
        e <- matrix(rbinom(K*N[i], size = 1, prob = p), ncol = K)
        this <- data.frame(e = e,
                           repl = rep(repl, N[i]),
                           Freq = rep(1, N[i]))
        if(aggregate) this <- sum_histories(this, eventcols = 1:K)
        y <- rbind(y, this)
    }
    y <- y[!rowSums(y[, (1:K)]) == 0, ]
    class(y) <- c("EventHistory", "data.frame")
    result <- list(call = match.call(), sample = y, type = type,
                   Ntrue = Ntrue, gen.modl = "Mht",
                   gen.parms = paste("betaparms=(", betaparms[1],
                                     ",", betaparms[2], ")", sep =""))
    class(result) <- c("PSEsimSample")
    invisible(result)
}

## Function trend_gen
#' Generate encounter histories from a temporal trend model
#'
#' Generate 5-source event encounter histories for which all individuals share
#' a common trend in encounter probabilities.  This is a special case of model
#' Mt.
#'
#' @param N Integer-valued true population size
#' @param reps  Number of replicate samples (integer)
#' @param aggregate  (Logical) return aggregated event histories if TRUE
#' @param pvec A length-5 vector of the temporal sequence of encounter probabilties
#' @param seed RNG seed (defaults to clock time)
#'
#' @return A list containing the elements:
#' \describe{
#' \item{\code{Call}}{The function call.}
#' \item{\code{sample}}{A dataframe containing the maximum sample.}
#' \item{\code{type}}{Envent-history type; "aggregated" or "individual".
#' \item{\code{Ntrue}}{True population size}
#' \item{\code{gen.modl}}{Name of the data-generating model.}
#' \item{\code{gen.parms}}{Parameters for the data-generating model.}
#' \item{\code{maxenc}}{The number of individual observed in the maximum sample.}
#' }
#' @export
trend_gen <- function(N = NULL, reps = NULL, pvec = NULL, aggregate = TRUE,
                      seed = NULL){
    if(aggregate){type <- "aggregated"
                 } else {
                  type <- "individual"
                 }
    repl <- 0
    Ntrue <- N
    N <- rep(N, reps)
    e <- matrix(rep(numeric(), 5), ncol = 5)
    y <- data.frame(e = e, repl = numeric(0), Freq = numeric())
    set.seed(seed)
    for(i in seq_along(N)){
        repl <- repl + 1
        e <- matrix(rep(0, 5*N[i]), ncol = 5)
        for(j in 1:N[i]){
            for(t in 1:5){
                p <- pvec[t]
                e[j, t] <- rbinom(1, 1, p)
            }
        }
        this <- data.frame(e = e,
                           repl = rep(repl, N[i]),
                           Freq = rep(1, N[i]))
        if(aggregate) this <- sum_histories(this, eventcols = 1:5)
        y <- rbind(y, this)
    }
    y <- y[!rowSums(y[, (1:5)]) == 0, ]
    class(y) <- c("EventHistory", "data.frame")
        result <- list(call = match.call(), sample = y, type = type,
                       Ntrue = Ntrue, gen.modl = "trend",
                       gen.parms = paste("trendparms=(", pvec[1], ",", pvec[2],
                                         ",", pvec[3], ",", pvec[4], ",", pvec[5],
                                         ")", sep = ""))
    class(result) <- c("PSEsimSample")
    invisible(result)
}

## Function llm.estimate
#' Estimate population size using log-linear models
#'
#' Fit loglinear models the simulated event-history data using the Rcapture
#' package.  All feasible models are fit and results from the AIC-best model
#' are returned.
#'
#' @param data a PSEsimSample data object
#' @param events a vector of numbers of encounter events to include
#' @param nreps number of replicates from data object to estimate from
#'
#' @return A data frame containing the variables:
#' \describe{
#' \item{Events}{The number of encounter events used in estimation.}
#' \item{Rep}{The replicate number.}
#' \item{Ntrue}{True population size.}
#' \item{Model}{The name of the best-fitting loglinear model.}
#' \item(Nest){Estimate of population size.}
#' \item{SE}{Standard error of Nest}
#' \item{Icoverage}{Indicator for confidence-interval coverage}
#' \item{CIwidth}{Confidence-interval width}
#' \item{lcl}{Lower confidence limit}
#' \item{ucl}{Upper confidence limit}
#' \item{ucl.truncated}{Indicator for CL truncation at 50 times the estimate}
#' \item{gen.modl}{Name of the data-generating model}
#' \item{gen.parms}{Parameters of the data-generating model}
#' }
#' @export
llm.estimate <- function(data = NULL, events = seq(2, 5), nreps = NULL,
                         progress = FALSE){
    if(any(min(events) < 2)) stop("events not greater than 1\n")
    if(!"PSEsimSample" %in% class(data))
        stop("data not a PSEsimSample object\n")
    if(is.null(nreps)) stop("Invalid nreps\n")
    Ntrue <- data$Ntrue
    if(data$type == "individual"){
        data <- sum_histories.PSEsimSample(data, max(events))
    }
    nreps_ <- min(nreps, length(unique(data[["sample"]][["repl"]])))
    dat <- data$sample[data$sample$repl <= nreps_, ]
    y <- data.frame(Events = numeric(), Repl = numeric(),
                    Ntrue = numeric(), Model = character(),
                    Nest = numeric(), Bias = numeric(),
                    Icoverage = numeric(), CIwidth = numeric(),
                    lcl = numeric(), ucl = numeric(),
                    ucl.truncated = logical(),
                    gen.modl = character(),
                    gen.parms = character())
    for(i in events){
        dati <- sum_histories(dat, 1:i)
        for(j in 1:nreps_){
            if(progress & (j %% 20 == 0)){
                cat(paste("Estimating event = ", i, ", rep = ", j, "\r", sep = ""))}
            datf <- dati[dati$repl == j, 1:i]
            if(is.BadList(datf)) next
            datf$Freq <- dati[dati$repl == j, "Freq"]
            if(i > 2){
                fit <- data.frame(closedp.t(datf, dfreq = TRUE)$results)
            } else {
                fit <- data.frame(closedp.0(datf, dfreq = TRUE)$results)
            }
            if(dim(fit)[1] >= 2) {
                model <- row.names(fit[-grep("(LB)", row.names(fit)), ])
                model <- model[-grep("Darroch", model)]
                model <- model[-grep("Gamma", model)]
                fit <- fit[row.names(fit) %in% model, ]
            }
            minAIC <- min(fit$AIC)
            fit <- fit[fit$AIC == minAIC, ][1, ]
            est <- fit$abundance
            tags <- strsplit(row.names(fit), split = " ")
            mm <- as.character(tags[[1]][1])
            if(!(mm %in% c("Mb", "Mbh"))){
                if(length(tags[[1]]) > 1){
                    hh <- as.character(tags[[1]][2])
                    if(substr(hh, 1, 3) == "Poi") hh <- gsub("2", "", hh)
                    if(i > 2){
                        CI <- closedpCI.t(datf, dfreq = TRUE, m = mm,
                                          h = hh, fmaxSupCL = 50)$CI
                    } else {
                        CI <- closedpCI.0(datf, dfreq = TRUE, m = mm,
                                          h = hh, fmaxSupCL = 50)$CI
                    }
                } else {
                    if(i > 2){
                        CI <- closedpCI.t(datf, dfreq = TRUE, m = mm,
                                          fmaxSupCL = 50)$CI
                    } else {
                        CI <- closedpCI.0(datf, dfreq = TRUE, m = mm,
                                          fmaxSupCL = 50)$CI
                    }
                }
                lcl <- as.numeric(CI[, "infCL"])
                if(substr(CI[, "supCL"], 1, 1) == ">"){
                    ucl <- 50*est
                    ucl.truncated <- TRUE
                } else {
                    ucl <- as.numeric( CI[, "supCL"])
                    ucl.truncated <- FALSE
                }
            } else {
                lcl <- est - 1.96 * fit$stderr
                ucl <- est + 1.96 * fit$stderr
            }
            Icov <- as.integer(Ntrue >= lcl & Ntrue <= ucl)
            CIwidth <- ucl - lcl
            y <- rbind(y, data.frame(Events = i, Rep = j,
                                     Ntrue = data$Ntrue,
                                     Model = row.names(fit),
                                     Nest = est,
                                     SE = fit$stderr,
                                     Icoverage = Icov,
                                     CIwidth = CIwidth,
                                     lcl = lcl, ucl = ucl,
                                     ucl.truncated = ucl.truncated,
                                     gen.modl = as.character(data$gen.modl),
                                     gen.parms = as.character(data$gen.parms)))
        }
    }
    row.names(y) <- 1:dim(y)[1]
    invisible(y)
}

## Function lcm_estimate
#' Estimate population size using Baysian nonparametric latent-class models
#'
#' Fit nonparametric latent-class models to the simulated event-history data
#' using the LCMCR package.  Highest posterior density credible sets are
#' computed using the HDinterval package.
#'
#' @param data a PSEsimSample data object
#' @param events a vector of numbers of encounter events to include
#' @param nreps number of replicates from data object to estimate from
#' @param seed  RNG seed for MCMC simulation (defaults to clock time)
#' @param buffer  Memory buffer (see help(lcmCR))
#' @param thinning  MC chain thinning
#' @param burnin  MCMC burn-in iterations
#' @param samples  Retained MCMC sample size after thinning
#' @param nreps  Number of replicates from the data object to estimate
#' @param a_alpha  Stick-breaking prior parameter (see help(lcmCR))
#' @param b_alpha  Stick-breaking prior parameter (see help(lcmCR))
#' @param progress print progress to teminal (logical)
#'
#' @return A data frame containing the variables:
#' \describe{
#' \item{Events}{The number of encounter events used in estimation.}
#' \item{Rep}{The replicate number.}
#' \item{Ntrue}{True population size.}
#' \item{Model}{The name of the best-fitting loglinear model.}
#' \item(Nest){Estimate of population size.}
#' \item{SE}{Standard error of Nest}
#' \item{Icoverage}{Indicator for confidence-interval coverage}
#' \item{CIwidth}{Credible set width}
#' \item{lcl}{Lower limit of credible set}
#' \item{ucl}{Upper limit of credible set}
#' \item{HPDIcov}{Indicator for highest posterior density (HPD) interval
#' coverage}
#' \item{HPDCIwidth}{HPD interval width}
#' \item{HPDlcl}{Lower limit of HPD credible set}
#' \item{HPDucl}{Upper limit of HPD credible set}
#' \item{ucl.truncated}{Indicator for CL truncation at 50 times the estimate}
#' \item{gen.modl}{Name of the data-generating model}
#' \item{gen.parms}{Parameters of the data-generating model}
#' }
#' @export
lcm.estimate <- function(data = NULL, events = seq(2, 5), seed = NULL,
                         buffer_size = 100000, thinning = 100,
                         burnin = 1000, samples = 1000, nreps = NULL,
                         a_alpha = 0.025, b_alpha = 0.025,
                         progress = FALSE){
    if(any(min(events) < 2)) stop("events not greater than 1\n")
    if(!"PSEsimSample" %in% class(data))
        stop("data not a PSEsimSample object\n")
    if(is.null(nreps)) stop("Invalid nreps\n")
    Ntrue <- data$Ntrue
    if(data$type == "individual"){
        data <- sum_histories.PSEsimSample(data, max(events))
    }
    nreps_ <- min(nreps, length(unique(data[["sample"]][["repl"]])))
    dat <- data$sample[data$sample$repl <= nreps_, ]
    factolog <- function(x){as.logical(as.integer(as.character(x)))}
    y <- data.frame(Events = numeric(), Repl = numeric(),
                    Ntrue = numeric(), Model = character(),
                    Nest = numeric(), Bias = numeric(),
                    Icoverage = numeric(), CIwidth = numeric(),
                    lcl = numeric(), ucl = numeric(),
                    gen.modl = character(),
                    gen.parms = character())
    for(i in events){
        dati <- sum_histories(dat, 1:i)
        dati[ , seq(1, i)] <- lapply(dati[ , seq(1, i)], factor)
        class(dati) <- "data.frame"
        for(j in 1:nreps_){
            if(progress & (j %% 20 == 0)){
                cat(paste("Estimating event = ", i, ", rep = ", j, "\r", sep = ""))}
            datf <- dati[dati$repl == j, 1:i]
            if(is.BadList(datf)) next
            datf$Freq <- dati[dati$repl == j, "Freq"]
            smplr <- lcmCR(datf, tabular = TRUE, K = 10,
                           a_alpha = a_alpha, b_alpha = b_alpha,
                           seed = seed, buffer_size = buffer_size,
                           thinning = thinning)
            post <- lcmCR_PostSampl(smplr, burnin = burnin,
                                    samples = samples,
                                    thinning = thinning,
                                    output = FALSE)
            est <- mean(post, na.rm = TRUE)
            SE <- sd(post, na.rm = TRUE)
            CS <- quantile(post, probs = c(0.025, 0.975))
            Icov <- as.integer(Ntrue >= CS[1] & Ntrue <= CS[2])
            CIwidth <- CS[2] - CS[1]
            HPDCS <- hdi(post)
            HPDIcov <- as.integer(Ntrue >= HPDCS[1] & Ntrue <= HPDCS[2])
            HPDCIwidth <- HPDCS[2] - HPDCS[1]
            y <- rbind(y, data.frame(Events = i, Rep = j,
                                     Ntrue = Ntrue,
                                     Model = "LCMCR",
                                     Nest = est,
                                     SE = SE,
                                     Icoverage = Icov,
                                     CIwidth = CIwidth,
                                     lcl = CS[1],
                                     ucl = CS[2],
                                     HPDIcov = HPDIcov,
                                     HPDCIwidth = HPDCIwidth,
                                     HPDlcl = HPDCS[1],
                                     HPDucl = HPDCS[2],
                                     gen.modl = as.character(data$gen.modl),
                                     gen.parms = as.character(data$gen.parms)))
        }
    }
    row.names(y) <- 1:dim(y)[1]
    invisible(y)
}

## Function bma_estimate
#' Estimate population size using Bayesian model averaging
#'
#' Estimate population size by Bayesian model averaging of log-linear model
#' estimates as implemented in the dga package. A non-informative prior equal
#' to -log(1:Nmax) is used.  The Dirichlet  hyperpior parameter is 2^-K, where
#' K is the number of observation events.
#'
#' @param data a PSEsimSample data object
#' @param events a vector of numbers of encounter events to include
#' @param priorNmax the maximum number of unobserved individuals
#' @param Nreps the number of samples from the data objec to estimate
#' @param progress print progress to teminal (logical)
#'
#' @return A data frame containing the variables:
#' \describe{
#' \item{Events}{The number of encounter events used in estimation.}
#' \item{Rep}{The replicate number.}
#' \item{Ntrue}{True population size.}
#' \item{Model}{The name of the best-fitting loglinear model.}
#' \item(Nest){Estimate of population size.}
#' \item{SE}{Standard error of Nest}
#' \item{Icoverage}{Indicator for confidence-interval coverage}
#' \item{CIwidth}{Credible set width}
#' \item{lcl}{Lower limit of credible set}
#' \item{ucl}{Upper limit of credible set}
#' \item{HPDIcov}{Indicator for highest posterior density (HPD) interval
#' coverage}
#' \item{HPDCIwidth}{HPD interval width}
#' \item{HPDlcl}{Lower limit of HPD credible set}
#' \item{HPDucl}{Upper limit of HPD credible set}
#' \item{ucl.truncated}{Indicator for CL truncation at 50 times the estimate}
#' \item{gen.modl}{Name of the data-generating model}
#' \item{gen.parms}{Parameters of the data-generating model}
#' }
#' @export
bma.estimate <- function(data = NULL, events = seq(3, 5),
                         priorNmax = NULL,
                         nreps = NULL, progress = FALSE){
    if(any(min(events) < 3)) stop("events not greater than 2\n")
    if(!"PSEsimSample" %in% class(data))
        stop("data not a PSEsimSample object\n")
    if(is.null(priorNmax)) stop("priorNmax must be specified\n")
    if(is.null(nreps)) stop("Invalid nreps\n")
    Ntrue <- data$Ntrue
    if(data$type == "individual"){
        data <- sum_histories.PSEsimSample(data, max(events))
    }
    nreps_ <- min(nreps, length(unique(data[["sample"]][["repl"]])))
    dats <- data$sample[data$sample$repl <= nreps_, ]
    y <- data.frame(Events = numeric(), Repl = numeric(),
                    Ntrue = numeric(), Model = character(),
                    Nest = numeric(), Bias = numeric(),
                    Icoverage = numeric(), CIwidth = numeric(),
                    lcl = numeric(), ucl = numeric(),
                    gen.modl = character(),
                    gen.parms = character())
    for(i in events){
        dati <- sum_histories(dats, 1:i)
        class(dati) <- "data.frame"
        if(i == 3){
            graph <- graphs3
        } else {
            if(i == 4){
                graph <- graphs4
                } else {
                  graph <- graphs5
                }}
        for(j in 1:nreps_){
            if(progress & (j %% 20 == 0)){
                cat(paste("Estimating event = ", i, ", rep = ", j, "\r", sep = ""))}
            datf <- dati[dati$repl == j, 1:i]
            if(is.BadList(datf)) next
            datf$Freq <- dati[dati$repl == j, "Freq"]
            delta <- 2^(-i)
            priorNmiss <- 1:priorNmax
            datij  <- tidyr::uncount(datf[,-ncol(datf)], datf[[ncol(datf)]])
            oc <- make.strata(datij,
                              locations=rep("a", nrow(datij)))$overlap.counts
            N.obs <- sum(oc)
            oca <- array(oc, dim = rep(2, ncol(datf[, -ncol(datf)])))
            post.dens <- bma.cr(oca, delta = delta,
                                Nmissing = priorNmiss,
                                graphs = graph)
            sat.p <- sum(post.dens[nrow(post.dens), ])
            if(sat.p >= 0.15){
                post.dens <- post.dens[-nrow(post.dens), , drop = FALSE]}
            post.dens <- colSums(post.dens)
            post.dens <- post.dens / sum(post.dens)
            res <- list(x = priorNmiss, y = post.dens)
            class(res) <- "density"
            Ex <- sum(res$x * res$y)
            SE <- (sum(((res$x - Ex)^2) * res$y))^0.5
            Nest <- Ex + N.obs
            CS <- N.obs + quantile((res$x + res$y), probs = c(0.025, 0.975))
            Icov <- as.integer(Ntrue >= CS[1] & Ntrue <= CS[2])
            CIwidth <- CS[2] - CS[1]
            HPDCS <- N.obs + hdi(res)
            HPDIcov <- as.integer(Ntrue >= HPDCS[1] & Ntrue <= HPDCS[2])
            HPDCIwidth <- HPDCS[2] - HPDCS[1]
            y <- rbind(y, data.frame(Events = i, Rep = j,
                                     Ntrue = Ntrue,
                                     Model = "BMA",
                                     Nest = Nest,
                                     SE = SE,
                                     Icoverage = Icov,
                                     CIwidth = CIwidth,
                                     lcl = CS[1], ucl = CS[2],
                                     HPDIcov = HPDIcov,
                                     HPDCIwidth = HPDCIwidth,
                                     HPDlcl = HPDCS[1],
                                     HPDucl = HPDCS[2],
                                     gen.modl = as.character(data$gen.modl),
                                     gen.parms = as.character(data$gen.parms)))
        }
    }
    row.names(y) <- 1:dim(y)[1]
    invisible(y)
}

## Function sum_histories.PSEsimSample
#' Sum encounter histories over events
#'
#' Sum encounter histories over eventcols for a PSEsimSample object.  The
#' function is used to create aggregated event-history data from smaller
#' subsets of events in a PSEsimSample data object.
#'
#' @param x a PSEsimSample data object
#'   eventcols  Columns retained for summation
#'
#' @return A list of class "PSEsimSample" containing the elements:
#' \describe{
#' \item{\code{call}}{The function call}
#' \item{\code{sample}}{The summed sample}
#' \item{\code{type}}{Encounter-history type ("individual" or "aggregated")}
#' \item{\code{Ntrue}}{True population size}
#'}
sum_histories.PSEsimSample <- function(obj, eventcols){
    if(!"PSEsimSample" %in% class(obj)) stop("obj must have class PSEsimSample\n")
    x <- obj$sample
    keepcols <- dim(x)[2] - c(1, 0)
    events <- names(x[ , eventcols])
    x <- x[ , c(eventcols, keepcols)]
    x <- x[(!rowSums(x[ , eventcols]) == 0), ]
    event.hist <- NULL
    for(i in eventcols){
        event.hist <- paste(event.hist, as.character(x[[i]]), sep = "")
    }
    x$event.hist <- event.hist
    sums <- aggregate(Freq ~ repl + event.hist, data = x, FUN = sum)
    y <- subset(x, select = -Freq)
    y <- unique(y)
    y <- merge(y, sums)
    y$event.hist <- NULL
    y <- y[ , c(events, "repl", "Freq")]
    y <- y[order(y$repl), ]
    class(y) <- c("EventHistory", "data.frame")
    result <- list(call = obj$call, sample = y, type = "aggregated",
                   Ntrue = obj$Ntrue)
    class(result) <- "PSEsimSample"
    result
}

## Function sum_histories
#' Sum encounter histories over eventcols for a data frame or EventHistory object
#'
#' Sum encounter histories over eventcols for a data frame or EventHistory object.
#' The first columns are in {0,1} and the rows of those columns code an event
#' history.  The next-to-last column is the replicate number and an optional final
#' column is the frequency.
#'
#' @param x a dataframe or EventHistory object
#' @param eventcols  Columns retained for summation
#'
#' @return An EventHistory data frame
sum_histories <- function(x, eventcols){
    keepcols <- dim(x)[2] - c(1, 0)
    events <- names(x[ , eventcols])
    x <- x[ , c(eventcols, keepcols)]
    x <- x[(!rowSums(x[ , eventcols]) == 0), ]
    event.hist <- NULL
    for(i in eventcols){
        event.hist <- paste(event.hist, as.character(x[[i]]), sep = "")
    }
    x$event.hist <- event.hist
    sums <- aggregate(Freq ~ repl + event.hist, data = x, FUN = sum)
    y <- subset(x, select = -Freq)
    y <- unique(y)
    y <- merge(y, sums)
    y$event.hist <- NULL
    y <- y[ , c(events, "repl", "Freq")]
    y <- y[order(y$repl), ]
    class(y) <- c("EventHistory", "data.frame")
    y
}

## Function disagg_PSEsimSample
#' Disaggregate a PSEsimSample object having type = "aggregated"
#'
#' Disaggregate a PSEsimSample object having type = "aggregated" to produce
#' an equivalent PSEsimSaple object having type = "individual".  This is an
#' inverse function to sum_histories.PSEsimSample.
#'
#' @param obj a PSEsimSample data object of type = "aggregated"
#' @param eventcols columns retained for summation
#'
#' @return A PSEsimSample data object
disagg_PSEsimSample <- function(obj){
    if(!"PSEsimSample" %in% class(obj)) stop("obj must have class PSEsimSample\n")
    if(!obj$type == "aggregated") stop('obj must have type = "aggregated"\n')
    x <- obj$sample %>%
        tidyr::uncount(Freq)
    x$Freq <- 1
    x <- x[order(x$repl), ]
    class(x) <- c("EventHistory", "data.frame")
    result <- list(call = obj$call, sample = x, type = "individual",
                   Ntrue = obj$Ntrue)
    class(result) <- "PSEsimSample"
    invisible(result)
}

## Function plotBeta
#' Plot the Beta density having parameters a and b
#'
#'
#' @param a Beta distribution shape1 parameter
#' @param b Beta distribution shape2 parameter
#' @param N The number of points on the x axis
#'
#' @return Nothing.  This function produces a a plot as a side effect.
plotBeta <- function(a = 1, b = 1, N = 200){
    if(!(a > 0 & b > 0)) stop("a and b must be greater than 0")
    x <- (1:N)/N
    y <- dbeta(x, a, b)
    plot(x, y, type = "l",
         main = paste("a = ", a, "\n b = ", b, sep = ""))
}

## Function BetaMoments
#' Compute the mean and variance of the Beta distribution
#'
#' @param a Beta distribution shape1 parameter
#' @param b Beta distribution shape2 parameter
#'
#' @return A data frame containing the mean and variance
BetaMoments <- function(a = 1, b = 2){
    if(!(a > 0 & b > 0)) stop("a and b must be greater than 0")
    mean <- a / (a + b)
    var <- a * b / (((a + b)^2) * (a + b + 1))
    data.frame(mean = mean, variance = var)
}

## Function BetaParms
#' Compute the Beta distribution parameters from the mean and variance
#'
#' @param mu distribution mean
#' @param var distribution variance
#'
#' @return
BetaParms <- function(mu, var) {
    if(!((mu < 1) | (mu > 0))) stop("mu must be in (0, 1)")
    if(var >= mu * (1 - mu)) stop("var must be < mu * (1 - mu)")
    a <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    b <- a * (1 / mu - 1)
    return(list(shape1 = a, shape2 = b))
}

## Function getPlotData
#' Compute and extract data for package Rcapture descriptive plots
#'
#' The purpose is to isolate the neccessary data so that multiple
#' descriptive plots for multiple reps from multiple models can be
#' created using package lattice.
#' See: \code{Rcapture:::plot.descriptive}
#'
#' @param x a PSEsimSample object
#' @param reps number of replicates to process
#' @param dfreq see \code{Rcapture:::plot.descriptive}
#' @param model model name
#'
#' @return An object of class "descriptivePlotData"
getPlotData <- function(x, reps = NULL, dfreq = TRUE, model = NULL){

    if(is.null(reps)) stop("Specify number of replicates\n")
    if(!class(x) == "PSEsimSample") {
        stop("x must have class PSEsimSample\n")}
    x <- x[["sample"]]
    x <- x[x$repl <= reps, ]
    g1 <- data.frame(model = NULL, repl = NULL, captures = NULL,
                     y1 = NULL)
    g2 <- data.frame(model = NULL, repl = NULL, observation = NULL,
                     ui = NULL, logui = NULL)
    for(i in 1:reps){
        dsc <- x[x$repl == i, ]
        dsc <- Rcapture::descriptive(dsc[, -(ncol(dsc) - 1)],
                                     dfreq = dfreq)
        tinf <- if(is.null(dsc$call$t)){
                    FALSE
                } else {
                    is.infinite(dsc$call$t)
                }
        t <- nrow(dsc$base.freq)
        fi <- dsc$base.freq[, "fi"]
        g1.c3 <- if(tinf){
                     log(fi * factorial(1:t))
                 } else {
                     log(fi / choose(t, 1:t))
                 }
        g1i <- cbind(1:t, fi, g1.c3)
        g1i <- g1i[g1i[, 2] != 0, , drop = FALSE]
        g1i <- data.frame(model = model, i, g1i)
        g1i$model <- as.character(g1i$model)
        names(g1i) = c("model", "repl", "captures", "fi", "y1")
        g1 <- rbind(g1, g1i)
        if (dim(dsc$base.freq)[2] == 4) {
            ui <- dsc$base.freq[, "ui"]
            g2i <- cbind(1:t, ui, log(ui))
            g2i <- g2i[g2i[, 2] != 0, , drop = FALSE]
            g2i <- data.frame(model = rep(model, nrow(g2i)),
                              repl = i, g2i)
            g2i$model <- as.character(g2i$model)
            names(g2i) <- c("model", "repl", "observation", "ui",
                            "logui")
            g2 <- rbind(g2, g2i)
        }
    }
    result <- list(graph1 = g1, graph2 = g2)
    class(result) <- "descriptivePlotData"
    result
}

## Function is.BadList
#' Flag encounter-history matrices which cannot be used for estimation
#'
#' A prophylactic function to detect encounter history  matrices which cannot
#' be used for estimation for lack of overlap between the first two encounter
#' events.  It ignores lack of overlap among three or more events based on
#' the premise that those cases are unlikely.
#'
#' @param x  An event history matrix
#'
#' @return Logical; TRUE if the event history matrix is "bad"
is.BadList <- function(x){
    result <- FALSE
    x1 <- x[!((x[ , 1] == 0) & (x[ , 2] == 0)), ]
    if(dim(x1)[1] > 1){
        if(var(as.numeric(x1[ , 2])) == 0 & var(as.numeric(x1[ , 1])) == 0)
            result <- TRUE
    } else {
        result <- TRUE
        }
    result
}

## Function countEncounters
#' Count the total number of encounters in 2-5 source sampling
#'
#' Count the total number of encounters (realized sample size) in each of
#' 2-, ..., 5-source sampling from PSEsimSample data objects.
#'
#' @param x a character vector of names of 'PSEsimSample' objects
#'
#' @return A dataframe containing the elements:
#' \describe{
#' \item{\code{gen.modl}}{The data-generating model}
#' \item{\code{gen.parms}}{The parameters of the data-generating model}
#' \item{\code{Ntrue}}{The true population size}
#' \item{\code{repl}}{The sample replicate number
#' \item{\code{src2.count}}{The average sample size from 2-source sampling}
#' \item{\code{src3.count}}{The average sample size from 3-source sampling}
#' \item{\code{src4.count}}{The average sample size from 4-source sampling}
#' \item{\code{src5.count}}{The average sample size from 5-source sampling}
#' \item{\code{src2.frac}}{The average fraction of the population encountered during 2-source sampling}
#' \item{\code{src3.frac}}{The average fraction of the population encountered during 3-source sampling}
#' \item{\code{src4.frac}}{The average fraction of the population encountered during 4-source sampling}
#' \item{\code{src5.frac}}{The average fraction of the population encountered during 5-source sampling}
#' }
#' @export
countEncounters <- function(x){
    EnCounts <- data.frame(rep(character(), 2), rep(numeric(0), 9))
    for(i in seq_along(x)){
        datobj <- get(x[i])
        if(!class(datobj) == "PSEsimSample")
            stop("x does not name a 'PSEsimSample' object")
        samp <- datobj[["sample"]]
        nreps <- length(unique(samp$repl))
        cvec <- rep(NA, 4)
        for(k in 2:5){
            cvec[k - 1]  <- sum(samp[rowSums(samp[ , 1:k]) > 0, "Freq"])
            }
            encount <- rbind(data.frame(numeric(0)), cvec/nreps)
            names(encount) <- c("src2.count", "src3.count", "src4.count",
                                "src5.count")
            EnCounts <- rbind(EnCounts,
                      data.frame(cbind(datobj[["gen.modl"]],
                                       datobj[["gen.parms"]],
                                       datobj[["Ntrue"]],
                                       round(encount),
                                       round((encount / datobj$Ntrue) , 4))))
            }
    names(EnCounts) <- c("gen.modl", "gen.parms", "Ntrue", "src2.count",
                         "src3.count", "src4.count", "src5.count",
                         "src2.frac", "src3.frac", "src4.frac", "src5.frac")
    EnCounts
}

## Function PSE_sample_counts
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
    dfrm <- data.frame(p = rep(p, K),
                         k = 1:K)
    dfrm  <- dfrm %>%
        mutate(N_unique = N * p * (1 - p)^(k - 1),
               N_undetected = N * (1 - p)^k,
               Cum_unique = cumsum(N_unique),
               Pct_unique = 100 * N_unique / N,
               Pct_cum_unique = 100 * Cum_unique / N,
               Pct_undetected = 100 * N_undetected / N) %>%
        select(p, k, N_unique, Pct_unique, Cum_unique, Pct_cum_unique,
               N_undetected, Pct_undetected)
    dfrm
}

## Function llm.estimate.hetero
#' Estimate population size using models Mh, Mth, Mtbh
#'
#' Fit loglinear models the simulated event-history data using the Rcapture
#' package.  All feasible models Mh, Mth, Mtbh are fit in order to assess
#' effects of alternate heterogeniety formulations.
#'
#' @param data a PSEsimSample data object
#' @param events a vector of numbers of encounter events to include
#' @param nreps number of replicates from data object to estimate from
#'
#' @return A data frame containing the variables:
#' \describe{
#' \item{Events}{The number of encounter events used in estimation.}
#' \item{Rep}{The replicate number.}
#' \item{Ntrue}{True population size.}
#' \item{Model}{The name of the best-fitting loglinear model.}
#' \item{Hetero}{Type of heterogeneity correction.}
#' \item(Nest){Estimate of population size.}
#' \item{Icoverage}{Indicator for confidence-interval coverage}
#' \item{CIwidth}{Confidence-interval width}
#' \item{lcl}{Lower confidence limit}
#' \item{ucl}{Upper confidence limit}
#' \item{ucl.truncated}{Indicator for CL truncation at 50 times the estimate}
#' \item{gen.modl}{Name of the data-generating model}
#' \item{gen.parms}{Parameters of the data-generating model}
#' }
#' @export
llm.estimate.hetero <- function(data = NULL, events = seq(3, 5), nreps = NULL,
                         progress = FALSE){
    if(any(min(events) < 2)) stop("events not greater than 1\n")
    if(!"PSEsimSample" %in% class(data))
        stop("data not a PSEsimSample object\n")
    if(is.null(nreps)) stop("Invalid nreps\n")
    Ntrue <- data$Ntrue
    if(data$type == "individual"){
        data <- sum_histories.PSEsimSample(data, max(events))
    }
    nreps_ <- min(nreps, length(unique(data[["sample"]][["repl"]])))
    dat <- data$sample[data$sample$repl <= nreps_, ]
    y <- data.frame(Events = numeric(), Repl = numeric(),
                    Ntrue = numeric(), Model = character(),
                    Hetero = character(), Nest = numeric(),
                    lcl = numeric(), ucl = numeric(),
                    ucl.truncated = logical(), Icoverage = integer(),
                    CIwidth = numeric(),
                    gen.modl = character(),
                    gen.parms = character())
    for(i in events){
        dati <- sum_histories(dat, 1:i)
        for(j in 1:nreps_){
            if(progress & (j %% 20 == 0)){
                cat(paste("Estimating event = ", i, ", rep = ", j, "\r", sep = ""))
            }
            datf <- dati[dati$repl == j, 1:i]
            if(is.BadList(datf)) next
            datf$Freq <- dati[dati$repl == j, "Freq"]
            Mhfit.P <- closedpCI.t(datf, dfreq = TRUE, m = "Mh", h = "Poisson",
                                   fmaxSupCL = 50)$CI
            Mhfit.D <- closedpCI.t(datf, dfreq = TRUE, m = "Mh", h = "Darroch",
                                   fmaxSupCL = 50)$CI
            Mhfit.G <- closedpCI.t(datf, dfreq = TRUE, m = "Mh", h = "Gamma",
                                   fmaxSupCL = 50)$CI
            Mhtfit.P <- closedpCI.t(datf, dfreq = TRUE, m = "Mth", h = "Poisson",
                                   fmaxSupCL = 50)$CI
            Mhtfit.D <- closedpCI.t(datf, dfreq = TRUE, m = "Mth", h = "Darroch",
                                   fmaxSupCL = 50)$CI
            Mhtfit.G <- closedpCI.t(datf, dfreq = TRUE, m = "Mth", h = "Gamma",
                                   fmaxSupCL = 50)$CI
            comb <- data.frame(rbind(Mhfit.P, Mhfit.D, Mhfit.G, Mhtfit.P,
                                     Mhtfit.D, Mhtfit.G))
            tags <- strsplit(row.names(comb), split = " ")
            tags <- matrix(unlist(tags), nrow = 6, byrow = TRUE)
            comb$Model <- tags[ , 1]
            comb$Hetero <- tags[ , 2]
            comb$Events <- i
            comb$Repl = j
            comb$Ntrue = Ntrue
            comb$lcl <- as.numeric(comb[, "infCL"])
            idx <- substr(comb[, "supCL"], 1, 1) == ">"
            tryCatch(comb$ucl <- ifelse(idx, 50*comb$est, as.numeric(comb[, "supCL"])),
                     error = function(e) { print(c(i, j)) })
            comb$ucl.truncated <- ifelse(idx, TRUE, FALSE)

            comb$Icoverage <- as.integer(comb$Ntrue >= comb$lcl &
                                         comb$Ntrue <= comb$ucl)
            comb$CIwidth <- comb$ucl - comb$lcl
            comb <- comb[c("Events", "Repl", "Ntrue", "Model", "Hetero",
                           "abundance", "lcl", "ucl", "CIwidth", "Icoverage",
                           "ucl.truncated", "infoCI")]
            names(comb) <- c("Events", "Repl", "Ntrue", "Model", "Hetero",
                             "Nest", "lcl", "ucl", "CIwidth", "Icoverage",
                             "ucl.truncated", "infoCI")
            comb$gen.modl  <-  data$gen.modl
            comb$gen.parms <- data$gen.parms
            y <- rbind(y, comb)
        }
    }
    row.names(y) <- 1:dim(y)[1]
    invisible(y)
}
