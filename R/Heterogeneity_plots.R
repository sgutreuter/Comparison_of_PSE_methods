################################################################################
##       R PROGRAM: Heterogeneity_plots
##
##         PROJECT:
##
## INVESTIGATOR(S):
##
##     DESCRIPTION: Construct heterogeneity plots for PSE data.
##
##       REFERENCE: Baillargeon S, Rivest L-P. Rcapture: Loglinear models for
##                    capture-recapture in R.  Journal of Statistical
##                    software 2007; 19(5):1-31.
##
##      WRITTEN BY: Steve Gutreuter, CDC/CGH/DGHT/HIDMS/SEM
##                  E-mail:  sgutreuter@cdc.gov
##
##            DATE: 2020-05-11
##
################################################################################
library(Rcapture)
library(ggplot2)
library(ggthemes)
library(cowplot)

basepath <- file.path(Sys.getenv("PROJ"), "PSE/PSEsim")
workpath <- file.path(basepath, "R")
datapath <- file.path(basepath, "data")
outpath <- file.path(basepath, "output")
graphs <- file.path(outpath, "graphs")
setwd(workpath)
source(file.path(workpath, "PSE_sim_functions.R"))

################################################################################
## Get the simulated data
################################################################################
dataname <- file.path(datapath, "PSEsimSamples_alt.Rdata")
attach(dataname)

################################################################################
## Extract the simulated data data from which to construct the graphs
################################################################################
addEp <- function(x, Ep){
    x$graph1$Ep <- Ep
    x$graph2$Ep <- Ep
    x
}
nreps <- 15
dpM0_025 <- getPlotData(M0.10K.025, reps = nreps, model = "M0")
dpM0_025 <- addEp(dpM0_025, "0.025")
dpM0_150 <- getPlotData(M0.10K.150, reps = nreps, model = "M0")
dpM0_150 <- addEp(dpM0_150, "0.150")
dpMh_025 <- getPlotData(Mh.10K.025, reps = nreps, model = "Mh")
dpMh_025 <- addEp(dpMh_025, "0.025")
dpMh_150 <- getPlotData(Mh.10K.150, reps = nreps, model = "Mh")
dpMh_150 <- addEp(dpMh_150, "0.150")
dpMb_025 <- getPlotData(Mb.10K.025, reps = nreps, model = "Mb")
dpMb_025 <- addEp(dpMb_025, "0.025")
dpMb_150 <- getPlotData(Mb.10K.150, reps = nreps, model = "Mb")
dpMb_150 <- addEp(dpMb_150, "0.150")
dpMt_025 <- getPlotData(Mt.10K.025, reps = nreps, model = "Mt")
dpMt_025 <- addEp(dpMt_025, "0.025")
dpMt_150 <- getPlotData(Mt.10K.150, reps = nreps, model = "Mt")
dpMt_150 <- addEp(dpMt_150, "0.150")
dpMtbh_025 <- getPlotData(Mbht.10K.025, reps = nreps, model = "Mbht")
dpMtbh_025 <- addEp(dpMtbh_025, "0.025")
dpMtbh_150 <- getPlotData(Mbht.10K.150, reps = nreps, model = "Mbht")
dpMtbh_150 <- addEp(dpMtbh_150, "0.150")

graph1 <- rbind(dpM0_025$graph1, dpM0_150$graph1,
                dpMh_025$graph1, dpMh_150$graph1,
                dpMb_025$graph1, dpMb_150$graph1,
                dpMt_025$graph1, dpMt_150$graph1,
                dpMtbh_025$graph1, dpMtbh_150$graph1)
graph2 <- rbind(dpM0_025$graph2, dpM0_150$graph2,
                dpMh_025$graph2, dpMh_150$graph2,
                dpMb_025$graph2, dpMb_150$graph2,
                dpMt_025$graph2, dpMt_150$graph2,
                dpMtbh_025$graph2, dpMtbh_150$graph2)
graph1$model <- ordered(graph1$model, levels = c("M0", "Mb", "Mh", "Mt", "Mbht"))
graph1$Ep <- factor(graph1$Ep)
graph2$model <- ordered(graph2$model, levels = c("M0", "Mb", "Mh", "Mt", "Mbht"))
graph2$Ep <- factor(graph2$Ep)

################################################################################
## Construct the graphs and align the top and bottom elements using
## cowplot::plotgrid
################################################################################
## Top half: The fi plot
lwd_ <- 0.25
ylbl <- expression("log" * bgroup("(", frac(italic(f)[italic(i)],
                                            bgroup("(",
                                            atop(italic(t), italic(i)), ")")),
                                  ")"))
g1 <- ggplot(graph1, aes(captures, y1, group = repl)) +
    geom_line(size = lwd_) + ylab(ylbl) +
    facet_grid(Ep ~ model) +
    coord_fixed(ratio = 0.5) +
    theme_base()
## Botttom half: The ui plot
g2 <- ggplot(graph2, aes(observation, logui, group = repl)) +
    geom_line(size = lwd_) + ylab(expression('log('*italic(u)[italic(i)]*')')) +
    facet_grid(Ep ~ model) +
    coord_fixed(ratio = 0.5) +
    theme_base()
## Create the aligned plot
png(file.path(graphs, "Heterogeneity_plots.png"), width = 6, height = 6,
    units = "in", res = 600)
plot_grid(g1, g2, align = 'hv', nrow = 2)
dev.off()
pdf(file.path(graphs, "Heterogeneity_plots.pdf"))
plot_grid(g1, g2, align = 'hv', nrow = 2)
dev.off()

## Alternate formatting of y-axis labels:
## Top half: The fi plot
ylblx <- expression("log" * bgroup("(", italic(f)[italic(i)]*"/"*bgroup("(",
                                            atop(scriptstyle(italic(t)), scriptstyle(italic(i))), ")"),
                                  ")"))
g1x <- ggplot(graph1, aes(captures, y1, group = repl)) +
    geom_line(size = lwd_) + ylab(ylblx) +
    facet_grid(Ep ~ model) +
    coord_fixed(ratio = 0.5) +
    theme_base()
## Create the aligned plot
png(file.path(graphs, "Heterogeneity_plots-alt.png"), width = 6, height = 6,
    units = "in", res = 600)
plot_grid(g1x, g2, align = 'hv', nrow = 2)
dev.off()
pdf(file.path(graphs, "Heterogeneity_plots-alt.pdf"))
plot_grid(g1x, g2, align = 'hv', nrow = 2)
dev.off()
