#############################
# MCMCTree.R
# Jessie Pelosi 
# Last modified Nov 10, 2021
#############################

library(ggplot2)
library(MCMCtreeR)
library(phytools)
library(geiger)
library(dplyr)

## Filter loci to generate ultrametric tree 

MARE_LB <- read.csv("MARE_LB.csv")

quantile(MARE_LB$LB_score_Heterogeneity, probs = c(0.1,0.25,0.5,0.75))
#       10%      25%      50%      75% 
#  22.01931 24.18048 26.11491 28.31592

quantile(MARE_LB$MARE, probs = c(0.1,0.25,0.5,0.75))
#       10%      25%      50%      75% 
#   0.353846 0.397203 0.447552 0.499301 

# Only use loci with LB < 28.31 and MARE > 0.397, this leaves 364 loci 

## Prepare tree for input in MCMCTree

species_tree <- read.tree(file = "species_tree.tre")
plot(species_tree,edge.width=1, no.margin=TRUE)
nodelabels(cex=0.5)
monophyleticgroups <- tipDes(species_tree,c(15))
min_ages <- 98.7
max_ages <- 150

N <- estimateSkewNormal(minAge = 98.7,maxAge = 150, monoGroups = monophyleticgroups, phy = species_tree, shape = 50, estimateScale = F, addMode = 0.5, scale= 15, writeMCMCtree = TRUE)

plotMCMCtree(N$parameters[1,], method = "skewNormal", upperTime = max(max_ages), title = "Blah")

dev.off()

cauchy_results <- estimateCauchy(minAge = 93.7, maxAge = 156.39, monoGroups = monophyleticgroups, 
                                 phy = species_tree, plot = FALSE, offset = 0.25, estimateScale = F, scale = 0.05, writeMCMCtree = T)

plotMCMCtree(cauchy_results, method = "cauchy", upperTime = max(max_ages), title = "Blah")
