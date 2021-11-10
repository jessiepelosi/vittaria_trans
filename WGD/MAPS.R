##################################
# MAPS.R
# Jessie Pelosi 
# Last modified Nov 10, 2021
##################################

library(ggplot2)
library(MCMCtreeR)
library(phytools)
library(geiger)
library(WGDgc)
library(dplyr)

## Prepare MAPS input ##

ultrametric_tree <- read.tree("../../mcmctree_aspl_out/cauchy_1/Cauchy_output_nwk.tre")
plot(ultrametric_tree)

taxa_to_keep <- read.csv("tips_to_keep_MAPS.csv", row.names = 2, header = T)
name_list <- name.check(ultrametric_tree, taxa_to_keep)
checked_names <- name_list$tree_not_data
names_to_drop <- as.vector(checked_names)
trimmed_tree <- drop.tip(ultrametric_tree, names_to_drop)
plot(trimmed_tree)
write.tree(phy = trimmed_tree, file = "trimmed_tree_for_MAPS_sim.tre")

gf_count <- read.delim("vittaria_trans/count_genefamilies_wagner.txt")

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
gm_mean(gf_count$X12)

gfs_present_at_root <- gf_count %>% 
  filter(X12 >= 1)

# Null simulations 
sim_tree <- read.simmap("vittaria_trans/trimmed_tree_for_MAPS_sim.tre", vers = 1.1)
plot(sim_tree)

genecounts_for_WGDgc <- dplyr::select(gfs_present_at_root, 'NDUV', 'SKYV', 'WCLG', 
                                      'DCDT', 'WQML','FLTD', 'KJZG')

subset1 <- sample_n(genecounts_for_WGDgc, 1000)
MLEGeneCount(sim_tree, geneCountData = subset1, geomMean = 1.186565, 
             conditioning = "none", fixedRetentionRates = T)
             
# Birth = 0.002513682; Death = 0.005460978 (may or may not be present at root)
# Birth = 0.002266459; Death = 0.002397053 (present at root)

subset2 <- sample_n(genecounts_for_WGDgc, 1000)
MLEGeneCount(sim_tree, geneCountData = subset2, geomMean = 1.186565, 
             conditioning = "none", fixedRetentionRates = T)
# Birth = 0.002417878; Death = 0.00551798 (may or may not be present at root)
# Birth = 0.002313318; Death = 0.002496395 (present at root)

subset3 <- sample_n(genecounts_for_WGDgc, 1000)
MLEGeneCount(sim_tree, geneCountData = subset3, geomMean = 1.186565, 
             conditioning = "none", fixedRetentionRates = T)
# Birth = 0.002558946; Death = 0.005221092 (may or may not be present at root)
# Birth = 0.002315056; Death = 0.002458696 (present at root)

# Average Birth = 0.002496835; Death = 0.005400017 may or may not be present at root
# Avearge Birth = 0.002298278; Death = 0.002450715 present at root 

# Positive simulations 

possim_tree <- read.simmap("vittaria_trans/trimmed_tree_for_MAPS_POSSIM.tre", vers = 1.1)

# Lambda and mu estimates for positive simulations with a single WGD at the MRCA
# of V. appalachiana and V. lineata 

subset1 <- sample_n(genecounts_for_WGDgc, 1000)
MLEGeneCount(possim_tree, geneCountData = subset1, geomMean = 1.186565, 
             conditioning = "none", fixedRetentionRates=TRUE, startingQ=0.2)
# Avearge Birth = 0.00216097; Death = 0.002604341 present at root 

subset2 <- sample_n(genecounts_for_WGDgc, 1000)
MLEGeneCount(possim_tree, geneCountData = subset2, geomMean = 1.186565, 
             conditioning = "none", fixedRetentionRates=TRUE, startingQ=0.2)
# Avearge Birth = 0.001991621; Death = 0.002581447 present at root 

subset3 <- sample_n(genecounts_for_WGDgc, 1000)
MLEGeneCount(possim_tree, geneCountData = subset3, geomMean = 1.186565, 
             conditioning = "none", fixedRetentionRates=TRUE, startingQ=0.2)
# Avearge Birth = 0.002033379; Death = 0.002769139 present at root 

mean(0.00216097, 0.001991621, 0.002033379)
mean(0.002604341, 0.002581447, 0.002769139)

## MAPS Results ##

MAPS_subtrees <- read.csv("vittaria_trans/pteridaceae_subtree.csv")
nullsim <- read.csv("vittaria_3000_nullsimulations.csv")
possim <- read.csv("Vittaria_3000_possimulations.csv")

ggplot() + 
  geom_line(data = pos, mapping = aes(x = Node, y = Duplicated.Per, group = Tree), color = "black", alpha = 0.5) +
  geom_line(data = sim, mapping = aes(x=Node, y = Value, group  = Tree), color = "grey", alpha = 0.5) +
  geom_point(data = MAPS_subtrees, mapping = aes(x = MRCA, y = Duplication.Per), color = "dodgerblue") +
  geom_line(data = MAPS_subtrees, mapping = aes(x = MRCA, y = Duplication.Per), size =1, color = "dodgerblue") +
  xlab("Node") + ylab("Proportion of Subtrees") + ylim(0,0.4) + theme_classic()

ggsave("MAPS_SimFig.png", dpi = 300, height = 5, width = 7)
