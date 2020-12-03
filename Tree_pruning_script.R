# Prune time-calibrated tree from Testo and Sundue 2016 (https://www.sciencedirect.com/science/article/pii/S1055790316302287#s0100)
# Written by W.L. Testo, modified by J.A. Pelosi
# Last modified 3 Decemeber 2020 

library("geiger")

#load tree
tree <- read.tree("4ktree.tre")
#tree <- ladderize(tree,right=F)

plot(tree, cex=0.3)

#read in taxon list
keep <- read.csv("tips_to_keep.csv", row.names = 2, header = T)

#check data
keep

#check names in list against tree
name_list<-name.check(tree,keep)

#make list of names not in data file
checked_names<-name_list$tree_not_data

#convert to a vector
names_to_drop<-as.vector(checked_names)

#trim tree to just taxa represented in list
trimmed_tree<-drop.tip(tree,names_to_drop)

#visualize trimmed tree
plot(trimmed_tree, cex=1)

#check trimmed tree against species list
name.check(trimmed_tree, keep)

write.tree(phy = trimmed_tree, file = "trimmed_tree.tre")
