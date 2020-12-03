##set working directory
setwd("~/Desktop/Wes/Jessie")

##load geiger and dependencies
library("geiger")


##load tree
tree<-read.tree("4ktree.tre")
tree<-ladderize(tree,right=F)

##check tree
plot(tree, cex=0.3)


##read in taxon list
data<-read.csv("taxa_for_pic.csv", header = T, row.names = 1)


##check data
data

###check names in list against tree
name_list<-name.check(tree,data)

##make list of names not in data file
checked_names<-name_list$tree_not_data


##convert to a vector
names_to_drop<-as.vector(checked_names)


##trim tree to just taxa represented in list
trimmed_tree<-drop.tip(tree,names_to_drop)

##visualize trimmed tree
plot(trimmed_tree, cex=0.25)


##check trimmed tree against species list
name.check(trimmed_tree, data)



##read in trimmed data (no missing allowed for PIC)
data2<-read.csv("taxa_for_pic2.csv", header = T, row.names = 1)

###check names in list against trimmed tree
name_list<-name.check(trimmed_tree,data2)

##make list of names not in data file
checked_names<-name_list$tree_not_data


##convert to a vector
names_to_drop<-as.vector(checked_names)


##trim tree to just taxa represented in list
trimmed_tree2<-drop.tip(trimmed_tree,names_to_drop)

##visualize trimmed tree2
plot(trimmed_tree2, cex=0.25)


##check trimmed tree2 against species list
name.check(trimmed_tree2, data2)



###perform regression without correcting for phylogeny
x<-data2$Mean_fst
y<-data2$speciation_rate
plot(x,y)
fit<-lm(y~x)
abline(fit)

##show summary of regression
summary(fit)

##perform ANOVA
anova(fit)


###correct for phylogeny using Phylogenetic Independent Contrasts
ix<-pic(x,trimmed_tree2)
iy<-pic(y,trimmed_tree2)
fit_PIC<-lm(iy~ix-1)

##show summary of fit with PIC
summary(fit_PIC)


##plot your PIC
plot(iy~ix)
abline(a=0, b=coef(fit_PIC))




