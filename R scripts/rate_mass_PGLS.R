#=============#
#load packages
#=============#

library(tidytree) #version 0.4.7; loaded for processing tree data
library(ape) #version 5.8-1; loaded for processing trees
library(phylolm) #version 2.6.5; loaded for performing phylogenetic generalized least squares

#------------------------------------------------------------------------------#

#================#
#process the data
#================#

#change working pathway

#setwd()

#load the unpartitioned trees

load("tree.RData")

#load the partitioned trees

load("tree_partitioned.RData")

#read the body mass data

mass.data <- read.csv("mass.csv", header = T, row.names = 1)

mass.data <- mass.data[,-1]

#log-transform the data

mass.data$mass <- log10(mass.data$mass)

#prune the trees

tree.p <- ape::keep.tip.multiPhylo(tree, row.names(mass.data))
tree.cd.p <- ape::keep.tip.multiPhylo(tree.cd, row.names(mass.data))
tree.ax.p <- ape::keep.tip.multiPhylo(tree.ax, row.names(mass.data))
tree.ap.p <- ape::keep.tip.multiPhylo(tree.ap, row.names(mass.data))

#plot to check

#note that the craniodental, axial, and appendicular trees
#should be the same

plot(ladderize(tree.p[[1]]), no.margin = T)
plot(ladderize(tree.cd.p[[1]]), no.margin = T)
plot(ladderize(tree.ax.p[[1]]), no.margin = T)
plot(ladderize(tree.ap.p[[1]]), no.margin = T)

#------------------------------------------------------------------------------#

#============#
#perform PGLS
#============#

#create lists to contain the regression results

RPR <- list()
RPR.cd <- list()
RPR.ax <- list()
RPR.ap <- list()

#Here I employ a family of temporary variables called "moa" as counters

moa <- 0
moa.cd <- 0
moa.ax <- 0
moa.ap <- 0

#unpartitioned

for(i in 1:100)
{
  #Here I employ a family of temporary variables called "auk"
  #and the family "dodo"
  
  auk <- as_tibble(TREE[[i]])
  
  #remove the internal branches
  #we shall just use the external branches (tips) here
  
  auk <- na.omit(auk)
  auk <- setNames(auk$branch_rates, auk$label)
  
  #a warning will be returned:
  #Setting row names on a tibble is deprecated
  #but this does not affect the result
  
  #extract the data
  dodo.mass <- data.frame(mass = mass.data[tree.p[[i]]$tip.label,], row.names = tree.p[[i]]$tip.label)
  dodo.morph <- data.frame(auk[tree.p[[i]]$tip.label])
  
  #log-transform the morphological rate
  dodo.morph[,1] <- log10(dodo.morph[,1])
  
  #combine the datasets
  
  dodo <- cbind(dodo.mass, dodo.morph)
  
  colnames(dodo) <- c("clade", "mass", "morph_rate")
  
  #perform robust phylogenetic regresion
  
  RPR[[i]] <- phylolm(formula = morph_rate ~ mass, 
                      data = dodo[tree.p[[i]]$tip.label,], 
                      phy = tree.p[[i]], model = "lambda")
  
  #count the number of cases which fail to reject the null hypothesis
  
  if(summary(RPR[[i]])$coefficients[2,4] > 0.05)
  {
    moa <- moa + 1
  }
}

#craniodental

for(i in 1:100)
{
  #Here I employ a temporary variable called "auk.cd."
  #and the family "dodo"
  
  auk.cd <- as_tibble(TREE.cd[[i]])
  
  #remove the internal branches
  #we shall just use the external branches (tips) here
  
  auk.cd <- na.omit(auk.cd)
  auk.cd <- setNames(auk.cd$branch_rates_cd, auk.cd$label)
  
  #a warning will be returned:
  #Setting row names on a tibble is deprecated
  #but this does not affect the result
  
  #extract the data
  dodo.mass <- data.frame(mass = mass.data[tree.cd.p[[i]]$tip.label,], row.names = tree.cd.p[[i]]$tip.label)
  dodo.morph <- data.frame(auk.cd[tree.cd.p[[i]]$tip.label])
  
  #log-transform the morphological rate
  dodo.morph[,1] <- log10(dodo.morph[,1])
  
  #combine the datasets
  
  dodo <- cbind(dodo.mass, dodo.morph)
  
  colnames(dodo) <- c("clade", "mass", "morph_rate")
  
  #perform robust phylogenetic regresion
  
  RPR.cd[[i]] <- phylolm(formula = morph_rate ~ mass,
                         data = dodo[tree.cd.p[[i]]$tip.label,],
                         phy = tree.cd.p[[i]], model = "lambda")
  
  #count the number of cases which fail to reject the null hypothesis
  
  if(summary(RPR.cd[[i]])$coefficients[2,4] > 0.05)
  {
    moa.cd <- moa.cd + 1
  }
}

#axial

for(i in 1:100)
{
  #Here I employ a temporary variable called "auk.ax.ax"
  #and the family "dodo"
  
  auk.ax <- as_tibble(TREE.ax[[i]])
  
  #remove the internal branches
  #we shall just use the external branches (tips) here
  
  auk.ax <- na.omit(auk.ax)
  auk.ax <- setNames(auk.ax$branch_rates_ax, auk.ax$label)
  
  #a warning will be returned:
  #Setting row names on a tibble is deprecated
  #but this does not affect the result
  
  #extract the data
  dodo.mass <- data.frame(mass = mass.data[tree.ax.p[[i]]$tip.label,], row.names = tree.ax.p[[i]]$tip.label)
  dodo.morph <- data.frame(auk.ax[tree.ax.p[[i]]$tip.label])
  
  #log-transform the morphological rate
  dodo.morph[,1] <- log10(dodo.morph[,1])
  
  #combine the datasets
  
  dodo <- cbind(dodo.mass, dodo.morph)
  
  colnames(dodo) <- c("clade", "mass", "morph_rate")
  
  #perform robust phylogenetic regresion
  
  RPR.ax[[i]] <- phylolm(formula = morph_rate ~ mass,
                         data = dodo[tree.ax.p[[i]]$tip.label,],
                         phy = tree.ax.p[[i]], model = "lambda")
  
  #count the number of cases which fail to reject the null hypothesis
  
  if(summary(RPR.ax[[i]])$coefficients[2,4] > 0.05)
  {
    moa.ax <- moa.ax + 1
  }
}

#appendicular

for(i in 1:100)
{
  #Here I employ a temporary variable called "auk.ap.ap"
  #and the family "dodo"
  
  auk.ap <- as_tibble(TREE.ap[[i]])
  
  #remove the internal branches
  #we shall just use the external branches (tips) here
  
  auk.ap <- na.omit(auk.ap)
  auk.ap <- setNames(auk.ap$branch_rates_ap, auk.ap$label)
  
  #a warning will be returned:
  #Setting row names on a tibble is deprecated
  #but this does not affect the result
  
  #extract the data
  dodo.mass <- data.frame(mass = mass.data[tree.ap.p[[i]]$tip.label,], row.names = tree.ap.p[[i]]$tip.label)
  dodo.morph <- data.frame(auk.ap[tree.ap.p[[i]]$tip.label])
  
  #log-transform the morphological rate
  dodo.morph[,1] <- log10(dodo.morph[,1])
  
  #combine the datasets
  
  dodo <- cbind(dodo.mass, dodo.morph)
  
  colnames(dodo) <- c("clade", "mass", "morph_rate")
  
  #perform robust phylogenetic regresion
  
  RPR.ap[[i]] <- phylolm(formula = morph_rate ~ mass,
                         data = dodo[tree.ap.p[[i]]$tip.label,],
                         phy = tree.ap.p[[i]], model = "lambda")
  
  #count the number of cases which fail to reject the null hypothesis
  
  if(summary(RPR.ap[[i]])$coefficients[2,4] > 0.05)
  {
    moa.ap <- moa.ap + 1
  }
}

