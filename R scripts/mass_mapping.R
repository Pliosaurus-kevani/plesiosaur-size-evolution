#==============================#
#load packages
#==============================#

library(RevGadgets) #version 1.2.1; loaded for reading the tree
library(ape) #version 5.8.1; loaded for processing the tree
library(mvMORPH) #version 1.2.1; loaded for ancestral state reconstruction
library(ggtree) #version 4.04; loaded for plotting the tree
library(tidytree) #version 0.4.7; loaded for plotting the tree
library(ggplot2) #version 4.0.2; loaded for plotting the tree
library(RColorBrewer) #version 1.1-3; loaded for loading color palette

#------------------------------------------------------------------------------#

#==============================#
#read and process data
#==============================#

#change working directory

#setwd()

#read the maximum a posteriori tree

MAPT <- readTrees("map1.tre")

#extract the tree topology

tree <- MAPT[[1]][[1]]@phylo

#read body mass data

data <- read.csv("mass.csv", row.names = 1, header = T)

#log-transform the mass

data$mass <- log10(data$mass)

#change the clade labels

for(i in 1:nrow(data))
{
  if(data[i,2] != "Elasmosauridae" & data[i,2] != "Thalassophonea")
  {
    data[i,2] <- "others"
  }
}

data <- data[,-1]

#prune the tree

tree.p <- ape::keep.tip(tree, row.names(data))

#extract body mass

mass <- setNames(data$mass, row.names(data))

#------------------------------------------------------------------------------#

#============================#
#stochastic character mapping
#============================#

#define the Mk model
#here I forbid the state to change among the three regimes
#since the regimes are defined based on monophyletic groups

mkmodel <- matrix(c(0,0,0,1,0,2,0,0,0), 3, 3, byrow = T, dimnames = list(0:2, 0:2))

#extract the information and convert it into factors

clade.info <- setNames(data$clade, row.names(data))

clade.info <- as.factor(clade.info)

#stochastic character mapping

tree.m <- make.simmap(tree.p, clade.info, model = mkmodel)

#plot to check

cols3 <- setNames(c("#F34C3F","#8498AB","#5061c5"), levels(clade.info))

plot(ladderize.simmap(tree.m), cols3, fsize = 0.4, ftype = "i", lwd = 2, offset = 0.4, ylim = c(-1, Ntip(tree.m)))

#------------------------------------------------------------------------------#

#==============================#
#ancestral state reconstruction
#==============================#

#fit an OU3 model to the data

OU3 <- mvOU(tree.m, data[tree.m$tip.label,2], model = "OUM", method = "pseudoinverse", param = list(vcv = "fixedRoot", root = "stationary"), scale.height = F, echo = F)

#ancestral state reconstruction

ans <- estim(tree.m, data[tree.m$tip.label,2], OU3, asr = T)

#=============================================================#
#construct a data frame that contains both tip and node values
#=============================================================#

#extract the estimated values for internal nodes

node_values <- as.numeric(ans$estim)

#generate node IDs for internal nodes

node_ids <- (Ntip(tree.m) + 1):(Ntip(tree.m) + tree.m$Nnode)
nodes <- data.frame(node = node_ids, trait = node_values)

#extract the actual observed mass values for the tips

tip_values <- as.numeric(data[tree.m$tip.label, 2])

#tip IDs are simply 1 to N

tip_ids <- 1:Ntip(tree.m)
tips <- data.frame(node = tip_ids, trait = tip_values)

#combine the data frames

data.p <- rbind(tips, nodes)

#combine the continuous data with the simmap tree object

tree.plot <- full_join(tree.m, data.p, by = 'node')

#define the clade node

Pliosauridae.node <- getMRCA(tree.p, c("Sachicasaurus_vitae", "Thalassiodracon_hawkinsii"))
Rhomaleosauridae.node <- getMRCA(tree.p, c("Macroplata_tenuiceps", "Rhomaleosaurus_cramptoni"))
Cryptoclididae.node <- getMRCA(tree.p, c("Abyssosaurus_nataliae", "Tricleidus_seeleyi"))
Polycotylidae.node <- getMRCA(tree.p, c("Dolichorhynchops_osborni", "Edgarosaurus_muddi"))
Elasmosauridae.node <- getMRCA(tree.p, c("Aristonectes_quiriquinensis", "Wapuskanectes_betsynichollsae"))

clade.node <- data.frame(node = c(Pliosauridae.node,
                                  Rhomaleosauridae.node,
                                  Cryptoclididae.node,
                                  Polycotylidae.node,
                                  Elasmosauridae.node), 
                         clade = c("Pliosauridae",
                                  "Rhomaleosauridae",
                                  "Cryptoclididae",
                                  "Polycotylidae",
                                  "Elasmosauridae"))

#------------------------------------------------------------------------------#

#=====#
#plot
#=====#

#extract the "Set3" color palette

set3_colors <- brewer.pal(8, "Set3")

#define the label size

label_size_mm <- 1.5

p <- ggtree(tree.plot, layout="circular", ladderize=FALSE) +
  geom_tree(aes(color=trait), linewidth = 0.3) +
  scale_color_viridis_c(name = expression(log[10](mass))) +
  geom_hilight(
    data = clade.node,
    aes(node=node, fill=clade),
    alpha = 0.35,
    colour = "lightgrey",
    linewidth = 0.1
  ) +
  scale_fill_manual(values=set3_colors, guide = "none")+
  geom_cladelab(node = Cryptoclididae.node,
                label = "Cryptoclididae",
                hjust = 1,
                barcolor = NA,
                geom = "label",
                offset.text = max(nodeHeights(tree.p))*0.05,
                fill = set3_colors[1],
                alpha = 0.35,
                fontsize = label_size_mm)+
  geom_cladelab(node = Elasmosauridae.node,
                label = "Elasmosauridae",
                hjust = 0,
                barcolor = NA,
                geom = "label",
                offset.text = max(nodeHeights(tree.p))*0.05,
                fill = set3_colors[2],
                alpha = 0.35,
                fontsize = label_size_mm)+
  geom_cladelab(node = Pliosauridae.node,
                label = "Pliosauridae",
                hjust = 0.5,
                vjust = 1,
                barcolor = NA,
                geom = "label",
                offset.text = max(nodeHeights(tree.p))*0.05,
                fill = set3_colors[3],
                alpha = 0.35,
                fontsize = label_size_mm)+
  geom_cladelab(node = Polycotylidae.node,
                label = "Polycotylidae",
                hjust = 1,
                barcolor = NA,
                geom = "label",
                offset.text = max(nodeHeights(tree.p))*0.05,
                fill = set3_colors[4],
                alpha = 0.35,
                fontsize = label_size_mm)+
  geom_cladelab(node = Rhomaleosauridae.node,
                label = "Rhomaleosauridae",
                hjust = 0,
                barcolor = NA,
                geom = "label",
                offset.text = max(nodeHeights(tree.p))*0.05,
                fill = set3_colors[5],
                alpha = 0.35,
                fontsize = label_size_mm)+
  theme_void() +
  theme(
    legend.text = element_text(size = 6),        
    legend.title = element_text(size = 7),
    plot.margin = margin(0, 0, 0, 0)
  )

#A warning "Ignoring unknown parameters: `label.size`"
#will be returned
#but this does not affect plotting

p

#save the figure

saveRDS(p, "mass.mapping.rds")
