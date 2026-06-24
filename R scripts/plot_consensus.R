#=============#
#load packages
#=============#

library(RevGadgets) #version 1.2.1; loaded for reading the tree
library(ggtree) #version 4.0.4; loaded for visualization
library(deeptime) #version 2.3.1; loaded for visualization
library(ggplot2) #version 4.0.2; loaded for visualization
 
#------------------------------------------------------------------------------#

#======================#
#read and plot the tree
#======================#

#change the working directory

#setwd()

#read the tree

consensus <- readTrees("consensus_combined.tre")

tree <- consensus[[1]][[1]]

#remove the "_" in tip labels

tree@phylo$tip.label <- gsub("_", " ", tree@phylo$tip.label)

#create an initial plot

p <- ggtree(tree, ladderize = T, right = T)

#calculate the tree height

tree_height <- max(p$data$x)

#change the root time of the tree

p$data$x <- p$data$x - tree_height - 66

#plot

max_y <- max(p$data$y)

ages_dat <- get_scale_data("ages")

p_geo <- p +
  geom_rect(
    data = ages_dat,
    aes(xmin = -max_age, xmax = -min_age, ymin = -Inf, ymax = Inf),
    fill = rep(c("gray92", "white"), length.out = nrow(ages_dat)),
    inherit.aes = F,
    color = NA 
  ) +
  geom_tree() +
  geom_tippoint(color = "black", size = 2) +
  geom_tiplab(size = 2.5, offset = 0.5, fontface = "italic") + 
  theme_tree2() +
  scale_y_continuous(limits = c(-1, max_y + 1)) +
  coord_geo(
    dat = list("ages", "periods"), 
    xlim = c(min(p$data$x), -45),
    pos = "bottom",
    neg = T,     
    abbrv = list(T, F),
    clip = "off",
    size = list(3, 4)
  ) +
  scale_x_continuous(labels = abs, name = "Time (Ma)") +
  theme(
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold")
  )

p_geo

#save the plot

#setwd()

ggsave("consensus.pdf", p_geo, width = 20, height = 15)
