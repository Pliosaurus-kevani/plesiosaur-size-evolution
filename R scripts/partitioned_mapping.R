#=============#
#load packages
#=============#

library(RevGadgets) #version 1.2.1; loaded for reading trees
library(ape) #version 5.8-1; loaded for calculating tree heights
library(ggplot2) #version 4.0.2; loaded for visualization
library(ggtree) #version 4.0.4; loaded for visualization

#------------------------------------------------------------------------------#

#==============#
#read the trees
#==============#

#change the working pathway

#setwd()

#read the trees
#cd, craniodental; ax, axial; ap, appendicular

map <- readTrees("map1.tre")
map.cd <- readTrees("map.craniodental.tre")
map.ax <- readTrees("map.axial.tre")
map.ap <- readTrees("map.appendicular.tre")

#extract the treedata object

treed <- map[[1]][[1]]
treed.cd <- map.cd[[1]][[1]]
treed.ax <- map.ax[[1]][[1]]
treed.ap <- map.ap[[1]][[1]]

#------------------------------------------------------------------------------#

#====================#
#plot (unpartitioned)
#====================#

#define the limit

global_limits <- c(-2.5, -1)

#find the root
#This is just for easier alignment in subsequent modification
#using Adobe Illustrator

root.treed <- Ntip(treed@phylo)+1

max <- max(node.height(treed@phylo))

#find some nodes with high evolutionary rates

Thalassophonea <- getMRCA(treed@phylo, c("Sachicasaurus_vitae", "Peloneustes_philarchus"))
Plesiosauroidea <- getMRCA(treed@phylo, c("Eoplesiosaurus_antiquior", "Plesiosaurus_dolichodeirus"))
Cryptoclidia <- getMRCA(treed@phylo, c("Leptocleidus_superstes", "Tricleidus_seeleyi"))
Polycotylidae <- getMRCA(treed@phylo, c("Edgarosaurus_muddi", "Trinacromerum_bentonianum"))

p_unpart <- 
  ggtree(treed, aes(color = log10(branch_rates)),
         layout = "fan",
         open.angle = 270, 
         linewidth = 0.5) +
  geom_hilight(
    node = root.treed,             
    fill = NA,         
    color = "grey",             
    linewidth = 0.5,
    extend = max*0.05                   
  ) +
  geom_nodepoint(
    aes(subset = (node %in% c(Thalassophonea, Plesiosauroidea, Cryptoclidia, Polycotylidae))), 
    shape = 21,                              
    fill = "black",                    
    color = "black",                         
    size = 0.8                             
  )+
  scale_x_continuous(limits = c(-max * 2, NA)) +
  scale_color_viridis_c(name = expression(log[10](rate:unpartitioned)),
                        limits = global_limits,
                        oob = scales::squish,
                        guide = guide_colorbar(direction = "horizontal",
                                               title.position = "top",         
                                               title.hjust = 0.5,              
                                               barwidth = unit(10, "lines"),   
                                               barheight = unit(1, "lines"))) +
  coord_polar(theta = "y", start = pi*3/2, direction = 1, clip = "off") +
  theme_void() + 
  theme(legend.position = "left",
        plot.margin = margin(10, 10, 10, 10))

p_unpart

#------------------------------------------------------------------------------#

#===================#
#plot (craniodental)
#===================#

#find the root
#This is just for easier alignment in subsequent modification
#using Adobe Illustrator

root.treecd <- Ntip(treed.cd@phylo)+1

#find some nodes with high evolutionary rates

Plesiosauria.cd <- getMRCA(treed.cd@phylo, c("Sachicasaurus_vitae", "Rhomaleosaurus_cramptoni"))
Thalassophonea.cd <- getMRCA(treed.cd@phylo, c("Sachicasaurus_vitae", "Peloneustes_philarchus"))
Pliosauridae.cd <- getMRCA(treed.cd@phylo, c("Sachicasaurus_vitae", "Thalassiodracon_hawkinsii"))
Plesiosauroidea.cd <- getMRCA(treed.cd@phylo, c("Styxosaurus_snowii", "Attenborosaurus_conybeari"))
Cryptoclididae.cd <- getMRCA(treed.cd@phylo, c("Muraenosaurus_leedsii", "Abyssosaurus_nataliae"))
Xenopsaria.cd <- getMRCA(treed.cd@phylo, c("Trinacromerum_bentonianum", "Styxosaurus_snowii"))
Aristonectinae.cd <- getMRCA(treed.cd@phylo, c("Aristonectes_quiriquinensis", "Wunyelfia_maulensis"))
Leptocleidia.cd <- getMRCA(treed.cd@phylo, c("Trinacromerum_bentonianum", "Brancasaurus_brancai"))
Elasmosauridae.cd <- getMRCA(treed.cd@phylo, c("Styxosaurus_snowii", "Speeton_Clay_Plesiosaurian"))

p_cd <- 
  ggtree(treed.cd, aes(color = log10(branch_rates_cd)),
         layout = "fan",
         open.angle = 270,
         linewidth = 0.5) +
  geom_hilight(
    node = root.treecd,             
    fill = NA,         
    color = "grey",             
    linewidth = 0.5,
    extend = max*0.05                   
  ) +
  geom_nodepoint(
    aes(subset = (node %in% c(Plesiosauria.cd, Thalassophonea.cd, 
                              Pliosauridae.cd, Plesiosauroidea.cd,
                              Cryptoclididae.cd, Xenopsaria.cd, 
                              Aristonectinae.cd, Leptocleidia.cd,
                              Elasmosauridae.cd))), 
    shape = 21,                              
    fill = "black",                    
    color = "black",                         
    size = 0.8                             
  )+
  scale_color_viridis_c(name = expression(log[10](rate:craniodental)),
                        limits = global_limits,
                        oob = scales::squish,
                        guide = guide_colorbar(direction = "horizontal",
                                               title.position = "top",         
                                               title.hjust = 0.5,              
                                               barwidth = unit(10, "lines"),   
                                               barheight = unit(1, "lines"))) +
  scale_x_continuous(limits = c(-max * 2, NA)) +
  coord_polar(theta = "y", start = 0, direction = 1, clip = "off") +
  theme_void() + 
  theme(legend.position = "right",
        plot.margin = margin(10, 10, 10, 10))

p_cd

#------------------------------------------------------------------------------#

#============#
#plot (axial)
#============#

#find the root
#This is just for easier alignment in subsequent modification
#using Adobe Illustrator

root.treeax <- Ntip(treed.ax@phylo)+1

#find some nodes with high evolutionary rates

Plesiosauroidea.ax <- getMRCA(treed.ax@phylo, c("Styxosaurus_snowii", "Attenborosaurus_conybeari"))
Leptocleidia.ax <- getMRCA(treed.ax@phylo, c("Trinacromerum_bentonianum", "Brancasaurus_brancai"))
Aristonectinae.ax <- getMRCA(treed.ax@phylo, c("Aristonectes_quiriquinensis", "Wunyelfia_maulensis"))

p_ax <- 
  ggtree(treed.ax, aes(color = log10(branch_rates_ax)),
         layout = "fan",
         open.angle = 270,
         linewidth = 0.5) +
  geom_hilight(
    node = root.treeax,             
    fill = NA,         
    color = "grey",             
    linewidth = 0.5,
    extend = max*0.05
  ) +
  geom_nodepoint(
    aes(subset = (node %in% c(Plesiosauroidea.ax, Leptocleidia.ax,
                              Aristonectinae.ax))), 
    shape = 21,                              
    fill = "black",                    
    color = "black",                         
    size = 0.8                             
  )+
  scale_color_viridis_c(name = expression(log[10](rate:axial)),
                        limits = global_limits,
                        oob = scales::squish,
                        guide = guide_colorbar(direction = "horizontal",
                                               title.position = "top",         
                                               title.hjust = 0.5,              
                                               barwidth = unit(10, "lines"),   
                                               barheight = unit(1, "lines"))) +
  scale_x_continuous(limits = c(-max * 2, NA)) +
  coord_polar(theta = "y", start = pi/2, direction = 1, clip = "off") +
  theme_void() + 
  theme(legend.position = "right",
        plot.margin = margin(10, 10, 10, 10))

p_ax

#------------------------------------------------------------------------------#

#===================#
#plot (appendicular)
#===================#

#find the root
#This is just for easier alignment in subsequent modification
#using Adobe Illustrator

root.treeap <- Ntip(treed.ap@phylo)+1

#find some nodes with high evolutionary rates

Thalassophonea.ap <- getMRCA(treed.ap@phylo, c("Sachicasaurus_vitae", "Peloneustes_philarchus"))
Cryptoclidia.ap <- getMRCA(treed.ap@phylo, c("Muraenosaurus_leedsii", "Trinacromerum_bentonianum"))

p_ap <- 
  ggtree(treed.ap, aes(color = log10(branch_rates_ap)),
         layout = "fan",
         open.angle = 270,
         linewidth = 0.5) +
  geom_hilight(
    node = root.treeap,             
    fill = NA,         
    color = "grey",             
    linewidth = 0.5,
    extend = max*0.05
  ) +
  geom_nodepoint(
    aes(subset = (node %in% c(Thalassophonea.ap, Cryptoclidia.ap))), 
    shape = 21,                              
    fill = "black",                    
    color = "black",                         
    size = 0.8                             
  )+
  scale_color_viridis_c(name = expression(log[10](rate:appendicular)),
                        limits = global_limits,
                        oob = scales::squish,
                        guide = guide_colorbar(direction = "horizontal",
                                               title.position = "top",         
                                               title.hjust = 0.5,              
                                               barwidth = unit(10, "lines"),   
                                               barheight = unit(1, "lines"))) +
  scale_x_continuous(limits = c(-max * 2, NA)) +
  coord_polar(theta = "y", start = pi, direction = 1, clip = "off") +
  theme_void() + 
  theme(legend.position = "left",
        plot.margin = margin(10, 10, 10, 10))

p_ap

#------------------------------------------------------------------------------#

#save the images

ggsave("unpartitioned.pdf", p_unpart, width = 15, height = 15)
ggsave("craniodental.pdf", p_cd, width = 15, height = 15)
ggsave("axial.pdf", p_ax, width = 15, height = 15)
ggsave("appendicular.pdf", p_ap, width = 15, height = 15)
