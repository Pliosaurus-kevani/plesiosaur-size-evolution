#=============#
#load packages
#=============#

library(ggplot2) #version 4.0.2; loaded for visualization
library(patchwork) #version 1.3.2; loaded for visualization

#------------------------------------------------------------------------------#

#change the working pathway

#setwd()

#read the plots

tree.plot <- readRDS("mass.mapping.rds")
aicc.plot <- readRDS("aicc.plot.rds")
bar.plot <- readRDS("bar.plot.rds")
mass.plot <- readRDS("mass.distribution.rds")
alpha.plot <- readRDS("alpha.plot.rds")
theta.plot <- readRDS("theta.plot.rds")

#------------------------------------------------------------------------------#

#====#
#plot
#====#

#remove the legend of some figures

mass.plot <- 
  mass.plot +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11))

aicc.plot <- 
  aicc.plot + 
  guides(fill = "none", color = "none")+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))

bar.plot <- 
  bar.plot + 
  guides(fill = "none", color = "none") + 
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))
  
  
alpha.plot <- 
  alpha.plot + 
  guides(fill = "none", color = "none") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))

theta.plot <-
  theta.plot + 
  guides(fill = "none", color = "none") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))


#partition the space

A = tree.plot
B = mass.plot 
C = aicc.plot  
D = bar.plot 
E = theta.plot
F = alpha.plot

layout_1_5 <- "
AAD
AAE
BCF
"

#plot

final_figure <- tree.plot + mass.plot + aicc.plot + bar.plot + theta.plot + alpha.plot + 
  
  plot_layout(
    design = layout_1_5,
    guides = "collect"   
  ) +
  plot_annotation(
    tag_levels = 'A'
  ) & 
  theme(
    legend.position = "bottom",    
    legend.box = "horizontal",
    legend.spacing.x = unit(0.5, "cm"), 
    plot.tag = element_text(size = 16, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA)
  )

final_figure

#save the plot

ggsave(
  filename = "figure2_new2.pdf", 
  plot = final_figure, 
  width = 12, 
  height = 8, 
  dpi = 300,
  device = cairo_pdf 
)

#some further modifications were carried out using Adobe Illustrator
