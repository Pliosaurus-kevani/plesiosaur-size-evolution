#==============#
#load packages
#==============#

library(tidytree) #version 0.4.7; loaded for loading trees
library(phytools) #version 2.5-2; loaded for stochastic character mapping
library(mvMORPH) #version 1.2.1; loaded for fitting multi-regime models
library(dplyr) #version 1.2.0; loaded for processing data
library(ggplot2) #version 4.0.2; loaded for visualization
library(hexbin) #version 1.28.5; loaded for visualization
library(viridis) #version 0.6.5; loaded for visualization

#------------------------------------------------------------------------------#

#=====================#
#read and process data
#=====================#

#change the working directory

#setwd()

#load the trees

load("tree_partitioned.RData")

#read the mass data (all 130 OTUs)

data <- read.csv("mass_all.csv", row.names = 1, header = T)

#log-transform the mass

data$mass <- log10(data$mass)

#change the clade labels

for(i in 1:nrow(data))
{
  if(data[i,1] != "Elasmosauridae" & data[i,1] != "Thalassophonea")
  {
    data[i,1] <- "others"
  }
}

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

tree.m <- list()

#The trees having the same index within tree.cd, tree.ax, and tree.ap
#are extactly the same
#so here I shall use tree.cd

for(i in 1:100)
{
  tree.m[[i]] <- make.simmap(tree.cd[[i]], clade.info, model = mkmodel)
  
  #report progress
  print(paste("working on tree", i, sep = " "))
}

#plot to check

cols3 <- setNames(c("#F34C3F","#8498AB","#5061c5"), levels(clade.info))

plot(ladderize.simmap(tree.m[[1]]), cols3, fsize = 0.4, ftype = "i", lwd = 2, offset = 0.4, ylim = c(-1, Ntip(tree.m[[1]])))
plot(ladderize.simmap(tree.m[[50]]), cols3, fsize = 0.4, ftype = "i", lwd = 2, offset = 0.4, ylim = c(-1, Ntip(tree.m[[1]])))
plot(ladderize.simmap(tree.m[[100]]), cols3, fsize = 0.4, ftype = "i", lwd = 2, offset = 0.4, ylim = c(-1, Ntip(tree.m[[1]])))

#------------------------------------------------------------------------------#

#=================================#
#fit 3-regime OU model to the data
#=================================#

OUs <- list()

for(i in 1:100)
{
  OUs[[i]] <- mvOU(tree.m[[i]], data[tree.m[[i]]$tip.label,2], model = "OUM", method = "pseudoinverse", param = list(vcv = "fixedRoot", root = "stationary"), scale.height = F, echo = F)
  
  #report progress
  print(paste("working on tree", i, sep = " "))
}

#------------------------------------------------------------------------------#

#==============================#
#ancestral state reconstruction
#==============================#

tree.ans <- list()

for(i in 1:100)
{
  tree.ans[[i]] <- estim(tree.m[[i]], data[tree.m[[i]]$tip.label,2], OUs[[i]], asr = T)
  
  #report progress
  print(paste("working on tree", i, sep = " "))
}

#In each case, a warning will be returned:
#Missing cases were first imputed before estimating the ancestral values!!
#However, this is not an error but a necessary step
#before ancestrak state reconstruction

#------------------------------------------------------------------------------#

#=======================================#
#calculate the mean mass for each branch
#=======================================#

branch_mass_list <- list()

for(i in 1:100)
{
  current_tree <- tree.m[[i]]
  current_ans <- tree.ans[[i]]
  
  N <- Ntip(current_tree)
  
  #extract the body mass for the tips (inlcuding NA values)
  tip_mass <- data[current_tree$tip.label, 2]
  
  #extract the body mass estimates for internal nodes
  node_mass <- as.numeric(current_ans$estimates)
  
  #combine the node indexes
  all_nodes_mass <- c(tip_mass, node_mass)
  
  #extract the tip labels
  #internal nodes are represented by node indexes
  all_labels <- c(current_tree$tip.label, as.character((N + 1):(2 * N - 1)))
  
  #constract the data frame
  edge_df <- data.frame(
    tree_id = i,                       
    parent = current_tree$edge[, 1],   
    node = current_tree$edge[, 2]    
  )
  
  #calculate the mean body mass of each branch
  edge_df <- edge_df %>%
    mutate(
      child_label = all_labels[node],           
      parent_mass = all_nodes_mass[parent],      
      child_mass  = all_nodes_mass[node],       
      branch_mass = (child_mass + parent_mass)/2 
    )
  
  #save the results
  branch_mass_list[[i]] <- edge_df
}

#------------------------------------------------------------------------------#

#========================================#
#combine the branch mass with branch rate
#========================================#

combined_data_cd <- list()

for(i in 1:100)
{
  mass_df <- branch_mass_list[[i]]
  
  rate_df <- as_tibble(TREE.cd[[i]])
  
  merged_df <- mass_df %>%
    inner_join(rate_df, by = c("parent", "node"))
  
  #also, log-transform the branch rate at this stage
  
  merged_df$branch_rates_cd <- log10(merged_df$branch_rates_cd)
  
  combined_data_cd[[i]] <- merged_df
}

combined_data_ax <- list()

for(i in 1:100)
{
  mass_df <- branch_mass_list[[i]]
  
  rate_df <- as_tibble(TREE.ax[[i]])
  
  merged_df <- mass_df %>%
    inner_join(rate_df, by = c("parent", "node"))
  
  #also, log-transform the branch rate at this stage
  
  merged_df$branch_rates_ax <- log10(merged_df$branch_rates_ax)
  
  combined_data_ax[[i]] <- merged_df
}

combined_data_ap <- list()

for(i in 1:100)
{
  mass_df <- branch_mass_list[[i]]
  
  rate_df <- as_tibble(TREE.ap[[i]])
  
  merged_df <- mass_df %>%
    inner_join(rate_df, by = c("parent", "node"))
  
  #also, log-transform the branch rate at this stage
  
  merged_df$branch_rates_ap <- log10(merged_df$branch_rates_ap)
  
  combined_data_ap[[i]] <- merged_df
}

#------------------------------------------------------------------------------#

#=============#
#plot: hexbins
#=============#

#merge the datasets

final_df_cd <- bind_rows(combined_data_cd) %>%
  filter(!is.na(branch_mass), !is.na(branch_rates_cd))

final_df_ax <- bind_rows(combined_data_ax) %>%
  filter(!is.na(branch_mass), !is.na(branch_rates_ax))

final_df_ap <- bind_rows(combined_data_ap) %>%
  filter(!is.na(branch_mass), !is.na(branch_rates_ap))

#plot: craniodental

hex_plot_cd <- 
  ggplot(final_df_cd, aes(x = branch_mass, y = branch_rates_cd)) +
  geom_hex(bins = 50, linewidth = 0.1) + 
  scale_fill_viridis_c(
    option = "magma",
    begin = 0.22,  
    end = 1.0,
    name = expression(N[branches])
  ) + 
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 11),
    legend.title = element_text(size = 9, face = "bold"),
    panel.grid.major = element_line(color = "grey95"), 
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "craniodental characters",
    x = expression(log[10]~(branch~mass)),
    y = expression(log[10]~(morphological~rate))
  )

print(hex_plot_cd)

#plot: axial

hex_plot_ax <- 
  ggplot(final_df_ax, aes(x = branch_mass, y = branch_rates_ax)) +
  geom_hex(bins = 50, linewidth = 0.1) + 
  scale_fill_viridis_c(
    option = "magma",
    begin = 0.22,  
    end = 1.0,
    name = expression(N[branches])
  ) + 
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 11),
    legend.title = element_text(size = 9, face = "bold"),
    panel.grid.major = element_line(color = "grey95"), 
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "axial characters",
    x = expression(log[10]~(branch~mass)),
    y = expression(log[10]~(morphological~rate))
  )

print(hex_plot_ax)

#plot: appendicular

hex_plot_ap <- 
  ggplot(final_df_ap, aes(x = branch_mass, y = branch_rates_ap)) +
  geom_hex(bins = 50, linewidth = 0.1) + 
  scale_fill_viridis_c(
    option = "magma",
    begin = 0.22,  
    end = 1.0,
    name = expression(N[branches])
  ) + 
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 11),
    legend.title = element_text(size = 9, face = "bold"),
    panel.grid.major = element_line(color = "grey95"), 
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "appendicular characters",
    x = expression(log[10]~(branch~mass)),
    y = expression(log[10]~(morphological~rate))
  )

print(hex_plot_ap)

#------------------------------------------------------------------------------#

#=========================================================#
#plot: distribution of the highest morphological rates
#=========================================================#

#extract the highest 5% of rates from each dataset

top_5_cd <- 
  final_df_cd %>%
  group_by(tree_id) %>%
  filter(branch_rates_cd >= quantile(branch_rates_cd, probs = 0.95, na.rm = TRUE)) %>%
  ungroup()

top_5_ax <- 
  final_df_ax %>%
  group_by(tree_id) %>%
  filter(branch_rates_ax >= quantile(branch_rates_ax, probs = 0.95, na.rm = TRUE)) %>%
  ungroup()

top_5_ap <- 
  final_df_ap %>%
  group_by(tree_id) %>%
  filter(branch_rates_ap >= quantile(branch_rates_ap, probs = 0.95, na.rm = TRUE)) %>%
  ungroup()

comp_cd <- ggplot() +
  geom_histogram(
    data = final_df_cd, aes(x = branch_mass, y = after_stat(density)),
    fill = "grey90", color = "white", bins = 40, alpha = 1
  ) +
  geom_density(
    data = top_5_cd, aes(x = branch_mass),
    fill = "#83D350",   
    color = "#519129",   
    alpha = 0.35,        
    linewidth = 1
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 15),
    legend.title = element_text(size = 9, face = "bold")
  ) +
  labs(
    title = "craniodental rate bursts across the body mass landscape",
    x = expression(log[10]~(branch~mass)),
    y = "density"
  )

comp_cd

comp_ax <- ggplot() +
  geom_histogram(
    data = final_df_ax, aes(x = branch_mass, y = after_stat(density)),
    fill = "grey90", color = "white", bins = 40, alpha = 1
  ) +
  geom_density(
    data = top_5_ax, aes(x = branch_mass),
    fill = "#83D350",   
    color = "#519129",   
    alpha = 0.35,        
    linewidth = 1
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 15),
    legend.title = element_text(size = 9, face = "bold")
  ) +
  labs(
    title = "axial rate bursts across the body mass landscape",
    x = expression(log[10]~(branch~mass)),
    y = "density"
  )

comp_ax

comp_ap <- ggplot() +
  geom_histogram(
    data = final_df_ap, aes(x = branch_mass, y = after_stat(density)),
    fill = "grey90", color = "white", bins = 40, alpha = 1
  ) +
  geom_density(
    data = top_5_ap, aes(x = branch_mass),
    fill = "#83D350",   
    color = "#519129",   
    alpha = 0.35,        
    linewidth = 1
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 15),
    legend.title = element_text(size = 9, face = "bold")
  ) +
  labs(
    title = "appendicular rate bursts across the body mass landscape",
    x = expression(log[10]~(branch~mass)),
    y = "density"
  )

comp_ap

#save the plots as RDS files

all_plots <- list(
  hex_cd  = hex_plot_cd,
  hex_ax  = hex_plot_ax,
  hex_ap  = hex_plot_ap,
  comp_cd = comp_cd,
  comp_ax = comp_ax,
  comp_ap = comp_ap
)

#setwd()

saveRDS(all_plots, file = "morpho_rate_partitioned_all branch.rds")

