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
library(ggdist) #version 3.3.3; loaded for visualization
library(patchwork) #version 1.3.2; loaded for visualization

#------------------------------------------------------------------------------#

#=====================#
#read and process data
#=====================#

#change the working directory

#setwd()

#load the trees

load("tree.RData")

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

for(i in 1:100)
{
  tree.m[[i]] <- make.simmap(tree[[i]], clade.info, model = mkmodel)
  
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
#before ancestral state reconstruction

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

#==================================================#
#plot: branch mass distribution across four subsets
#==================================================#

#define the clades

Cryptoclididae_tips <- c("Cryptoclidus_eurymerus","Tricleidus_seeleyi",
                         "Muraenosaurus_leedsii","Kimmerosaurus_langhami",
                         "Pantosaurus_striatus","Picrocleidus_beloclis",
                         "Tatenectes_laramiensis","Plesiosaurus_mansellii",
                         "Colymbosaurus_megadeirus","Djupedalia_engeri",
                         "Spitrasaurus_spp","Ophthalmothule_cryostea","Abyssosaurus_nataliae")  

Thalassophonea_tips <- c("Peloneustes_philarchus","Simolestes_vorax",
                         "Pliosaurus_funkei","Pliosaurus_westburyensis",
                         "Pliosaurus_brachydeirus","Pliosaurus_almanzaensis",
                         "Gallardosaurus_iturraldei","Pliosaurus_rossicus",
                         "Pliosaurus_irgisensis","Pliosaurus_andrewsi",
                         "Luskhan_itilensis","Makhaira_rossica",
                         "Liopleurodon_ferox","Kronosaurus_MCZ_1285",
                         "Brachauchenius_lucasi","Brachauchenius_MNA_V9433",
                         "QM_F51291","Sachicasaurus_vitae",
                         "Acostasaurus_pavachoquensis","Stenorhynchosaurus_munozi")

Elasmosauridae_tips <- c("Speeton_Clay_Plesiosaurian","Wapuskanectes_betsynichollsae",
                         "Eromangasaurus_australis","Tuarangisaurus_keyesi",
                         "Elasmosaurus_platyurus","Terminonatator_ponteixensis",
                         "Libonectes_morgani","Kaiwhekea_katiki",
                         "Aristonectes_quiriquinensis","Aristonectes_parvidens",
                         "Morturneria_seymourensis","Albertonectes_vanderveldei",
                         "Kawanectes_lafquenianum","Vegasaurus_molyi",
                         "Morenosaurus_stocki","Hydrotherosaurus_alexandrae",
                         "Wunyelfia_maulensis","Aphrosaurus_furlongi",
                         "Zarafasaura_oceanis","Futabasaurus_suzukii",
                         "Callawayasaurus_colombiensis","MLP_99_XII_1_5",
                         "Alexandronectes_zealandiensis","Fluvionectes_sagecrensis",
                         "Cardiocorax_mukulu","Jucha_squalea",
                         "Chubutinectes_carmeloi","Marambionectes_molinai",
                         "Lagenanectes_richterae","Thalassomedon_haningtoni",
                         "Nakonanectes_bradti","Styxosaurus_snowii",
                         "Styxosaurus_rezaci","Styxosaurus_SDSM_451",
                         "Styxosaurus_browni","Traskasaura_sandrae")

#extract the branch masses across the 100 trees

all_mass_df <- data.frame()

for(i in 1:100) {
  current_tree <- tree.m[[i]]
  edge_df <- branch_mass_list[[i]]
  
  edge_df$clade <- "Remaining taxa"
  
  get_clade_nodes <- function(tip_list) {
    
    valid_tips <- intersect(tip_list, current_tree$tip.label)
    
    if(length(valid_tips) > 1) {

      mrca_node <- getMRCA(current_tree, valid_tips)

      descendants <- getDescendants(current_tree, mrca_node)

      return(c(mrca_node, descendants))
    }
    return(numeric(0))
  }
  
  crypto_nodes <- get_clade_nodes(Cryptoclididae_tips)
  thalasso_nodes <- get_clade_nodes(Thalassophonea_tips)
  elasmo_nodes <- get_clade_nodes(Elasmosauridae_tips)
  
  edge_df$clade[edge_df$node %in% crypto_nodes] <- "Cryptoclididae"
  edge_df$clade[edge_df$node %in% thalasso_nodes] <- "Thalassophonea"
  edge_df$clade[edge_df$node %in% elasmo_nodes] <- "Elasmosauridae"
  
  all_mass_df <- rbind(all_mass_df, edge_df)
}

#remove the NA values (if present)

all_mass_clean <- all_mass_df %>%
  filter(!is.na(branch_mass))

#order the clades

clade_levels <- c("Cryptoclididae", "Thalassophonea", "Elasmosauridae", "Remaining taxa")
all_mass_clean$clade <- factor(all_mass_clean$clade, levels = rev(clade_levels))

#define the colors

my_colors <- c("Cryptoclididae" = "#D2A16E",   
               "Thalassophonea" = "#5061c5",   
               "Elasmosauridae" = "#F34C3F",   
               "Remaining taxa" = "#8498AB")

#raincloud plot for branch mass distribution

p_mass_raincloud <- ggplot(all_mass_clean, aes(x = branch_mass, y = clade, fill = clade)) +
  
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, scale = 0.45, justification = -0.35, point_colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5, color = "grey30") +
  geom_jitter(data = all_mass_clean[sample(nrow(all_mass_clean), min(5000, nrow(all_mass_clean))), ], 
              aes(color = clade), width = 0, height = 0.08, alpha = 0.15, size = 0.5) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  coord_cartesian(ylim = c(1.2, 4.1)) +
  theme_bw() + 
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12)) +
  labs(title = "body mass variations across clades", 
       x = expression(log[10]~(branch~mass)))

print(p_mass_raincloud)

#setwd()

saveRDS(p_mass_raincloud, file = "mass.raincloud.rds")

#------------------------------------------------------------------------------#

#========================================#
#combine the branch mass with branch rate
#========================================#

combined_data <- list()

for(i in 1:100)
{
  mass_df <- branch_mass_list[[i]]
  
  rate_df <- as_tibble(TREE[[i]])
  
  merged_df <- mass_df %>%
    inner_join(rate_df, by = c("parent", "node"))
  
  #also, log-transform the branch rate at this stage
  
  merged_df$branch_rates <- log10(merged_df$branch_rates)
  
  combined_data[[i]] <- merged_df
}

#------------------------------------------------------------------------------#

#=============#
#plot: hexbins
#=============#

#merge the datasets

final_df <- bind_rows(combined_data) %>%
  filter(!is.na(branch_mass), !is.na(branch_rates))

#plot

hex_plot <- 
  ggplot(final_df, aes(x = branch_mass, y = branch_rates)) +
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
    title = "all characters",
    x = expression(log[10]~(branch~mass)),
    y = expression(log[10]~(morphological~rate))
  )

print(hex_plot)

#------------------------------------------------------------------------------#

#=========================================================#
#plot: distribution of the highest morphological rates
#=========================================================#

#extract the highest 5% of rates from each dataset

top_5 <- 
  final_df %>%
  group_by(tree_id) %>%
  filter(branch_rates >= quantile(branch_rates, probs = 0.95, na.rm = TRUE)) %>%
  ungroup()

comp <- 
  ggplot() +
  geom_histogram(
    data = final_df, aes(x = branch_mass, y = after_stat(density)),
    fill = "grey90", color = "white", bins = 40, alpha = 1
  ) +
  geom_density(
    data = top_5, aes(x = branch_mass),
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
    title = "morphological rate bursts across the body mass landscape",
    x = expression(log[10]~(branch~mass)),
    y = "density"
  )

comp

#------------------------------------------------------------------------------#

#=================#
#combine the plots
#=================#

#read the plots

#setwd()

loaded_plots <- readRDS("morpho_rate_partitioned_all branch.rds")

#extract the data for some minor adjustments

final_df_cd  <- loaded_plots$hex_cd$data
final_df_ax  <- loaded_plots$hex_ax$data
final_df_ap  <- loaded_plots$hex_ap$data

top_5_cd  <- loaded_plots$comp_cd$layers[[2]]$data
top_5_ax  <- loaded_plots$comp_ax$layers[[2]]$data
top_5_ap  <- loaded_plots$comp_ap$layers[[2]]$data

#adjust the plots minorly

hex_plot <-
  ggplot(final_df, aes(x = branch_mass, y = branch_rates)) +
  geom_hex(bins = 50, linewidth = 0.1, show.legend = FALSE) + 
  scale_fill_viridis_c(option = "magma", begin = 0.22, end = 1.0, name = expression(N[branches])) + 
  #geom_smooth(method = "loess", color = "#00FFCC", se = FALSE, linewidth = 1.2, linetype = "dashed") +
  theme_bw() + theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5), axis.title = element_text(face = "bold", size = 11), legend.title = element_text(size = 10, face = "bold"), panel.grid.major = element_line(color = "grey95"), panel.grid.minor = element_blank()) +
  labs(title = "all characters", x = expression(log[10]~(branch~mass)), y = expression(log[10]~(morphological~rate)))

hex_plot_cd <- 
  ggplot(final_df_cd, aes(x = branch_mass, y = branch_rates_cd)) +
  geom_hex(bins = 50, linewidth = 0.1, show.legend = FALSE) + 
  scale_fill_viridis_c(option = "magma", begin = 0.22, end = 1.0) + 
  #geom_smooth(method = "loess", color = "#00FFCC", se = FALSE, linewidth = 1.2, linetype = "dashed") +
  theme_bw() + theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5), axis.title = element_text(face = "bold", size = 11), panel.grid.major = element_line(color = "grey95"), panel.grid.minor = element_blank()) +
  labs(title = "craniodental characters", x = expression(log[10]~(branch~mass)), y = expression(log[10]~(morphological~rate)))

hex_plot_ax <- 
  ggplot(final_df_ax, aes(x = branch_mass, y = branch_rates_ax)) +
  geom_hex(bins = 50, linewidth = 0.1, show.legend = FALSE) + 
  scale_fill_viridis_c(option = "magma", begin = 0.22, end = 1.0) +
  #geom_smooth(method = "loess", color = "#00FFCC", se = FALSE, linewidth = 1.2, linetype = "dashed") +
  theme_bw() + theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5), axis.title = element_text(face = "bold", size = 11), panel.grid.major = element_line(color = "grey95"), panel.grid.minor = element_blank()) +
  labs(title = "axial characters", x = expression(log[10]~(branch~mass)), y = expression(log[10]~(morphological~rate)))

hex_plot_ap <- 
  ggplot(final_df_ap, aes(x = branch_mass, y = branch_rates_ap)) +
  geom_hex(bins = 50, linewidth = 0.1, show.legend = FALSE) + 
  scale_fill_viridis_c(option = "magma", begin = 0.22, end = 1.0) + 
  #geom_smooth(method = "loess", color = "#00FFCC", se = FALSE, linewidth = 1.2, linetype = "dashed") +
  theme_bw() + theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5), axis.title = element_text(face = "bold", size = 11), panel.grid.major = element_line(color = "grey95"), panel.grid.minor = element_blank()) +
  labs(title = "appendicular characters", x = expression(log[10]~(branch~mass)), y = expression(log[10]~(morphological~rate)))

comp <- ggplot() +
  geom_histogram(data = final_df, aes(x = branch_mass, y = after_stat(density)), fill = "grey90", color = "white", bins = 40, alpha = 1) +
  geom_density(data = top_5, aes(x = branch_mass), fill = "#83D350", color = "#519129", alpha = 0.35, linewidth = 1) +
  theme_bw() + theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5), axis.title = element_text(size = 11), panel.grid.major = element_line(color = "grey95"), panel.grid.minor = element_blank(), plot.margin = margin(5, 5, 5, 15)) +
  labs(title = "all characters rate bursts", x = expression(log[10]~(branch~mass)), y = "density")

comp_cd <- ggplot() +
  geom_histogram(data = final_df_cd, aes(x = branch_mass, y = after_stat(density)), fill = "grey90", color = "white", bins = 40, alpha = 1) +
  geom_density(data = top_5_cd, aes(x = branch_mass), fill = "#83D350", color = "#519129", alpha = 0.35, linewidth = 1) +
  theme_bw() + theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5), axis.title = element_text(size = 11), panel.grid.major = element_line(color = "grey95"), panel.grid.minor = element_blank(), plot.margin = margin(5, 5, 5, 15)) +
  labs(title = "craniodental rate bursts", x = expression(log[10]~(branch~mass)), y = "density")

comp_ax <- ggplot() +
  geom_histogram(data = final_df_ax, aes(x = branch_mass, y = after_stat(density)), fill = "grey90", color = "white", bins = 40, alpha = 1) +
  geom_density(data = top_5_ax, aes(x = branch_mass), fill = "#83D350", color = "#519129", alpha = 0.35, linewidth = 1) +
  theme_bw() + theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5), axis.title = element_text(size = 11), panel.grid.major = element_line(color = "grey95"), panel.grid.minor = element_blank(), plot.margin = margin(5, 5, 5, 15)) +
  labs(title = "axial rate bursts", x = expression(log[10]~(branch~mass)), y = "density")

comp_ap <- ggplot() +
  geom_histogram(data = final_df_ap, aes(x = branch_mass, y = after_stat(density)), fill = "grey90", color = "white", bins = 40, alpha = 1) +
  geom_density(data = top_5_ap, aes(x = branch_mass), fill = "#83D350", color = "#519129", alpha = 0.35, linewidth = 1) +
  theme_bw() + theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5), axis.title = element_text(size = 11), panel.grid.major = element_line(color = "grey95"), panel.grid.minor = element_blank(), plot.margin = margin(5, 5, 5, 15)) +
  labs(title = "appendicular rate bursts", x = expression(log[10]~(branch~mass)), y = "density")

#combine the plots

combined_figure <- 
  hex_plot + hex_plot_cd + hex_plot_ax + hex_plot_ap +  
  comp + comp_cd + comp_ax + comp_ap +  
  plot_layout(
    ncol = 2, 
    byrow = FALSE,            
    guides = "collect"        
  ) + 
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )
  ) & 
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    text = element_text(family = "sans")
  )

ggsave(
  filename = "fig4.pdf", 
  plot = combined_figure,
  width = 261,                            
  height = 353,                           
  units = "mm",
  dpi = 600
)
