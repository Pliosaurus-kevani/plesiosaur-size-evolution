#=============#
#load packages
#=============#

library(cluster) #version 2.1.8.2; loaded for performing PCoA
library(tidytree) #version 0.4.7; loaded for processing trees
library(ape) #version 5.8-1; loaded for processing trees
library(coda) #version 0.19-4.1; loaded for checking convergence
library(phytools) #version 2.5-2; loaded for plotting trees
library(phylolm) #version 2.6.5; loaded for phylogenetic generalized least squares
library(dplyr) #version 1.1.2; loaded for visualization
library(ggplot2) #version 4.0.2; loaded for visualization
library(ggtree) #version 4.0.4; loaded for visualization
library(ggdist) #version 3.3.3; loaded for visualization
library(patchwork) #version 1.3.2; loaded for visualization

#------------------------------------------------------------------------------#

#==========================#
#read data and perform PCoA
#==========================#

#change working directory

#setwd()

#read the neck data

neck.data <- read.csv("neck data.csv", row.names = 1, header = T)

neck <- data.frame(neck.data[,-1], row.names = row.names(neck.data))

#extract and log-transform the continuous variables

cont <- data.frame(neck[,-3], row.names = row.names(neck))

cont <- log10(cont)

#z-transform the continuous variables

cont <- scale(cont)

#combine the data

neck.p <- data.frame(SKL.neck = cont[,1], neck.trunk = cont[,2], CN = neck$cervical, row.names = row.names(neck))

#set the type of variable for the cervical number

neck.p$CN <- factor(neck.p$CN, ordered = T)

#compute the Gower distance

dist <- daisy(neck.p, metric = "gower", stand = F)

#perform principal coordinate analysis

pcoa <- cmdscale(dist, eig = T)

#extract PCo1

pco1 <- data.frame(neck = pcoa$points[,1], row.names = row.names(pcoa$points))

#load the trees

load("tree.RData")

#find the taxa both present in the data set and the trees

name.int <- intersect(row.names(neck.data), tree[[1]]$tip.label)

#prune the trees

tree.p <- ape::keep.tip.multiPhylo(tree, name.int)

#plot to check

plot(ladderize(tree.p[[1]]), no.margin = T)

pco1.p <- setNames(pco1[tree.p[[1]]$tip.label,], tree.p[[1]]$tip.label)

#save the neck length data

#setwd()

write.table(pco1.p, file = "neck.txt", sep = "\t", col.names = F, quote = F)

#save the trees

for(i in 1:100)
{
  write.nexus(tree.p[[i]], file = paste("tree_",i,".nex", sep = ""), translate = T)
}

#------------------------------------------------------------------------------#

#============================#
#create files for BayesTraits
#============================#

#create lists to contain the cmd files

#create cmd files
neck.cmd <- list()

#Use global transformation (Lambda) here
#Use the threaded version (using Cores command to set the number of cores)

for(i in 1:100)
{
  neck.cmd[[i]] <- paste("7\n",
                         "2\n",
                         "VarRates\n",
                         "iterations 12000000\n",
                         "sample 10000\n",
                         "burnin 2000000\n",
                         "Lambda\n",
                         "stones 100 1000\n",
                         "Cores 8\n",
                         "Logfile neck_",i,".log.txt\n", sep="",
                         "Info\n",
                         "Run")
}

#export the cmd files

for(i in 1:100)
{
  write(neck.cmd[[i]], file=paste("neck_script_",i,".cmd", sep=""))
}

#run BayesTraits
#run the following for-loop in cmd

#for /l %i in (1,1,100) do (BayesTraitsV4.exe tree_%i.nex neck.txt < neck_script_%i.cmd)

#------------------------------------------------------------------------------#

#=================#
#check convergence
#=================#

skip_neck_convergence <- list()
neck_convergence <- list()

#read data

for(i in 1:100)
{
  skip_neck_convergence[[i]] <- grep("Tree No", scan(file = paste("neck_",i,".log.txt.Log.txt", sep=""), what="c", quiet=T, sep="\n", blank.lines.skip=FALSE)) - 1
  neck_convergence[[i]] = read.table(paste("neck_",i,".log.txt.Log.txt", sep=""), skip = skip_neck_convergence[[i]], sep = "\t",  quote="\"", header = TRUE)
  neck_convergence[[i]] = neck_convergence[[i]][,-ncol(neck_convergence[[i]])]
}

#convert the results to coda format
res_neck <- list()

for(i in 1:100)
{
  res_neck[[i]] <- mcmc(subset(neck_convergence[[i]], select=-c(Iteration, Tree.No)),
                            start=min(neck_convergence[[i]]$Iteration),
                            end=max(neck_convergence[[i]]$Iteration),thin=10000)
}

#get effective size(should be >200)
ess_neck_list <- list()

for(i in 1:100)
{
  ess_neck_list[[i]] <- effectiveSize(res_neck[[i]])
}

#find minimum effective size

tmp_ess_neck_min <- do.call(rbind, ess_neck_list)
ess_neck_min <- min(tmp_ess_neck_min)
ess_neck_min

#------------------------------------------------------------------------------#

#===================#
#posterior processor
#===================#

#run the following code in cmd:

#for /l %i in (1,1,100) do (PPPostProcess.exe neck_%i.log.txt.VarRates.txt > PP_neck_%i.txt)

#------------------------------------------------------------------------------#

#==========================#
#summarize the branch rates
#==========================#

#read the results
bayes.postproc.neck <- list()

#As we checked above, the effective sample size of all the cases are larger than 200.
#Thus we don't need to prune any case.

for(i in 1:100)
{
  bayes.postproc.neck[[i]] <- read.delim(paste('PP_neck_',i,'.txt', sep = ""), header = TRUE, check.names = FALSE, sep = "\t", colClasses = c(rep(NA, 22), rep("NULL", 1)))
}

for(i in 1:100)
{
  #convert the spaces in the column names to underscores "_"
  colnames(bayes.postproc.neck[[i]]) <- gsub(" ", "_", colnames(bayes.postproc.neck[[i]])) 
  
  #Split reported strings/names, so that we can use them with getMRCA
  bayes.postproc.neck[[i]]$Taxa_List <- sapply(strsplit(as.character(bayes.postproc.neck[[i]]$Taxa_List), ","), "[") 
  
  #create new ID columns
  bayes.postproc.neck[[i]]["NEW_Edge_ID"] <- NA 
  bayes.postproc.neck[[i]]["NEW_Node_ID"] <- NA 
}

#Find new node and edge IDs of our data from the trees

for(i in 1:100)
{
  for(j in 1:length(bayes.postproc.neck[[i]]$Taxa_List))
  {
    if(length(bayes.postproc.neck[[i]]$Taxa_List[[j]]) == 1){
      #Single tips: find tip nodes and edges
      bayes.postproc.neck[[i]]$NEW_Node_ID[[j]]  <- which(tree.p[[i]]$tip.label == bayes.postproc.neck[[i]]$Taxa_List[[j]]) #find tip nodes
      bayes.postproc.neck[[i]]$NEW_Edge_ID[[j]] <- which.edge(tree.p[[i]], bayes.postproc.neck[[i]]$Taxa_List[[j]]) #find tip edges
    }  else{
      #Multiple tips: find internal nodes and edges of MRCA
      #Notice that the root edge is not assigned an edge ID, i.e. it stays NA
      bayes.postproc.neck[[i]]$NEW_Node_ID[[j]]  <- getMRCA(tree.p[[i]], bayes.postproc.neck[[i]]$Taxa_List[[j]]) #find internal node
      bayes.postproc.neck[[i]]$NEW_Edge_ID[[j]] <- which.edge(tree.p[[i]], getMRCA(tree.p[[i]], bayes.postproc.neck[[i]]$Taxa_List[[j]])) #find internal edges
    }
  }
}

#re-order the data frame to match the trees
#since the edge ID the root is placed on the first row in the data frame

for(i in 1:100)
{
  bayes.postproc.neck[[i]]$NEW_Edge_ID[is.na(bayes.postproc.neck[[i]]$NEW_Edge_ID)] <- length(bayes.postproc.neck[[i]]$NEW_Edge_ID)
  
  #Replace rownames with NEW_Edge_ID
  rownames(bayes.postproc.neck[[i]]) <- bayes.postproc.neck[[i]]$NEW_Edge_ID 
  
  #Order rows according to row numbers
  bayes.postproc.neck[[i]] <- bayes.postproc.neck[[i]][ order(as.numeric(row.names(bayes.postproc.neck[[i]]))),] 
}

myvarpar <- "Mean_Scalar"

#plot to check
for(i in 51:60)
{
  plotBranchbyTrait(tree.p[[i]],bayes.postproc.neck[[i]][[myvarpar]], mode="edges", type="phylogram", show.tip.label=TRUE, show.node.label=TRUE, cex=0.5, palette = colorRampPalette(c("#2C1B47", "#724C9D", "#EDC5AB", "#FF6933", "#F90627")))
}

#------------------------------------------------------------------------------#

#================================#
#check correlation with body mass
#================================#

#change working directory

#setwd()

#read body mass data

mass.data <- read.csv("mass.csv", header = T, row.names = 1)

mass.data <- mass.data[,-1]

#log-transform the data

mass.data$mass <- log10(mass.data$mass)

#prune the trees

name.int <- intersect(tree.p[[1]]$tip.label, row.names(mass.data))

tree.p2 <- keep.tip.multiPhylo(tree.p, name.int)

#combine the mass and neck rate data
#then perform phylogenetic generalized least squares

PGLS <- list()

#Here I employ a temporary variable called "moa" as counter

moa <- 0

for(i in 1:100)
{
  #here I employ some temporary variables: the dodo family
  #and auk
  
  #extract the body masses of the species within the tree
  
  dodo.mass <- mass.data[tree.p2[[i]]$tip.label,]
  
  #extract the neck rates of the species within the tree
  
  auk <- setNames(bayes.postproc.neck[[i]][[myvarpar]], bayes.postproc.neck[[i]]$Taxa_List)
  dodo.rate <- data.frame(auk[tree.p2[[i]]$tip.label], row.names = tree.p2[[i]]$tip.label)  
  colnames(dodo.rate) <- c("neck_rate")
  
  #log-transform the neck rates
  dodo.rate$neck_rate <- log10(dodo.rate$neck_rate)
  
  #combine the datasets
  dodo <- cbind(dodo.mass, dodo.rate)
  
  #perform robust phylogenetic regression
  
  PGLS[[i]] <- phylolm(formula = neck_rate ~ mass,
                      data = dodo[tree.p2[[i]]$tip.label,],
                      phy = tree.p2[[i]],
                      model = "lambda")
  
  #count the number of cases which fail to reject the null hypothesis
  
  if(summary(PGLS[[i]])$coefficients[2,4] > 0.05)
  {
    moa <- moa + 1
  }
}

#------------------------------------------------------------------------------#

#======================================#
#plot (neck rates mapped onto branches)
#======================================#

#select a tree randomly

my_tree <- tree.p[[20]]
my_data <- bayes.postproc.neck[[20]]

#summarize the data

plot_data <- data.frame(
  node = my_data$NEW_Node_ID,
  rate = my_data[[myvarpar]] 
)

plot_data$rate[103] <- NA

#correct the species name and
#remove the dash in species names

my_tree$tip.label <- gsub("_", " ", my_tree$tip.label)


p_bayes <- 
  ggtree(my_tree, ladderize = T, right = T, color = "grey20", linewidth = 2.5) %<+% plot_data +
  geom_tree(aes(color = rate), linewidth = 1.2) +
  geom_tiplab(fontface = "italic", size = 3, offset = max(nodeHeights(my_tree))*0.01) +
  scale_color_gradientn(
    name = "neck rate", 
    colors = c("#73A9AD", "#90C8AC", "#FEEB7D", "#FF8787", "#F03E3E"),
    values = scales::rescale(c(
      1.5,                            
      3,                          
      4.5,                          
      6,                          
      max(plot_data$rate, na.rm = TRUE) 
    )),
    na.value = "grey50"
  )+
  xlim(0, max(nodeHeights(my_tree))*1.3) +
  theme(
    legend.position = c(0.1, 0.2), 
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(1, "cm"),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(0, 0, 20, 0)
  )

p_bayes

#save the plot

#setwd()

saveRDS(p_bayes, file = "neck.rate.rds")

#------------------------------------------------------------------------------#

#=====================#
#plot: raincloud plot
#=====================#

Cryptoclididae_tips <- c("Muraenosaurus_leedsii", "Tricleidus_seeleyi", 
                         "Cryptoclidus_eurymerus", "Ophthalmothule_cryostea",
                         "Colymbosaurus_megadeirus", "Abyssosaurus_nataliae")  

Thalassophonea_tips <- c("Peloneustes_philarchus", "Simolestes_vorax",
                         "Liopleurodon_ferox", "Luskhan_itilensis",
                         "Brachauchenius_lucasi", "Kronosaurus_MCZ_1285",
                         "Sachicasaurus_vitae", "Stenorhynchosaurus_munozi")

Elasmosauridae_tips <- c("Callawayasaurus_colombiensis", "Libonectes_morgani",
                         "Albertonectes_vanderveldei", "Hydrotherosaurus_alexandrae",
                         "Vegasaurus_molyi", "Aristonectes_quiriquinensis",
                         "Kaiwhekea_katiki", "Nakonanectes_bradti",
                         "Thalassomedon_haningtoni", "Styxosaurus_rezaci",
                         "Styxosaurus_SDSM_451", "Styxosaurus_browni")


#create a data frame to extract the data

all_rates_df <- data.frame()


for(i in 1:100) {
  current_data <- bayes.postproc.neck[[i]]
  
  clade_assignment <- rep("Remaining taxa", nrow(current_data))
  
  for(j in 1:nrow(current_data)) {
    branch_taxa <- current_data$Taxa_List[[j]]
    
    if(all(branch_taxa %in% Cryptoclididae_tips)) {
      clade_assignment[j] <- "Cryptoclididae"
    } else if(all(branch_taxa %in% Thalassophonea_tips)) {
      clade_assignment[j] <- "Thalassophonea"
    } else if(all(branch_taxa %in% Elasmosauridae_tips)) {
      clade_assignment[j] <- "Elasmosauridae"
    }
  }

  temp_df <- data.frame(
    tree_id = i,
    rate = current_data[[myvarpar]],
    clade = clade_assignment
  )
  
  all_rates_df <- rbind(all_rates_df, temp_df)
}

all_rates_clean <- all_rates_df %>%
  filter(!is.na(rate) & rate > 0) %>%
  mutate(log_rate = log10(rate))

#order the subsets

clade_levels <- c("Cryptoclididae", "Thalassophonea", "Elasmosauridae", "Remaining taxa")
all_rates_clean$clade <- factor(all_rates_clean$clade, levels = rev(clade_levels))

#define the colors

my_colors <- c("Cryptoclididae" = "#D2A16E",   
               "Thalassophonea" = "#5061c5",   
               "Elasmosauridae" = "#F34C3F",   
               "Remaining taxa" = "#8498AB")   

p_raincloud <- 
  ggplot(all_rates_clean, aes(x = log_rate, y = clade, fill = clade)) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, scale = 0.5, justification = -0.35, point_colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5, color = "grey30") +
  geom_jitter(data = all_rates_clean[sample(nrow(all_rates_clean), min(5000, nrow(all_rates_clean))), ], 
              aes(color = clade), width = 0, height = 0.1, alpha = 0.1, size = 0.5) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  coord_cartesian(ylim = c(1.2, 4.1)) +
  theme_bw() + 
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12)) +
  labs(title = "neck evolution rate variations across clades", 
       x = expression(log[10]~(branch~rate)))

print(p_raincloud)

saveRDS(p_raincloud, file = "neck.raincloud.rds")

#------------------------------------------------------------------------------#

#======================#
#combine the subfigures
#======================#

scatter.plot <- readRDS("mass_neck_scatter.rds")
bin.plot <- readRDS("neck_binplot.rds")
mass.rain.plot <- readRDS("mass.raincloud.rds")

all_mass_clean <- mass.rain.plot$data

mass.rain.plot <- 
  ggplot(all_mass_clean, aes(x = branch_mass, y = clade, fill = clade)) +
  
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, scale = 0.5, justification = -0.35, point_colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5, color = "grey30") +
  geom_jitter(data = all_mass_clean[sample(nrow(all_mass_clean), min(5000, nrow(all_mass_clean))), ], 
              aes(color = clade), width = 0, height = 0.1, alpha = 0.1, size = 0.5) +
  
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  coord_cartesian(ylim = c(1.2, 4.1)) +
  theme_bw() + 
  theme(legend.position = "none", 
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12)) +
  labs(title = "body mass variations across clades", 
       x = expression(log[10]~(branch~mass)))


layout_matrix <- "
AABBBB
CCCCCC
CCCCCC
CCCCCC
DDDEEE
"

final_figure <- 
  scatter.plot + bin.plot + p_bayes + p_raincloud + mass.rain.plot +
  plot_layout(
    design = layout_matrix
  ) +
  plot_annotation(tag_levels = 'A') & 
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA)
  )

final_figure

ggsave(
  filename = "fig3.pdf", 
  plot = final_figure, 
  width = 10, 
  height = 13, 
  dpi = 300,
  device = cairo_pdf 
)
