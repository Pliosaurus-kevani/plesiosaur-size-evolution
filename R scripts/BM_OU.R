#=============#
#load packages
#=============#

library(tidytree) #version 0.4.7; loaded for loading trees
library(ape) #version 5.8-1; loaded for processing trees
library(RevGadgets) #version 1.2.1; loaded for reading trace samples
library(coda) #version 0.19-4.1; loaded for checking convergence
library(ggplot2) #version 4.0.2; loaded for visualization
library(scales) #version 1.4.0; loaded for visualization

#------------------------------------------------------------------------------#

#=====================#
#read and process data
#=====================#

#change the working directory

#setwd()

#load the trees

load("tree.RData")

#read the body mass data

data <- read.csv("mass.csv", header = T, row.names = 1)

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

#split the data sets 

Elasmosauridae <- data[data$clade == "Elasmosauridae",]
Thalassophonea <- data[data$clade == "Thalassophonea",]
others <- data[data$clade == "others",]

Elas.mass <- data.frame(taxon = row.names(Elasmosauridae), mass = Elasmosauridae$mass)
Tha.mass <- data.frame(taxon = row.names(Thalassophonea), mass = Thalassophonea$mass)
others.mass <- data.frame(taxon = row.names(others), mass = others$mass)

#prune the trees

tree.Elas <- ape::keep.tip.multiPhylo(tree, row.names(Elasmosauridae))
tree.Tha <- ape::keep.tip.multiPhylo(tree, row.names(Thalassophonea))
tree.others <- ape::keep.tip.multiPhylo(tree, row.names(others))

#plot to check the pruned trees

plot(ladderize(tree.Elas[[1]]))
plot(ladderize(tree.Tha[[1]]))
plot(ladderize(tree.others[[1]]))

#change working directory and save the data

#setwd()

write.nexus(tree.Elas, file = "trees_Elas.nex", translate = T)
write.table(Elas.mass, file = "mass_Elas.txt", sep = "\t", row.names = F, quote = F)

#setwd()

write.nexus(tree.Tha, file = "trees_Tha.nex", translate = T)
write.table(Tha.mass, file = "mass_Tha.txt", sep = "\t", row.names = F, quote = F)

#setwd()

write.nexus(tree.others, file = "trees_others.nex", translate = T)
write.table(others.mass, file = "mass_others.txt", sep = "\t", row.names = F, quote = F)

#------------------------------------------------------------------------------#

#=================#
#check convergence
#=================#

#create lists to save the samples

Elas.s <- list()
Tha.s <- list()
others.s <- list()

#read the samples

#change the working directory to the folders 
#where you store the RevBayes outputs 

#setwd()

for(i in 1:100)
{
  Elas.s[[i]] <- readTrace(path = paste("simple_OU_RJ_", i, ".log", sep = ""), burnin = 0.25)
}

#setwd()

for(i in 1:100)
{
  Tha.s[[i]] <- readTrace(path = paste("simple_OU_RJ_", i, ".log", sep = ""), burnin = 0.25)
}

#setwd()

for(i in 1:100)
{
  others.s[[i]] <- readTrace(path = paste("simple_OU_RJ_", i, ".log", sep = ""), burnin = 0.25)
}

#check convergence: Elasmosauridae

for(i in 1:100)
{
  for(j in 2:9)
  {
    if(effectiveSize(Elas.s[[i]][[1]][[j]]) < 200)
    {
      print(paste("Elasmosauridae: case", i, "fails to reach convergence"))
    }
  }
  #if nothing is returned, then all cases reach convergence
}

#check convergence: Thalassophonea

for(i in 1:100)
{
  for(j in 2:9)
  {
    if(effectiveSize(Tha.s[[i]][[1]][[j]]) < 200)
    {
      print(paste("Thalassophonea: case", i, "fails to reach convergence"))
    }
  }
  #if nothing is returned, then all cases reach convergence
}

#check convergence: other plesiosaurs

for(i in 1:100)
{
  for(j in 2:9)
  {
    if(effectiveSize(others.s[[i]][[1]][[j]]) < 200)
    {
      print(paste("others: case", i, "fails to reach convergence"))
    }
  }
  #if nothing is returned, then all cases reach convergence
}

#------------------------------------------------------------------------------#

#==================#
#combine the traces
#==================#

#initialize

Elas.s.com <- data.frame(Elas.s[[1]][[1]][,5:9], row.names = NULL)
Tha.s.com <- data.frame(Tha.s[[1]][[1]][,5:9], row.names = NULL)
others.s.com <- data.frame(others.s[[1]][[1]][,5:9], row.names = NULL)

for(i in 2:100)
{
  #Here I employ some temporary variables: the auk family
  
  auk.Elas <- data.frame(Elas.s[[i]][[1]][,5:9], row.names = NULL)
  auk.Tha <- data.frame(Tha.s[[i]][[1]][,5:9], row.names = NULL)
  auk.others <- data.frame(others.s[[i]][[1]][,5:9], row.names = NULL)
  
  Elas.s.com <- rbind(Elas.s.com, auk.Elas)
  Tha.s.com <- rbind(Tha.s.com, auk.Tha)
  others.s.com <- rbind(others.s.com, auk.others)
}

#------------------------------------------------------------------------------#

#=====================================#
#mass distribution of the three groups
#=====================================#

plot.mass <- data.frame(clade = data$clade, mass = 10^(data$mass), row.names = row.names(data))

plot.mass$clade <- factor(plot.mass$clade, level = c("Elasmosauridae", "Thalassophonea", "others"))

mass.distribution <- 
  ggplot(plot.mass, aes(x = mass, fill = clade))+
  geom_density(alpha = 0.6, color = "black",          
               linewidth = 0.3,               
               adjust = 1.2 )+
  scale_fill_manual(
    values = c(
      "Elasmosauridae"   = "#F34C3FCC",  
      "Thalassophonea"    = "#5061c5CC",  
      "others" = "#8498ABCC"
    )
  )+
  scale_x_continuous(
    name = "body mass (kg)", 
    expand = expansion(mult = c(0.05, 0.05))
  )+
  ggtitle("body mass distribution")+
  theme_bw(base_size = 7)+
  theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 12),         
        axis.text = element_text(size = 10),          
        legend.text = element_text(size = 10),        
        legend.title = element_text(size = 12))

mass.distribution

#=================================#
#posterior probability of BM vs OU
#=================================#

#summarize the posterior probability of the models 

isOU <- data.frame(matrix(data = NA, nrow = 6, ncol = 3))

isOU[,1] <- c("Elasmosauridae", "Elasmosauridae", "Thalassophonea", "Thalassophonea", "others", "others")
isOU[,2] <- c("Elasmosauridae (BM)", "Elasmosauridae (OU)", "Thalassophonea (BM)", "Thalassophonea (OU)", "others (BM)", "others (OU)")
isOU[1,3] <- sum(Elas.s.com$is_BM)/375000
isOU[2,3] <- sum(Elas.s.com$is_OU)/375000
isOU[3,3] <- sum(Tha.s.com$is_BM)/375000
isOU[4,3] <- sum(Tha.s.com$is_OU)/375000
isOU[5,3] <- sum(others.s.com$is_BM)/375000
isOU[6,3] <- sum(others.s.com$is_OU)/375000

colnames(isOU) <- c("clade", "model", "probability")

isOU$model <- factor(isOU$model, level = c("others (OU)", "others (BM)", "Thalassophonea (OU)", "Thalassophonea (BM)", "Elasmosauridae (OU)", "Elasmosauridae (BM)"))

BF <- c(1, 3, 10)
prob = BF/(1+BF)

bar.plot <- 
  ggplot(isOU, aes(x = model, y = probability, fill = clade)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.6) +
  scale_fill_manual(
    values = c(
      "Elasmosauridae"   = "#F34C3FCC",  
      "Thalassophonea"    = "#5061c5CC",  
      "others" = "#8498ABCC"   
    )) +
  labs(
    title = "BM vs OU",
    y = "posterior probability",
    x = NULL
  ) +
  theme_bw(base_size = 7)+
  theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 12),         
        axis.text = element_text(size = 10),          
        legend.position = "none")+
  geom_hline(yintercept=prob, linetype="dashed", color = "grey", linewidth = 0.5)+
  coord_flip()

bar.plot

#=====================#
#distribution of theta
#=====================#

#create a data frame containing the theta values

thetas <- data.frame(matrix(data = NA, nrow = 1125000, ncol = 2))
colnames(thetas) <- c("clade", "theta")

thetas[1:375000,1] <- "Elasmosauridae"
thetas[1:375000,2] <- Elas.s.com$theta
thetas[375001:750000,1] <- "Thalassophonea"
thetas[375001:750000,2] <- Tha.s.com$theta
thetas[750001:1125000,1] <- "others"
thetas[750001:1125000,2] <- others.s.com$theta

thetas$clade <- factor(thetas$clade, level = c("Elasmosauridae", "Thalassophonea", "others"))

theta.plot <-
  ggplot(thetas, aes(x = theta, fill = clade))+
  geom_density(alpha = 0.6, color = "black", linewidth =0.3, adjust = 1.2)+
  scale_fill_manual(
    values = c(
      "Elasmosauridae"   = "#F34C3FCC",  
      "Thalassophonea"    = "#5061c5CC",  
      "others" = "#8498ABCC"
    )
  )+
  scale_x_continuous(
    name = expression(theta~"(log"["10"]*" mass)"), 
    expand = expansion(mult = c(0.05, 0.05))
  )+
  ggtitle("selective optima")+
  theme_bw(base_size = 7)+
  theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 12),         
      axis.text = element_text(size = 10),          
      legend.position = "none")

theta.plot

#=====================#
#distribution of alpha
#=====================#

#create a data frame containing the theta values

alphas <- data.frame(matrix(data = NA, nrow = 1125000, ncol = 2))
colnames(alphas) <- c("clade", "alpha")

alphas[1:375000,1] <- "Elasmosauridae"
alphas[1:375000,2] <- Elas.s.com$alpha
alphas[375001:750000,1] <- "Thalassophonea"
alphas[375001:750000,2] <- Tha.s.com$alpha
alphas[750001:1125000,1] <- "others"
alphas[750001:1125000,2] <- others.s.com$alpha

alphas$clade <- factor(alphas$clade, level = c("others", "Thalassophonea", "Elasmosauridae"))

alpha.plot <- 
  ggplot(alphas, aes(x = clade, y = alpha, fill = clade))+
  geom_violin(
    alpha = 0.6,          
    color = "black",      
    linewidth = 0.3,           
    scale = "width",      
    trim = FALSE,         
    adjust = 1.2          
  ) +
  geom_boxplot(
    width = 0.1,         
    fill = "white",       
    color = "black",      
    alpha = 0.8,          
    outlier.alpha = 0,    
    outlier.size = 1.5,
    linewidth = 0.3             
  ) +
  scale_fill_manual(
    name = "clade",
    values = c(
      "Elasmosauridae" = "#F34C3FCC",  
      "Thalassophonea" = "#5061c5CC",  
      "others"         = "#8498ABCC"
    )
  ) +
  scale_y_continuous(   
    name = expression(alpha), 
    expand = expansion(mult = c(0.05, 0.1)), 
    breaks = waiver()   
  ) +
  scale_x_discrete(     
    name = NULL
  ) +
  labs(title = "selection strength") +
  theme_bw(base_size = 7) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(size = 12),         
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )+
  coord_flip()

alpha.plot

#save the plots

#setwd()

saveRDS(mass.distribution, file = "mass.distribution.rds")
saveRDS(alpha.plot, file = "alpha.plot.rds")
saveRDS(theta.plot, file = "theta.plot.rds")
saveRDS(bar.plot, file = "bar.plot.rds")
