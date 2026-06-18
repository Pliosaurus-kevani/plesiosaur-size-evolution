#=============#
#load packages
#=============#

library(tidytree) #version 0.4.7; loaded for loading trees
library(phytools) #version 2.5-2; loaded for stochastic character mapping
library(geiger) #version 2.0.11; loaded for fitting single-regime models
library(mvMORPH) #version 1.2.1; loaded for fitting multi-regime models
library(ggplot2) #version 4.0.2; loaded for visualization

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

data <- data[,-1]

#save the data in another way

data.v <- setNames(data$mass, row.names(data))

#------------------------------------------------------------------------------#

#============================#
#stochastic character mapping
#============================#

#prune the trees

tree.p <- keep.tip.multiPhylo(tree, row.names(data))

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
  tree.m[[i]] <- make.simmap(tree.p[[i]], clade.info, model = mkmodel)
  
  #report progress
  paste("working on tree", i, sep = " ")
}

#plot to check

cols3 <- setNames(c("#F34C3F","#8498AB","#5061c5"), levels(clade.info))

plot(ladderize.simmap(tree.m[[1]]), cols3, fsize = 0.4, ftype = "i", lwd = 2, offset = 0.4, ylim = c(-1, Ntip(tree.m[[1]])))
plot(ladderize.simmap(tree.m[[50]]), cols3, fsize = 0.4, ftype = "i", lwd = 2, offset = 0.4, ylim = c(-1, Ntip(tree.m[[1]])))
plot(ladderize.simmap(tree.m[[100]]), cols3, fsize = 0.4, ftype = "i", lwd = 2, offset = 0.4, ylim = c(-1, Ntip(tree.m[[1]])))

#------------------------------------------------------------------------------#

#===============================#
#fit multi-regime models to data
#===============================#

#create lists to save the results

BM1 <- list()
OU1 <- list()
BMs <- list()
OUs <- list()
EB <- list()
trend <- list()
drift <- list()

#fit the models to data

for(i in 1:100)
{
  BMs[[i]] <- mvBM(tree.m[[i]], data[tree.m[[i]]$tip.label,2], model = "BMM", method = "pseudoinverse", scale.height = F, echo = F)
  OUs[[i]] <- mvOU(tree.m[[i]], data[tree.m[[i]]$tip.label,2], model = "OUM", method = "pseudoinverse", param = list(vcv = "fixedRoot", root = "stationary"), scale.height = F, echo = F)
  BM1[[i]] <- mvBM(tree.m[[i]], data[tree.m[[i]]$tip.label,2], model = "BM1", method = "pseudoinverse", scale.height = F, echo = F)
  OU1[[i]] <- mvOU(tree.m[[i]], data[tree.m[[i]]$tip.label,2], model = "OU1", method = "pseudoinverse", param = list(vcv = "fixedRoot", root = "stationary"), scale.height = F, echo = F)
  EB[[i]] <- fitContinuous(tree.p[[i]], data.v[tree.p[[i]]$tip.label], model = "EB")
  trend[[i]] <- fitContinuous(tree.p[[i]], data.v[tree.p[[i]]$tip.label], model = "rate_trend")
  drift[[i]] <- fitContinuous(tree.p[[i]], data.v[tree.p[[i]]$tip.label], model = "mean_trend")
  
  #report progress
  print(paste("working on tree", i, sep = " "))
}

#In 4 cases, the parameter a of the early burst model appears at bounds
#This suggests that the EB models reduced to BM model

#summarize the results into a data frame

multi.AICc <- data.frame(matrix(ncol = 7, nrow = 100))
colnames(multi.AICc) <- c("BM1", "OU1", "EB", 
                          "rate_trend", "mean_trend",  
                          "BM3", "OU3")

for(i in 1:100)
{
  multi.AICc[i,1] <- BM1[[i]]$AICc
  multi.AICc[i,2] <- OU1[[i]]$AICc
  multi.AICc[i,3] <- EB[[i]]$opt$aicc
  multi.AICc[i,4] <- trend[[i]]$opt$aicc
  multi.AICc[i,5] <- drift[[i]]$opt$aicc
  multi.AICc[i,6] <- BMs[[i]]$AICc
  multi.AICc[i,7] <- OUs[[i]]$AICc
}

#prune the cases where convergence is not reached

for(i in 1:100)
{
  if(BM1[[i]]$convergence != 0 |
     BMs[[i]]$convergence != 0 |
     OU1[[i]]$convergence != 0 |
     OUs[[i]]$convergence != 0)
  {
    multi.AICc[i,] <- NA
  }
}

#remove these failed cases

multi.AICc <- na.omit(multi.AICc)

#------------------------------------------------------------------------------#

#====#
#plot
#====#

#change the way how data is saved

plot.data <- data.frame(matrix(ncol = 2, nrow = nrow(multi.AICc)*7))
colnames(plot.data) <- c("AICc", "model")

for(i in 1:nrow(multi.AICc))
{
  plot.data[i,1] <- multi.AICc[i,1]
  plot.data[i,2] <- "BM1"
  
  plot.data[i+nrow(multi.AICc),1] <- multi.AICc[i,2]
  plot.data[i+nrow(multi.AICc),2] <- "OU1"
  
  plot.data[i+2*nrow(multi.AICc),1] <- multi.AICc[i,3]
  plot.data[i+2*nrow(multi.AICc),2] <- "EB"
  
  plot.data[i+3*nrow(multi.AICc),1] <- multi.AICc[i,4]
  plot.data[i+3*nrow(multi.AICc),2] <- "rate trend"
  
  plot.data[i+4*nrow(multi.AICc),1] <- multi.AICc[i,5]
  plot.data[i+4*nrow(multi.AICc),2] <- "mean trend"
  
  plot.data[i+5*nrow(multi.AICc),1] <- multi.AICc[i,6]
  plot.data[i+5*nrow(multi.AICc),2] <- "BM3"
  
  plot.data[i+6*nrow(multi.AICc),1] <- multi.AICc[i,7]
  plot.data[i+6*nrow(multi.AICc),2] <- "OU3"
}

#order the models

plot.data$model <- factor(plot.data$model, levels = c("BM1", "OU1", "EB", "rate trend", "mean trend", "BM3", "OU3"))


aicc.plot <- ggplot(plot.data, aes(
  x = model, 
  y = AICc, 
  fill = model
)) +
  geom_boxplot(
    width = 0.5,             
    color = "black",         
    linewidth = 0.3,         
    alpha = 0.8,             
    outlier.size = 1.5,      
    outlier.shape = 21,      
    outlier.fill = NA,  
    outlier.color = NA, 
    outlier.stroke = 0.3     
  ) +
  labs(
    title = "model comparison (AICc)",
    y = "AICc score", 
    x = NULL                 
  ) +
  theme_bw(base_size = 7) +
  theme(
    plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 8),         
    axis.text = element_text(size = 6),          
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "none"
  )

aicc.plot

#save the plot

saveRDS(aicc.plot, "aicc.plot.rds")

