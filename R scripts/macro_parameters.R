#=============#
#load packages
#=============#

library(RevGadgets) #version 1.2.1; loaded for reading trace samples
library(ggplot2) #version 4.0.2; loaded for visualization

#------------------------------------------------------------------------------#

#=========================#
#read and process the data
#=========================#

#change the working pathway

#setwd()

#read the trace file

trace <- readTrace("results_run_1.log", burnin = 0.1)

#read the speciation rate

lambda_list <- list()

for(bin_num in 1:5) {
  
  lambda_index <- 6 - bin_num 
  
  col_name <- paste0("lambda[", lambda_index, "]")
  
  rates <- trace[[1]][[col_name]]
  
  bins <- rep(paste("bin", bin_num), 900)
  
  lambda_list[[bin_num]] <- data.frame(Rate = rates, Bin = bins)
}

lambda.d <- do.call(rbind, lambda_list)

lambda.d$Bin <- factor(lambda.d$Bin, levels = paste("bin", 1:5))

#read the extinction rate

mu_list <- list()

for(bin_num in 1:5) {
  
  mu_index <- 6 - bin_num 
  
  col_name <- paste0("mu[", mu_index, "]")
  
  rates <- trace[[1]][[col_name]]
  
  bins <- rep(paste("bin", bin_num), 900)
  
  mu_list[[bin_num]] <- data.frame(Rate = rates, Bin = bins)
}

mu.d <- do.call(rbind, mu_list)

mu.d$Bin <- factor(mu.d$Bin, levels = paste("bin", 1:5))

#read the net diversification rate

div_list <- list()

for(bin_num in 1:5) {
  div_index <- 6 - bin_num
  
  col_name <- paste0("div[", div_index, "]")
  
  rates <- trace[[1]][[col_name]]
  
  bins <- rep(paste("bin", bin_num), 900)
  
  div_list[[bin_num]] <- data.frame(Rate = rates, Bin = bins)
}

div.d <- do.call(rbind, div_list)

div.d$Bin <- factor(div.d$Bin, levels = paste("bin", 1:5))

#read the fossil sampling rate

psi_list <- list()

for(bin_num in 1:5) {
  psi_index <- 6 - bin_num
  
  col_name <- paste0("psi[", psi_index, "]")
  
  rates <- trace[[1]][[col_name]]
  
  bins <- rep(paste("bin", bin_num), 900)
  
  psi_list[[bin_num]] <- data.frame(Rate = rates, Bin = bins)
}

psi.d <- do.call(rbind, psi_list)

psi.d$Bin <- factor(psi.d$Bin, levels = paste("bin", 1:5))

#------------------------------------------------------------------------------#

#====#
#plot
#====#

p_lambda <- ggplot(lambda.d, aes(x = Bin, y = Rate, fill = Bin)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  theme_bw() +
  labs(
    x = "Time Bin",
    y = "Speciation Rate"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.major.x = element_blank()
  )

p_lambda

p_mu <- ggplot(mu.d, aes(x = Bin, y = Rate, fill = Bin)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  theme_bw() +
  labs(
    x = "Time Bin",
    y = "Extinction Rate"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.major.x = element_blank()
  )

p_mu

p_div <- ggplot(div.d, aes(x = Bin, y = Rate, fill = Bin)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  theme_bw() +
  labs(
    x = "Time Bin",
    y = "Net Diversification Rate"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.major.x = element_blank()
  )

p_div

p_psi <- ggplot(psi.d, aes(x = Bin, y = Rate, fill = Bin)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  theme_bw() +
  labs(
    x = "Time Bin",
    y = "Fossil Sampling Rate"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.major.x = element_blank()
  )

p_psi

#combine the plots

lambda.d$Parameter <- "speciation"
mu.d$Parameter     <- "extinction"
div.d$Parameter    <- "net diversification"
psi.d$Parameter    <- "fossil sampling"

all_rates <- rbind(lambda.d, mu.d, div.d, psi.d)

all_rates$Parameter <- factor(all_rates$Parameter, 
                              levels = c("speciation", 
                                         "extinction", 
                                         "net diversification", 
                                         "fossil sampling"))

p_all <- ggplot(all_rates, aes(x = Bin, y = Rate, fill = Bin)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  scale_fill_viridis_d(option = "mako", end = 0.8) + 
  facet_wrap(~ Parameter, scales = "free", ncol = 2) + 
  theme_bw() +
  labs(
    x = "Time Bin",
    y = "Rate (Events per Million Years)"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25, face = "bold"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    strip.text = element_text(size = 30, face = "bold"), 
    strip.background = element_rect(fill = "gray95"),    
    panel.grid.major.x = element_blank()
  )

p_all

#save the plot

ggsave("evo_parameters.pdf", p_all, width = 8, height = 6)
