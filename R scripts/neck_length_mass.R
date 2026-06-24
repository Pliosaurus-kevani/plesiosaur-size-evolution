#=============#
#load packages
#=============#

library(cluster) #version 2.1.8.2; loaded for performing cluster analysis
library(fpc) #version 2.2-14; loaded for performing cluster analysis
library(tidytree) #version 0.4.7; loaded for loading trees
library(ape) #version 5.8-1; loaded for processing trees
library(phylolm) #version 2.6.5; loaded for phylogenetic generalize least squares
library(dispRity) #version 1.9; loaded for bootstrap of disparity computation
library(ggplot2) #version 4.0.2; loaded for visualization
library(dplyr) #version 1.1.2; loaded for visualization
library(ggtree) #version 4.0.4; loaded for visualization
library(patchwork) #verison 1.3.2; loaded for visualization
library(ggdist) #version 3.3.3; loaded for visualization

#------------------------------------------------------------------------------#

#=====================#
#read and process data
#=====================#

#change working pathway

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

#------------------------------------------------------------------------------#

#========================#
#perform cluster analysis
#========================#

#compute the Gower distance

dist <- daisy(neck.p, metric = "gower", stand = F)

set.seed(123)

#perform cluster analysis

cboot <- clusterboot(dist,
                     B = 10000,  
                     distances = T,
                     clustermethod = hclustCBI,
                     method = "ward.D2",
                     k = 2) 

#check level of support
#values >0.85 indicate very strong support

cboot$bootmean

#------------------------------------------------------------------------------#

#=====================================#
#perform principal coordinate analysis
#=====================================#

pcoa <- cmdscale(dist, eig = T)

#check the amount of variance explained by PCo1

total_variance <- sum(pcoa$eig)
var_explained_gower <- pcoa$eig[1] / total_variance

#extract PCo1

pco1 <- data.frame(neck = pcoa$points[,1], row.names = row.names(pcoa$points))

#------------------------------------------------------------------------------#

#read the mass data

mass.data <- read.csv("mass.csv", header = T, row.names = 1)

mass <- data.frame(mass = mass.data$mass, row.names = row.names(mass.data))

#log-transform the mass

mass$mass <- log10(mass$mass)

#combine the neck PCo1 and mass data

name.int <- base::intersect(row.names(mass), row.names(pco1))

mass.neck <- cbind(data.frame(mass[name.int,], row.names = name.int), data.frame(pco1[name.int,], row.names = name.int), data.frame(cboot$partition[name.int], row.names = name.int))

colnames(mass.neck) <- c("mass", "neck", "type")

#converte the cluster information to factors

mass.neck$type <- factor(mass.neck$type, levels = c("1","2"))
 
#------------------------------------------------------------------------------#

#================================================#
#whether neck length is correlated with body mass
#================================================#

#load the trees

load("tree.RData")

#separate the datasets

long <- mass.neck[mass.neck$type == "2",]
short <- mass.neck[mass.neck$type == "1",]

#prune the trees

tree.long <- ape::keep.tip.multiPhylo(tree, row.names(long))
tree.short <- ape::keep.tip.multiPhylo(tree, row.names(short))
tree.p <- ape::keep.tip.multiPhylo(tree, row.names(mass.neck))

#plot to check

plot(ladderize(tree.p[[1]]), no.margin = T, cex = 0.5)
plot(ladderize(tree.long[[1]]), no.margin = T)
plot(ladderize(tree.short[[1]]), no.margin = T)

#perform phylogenetic generalized least squares

#create lists to contain the results

fit <- list()
fit.long <- list()
fit.short <- list()

#here I employ the "moa" family as counters

moa <- 0
moa.long <- 0
moa.short <- 0

for(i in 1:100)
{
  fit[[i]] <- phylolm(formula = neck ~ mass, data = mass.neck[tree.p[[i]]$tip.label,], phy = tree.p[[i]], model = "lambda")
  fit.long[[i]] <- phylolm(formula = neck ~ mass, data = long[tree.long[[i]]$tip.label,], phy = tree.long[[i]], model = "lambda")
  fit.short[[i]] <- phylolm(formula = neck ~ mass, data = short[tree.short[[i]]$tip.label,], phy = tree.short[[i]], model = "lambda")

  #count the cases which fail to reject the null hopothesis
  
  if(summary(fit[[i]])$coefficients[2,4] > 0.05)
  {
    moa <- moa + 1
  }
  if(summary(fit.long[[i]])$coefficients[2,4] > 0.05)
  {
    moa.long <- moa.long + 1
  }
  if(summary(fit.short[[i]])$coefficients[2,4] > 0.05)
  {
    moa.short <- moa.short + 1
  }
}

#extract the residuals

short.residuals <- data.frame(col1 = numeric())
long.residuals <- data.frame(col1 = numeric())

for(i in 1:100)
{
    #here we summarize the residuals
    #in the same order as the first tree
    
    #I shall employ a temporary variable called "dodo"
    
    dodo <- fit.short[[i]]$residuals[tree.short[[1]]$tip.label]
    
    if(i == 1)
    {
      short.residuals <- dodo
    }
    else
    {
      short.residuals <- cbind(short.residuals, dodo)
    }
}

for(i in 1:100)
{
  #here we summarize the residuals
  #in the same order as the first tree
  
  #I shall employ a temporary variable called "dodo"
  
  dodo <- fit.long[[i]]$residuals[tree.long[[1]]$tip.label]
  
  if(i == 1)
  {
    long.residuals <- dodo
  }
  else
  {
    long.residuals <- cbind(long.residuals, dodo)
  }
}

#take the mean value of the residuals

short.p <- data.frame(matrix(nrow = 24, ncol = 3))
long.p <- data.frame(matrix(nrow = 24, ncol = 3))
colnames(short.p) <- c("mass", "neck", "type")
colnames(long.p) <- c("mass", "neck", "type")

short.p$mass <- short[tree.short[[1]]$tip.label,]$mass
long.p$mass <- long[tree.long[[1]]$tip.label,]$mass

#combine the datasets

for(i in 1:nrow(short.residuals))
{
  short.p[i,2] <- mean(short.residuals[i,1:100])
  
  short.p[i,3] <- "1"
}

for(i in 1:nrow(long.residuals))
{
  long.p[i,2] <- mean(long.residuals[i,1:100])
  
  long.p[i,3] <- "2"
}

row.names(short.p) <- row.names(short.residuals)
short.p$type <- factor(short.p$type)

row.names(long.p) <- row.names(long.residuals)
long.p$type <- factor(long.p$type)

#combine the long-necked and short-necked datasets

mass.neck.p <- rbind(long.p, short.p)

#------------------------------------------------------------------------------#

#===================================================#
#separate body mass into bins
#and investigate variance of neck length in each bin
#===================================================#

make_bins <- function(data, mass_col, trait_col, nbin, boot_reps = 100) {
  
  #step 1: sort the entire dataset in ascending order based on the body mass column
  data <- data[order(data[[mass_col]]), ]
  
  #safety check
  if (nrow(data) < nbin) {
    stop("Number of data rows is fewer than the number of bins (nbin).")
  }
  
  #step 2: create row indices for "fixed bins"
  fixed_indices <- split(seq_len(nrow(data)), cut(seq_len(nrow(data)), nbin, labels = FALSE))
  
  #step 3: create row indices for "overlapping bins"
  overlap_indices <- list()
  for (i in seq_len(nbin)) {
    current <- fixed_indices[[i]]
    if (nbin == 1) {
      overlap_indices[[i]] <- current
      next
    }
    if (i == 1) {
      next_idx <- fixed_indices[[i + 1]]
      overlap_indices[[i]] <- c(current, head(next_idx, floor(length(next_idx) / 2)))
    } else if (i == nbin) {
      prev_idx <- fixed_indices[[i - 1]]
      overlap_indices[[i]] <- c(tail(prev_idx, floor(length(prev_idx) / 2)), current)
    } else {
      prev_idx <- fixed_indices[[i - 1]]
      next_idx <- fixed_indices[[i + 1]]
      overlap_indices[[i]] <- c(tail(prev_idx, floor(length(prev_idx) / 4)), 
                                current, 
                                head(next_idx, ceiling(length(next_idx) / 4)))
    }
  }
  
  #define internal helper: dispaRity core workflow
  run_disparity_logic <- function(indices, trait_data, reps) {
    mat <- as.matrix(trait_data)
    rownames(mat) <- paste0("sp", seq_len(nrow(mat))) 
    
    names(indices) <- paste0("bin", seq_len(length(indices)))
    disp_obj <- custom.subsets(mat, group = indices)
    
    disp_boot <- boot.matrix(disp_obj, bootstraps = reps)
    disp_res <- dispRity(disp_boot, metric = c(sd))
    
    return(summary(disp_res))
  }
  # ------------------------------------------------
  
  trait_vector <- data[[trait_col]]
  
  fixed_boot_summary <- run_disparity_logic(fixed_indices, trait_vector, boot_reps)
  overlap_boot_summary <- run_disparity_logic(overlap_indices, trait_vector, boot_reps)
  
  #step 4: finalize data frames
  finalize_summary <- function(summ_df, idx_list, type_label) {
    
    mass_medians <- sapply(idx_list, function(idx) median(data[[mass_col]][idx], na.rm = TRUE))
    
    base_df <- data.frame(
      subset_name = paste0("bin", seq_along(idx_list)),
      bin = seq_along(idx_list),
      mass_median = mass_medians, 
      type = type_label,
      stringsAsFactors = FALSE
    )
    
    merged_df <- merge(base_df, summ_df, by.x = "subset_name", by.y = "subsets", all.x = TRUE)
    merged_df <- merged_df[order(merged_df$bin), ] 
    
    med_val <- if ("bs.median" %in% names(merged_df)) merged_df$bs.median else merged_df$obs
    lower_val <- if ("2.5%" %in% names(merged_df)) merged_df[["2.5%"]] else NA
    upper_val <- if ("97.5%" %in% names(merged_df)) merged_df[["97.5%"]] else NA
    
    final_df <- data.frame(
      bin = merged_df$bin,
      mass_median = merged_df$mass_median,
      trait_disparity = med_val,
      lower = lower_val,
      upper = upper_val,
      n_species = ifelse(is.na(merged_df$n), 0, merged_df$n),
      type = merged_df$type
    )
    
    return(final_df)
  }
  
  fixed_summary <- finalize_summary(fixed_boot_summary, fixed_indices, "Fixed")
  overlap_summary <- finalize_summary(overlap_boot_summary, overlap_indices, "Overlap")
  
  list(
    fixed_summary = fixed_summary,
    overlap_summary = overlap_summary
  )
}

#run the analysis

result <- make_bins(mass.neck.p, "mass", "neck", nbin = 5)

result$fixed_summary
result$overlap_summary

#------------------------------------------------------------------------------#

#===================================#
#plot (neck length across mass bins)
#===================================#

#summarize the data into a data frame

fixed_df <- result$fixed_summary %>% mutate(type = "fixed")
overlap_df <- result$overlap_summary %>% mutate(type = "overlap")

combined <- bind_rows(fixed_df, overlap_df) %>%
  filter(!is.na(trait_disparity)) %>%                
  arrange(mass_median) %>%                     
  mutate(seq = row_number())                  

#plot

bin.plot <- 
  ggplot(combined, aes(x = seq, y = trait_disparity)) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper), 
    width = 0.15,             
    linewidth = 1,          
    color = "#16a085",
    alpha = 0.7               
  ) +
  geom_line(linewidth = 1, color = "#16a085") +        
  geom_point(size = 3, color = "#16a085", alpha = 0.8) +       
  scale_x_continuous(
    breaks = combined$seq,                           
    labels = paste0(combined$type, combined$bin, "\n", 
                    round(10^(combined$mass_median), 2), "kg")  
  ) +
  labs(
    x = "bin order (sorted by median mass)",
    y = "standard deviation",
    title = "standard deviation of neck length"
  ) +
  theme_bw(base_size = 7)+
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(5, 5, 5, 5)
  )

bin.plot

#------------------------------------------------------------------------------#

#=========================#
#plot (cluster dendrogram)
#=========================#

#extract the clustering result

clus.res <- cboot$result$result

#get the grouping information

groups <- cboot$partition

#treat the clustering result as a phylogeny

clus.tree <- as.phylo(clus.res)

#change some names within the clustering tree and the grouping information

clus.tree$tip.label[1] <- '"Monquirasaurus"_boyacensis'
clus.tree$tip.label[2] <- "Serpentisuchops_pfisterae"
clus.tree$tip.label[3] <- "Plesionectes_longicollum"
clus.tree$tip.label <- gsub("_", " ", clus.tree$tip.label)

names(groups)[[1]] <- '"Monquirasaurus"_boyacensis'
names(groups)[[2]] <- "Serpentisuchops_pfisterae"
names(groups)[[3]] <- "Plesionectes_longicollum"
names(groups) <- gsub("_", " ", names(groups))

#mapping the clustering information onto the branches

cluster_list <- split(names(groups), groups)

clus_grouped <- groupOTU(clus.tree, cluster_list)

#calculate the tree height

max_step <- max(node.depth(clus.tree))

#find the "ancestor" node of the two clusters

tips_short <- names(groups[groups == "1"])
tips_long <- names(groups[groups == "2"])

node_short <- getMRCA(clus.tree, tips_short)
node_long <- getMRCA(clus.tree, tips_long)

#extract the bootstrap value

bs_short <- round(cboot$bootmean[1], 3)
bs_long <- round(cboot$bootmean[2], 3)

#plot

cluster.plot <- 
  ggtree(clus_grouped, 
         ladderize = T, 
         linewidth = 0.3, 
         branch.length = "none",
         aes(color = group)) +
  geom_tippoint(
    aes(fill = group),
    size = 1,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  geom_tiplab(
    fontface = "italic", 
    size = 2.2,
    offset = 0.2
  ) +
  geom_text2(
    aes(subset = (node == node_short)),
    label = "short",  
    size = 3,
    fontface = "bold",
    hjust = 1.1,
    vjust = -0.5
  ) +
  geom_text2(
    aes(subset = (node == node_short)),
    label = bs_short,  
    size = 3,        
    fontface = "plain",
    hjust = 1.1,
    vjust = 1.5   
  ) +
  geom_text2(
    aes(subset = (node == node_long)),
    label = "long",
    size = 3,
    fontface = "bold",
    hjust = 1.1,
    vjust = -0.5
  ) +
  geom_text2(
    aes(subset = (node == node_long)),
    label = bs_long,  
    size = 3,
    fontface = "plain",
    hjust = 1.1,
    vjust = 1.5
  ) +
  scale_color_manual(
    values = c("1" = "#5061c5CC", "2" = "#F34C3FCC")
  ) +
  scale_fill_manual(
    values = c("1" = "#5061c5CC", "2" = "#F34C3FCC"),
    guide = "none"
  ) +
  xlim(0, 12)+
  theme(
    legend.position = "none", 
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(5, 20, 5, 0)
  )

cluster.plot

#------------------------------------------------------------------------------#

#=============================#
#plot (mass-neck scatter plot)
#=============================#

#summarize the ranges of the PGLS models

long.params <- 
  do.call(rbind, lapply(fit.long, function(mod) {
  coefs <- mod$coefficients
  data.frame(
    intercept = coefs[1],
    slope     = coefs[2]
  )
}))

short.params <- 
  do.call(rbind, lapply(fit.short, function(mod) {
    coefs <- mod$coefficients
    data.frame(
      intercept = coefs[1],
      slope     = coefs[2]
    )
  }))

#define the range of body mass

x_range <- seq(min(mass.neck$mass), max(mass.neck$mass), length.out = 200)

#calculate the corresponding neck values

lines_long <- sapply(1:nrow(long.params), function(i) {
  long.params$intercept[i] + long.params$slope[i] * x_range
})

lines_short <- sapply(1:nrow(short.params), function(i) {
  short.params$intercept[i] + short.params$slope[i] * x_range
})

#take the median and 95% intervals
long_median <- apply(lines_long, 1, median)
long_lower  <- apply(lines_long, 1, quantile, probs = 0.025)
long_upper  <- apply(lines_long, 1, quantile, probs = 0.975)

short_median <- apply(lines_short, 1, median)
short_lower  <- apply(lines_short, 1, quantile, probs = 0.025)
short_upper  <- apply(lines_short, 1, quantile, probs = 0.975)

#create the data frames for plotting

ribbon_long <- data.frame(x = x_range, ymin = long_lower, ymax = long_upper)
line_long <- data.frame(x = x_range, y = long_median)

ribbon_short <- data.frame(x = x_range, ymin = short_lower, ymax = short_upper)
line_short <- data.frame(x = x_range, y = short_median)

#fit the OLS models

long.ols <- lm(neck ~ mass, data = mass.neck[mass.neck$type == "2",])
short.ols <- lm(neck ~ mass, data = mass.neck[mass.neck$type == "1",])

#create data frames for plotting

line_long_ols <- data.frame(
  x = x_range,
  y = long.ols$coefficients[1] + long.ols$coefficients[2] * x_range
)

line_short_ols <- data.frame(
  x = x_range,
  y = short.ols$coefficients[1] + short.ols$coefficients[2] * x_range
)

scatter.plot <-
  ggplot(data = mass.neck, aes(x = mass, y = neck,
                               color = factor(type))) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(
    values = c("1" = "#5061c5CC", "2" = "#F34C3FCC"),          
    labels = c("1" = "short", "2" = "long"),    
    name = "neck type"                                   
  )+
  labs(title = "mass-neck scatter plot",
       x = expression("log"["10"]*"(body mass)"),
       y = "neck PCo1 (98.8% of variance)",
       color = "neck type") +
  geom_ribbon(data = ribbon_long, 
              aes(x = x, ymin = ymin, ymax = ymax),
              fill = "#F34C3FCC", alpha = 0.35, inherit.aes = F) +
  geom_line(data = line_long, 
            aes(x = x, y = y), 
            color = "#F34C3FCC", linewidth = 1, inherit.aes = F) +
  geom_line(data = line_long_ols,
            aes(x = x, y = y),
            color = "#F34C3FCC", linewidth = 1, inherit.aes = F,
            linetype = "dashed")+
  geom_ribbon(data = ribbon_short, 
              aes(x = x, ymin = ymin, ymax = ymax),
              fill = "#5061c5CC", alpha = 0.35, inherit.aes = F) +
  geom_line(data = line_short, 
            aes(x = x, y = y), 
            color = "#5061c5CC", linewidth = 1, inherit.aes = F) +
  geom_line(data = line_short_ols,
            aes(x = x, y = y),
            color = "#5061c5CC", linewidth = 1, inherit.aes = F,
            linetype = "dashed")+
  theme_bw(base_size = 7) +
  theme(
    legend.position = "top",
    legend.key = element_rect(fill = "transparent"),
    legend.direction = "horizontal",
    legend.text = element_text(size = 11),     
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(5, 5, 5, 5)
  )


scatter.plot

#save the plots

#setwd()

saveRDS(scatter.plot, file = "mass_neck_scatter.rds")
saveRDS(bin.plot, file = "neck_binplot.rds")
