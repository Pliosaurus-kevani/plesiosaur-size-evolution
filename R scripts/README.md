# R scripts

Commented R scripts to enable readers to reproduce the analyses. All scripts were run based on [R](https://www.r-project.org/) 4.5.2.

## figure 1

**partitioned_mapping.R**:
>Script for subfigures A-D of figure 1 presented in the manuscript. It maps the branch-specific rates of morphological evolution 
>onto the branches of the maximum a posteriori trees of the unpartitioned and partitioned analysis, respectively.

**macro_parameters.R**:
>Script for subfigures E-H of figure 1 presented in the manuscript. It creates violin plots of the speciation, extinction, net diversification, and fossil sampling rates across the five time bins.

## figure 2

**mass_mapping.R**:
>Script for subfigure A of figure 2 presented in the manuscript. It performs ancestral state reconstruction on the pruned maximum a posteriori tree, and maps body mass values onto the branches.

**multi_mode_mass_evolution.R**:
>Script for subfigure C of figure 2 presented in the manuscript. It fits multiple models to the body mass evolution of plesiosaurs. A figure file named **aicc.plot.rds** is generated in this R script. This script may take more than twenty minutes to run.

**BM_OU.R**:
>Script for subfigures B and D-F of figure 2 presented in the manuscript. It generates the body mass distribution plot of the three subsets, and summarizes the results of model comparison (Brownian motion vs Ornstein-Uhlenbeck process) in each subset. Before running this script, model fitting results performed in [RevBayes](https://revbayes.github.io/) are required (see the [BM_OU](https://github.com/Pliosaurus-kevani/plesiosaur-size-evolution/tree/main/BM_OU) folder in this repository). Alternatively, I have also provided my RevBayes results in this repository (see the [BM_OU_output](https://github.com/Pliosaurus-kevani/plesiosaur-size-evolution/tree/main/BM_OU_output) folder). Figure files including  **mass.distribution.rds**, **alpha.plot.rds**, **theta.plot.rds**, **bar.plot.rds** are generated in this R script.

**figure2_patch.R**:
>Script for figure 2 presented in the manuscript. Before running this file, you need to run the R files **mass_mapping.R**, **multi_mode_mass_evolution.R**, and **BM_OU.R**, as this script aims to combine the subfigures generated therein to create figure 2.

## figure 4

**rate_distribution_par.R**:
>Script for subfigures B-D and F-H of figure 4 presented in the manuscript. It generates the hexbin plots and mass distribution of branches that exhibit rate bursts in craniodental, axial, or appendicular characters. It generates a figure file named as **morpho_rate_partitioned_all branch.rds**.

**rate_distribution_unpar.R**:
>Script for subfigures A and E of figure 4, and subfigure E of figure 3. It generates the hexbin plot and mass distribution of branches that exhibit rate bursts in the unpartitioned dataset, as well as the raincloud plots of branch mass of four subsets. This script also combines the subfigures **morpho_rate_partitioned_all branch.rds** to create figure 4, so **rate_distribution_par.R** should be run before this script.

## figure 3

**neck_length_mass.R**:
>Script for subfigures A and B of figure 3 presented in the manuscript, as well as the cluster dendrogram presented in figure S2. This R script performs phylogenetic generalized least squares between body mass and neck principal coordinate, and generates figure files **mass_neck_scatter.rds** and **neck_binplot.rds**.

**neck_rate_mass.R**:
>Script for subfigures C and D of figure 3 presented in the manuscript. It creates files for the variable-rates model, which was run on [BayesTraits](https://www.evolution.reading.ac.uk/BayesTraitsV5.0.0/BayesTraitsV5.0.0.html) 4.0.0. To run BayesTraits and Posterior Processor, you need to download the software from the [official website](https://www.evolution.reading.ac.uk/SoftwareMain.html), and include the two tools in the folder where the scripts and data files are contained. I have uploaded my results in the folder [BayesTraits_output](https://github.com/Pliosaurus-kevani/plesiosaur-size-evolution/tree/main/BayesTraits_output). This R script also combines the subfigures to create figure 3. The .rds files generated in R scripts **rate_distribution_unpar.R** and **neck_length_mass.R** are required to successfully run this script. 

## figure S1

**plot_consensus.R**:
>Script for figure S1. It plots the 50% majority-rule consensus tree.

## other analyses

**rate_mass_PGLS.R**:
>Script for phylogenetic generalized least squares between body mass and morphological rates in unpartitioned, craniodental, axial, or appendicular characters.
