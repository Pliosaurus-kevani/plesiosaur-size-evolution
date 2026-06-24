# Data for R

Data files required to run the [R scripts](https://github.com/Pliosaurus-kevani/plesiosaur-size-evolution/tree/main/R%20scripts).

**consensus_combined.tre**:
>50% majority-rule consensus tree of the unpartitioned analysis.

**map_combined.tre**:
>Maximum a posteriori tree of the unpartitioned analysis.

**map.craniodental.tre**, **map.axial.tre**, and **map.appendicular.tre**:
>Maximum a posteriori tree of the partitioned analysis. The three files should be completely the same tree, and they were created to ensure the results of the partitioned analysis were exported correctly.

**mass.csv**:
>Body mass file of 89 plesiosaur species.

**mass_all.csv**:
>Body mass file of all 130 plesiosaur taxa presented in the trees. Taxa without body mass estimates are marked as NA.

**neck data.csv**:
>Neck length file of plesiosaurs

**results_run_1.log**:
>Parameters of the first run of the unpartitioned phylogenetic analysis.

**tree.RData**:
>A hundred randomly selected posterior trees of the unpartitioned analysis. The variable **TREE** is a list containing *treedata* files (see, for example, R package [tidytree](https://cran.r-project.org/web/packages/tidytree/index.html) for introduction). The variable **tree** is a list containing *phylo* structures.

**tree_partitioned.RData**
>A hundred randomly selected posterior trees of the partitioned analysis.
