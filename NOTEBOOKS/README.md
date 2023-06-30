# Notebooks and pythonic libraries

## Python scripts

+ simulation.py - python code for generating synthetic replicates from input bam file 

+ submitsim.py - python code for formating and submitting calls to simulation.py to slurm 

+ bootlegger.py - Library of python code for conducting bootstrapping and statistical testing

+ epigenomievisulaiztion.py - Library of python code for plotting epigenoimc datasets

+ mystatslib.py - Library of python code for conducting correlation analysis between ATAC-seq experiment

+ mystatslib_wcz.py - Same as above with co-zeros retained between comparisons. 

## Notebooks

Notebooks should be run from top to bottom as presented here. 

+ Mapping_stats_and_counts_Table1_STable1.ipynb - Formats tables for publication of mapping statistics from output logs. 

+ Plot_Simulation_Results_Figure3_SFigure1.ipynb - Plots results from statistical simulations on ATAC-seq data. Generates Figure 3 and Supplementary Figure 1. 

+ Real_Correlation_Counts_Figure1_Figure2_Figure4_SFigure2_SFigure3.ipynb - Calculates the real correlation and association metrics between ATAC-seq data from A549 experiments and ENCODE samples. Cozeros are removed from analysis. Generates Figures 1, 2, and 4 as well as supplementary Figures 2 and 3. 

+ Real_Correlation_Counts_with_cozeros.ipynb - Same as above, but co-zeros are retained. 

+ Random_Forest_Correlation_Replication_Figure5_SFigure7.ipynb - Conducts classification on ATAC-seq replicate class using a random forest. Generates Figure 5 and supplementary Figure 7. 

+ Diff_Real_Corr_Co-zero_Removed_SFigure6.ipynb - Calculates the difference in correlation and assoicaiton statistics run on real ATAC-seq experiments. Generates Supplementary Figure 6 in Roth et al 2023. 
