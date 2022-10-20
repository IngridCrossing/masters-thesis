# masters-thesis
Code for the simulation and inference model from Crossing et al. 2022

The code for running the simulation and inference model for each design space is provided here. The code was designed to run on Spartan, the University of Melbourne's HPC System for parallel processing, therefore it will require alterations to run on R. 

The following packages must be loaded into R:

readxl
tidyr
rjags
dplyr
jagsUI
magrittr
ggplot2


The main script for each Design Space i is found in Design_Spaces/Design_Space_i/scripts. Each directory also contains a job submission and batch submission slurm script for running on Spartan. By manually assigning values to design_index, detProb_index, and path_index, single iterations of the model can be run. 

In each directory Design_Spaces/Design_Space_i/results, the required folders for saving results are provided. 