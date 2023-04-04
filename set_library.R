.libPaths("../R/lib")

all_packages <- c("readxl", "tidyr", "rjags", "dplyr", "jagsUI", "magrittr", "ggplot2")
lapply(all_packages, require, character.only = TRUE) 