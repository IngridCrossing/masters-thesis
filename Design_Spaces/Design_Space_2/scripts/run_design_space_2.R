

command_args <- commandArgs(trailingOnly = TRUE) #Takes arguments for the .slurm submission files for parallel processing in Spartan



design_index <- as.numeric(command_args[1]) #Index for survey design for this iteration
detProb_index <- as.numeric(command_args[2]) #Index for detection probability for this iteration
path_index <- as.numeric(command_args[1]) #Index for path to save outputs for this iteration


design_index <- 1
detProb_index <- 1
path_index <- 1



path_options<-rep(NA, 10)
for(i in 1:length(path_options)){
  path_options[i] <- paste0("Design_Spaces/Design_Space_2/results/Design_", as.character(i), "/")
}

source("General_Scripts/functions.R")



designs_list <- vector(mode = "list", length = 20) #Creates survey designs with increasing intervals between surveys
for (i in 1:length(designs_list)){
  designs_list[[i]] <- seq(1, 61, by = i)
}


det_probs <- seq(0.1, 1, 0.1) #Detection probabilities

survey_design <- designs_list[[design_index]] #Survey design for this iteration
det_prob <- det_probs[detProb_index] #Detection probability for this iteration
path <- path_options[path_index] #Path for this iteration
csv.name <- paste0("mark_recapture", as.character(det_prob), ".csv") #mark recapture data csv name for this iteration

print("Iteration:") 
print(design_index)
print(det_prob)


n = 100 #Initial size
uninf = 0.05 #Daily probability of becoming infected
inf = 0.85 #Daily probability of remaining infected
fecundity = 720 #Average number of offspring per breeding-age female
adult.age = 3 #age at which individuals can breed and be captured in CMR
off.mortality = 0.989 #Probability of tadpoles <1 year old dying following a breeding event
ovi.sites = 50 #Number of oviposition sites available (density dependence)
age.classes = 4 #Number of age classes
alpha = c(6.374251, 6.428781, 6.706079, 6.706079) #Original alphas
beta = c(-0.5162382, -0.5408142, -0.6679783, -0.6679783) # Original beta
max.age = 12 #Maximum age individuals can reach
init.years=20 #Number of years the simulation runs for before survey begins
survey.years=1 #Number of years to be surveyed 


loop.recaptures.inf(n = n, uninf = uninf, fecundity = fecundity, adult.age = adult.age, off.mortality = off.mortality, ovi.sites = ovi.sites, age.classes = age.classes, alpha = alpha, max.age = max.age, det.prob = det_prob, init.years = init.years, survey.years = survey.years, path = path, survey.days = survey_design, inf = inf) #Run simulation and inference

