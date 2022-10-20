

command_args <- commandArgs(trailingOnly = TRUE) #Takes arguments for the .slurm submission files for parallel processing in Spartan


design_index <- as.numeric(command_args[1]) #Index for survey design for this iteration
detProb_index <- as.numeric(command_args[2]) #Index for detection probability for this iteration
path_index <- as.numeric(command_args[1]) #Index for path to save outputs for this iteration



path_options<-rep(NA, 10)
for(i in 1:length(path_options)){
  path_options[i] <- paste0("Design_Spaces/Design_Space_3/results/Design_", as.character(i), "/")
}
source("General_Scripts/functions.R")

designs_list <- vector(mode = "list", length = 29)
for(i in 1:29){
  d = i+1 #
    trips <- ifelse(d %% 3 > 0,  ((d - (d%%3))/3 + 1), (d - (d%%3))/3) #Finds the number of max  3-day field trips for the given survey number
    trip.start <- seq(from = 1, to = 61, length.out = trips) #Find the start day for eacg 3-day trip
    trips.list <- vector(mode = "list", length = length(trip.start)) #creates empty list
  for(j in 1:length(trip.start)){
    trips.list[[j]] <- seq(from = trip.start[j], by = 1, length.out = 3) #Sequences 3 days from each start day and stores them in list
  }

survey.days <- c(trips.list %>% unlist()) #joins all the survey days back together from the list

if(d%%3 > 0 ){
  survey.days %<>% head(-(3-(d%%3))) #Subtracts surveys if there are remainders
}

designs_list[[i]] <- survey.days %>% round() #Rounds any decimal numbers
}

keep <- seq(2, 29, 3) #Index for choosing the correct designs to keep

designs_list <- designs_list[keep] 

#Corrections to the designs list
designs_list[[9]][5] <- 9
designs_list[[9]][11]<- 25
designs_list[[9]][[17]] <- 39
designs_list[[9]][[23]] <- 55


det_probs <- seq(0.1, 1, 0.1) #Detection probabilities

survey_design <- designs_list[[design_index]] #Survey design for this iteration
det_prob <- det_probs[detProb_index] #Detection probabilities for this iteration
path <- path_options[path_index] #Path options for this iteration
csv.name <- paste0("mark_recapture", as.character(det_prob), ".csv")

print("Iteration:")
print(design_index)
print(det_prob)



#Infection values 
#uninf<-year.to.day(0.694) #Remains constant at rate taken from West et al. 2019

n = 100 #Initial size
uninf = 0.05 #Daily probability of infection
inf = 0.85 #Daily probability of remaining infected
fecundity = 720 #Average number of offspring per breeding-age female
adult.age = 3 #age at which individuals can breed and be captured in CMR
off.mortality = 0.989 #Probability of tadpoles <1 year old dying following a breeding event
ovi.sites = 50 #Number of oviposition sites available (density dependence)
age.classes = 4 #Number of age classes
alpha = c(6.374251, 6.428781, 6.706079, 6.706079)
beta = c(-0.5162382, -0.5408142, -0.6679783, -0.6679783) # Original beta
#beta = c(-0.5162382, -0.5408142, -0.993452, -0.993452) # new betas 
max.age = 12 #Maximum age individuals can reach
init.years=20 #Number of years the simulation runs for before survey begins
survey.years=1 #Number of years to be surveyed 




loop.recaptures.inf(n = n, uninf = uninf, fecundity = fecundity, adult.age = adult.age, off.mortality = off.mortality, ovi.sites = ovi.sites, age.classes = age.classes, alpha = alpha, max.age = max.age, det.prob = det_prob, init.years = init.years, survey.years = survey.years, path = path, survey.days = survey_design, inf = inf) #Run simulation and inference
