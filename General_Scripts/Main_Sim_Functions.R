

####################---IBM Functions---#####################

#-----Initialise a population--------#
#Uses random draws from a binomial distribution to determine sex
initialise<-function(n=100, age.classes=4, adult.age=3){
  inds<-array(data=0, dim = c(n, 7))   #Creates array with n rows and 7 columns
  colnames(inds)<-c("sex", "infected", "repro", "death", "age", "ID", "marked") #Adds column names
  inds[, "sex"]<-rbinom(n=n, size= 1, prob = 0.5) + 1  # Female = 2, Male = 1
  inds[, "infected"]<-1  # Infected = 1, uninfected = 2
  inds[, "age"]<- sample(adult.age:age.classes, size=n, replace=TRUE) #Population of all adults
  inds[, "ID"]<-c(1:n)  # Assigns ID 
  ID.bank <<-as.integer(inds[, "ID"])
  return(inds)
}


#------Infection------#
infection<-function(inds, uninf, inf){ # maps to B.I.c <- uninf; b.S.c <- 1-inf
  total_inds<-dim(inds)[1]  #Total number of individuals (rows in column 1 of inds)
  pInf<-c(uninf, inf) #Probability of being infected in the next timestep if uninfected, infected
  pI<-pInf[inds[, "infected"]] #Assigns pInf according to whether infection status is 1 (uninfected) or 2 (infected)
  inds[, "infected"]<-rbinom(total_inds, 1, pI) +1
  return(inds)
}

#-----Reproduction------#
#each adult produces a mean of 720 offspring
#Survival probability of eggs is 0.011
reproduction<-function(inds, n=100, fecundity = 720, adult.age = 4, off.mortality = 0.980, ovi.sites = 20){
  young_females<-which(inds[, "sex"]==2 & inds[, "age"]==adult.age)
  ovi.sites <- rpois(1, ovi.sites)
  female_sample_size<-sum(rbinom(length(young_females), 1, 0.5))
  breeding_females<-sample(young_females, female_sample_size, replace = FALSE)
  female_row<-which(inds[, "sex"]==2 & inds[, "age"]> adult.age) #Identifies adult individuals not in first year stores row numbers
  female_row <- append(female_row, breeding_females)
  if(length(female_row)>=ovi.sites){
    female_row<-sample(female_row, ovi.sites, replace = FALSE) #Takes ovi.sites number random draws from female_row without replacement
  }
  inds[female_row, "repro"]<-rpois(length(female_row), fecundity) #Every female produces roughly 720 offspring (drawn from a poisson distribution with mean fecundity)
  total_off<-sum(inds[, "repro"]) #Calculates total number of offspring produced
  total_off<-sum(rbinom(total_off, 1, (1-off.mortality))) #Total number of offspring that will survive
  
  inds_cols<-dim(inds)[2] #Number of columns in inds to be replicated in data aeeay
  new_inds<-array(data = 0, dim = c(total_off, inds_cols)) #Creates an empty data frame with total number of surviving offspring as rows, and columns same as original data frame
  new_inds[, 1]<-rbinom(total_off, size = 1, prob = 0.5) + 1 #Assigns sex with a probability of 0.5 
  new_inds[, 2]<-1 #Assigns infection status as 1 (uninfected)
  
  #Assigning ID
  first.ID <- max(ID.bank) + 1 #Finds the last ID given to an individual, then adds 1
  new.IDs <- seq(from = first.ID, by = 1, length.out = dim(new_inds)[1])
  ID.bank <<-append(ID.bank, new.IDs)
  new_inds[, 6]<-new.IDs
  
  new_inds[, 7]<-0 
  new_inds[, 5]<-1 #Puts them all in the first age class
  
  #Attach new offspring to inds array
  inds<-rbind(inds, new_inds)
  return(inds)
}

#----------Death-------------#
#Includes create mortality matrix function
#First, must create matrix for mortality for each age class, and added mortality due to infection

create.mortality.matrix<-function(death.params, age.classes){
  if(((length(death.params))/age.classes)==2){
    mortality<<-matrix(death.params, nrow = age.classes, ncol = 2, byrow = TRUE)
    age.classes<<-age.classes
    print(mortality)
  }
  else{ 
    stop("Number of death parameters not equal to age classes * 2")}
}   

logit.function <- function(a, b, inf = FALSE){ #Takes the alpha and beta values, and infection status, and calculates the log-odds of survival 
  if(inf == TRUE){
    x <- 1
  } else {
    x <- 0
  }
  log.odds <- a + b*x 
  return(log.odds)
}


logit.to.prob <- function(logodds){
  p <- exp(logodds)/(1+exp(logodds))
  return(p)
}

#then, function for death 

death<-function(inds, age.classes = 4, mortality.matrix = mortality, max.age = 13, alpha = c(6.374251, 6.428781, 6.706079, 6.706079), beta = c(-0.5162382, -0.5408142, -0.6679783, -0.6679783)){
  total_inds<-dim(inds)[1] #Total number of individuals in current population
  #Assign probability of death due to age 
  age<-inds[, "age"]  #Creates a vector of the ages of individuals
  prob_death_age <- 1-(logit.to.prob(alpha[-age.classes]))#Probability of death for individual ages (eg. 1, 2, 3 but not 4+ year olds)
  prob_death<-prob_death_age[age] #Assigns probability of death based on age for 1, 2 and 3 year olds. Adults currently NA 
  prob_death[is.na(prob_death)] = 1-(logit.to.prob(alpha[age.classes])) #Replaces all NAs with mortality for individuals in "final" age and older
  
  inf.alpha.beta <- logit.function(a = alpha, b = beta, inf = TRUE) #Find log odds of survival for infected individuals
  
  
  #Alters probability based on infection status:
  prob_death_inf<- 1-(logit.to.prob(inf.alpha.beta[-age.classes])) #Additional additional probability of death with infection for all ages <4
  inf_prob_death<-prob_death_inf[age] #Assign probability of death based on age
  inf_prob_death[is.na(inf_prob_death)] <- 1-(logit.to.prob(inf.alpha.beta[age.classes])) #Replaces all NAs with mortality for individuals in "final" age or older
  inf_rows<- which(inds[, "infected"]==2)  #Finds all the individuals who are infected
  
  prob_death[inf_rows] <- inf_prob_death[inf_rows]
  
  
  #Now have a vector of probability of dying for every age class + infection status
  inds[, "death"]<-rbinom(total_inds, 1, prob_death)
  elderly <- which(inds[, "age"] >= max.age) #Finds individuals who are 13 or older
  inds[elderly, "death"] <- 1
  inds<-inds[inds[, "death"]==0, , drop = FALSE]
  return(inds)
}


#-----------Mark-Recapture Function--------------#

mark_recapture<-function(
  n = 100, age.classes = 4, adult.age = 3,
  fecundity = 720, off.mortality = 0.98, ovi.sites = 50,
  uninf, inf, 
  alpha = c(6.374251, 6.428781, 6.706079, 6.706079), 
  beta = c(-0.5162382, -0.5408142, -0.6679783, -0.6679783), 
  max.age = 13,
  survey.days=c(30, 60, 90), det.prob = 0.7,
  init.years = 10,
  survey.years=3,
  save.csv = TRUE, csv.name = "cmr_data.csv"){
  years = init.years + survey.years
  pop.history=matrix(0, nrow = 1, ncol = years*365)  #Empty matrix to store population information
  inds<-initialise(n = n, age.classes = age.classes, adult.age = adult.age) #initialise the population
  adult.history = matrix(0, nrow = 1, ncol = years*365) #Empty matrix for storing information about adults
  
  
  for(i in 1:years){
    
    if(i<=init.years){
      inds[, "age"]<-(inds[, "age"]+1) #Age up individuals by 1 year  
      inds<-reproduction(inds, n=n, fecundity = fecundity, adult.age = adult.age, off.mortality = off.mortality, ovi.sites = ovi.sites) #Reproduce at the start of the year 
      for (j in 1:365){
        inds %<>% 
          infection(uninf = uninf, inf = inf) %>% #infection dynamics
          death(age.classes = age.classes, alpha = alpha, beta = beta, max.age = max.age) #death
        total_pop<-dim(inds)[1] #sums total population
        pop.history[, i*365-365+j]<-total_pop #Saves total population number to population matrix
        
        
      }
    }
    
    else{
      inds[, "age"]<-(inds[, "age"]+1) #Age up individuals by 1 year
      inds<-reproduction(inds, n=n, fecundity = fecundity, adult.age = adult.age, off.mortality = off.mortality, ovi.sites = ovi.sites) #Reproduce at the start of the year 
      
      for (j in 1:365){
        inds %<>% 
          infection(uninf = uninf, inf = inf) %>% #Infection dynamics
          death(age.classes = age.classes, alpha = alpha, beta = beta, max.age = max.age) #Death
        total_pop<-dim(inds)[1] #record total population number
        pop.history[, i*365-365+j]<-total_pop #Save total population number
        
        if(j %in% survey.days){
          adult_inds<-which(inds[, "age"] >= adult.age) #Finds rows where individuals are adults 
          
          total_adults<-length(adult_inds) #Total number of adults
          total_marked<-sum(rbinom(total_adults, 1, det.prob)) #Chooses teh number of adult invidivuals that will be sampled based on detection probability 
          capture_rows<-sample(adult_inds, total_marked, replace = FALSE) #Samples adult population
          capture_rows %<>% sort()
          capture_sample<-inds[capture_rows, , drop = FALSE] #Subsets individuals that are captured
          
          matrix_position<-match(j, survey.days) 
          #Individuals that are being re-captured (have been seen before)
          #Uninfected 
          
          if(i == init.years+1 & matrix_position==1){ #If it's the first year of capture (init.years+1)
            
            cmr_data<-matrix(0, nrow = length(capture_rows), ncol = survey.years*length(survey.days))  #Create empty matrix for cmr data
            rownames(cmr_data)<-as.character(capture_sample[, "ID"])
            colnames(cmr_data)<-as.character(rep(c(1:survey.years), each = length(survey.days))) #Label columns of cmr data by year
            uninf_rows<-which(capture_sample[, "infected"]==1) #Identify uninfected individuals
            uninf_id<-as.character(capture_sample[uninf_rows, "ID"]) #Identify ID of uninfected individuals
            
            inf_rows<-which(capture_sample[, "infected"]==2) #Identify infected individuals
            inf_id<-as.character(capture_sample[inf_rows, "ID"]) #Identify ID of infected individuals
            
            cmr_data[uninf_rows, 1]<-1 
            cmr_data[inf_rows, 1]<-2 
            
            inds[capture_rows, "marked"] <-1
            
          }
          
          else{
            
            recapture_marked_uninf_rows<-which(capture_sample[, "marked"]==1 & capture_sample[, "infected"]==1)
            recapture_marked_uninf_id<-as.character(capture_sample[recapture_marked_uninf_rows, "ID"])
            
            
            #Infected
            recapture_marked_inf_rows<-which(capture_sample[, "marked"]==1 & capture_sample[, "infected"]==2)
            recapture_marked_inf_id<-as.character(capture_sample[recapture_marked_inf_rows, "ID"])
            
            
            
            cmr_data[recapture_marked_uninf_id, ((length(survey.days)*(i-init.years-1))+matrix_position)]<-1 #Assign 1 to cmr matrix for uninfected individuals
            cmr_data[recapture_marked_inf_id, ((length(survey.days)*(i-init.years-1))+matrix_position)]<-2 #Assign 2 to cmr matrix for infected individuals
            
            #individuals being capture for the first time 
            
            #Uninfected 
            recapture_unmarked_uninf_rows<-which(capture_sample[, "marked"]==0 & capture_sample[, "infected"]==1)
            recapture_unmarked_uninf_id<-as.character(capture_sample[recapture_unmarked_uninf_rows, "ID"])
            new_matrix_rows<-matrix(0, nrow = length(recapture_unmarked_uninf_id), ncol = dim(cmr_data)[2])
            rownames(new_matrix_rows)<-recapture_unmarked_uninf_id
            new_matrix_rows[, ((length(survey.days)*(i-init.years-1))+matrix_position)]<-1
            cmr_data %<>% rbind(new_matrix_rows)
            
            #Infected 
            recapture_unmarked_inf_rows<-which(capture_sample[, "marked"]==0 & capture_sample[, "infected"]==2)
            recapture_unmarked_inf_id<-as.character(capture_sample[recapture_unmarked_inf_rows, "ID"])
            new_matrix_rows<-matrix(0, nrow = length(recapture_unmarked_inf_id), ncol = dim(cmr_data)[2])
            rownames(new_matrix_rows)<-recapture_unmarked_inf_id
            new_matrix_rows[, ((length(survey.days)*(i-init.years-1))+matrix_position)]<-2 #Was 1, now correct to 2
            cmr_data %<>% rbind(new_matrix_rows)
            
            #Finally, mark new individuals 
            inds[capture_rows, "marked"]<-1
            
          }
          
        }
        
      }
    }
  }
  inds<<-inds #Save final inds to global environment
  cmr_data<<-cmr_data
  pop.history<<-pop.history #Save population history to global environment
  matplot(1:(years*365), t(pop.history), type = "l", lty = 1, ylab = "Number of Individuals", xlab = "Days")
  
  
  if(save.csv == TRUE){
    write.csv(cmr_data, csv.name)}
}


###-----------Functions for Inference Model-------------#

#From Matt's Model, line 38
########################################################
###     define functions to prep data for model   ######
########################################################

# Function to collapse capture history matrix into unique histories and the
# frequency of these unique histories. Returns a list: first element is the
# collapsed capture histories, second element is the frequency.


collapse.ch <- function(ch){
  ch.char = apply(ch, 1, function(x) paste(x, collapse = ","))
  
  ch.sum.out = t(sapply(strsplit(names(table(ch.char)), split = ","), as.numeric))
  fr.out = as.numeric(as.vector(table(ch.char)))
  
  return(list(ch.sum.out, fr.out))
}

# Function to create the matrix of latent state z 
known.state.cjs <- function(ch){
  state <- ch
  for(i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,] != 5))
    n2 <- max(which(ch[i,] != 5))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state == 0] <- NA
  return(state)
}

# Create vector with occasion of marking
get.first <- function(x) min(which(x != 5)) #find the first value that is not 5 i.e. the first detection of a frog



# Function to create vector of suml

find.suml <- function(x){which(x !=0) %>%
max()}


##--------------- csv naming function-------------------#


### Function for creating names to save CMR csv's (Works from 1-99)
create.csv.names.det = function(min.detect= 0, max.detect = 1, detect.interval = 0.1){             #Recaptures specifies the number of recaptures per year
  suffix.list<-seq(min.detect, max.detect, detect.interval)    #Sequences from min det prob to max det prob by a specified interval
  suffix.list<-format(suffix.list, nsmall = 1)                 #Changes 0 and 1 to 0.0 and 1.0 for consistency (nsmall = the minimum number of digits to the right of the decimal place)
  total.csvs<-length(suffix.list)                             #Saves the total number of possibly det probabilitys ie. the total number of required csvs
  csv.names<-(rep("mark.recapture.xxx.det.prob.csv", total.csvs))   #Creates dataframe of strings of "mark.recapture.xx.csv"
  for (i in 1:total.csvs){
    substr(csv.names[i], 16, 18)<-as.character(suffix.list[i])} #Replaces xxx (characters 16 to 19) in each csv name with the detection probability 
  csv.names<<-csv.names                               #Saves csv.names to the global environment
}





