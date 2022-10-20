
source("General_Scripts/Main_Sim_Functions.R")


create.survey.days<-function(design_total){## Function for creating names to save CMR csv's (Works from 1-99)
  designs_list<-vector(mode = "list", length = design_total) #Creates an empty list for storing all the design optios
  designs_list<<-list(rep(NA, design_total)) #...Not sure if this has a purpose?
  for (i in 2:(design_total+1)){ #Starts as 2 as everything breaks if there's only 1 survey day and I don't have time to fix that
    designs_list[[i-1]]<<-
      seq(from = 1, to = 61, length.out = i) %>% #Evenly spaces survey days according to length.out (number of desired days)
      round() #Rounds to whole numbers, as days are whole numbers
  }
}







##############################################################
##############################################################
#########-------- Simulation -> Inference -----------#########
##############################################################
##############################################################

loop.recaptures.inf = function(n = 100, #number of initial individuals
                               uninf = 0.003, 
                               fecundity = 780, adult.age = 3, off.mortality = 0.980, ovi.sites = 20,
                               age.classes = 4, alpha = c(6.374251, 6.428781, 6.706079, 6.706079), 
                               beta = c(-0.5162382, -0.5408142, -0.6679783, -0.6679783),
                               max.age = 13,
                               years = 3, 
                               det.prob = 0.6, survey.days = c(30, 60, 90),
                               inf = 0.007,
                               init.years = 5,
                               survey.years=2,
                               save.csv = TRUE,
                               path = "Chapter_1_Inf_Dyn/results/") {
  
  
  init.list<-list(c(phi.alpha = 6, phi.beta = 0, p.I.c = 0.5, p.S.c = 0.5,
                    p.I.D.c = 0.5, p.S.D.c = 0.5, b.I.c = 0.003, b.S.c = 0.05), 
                  c(pphi.alpha = 6, phi.beta = 0, p.I.c = 0.5, p.S.c = 0.5,
                    p.I.D.c = 0.5, p.S.D.c = 0.5, b.I.c = 0.003, b.S.c = 0.05),
                  c(phi.alpha = 6, phi.beta = 0, p.I.c = 0.5, p.S.c = 0.5,
                    p.I.D.c = 0.5, p.S.D.c = 0.5, b.I.c = 0.003, b.S.c = 0.05)) #Setting priors for the inference model. General.
  
  
  
  #Create .csv  
  mark_recapture(n, 
                 uninf = uninf, inf = inf, #Infection parameters
                 fecundity = fecundity, adult.age = adult.age, alpha = alpha, beta = beta, ovi.sites = ovi.sites, #Reproduction parameters
                 age.classes = age.classes,  max.age = max.age,  #Death parameters
                 init.years = init.years, survey.years = survey.years, save.csv = TRUE, survey.days = survey.days,
                 det.prob = det.prob,
                 csv.name = paste0(path, "CMR_csvs/", csv.name)) #Mark-recapture parameters
  
  
  pop_path<-paste(path, "population_summary/pop_summary_", as.character(det_prob), ".csv", sep = "")
  write.csv(pop.history, file = pop_path) #Saves the population history as a csv file
  
  ###---Matt's Model---###
  
  
  ###############################################################
  ###   load in capture history (incl SI) data              ######
  ###############################################################
  
  data <- paste0(path, "CMR_csvs/", csv.name)
  
  CH_SI_data.df <- read.csv(data) # CH = capture history, SI=susceptible infected
  
  head(CH_SI_data.df)
  dim(CH_SI_data.df[-1])
  
  # Note in this data frame 0 = not captured, 1 = captured uninfected, 2 = captured infected
  
  
  # convert all not captured histories to equal 5
  
  CH_SI_data.df[CH_SI_data.df==0] <- 5
  
  
  ########################################################
  ###     get data ready for model                  ######
  ########################################################
  
  # format Tap data for model fitting
  tmpCH = collapse.ch(CH_SI_data.df[-1])[[1]] #creates a data frame of each of the unique capture histories
  
  FR = collapse.ch(CH_SI_data.df[-1])[[2]] #creates a vector of the frequency of each unique capture history
  
  sumf <- apply(tmpCH, 1, get.first) # creates a vector of the first time (occasion) that individuals were captured for each unique capture history 
  
  NCH = nrow(tmpCH)          # number of capture histories 
  Nint = ncol(tmpCH) - 1    # number of sampling intervals
  
  ones <- FR
  
  State.index <- tmpCH 
  
  #############################################################
  ###  create supporting data vectors for use in model   ######
  #############################################################
  
  
  
  
  
  
  ########################################################################################
  ## One final data check /manipulation                                                 ##
  ##  as the model does not allow individuals to only be observed on the last occasion  ##
  ########################################################################################
  
  #so need to do the following:
  
  sumf <- ifelse(sumf==(Nint+1),Nint,sumf) 
  ##################################################################################################################################
  ##### load Model Script constant monthly survival, infection and recovery, detection does not vary with time ###################
  #################################################################################################################################
  
  model_file="General_Scripts/jags_model.txt"
  
  ################################################
  ##### Fit Jolly-Seber model using JagsUI  ######
  ################################################
  ##---define the data for the marginalized model -------
  
  # year.index.data.m <- year.index.data[-1]
  # season <- unique(year.index.data.m)
  
  #------------------------------------- This code block can be used to alter intervalLength for multiple years----#
  # generate interval length vector
  # survey.days <- c(survey.days, 366)
  # intervalLength <- diff(survey.days)
  # intervalLength <- rep(intervalLength, survey.years) 
#---------------------------------------------#
  
  intervalLength <- diff(survey.days)
  
  
  
  jagsM.data <- list(NCH = NCH, Nint = Nint, State.index = State.index,
                     sumf = sumf, FR = FR, ones = ones, 
                     intervalLength = intervalLength)
  
  
  #----identify the parameters to monitor-------------------
  
  params=c('phi.alpha','phi.beta', 'b.I.c', 'b.S.c', 'p.I.c', 'p.S.c')#
  
  #-----fit the model-------------------------
  
  #JM.inits <- function(){list(phi.I.c = c(.7, .7, .7, .7), p.I.D.c= 0, p.S.D.c= 0)}
  
  jags_mod<-jags(data=jagsM.data,               #data object
                 parameters.to.save = params, #names of parameters to monitor
                 inits=init.list,              #function to generate initial values #initfunc #JM.inits,#
                 seed=1492,                   #random seed
                 model.file=model_file,       #file containing JAGS model
                 n.chains=3,                  
                 n.adapt=100, 
                 n.iter = 80000, #burnin for Question 3 test 4
                 n.burnin = 40000, #burnin for Question 3 test 4
                 n.thin=10,
                 DIC = FALSE,
                 #DIC = TRUE,#defualt = True
                 parallel=TRUE)               
  
 
  
  jags_mod # print(jags_mod, digits = 3)
  
  
  filepath<-paste(path, "traceplots/traceplot_det_prob_", as.character(det_prob), ".jpeg", sep = "")
  jpeg(file = filepath,  width = 10, height = 7.5, units = 'in', res = 300)
  #check for convergence
  try(traceplot(jags_mod))
  dev.off()
  
  filepath<-paste(path, "densityplots/density_det_prob_", as.character(det_prob), ".jpeg", sep = "")
  jpeg(file=filepath, width = 10, height = 7.5, units = 'in', res = 300)
  #examine density plots
  try(densityplot(jags_mod)) #parameters=NULL, layout=NULL, ask=
  dev.off()
  
  jags_mod<<-jags_mod #Save jags_mod to global environment
  
  #plot parameter estimates mean and 95% credible interval
  jpeg(paste(path, "boxplots/", paste(as.character(det_prob), "Det Prob, phi.alpha, phi.beta.jpeg"), sep = "")) 
  whiskerplot(jags_mod, parameters=c('phi.alpha','phi.beta'), quantiles=c(0.025,0.975), zeroline=TRUE) #plot survival estimates 
  dev.off()
  jpeg(paste(path, "boxplots/", paste(as.character(det_prob), "Det Prob, b.I.c, b.S.c.jpeg"), sep = "")) 
  whiskerplot(jags_mod, parameters=c('b.I.c', 'b.S.c'), quantiles=c(0.025,0.975), zeroline=TRUE) #plot infection and recovery estimates 
  dev.off()
  jpeg(paste(path, "boxplots/", paste(as.character(det_prob), "Det Prob, p.I.c, p.S.c.jpeg"), sep = "")) 
  whiskerplot(jags_mod, parameters=c('p.I.c', 'p.S.c'), quantiles=c(0.025,0.975), zeroline=TRUE) #plot detection estimates
  dev.off()
  
  
  
  ####################################################################
  ##### Define model output name and location  as an Rdata file ######
  ####################################################################
  fitted_model_filename<-paste(paste(path, "estimates_summaries/", as.character(det_prob), sep = ""), "Det Probability")
  
  
  
  save.image(paste(fitted_model_filename)) #change name of fitted model object as appropriate
  
  
  
  fitted_model_csv_filename<-paste(paste(path, "estimates_summaries/", as.character(det_prob), sep = ""), "det probability.csv")
  filename<-paste(fitted_model_csv_filename)
  write.table(jags_mod$summary, file = filename, sep = ",", col.names = NA, qmethod = "double") 
  
  
  save.image(paste(fitted_model_filename)) #change name of fitted model object as appropriate
  
  r.hat <- (jags_mod$summary)[1:6, 'Rhat']  #Prints a message if R.hat is greater than 1.1
  if(any(r.hat>1.1)==TRUE){
    print(paste("design index", as.character(design_index), "det prob", as.character(inf), "R.hat greater than 1.1"))
  } else {
    print("Convergence achieved")
  }
  
}




#Extra: Convert yearly probabilities to daily

year.to.day<-function(probs){
  daily.probs=1-((1-probs)^(1/365))
  return(daily.probs)
}






