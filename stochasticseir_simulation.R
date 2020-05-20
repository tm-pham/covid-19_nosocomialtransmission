# =================================================================== #
# Simulated Data 
# =================================================================== #
# Discrete time stochastic SEIR model
# Note: 
# Infection of HCWs from community is still missing
# Death of patients is still missing

# Necessary packages
library(dplyr)
library(tidyr)
library(ggplot2)
# File paths
setwd("~/covid-19/nosocomialtransmission/stochasticepidemicmodel/stan/")
source("stochasticseir_functions.R")


set.seed(12345)
# Time frame
t <- 20

# Epidemiological parameters
gamma_A <- 1/7  # recovery rate for asymptomatic infected
gamma_S <- 1/14  # recovery rate for symptomatic infected 
alpha_1 <- 1/4  # 1/duration of latent period
alpha_2 <- 1/2  # Rate from presymptomatic to symptomatic
eta <- 0.0011   # eta=gamma_S*f/(1-f) where f is the case-fatality rate 1.6%
pA <- 0.4       # proportion of asymptomatics
f <- 0.5        # rate at which symptomatic covid cases arrive at the hospital
  
# Infectiousness function (from stochasticseir_functions.R)
# Weibull with shape=2.826, scale=5.665 according to Ferretti et al
max_gen <- 14
x  <- seq(1,max_gen,by=1)
gen_shape <- 2.826; gen_scale <- 5.665
R <- 2
beta <- R*gen.time(x,shape=gen_shape,scale=gen_scale)
# Contact rates per day
c_pat_pat <- 5  # contact between patients
c_pat_hcw <- 5  # contact between patient and hcw
c_hcw_hcw <- 20
c_hcw_pat <- 20

# Dispersion parameter for beta-binomial distribution (infection process)
# see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6027774/
pDshape=1; pDrate=1
pDis <- rgamma(1,shape=pDshape,rate=pDrate)

# ============================= #
# Patient population
# ============================= #
N_p <- 600 # Maximum number of hospital beds
p_b <- 0.8 # Percentage of occupied beds
# Initial values
E0_p <- 10
S0_p <- p_b*N_p-E0_p

pat_dataset <- create.data.template(N_p,S0_p,E0_p,0,0,0,0,EE=ceiling((1-p_b)*N_p),t=t)
pat_data <- pat_dataset$data
pat_wdata <- pat_dataset$wdata

# Store number of patients that were infected 1,2,...,max_gen days ago
pat_past_inf <- c(E0_p,rep(0,max_gen))
pat_data[1:E0_p,3] <- 1 # Time of infection for first infected individuals (assume=1)


# ============================= #
# HCW population
# ============================= #
N_hcw <- 200
E0_hcw <- 1
S0_hcw <- N_hcw-E0_hcw
hcw_dataset <- create.data.template(N_hcw,S0_hcw,E0_hcw,0,0,0,0,EE=0,t=t)
hcw_data <- hcw_dataset$data
hcw_wdata <- hcw_dataset$wdata

# Store number of HCWs infected 1,2,...,max_gen days ago
# Note that symptomatic HCWs get isolated and are not included here
hcw_past_inf <- c(E0_hcw,rep(0,max_gen))
pat_data[1:E0_hcw,3] <- 1


# ============================= #
# Run through the simulation:
# ============================= #
for(i in 1:t){
  print(paste("i=",i))
  # Arrival of symptomatic covid patients
  ind_free_beds <- which(pat_wdata[i,]=='EE')
  num_free_beds <- length(ind_free_beds)
  new_covid_arrivals <- rpois(1,f)
  if(num_free_beds>0 && new_covid_arrivals>0){
    pat_wdata[i,ind_free_beds[1:min(new_covid_arrivals,num_free_beds)]] <- 'I_S'
  }
  
  # Current state of individuals
  ind_susc_pat <- which(pat_wdata[i,]=='S')
  ind_exp_pat <- which(pat_wdata[i,]=='E')
  ind_presymptomatic_pat <- which(pat_wdata[i,]=='I_P')
  ind_symptomatic_pat <- which(pat_wdata[i,]=='I_S')
  ind_asymptomatic_pat <- which(pat_wdata[i,]=='I_A')
  
  S_p <- length(ind_susc_pat)
  E_p <- length(ind_exp_pat)
  I_pP <- length(ind_presymptomatic_pat)
  I_pS <- length(ind_symptomatic_pat)
  I_pA <- length(ind_asymptomatic_pat)
  R_p <- sum(as.numeric(pat_wdata[i,]=='R',T))
  
  ind_susc_hcw <- which(hcw_wdata[i,]=='S')
  ind_exp_hcw <- which(hcw_wdata[i,]=='E')
  ind_presymptomatic_hcw <- which(hcw_wdata[i,]=='I_P')
  ind_symptomatic_hcw <- which(hcw_wdata[i,]=='I_S')
  ind_asymptomatic_hcw <- which(hcw_wdata[i,]=='I_A')
  
  S_hcw <- length(ind_susc_hcw)
  E_hcw <- length(ind_exp_hcw)
  I_hcwP <- length(ind_presymptomatic_hcw)
  I_hcwS <- length(ind_symptomatic_hcw)
  I_hcwA <- length(ind_asymptomatic_hcw)
  R_hcw <- sum(as.numeric(hcw_wdata[i,]=='R',T))
  

  # Number of patients participating in the contact process 
  N_p <- S_p + E_p + I_pP + I_pA
  # Number of HCWs participating in the contact process 
  N_hcw <- S_hcw + E_hcw + I_hcwP + I_hcwA + R_hcw
       
  # How many of the patients that are present and infectious were infected 1,2,..,14 days ago
  days <- max(i-max_gen+1,1):i
  present_pat <- which(pat_wdata[i,]%in%c('I_A','I_P','I_S'))
  pat_past_inf <- c(sapply(days, function(x) length(intersect(which(pat_data[,3]==x),present_pat))),rep(0,max_gen-length(days)))
  
  # How many of the HCWs that are present and infectious were infected 1,2,..,14 days ago
  present_hcw <- which(hcw_wdata[i,]%in%c('I_A','I_P','I_S'))
  hcw_past_inf <- c(sapply(days, function(x) length(intersect(which(hcw_data[,3]==x),present_hcw))),rep(0,max_gen-length(days)))
  
  
  # ============================ #
  # Force of infection
  # ============================ #
  # Infectiousness from infected patients who were infected 1,2,...,max_gen days ago
  beta_by_pat <- as.numeric(pat_past_inf%*%beta)/N_p
  # Infectiousness from infected HCW who were infected 1,2,...,max_gen days ago
  beta_by_hcw <- as.numeric(hcw_past_inf%*%beta)/N_hcw
  
  # Force of infection experienced by susceptible patients
  foi_pat <- c_pat_pat*beta_by_pat + c_pat_hcw*beta_by_hcw
  prob_infected_pat <- 1-exp(-foi_pat)
  # Force of infection experienced by HCWs
  foi_hcw <- c_hcw_pat*beta_by_pat + c_hcw_hcw*beta_by_hcw
  prob_infected_hcw <- 1-exp(-foi_hcw)
  
  # ============================ #
  # TRANSITIONS FOR PATIENTS
  # ============================ #
  # Newly exposed (latent) patients
  print("Infecting patients.")
  new_exposed_pat <- rbbinom(1,prob=prob_infected_pat,
                             k=pDis,
                             size=S_p)
  # Update wdata (change respective entries to 'E")
  ind_new_exposed_pat <- ind_susc_pat[sample(length(ind_susc_pat), new_exposed_pat)]
  pat_wdata[(i+1):t,ind_new_exposed_pat] <- 'E'
  pat_data[ind_new_exposed_pat,3] <- i+1 # Time of infection
  
  print("Exposed to presymptomatics and asymptomatics.")
  # Transition from exposed to pre/asymptomatic patients
  new_presymptomatic_pat <- sum(rbinom(floor((1-pA)*E_p), 1, 1-exp(-alpha_1)))
  new_asymptomatic_pat <- sum(rbinom(ceiling(pA*E_p), 1, 1-exp(-alpha_1)))
  # Update wdata (change respective entries to 'I_A', 'I_P')
  len_exp_pat_s <- floor((1-pA)*length(ind_exp_pat))
  len_exp_pat_a <- ceiling(pA*length(ind_exp_pat))
  ind_exp_pat_s <- ind_exp_pat[min(1,len_exp_pat_s):len_exp_pat_s]
  ind_exp_pat_a <- setdiff(ind_exp_pat,ind_exp_pat_s)
  ind_new_presymptomatic_pat <- ind_exp_pat_s[sample(len_exp_pat_s, max(new_presymptomatic_pat,0))]
  ind_new_asymptomatic_pat <- ind_exp_pat_a[sample(len_exp_pat_a, max(new_asymptomatic_pat,0))]
  pat_wdata[(i+1):t,ind_new_presymptomatic_pat] <- 'I_P'
  pat_data[ind_new_presymptomatic_pat,2] <- i+1 # Time of infectiousness 
  pat_wdata[(i+1):t,ind_new_asymptomatic_pat] <- 'I_A'
  pat_data[ind_new_asymptomatic_pat,2] <- i+1 # Time of infectiousness
  
  print("Presymptomatic to symptomatic.")
  # Transition from presymptomatic to symptomatic patients
  new_symptomatic_pat <- sum(rbinom(I_pP, 1, 1-exp(-alpha_2)))
  ind_new_presymptomatic_pat <- ind_presymptomatic_pat[sample(length(ind_presymptomatic_pat), max(new_symptomatic_pat,0))]
  pat_wdata[(i+1):t,ind_new_presymptomatic_pat] <- 'I_S'
  
  print("Symptomatic and asymptomatic to recovered.")
  # Transition from infectious  to recovered patients
  new_recovered_pat_a <-sum(rbinom(I_pA, 1, 1-exp(-gamma_A)))
  new_recovered_pat_s <- sum(rbinom(I_pS, 1, 1-exp(-gamma_S)))
  ind_new_recovered_pat_a <- ind_asymptomatic_pat[sample(length(ind_asymptomatic_pat),max(new_recovered_pat_a,0))]
  ind_new_recovered_pat_s <- ind_symptomatic_pat[sample(length(ind_symptomatic_pat),max(new_recovered_pat_s,0))]
  # Recovered patients are discharged and immediately replaced by susceptibles
  pat_wdata[(i+1):t,ind_new_recovered_pat_a] <- 'S'
  pat_wdata[(i+1):t,ind_new_recovered_pat_s] <- 'S'

  # ============================ #
  # TRANSITIONS FOR HCWS
  # ============================ #
  # Newly exposed (latent) HCWs
  print("Infecting HCWs.")
  new_exposed_hcw <- rbbinom(1,prob=prob_infected_hcw,
                             k=pDis,
                             size=S_hcw)
  # Update wdata
  ind_new_exposed_hcw <- ind_susc_hcw[sample(length(ind_susc_hcw), max(new_exposed_hcw,0))]
  hcw_wdata[(i+1):t,ind_new_exposed_hcw] <- 'E'
  hcw_data[ind_new_exposed_hcw,3] <- i+1 # Time of infection
  
  print("Exposed to presymptomatics and asymptomatics.")
  # Transition from exposed to asymptomatic/presymptomatic HCW
  new_asymptomatic_hcw <- sum(rbinom(floor(pA*E_hcw), 1, 1-exp(-alpha_1)))
  new_presymptomatic_hcw <- sum(rbinom(floor((1-pA)*E_hcw), 1, 1-exp(-alpha_1)))
  # Update wdata
  len_exp_hcw_s <- floor((1-pA)*length(ind_exp_hcw))
  len_exp_hcw_a <- ceiling(pA*length(ind_exp_hcw))
  ind_exp_hcw_s <- ind_exp_hcw[min(1,len_exp_hcw_s):len_exp_hcw_s]
  ind_exp_hcw_a <- setdiff(ind_exp_hcw,ind_exp_hcw_s)
  ind_new_presymptomatic_hcw <- ind_exp_hcw_s[sample(len_exp_hcw_s, max(new_presymptomatic_hcw,0))]
  ind_new_asymptomatic_hcw <- ind_exp_hcw_a[sample(len_exp_hcw_a, max(new_asymptomatic_hcw,0))]
  hcw_wdata[(i+1):t,ind_new_presymptomatic_hcw] <- 'I_P'
  hcw_data[ind_new_presymptomatic_hcw,2] <- i+1 # Time of infectiousness 
  hcw_wdata[(i+1):t,ind_new_asymptomatic_hcw] <- 'I_A'
  hcw_data[ind_new_asymptomatic_hcw,2] <- i+1 # Time of infectiousness
  
  print("Presymptomatic to symptomatic.")
  # Transition from presymptomatic to symptomatic HCW
  new_symptomatic_hcw <- sum(rbinom(I_hcwP, 1, 1-exp(-alpha_2)))
  ind_new_presymptomatic_hcw <- ind_presymptomatic_hcw[sample(length(ind_presymptomatic_hcw), max(new_symptomatic_hcw,0))]
  hcw_wdata[(i+1):t,ind_new_presymptomatic_hcw] <- 'I_S'
  
  print("Symptomatic and asymptomatic to recovered.")
  # Transition from aasymptomatic/symptomatic to recovered HCW
  new_recovered_hcw_a <- sum(rbinom(I_hcwA, 1, 1-exp(-gamma_A)))
  new_recovered_hcw_s <- sum(rbinom(I_hcwS, 1, 1-exp(-gamma_S)))
  ind_new_recovered_hcw_a <- ind_asymptomatic_hcw[sample(length(ind_asymptomatic_hcw),max(new_recovered_hcw_a,0))]
  ind_new_recovered_hcw_s <- ind_symptomatic_hcw[sample(length(ind_symptomatic_hcw),max(new_recovered_hcw_s,0))]
  hcw_wdata[(i+1):t,ind_new_recovered_hcw_a] <- 'R'
  hcw_wdata[(i+1):t,ind_new_recovered_hcw_s] <- 'R'
}
pat_wdata <- pat_wdata[1:t,]
hcw_wdata <- hcw_wdata[1:t,]

# Plot detected infected patients
no_symp_pat_per_day <- apply(pat_wdata[,-1],1,function(x) sum(as.numeric(x=='I_S')))
data_symp_pat_per_day <- as.data.frame(cbind(time=1:t,I_S=no_symp_pat_per_day))

ggplot(data_symp_pat_per_day, aes(x=time,y=I_S))+
  geom_point() + 
  theme_bw() + 
  ylab("Number of symptomatic patients")

# Plot detected infected HCWs
no_symp_hcw_per_day <- apply(hcw_wdata[,-1],1,function(x) sum(as.numeric(x=='I_S')))
data_symp_hcw_per_day <- as.data.frame(cbind(time=1:t,I_S=no_symp_hcw_per_day))

ggplot(data_symp_hcw_per_day, aes(x=time,y=I_S))+
  geom_point() + 
  theme_bw() + 
  ylab("Number of symptomatic HCWs")
