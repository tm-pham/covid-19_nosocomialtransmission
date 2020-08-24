# =============================================================================#
# NOSOCOMIAL TRANSMISSION
# Simulate data
# 19th August 2020
# =============================================================================#

# In this simplified version, we assume
# - only one population (namely the patient population)
# - constant infecitous period
# - constant infectiousness over infectious period
# - prob_symptom_onset = fixed daily probability that asymptomatic infections develop symptoms

# setwd("/Users/tm-pham/PhD/covid-19/nosocomialtransmission/stochasticepidemicmodel/stan")
# Libraries

# rbbinom<-function(n, prob, k, size){
#   mtilde <- rbeta(n, k/(1-prob), k/prob)
#   return(rbinom(n, prob=mtilde, size=size))
# }

rbbinom<-function(n, size, alpha, beta){
  mtilde <- rbeta(n, alpha, beta)
  return(rbinom(n, prob=mtilde, size=size))
}


# ============================#
# Parameters
# ============================#
max_time <- 3*30                                                                # Study period (days)
delta <- 0.1                                                                    # Probability of nosocomial infections being discharged or dying in a given time step 
infectious_period <- 7                                                          # Infectious period (assumed to be constant for now)

# For now assume a fixed isolation duratio
disp_inf <- 1000                                                                # Dispersion parameter for beta-binomial distribution of infection process
disp_obs <- 1000                                                                # Dispersion parameter for beta-binomial distribution of observation process
p_p_obs <- 1/3                                                                  # Proportion of observed patient infections (CO-CIN)
alpha_obs <- disp_obs/(1-p_p_obs)                                               # Parameter for beta-binomial distribution of observation process
beta_obs <- disp_obs/p_p_obs                                                    # Parameter for beta-binomial distribution of observation process

# Transmission parameters (assumed to be known in simulation)
# In parameterisation below these represent components of the next generation matrix
# i.e. umber of secondary cases in a given group caused by one infection in a 
# given group if all susceptible
R_pU_p <- 0.8                                                                   # from unkown infected patient to susceptible patient
R_p_p <- 0.4                                                                    # from known infected patient to susceptible patient

# In this simplified model we assume prevalent presymptomatic cases have a 
# fixed daily probablity of developing symptoms
prob_symptom_onset<- 0.1 

# Probability distribution for length-of-stay
# Assume exponential distribution 
meanlos <- 7
cum_prob_los <- pexp(1:max_time,1/meanlos)
prob_los <- cum_prob_los-c(0,cum_prob_los[1:(max_time-1)])


# ============================#
# Data
# ============================#
# Patients
N_ncp <- rep(900, max_time)                                                     # Number of non-cohorted patients at time t
I_p <- c(5, rep(10, max_time-1))                                                # Number of isolated infected patients at time t
obs_nosocomial<-rep(0,max_time)                                                 # Observed nosocomial infections
sum_symp_pat <- c(0,rep(NA, max_time-1))                                        # Number of patients that develop symptoms at time t
discharged_dead_pat <- rep(0, max_time)                                         # Number of patients that are discharged or die at time t

I_pU <- c(5,rep(NA, max_time-1))                                                # Incidence of unknown (unisolated) infected patients at time t 
Prev_pU <- c(I_pU[1],rep(NA, max_time-1))                                       # Prevalence of unknown (unisolated) infected patients 
S_p <- c(N_ncp[1]-Prev_pU[1], rep(NA,max_time-1))                               # Number of susceptible patients at time t

# Initialize vector for probability of infection
p_p <- rep(NA,max_time-1)                                                       # Probability of infection for patients at time t

# =============================================================================#
# START SIMULATION
# =============================================================================#
for(t in 2:max_time){
  # ======================================= #
  # Transmission from nosocomial infections
  # Here assuming a constant infectious period
  # Probability of infection for patients at time t
  p_p[t-1] = 1-exp(-R_pU_p*Prev_pU[t-1]/(N_ncp[t-1]*infectious_period) - R_p_p*I_p[t-1]/(N_ncp[t-1]*infectious_period))
  
  # ======================================= #
  # Infection process
  # Number of newly infected patients at time t (beta-binomial)
  alpha_p <- disp_inf/(1-p_p[t-1])                                              # Parameter for beta-binomial distribution
  beta_p <- disp_inf/p_p[t-1]                                                   # Parameter for beta-binomial distribution
  I_pU[t] <- rbbinom(1, size=S_p[t-1], alpha=alpha_p, beta=beta_p)
  
  # Number of non-cohorted patients
  S_p[t] <- max(S_p[t-1] - I_pU[t], 0)
  
  # The following three things can happen to these patients at the next time step
  # 1) they may be discharged or die (or recover). Probability of this is delta
  # 2) they may stay on the ward and develop symptoms. Probability of this is (1-delta)*prob_symptom_onset
  # 3) they may stay on the ward and not develop symtpoms. Probability of this is (1-delta)*(1-prob_symptom_onset)
  outcomes<-rmultinom(1, Prev_pU[t-1], c(delta,(1-delta)*prob_symptom_onset,(1-delta)*(1-prob_symptom_onset)))
  discharged_dead_pat[t] <- outcomes[1,1]
  sum_symp_pat[t] <- outcomes[2,1]
  Prev_pU[t] <- outcomes[3,1] + I_pU[t]
  
  # Number of newly observed infected patients at time t
  # Beta-binomial distribution for observation process
  obs_nosocomial[t] <- rbbinom(1, size=sum_symp_pat[t], alpha=alpha_obs, beta=beta_obs)
}

sim_data <- list(T=max_time, 
                 N_ncp=N_ncp, 
                 I_p=I_p, 
                 i_pU0 = I_pU[1], 
                 sum_symp_pat0 = sum_symp_pat[1],
                 obs_nosocomial=obs_nosocomial,
                 delta=delta, 
                 infectious_period=infectious_period,
                 prob_symptom_onset = prob_symptom_onset,
                 p_p_obs=p_p_obs)
saveRDS(sim_data, file="sim_data_only_patients_1.RDS")

parameters <- list(R_pU_p = R_pU_p,                                             # from unkown infected patient to susceptible patient
                    R_p_p = R_p_p)                                              # from unknown infected HCW to susceptible patient
saveRDS(parameters, file="parameters_only_patients_1.RDS")
