# =============================================================================#
# NOSOCOMIAL TRANSMISSION
# Simulate data
# 8th June 2020
# adapted from https://github.com/tm-pham/covid-19_nosocomialtransmission/blob/master/simulation.R
# =============================================================================#

setwd("/Users/tm-pham/PhD/covid-19/nosocomialtransmission/stochasticepidemicmodel/stan")
# Libraries
# For beta binomial distribution
if("TailRank"%in%installed.packages()) library("TailRank")
if("extraDistr"%in%instsalled.packages()) library("extraDistr") 

# ============================#
# Parameters
# ============================#
max_time <- 3*30                   # Study period (days)
S <- 14                            # Maximum number of days to be infectious
delta <- rep(0.2, max_time)        # Proportion/Probability of nosocomial infections who are discharged or died 

hcw_isolation_period <- 7          # Number of days HCWs isolate themselves on developing symptoms consistent with COVID
                                   # for now assume a fixed isolation duration
gen_shape <- 2.826                 # shape parameter for generation time distribution (Feretti et al)
gen_scale <- 5.665                 # scale parameter for generation time distribution (Feretti et al)
disp_inf <- 0.01                   # Dispersion parameter for beta-binomial distribution of infection process
disp_obs <- 0.01                   # Dispersion parameter for beta-binomial distribution of observation process
p_p_obs <- 1/3                     # Proportion of observed patient infections (CO-CIN)
alpha_obs <- disp_obs/(1-p_p_obs)  # Parameter for beta-binomial distribution of observation process
beta_obs <- disp_obs/p_p_obs       # Parameter for beta-binomial distribution of observation process

# Transmission parameters (assume to be known in simulation)
f_hcw_pp <- 0.001                  # unknown infected patient to HCW
f_hcw_hcw <- 0.001                 # HCW to HCW
f_p_pp <- 0.001                    # Unkown infected patient to patient
f_p_hcw <- 0.001                   # HCW to susceptible patient

# Probability distribution for incubation period
# Source: Ferretti et al (2020): https://science.sciencemag.org/content/368/6491/eabb6936.full
prob_asymptomatic <- 0.4
meanlog<-1.644; sdlog<-0.363
inc_distr <- (1-prob_asymptomatic)*dlnorm(1:max_time,meanlog,sdlog)
# inc_distr should give probablity of becoming symptomatic on each day after infection (given still infected)
# these probabilities don’t need to sum to one since not all people will develop symptom

# Probability distribution for LOS
meanlos <- 7
cum_prob_los <- pexp(1:max_time,1/meanlos)
prob_los <- cum_prob_los-c(0,cum_prob_los[1:(max_time-1)])


# ============================#
# Data
# ============================#
# Health-care workers
S_hcw <- rep(30, max_time)                                                      # Number of susceptible HCWs at time t
I_hcwU <- matrix(c(10,rep(0,max_time*max_time-1)), ncol=max_time)               # Number of unknown infected HCWs at time t who got infected s-1 days ago
new_symptomatic_hcw <- matrix(c(1,rep(0,max_time*max_time-1)), ncol=max_time)   # Number of unknown infected HCWs who got infected s-1 days ago and develop symptoms at time t
R_hcw <- rep(0,max_time)                                                        # Number of immune HCWs at time t
isolated_hcw <- rep(10,max_time)                                                 # Number of isolated HCWs at time t
obs_nosocomial<-rep(10,max_time)                                                 # Observed nosocomial infections

# Patients
S_p <- rep(100, max_time)                                                       # Number of susceptible patients at time t
I_pU <- matrix(c(20,rep(0, max_time*max_time-1)), ncol=max_time)                # Number of unknown (unisolated) infected patients at time t who were infected s days ago
new_symptomatic_pat <- matrix(rep(1, max_time*max_time), ncol=max_time)         # Number of unisolated infected patients who were infected s days ago and developed symptoms
N_ncp <- rep(S_p[1]+I_pU[1,1], max_time)                                    # Number of non-cohorted patients at time t
# Assume the following is observed
I_p <- rep(10, max_time)                                                        # Number of isolated infected patients at time t

# Probability of infection
p_hcw <- rep(0,max_time)        # Probability of infection for HCWs at time t
p_p <- rep(0,max_time)          # Probability of infection for patients at time t
inf_hcw <- rep(0,max_time)      # Infectivity from HCWs (densitiy dependent)
inf_p <- rep(0,max_time)        # Infectivity from patients (densitiy dependent)

# Generation time distribution (Ferretti et al, 2020)
gen_time <- dweibull(seq(1,max_time,by=1),shape=gen_shape, scale=gen_scale)

# =============================================================================#
# START SIMULATION
# =============================================================================#s
for(t in 2:max_time){
  # Tranmission from nosocomial infections
  # Assume that only unknown (unisolated) infected HCWs and patients are infectious
  # Cumulative infectivity from HCWs
  inf_hcw[t] <- sum(gen_time*I_hcwU[,t-1])
  # Cumulative infectivity from patients
  inf_p[t] <- sum(gen_time*I_pU[,t-1])
  # Probability of infection for HCWs at time t
  p_hcw[t] = 1-exp(-f_hcw_hcw*inf_hcw[t] - f_hcw_pp*inf_p[t] - f_p_p*I_p[t])
  # Probability of infection for patients at time t
  p_p[t] = 1-exp(-f_p_hcw*inf_hcw[t] - f_p_pp*inf_p[t] - f_hcw_p*I_p[t])
  
  # Number of newly infected HCWs at time t (beta-binomial)
  alpha_hcw <- disp_inf/(1-p_hcw)  # Parameter for beta-binomial distribution
  beta_hcw <- disp_inf/(1-p_hcw)   # Parameter for beta-binomial distribution
  
  # Number of newly infected patients at time t (beta-binomial)
  alpha_p <- disp_inf/(1-p_p)     # Parameter for beta-binomial distribution
  beta_p <- disp_inf/(1-p_p)      # Parameter for beta-binomial distribution
  
  # For beta binomial distribution
  if("TailRank"%in%installed.packages()){
    I_hcwU[1,t] <- rbb(1, N=S_hcw[t], u=alpha_hcw, v=beta_hcw)
    I_pU[1,t] <- rbb(1, N=S_p[t], u=alpha_p, v=beta_p)
  }else{
    if("extraDistr"%in%installed.packages()){
      I_hcwU[1,t] <- rbbinom(1, size=S_hcw[t], alpha = alpha_hcw, beta = beta_hcw)
      I_pU[1,t] <- rbbinom(1, size=S_p[t], alpha=alpha_p, beta=beta_p)
    }else{
      print("No package for beta-binomial distribution installed. Aborting simulation...")
      break;
    }
  } 
  

  for(s in 2:S){
    # (note that R starts counting at 1)
    # HEALTH-CARE WORKERS
    # Number of unknown infected HCWs who got infected s-1 days ago and develop symptoms at time t
    hcw_symp <- rbinom(1, size=I_hcwU[s-1,t-1], prob=inc_distr[s])
    new_symptomatic_hcw[s,t] <- hcw_symp
    # Remaining unknown infected HCWs at time t who got infected s-1 days ago
    I_hcwU[s,t] <- I_hcwU[s-1,t-1] - hcw_symp

    # PATIENTS
    # I_pU[s-1,t-1] = number of unknown infected patients at time t-1 who got infected s-1 days ago
    # The following three things can happen to these patients at the next time step
    # 1) they may be discharged or die (or recover). Probability of this is delta[s-1]
    # 2) they may stay on the ward and develop symptoms. Probability of this is (1-delta[s-1])*inc_distr[s]
    # 3) they may stay on the ward and not develop symtpoms. Probability of this is (1-delta[s-1])*(1-inc_distr[s])
    outcomes<-rmultinom(1,I_pU[s-1,t-1],c(delta[s-1],(1-delta[s-1])*inc_distr[s],(1-delta[s-1])*(1-inc_distr[s])))
    I_pU[s,t]<-outcomes[3,1]
    # Number of unisolated infected patients who were infected s days ago and developed symptoms at time t
    new_symptomatic_pat[s,t]<-outcomes[2,1]
  }
  
  # Number of known (isolated) infected patients at time t
  # Beta-binomial distribution for observation process
  if("TailRank"%in%installed.packages()){
    obs_nosocomial[t] <- rbb(1, N=sum(new_symptomatic_pat[1:t,t]), u=alpha_obs, v=beta_obs)
  }else{
    if("extraDistr"%in%installed.packages()) library("extraDistr"){
      obs_nosocomial[t] <- rbbinom(1, size=sum(new_symptomatic_pat[1:t,t]), alpha=alpha_obs, beta=beta_obs)
    }else{
      print("No package for beta-binomial distribution installed. Aborting simulation...")
      break;
    }
  }
  
  # Assume that symptomatic HCWs self-isolate and return to work after 7 days with immunity
  hcw_recover <- ifelse(t>hcw_isolation_period, sum(new_symptomatic_hcw[,t-hcw_isolation_period]), 0)

  # Number of isolated HCWs
  # = Number of isolated HCWs the day before - those that recover and those that are still infected at time t (got infected 1,2,...,t days ago)
  # Assume HCWs are isolated immediately
  isolated_hcw[t] <- isolated_hcw[t-1] - hcw_recover + sum(new_symptomatic_hcw[1:t,t])
  # Number of recovered HCWs
  R_hcw[t] = R_hcw[t-1] + hcw_recover
  # Number of non-cohorted patients
  N_ncp[t] <- S_p[t] + I_pU[1,t]
}


sim_data <- list(T=max_time, 
                 S=S, 
                 N_ncp=N_ncp, 
                 I_p=I_p, 
                 obs_nosocomial=obs_nosocomial, 
                 S_hcw=S_hcw, 
                 hcw_isolation_period=hcw_isolation_period, 
                 delta=delta, 
                 gen_shape=gen_shape, 
                 meanlog=meanlog, 
                 sdlog=sdlog, 
                 prob_asymptomatic=prob_asymptomatic, 
                 p_p_obs=p_p_obs)


saveRDS(sim_data, file="sim_data.RDS")


