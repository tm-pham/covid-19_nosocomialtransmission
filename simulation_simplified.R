# =============================================================================#
# NOSOCOMIAL TRANSMISSION
# Simulate data
# 8th June 2020
# adapted from https://github.com/tm-pham/covid-19_nosocomialtransmission/blob/master/simulation.R
# =============================================================================#

# setwd("/Users/tm-pham/PhD/covid-19/nosocomialtransmission/stochasticepidemicmodel/stan")
# Libraries
# For beta binomial distribution
if("extraDistr"%in%installed.packages()) library("extraDistr") 

# ============================#
# Parameters
# ============================#
max_time <- 3*30                         # Study period (days)
S <- 14                                  # Maximum number of days to be infectious
delta <- rep(0.01, max_time)             # Proportion/Probability of nosocomial infections who are discharged or died 

hcw_isolation_period <- 7                # Number of days HCWs isolate themselves on developing symptoms consistent with COVID
hcw_recovery_rate <- 1/7
# for now assume a fixed isolation duration
gen_shape <- 2.826                       # shape parameter for generation time distribution (Feretti et al)
gen_scale <- 5.665                       # scale parameter for generation time distribution (Feretti et al)
disp_inf <- 10                           # Dispersion parameter for beta-binomial distribution of infection process
disp_obs <- 10                           # Dispersion parameter for beta-binomial distribution of observation process
p_p_obs <- 1/3                           # Proportion of observed patient infections (CO-CIN)
alpha_obs <- disp_obs/(1-p_p_obs)        # Parameter for beta-binomial distribution of observation process
beta_obs <- disp_obs/p_p_obs             # Parameter for beta-binomial distribution of observation process
p_hcw_obs <- 0.1                         # Proportion of observed hcw infections (assumed to be observed when hcws admitted)
alpha_obs_hcw <- disp_obs/(1-p_hcw_obs)  # Parameter for beta-binomial distribution of observation process
beta_obs_hcw <- disp_obs/p_hcw_obs       # Parameter for beta-binomial distribution of observation process

# Transmission parameters (assume to be known in simulation)
f_pU_hcw <- 0.007                        # from unknown infected patient to susceptible HCW
f_pU_p <- 0.005                          # from unkown infected patient to susceptible patient
f_p_hcw <- 0.0005                        # from known infected patient to susceptible HCW
f_p_p <- 0.0001                          # from known infected patient to susceptible patient
f_hcw_hcw <- 0.005                       # from unknown infected HCW to susceptible HCW
f_hcw_p <- 0.001                         # from unknown infected HCW to susceptible patient

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
N_hcw <- rep(500, max_time)                                                     # Total number of HCWs at time t
I_hcwU <- c(30,rep(0,max_time-1))                                               # Number of unknown infected HCWs at time t who got infected s-1 days ago
R_hcw <- rep(0,max_time)                                                        # Number of immune HCWs at time t
isolated_hcw <- rep(10,max_time)                                                # Number of isolated HCWs at time t
Prev_hcwU <- rep(I_hcwU[1], max_time)                                           # Prevalene of unknown (unisolated) infected HCWs
sum_symp_hcw <- rep(1, max_time)                                                # Number of patients that develop symptoms at time t
S_hcw <- rep(N_hcw[1]-I_hcwU[1]-R_hcw[1]-isolated_hcw[1]-sum_symp_hcw[1], max_time) # Number of susceptible HCWs at time t


# Patients
N_ncp <- rep(900, max_time)                                                     # Number of non-cohorted patients at time t
I_pU <- c(100,rep(0, max_time-1))                                               # Number of unknown (unisolated) infected patients at time t who were infected s days ago
Prev_pU <- rep(I_pU[1], max_time)                                               # Prevalene of unknown (unisolated) infected HCWs 
discharged_dead_pat <- rep(1, max_time)                                         # Number of patients that are discharged or die at time t
obs_nosocomial<-rep(0,max_time)                                                 # Observed nosocomial infections
obs_hcw_infections<-rep(0,max_time)                                             # Observed  infections in hcws
sum_symp_pat <- rep(0, max_time)                                                # Number of patients that develop symptoms at time t
S_p <- rep(N_ncp[1]-I_pU[1], max_time)                                          # Number of susceptible patients at time t

# Assume the following is observed
I_p <- rep(50, max_time)                                                        # Number of isolated infected patients at time t

# Probability of infection
p_hcw <- rep(0,max_time)        # Probability of infection for HCWs at time t
p_p <- rep(0,max_time)          # Probability of infection for patients at time t
inf_hcw <- rep(0,max_time)      # Infectivity from HCWs (density dependent)
inf_p <- rep(0,max_time)        # Infectivity from patients (density dependent)

# Generation time distribution (Ferretti et al, 2020)
gen_time <- dweibull(seq(1,max_time,by=1),shape=gen_shape, scale=gen_scale)

# =============================================================================#
# START SIMULATION
# =============================================================================#s
for(t in 2:max_time){
  # Tranmission from nosocomial infections
  # Assume that only unknown (unisolated) infected HCWs and patients are infectious
  # Cumulative infectivity from HCWs
  inf_hcw[t-1] <- sum(gen_time[1:(t-1)]*I_hcwU[1:(t-1)])
  # Cumulative infectivity from patients
  inf_p[t-1] <- sum(gen_time[1:(t-1)]*I_pU[1:(t-1)])
  # Probability of infection for HCWs at time t
  p_hcw[t-1] = 1-exp(-f_hcw_hcw*inf_hcw[t-1] - f_pU_hcw*inf_p[t-1] - f_p_hcw*I_p[t-1])
  # Probability of infection for patients at time t
  p_p[t-1] = 1-exp(-f_hcw_p*inf_hcw[t-1] - f_pU_p*inf_p[t-1] - f_p_p*I_p[t-1])
  
  # Number of newly infected HCWs at time t (beta-binomial)
  alpha_hcw <- disp_inf/(1-p_hcw)  # Parameter for beta-binomial distribution
  beta_hcw <- disp_inf/p_hcw       # Parameter for beta-binomial distribution
  
  # Number of newly infected patients at time t (beta-binomial)
  alpha_p <- disp_inf/(1-p_p)      # Parameter for beta-binomial distribution
  beta_p <- disp_inf/p_p           # Parameter for beta-binomial distribution
  
  
  if("extraDistr"%in%installed.packages()){
    I_hcwU[t] <- rbbinom(1, size=S_hcw[t-1], alpha = alpha_hcw, beta = beta_hcw)
    I_pU[t] <- rbbinom(1, size=S_p[t-1], alpha=alpha_p, beta=beta_p)
  }else{
    print("No package for beta-binomial distribution installed. Aborting simulation...")
    break;
  }
  
  # Number of non-cohorted patients
  S_p[t] <- max(S_p[t-1] - I_pU[t], 0)
  
  # HEALTH-CARE WORKERS
  # Number of unknown infected HCWs who develop symptoms at time t
  # Sample from HCW prevalence at t
  sum_symp_hcw[t] <- rbinom(1, size=Prev_hcwU[1:(t-1)], prob=sum(I_hcwU[1:(t-1)]*rev(inc_distr[1:(t-1)])/sum(I_hcwU[1:(t-1)])))

  # Update the prevalence of HCWs at time t
  Prev_hcwU[t] <- max(Prev_hcwU[t-1] - sum_symp_hcw[t] + I_hcwU[t], 0)
  
  
  
  # PATIENTS
  # The following three things can happen to these patients at the next time step
  # 1) they may be discharged or die (or recover). Probability of this is delta[s-1]
  # 2) they may stay on the ward and develop symptoms. Probability of this is (1-delta[s-1])*inc_distr[s]
  # 3) they may stay on the ward and not develop symtpoms. Probability of this is (1-delta[s-1])*(1-inc_distr[s])
  outcomes<-rmultinom(1, Prev_pU[t-1], c(delta[t],(1-delta[t])*sum(I_pU[1:(t-1)]*rev(inc_distr[1:(t-1)])/Prev_pU[t-1]),(1-delta[t])*(1-sum(I_pU[1:(t-1)]*rev(inc_distr[1:(t-1)])/Prev_pU[t-1]))))
  discharged_dead_pat[t] <- outcomes[1,1]
  sum_symp_pat[t] <- outcomes[2,1]
  Prev_pU[t] <- outcomes[3,1] + I_pU[t]
  
  # Number of newly observed infected patients at time t
  # Beta-binomial distribution for observation process
  if("extraDistr"%in%installed.packages()){
    obs_nosocomial[t] <- rbbinom(1, size=sum_symp_pat[t], alpha=alpha_obs, beta=beta_obs)
    obs_hcw_infections[t] <- rbbinom(1, size=sum_symp_hcw[t], alpha=alpha_obs_hcw, beta=beta_obs_hcw)
  }else{
    print("No package for beta-binomial distribution installed. Aborting simulation...")
    break;
  }
  # Symptomatic HCWs that isolated (at home) upon symptom onset and return to work after 7 days with immunity
  hcw_isolated_recover <- ifelse(t>hcw_isolation_period, sum_symp_hcw[t-hcw_isolation_period], 0)
  
  # Number of isolated HCWs
  # = Number of isolated HCWs the day before - those that recover and those that are still infected at time t (got infected 1,2,...,t days ago)
  # Assume HCWs are isolated immediately
  isolated_hcw[t] <- max(isolated_hcw[t-1] - hcw_isolated_recover + sum(sum_symp_hcw[1:t]), 0)
  

  # Number of recovered HCWs
  R_hcw[t] = R_hcw[t-1] + hcw_recover
  # Number of susceptible HCWs
  S_hcw[t] = max(N_hcw[t] - R_hcw[t] - I_hcwU[t] - isolated_hcw[t], 0)
}

sim_data <- list(T=max_time, 
                 S=S, 
                 N_ncp=N_ncp, 
                 I_p=I_p, 
                 i_pU0 = I_pU[1], 
                 sum_symp_pat0 = sum_symp_pat[1],
                 obs_nosocomial=obs_nosocomial,
                 N_hcw=N_hcw, 
                 hcw_isolation_period=hcw_isolation_period, 
                 i_hcwU0 = I_hcwU[1],
                 isolated_hcw0 = isolated_hcw[1],
                 sum_symp_hcw0 = sum_symp_hcw[1],
                 r_hcw0=R_hcw[1],
                 obs_hcw_infections=obs_hcw_infections,
                 delta=delta, 
                 gen_shape=gen_shape,
                 gen_scale=gen_scale,
                 meanlog=meanlog, 
                 sdlog=sdlog, 
                 prob_asymptomatic=prob_asymptomatic, 
                 p_p_obs=p_p_obs,
                 p_hcw_obs=p_hcw_obs)
saveRDS(sim_data, file="sim_data_simplified.RDS")

parameters <- list(f_pU_hcw = 0.01,                   # from unknown infected patient to susceptible HCW
                   f_pU_p = 0.005,                    # from unkown infected patient to susceptible patient
                   f_p_hcw = 0.0005,                  # from known infected patient to susceptible HCW
                   f_p_p = 0.0001,                    # from known infected patient to susceptible patient
                   f_hcw_hcw = 0.01,                  # from unknown infected HCW to susceptible HCW
                   f_hcw_p = 0.001)                   # from unknown infected HCW to susceptible patient
saveRDS(parameters, file="parameters_simplified.RDS")




