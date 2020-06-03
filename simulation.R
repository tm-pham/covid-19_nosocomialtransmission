# =============================================================================#
# NOSOCOMIAL TRANSMISSION
# Simulate data
# =============================================================================#

# Libraries
library("TailRank") # for beta binomial distribution


# ============================#
# Parameters/Data
# ============================#
T <- 3*30                      # Study period (days)

S_p <- rep(100, T)             # Number of susceptible patients at time t
S_hcw <- rep(30, T)            # Number of susceptible HCWs at time t

delta <- rep(0.2, T)           # Proportion/Probability of nosocomial infections, infected s days ago who are discharged, isolated or died
gamma <- rep(0.1, T)           # Proportion/Probability of nosocomial HCW infections, infected s days ago who recover (and thus are immune)

gen_shape <- 2.826             # shape parameter for generation time distribution (Feretti et al)
gen_scale <- 5.665             # scale parameter for generation time distribution (Feretti et al)

disp_inf <- 0.01                 # Dispersion parameter for beta-binomial distribution of infection process
disp_obs <- 0.01                 # Dispersion parameter for beta-binomial distribution of observation process
p_p_obs <- 1/3                   # Proportion of observed patient infections (CO-CIN)
alpha_obs <- disp_obs/(1-p_p_obs)  # Parameter for beta-binomial distribution of observation process
beta_obs <- disp_obs/p_p_obs       # Parameter for beta-binomial distribution of observation process

p_hcw_obs <- 0.6                 # Parameter of observed HCW infections 


# Transmission parameters (assume to be known in simulation)
f_hcw_p <- 0.001               # known infected patient to HCW
f_hcw_pp <- 0.001              # unknown infected patient to HCW
f_hcw_hcw <- 0.001             # HCW to HCW
f_p_p <- 0.001                 # Known infected patient to patient
f_p_pp <- 0.001                # Unkown infected patient to patient
f_p_hcw <- 0.001               # HCW to susceptible patient

# Delay distribution
# First entry in prob_delay corresponds to delay=0
prob_delay <- c(1,rep(0,19))
n <- 100
cum_prob_delay <- pgamma(1:n, shape = 0.811, rate = 0.064)
prob_delay <- cum_prob_delay-c(0,cum_prob_delay[1:(n-1)])
onset_to_discharge_cum_distr <- distr.onset.to.discharge(prob_los,prob_inc)$cum_distr
len_delay <- min(length(prob_delay),length(onset_to_discharge_cum_distr)-1)
# Counterfactual delay distribution (delay from symptom onset to study enrolment)
# that would have occurred if there was no discharge
cf_delay_distr <- prob_delay[1:len_delay]/(1-onset_to_discharge_cum_distr[1:len_delay])
cf_delay_distr <- cf_delay_distr/sum(cf_delay_distr)


# ============================#
# Simulated data
# ============================#
# Health-care workers
I_hcw <- matrix(c(10,rep(0,T*T-1)), ncol=T)    # Number of infected HCWs at time t who were infected (s-1) days ago and are still working
I_hcwU <- matrix(c(10,rep(0,T*T-1)), ncol=T)   # Number of unknown infected HCWs at time t who got infected s-1 days ago
I_hcwS <- matrix(c(1,rep(0,T*T-1)), ncol=T)    # Number of unknown infected HCWs who got infected s-1 days ago and develop symptoms at time t 
I_hcwR <- matrix(c(2,rep(0,T*T-1)), ncol=T)    # Number of symptomatic HCWs who got infected s-1 days ago and recover at time t
R_hcw <- rep(2,T)                              # Number of immune HCWs at time t
w <- rep(0, T)                                 # Number of immune HCWs returning to work at time t

# Patients
I_pU <- matrix(c(20,rep(0, T*T-1)), ncol=T)    # Number of unknown (unisolated) infected patients at time t who were infected s days ago
I_pUS <- matrix(rep(0, T*T), ncol=T)   # Number of unisolated infected patients who were infected s days ago and developed symptoms
I_pUD <- matrix(rep(0, T*T), ncol=T)   # Number of unknown infected patients eligible for detection at time t who got infected s-1 days ago
                                       # This accounts for the delay from symptom onset till detection
I_p <- rep(10, T)                      # Number of (isolated) infected patients in hospital at time t
                                       # initialize with number of severly infected patients arriving from community 
N_ncp <- rep(0, T)                     # Number of non-cohorted patients at time t

# Probability of infection
p_hcw <- rep(0,T)                      # Probability of infection for HCWs at time t
p_p <- rep(0,T)                        # Probability of infection for patients at time t
inf_hcw <- rep(0,T)                    # Infectivity from HCWs (densitiy dependent)
inf_p <- rep(0,T)                      # Infectivity from patients (densitiy dependent)

# Generation time distribution (Ferretti et al, 2020)
gen_time <- dweibull(seq(1,T,by=1),shape=gen_shape, scale=gen_scale)

for(t in 2:T){
  # cumulative infectivity from HCWs
  inf_hcw[t] <- sum(gen_time*I_hcwU[,t-1])
  # Infectivity from patients
  inf_p[t] <- sum(gen_time*I_pU[,t-1])
  # Probability of infection for HCWs at time t
  p_hcw[t] = 1-exp(-f_hcw_hcw*inf_hcw[t] - f_hcw_pp*inf_p[t] - f_hcw_p*I_p[t])
  # Probability of infection for patients at time t
  p_p[t] = 1-exp(-f_p_hcw*inf_hcw[t] - f_p_pp*inf_p[t] - f_p_p*I_p[t])
  
  # Number of newly infected HCWs at time t
  # Beta binomial
  alpha_hcw <- disp_inf/(1-p_hcw)  # Parameter for beta-binomial distribution
  beta_hcw <- disp_inf/(1-p_hcw)   # Parameter for beta-binomial distribution
  I_hcwU[1,t] <- rbb(1, N=S_hcw[t], u=alpha_hcw, v=beta_hcw)
  # Alternative: Binomial: I_hcw[1,t] = rbinom(1, size=S_hcw[t], prob=p_hcw[t]) 

  # Number of newly infected patients at time t
  # Beta binomial
  alpha_p <- disp_inf/(1-p_p)    # Parameter for beta-binomial distribution
  beta_p <- disp_inf/(1-p_p)     # Parameter for beta-binomial distribution
  I_pU[1,t] <- rbb(1, N=S_p[t], u=alpha_p, v=beta_p) 
  # Alternative: Binomial: I_pU[1,t] = rbinom(1, size=S_p[t], prob=p_p[t]) 
  
  
  for(s in 2:t){
    # (note that R starts counting at 1)
    # Number of unknown infected HCWs who got infected s-1 days ago and develop symptoms at time t 
    hcw_symp <- rbinom(1, size=I_hcwU[s-1,t-1], prob=inc_distr[s])
    I_hcwS[s,t] <- I_hcwS[s-1,t-1] + hcw_symp
    # Remaining unknown infected HCWs at time t who got infected s-1 days ago
    I_hcwU[s,t] <- I_hcwU[s-1,t-1]-hcw_symp
    # Number of symptomatic HCWs who got infected s-1 days ago and recover at time t
    # Assume that they are immediately isolated
    hcw_recover <- rbinom(1, size=I_hcwS[s-1,t-1], prob=1-gamma[s-1])
    I_hcwR[s,t] <- hcw_recover
    I_hcwS[s,t] <- I_hcwS[s,t] - hcw_recover

    # Number of unknown infected patients that develop symptoms at time t who got infected s-1 days ago 
    I_pUS[s,t] <- rbinom(1, size=I_pU[s-1,t-1], prob=inc_distr[s])
    # Number of unknown infected patients eligible for detection at time t who got infected s-1 days ago
    # Delay from symptom onset till detection
    I_pUD[s,t] <- rbinom(1, size=I_pUS[s-1,t-1], prob=cf_delay_distr[s])
    
    # Number of known infected patients at time t who got infected s-1 days ago
    # Beta-binomial distribution for observation process
    I_p[t] <- I_p[t] +  rbb(1, N=I_pUD[s-1,t-1], u=alpha_obs, v=beta_obs)
    # Number of unknown infected patients at time t who got infected s-1 days ago 
    # = Number of unknown infected patients a day before - those that are discharged, isolated, or die
    I_pU[s,t] <- rbinom(1, size=I_pU[s-1,t-1], prob=1-delta[s-1])
  }
  # Number of immune HCWs at time t
  R_hcw[t] = R_hcw[t-1] + sum(I_hcwR[2:t,t])
  
  # Number of non-cohorted patients
  N_ncp[t] <- S_p[t] + I_pU[1,t]
}


