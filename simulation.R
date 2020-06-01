# =============================================================================#
# NOSOCOMIAL TRANSMISSION
# Simulate data
# =============================================================================#

# ============================#
# Parameters/Data
# ============================#
T <- 3*30                      # Study period (days)

S_p <- rep(100, T)             # Number of susceptible patients at time t
S_hcw <- rep(30, T)            # Number of susceptible HCWs at time t

delta <- rep(0.2, T)           # Proportion of nosocomial infections, infected s days ago who are discharged, isolated or died
gamma <- rep(0.1, T)           # Proportion of nosocomial HCW infections, infected s days ago who recover (and thus are immune)

gen_shape <- 2.826             # shape parameter for generation time distribution (Feretti et al)
gen_scale <- 5.665             # scale parameter for generation time distribution (Feretti et al)

# Transmission parameters (assume to be known in simulation)
f_hcw_p <- 0.001               # known infected patient to HCW
f_hcw_pp <- 0.001              # unknown infected patient to HCW
f_hcw_hcw <- 0.001             # HCW to HCW
f_p_p <- 0.001                 # Known infected patient to patient
f_p_pp <- 0.001                # Unkown infected patient to patient
f_p_hcw <- 0.001               # HCW to susceptible patient

# ============================#
# Simulated data
# ============================#
# Health-care workers
I_hcw0 <- 10                           # Intial number of infected HCWs at t = 1
I_hcw <- matrix(rep(0,T*T), ncol=T)    # Number of infected HCWs at time t who were infected (s-1) days ago and are still working
I_hcw[1,1] <- I_hcw0
R_hcw <- rep(0,T)                      # Number of immune HCWs at time t
w <- rep(0, T)                         # Number of immune HCWs returning to work at time t

# Patients
I_pU0 <- 20
I_pU <- matrix(rep(0, T*T), ncol=T)    # Number of unknown (unisolated) infected patients at time t who were infected s days ago
I_pU[1,1] <- I_pU0
I_p <- rep(10, T)                       # Number of (isolated) infected patients in hospital at time t
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
  inf_hcw[t] <- sum(gen_time*I_hcw[,t-1])
  # Infectivity from patients
  inf_p[t] <- sum(gen_time*I_pU[,t-1])
  # Probability of infection for HCWs at time t
  p_hcw[t] = f_hcw_hcw*inf_hcw[t] + f_hcw_pp*inf_p[t] + f_hcw_p*I_p[t]
  # Probability of infection for patients at time t
  p_p[t] = f_p_hcw*inf_hcw[t] + f_p_pp*inf_p[t] + f_p_p*I_p[t]
  
  
  # Number of newly infected HCWs at time t
  I_hcw[1,t] = rbinom(1, size=S_hcw[t], prob=p_hcw[t])
  # Number of newly infected patients at time t
  I_pU[1,t] = rbinom(1, size=S_p[t], prob=p_p[t])
  
  for(s in 2:t){
    # Number of infected HCWs at time t who got infected s-1 days ago (note that stan starts counting at 1)
    I_hcw[s,t] <- (1-gamma[s-1])*I_hcw[s-1,t-1]
    # Number of unknown infected patients at time t who got infected s-1 days ago 
    # = Number of unknown infected patients a day before - those that are discharged, isolated, or die
    I_pU[s,t] = (1-delta[s-1])*I_pU[s-1,t-1]
  }
  # Number of immune HCWs at time t
  w[t] <- sum(gamma[1:(t-1)]*I_hcw[2:t,t])
  R_hcw[t] = R_hcw[t-1] + w[t]
  
  # Number of non-cohorted patients
  N_ncp[t] <- S_p[t] + I_pU[1,t]
}



