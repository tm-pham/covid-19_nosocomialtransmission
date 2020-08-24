# =============================================================================#
# NOSOCOMIAL TRANSMISSION
# Run file for STAN code
# Simplified version - v5
# Correspondonds to simulation_only_patients.R and 
# nosocomialtransmission_only_patients.stan
# =============================================================================#
library(rstan)
setwd("/Users/tm-pham/PhD/covid-19/nosocomialtransmission/stochasticepidemicmodel/stan")
sim_data <- readRDS("sim_data_only_patients_1.RDS")
parameters <- readRDS("parameters_only_patients_1.RDS")

# Data
eps <- 10^(-7)

pDshape <- 10                      # Shape for gamma distribution for dispersion parameter for transmission process
pDrate <- 10                       # Rate for gamma distribution for dispersion parameter for transmission process
pDshape_obs <- 10                  # Shape for gamma distribution for dispersion parameter for observation process
pDrate_obs <- 10                   # Rate for gamma distribution for dispersion parameter for observation process

f_mu <- unlist(parameters)         # Mean for normal distribution for transmission parameters
f_sigma <-  rep(0.5, length(f_mu)) # Sigma for normal distribution for transmission parameters

sim_data <- append(sim_data, list(eps=eps, 
                                  pDshape=pDshape, pDrate=pDrate, 
                                  pDshape_obs=pDshape_obs, pDrate_obs=pDrate_obs, 
                                  f_mu=f_mu, f_sigma=f_sigma))

fit <- stan(
  file = "nosocomialtransmission_only_patients.stan",        # Input model version here
  data = sim_data,                                        # named list of data
  chains = 1,                                             # number of Markov chains
  warmup = 100,                                           # number of warmup iterations per chain
  iter = 1000,                                            # total number of iterations per chain
  cores = 1,                                              # number of cores (use one per chain)
  refresh = 1,                                            # no of runs at which progress is shown
  control = list(max_treedepth = 15, adapt_delta=0.99)
)

save(fit, file="fit_simplified_only_patients_1.RDS")

p<-c("R_pU_p", 
     "R_p_p")

summary(fit,p)
stan_plot(fit,p)
stan_trace(fit,p)
check_divergences(fit)

