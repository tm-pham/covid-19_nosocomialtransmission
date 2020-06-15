# =============================================================================#
# NOSOCOMIAL TRANSMISSION
# Run file for STAN code
# =============================================================================#
library(rstan)
setwd("/Users/tm-pham/PhD/covid-19/nosocomialtransmission/stochasticepidemicmodel/stan")
sim_data <- readRDS("sim_data.RDS")

# Data
pDis <- 10
pDis_obs <- 10

pDshape <- 10                      # Shape for gamma distribution for dispersion parameter for transmission process
pDrate <- 10                       # Rate for gamma distribution for dispersion parameter for transmission process
pDshape_obs <- 10                  # Shape for gamma distribution for dispersion parameter for observation process
pDrate_obs <- 10                   # Rate for gamma distribution for dispersion parameter for observation process

f_mu <- c(0.001, 0.001, 0.0001, 0.0001, 0.0005, 0.001)  # Mean for normal distribution for transmission parameters
f_sigma <-  rep(0.001, length(f_mu))   # Sigma for normal distribution for transmission parameters

sim_data <- append(sim_data, list(pDis=pDis, pDis_obs=pDis_obs, pDshape=pDshape, pDrate=pDrate, pDshape_obs=pDshape_obs, pDrate_obs=pDrate_obs, 
                                  f_mu=f_mu, f_sigma=f_sigma))

fit <- stan(
  file = "nosocomialtransmission.stan",  # Input model version here 
  data = sim_data,                       # named list of data defined in metareg_define_data.R
  chains = 1,                            # number of Markov chains
  warmup = 100,                         # number of warmup iterations per chain
  iter = 1000,                         # total number of iterations per chain
  cores = 1,                             # number of cores (use one per chain)
  refresh = 1000,                        # no of runs at which progress is shown
  control = list(max_treedepth = 15, adapt_delta=0.99)
)

