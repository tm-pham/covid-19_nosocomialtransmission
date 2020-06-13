# =============================================================================#
# NOSOCOMIAL TRANSMISSION
# Run file for STAN code
# =============================================================================#
library(rstan)
setwd("/Users/tm-pham/PhD/covid-19/nosocomialtransmission/stochasticepidemicmodel/stan")
sim_data <- readRDS("sim_data.RDS")

# Data
pDshape <- 0.1                      # Shape for gamma distribution for dispersion parameter for transmission process
pDrate <- 0.5                       # Rate for gamma distribution for dispersion parameter for transmission process

pDshape_obs <- 0.1                  # Shape for gamma distribution for dispersion parameter for observation process
pDrate_obs <- 0.5                   # Rate for gamma distribution for dispersion parameter for observation process

fshape <- 1                         # Shape for gamma distribution for transmission parameters
frate <-  0.5                       # Rate for gamma distribution for transmission parameters

sim_data <- append(sim_data, list(pDshape=pDshape, pDrate=pDrate, pDshape_obs=pDshape_obs, pDrate_obs=pDrate_obs, fshape=fshape, frate=frate)) 

fit <- stan(
  file = "nosocomialtransmission.stan",  # Input model version here 
  data = sim_data,                       # named list of data defined in metareg_define_data.R
  chains = 2,                            # number of Markov chains
  warmup = 1000,                         # number of warmup iterations per chain
  iter = 100000,                         # total number of iterations per chain
  cores = 1,                             # number of cores (use one per chain)
  refresh = 1000,                        # no of runs at which progress is shown
  control = list(max_treedepth = 15, adapt_delta=0.99)
)

