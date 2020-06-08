// ========================================================================== //
// NOSOCOMIAL TRANSMISSION
// STAN code
// ========================================================================== //
data {
  int T;                              // Study period
  int S;                              // Maximum number of days to be infectious 
  int<lower=0> N_ncp[T];              // Number of non-cohorted patients at time t
  int<lower=0> I_p[T];                // Number of (isolated) infected patients in hospital at time t
  int obs_nosocomial[T];              // Observed number of nosocomial cases at time t
  
  int<lower=0> S_hcw[T];              // Number of susceptible HCWs at time t
  int hcw_isolation_period;           // Number of days HCWs are isolated

  real delta[T];                      // Probability to be discharged or die at time t

  real gen_shape;                     // shape parameter for generation time distribution (Feretti et al)
  real gen_scale;                     // scale parameter for generation time distribution (Feretti et al)
  
  real meanlog;                       // log mean for log normal distribution for probability of developing symptoms after infection
  real sdlog;                         // log sd for log normal distribution for probability of developing symptoms after infection

  real prob_asymptomatic;             // Proportion of infected individuals who are asymptomatic
  
  real p_p_obs;                       // Proportion of observed patient infections (CO-CIN)
  
  real pDshape;                       // Shape for gamma distribution for dispersion parameter for transmission process
  real pDrate;                        // Rate for gamma distribution for dispersion parameter for transmission process

  real pDshape_obs;                   // Shape for gamma distribution for dispersion parameter for observation process
  real pDrate_obs;                    // Rate for gamma distribution for dispersion parameter for observation process

  real fshape;                        // Shape for gamma distribution for transmission parameters
  real frate;                         // Rate for gamma distribution for transmission parameters
}

transformed data { 
  real inc_distr[T];                 // Probability of developing symptoms after s days of infection
  real gen_time[T];                  // Generation time distribution 
  simplex[3] p_multi[T-1];           // Array of simplex for multinomial distribution
  vector[3] temp;
  
  for(t in 1:T){
    inc_distr[t] = (1-prob_asymptomatic)*exp(lognormal_lpdf(t|meanlog, sdlog));
    gen_time[t] = weibull_lpdf(t| gen_shape, gen_scale);
  }
  
  for(s in 2:S){
    temp[1] = delta[s-1];
    temp[2] = (1-delta[s-1])*inc_distr[s];
    temp[3] = (1-delta[s-1])*(1-inc_distr[s]);
    p_multi[s-1]= temp;  
  }
}

parameters {
  real f_pU_hcw;                      // from unknown infected patient to HCW
  real f_pU_p;                        // from unkown infected patient to patient
  real f_p_hcw;                       // from known infected patient to HCW
  real f_p_p;                         // from known infected patient to patient
  real f_hcw_hcw;                     // from HCW to HCW
  real f_hcw_p;                       // from HCW to susceptible patient
  real pDis;                          // Dispersion parameter for gamma distribution of transmission process 
  real pDis_obs;                      // Dispersion parameter for gamma distribution of observation process 
}

model{

  real S_p[T];                       // Number of susceptible patients at time t
  real I_pU[S,T];                    // Number of unknown (unisolated) infected patients at time t who were infected s days ago
  real new_symptomatic_pat[S,T];     // Number of unknown infected HCWs who got infected s-1 days ago and develop symptoms at time t
  int pat_outcomes[T,3];             // Array for multinomial distribution for change of I_pU
  
  real I_hcwU[S,T];                  // Number of unknown infected HCWs at time t who got infected s-1 days ago
  real hcw_symp;                     // Number of newly symptomatic HCWs
  real hcw_recover;                  // Number of newly recovered HCWs
  real new_symptomatic_hcw[S,T];     // Number of unknown infected HCWs who got infected s-1 days ago and develop symptoms at time t
  real isolated_hcw[T];              // Number of isolated HCWs at time t
  real R_hcw[T];                     // Number of immune HCWs at time t
  
  real p_gamma_rate[T];              // Rate for gamma distribution for I_pU
  real p_gamma_shape[T];             // Shape for gamma distribution for I_pU
  
  real hcw_gamma_rate[T];            // Rate for gamma distribution for I_hcwU
  real hcw_gamma_shape[T];           // Shape for gamma distribution for I_hcwU
  
  real inf_p[T];                     // Infectivity of (unknown) infected patients
  real inf_hcw[T];                   // Infectivity of (unknown) infected HCWs
  
  real p_p[T];                       // Probability of infection for patients at time t
  real p_hcw[T];                     // Probability of infection for HCWs at time t

  // Priors for dispersion parameters
  pDis ~ gamma(pDshape,pDrate);                // Dispersion parameter for transmission process
  pDis_obs ~ gamma(pDshape_obs, pDrate_obs);   // Dispersion parameter for observation process
  
  // Priors for transmission parameters
  f_pU_p ~ gamma(fshape, frate);
  f_pU_hcw ~ gamma(fshape, frate);
  f_p_hcw ~ gamma(fshape, frate);
  f_p_p ~ gamma(fshape, frate);
  f_hcw_hcw ~ gamma(fshape, frate);
  f_hcw_p ~ gamma(fshape, frate);
  
  
  for(t in 2:T){
    inf_hcw[t] = 0;
    inf_p[t] = 0;
    for(s in 1:(t-1)){
      // Cumulative infectivity from HCWs
      inf_hcw[t] += gen_time[s]*I_hcwU[s,t-1];
      // Cumulative infectivity from unknown infected patients
      inf_p[t] += gen_time[s]*I_pU[s,t-1];
    }
    
    // Probability of infection for HCWs at time t
    p_hcw[t] = 1-exp(-f_hcw_hcw*inf_hcw[t] - f_pU_hcw*inf_p[t] - f_p_p*I_p[t]);
    // Probability of infection for patients at time t
    p_p[t] = 1-exp(-f_hcw_p*inf_hcw[t] - f_pU_p*inf_p[t] - f_p_hcw*I_p[t]);
    
    
    hcw_gamma_shape[t] = p_hcw[t]*hcw_gamma_rate[t]*S_hcw[t];
    hcw_gamma_rate[t] = (pDis + p_hcw[t]*(1-p_hcw[t]))/((1-p_hcw[t])*(pDis + S_hcw[t]*p_hcw[t]*(1-p_hcw[t]))); 
    
    p_gamma_shape[t] = p_p[t]*p_gamma_rate[t]*S_p[t];
    p_gamma_rate[t] = (pDis + p_p[t]*(1-p_p[t]))/((1-p_p[t])*(pDis + S_p[t]*p_p[t]*(1-p_p[t]))); 
    
    // Number of newly infected HCWs at time t (Continuous gamma approximation)
    I_hcwU[1,t] ~ gamma(hcw_gamma_shape[t],hcw_gamma_rate[t]/p_p_obs);
    
    // Number of newly infected patients at time t (Continuous gamma approximation)
    I_pU[1,t] ~ gamma(p_gamma_shape[t],p_gamma_rate[t]/p_p_obs);
    
    // Number of susceptible patients
    S_p[t] = N_ncp[t] - I_pU[1,t];
    
    for(s in 2:S){
      // HEALTH-CARE WORKERS
      // Number of unknown infected HCWs who got infected s-1 days ago and develop symptoms at time t
      real hcw_symp_gamma_rate = 1.0/(1.0-inc_distr[s]);
      real hcw_symp_gamma_shape = I_hcwU[s-1,t-1]*inc_distr[s]*hcw_symp_gamma_rate;
      hcw_symp ~ gamma(hcw_symp_gamma_shape,hcw_symp_gamma_rate);
      new_symptomatic_hcw[s,t] = hcw_symp;
      // Remaining unknown infected HCWs at time t who got infected s-1 days ago
      I_hcwU[s,t] = I_hcwU[s-1,t-1] - hcw_symp;
      
      // PATIENTS
      pat_outcomes[s,] ~ multinomial(p_multi[s-1]);
      // Number of unknown (unisolated) infected patients at time t who got infected s days ago
      I_pU[s,t] = pat_outcomes[s,3];
      // Number of unisolated infected patients who were infected s days ago and developed symptoms at time t
      new_symptomatic_pat[s,t] = pat_outcomes[s,2];
    }
    // Assume that symptomatic HCWs self-isolate and return to work after 7 days with immunity
    if(t > hcw_isolation_period){
      hcw_recover = sum(new_symptomatic_pat[1:(t-hcw_isolation_period),t-hcw_isolation_period]);
    }else{
      hcw_recover = 0;
    }
  // Number of known (isolated) infected patients at time t
  obs_nosocomial ~ neg_binomial(I_pU[1,t], pDis_obs);
  
  // Number of isolated HCWs
  // = Number of isolated HCWs the day before - those that recover and those that are still infected at time t (got infected 1,2,...,t days ago)
  // Assume HCWs are isolated immediately
  isolated_hcw[t] = isolated_hcw[t-1] - hcw_recover + sum(new_symptomatic_hcw[1:t,t]);
  // Number of recovered HCWs
  R_hcw[t] = R_hcw[t-1] + hcw_recover;
    
  }
}
