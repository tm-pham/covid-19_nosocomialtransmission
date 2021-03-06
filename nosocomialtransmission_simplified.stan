// ========================================================================== //
  // NOSOCOMIAL TRANSMISSION
// STAN code
// ========================================================================== //
  data {
    int T;                              // Study period
    int S;                              // Not used for now. Maximum number of days to be infectious 
    int <lower=0> N_ncp[T];              // Number of non-cohorted patients at time t
    int <lower=0> I_p[T];                // Number of (isolated) infected patients in hospital at time t
    real <lower=0> i_pU0;                // Initial number of unknown (unisolated) infected patients at time 1 (start of time period)
    real <lower=0> sum_symp_pat0;        // Initial number of patients that develop symptoms at time 1 (start of time period)
    int <lower=0> obs_nosocomial[T];     // Observed number of nosocomial cases at time t
    
    int <lower=0> N_hcw[T];              // Total number of HCWs at time t
    int <lower=0> hcw_isolation_period;  // Number of days HCWs are isolated
    real <lower=0> i_hcwU0;              // Initial number of unknown (unisolated) HCWs at time 1
    real <lower=0> isolated_hcw0;        // Initial number of isolated HCWs at time 1
    real sum_symp_hcw0;                 // Initial number of patients that develop symptoms at time 1 (start of time period)
    real r_hcw0;                        // Initial number of recovered HCWs at time 1
    int <lower=0> obs_hcw_infections[T]; // Observed number of infections in hcws at time t
    
    real <lower=0> delta[T];             // Probability to be discharged or die at time t
    
    real gen_shape;                     // Shape parameter for generation time distribution (Feretti et al)
    real gen_scale;                     // Scale parameter for generation time distribution (Feretti et al)
    
    real meanlog;                       // Log mean for log normal distribution for probability of developing symptoms after infection
    real sdlog;                         // Log sd for log normal distribution for probability of developing symptoms after infection
    
    real <lower=0> prob_asymptomatic;    // Proportion of infected individuals who are asymptomatic
    
    real <lower=0> p_p_obs;              // Proportion of observed patient infections (CO-CIN)
    real <lower=0> p_hcw_obs;            // Proportion of observed patient infections (CO-CIN)
    
    real eps;                           // Minimum value for S_p, S_hcw, ...
    
    real pDis;                          // Dispersion parameter for infection process
    real pDis_obs;                      // Dispersion parameter for observation process
    
    real pDshape;                       // Shape for gamma distribution for dispersion parameter for transmission process
    real pDrate;                        // Rate for gamma distribution for dispersion parameter for transmission process
    
    real pDshape_obs;                   // Shape for gamma distribution for dispersion parameter for observation process
    real pDrate_obs;                    // Rate for gamma distribution for dispersion parameter for observation process
    
    real f_mu[6];                       // Mean values for normal distribution for transmission parameters
    real f_sigma[6];                    // Sd values for normal distribution for transmission parameters
  }

transformed data { 
  real <lower=0> inc_distr[T];        // Probability of developing symptoms after s days of infection
  real gen_time[T];                  // Generation time distribution 
  simplex[3] p_multi[T-1];           // Array of simplex for multinomial distribution
  vector[3] temp;
  real <lower=0> rev_inc_distr[T];
  
  for(t in 1:T){
    inc_distr[t] = (1.0-prob_asymptomatic)*exp(lognormal_lpdf(t|meanlog, sdlog));
    gen_time[t] = exp(weibull_lpdf(t| gen_shape, gen_scale));
  }
  
  for(t in 1:T){
    rev_inc_distr[t] = inc_distr[T-t+1];
  }
  
  for(s in 2:T){
    temp[1] = delta[s-1];
    temp[2] = (1-delta[s-1])*inc_distr[s];
    temp[3] = (1-delta[s-1])*(1-inc_distr[s]);
    p_multi[s-1]= temp;  
  }
}

parameters {
  // For uniform distribution
  // real <lower=0.00001, upper=f_mu[1]+0.00015>f_pU_hcw;            // from unknown infected patient to susceptible HCW
  // real <lower=0.00001, upper=f_mu[2]+0.00015>f_pU_p;                // from unkown infected patient to susceptible patient
  // real <lower=0.00001, upper=f_mu[3]+0.00015>f_p_hcw;              // from known infected patient to susceptible HCW
  // real <lower=0.00001, upper=f_mu[4]+0.00015>f_p_p;                  // from known infected patient to susceptible patient
  // real <lower=0.00001, upper=f_mu[5]+0.00015>f_hcw_hcw;          // from unknown infected HCW to susceptible HCW
  // real <lower=0.00001, upper=f_mu[6]+0.00015>f_hcw_p;              // from unknown infected HCW to susceptible patient
  
  // For normal distribution
  real <lower=0>f_pU_hcw;                      // from unknown infected patient to HCW
  real <lower=0>f_pU_p;                        // from unkown infected patient to patient
  real <lower=0>f_p_hcw;                       // from known infected patient to HCW
  real <lower=0>f_p_p;                         // from known infected patient to patient
  real <lower=0>f_hcw_hcw;                     // from HCW to HCW
  real <lower=0>f_hcw_p;                       // from HCW to susceptible patient
  
  //real <lower=0>pDis;                        // Dispersion parameter for gamma distribution of transmission process 
  //real <lower=0>pDis_obs;                    // Dispersion parameter for gamma distribution of observation process 
  
  real <lower=0>I_hcwU[T];                    // Number of newly infected HCWs at time t (Continuous gamma approximation)
  real <lower=0>I_pU[T];                      // Number of newly infected patients at time t (Continuous gamma approximation)
  
  real <lower=0>sum_symp_hcw[T];               // Number of HCWs that develop symptoms at time t
  real <lower=0>sum_symp_pat[T];               // Number of patients that develop symptoms at time t
  real <lower=0>discharge_dead_pat[T];
}

model{
  real S_p[T];                       // Number of susceptible patients at time t
  real Prev_pU[T];                   // Prevalence of unknown (unisolated) infected patients at time t
  real temp_prev_pU;                 // Temporary variable for prevalence of unknown infected patients
  
  real S_hcw[T];                     // Number of susceptible HCWs
  real hcw_recover;                  // Number of newly recovered HCWs
  real isolated_hcw[T];              // Number of isolated HCWs at time t
  real R_hcw[T];                     // Number of immune HCWs at time t

  real Prev_hcwU[T];                 // Prevelance of unknown (unisolated) infected patients at time t
  
  real p_gamma_rate[T-1];            // Rate for gamma distribution for I_pU
  real p_gamma_shape[T-1];           // Shape for gamma distribution for I_pU
  
  real p_d_gamma_rate;               // Rate for gamma distribution for I_pU that are discharged or die
  real p_d_gamma_shape;              // Shape for gamma distribution for I_pU that are discharged or die
  
  real p_symp_gamma_rate;            // Rate for gamma distribution for I_pU[t] that become symptomatic
  real p_symp_gamma_shape;           // Shape for gamma distribution for I_pU[t] that become symptomatic
  
  real hcw_gamma_rate[T-1];          // Rate for gamma distribution for I_hcwU[t]
  real hcw_gamma_shape[T-1];         // Shape for gamma distribution for I_hcwU[t]
  
  real hcw_symp_gamma_rate;          // Rate for gamma distribution for I_hcwU[s,t]
  real hcw_symp_gamma_shape;         // Shape for gamma distribution for I_hcwU[s,t]
  
  real inf_p[T-1];                   // Infectivity of (unknown) infected patients
  real inf_hcw[T-1];                 // Infectivity of (unknown) infected HCWs
  
  real p_p[T-1];                     // Probability of infection for patients at time t
  real p_hcw[T-1];                   // Probability of infection for HCWs at time t
  
  real prob_symp;                    // Probability for patient to stay on the ward and develop symptoms
  
  int ind;                           // Temporary variable for indices
  
  real prob_inc;
  
  // Priors for dispersion parameters
  // pDis ~ gamma(pDshape,pDrate);               // Dispersion parameter for transmission process
  // pDis_obs ~ gamma(pDshape_obs, pDrate_obs);  // Dispersion parameter for observation process
  
  // Priors for transmission parameters
  // Half-normal distribution
  f_pU_hcw ~ normal(f_mu[1], f_sigma[1]);
  f_pU_p ~ normal(f_mu[2], f_sigma[2]);
  f_p_hcw ~ normal(f_mu[3], f_sigma[3]);
  f_p_p ~ normal(f_mu[4], f_sigma[4]);
  f_hcw_hcw ~ normal(f_mu[5], f_sigma[5]);
  f_hcw_p ~ normal(f_mu[6], f_sigma[6]);
  
  // Uniform distribution 
  //f_pU_hcw ~ uniform(0.00001, f_mu[1]+0.00015);
  //f_pU_p ~ uniform(0.00001, f_mu[2] +0.00015);
  //f_p_hcw ~ uniform(0.00001, f_mu[3] +0.00015);
  //f_p_p ~ uniform(0.00001, f_mu[4] +0.00015);
  //f_hcw_hcw ~ uniform(0.00001, f_mu[5] +0.00015);
  //f_hcw_p ~ uniform(0.00001, f_mu[6] +0.00015);
  
  // ================ //
    // Initial values
  // ================ //
  // Initial number of susceptible patients
  S_p[1] = N_ncp[1] - i_pU0;
  // Intitial values for HCWs
  isolated_hcw[1] = isolated_hcw0;
  R_hcw[1] = r_hcw0;
  S_hcw[1] = fmax(N_hcw[1] - R_hcw[1] - i_hcwU0 - isolated_hcw[1], eps);
  
  Prev_hcwU[1] = i_hcwU0; 
  Prev_pU[1] = i_pU0; 
  // Start loop
  for(t in 2:T){
    // Cumulative infectivity 
    if(t==2){
      inf_hcw[t-1] = gen_time[t-1]*i_hcwU0;
      inf_p[t-1] = gen_time[t-1]*i_pU0;
    }else{
      inf_hcw[t-1] = dot_product(gen_time[1:(t-1)], I_hcwU[1:(t-1)]);
      inf_p[t-1] = dot_product(gen_time[1:(t-1)], I_pU[1:(t-1)]);
    }

    
    // =================================== //
      // Probability of infection at time t
    // ================================== //
      // HCWs
    p_hcw[t-1] = 1-exp(-f_hcw_hcw*inf_hcw[t-1] - f_pU_hcw*inf_p[t-1] - f_p_hcw*I_p[t-1]); // Should I_p be divided by N_ncp?
      // print("inf_hcw[", t-1, "]=", inf_hcw[t-1]);
    // print("inf_p[", t-1, "]=", inf_p[t-1]);
    // print("I_p[", t-1, "]=", I_p[t-1]);
    // print("N_ncp[", t-1, "]=", N_ncp[t-1]);
    // Patients
    p_p[t-1] = 1-exp(-f_hcw_p*inf_hcw[t-1] - f_pU_p*inf_p[t-1] - f_p_p*I_p[t-1]); // Should I_p be divided by N_ncp?
      
      // Shape and rate parameters for gamma distribution for new infections
    hcw_gamma_rate[t-1] = (pDis + p_hcw[t-1]*(1-p_hcw[t-1]))/((1-p_hcw[t-1])*(pDis + S_hcw[t-1]*p_hcw[t-1]*(1-p_hcw[t-1]))); 
    hcw_gamma_shape[t-1] = fmax(p_hcw[t-1]*hcw_gamma_rate[t-1]*S_hcw[t-1],eps);
    
    p_gamma_rate[t-1] = (pDis + p_p[t-1]*(1-p_p[t-1]))/((1-p_p[t-1])*(pDis + S_p[t-1]*p_p[t-1]*(1-p_p[t-1]))); 
    p_gamma_shape[t-1] = fmax(p_p[t-1]*p_gamma_rate[t-1]*S_p[t-1], eps);
    
    // ========================================= //
      // Number of newly infected HCWs at time t 
    // (Continuous gamma approximation)
    // ========================================= //
      // HCWs
    I_hcwU[t] ~ gamma(hcw_gamma_shape[t-1],hcw_gamma_rate[t-1]);
    // print("I_hcwU[", t, "]=", I_hcwU[1,t]);
    // print("S_hcw[", t-1, "]=", S_hcw[t-1]);
    // print("p_hcw[", t-1, "]=", p_hcw[t-1]);
    // print("hcw_gamma_shape[", t-1, "]=", hcw_gamma_shape[t-1]);
    // Patients
    I_pU[t] ~ gamma(p_gamma_shape[t-1],p_gamma_rate[t-1]/p_p_obs);
    // print("I_pU[", t, "]=", I_pU[1,t]);
    // print("S_p[", t-1, "]=", S_p[t-1]);
    // print("p_p[", t-1, "]=", p_p[t-1]);
    // print("p_gamma_shape[", t-1, "]=", p_gamma_shape[t-1]);
    
    
    // Remaining number of susceptible patients
    S_p[t] = fmax(N_ncp[t] - I_pU[t], eps);
    
    // ========================================================================== //
    // Number of individuals that develop symptoms at time t  
    // ========================================================================= //
    // HCWs
    ind = max(t-S, 1);
    prob_inc = sum(inc_distr[ind:(t-1)]);
    //prob_inc = dot_product(I_hcwU[1:(t-1)], rev_inc_distr[1:(t-1)]); 
    hcw_symp_gamma_rate = 1.0/(1.0-prob_inc);
    hcw_symp_gamma_shape = fmax(Prev_hcwU[t-1]*prob_inc*hcw_symp_gamma_rate, eps);
    sum_symp_hcw[t] ~ gamma(hcw_symp_gamma_shape,hcw_symp_gamma_rate);
    Prev_hcwU[t] = fmax(Prev_hcwU[t-1] - sum_symp_hcw[t] + I_hcwU[t], eps);
    
    // Patients
    // Patients that get discharged or die at time t
    p_d_gamma_rate = (pDis + delta[t]*(1-delta[t]))/((1-delta[t])*(pDis + Prev_pU[t-1]*delta[t]*(1-delta[t]))); 
    p_d_gamma_shape = fmax(Prev_pU[t-1]*delta[t]*p_d_gamma_rate,eps);
    discharge_dead_pat[t] ~ gamma(p_d_gamma_shape,p_d_gamma_rate);
    // Unknown infected patients that develop symptoms
    temp_prev_pU = fmax(Prev_pU[t-1] - discharge_dead_pat[t], eps);
    prob_symp = (1-delta[t])*prob_inc; 
    //prob_symp = (1-delta[t])*dot_product(I_pU[1:(t-1)], rev_inc_distr[1:(t-1)]);
    p_symp_gamma_rate = (pDis + prob_symp*(1-prob_symp))/((1-prob_symp)*(pDis + Prev_pU[t-1]*prob_symp*(1-prob_symp))); 
    p_symp_gamma_shape = fmax(temp_prev_pU*prob_symp*p_symp_gamma_rate, eps);
    sum_symp_pat[t] ~ gamma(p_symp_gamma_shape,p_symp_gamma_rate);
    // Remaining number of unknonwn unisolated patient at time t (prevalence)
    Prev_pU[t] = fmax(temp_prev_pU - sum_symp_pat[t] + I_pU[t], eps);
      
    
    // ========================================================= //
    // Number of known (isolated) infected patients at time t
    // ========================================================= //
    obs_nosocomial[t] ~ neg_binomial(p_p_obs*fmax(sum_symp_pat[t],eps), pDis_obs);
    obs_hcw_infections[t] ~ neg_binomial(p_hcw_obs*fmax(sum_symp_hcw[t],eps), pDis_obs);
    
    // ========================================================================================= //
    // Number of isolated HCWs
    // Assume that symptomatic HCWs self-isolate and return to work after 7 days with immunity
    // Assume HCWs are isolated immediately
    // ========================================================================================= //
    
    if(t > hcw_isolation_period){
      hcw_recover = sum_symp_hcw[t-hcw_isolation_period];
    }else{
      hcw_recover = 0;
    }
    // Number of isolated HCWs
    // = Number of isolated HCWs the day before - those that recover and those that are still infected at time t (got infected 1,2,...,t days ago)
    isolated_hcw[t] = fmax(isolated_hcw[t-1] - hcw_recover + sum(sum_symp_hcw[1:t]), eps);

    // Number of recovered HCWs
    R_hcw[t] = R_hcw[t-1] + hcw_recover;
    // Number of remaining susceptible HCWs
    S_hcw[t] = fmax(N_hcw[t] - R_hcw[t] - I_hcwU[t] - isolated_hcw[t], eps);
  }
}
