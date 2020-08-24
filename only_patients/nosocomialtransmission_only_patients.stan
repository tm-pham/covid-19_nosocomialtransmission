// ========================================================================== //
// NOSOCOMIAL TRANSMISSION
// STAN code 
// Highly simplified version - v5
// This version corresponds to simulation_only_patients.R
// We assume 
//- only one population (namely the patient population)
// - constant infecitous period
// - constant infectiousness over infectious period
// - prob_symptom_onset = fixed daily probability that asymptomatic infections develop symptoms

// The discrete beta binomial distribution is approximated by a 
// moment matched scaled beta distribution (instead of a gamma distribution
// used in Li et al, 2018)
// Reason is that when we used a gamma distribution, we had to use fmax
// because 
// 21st August 2020
// ========================================================================== //
  
  
functions {
    real scaledbeta_rng(real n, real alpha, real beta){
      real x;
      if (n <=0) reject("n must be a positive number; n=", n );
      if(alpha <= 0) reject("alpha must be >0; alpha=", alpha);
      if(beta <= 0) reject("beta must be >0; beta=", beta);  
      x = beta_rng(alpha, beta)* n;
      return x;
    }
    
    real scaledbeta_lpdf(real y, real n, real alpha, real beta){
      real x;
      real y2;
      print("y is ", y);
      print("n is ", n);
      print("alpha is ", alpha);
      print("beta is ", beta);
      
      
      if(y <= 0 || y>n ) reject("y must be between 0 and n; y=", y); 
      if (n <= 0) reject("n must be  positive ; n=", n );
      if(alpha <= 0) reject("alpha must be >0; alpha=", alpha);
      if(beta <= 0) reject("beta must be >0; beta=", beta); 
      y2 = y/n;
      x =  beta_lpdf(y2| alpha, beta);
      return x;
    }  
    
    real[] getscaledbetaparams(real m, real v, real n){
      // Arguments required: mean of distribution (m), variance (v) and max (n)
      // Functions returns alpha and beta of scaled beta distribution with scaling factor (max) n
      // that has the given mean and var (used for choosing a moment matched scaled beta distribution)
      real alpha;
      real beta;
      real x[2];
      alpha =  ((m^2)*n-m^3 -m*v)/(n*v);
      beta=alpha*((n/m)-1);
      x[1]=alpha;
      x[2]=beta;
      return x;
    }
  }  

data {
  int T;                               // Study period
  int <lower=0> N_ncp[T];              // Number of non-cohorted patients at time t
  int <lower=0> I_p[T];             // Number of (isolated) infected patients in hospital at time t
  real <lower=0> i_pU0;                // Initial number of unknown (unisolated) infected patients at time 1 (start of time period)
  real <lower=0> sum_symp_pat0;        // Initial number of patients that develop symptoms at time 1 (start of time period)
  int <lower=0> obs_nosocomial[T];     // Observed number of nosocomial cases at time t

  real <lower=0> delta;                // Probability to be discharged or die at time t
  
  int <lower=0> infectious_period;     // Infectious period (assumed to be constant for now)
  real <lower=0> prob_symptom_onset;   // Presymptomatic cases have a fixed daily probablity of developing symptoms
  
  real <lower=0> p_p_obs;              // Proportion of observed patient infections (CO-CIN)
  
  real eps;                            // Minimum value to avoid zero values
  
  real pDshape;                        // Shape for gamma distribution for dispersion parameter for transmission process
  real pDrate;                         // Rate for gamma distribution for dispersion parameter for transmission process
  
  real pDshape_obs;                    // Shape for gamma distribution for dispersion parameter for observation process
  real pDrate_obs;                     // Rate for gamma distribution for dispersion parameter for observation process

  real f_mu[2];                        // Mean for normal distribution for transmission reproduction numbers
  real f_sigma[2];                     // Sigma for normal distribution for transmission reproduction numbers
}

parameters {
  // Components of the next generation matrix (expected number of secondary infections from x to y)
  real <lower=0>R_pU_p;                        // from unkown infected patient to patient
  real <lower=0>R_p_p;                         // from known infected patient to patient

  real <lower=0, upper=100>I_pU[T];           // Number of newly infected patients at time t (Continuous gamma approximation)
  
  real <lower=0> obsMean_pat[T];              // Gamma distributed mean which when mixed with Poisson gives NB observation model for patients
  
  real <lower=0, upper=100>sum_symp_pat[T];   // Number of patients that develop symptoms at time t
  real <lower=0>discharged_dead_pat[T];
  
  real pDis;                                  // Dispersion parameter for alpha and beta
  real pDis_obs;                              // Dispersion parameter for observed number of infected patients
}

model{
  real S_p[T];                       // Number of susceptible patients at time t
  real Prev_pU[T];                   // Prevalence of unknown (unisolated) infected patients at time t
  real temp_prev_pU;                 // Temporary variable for prevalence of unknown infected patients
  
  real p_gamma_rate[T-1];            // Rate for gamma distribution for I_pU
  real p_gamma_shape[T-1];           // Shape for gamma distribution for I_pU
  
  real p_d_gamma_rate;               // Rate for gamma distribution for I_pU that are discharged or die
  real p_d_gamma_shape;              // Shape for gamma distribution for I_pU that are discharged or die
  
  real beta_params[2];               // Parameters for scaled beta distribution used as continuous moment matched version of discrete values
  
  real p_symp_gamma_rate;            // Rate for gamma distribution for I_pU[t] that become symptomatic
  real p_symp_gamma_shape;           // Shape for gamma distribution for I_pU[t] that become symptomatic
  
  real inf_p[T-1];                   // Infectivity of (unknown) infected patients

  real p_p[T-1];                     // Probability of infection for patients at time t

  real prob_symp;                    // Probability for patient to stay on the ward and develop symptoms
  
  int ind;                           // Temporary variable for indices
  
  real prob_inc;
  
  real alpha;                       //alpha for beta - binom distribution
  real beta;                        //beta for beta binom distribution
  real mean_beta_binom; 
  real var_beta_binom;
  
  // Priors for dispersion parameters
  pDis ~ uniform(0, 2);
  pDis_obs ~ normal(1000,20);       //  might to experiment with this prior - data might not be very informative about it 
  
  
  // Priors for transmission parameters
  R_pU_p ~ normal(f_mu[1], f_sigma[1]);
  R_p_p ~ normal(f_mu[2], f_sigma[2]);
  
  
  // ================ //
  // Initial values
  // ================ //
  // Initial number of susceptible patients
  S_p[1] = N_ncp[1] - i_pU0;
  Prev_pU[1] = i_pU0; 
  I_pU[1]  ~ normal(0,20);
  
  for(t in 2:T){
    // =================================== //
    // Probability of infection at time t
    // ================================== //
    // Patients
    p_p[t-1] = 1-exp(- R_pU_p*Prev_pU[t-1]/(N_ncp[t-1]*infectious_period) - R_p_p*I_p[t-1]/(N_ncp[t-1]*infectious_period));

        
    // ======================================================== //
    // Number of newly infected patients at time t 
    // (Continuous scaled beta distribution used as continuous 
    // moment matched approximation of beta binomial model)
    // ======================================================== //
    alpha = pDis/(1- p_p[t-1]);
    beta =  pDis/p_p[t-1];
    mean_beta_binom = (S_p[t-1]*alpha)/(alpha + beta);
    var_beta_binom = S_p[t-1] * alpha * beta * (alpha + beta + S_p[t-1])/((alpha + beta)^2 + (alpha + beta + 1));
    beta_params = getscaledbetaparams(mean_beta_binom, var_beta_binom ,S_p[t-1]);
        
    I_pU[t] ~ scaledbeta( S_p[t-1], beta_params[1], beta_params[2]);
        
    // Remaining number of susceptible patients
    S_p[t] = N_ncp[t] - I_pU[t];
        
    // ========================================================================== //
    // Number of individuals that develop symptoms at time t  
    // ========================================================================= //
    // Patients
    beta_params = getscaledbetaparams(delta*Prev_pU[t-1], delta*(1-delta)*Prev_pU[t-1], Prev_pU[t-1]);
    discharged_dead_pat[t] ~ scaledbeta(Prev_pU[t-1], beta_params[1], beta_params[2]);
        
    // Unknown infected patients that develop symptoms - now modelled without extra binomial variation)
    temp_prev_pU = Prev_pU[t-1] - discharged_dead_pat[t];
    prob_symp = (1-delta)*prob_symptom_onset; 

    beta_params = getscaledbetaparams(prob_symp*temp_prev_pU,prob_symp*(1-prob_symp)*temp_prev_pU ,temp_prev_pU);
    sum_symp_pat[t] ~ scaledbeta( temp_prev_pU, beta_params[1], beta_params[2]);

    // Remaining number of unknonwn unisolated patient at time t (prevalence)
    Prev_pU[t] = temp_prev_pU - sum_symp_pat[t] + I_pU[t];

    // ========================================================= //
   // Number of known (isolated) infected patients at time t
   // ========================================================= //
    obsMean_pat[t] ~ gamma(pDis_obs,pDis_obs/(sum_symp_pat[t]+eps));
    obs_nosocomial[t] ~ poisson(obsMean_pat[t]); 
  }
}

generated quantities {  //testing
  real<lower=0, upper=20> y_new[20];
  real testtest[2];
  for (i in 1:20){
    y_new[i] = scaledbeta_rng(20.0, 20.0, 20.0);
  }
  testtest=getscaledbetaparams(5,4,10);//returns alpha and beta so that n*rbeta(alpha,beta) has mean 5 and var 4
}
