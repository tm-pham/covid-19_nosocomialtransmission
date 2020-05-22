# SEIR model
data {
  int<lower=0> N; // total number of individuals
  int<lower=0> T; // number of time points
  real<upper=1> alpha_1; // rate for latent period
  real<upper=1> alpha_2; // rate from pre-symptomatic to symptomatic
  real<upper=1> rep_delay; // rate for delay between developing symptoms and 
  real<upper=1> gamma_A; // recovery rate for asymptomatics
  real<upper=1> gamma_S; // recovery rate for symptomatics
  // Intitial values for patient compartments
  int<lower = 0, upper =N> S0_p[T]; // Number of initial susceptible patients
  int<lower = 0, upper =N> E0_p[T]; 
  int<lower = 0, upper =N> I0_pP[T];
  int<lower = 0, upper =N> I0_pS[T];
  int<lower = 0, upper =N> I0_pA[T];
  int<lower = 0, upper =N> R0_p[T];
    // Intitial values for HCW compartments
  int<lower = 0, upper =N> S0_hcw[T]; // Number of initial susceptible patients
  int<lower = 0, upper =N> E0_hcw[T]; 
  int<lower = 0, upper =N> I0_hcwP[T];
  int<lower = 0, upper =N> I0_hcwS[T];
  int<lower = 0, upper =N> I0_hcwA[T];
  int<lower = 0, upper =N> R0_hcw[T];
  
  
  int<lower = 0, upper =N> I_pCS[T]; // Number of symptomatic cohorted patients
  int<lower = 0, upper =N> I_hS[T];  // Number of symptomatic HCWs
  int<lower = 0, upper =N> I_hCS[T]; // Number of symptomatic cohorted HCWs (admitted as patients)
}

transformed data { 
}

parameters {
  real<lower=0> R0;
  real ker_shape; // shape parameter for transmission kernel 
  real ker_scale; // scale parameter for transmission kernel
  // Other paraemters
  real<upper=1> delta_P;

}


model{}

generated quantities {
  int N_p[T]; // Patients participating in the contact process
  int N_hcw[T]; // HCWs participating in the contact process
  int<lower = 0, upper =N> N_ncp[T]; // Number of non-covid patients
  int<lower = 0, upper =N> S_p[T]; // Number of initial susceptible patients
  int<lower = 0, upper =N> E_p[T]; 
  int<lower = 0, upper =N> I_pP[T];
  int<lower = 0, upper =N> I_pS[T];
  int<lower = 0, upper =N> I_pA[T];
  int<lower = 0, upper =N> R_p[T];
    // Intitial values for HCW compartments
  int<lower = 0, upper =N> S_hcw[T]; // Number of initial susceptible patients
  int<lower = 0, upper =N> E_hcw[T]; 
  int<lower = 0, upper =N> I_hcwP[T];
  int<lower = 0, upper =N> I_hcwS[T];
  int<lower = 0, upper =N> I_hcwA[T];
  int<lower = 0, upper =N> R_hcw[T];
  
  
  
  for(i in 1:(T-1)){
    N_p[i] = S_p[i] + E_p[i] + I_pA[i] + I_pP[i] + I_pS[i];
    N_hcw[i] = S_hcw[i] + E_hcw[i] + I_hcwP[i] + I_hcwA[i] + R_hcw[i];
      
      
      
  }
}
