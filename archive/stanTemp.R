# SEIR model
data {
  int<lower=0> N; // total number of individuals
  int<lower=0> T; // number of time points
  real<upper=1> alpha_1; // rate for latent period
  real<upper=1> alpha_2; // rate from pre-symptomatic to symptomatic
  real<upper=1> rep_delay; // rate for delay between developing symptoms and 
  real<upper=1> gamma_A; // recovery rate for asymptomatics
  real<upper=1> gamma_S; // recovery rate for symptomatics
  real<upper=1> P_A; // proportion of asymptomatics
  real discharge_IpP; // discharge rate of pre-symptomatic (but will return)
  real discharge_ncp; // discharge rate of non-covid patients (except pre-symptomatic)
  int max_gen;
  real gen_shape;
  real gen_scale;
  int r; // reduction factor of infectiousness of cohorted symptomatic patients
  
  // Intitial values for patient compartments
  int<lower = 0, upper =N> S0_p; // Number of initial susceptible patients
  int<lower = 0, upper =N> E0_p; 
  int<lower = 0, upper =N> I0_pP;
  int<lower = 0, upper =N> I0_pS;
  int<lower = 0, upper =N> I0_pA;
  int<lower = 0, upper =N> R0_p;
  // Intitial values for HCW compartments
  int<lower = 0, upper =N> S0_hcw; // Number of initial susceptible patients
  int<lower = 0, upper =N> E0_hcw; 
  int<lower = 0, upper =N> I0_hcwP;
  int<lower = 0, upper =N> I0_hcwS;
  int<lower = 0, upper =N> I0_hcwA;
  int<lower = 0, upper =N> R0_hcw;
  
  
  int<lower = 0, upper =N> I_pCS[T]; // Number of symptomatic cohorted patients
  int<lower = 0, upper =N> I_pDS[T]; // Number of readmitted patients (into I_pCS) who were discharged without symptoms
  int<lower = 0, upper =N> I_pNS[T]; // Number of newly admitted symptomatic patients
  int<lower = 0, upper =N> I_hS[T];  // Number of symptomatic HCWs (take from absence data)
  int<lower = 0, upper =N> I_hCS[T]; // Number of newly  symptomatic cohorted HCWs (admitted as patients)
  int<lower = 0, upper =N> N_NCP[T]; // Non-cohorted patients
}

transformed data { 
}

parameters {
  real<lower=0> R0;
  real ker_shape; // shape parameter for transmission kernel 
  real ker_scale; // scale parameter for transmission kernel
  // Other paraemters
  real<upper=1> dispersion;
  
}


model{}


functions {
  matrix[] create_data_template(S,E,I_A,I_P,I_S,R,EE=0,T){
    N = S+E+I_A+I_P+I_S+R+EE;
    vector[N] data_vec;
    if(S>0){
      for(i in 1:S){
        data_vec[i] = 1;
      }
    }
    if(E)
  }
  
  // returns the number of times the value x was found in vector y
  int num_val(int x, vector y){
    int n = 0;
    for(i in 1:rows(y)){
      if(y[i]==x){
        n = n + 1;
      }
    }
    return n;
  }
  
  // returns integer vactor of indexes of elements equal to x
  int[] find(int x, vector y){
    vector[num_val(x,y)] result;
    in n = 1;
    for(i in 1:rows(y)){
      if(y[i]==x){
        result[n] =i;
        n = n+1;
      }
    }
    return result;
  }
  
}

// 1 = S, 2 = E, 3 = I_A, 4 = I_P, 5 = I_S, 6 = I_SC, 7 = I_PD, 8 = R, 
// 9 = Discharged but will return, 10 = Discharged
generated quantities {
  int N_p[T]; // Patients participating in the contact process
  int N_hcw[T]; // HCWs participating in the contact process
  int<lower = 0, upper =N> N_ncp[T]; // Number of non-covid patients
  int<lower = 0, upper =N> dS_p[T]; // Number of initial susceptible patients
  int<lower = 0, upper =N> dE_p[T]; 
  int<lower = 0, upper =N> dI_pP[T];
  int<lower = 0, upper =N> dI_pPD[T];
  int<lower = 0, upper =N> dI_pS[T];
  int<lower = 0, upper =N> dI_pSC[T]; // Cohorted symptomatic patients 
  int<lower = 0, upper =N> dI_pA[T];
  int<lower = 0, upper =N> dR_p[T];
  
  int<lower = 0, upper =N> S_p[T]; // Number of initial susceptible patients
  int<lower = 0, upper =N> E_p[T]; 
  int<lower = 0, upper =N> I_pP[T];
  int<lower = 0, upper =N> I_pS[T];
  int<lower = 0, upper =N> I_pA[T];
  int<lower = 0, upper =N> R_p[T];
  // Intitial values for HCW compartments
  int<lower = 0, upper =N> dS_hcw[T]; // Number of initial susceptible patients
  int<lower = 0, upper =N> dE_hcw[T]; 
  int<lower = 0, upper =N> dI_hcwP[T];
  int<lower = 0, upper =N> dI_hcwS[T];
  int<lower = 0, upper =N> dI_hcwA[T];
  int<lower = 0, upper =N> dR_hcw[T];
  
  int<lower = 0, upper =N> S_hcw[T]; // Number of initial susceptible patients
  int<lower = 0, upper =N> E_hcw[T]; 
  int<lower = 0, upper =N> I_hcwP[T];
  int<lower = 0, upper =N> I_hcwS[T];
  int<lower = 0, upper =N> I_hcwA[T];
  int<lower = 0, upper =N> R_hcw[T];
  
  int pat_past_inf[max_gen];
  
  
  int num_E=0;
  int num_IA=0;
  int num_IP=0;
  int num_IS=0;
  
  real inf_ncp; // infectivity from non-cohorted patients
  real inf_cp;  // infectivity from cohorted patients
  
  matrix[T,N] pata_wdata; 
  matrix[E0_p,3] pat_data; // number of days since infection, state, flag for infection
  pat_data = append_col(rep_col_vector(1,E0_p),rep_col_vector(2,E0_p), rep_col_vector(1,E0_p));
  vector[T] inf_on_day =rep_row_vector(0,T);
  vector[T] inf_on_day_a = rep_row_vector(0,T);
  vector[T] inf_on_day_s = rep_row_vector(0,T);
  
  // Infectiousness function
  vector[max_gen] beta = weibull_lpdf(1:max_gen, gen_shape, gen_scale);
  
  
  vector[N] inf_time;
  //pat_past_inf = append_col(E0_p,rep_row_vector(0,max_gen));
  // Initialize:
    real prob_inf_pat;
  real prob_inf_hcw;
  
  // Initialize compartments
  // TBD
  
  for(i in 1:(T-1)){
    pat_past_inf = append_array(0,pat_past_inf);
    N_p[i] = S_p[i] + E_p[i] + I_pA[i] + I_pP[i] + I_pS[i];
    N_hcw[i] = S_hcw[i] + E_hcw[i] + I_hcwP[i] + I_hcwA[i] + R_hcw[i];
    
    // Force of infection
    // Days since infection for non-cohorted patients
    num_E = num_val(2, pat_data[,2]);
    // Days since infecttion for exposed patients
    vector[num_E] num_days_since_inf_NCP = pat_data[find(2, pat_data[,2]),1]; 
    // Days since infection for asymptomatic patients
    num_days_since_inf_NCP = append_array(num_days_since_inf_NCP, pat_data[find(3, pat_data[,2]),1]);
    // Days since infection for pre-symptomatic patients
    num_days_since_inf_NCP = append_array(num_days_since_inf_NCP, pat_data[find(4, pat_data[,2]),1]);
    // Days since infection for symptomatic patients
    num_days_since_inf_NCP = append_array(num_days_since_inf_NCP, pat_data[find(5, pat_data[,2]),1]);
    // Total number of non-cohorted patients
    total_num_ncp = num_E + num_IA + num_IP + num_IS;
    
    // Infectivity for non-cohorted patients
    inf_ncp = 0;
    for(d in 1:total_num_cp){
      inf_ncp += weibull_lpdf(num_days_since_inf_NCP[d], gen_shape, gen_scale);
    }
    
    // Days since infection for cohorted patients
    num_ISC = num_val(6, pat_data[,2]);
    vector[num_ISC] num_days_since_inf_CP = pat_data[find(6, pat_data[,2]),1]; 
    // Infectivity for non-cohorted patients 
    inf_cp = 0;
    for(d in 1:num_ISC){
      inf_cp += weibull_lpdf(num_days_since_inf_CP[d], gen_shape, gen_scale);
    }
    
    // Force of infection experienced by patients
    foi_pat = R_ncp*inf_ncp + R_cp*inf_cp;
    prob_inf_pat = 1-exp(-foi_pat);
    
    // ================================================== //
      // Patient model
    // ================================================== //
      // Number of non-covid patients (susceptibles + undetected infected)
    N_ncp[i] = S_p[i] + E_p[i] + I_pA[i] + I_pP[i] + I_pS[i];
    
    // Transitions
    // Discharge of susceptibles
    discharged_Sp = binomial_rng(S_p[i],prob_discharge_per_day);
    ind_S = find(1,pat_data[,2]); 
    // Change the state from S to Discharged
    pat_data[ind_S[1:discharged_Sp]] = 10;
    // Infection of susceptibles
    inf_on_day[i+1] = beta_binomial_rng(S_p[i]-discharged_Sp,prob_inf_pat,dispersion);
    new_infections = inf_on_day[i+1];
    S[i+1] = S[i]-inf_on_day[i+1];
    ind_S = find(1,pat_data[,2]); 
    
    pat_past_inf[i] = pat_past_inf + new_infections; 
    // Change the state from S to E for infected patients
    pat_data[ind_S[1:inf_on_day[i+1]],2] = 2; 
    
    // Discharge of exposed patients
    discharged_Ep = binomial_rng(E_p[i],prob_discharge_per_day);
    // Exposed to asymptomatic
    dI_pA[i+1] =  binomial_rng(E_p[i]-discharged_Ep, 1-exp(-P_A*alpha_1));
    // Exposed to pre-symptomatic
    dI_pP[i+1] = binomial_rng(E_p[i], 1-exp(-(1-P_A)*alpha_1));
    // Discharge of pre-symptomatic patients
    dI_pPD[i+1] = binomial_rng(I_pP[i], 1-exp(-discharge_IpP));
    // Transition from presymptomatic to symptomatic (not detected)
    dI_pS[i+1] = binomial_rng(I_pP[i], 1-exp(-alpha_2));
    // Transition from asymptomatic to recovered
    R_p[i+1] = binomial_rng(I_pA[i], 1-exp(-gamma_A)) + binomial_rng(I_pS[i], 1-exp(-gamma_S));
    
    // Transition from undetected symptomatics to detected symptomatic 
    E_p[i+1] = E_p[i] + dE_p[i+1] - dI_pA[i+1] - dI_pP[i+1];
    
    // Update pat_data
    ind_inf = find(1, pat_data[,3]);
    pat_data[ind_inf,1] +=1;
    
    
    // ================================================== //
      // HCW model
    // ================================================== //
      // Infection of susceptibles
    dE_hcw[i+1] = beta_binomial_rng(S_hcw[i],prob_inf_pat,dispersion);
    // Transition from exposed to asymptomatic
    dI_hcwA[i+1] = binomial_rng(E_hcw[i], 1-exp(-P_A*alpha_1));
    // Transition from exposed to pre-symptomatic
    dI_pP[i+1] = binomial_rng(E_hcw[i], 1-exp(-(1-P_A)*alpha_1));
    // Discharge of pre-symptomatic patients
    dI_pPD[i+1] = binomial_rng(I_pP[i], 1-exp(-discharge_IpP));
    // Transition from presymptomatic to symptomatic (not detected)
    dI_pS[i+1] = binomial_rng(I_pP[i], 1-exp(-alpha_2));
    // Transition from asymptomatic to recovered
    R_p[i+1] = binomial_rng(I_pA[i], 1-exp(-gamma_A)) + binomial_rng(I_pS[i], 1-exp(-gamma_S));
    
    // Transition from undetected symptomatics to detected symptomatic 
    E_p[i+1] = E_p[i] + dE_p[i+1] - dI_pA[i+1] - dI_pP[i+1];
    
    
    
    
    
    
    
  }
}
