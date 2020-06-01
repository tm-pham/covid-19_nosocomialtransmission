# SEIR model
data {
  int N;                              // Total number of individuals in hospital
  int T;                              // Study period

  int I_pU0;                          // Initial number of unknown infected patients
  int I_hcw0;                         // Initial number of unknown infected HCWs
  
  int<lower=0, upper=N> S_p[T];       // Number of susceptible patients at time t
  int<lower=0, upper=N> S_hcw[T];     // Number of susceptible HCWs at time t

  real delta[T];                      // Proportion of nosocomial infections, infected s days ago who are discharged, isolated or died
  real gamma[T];                      // Proportion of nosocomial HCW infections, infected s days ago who recover (and thus are immune)
  
  real gen_shape;                     // shape parameter for generation time distribution (Feretti et al)
  real gen_scale;                     // scale parameter for generation time distribution (Feretti et al)
  
  //real c_hcw_p;                     // Contact rate of HCWs with patients
  //real c_hcw_hcw;                   // Contact rate of HCWs with HCWs
  //real c_p_p;                       // Contact rate of patients with patients
  //real c_p_hcw;                     // Contact rate of patients with HCWs
  
  // Transmission parameters (assume to be known in simulation)
  real f_hcw_p;                       // known infected patient to HCW
  real f_hcw_pp;                      // unknown infected patient to HCW
  real f_hcw_hcw;                     // HCW to HCW
  real f_p_p;                         // Known infected patient to patient
  real f_p_pp;                        // Unkown infected patient to patient
  real f_p_hcw;                       // HCW to susceptible patient
}

transformed data { 
}

parameters {}


model{}

functions {
   // Returns the largest integer less than or equal to value
   // Return value is an integer (the built-in "floor" function by stan returns a real value)
   int floorInt(real value){
     int res=0; 
     // Artificially convert integer to real in order to compare (not sure if necessary)
     while(floor(res)<value){
       res +=1;
     }
     return(res) 
   }
}


generated quantities {
  // Usually observed but here for simulation?
  real II_hcw[T,T];    // (real) Number of infected HCWs at time t who were infected (s-1) days ago and are still working
  int I_hcw[T,T];      // (int) Number of infected HCWs at time t who were infected (s-1) days ago and are still working
  int R_hcw[T];        // Number of immune HCWs at time t

  real ww[T];         // (real) Number of (immune) HCWs returning to work at time t
  int w[T];           // (int) Number of (immune) HCWs returning to work at time t

  int II_pU[T,T];     // (real) Number of unknown (unisolated) infected patients at time t who were infected s days ago
  int I_pU[T,T];      // (int) Number of unknown (unisolated) infected patients at time t who were infected s days ago
  
  // Usually observed but here for simulation
  int I_p[T];         // Number of (isolated) infected patients in hospital at time t
  int N_ncp[T];       // Number of non-cohorted patients at time t

  real p_hcw[T];      // Probability of infection for HCWs at time t
  real p_p[T];        // Probability of infection for patients at time t
  real inf_hcw[T];    // Infectivity from HCWs (densitiy dependent)
  real inf_p[T];      // Infectivity from patients (densitiy dependent)
  
  // Generation time distribution used to determine the infectivity dependend on the time of infection
  real gen_time[T];
  for(t in 1:T){
    gen_time[t] = exp(weibull_lpdf(t | gen_shape, gen_scale)); 
  }

  I_pU[1,1] = I_pu0;   // Initial number of unknown infected unisolated patients
  I_hcw[1,1] = I_hcw0; // Intial number of infected HCWs 
  
  for(t in 2:T){
    for(s in 1:(t-1)){
      // cumualtive infectivity from HCWs
      inf_hcw[t] += gen_time[s]*I_hcw[s,t-1];
      // cumualtive infectivity from unknown infected patients
      inf_p[t] += gen_time[s]*I_pU[s,t-1];
    }
    // Probability of infection for HCWs at time t
    p_hcw[t] = f_hcw_hcw*inf_hcw[t] + f_hcw_pp*inf_p[t] + f_hcw_p*I_p[t];
    // Probability of infection for patients at time t
    p_p[t] = f_p_hcw*inf_hcw[t] + f_p_pp*inf_p[t] + f_p_p*I_p[t];
    
    // Number of newly infected HCWs at time t
    I_hcw[1,t] = binomial_rng(S_hcw[t], p_hcw[t]);
    // Number of newly infected patients at time t
    I_pU[1,t] = binomial_rng(S_p[t], p_p[t]);
  
    // Recoveries, discharges, isolations, deaths
    for(s in 2:t){
      // Number of infected HCWs at time t who got infected s-1 days ago (note that stan starts counting at 1)
      II_hcw[s,t] = (1-gamma[s])*I_hcw[s-1,t-1]); 
      I_hcw[s,t] = floorInt(II_hcw[s,t]);       // convert to int

      // Number of unknown infected patients at time t who got infected s-1 days ago 
      // = Number of unknown infected patients a day before - those that are discharged, isolated, or die
      II_pU[s,t] = (1-delta[s])*I_pU[s-1,t-1];
      I_pU[s,t] = floorInt(II_pU[s,t]); // convert to int
      
      // Returning (immune) HCWs at time t (here: real value)
      ww[t] += gamma[s]*I_hcw[s,t-1];
    }
    // Returning (immune) HCWs at time t, converted from real to int
    w[t] = floorInt(ww[t]);
    // Number of immune HCWs at time t
    R_hcw[t] = R_hcw[t-1] + w[t];

    N_ncp[t] = S_p[t] + I_pU[1,t]; // generate this for simulation 
  }

}
