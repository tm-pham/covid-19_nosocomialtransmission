# =================================================================== #
# Functions
# =================================================================== #
# N = total number
# S = Susceptible
# E = Exposed
# I_A = asymptomatic
# I_P = pre-symptomatic
# I_S = symptomatic
# R = Recovered
# EE = empty beds
# t = time frame
create.data.template <- function(N,S,E,I_A,I_P,I_S,R,EE=0,t){
  data <- expand.grid(subject = 1:N, time = 1:t) %>% 
    tbl_df()
  # 'wide' version of the data frame
  wdata <- data %>% mutate(shedding = NA_character_) %>% 
    spread(subject,shedding) %>% setNames(names(.) %>% make.names) 
  
  # Intitial values for different compartments
  exposed <- c(rep(1,E),rep(0,N-E)) 
  susceptible <- c(rep(0,E),rep(1,S),rep(0,N-E-S))
  asymptomatic <- presymptomatic <- recovered <- rep(0,N)
  symptomatic <- c(rep(0,E),rep(0,S),rep(1,I_S),rep(0,N-S-E-I_S))
  nopat <- c(rep(0,N-EE),rep(1,EE))
  
  wdata[1,-1][susceptible ==1] = 'S'
  wdata[1,-1][exposed == 1] = 'E'
  wdata[1,-1][asymptomatic == 1] = 'I_A'
  wdata[1,-1][presymptomatic == 1] = 'I_P'
  wdata[1,-1][symptomatic == 1] = 'I_S'
  wdata[1,-1][recovered == 1] = 'R'
  wdata[1,-1][nopat == 1] = 'EE'
  
  wdata[2:nrow(wdata),-1] <- wdata[1,-1]
  
  # data frame to save symptom onset and infection time
  ddata <- as.data.frame(cbind(subject=1:N,symptom_time=rep(0,N), inf_time=rep(0,N)))
  
  return(list(data=ddata,wdata=wdata))
}


# Generation time according to Ferretti et al (2020)
gen.time <- function(x,shape=2.826, scale=5.665){
  return(dweibull(x,shape=shape, scale=scale))
} 

# Beta binomial distribution 
rbbinom <- function(n, prob, k, size){
  mtilde <- rbeta(n, k/(1-prob), k/prob)
  return(rbinom(n, prob=mtilde, size=size))
}
