
## Cluster Randomised Trial for a Typhoid vaccine in Vellore ##
    ## MRes Biomedical Research - EECID - Project 2 ##

                ## SIR Model Exploration ##



# SET UP ------------------------------------------------------------------

library(tidyverse)
library(ggplot2)



# 1: ONE CLUSTER ----------------------------------------------------------


# A simple SIR model for one cluster

rm(list = ls())


#### Time ####

time_step <- 1                                          # One day time step
time_seq <- seq(from = 1, to = 2*365, by = time_step)   # Total time: two years


#### Parameters ####

# beta = R0/Duration of infectiousness - infections/time

R0 <- 2       # Basic reproduction number - CHECK 
dur_inf <- 3 # Duration of infectiousness (days) - CHECK

beta <- R0/dur_inf
# Hazard of infection: beta x prevalence; prevalence = I/N
# Probability of infection: 1 - exp(-hazard*time)

# Recovery

1/dur_inf
# Hazard of recovery: 1/duration of infectiousness
# Prob of recovery: 1 - exp(-hazard*time)

# Vaccination

p_vax <- 0.5   # Vaccination coverage - THIS WE WILL EXPLORE
vax_eff <- 0.8     # Vaccine effect in infection reduction - CHECK

# Infected observed
# 
# p_sym <- 0.55       # Probability of being symptomatic - CHECK
# p_test <- 0.75      # Probability of seeking a test - CHECK
# p_positive <- 0.85  # Probability of test being positive
# 
# mu <- p_sym*p_test*p_positive


#### Population ####

number_infected <- 10

no_N <- as.integer(rep(2000, times = length(time_seq))) # Total cluster population

no_I <- as.integer(rep(10, times = length(time_seq)))   # Infected
no_R <- as.integer(rep(0, times = length(time_seq)))    # Recovered

no_S <- as.integer(rep((no_N[1]-no_I[1])*(1 - p_vax),
                       times = length(time_seq)))       # Susceptible (non-vax)
no_V <- as.integer(rep((no_N[1]-no_I[1])*p_vax,
                       times = length(time_seq)))       # Vaccinated


#### Loop ####

# Initialize state vars

# Hazards
haz_inf <- vector(mode = "numeric", length = length(time_seq))
haz_rec <- vector(mode = "numeric", length = length(time_seq))

# Probabilities: 1 - exp(hazard)
lambda <- vector(mode = "numeric", length = length(time_seq))
lambda_v <- vector(mode = "numeric", length = length(time_seq))

sigma <- vector(mode = "numeric", length = length(time_seq))

# State incidence/recovery: add stochastic
inc_SI <- vector(mode = "integer", length = length(time_seq))
inc_VI <- vector(mode = "integer", length = length(time_seq))

rec_IR <- vector(mode = "integer", length = length(time_seq))


for (i in 2:length(time_seq)) {
  
  # Hazards
  haz_inf[i] = beta*no_I[i-1]/no_N[i-1] 
  haz_rec[i] = 1/dur_inf
  
  # Probabilities
  lambda[i] = (1 - exp(-haz_inf[i]*time_step))
  lambda_v[i] = (1 - exp(-haz_inf[i]*(1 - vax_eff)*time_step))
  sigma[i] = (1 - exp(-haz_rec[i]*time_step))  
  
  # State variables
  inc_SI[i] = rbinom(n = 1, size = no_S[i-1], prob = lambda[i])
  inc_VI[i] = rbinom(n = 1, size = no_V[i-1], prob = lambda_v[i])
  rec_IR[i] = rbinom(n = 1, size = no_I[i-1], prob = sigma[i])  

  # Model equations
  no_S[i] = no_S[i-1] - inc_SI[i]
  no_V[i] = no_V[i-1] - inc_VI[i]
  
  no_I[i] = no_I[i-1] + inc_SI[i] + inc_VI[i] - rec_IR[i]
  
  no_R[i] = no_R[i-1] + rec_IR[i]
  
}

# Simple graph

sir <- data.frame(time_seq,
                  no_N, no_S, no_V, no_I, no_R)

ggplot(data = sir, mapping = aes(x = time_seq, y = no_I)) +
  geom_line(size = rel(1.2), color = "firebrick") +
  theme_classic() +
  labs(title = "Incidence over time",
       x = "Time over two years (days)",
       y = "Number of infections") +
  theme(
    plot.title = element_text(size = rel(1.2), face="bold", hjust = 0.5),
    axis.title.x = element_text(size = rel(1.1), face="bold"),
    axis.title.y = element_text(size = rel(1.1), face="bold"),
    axis.text = element_text(size=rel(1)))



# 2: VARIOUS CLUSTERS -----------------------------------------------------


# Add all the clusters in here

# Basic structure: an array with three dim
# 1: Time (no of rows)
# 2: Variables S, I, etc. (no of cols)
# 3: Cluster number (no of mini matrices)


rm(list = ls())


#### Parameters ####

# Starting seed

number_infected <- as.integer(10)

# beta = R0/Duration of infectiousness - infections/time

R0 <- 2       # Basic reproduction number - CHECK 
dur_inf <- 3 # Duration of infectiousness (days) - CHECK

beta <- R0/dur_inf
# Hazard of infection: beta x prevalence; prevalence = I/N
# Probability of infection: 1 - exp(-hazard*time)

# Recovery

1/dur_inf
# Hazard of recovery: 1/duration of infectiousness
# Prob of recovery: 1 - exp(-hazard*time)

# Vaccination

p_vax <- 0.5   # Vaccination coverage - THIS WE WILL EXPLORE
vax_eff <- 0.8     # Vaccine effect in infection reduction - CHECK

# Infected observed
# 
# p_sym <- 0.55       # Probability of being symptomatic - CHECK
# p_test <- 0.75      # Probability of seeking a test - CHECK
# p_positive <- 0.85  # Probability of test being positive
# 
# mu <- p_sym*p_test*p_positive


#### 1D: Time ####

time_step <- 1                                          # One day time step
time_seq <- seq(from = 1, to = 365*2, by = time_step)   # Total time: two years


#### 2D: Variables ####

# First indications

cluster <- rep(0, times = length(time_seq))

v_cluster <- rep(0, times = length(time_seq))

time_seq

# Empty vectors

no_N <- as.integer(rep(0, times = length(time_seq)))     # Total n per cluster

no_S <- as.integer(rep(0, times = length(time_seq)))     # Susceptible 
no_V <- as.integer(rep(0, times = length(time_seq)))     # Vaccinated  

no_I <- as.integer(rep(0, times = length(time_seq)))     # Infected
no_R <- as.integer(rep(0, times = length(time_seq)))     # Recovered

# Create state vars 

# Hazards
haz_inf <- rep(0, times = length(time_seq))
haz_rec <- rep(0, times = length(time_seq))

# Probabilities: 1 - exp(hazard)
lambda <- rep(0, times = length(time_seq))
lambda_v <- rep(0, times = length(time_seq))

sigma <- rep(0, times = length(time_seq))

# State incidence/recovery: add stochastic
inc_SI <- as.integer(rep(0, times = length(time_seq)))
inc_VI <- as.integer(rep(0, times = length(time_seq)))

rec_IR <- as.integer(rep(0, times = length(time_seq)))


#### 3D: Clusters ####

N <- 200000                          # Total population in the study: FIXED

C <- 100                              # Number of clusters
cluster_no <- seq(1:C)                # Vector of clusters

                                      # Population in each cluster
n <- N/C                              # (FOR THE MOMENT IT'S SAME - TO CHANGE)
cluster_n <- as.integer(rep(n, times = C))        # Vector of cluster populations


V <- C/2                              # Number of clusters in the vaccine group
cluster_vstatus <- c(rep(1, times = V), rep(0, times = V)) # Flag for vax clusters


#### Array ####

# Build

names_row <- paste0("t_", time_seq)

names_column <- c("cluster", "vaccine", "time_seq",
                  "no_N", "no_S", "no_V", "no_I", "no_R",
                  "haz_inf", "haz_rec",
                  "lambda", "lambda_v", "sigma",
                  "inc_SI", "inc_VI", "rec_IR")

names_matrix <- paste0("cluster_", cluster_no)

sir <- array(c(cluster, v_cluster, time_seq, no_N, no_S, no_V, no_I, no_R,
               haz_inf, haz_rec,lambda, lambda_v, sigma, inc_SI, inc_VI, rec_IR),
             dim = c(length(time_seq), 16, C),
             dimnames = list(names_row, names_column, names_matrix))
sir

# Assign the correct values for each column

# Cluster number, vaccine status and no_N

for (i in 1:C) {

  sir[, 1, i] = cluster_no[i]

  sir[, 2, i] = cluster_vstatus[i]

  sir[, 4, i] = cluster_n[i]

}

# I an R

for (i in 1:C) {
  
  sir[, 7, i] = number_infected
  
  sir[, 8, i] = 0

}

# S and V population - depends on how many are infected

for (i in 1:C) {
  
  if (sir[, 2, i] == 1) {                                                  # In vaccinated clusters
    sir[, 5, i] = round((sir[, 4 ,i]-sir[, 7, i])*(1 - p_vax), digits = 0) # Susceptible (non-vax)
    sir[, 6, i] = round((sir[, 4 ,i]-sir[, 7, i])*p_vax, digits = 0)       # Vaccinated
  
  } else {                                                                 # In non-vax clusters
    sir[, 5, i] = round((sir[, 4 ,i]-sir[, 7, i]), digits = 0)             # Susceptible (non-vax)
    sir[, 6, i] = as.integer(0)                                            # Vaccinated
  
  }
}


#### Loop over time ####

for (j in 1:C) {
  
  for (i in 2:length(time_seq)) {
      
    # Hazards
    sir[i, 9, j] = beta*sir[i-1, 7, j]/sir[i-1, 4 ,j] 
    sir[i, 10, j] = 1/dur_inf
      
    # Probabilities
    sir[i, 11, j] = (1 - exp(-sir[i, 9, j]*time_step))
    sir[i, 12, j] = (1 - exp(-sir[i, 9, j]*(1 - vax_eff)*time_step))
    sir[i, 13, j] = (1 - exp(-sir[i, 10, j]*time_step))  
      
    # State variables
    sir[i, 14, j] = round(rbinom(n = 1, size = sir[i-1, 5, j], prob = sir[i, 11, j]), digits = 0)
    sir[i, 15, j] = round(rbinom(n = 1, size = sir[i-1, 6, j], prob = sir[i, 12, j]), digits = 0)
    sir[i, 16, j] = round(rbinom(n = 1, size = sir[i-1, 7, j], prob = sir[i, 13, j]), digits = 0)  
      
    # Model equations
    sir[i, 5, j] = sir[i-1, 5, j] - sir[i, 14, j]
    sir[i, 6, j] = sir[i-1, 6, j] - sir[i, 15, j]
  
    sir[i, 7, j] = sir[i-1, 7, j] + sir[i, 14, j] + sir[i, 15, j] - sir[i, 16, j]
      
    sir[i, 8, j] = sir[i-1, 8, j] + sir[i, 16, j]
      
  }
}

# Store results

sir_result <- as.data.frame(matrix(0, nrow = length(time_seq), ncol = C))
colnames(sir_result) <- names_matrix

for (i in 1:C) {
  sir_result[,i] <- sir[,7,i]
}

# Simple graph!

sir_result <- sir_result %>%
  mutate(time_seq = time_seq) %>%
  pivot_longer(-time_seq, names_to = "cluster", values_to = "infected") %>%
  mutate(num_cluster = readr::parse_number(cluster)) %>%
  arrange(num_cluster)

ggplot() +
  geom_line(data = sir_result[1:36500,],
            mapping = aes(x = time_seq, y = infected, group = num_cluster),
            color = "steelblue") +
  geom_line(data = sir_result[36501:73000,],
            mapping = aes(x = time_seq, y = infected, group = num_cluster),
            color = "firebrick") +
  xlim(c(1, 150)) +
  theme_classic() +
  labs(title = "Incidence over time",
       x = "Time over two years (days)",
       y = "Number of infections") +
  theme(
    plot.title = element_text(size = rel(1.2), face="bold", hjust = 0.5),
    axis.title.x = element_text(size = rel(1.1), face="bold"),
    axis.title.y = element_text(size = rel(1.1), face="bold"),
    axis.text = element_text(size=rel(1)))



# 3: UNEQUAL CLUSTERS -----------------------------------------------------


# Add all the clusters in here

# Basic structure: an array with three dim
# 1: Time (no of rows)
# 2: Variables S, I, etc. (no of cols)
# 3: Cluster number (no of mini matrices)


rm(list = ls())


#### Parameters ####

# Starting seed

number_infected <- as.integer(10)

# beta = R0/Duration of infectiousness - infections/time

R0 <- 2       # Basic reproduction number - CHECK 
dur_inf <- 3 # Duration of infectiousness (days) - CHECK

beta <- R0/dur_inf
# Hazard of infection: beta x prevalence; prevalence = I/N
# Probability of infection: 1 - exp(-hazard*time)

# Recovery

1/dur_inf
# Hazard of recovery: 1/duration of infectiousness
# Prob of recovery: 1 - exp(-hazard*time)

# Vaccination

p_vax <- 0.5   # Vaccination coverage - THIS WE WILL EXPLORE
vax_eff <- 0.8     # Vaccine effect in infection reduction - CHECK

# Infected observed

p_sym <- 0.55       # Probability of being symptomatic - CHECK
p_test <- 0.50      # Probability of seeking a test - CHECK
p_positive <- 0.80  # Probability of test being positive

mu <- p_sym*p_test*p_positive


#### 1D: Time ####

time_step <- 1                                          # One day time step
time_seq <- seq(from = 1, to = 365*2, by = time_step)   # Total time: two years


#### 2D: Variables ####

# First indications

cluster <- rep(0, times = length(time_seq))

v_cluster <- rep(0, times = length(time_seq))

time_seq

# Empty vectors

no_N <- as.integer(rep(0, times = length(time_seq)))     # Total n per cluster

no_S <- as.integer(rep(0, times = length(time_seq)))     # Susceptible 
no_V <- as.integer(rep(0, times = length(time_seq)))     # Vaccinated  

no_I <- as.integer(rep(0, times = length(time_seq)))     # Infected
no_R <- as.integer(rep(0, times = length(time_seq)))     # Recovered

# Create state vars 

# Hazards
haz_inf <- rep(0, times = length(time_seq))
haz_rec <- rep(0, times = length(time_seq))

# Probabilities: 1 - exp(hazard)
lambda <- rep(0, times = length(time_seq))
lambda_v <- rep(0, times = length(time_seq))

sigma <- rep(0, times = length(time_seq))

# State incidence/recovery: add stochastic
inc_SI <- as.integer(rep(0, times = length(time_seq)))
inc_VI <- as.integer(rep(0, times = length(time_seq)))

rec_IR <- as.integer(rep(0, times = length(time_seq)))


#### 3D: Clusters ####

N <- 200000                          # Total population in the study: FIXED

C <- 100                              # Number of clusters
cluster_no <- seq(1:C)                # Vector of clusters

                                                            # Population in each cluster: normal distribution
n <- round(rnorm(n = C, mean = N/C, sd = 100), digits = 0)  #  with mean total people/number clusters
cluster_n <- abs(n)                                         # Vector of cluster populations

V <- C/2                                                   # Number of clusters in the vaccine group
cluster_vstatus <- c(rep(1, times = V), rep(0, times = V)) # Flag for vax clusters

# Data frame for reference

cluster_data <- data.frame(cluster = cluster_no,
                           vaccine = cluster_vstatus,
                           pop = cluster_n)


#### Array ####

# Build

names_row <- paste0("t_", time_seq)

names_column <- c("cluster", "vaccine", "time_seq",
                  "no_N", "no_S", "no_V", "no_I", "no_R",
                  "haz_inf", "haz_rec",
                  "lambda", "lambda_v", "sigma",
                  "inc_SI", "inc_VI", "rec_IR")

names_matrix <- paste0("cluster_", cluster_no)

sir <- array(c(cluster, v_cluster, time_seq, no_N, no_S, no_V, no_I, no_R,
               haz_inf, haz_rec,lambda, lambda_v, sigma, inc_SI, inc_VI, rec_IR),
             dim = c(length(time_seq), 16, C),
             dimnames = list(names_row, names_column, names_matrix))
sir

# Assign the correct values for each column

# Cluster number, vaccine status and no_N

for (i in 1:C) {
  
  sir[, 1, i] = cluster_no[i]
  
  sir[, 2, i] = cluster_vstatus[i]
  
  sir[, 4, i] = cluster_n[i]
  
}

# I an R

for (i in 1:C) {
  
  sir[, 7, i] = number_infected
  
  sir[, 8, i] = 0
  
}

# S and V population - depends on how many are infected

for (i in 1:C) {
  
  if (sir[, 2, i] == 1) {                                                  # In vaccinated clusters
    sir[, 5, i] = round((sir[, 4 ,i]-sir[, 7, i])*(1 - p_vax), digits = 0) # Susceptible (non-vax)
    sir[, 6, i] = round((sir[, 4 ,i]-sir[, 7, i])*p_vax, digits = 0)       # Vaccinated
    
  } else {                                                       # In non-vax clusters
    sir[, 5, i] = round((sir[, 4 ,i]-sir[, 7, i]), digits = 0)   # Susceptible (non-vax)
    sir[, 6, i] = as.integer(0)                                  # Vaccinated
    
  }
}


#### Loop over time ####

for (j in 1:C) {
  
  for (i in 2:length(time_seq)) {
    
    # Hazards
    sir[i, 9, j] = beta*sir[i-1, 7, j]/sir[i-1, 4 ,j] 
    sir[i, 10, j] = 1/dur_inf
    
    # Probabilities
    sir[i, 11, j] = (1 - exp(-sir[i, 9, j]*time_step))
    sir[i, 12, j] = (1 - exp(-sir[i, 9, j]*(1 - vax_eff)*time_step))
    sir[i, 13, j] = (1 - exp(-sir[i, 10, j])*time_step)  
    
    # State variables
    sir[i, 14, j] = round(rbinom(n = 1, size = sir[i-1, 5, j], prob = sir[i, 11, j]), digits = 0)
    sir[i, 15, j] = round(rbinom(n = 1, size = sir[i-1, 6, j], prob = sir[i, 12, j]), digits = 0)
    sir[i, 16, j] = round(rbinom(n = 1, size = sir[i-1, 7, j], prob = sir[i, 13, j]), digits = 0)  
    
    # Model equations
    sir[i, 5, j] = sir[i-1, 5, j] - sir[i, 14, j]
    sir[i, 6, j] = sir[i-1, 6, j] - sir[i, 15, j]
    
    sir[i, 7, j] = sir[i-1, 7, j] + sir[i, 14, j] + sir[i, 15, j] - sir[i, 16, j]
    
    sir[i, 8, j] = sir[i-1, 8, j] + sir[i, 16, j]
    
  }
}


#### Results: all ####

# Store results

sir_res_infected <- as.data.frame(matrix(0, # Empty matrix
                                         nrow = length(time_seq), ncol = C))
colnames(sir_res_infected) <- names_matrix

for (i in 1:C) {
  sir_res_infected[,i] <- sir[,7,i]         # Get results from run
}

sir_res_observed <- as.data.frame(matrix(0, # Empty matrix
                                         nrow = length(time_seq), ncol = C))
colnames(sir_res_observed) <- names_matrix

for (i in 1:C) {
  sir_res_observed[,i] <- round(sir[,7,i]*mu, digits = 0) # Results from run x prob of detecting that infected
}

# Pivot the results and merge all together (infected and observed)

sir_res_infected <- sir_res_infected %>%
  mutate(time_seq = time_seq) %>%
  pivot_longer(-time_seq, names_to = "cluster", values_to = "infected") %>%
  mutate(cluster = readr::parse_number(cluster)) %>%
  arrange(cluster) %>%
  merge(y = cluster_data[, c("cluster", "vaccine")], by = "cluster", all = TRUE) %>%
  select(c("cluster", "vaccine", "time_seq", "infected"))

sir_res_observed <- sir_res_observed %>%
  mutate(time_seq = time_seq) %>%
  pivot_longer(-time_seq, names_to = "cluster", values_to = "observed") %>%
  mutate(cluster = readr::parse_number(cluster)) %>%
  arrange(cluster) %>%
  merge(y = cluster_data[, c("cluster", "vaccine")], by = "cluster", all = TRUE) %>%
  select(c("cluster", "vaccine", "time_seq", "observed"))

sir_result <- merge(sir_res_infected, sir_res_observed,
                    by = c("cluster", "vaccine", "time_seq"), all = TRUE)
sir_result <- arrange(sir_result, cluster, time_seq)

# Simple graph (mainly to check if model was OK)

ggplot() +
  geom_line(data = filter(sir_result, vaccine == 1),
            mapping = aes(x = time_seq, y = infected, group = cluster,
            color = "Inf_Vax")) +
  geom_line(data = filter(sir_result, vaccine == 0),
            mapping = aes(x = time_seq, y = infected, group = cluster,
            color = "Inf_No")) +
  geom_line(data = filter(sir_result, vaccine == 1),
            mapping = aes(x = time_seq, y = observed, group = cluster,
            color = "Obs_Vax")) +
  geom_line(data = filter(sir_result, vaccine == 0),
            mapping = aes(x = time_seq, y = observed, group = cluster,
            color = "Obs_No")) +
  scale_color_manual(name = NULL,
                     breaks = c("Inf_Vax", "Inf_No", "Obs_Vax", "Obs_No"),
                     values = c("Inf_Vax" = "steelblue4",
                                "Inf_No" = "steelblue1",
                                "Obs_Vax" = "palegreen4",
                                "Obs_No" = "palegreen1"),
                     labels = c("Infections in vaccine clusters",
                                "Infections in non-vaccine clusters",
                                "Detected infections in vaccine clusters",
                                "Detected infections in non-vaccine cluster")) +
  xlim(c(1, 150)) +
  theme_classic() +
  labs(title = "Incidence over time",
       x = "Time over two years (days)",
       y = "Number of infections/detected infections") +
  theme(
      plot.title = element_text(size = rel(1.2), face="bold", hjust = 0.5),
      axis.title.x = element_text(size = rel(1.1), face="bold"),
      axis.title.y = element_text(size = rel(1.1), face="bold"),
      axis.text = element_text(size=rel(1)),
      legend.position = "bottom",
      legend.text = element_text(size=rel(1)))


#### Results: summary ####

# Get summary stats

sir_summary <- sir_result %>%
  group_by(cluster) %>%
  mutate(mean_inf = mean(infected)) %>%
  mutate(mena_obs = mean(observed)) %>%
  mutate(sum_inf = sum(infected)) %>%
  mutate(sum_obs = sum(observed)) %>%
  filter(row_number() == 1) %>%
  select(-time_seq, -infected, -observed) %>%
  ungroup()

# Test

model <- glm(formula = sum_obs ~ vaccine,
             family = "poisson", data = sir_summary)
x <- exp(summary(model)$coef)
y <- exp(confint(model))



# 4: MORE CLUSTERS -------------------------------------------------------


# Array for all clusters and timepoints

# 1: Time (no of rows)
# 2: Variables S, I, etc. (no of cols)
# 3: Cluster number (no of mini matrices)

# Over three loops

# 1: Loop over all clusters
# 2: Inside each cluster, loop over time
# 3: At every time point, when calculating FOI, add external all other cluster FOI 


rm(list = ls())


#### Parameters ####

# Starting seed

number_infected <- as.integer(1)

# FOI: beta = R0/Duration of infectiousness - infections/time

R0 <- 2       # Basic reproduction number - CHECK 
dur_inf <- 3 # Duration of infectiousness (days) - CHECK

beta <- R0/dur_inf
# Hazard of infection: beta x prevalence; prevalence = I/N
# Probability of infection: 1 - exp(-hazard*time)

imp_rate <- 0.25
# Importation rate from other cluster, i.e.
# Even if clusters are right next to each other,
# we assume not everything comes in.


# Recovery

1/dur_inf
# Hazard of recovery: 1/duration of infectiousness
# Prob of recovery: 1 - exp(-hazard*time)

# Vaccination

p_vax <- 0.5       # Vaccination coverage - THIS WE WILL EXPLORE
p_clusvax <- 0.5   # Proportion of clusters assigned to vaccine
vax_eff <- 0.8     # Vaccine effect in infection reduction - CHECK

# Infected observed

p_sym <- 0.55       # Probability of being symptomatic - CHECK
p_test <- 0.50      # Probability of seeking a test - CHECK
p_positive <- 0.80  # Probability of test being positive

mu <- p_sym*p_test*p_positive


#### 1D: Time ####

time_step <- 1                                          # One day time step
years <- 1                                              # Total time in years
time_seq <- seq(from = 1, to = 365*years, by = time_step)   # Sequence


#### 2D: Variables ####

# First indications

cluster <- rep(0, times = length(time_seq))

v_cluster <- rep(0, times = length(time_seq))

time_seq

# Empty vectors

no_N <- as.integer(rep(0, times = length(time_seq)))     # Total n per cluster

no_S <- as.integer(rep(0, times = length(time_seq)))     # Susceptible 
no_V <- as.integer(rep(0, times = length(time_seq)))     # Vaccinated  

no_I <- as.integer(rep(0, times = length(time_seq)))     # Infected
no_R <- as.integer(rep(0, times = length(time_seq)))     # Recovered

# Create state vars 

# Hazards
haz_inf <- rep(0, times = length(time_seq))
haz_rec <- rep(0, times = length(time_seq))

# Probabilities: 1 - exp(hazard)
lambda <- rep(0, times = length(time_seq))
lambda_v <- rep(0, times = length(time_seq))

sigma <- rep(0, times = length(time_seq))

# State incidence/recovery: add stochastic
inc_SI <- as.integer(rep(0, times = length(time_seq)))
inc_VI <- as.integer(rep(0, times = length(time_seq)))

rec_IR <- as.integer(rep(0, times = length(time_seq)))


#### 3D: Clusters ####

N <- 20000                            # Total population in the study: FIXED

C <- 10                               # Number of clusters
cluster_no <- seq(1:C)                # Vector of clusters

# Population in each cluster: normal distribution
n <- round(rnorm(n = C, mean = N/C, sd = 100), digits = 0)  #  with mean total people/number clusters
cluster_n <- abs(n)                                         # Vector of cluster populations

V <- C*p_clusvax                                           # Number of clusters in the vaccine group
cluster_vstatus <- c(rep(1, times = V), rep(0, times = V)) # Flag for vax clusters

#  Distance matrix

cluster_dis <- matrix(1, nrow = C, ncol = C,
                      dimnames = list(cluster_no, cluster_no))
# Created a matrix that states all distance between each pair of clusters
# From 1: no distance (same cluster) to 100 (far away)

cluster_dis[lower.tri(cluster_dis, diag = FALSE)] <- runif(
  n = (C^2 - C)/2, min = 1, max = 100)
cluster_dis[upper.tri(cluster_dis, diag = FALSE)] <- t(cluster_dis)[upper.tri(cluster_dis)]
cluster_dis <- round(cluster_dis, digits = 1)
# Filled with random numbers in a symmetrical way

# Data frame for reference

cluster_data <- data.frame(cluster = cluster_no,
                           vaccine = cluster_vstatus,
                           pop = cluster_n,
                           cluster_dis)
colnames(cluster_data) <- c("cluster", "vaccine", "pop",
                            paste0("dis_", cluster_no))


#### Array ####

# Build

names_row <- paste0("t_", time_seq)

names_column <- c("cluster", "vaccine", "time_seq",
                  "no_N", "no_S", "no_V", "no_I", "no_R",
                  "haz_inf", "haz_rec",
                  "lambda", "lambda_v", "sigma",
                  "inc_SI", "inc_VI", "rec_IR")

names_matrix <- paste0("cluster_", cluster_no)

sir <- array(c(cluster, v_cluster, time_seq, no_N, no_S, no_V, no_I, no_R,
               haz_inf, haz_rec,lambda, lambda_v, sigma, inc_SI, inc_VI, rec_IR),
             dim = c(length(time_seq), 16, C),
             dimnames = list(names_row, names_column, names_matrix))
sir

# Assign the correct values for each column

# Cluster number, vaccine status and no_N

for (i in 1:C) {
  
  sir[, 1, i] = cluster_no[i]
  
  sir[, 2, i] = cluster_vstatus[i]
  
  sir[, 4, i] = cluster_n[i]
  
}

# I an R

for (i in 1:C) {
  
  sir[, 7, i] = number_infected
  
  sir[, 8, i] = 0
  
}

# S and V population - depends on how many are infected

for (i in 1:C) {
  
  if (sir[, 2, i] == 1) {                                                  # In vaccinated clusters
    sir[, 5, i] = round((sir[, 4 ,i]-sir[, 7, i])*(1 - p_vax), digits = 0) # Susceptible (non-vax)
    sir[, 6, i] = round((sir[, 4 ,i]-sir[, 7, i])*p_vax, digits = 0)       # Vaccinated
    
  } else {                                                       # In non-vax clusters
    sir[, 5, i] = round((sir[, 4 ,i]-sir[, 7, i]), digits = 0)   # Susceptible (non-vax)
    sir[, 6, i] = as.integer(0)                                  # Vaccinated
    
  }
}
# Ignore warnings at the end


#### Loop over time ####

for (j in 1:C) {
  
  for (i in 2:length(time_seq)) {
    
    # Hazards
    sir[i, 9, j] = beta*sir[i-1, 7, j]/sir[i-1, 4 ,j]
    
    for (k in 1:C) { # Loop to add external FOI
      
      if (k != j) {
        sir[i, 9, j] = sir[i, 9, j] + (imp_rate/cluster_dis[k,j])*beta*sir[i-1, 7, k]/sir[i-1, 4 ,k]
      }
      
    }
    
    sir[i, 10, j] = 1/dur_inf
    
    # Probabilities
    sir[i, 11, j] = (1 - exp(-sir[i, 9, j]*time_step))
    sir[i, 12, j] = (1 - exp(-sir[i, 9, j]*(1 - vax_eff)*time_step))
    sir[i, 13, j] = (1 - exp(-sir[i, 10, j])*time_step)  
    
    # State variables
    sir[i, 14, j] = round(rbinom(n = 1, size = sir[i-1, 5, j], prob = sir[i, 11, j]), digits = 0)
    sir[i, 15, j] = round(rbinom(n = 1, size = sir[i-1, 6, j], prob = sir[i, 12, j]), digits = 0)
    sir[i, 16, j] = round(rbinom(n = 1, size = sir[i-1, 7, j], prob = sir[i, 13, j]), digits = 0)  
    
    # Model equations
    sir[i, 5, j] = sir[i-1, 5, j] - sir[i, 14, j]
    sir[i, 6, j] = sir[i-1, 6, j] - sir[i, 15, j]
    
    sir[i, 7, j] = sir[i-1, 7, j] + sir[i, 14, j] + sir[i, 15, j] - sir[i, 16, j]
    
    sir[i, 8, j] = sir[i-1, 8, j] + sir[i, 16, j]
    
  }
}


#### Results: all ####

# Store results

sir_res_infected <- as.data.frame(matrix(0, # Empty matrix
                                         nrow = length(time_seq), ncol = C))
colnames(sir_res_infected) <- names_matrix

for (i in 1:C) {
  sir_res_infected[,i] <- sir[,7,i]         # Get results from run
}

sir_res_observed <- as.data.frame(matrix(0, # Empty matrix
                                         nrow = length(time_seq), ncol = C))
colnames(sir_res_observed) <- names_matrix

for (i in 1:C) {
  sir_res_observed[,i] <- round(sir[,7,i]*mu, digits = 0) # Results from run x prob of detecting that infected
}

# Pivot the results and merge all together (infected and observed)

sir_res_infected <- sir_res_infected %>%
  mutate(time_seq = time_seq) %>%
  pivot_longer(-time_seq, names_to = "cluster", values_to = "infected") %>%
  mutate(cluster = readr::parse_number(cluster)) %>%
  arrange(cluster) %>%
  merge(y = cluster_data[, c("cluster", "vaccine")], by = "cluster", all = TRUE) %>%
  select(c("cluster", "vaccine", "time_seq", "infected"))

sir_res_observed <- sir_res_observed %>%
  mutate(time_seq = time_seq) %>%
  pivot_longer(-time_seq, names_to = "cluster", values_to = "observed") %>%
  mutate(cluster = readr::parse_number(cluster)) %>%
  arrange(cluster) %>%
  merge(y = cluster_data[, c("cluster", "vaccine")], by = "cluster", all = TRUE) %>%
  select(c("cluster", "vaccine", "time_seq", "observed"))

sir_result <- merge(sir_res_infected, sir_res_observed,
                    by = c("cluster", "vaccine", "time_seq"), all = TRUE)
sir_result <- arrange(sir_result, cluster, time_seq)

# Simple graph (mainly to check if model was OK)

ggplot() +
  geom_line(data = filter(sir_result, vaccine == 1),
            mapping = aes(x = time_seq, y = infected, group = cluster,
                          color = "Inf_Vax"), size = rel(1.1)) +
  geom_line(data = filter(sir_result, vaccine == 0),
            mapping = aes(x = time_seq, y = infected, group = cluster,
                          color = "Inf_No"), size = rel(1.1)) +
  geom_line(data = filter(sir_result, vaccine == 1),
            mapping = aes(x = time_seq, y = observed, group = cluster,
                          color = "Obs_Vax"), size = rel(1.1)) +
  geom_line(data = filter(sir_result, vaccine == 0),
            mapping = aes(x = time_seq, y = observed, group = cluster,
                          color = "Obs_No"), size = rel(1.1)) +
  scale_color_manual(name = NULL,
                     breaks = c("Inf_Vax", "Inf_No", "Obs_Vax", "Obs_No"),
                     values = c("Inf_Vax" = "steelblue4",
                                "Inf_No" = "steelblue1",
                                "Obs_Vax" = "palegreen4",
                                "Obs_No" = "palegreen1"),
                     labels = c("Infections in vaccine clusters",
                                "Infections in non-vaccine clusters",
                                "Detected infections in vaccine clusters",
                                "Detected infections in non-vaccine cluster")) +
  xlim(c(1, 150)) +
  theme_classic() +
  labs(title = "Incidence over time",
       x = "Time over two years (days)",
       y = "Number of infections/detected infections") +
  theme(
    plot.title = element_text(size = rel(1.2), face="bold", hjust = 0.5),
    axis.title.x = element_text(size = rel(1.1), face="bold"),
    axis.title.y = element_text(size = rel(1.1), face="bold"),
    axis.text = element_text(size=rel(1)),
    legend.position = "bottom",
    legend.text = element_text(size=rel(1)))


#### Results: summary ####

# Get summary stats

sir_summary <- sir_result %>%
  group_by(cluster) %>%
  mutate(mean_inf = mean(infected)) %>%
  mutate(mena_obs = mean(observed)) %>%
  mutate(sum_inf = sum(infected)) %>%
  mutate(sum_obs = sum(observed)) %>%
  filter(row_number() == 1) %>%
  select(-time_seq, -infected, -observed) %>%
  ungroup()

# Test

model <- glm(formula = sum_obs ~ vaccine,
             family = "poisson", data = sir_summary)
x <- exp(summary(model)$coef)
y <- exp(confint(model))