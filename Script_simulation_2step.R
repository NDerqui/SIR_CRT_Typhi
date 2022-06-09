
## Cluster Randomised Trial for a Typhoid vaccine in Vellore ##
      ## MRes Biomedical Research - EECID - Project 2 ##

                  ## Two step simulation ##

## Run SIR open model without vaccination to reach infection equilibrium ##
    ## Then, allocate vaccine to some (1/2) clusters and follow  ##



# SET UP ------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(sandwich)
library(lmtest)
library(statnet)
library(ggnet)



# PARAMETERS --------------------------------------------------------------


# Parameters to be fixed across all runs


#### Incidence, birth and death ####

# Prevalence of infection in India
# Incidence according to GDB 2017

incidence <- 0.005492 # 549.2 cases per 100,000

# Birth rate
# The World Bank Data 2020

# birth <- 0.017      # 17 per 1,000
birth <- 0.007

# Death rate
# The World Bank Data 2020

# death <- 0.007        # 7 per 1,000
death <- 0.007        # 7 per 1,000


#### Force of infection ####

# (infections/time): beta = R0/Duration of infectiousness

R0 <- 2            # Basic reproduction number
dur_inf <- 7       # Duration of infectiousness (days)

# Importation rate: from external clusters to a given one (??)

imp_rate <- 0.5
imp_rate <- 0.25


#### Vaccination ####

p_vax <- 0.5       # Proportion of vaccinated population in vaccine clusters
p_clusvax <- 0.5   # Proportion of clusters assigned to vaccine
vax_eff <- 0.8     # Vaccine effectiveness (infection)


#### Detected infections ####

p_sym <- 0.55       # Probability of being symptomatic - CHECK
p_test <- 0.50      # Probability of seeking a test - CHECK
p_positive <- 0.60  # Probability of test being positive


#### Time frame ####

time_step <- 1      # Time step change (days)


#### Population ####

# Total population in the study and number of clusters

N <- 200000                       # Total population in the study
C <- 100                          # Number of clusters

n <- round(rnorm(n = C, mean = N/C, sd = 100), digits = 0)  # Pop in each cluster
cluster_n <- abs(n)                                         # Vector of cluster populations
hist(cluster_n, main = "Histogram of clusters' population",
     xlab = "Clusters' population", ylab = "Frequency of clusters")


#### Clusters ####

# Vector of clusters

cluster_no <- seq(1:C)

# Vaccination status of clusters

V <- C*p_clusvax                                           # Number of clusters in the vaccine group
cluster_vstatus <- c(rep(1, times = V), rep(0, times = V)) # Flag for vax clusters


# Cluster map
# Random location of clusters (cannot be in order, we vax the first proportion)

cluster_map <- matrix(sample(cluster_no), ncol = sqrt(C), nrow = sqrt(C),
                      dimnames = list(seq(1:sqrt(C)), seq(1:sqrt(C))))

# Cluster distance matrix

cluster_dis <- matrix(1, nrow = C, ncol = C,
                      dimnames = list(seq(1:C), seq(1:C)))

for (i in 1:C) {                          # Each cluster distance
  for (j in 1:C) {                        # with each of the others
    
    for (k in 1:sqrt(C)) {                # Look for location of first cluster
      for (l in 1:sqrt(C)) {
        if (i == cluster_map[k, l]) {
          
          for (m in 1:sqrt(C)) {          # Look for location of second cluster
            for (n in 1:sqrt(C)) {
              if (j == cluster_map[m, n]) {
                
                # Function to get the distance
                
                vertical <- abs(m - k)
                horizontal <- abs(n - l)
                distance <- sqrt(horizontal^2 + vertical^2)
                
                cluster_dis[i, j] <- round(distance, digits = 3)
              }
            }
          }
          
        }
      }
    }
    
  }
}
diag(cluster_dis) <- 1

# Network object

cluster_network <- as.network(x = cluster_dis,
                              directed = FALSE,
                              loops = FALSE,
                              matrix.type = "adjacency"
)
set.vertex.attribute(cluster_network,"Number", cluster_no)
set.vertex.attribute(cluster_network,"Vaccine", cluster_vstatus)
set.vertex.attribute(cluster_network,"Pop", cluster_n)
# Cluster(node) level attribute (let's put the number, pop and vax status)
edge_values <- as.numeric(1/cluster_dis)
set.edge.value(cluster_network, "Distance", edge_values)
# Edge(connection) level attribute (inversely proportion to distance)

summary.network(cluster_network,print.adj = FALSE)



# EQUILIBRIUM -------------------------------------------------------------


# Run open cohort SIR model for x years without the vaccine to reach equilibrium
# Birth rate will increase the number of S
# Death will be a competing risk
# Results will be used as starting point for the CRT simulation

equilibrium <- function(N, C, cluster_no, cluster_n, cluster_vstatus, cluster_dis,
                        incidence, R0, dur_inf, imp_rate,
                        time_step, years) {
  
  
  time_seq <- seq(from = 1, to = 365*years, by = time_step) # Total time   
  
  
  # Calculated parameters
  
  beta <- R0/dur_inf                        # Force of infection
  

  # Cluster vectors
  
  cluster_no <- cluster_no                  # Vector of clusters
  cluster_n <- cluster_n                    # Vector of cluster populations
  cluster_vstatus <- cluster_vstatus    # Flag for vax clusters
  

  # Empty columns for the array
  
  ## Basic cluster indications
  
  cluster <- rep(0, times = length(time_seq))
  v_cluster <- rep(0, times = length(time_seq))
  #time_seq is the second col
  
  ## SIR model compartments
  
  no_N <- as.integer(rep(0, times = length(time_seq)))     # Total n per cluster
  no_S <- as.integer(rep(0, times = length(time_seq)))     # Susceptible 
  no_I <- as.integer(rep(0, times = length(time_seq)))     # Infected 
  no_R <- as.integer(rep(0, times = length(time_seq)))     # Recovered
  
  ## From S to I
  
  # Hazard
  haz_inf <- rep(0, times = length(time_seq))
  # Probability: 1 - exp(hazard)
  prob_infec <- rep(0, times = length(time_seq))
  # State variables
  inc_SI <- as.integer(rep(0, times = length(time_seq)))
  
  ## From I to R
  
  # Hazard
  haz_rec <- rep(0, times = length(time_seq))
  # Probability: 1 - exp(hazard)
  prob_recov <- rep(0, times = length(time_seq))
  # State variable
  inc_IR <- as.integer(rep(0, times = length(time_seq)))
  
  ## Deaths
  
  # Hazard
  haz_dea <- rep(0, times = length(time_seq))
  # Probability: 1 - exp(hazard)
  prob_death <- rep(0, times = length(time_seq))
  # State variable
  inc_SD <- as.integer(rep(0, times = length(time_seq)))
  
  ## Births
  
  # Hazard
  haz_bir <- rep(0, times = length(time_seq))
  # Probability: 1 -exp(hazard)
  prob_birth <- rep(0, times = length(time_seq))
  # State variable
  inc_BS <- as.integer(rep(0, times = length(time_seq)))
  
  
  # Build the array
  
  ## Basic structure
  
  names_row <- paste0("t_", time_seq)
  
  names_column <- c("cluster", "vax_status", "time_seq", 
                    "no_N", "no_S", "no_I", "no_R",
                    "haz_inf","prob_infec", "inc_SI",
                    "haz_rec", "prob_recov", "inc_IR",
                    "haz_dea", "prob_death", "inc_SD",
                    "haz_bir", "prob_birth", "inc_BS")
  
  names_matrix <- paste0("cluster_", cluster_no)
  
  sir <- array(c(cluster, v_cluster, time_seq,
                 no_N, no_S, no_I, no_R,
                 haz_inf, prob_infec, inc_SI,
                 haz_rec, prob_recov, inc_IR,
                 haz_dea, prob_death, inc_SD,
                 haz_bir, prob_birth, inc_BS),
               dim = c(length(time_seq), 19, C),
               dimnames = list(names_row, names_column, names_matrix))
  
  ## Assign initial values
  
  # Cluster number and no_N
  for (i in 1:C) {
    sir[, 1, i] = cluster_no[i]
    sir[, 2, i] = cluster_vstatus[i]
    sir[, 4, i] = cluster_n[i]
  }
  # S, I, D and R
  for (i in 1:C) {
    sir[, 6, i] = round(incidence*sir[, 4, i], digits = 0)     # Initial I depends on incidence rate
    sir[, 5, i] = round((sir[, 4 ,i]-sir[, 6, i]), digits = 0) # Initial S depends on N - I
    sir[, 7, i] = 0
  }

  
  # Put the model to run
  
  for (j in 1:C) {
    
    for (i in 2:length(time_seq)) {
      
      # From S to I
      
      sir[i, 8, j] = beta*sir[i-1, 6, j]/sir[i-1, 4, j]     # Initial FOI: from cluster
      for (k in 1:C) {                                      # Loop to add external FOI
        if (k != j) {
          sir[i, 8, j] = sir[i, 8, j] + (imp_rate/cluster_dis[k,j])*beta*sir[i-1, 6, k]/sir[i-1, 4 ,k] 
        }
      }
      
      sir[i, 9, j] = (1 - exp(-sir[i, 8, j]*time_step))
      sir[i, 10, j] = round(rbinom(n = 1, size = sir[i-1, 5, j], prob = sir[i, 9, j]), digits = 0)
      
      # From I to R
      
      sir[i, 11, j] = 1/dur_inf
      sir[i, 12, j] = (1 - exp(-sir[i, 11, j])*time_step)  
      sir[i, 13, j] = round(rbinom(n = 1, size = sir[i-1, 6, j], prob = sir[i, 12, j]), digits = 0) 
      
      # Deaths
      
      sir[i, 14, j] = death
      sir[i, 15, j] = (1 - exp(-sir[i, 14, j])*time_step)  
      sir[i, 16, j] = round(rbinom(n = 1, size = sir[i-1, 5, j], prob = sir[i, 15, j]), digits = 0) 
      
      # Births
      
      sir[i, 17, j] = birth
      sir[i, 18, j] = (1 - exp(-sir[i, 17, j])*time_step)  
      sir[i, 19, j] = round(rbinom(n = 1, size = sir[i-1, 4, j], prob = sir[i, 18, j]), digits = 0) 
      
      # Model equations
      sir[i, 5, j] = sir[i-1, 5, j] - sir[i, 10, j] - sir[i, 16, j] + sir[i, 19, j]
      sir[i, 6, j] = sir[i-1, 6, j] + sir[i, 10, j] - sir[i, 13, j]
      sir[i, 7, j] = sir[i-1, 7, j] + sir[i, 13, j]
      #sir[i, 4, j] = sir[i-1, 4, j] - sir[i, 16, j] + sir[i, 19, j]
      
    }
  }
  
  
  # Store the results
  
  sir_result <- as.data.frame(matrix(0, nrow = C, ncol = 19))
  
  ## Last row of array
  
  for (i in 1:C) {
    sir_result[i,] <- sir[length(time_seq),,i]         
  }
  
  ## Returned object
  
  sir_result
}



# CRT: FUNCTIONS ----------------------------------------------------------


#### Main functions ####

sis_model <- function(N, C, cluster_n, cluster_dis,
                      R0, dur_inf, imp_rate,
                      p_vax, p_clusvax, vax_eff,
                      p_sym, p_test, p_positive,
                      time_step, years, equilibrium) {

  
  time_seq <- seq(from = 1, to = 365*years, by = time_step) # Total time   
    
  
  # Calculated parameters
  
  beta <- R0/dur_inf                                    # Force of infection
  mu <- p_sym*p_test*p_positive                         # Prob of detecting I
  
  
  # Cluster vectors

  cluster_no <- seq(1:C)                                     # Vector of clusters
  cluster_n <- cluster_n                                     # Vector of cluster populations
  V <- C*p_clusvax                                           # Number of clusters in the vaccine group
  cluster_vstatus <- c(rep(1, times = V), rep(0, times = V)) # Flag for vax clusters
  
  ## Data frame for reference
  
  cluster_data <- data.frame(cluster = cluster_no,
                             vaccine = cluster_vstatus,
                             pop = cluster_n,
                             cluster_dis)
  colnames(cluster_data) <- c("cluster", "vaccine", "pop",
                              paste0("dis_", cluster_no))
  

  # Empty columns for the array
  
  ## Basic cluster indications
  
  cluster <- rep(0, times = length(time_seq))
  v_cluster <- rep(0, times = length(time_seq))
  #time_seq is the third col
  
  ## SIS model compartments
  
  no_N <- as.integer(rep(0, times = length(time_seq)))     # Total n per cluster
  no_S <- as.integer(rep(0, times = length(time_seq)))     # Susceptible 
  no_V <- as.integer(rep(0, times = length(time_seq)))     # Vaccinated  
  no_I <- as.integer(rep(0, times = length(time_seq)))     # Infected
  no_R <- as.integer(rep(0, times = length(time_seq)))     # Recovered
  
  ## Hazards 
  
  haz_inf <- rep(0, times = length(time_seq))
  haz_rec <- rep(0, times = length(time_seq))
  
  ##  Probabilities: 1 - exp(hazard)
  
  lambda <- rep(0, times = length(time_seq))
  lambda_v <- rep(0, times = length(time_seq))
  sigma <- rep(0, times = length(time_seq))
  
  ## State variables: incidence/recovery
  
  inc_SI <- as.integer(rep(0, times = length(time_seq)))
  inc_VI <- as.integer(rep(0, times = length(time_seq)))
  rec_IR <- as.integer(rep(0, times = length(time_seq)))
  
  
  # Build the array
  
  ## Basic structure
  
  names_row <- paste0("t_", time_seq)
  names_column <- c("cluster", "vaccine", "time_seq",
                    "no_N", "no_S", "no_V", "no_I", "no_R",
                    "haz_inf", "haz_rec",
                    "lambda", "lambda_v", "sigma",
                    "inc_SI", "inc_VI", "rec_IR")
  names_matrix <- paste0("cluster_", cluster_no)
  sis <- array(c(cluster, v_cluster, time_seq, no_N, no_S, no_V, no_I, no_R,
                 haz_inf, haz_rec,lambda, lambda_v, sigma, inc_SI, inc_VI, rec_IR),
               dim = c(length(time_seq), 16, C),
               dimnames = list(names_row, names_column, names_matrix))
  
  ## Assign initial values
  
  # Cluster number, vaccine status and no_N
  for (i in 1:C) {
    sis[, 1, i] = cluster_no[i]
    sis[, 2, i] = cluster_vstatus[i]
    sis[, 4, i] = cluster_n[i]
  }
  # I an R
  for (i in 1:C) {
    sis[, 7, i] = equilibrium[i, 7]
    sis[, 8, i] = 0
  }
  # S and V population - depends on vaccination cluster status and coverage
  for (i in 1:C) {
    if (sis[1, 2, i] == 1) {                                                 # In vaccinated clusters
      sis[, 5, i] = round((sis[, 4 ,i]-sis[, 7, i])*(1 - p_vax), digits = 0) # Susceptible (non-vax)
      sis[, 6, i] = round((sis[, 4 ,i]-sis[, 7, i])*p_vax, digits = 0)       # Vaccinated
      
    } else {                                                       # In non-vax clusters
      sis[, 5, i] = round((sis[, 4 ,i]-sis[, 7, i]), digits = 0)   # Susceptible (non-vax)
      sis[, 6, i] = as.integer(0)                                  # Vaccinated
    }
  }
  

  # Put the model to run
    
  for (j in 1:C) {
    
    for (i in 2:length(time_seq)) {
      
      # Hazards
      sis[i, 9, j] = beta*sis[i-1, 7, j]/sis[i-1, 4 ,j]
      for (k in 1:C) { # Loop to add external FOI
        
        if (k != j) {
          sis[i, 9, j] = sis[i, 9, j] + (imp_rate/cluster_dis[k,j])*beta*sis[i-1, 7, k]/sis[i-1, 4 ,k]
        }
      }
      sis[i, 10, j] = 1/dur_inf
      
      # Probabilities
      sis[i, 11, j] = (1 - exp(-sis[i, 9, j]*time_step))
      sis[i, 12, j] = (1 - exp(-sis[i, 9, j]*(1 - vax_eff)*time_step))
      sis[i, 13, j] = (1 - exp(-sis[i, 10, j])*time_step)  
      
      # State variables
      sis[i, 14, j] = round(rbinom(n = 1, size = sis[i-1, 5, j], prob = sis[i, 11, j]), digits = 0)
      sis[i, 15, j] = round(rbinom(n = 1, size = sis[i-1, 6, j], prob = sis[i, 12, j]), digits = 0)
      sis[i, 16, j] = round(rbinom(n = 1, size = sis[i-1, 7, j], prob = sis[i, 13, j]), digits = 0)  
      
      # Model equations
      sis[i, 4, j] = sis[i, 5, j] + sis[i, 6, j] + sis[i, 7, j] + sis[i, 8, j]
      sis[i, 5, j] = sis[i-1, 5, j] - sis[i, 14, j] + round(sis[i-1, 8, j]*(1 - p_vax), digits = 0)
      sis[i, 6, j] = sis[i-1, 6, j] - sis[i, 15, j] + round(sis[i-1, 8, j]*p_vax, digits = 0)
      sis[i, 7, j] = sis[i-1, 7, j] + sis[i, 14, j] + sis[i, 15, j] - sis[i, 16, j]
      sis[i, 8, j] = sis[i-1, 8, j] + sis[i, 16, j] - sis[i-1, 8, j]
    }
  }
  
  
  # Store the results
  
  ## Infected
  
  sis_res_infected <- as.data.frame(matrix(0, nrow = length(time_seq), ncol = C))
  colnames(sis_res_infected) <- names_matrix
  for (i in 1:C) {
    sis_res_infected[,i] <- sis[,7,i]         
  }
  
  ## Detected infections
  
  sis_res_observed <- as.data.frame(matrix(0, nrow = length(time_seq), ncol = C))
  colnames(sis_res_observed) <- names_matrix
  for (i in 1:C) {
    sis_res_observed[,i] <- round(sis[,7,i]*mu, digits = 0)
  }
  
  ## Incidence of infection from S
  
  sis_res_inc_SI <- as.data.frame(matrix(0, nrow = length(time_seq), ncol = C))
  colnames(sis_res_inc_SI) <- names_matrix
  for (i in 1:C) {
    sis_res_inc_SI[,i] <- sis[,14,i]
  }
  
  ## Incidence of infection from V
  
  sis_res_inc_VI <- as.data.frame(matrix(0, nrow = length(time_seq), ncol = C))
  colnames(sis_res_inc_VI) <- names_matrix
  for (i in 1:C) {
    sis_res_inc_VI[,i] <- sis[,15,i]
  }
  
  ## Pivot the results and merge all together
  
  sis_res_infected <- sis_res_infected %>%
    mutate(time_seq = time_seq) %>%
    pivot_longer(-time_seq, names_to = "cluster", values_to = "infected") %>%
    mutate(cluster = readr::parse_number(cluster)) %>%
    arrange(cluster) %>%
    merge(y = cluster_data[, c("cluster", "vaccine")], by = "cluster", all = TRUE) %>%
    select(c("cluster", "vaccine", "time_seq", "infected"))
  
  sis_res_observed <- sis_res_observed %>%
    mutate(time_seq = time_seq) %>%
    pivot_longer(-time_seq, names_to = "cluster", values_to = "observed") %>%
    mutate(cluster = readr::parse_number(cluster)) %>%
    arrange(cluster) %>%
    merge(y = cluster_data[, c("cluster", "vaccine")], by = "cluster", all = TRUE) %>%
    select(c("cluster", "vaccine", "time_seq", "observed"))
  
  sis_res_inc_SI <- sis_res_inc_SI %>%
    mutate(time_seq = time_seq) %>%
    pivot_longer(-time_seq, names_to = "cluster", values_to = "inc_sus_inf") %>%
    mutate(cluster = readr::parse_number(cluster)) %>%
    arrange(cluster) %>%
    merge(y = cluster_data[, c("cluster", "vaccine")], by = "cluster", all = TRUE) %>%
    select(c("cluster", "vaccine", "time_seq", "inc_sus_inf"))
  
  sis_res_inc_VI <- sis_res_inc_VI %>%
    mutate(time_seq = time_seq) %>%
    pivot_longer(-time_seq, names_to = "cluster", values_to = "inc_vax_inf") %>%
    mutate(cluster = readr::parse_number(cluster)) %>%
    arrange(cluster) %>%
    merge(y = cluster_data[, c("cluster", "vaccine")], by = "cluster", all = TRUE) %>%
    select(c("cluster", "vaccine", "time_seq", "inc_vax_inf"))
  
  sis_result <- merge(sis_res_infected, sis_res_observed,
                      by = c("cluster", "vaccine", "time_seq"), all = TRUE)
  sis_result <- merge(sis_result, sis_res_inc_SI,
                      by = c("cluster", "vaccine", "time_seq"), all = TRUE)
  sis_result <- merge(sis_result, sis_res_inc_VI,
                      by = c("cluster", "vaccine", "time_seq"), all = TRUE)
  sis_result <- arrange(sis_result, cluster, time_seq)
  
  ## Returned object

  sis_result
}


# Function to run the SIS several times
# Append the results of each run to each other

sis_many <- function(..., n_runs) {
  
  # Generate output vector
  
  sis_result <- sis_model(...)
  sis_result <- sis_result %>%
    mutate(run = 1)
  sis_output <- sis_result
  
  for (i in 2:n_runs) {
    
    # Get result of one run
    
    sis_result <- sis_model(...)
    sis_result <- sis_result %>%
      mutate(run = i)
    
    # Merge to previous
    
    sis_output <- merge(sis_output, sis_result,
                        by = c("cluster", "vaccine", "time_seq",
                               "infected", "observed", "inc_sus_inf", "inc_vax_inf",
                               "run"), all = TRUE)
  }
  
  # Clean the result
  
  sis_output <- sis_output %>%
    select("run", "cluster", "vaccine", "time_seq",
           "infected", "observed", "inc_sus_inf", "inc_vax_inf") %>%
    arrange(run, cluster, time_seq)
  
  sis_output
}


#### Results functions ####


# Graph: mainly to check if run was okay

sis_graph <- function(sis_many_result) {
  
  ggplot() +
    geom_line(data = filter(sis_many_result, vaccine == 1),
              mapping = aes(x = time_seq, y = infected, group = run,
                            color = "Inf_Vax")) +
    geom_line(data = filter(sis_many_result, vaccine == 0),
              mapping = aes(x = time_seq, y = infected, group = run,
                            color = "Inf_No")) +
    geom_line(data = filter(sis_many_result, vaccine == 1),
              mapping = aes(x = time_seq, y = observed, group = run,
                            color = "Obs_Vax")) +
    geom_line(data = filter(sis_many_result, vaccine == 0),
              mapping = aes(x = time_seq, y = observed, group = run,
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
    xlim(c(1, 100)) +
    theme_classic() +
    labs(title = "Incidence over time (SIS)",
         x = "Time (days)",
         y = "Number of infections/detected infections") +
    theme(
      plot.title = element_text(size = rel(1.2), face="bold", hjust = 0.5),
      axis.title.x = element_text(size = rel(1.1), face="bold"),
      axis.title.y = element_text(size = rel(1.1), face="bold"),
      axis.text = element_text(size=rel(1)),
      legend.position = "bottom",
      legend.text = element_text(size=rel(1))) +
    facet_wrap(~cluster, ncol = 10, nrow = 10)
  
}


# Summary numbers

sis_summary <- function(sis_many_result, n_runs) {
  
  summary <- sis_many_result %>%
    group_by(run, cluster) %>%
    mutate(mean_inf = mean(infected)) %>%
    mutate(mean_obs = mean(observed)) %>%
    mutate(mean_SI = mean(inc_sus_inf)) %>%
    mutate(mean_VI = mean(inc_vax_inf)) %>%
    mutate(sum_inf = sum(infected)) %>%
    mutate(sum_obs = sum(observed)) %>%
    mutate(sum_SI = sum(inc_sus_inf)) %>%
    mutate(sum_VI = sum(inc_vax_inf)) %>%
    mutate(sum_all_inc = sum_SI + sum_VI) %>%
    filter(row_number() == 1) %>%
    select(-time_seq, -infected, -observed, -inc_sus_inf, -inc_vax_inf) %>%
    ungroup()

}


# Poisson regression : direct effects

sis_direct <- function(sis_summary_result, n_runs) {
  
  # Direct effect: vax vs unvax in vaccine cluster

  sis_summary_result <- sis_summary_result %>%
    
    # Only measured in vaccine clusters
    filter(vaccine == 1) %>%
    select(run, cluster, vaccine, sum_SI, sum_VI) %>%
    
    # Pivot to get the incidence SI vs VI only
    pivot_longer(-c(run, cluster, vaccine),
                 names_to = "vax_incluster", values_to = "incidence")
  
  # Output vector to store the results of the test
  
  output <- matrix(0, ncol = 6, nrow = n_runs)

  for (i in 1:n_runs) {
    
    model <- glm(formula = incidence ~ vax_incluster,
                 family = "poisson", data = filter(sis_summary_result, run == i))
    
    # Get the coefficients of the model
    
    x <- exp(summary(model)$coef)
    y <- exp(confint(model))
    
    # Store them
    
    output[i, 1:4] <- x[2,]
    output[i, 5:6] <- y[2,]
  }
  
  # Clean the result vector
  
  rownames(output) <- paste0("run_", seq(1:10))
  colnames(output) <- c("exp(Estimate)", "SE", "z Val", "Pr(>|z|)",
                        "2.5%", "97.5%")
  
  output
}


# Poisson regression : indirect effects

sis_indirect <- function(sis_summary_result, n_runs) {
  
  # Indirect effect: unvax in vaccine clusters vs non-vaccine clusters

  sis_summary_result <- sis_summary_result %>%
    
    # Focus therefore only on incidence S to I
    select(run, cluster, vaccine, sum_SI) 
  
  # Output vector to store the results of the test
  
  output <- matrix(0, ncol = 6, nrow = n_runs)
  
  for (i in 1:n_runs) {
    
    model <- glm(formula = sum_SI ~ vaccine,
                 family = "poisson", data = filter(sis_summary_result, run == i))
    
    # Get the coefficients of the model
    
    x <- exp(summary(model)$coef)
    y <- exp(confint(model))
    
    # Store them
    
    output[i, 1:4] <- x[2,]
    output[i, 5:6] <- y[2,]
  }
  
  # Clean the result vector
  
  rownames(output) <- paste0("run_", seq(1:10))
  colnames(output) <- c("exp(Estimate)", "SE", "z Val", "Pr(>|z|)",
                        "2.5%", "97.5%")
  
  output
}


# Poisson regression : total effects

sis_total <- function(sis_summary_result, n_runs) {
  
  # Total effect: vax in vaccine cluster vs unvax in non-vaccine
  
  sis_summary_result <- sis_summary_result %>%
    select(run, cluster, vaccine, sum_SI, sum_VI, sum_all_inc) %>%
    
    # Get one var with VI only from vaccine clusters,
    # and all incidence from non-vaccine
    mutate(incidence = case_when(vaccine == 1 ~ sum_VI,
                                 vaccine == 0 ~ sum_all_inc))
  
  # Output vector to store the results of the test
  
  output <- matrix(0, ncol = 6, nrow = n_runs)
  
  for (i in 1:n_runs) {
    
    model <- glm(formula = incidence ~ vaccine,
                 family = "poisson", data = filter(sis_summary_result, run == i))
    
    # Get the coefficients of the model
    
    x <- exp(summary(model)$coef)
    y <- exp(confint(model))
    
    # Store them
    
    output[i, 1:4] <- x[2,]
    output[i, 5:6] <- y[2,]
  }
  
  # Clean the result vector
  
  rownames(output) <- paste0("run_", seq(1:10))
  colnames(output) <- c("exp(Estimate)", "SE", "z Val", "Pr(>|z|)",
                        "2.5%", "97.5%")
  
  output
}


# Poisson regression : overall effects

sis_overall <- function(sis_summary_result, n_runs) {
  
  # Overall effect: all in vaccine cluster vs all in non-vaccine
  
  sis_summary_result <- sis_summary_result %>%
    
    # We are interested in overall incidence
    select(run, cluster, vaccine, sum_all_inc)

  # Output vector to store the results of the test
  
  output <- matrix(0, ncol = 6, nrow = n_runs)
  
  for (i in 1:n_runs) {
    
    model <- glm(formula = sum_all_inc ~ vaccine,
                 family = "poisson", data = filter(sis_summary_result, run == i))
    
    # Get the coefficients of the model
    
    x <- exp(summary(model)$coef)
    y <- exp(confint(model))
    
    # Store them
    
    output[i, 1:4] <- x[2,]
    output[i, 5:6] <- y[2,]
  }
  
  # Clean the result vector
  
  rownames(output) <- paste0("run_", seq(1:10))
  colnames(output) <- c("exp(Estimate)", "SE", "z Val", "Pr(>|z|)",
                        "2.5%", "97.5%")
  
  output
}



# CRT: SIMULATIONS --------------------------------------------------------


#### 1: Reach equilibrium ####

equilibrium <- sis_equi(N = N, C = C, cluster_n = cluster_n,
                        cluster_dis = cluster_dis,
                        number_infected = number_infected,
                        R0 = R0, dur_inf = dur_inf, imp_rate = imp_rate,
                        time_step = time_step, years = 5)


#### 2: Run several simulations ####

# Run!

test_sis <- sis_many(N = N, C = C, cluster_n = cluster_n, cluster_dis = cluster_dis,
                     R0 = R0, dur_inf = dur_inf, imp_rate = imp_rate,
                     p_vax = p_vax, p_clusvax = p_clusvax, vax_eff = vax_eff,
                     p_sym = p_sym, p_test = p_test, p_positive = p_positive,
                     time_step = time_step, years = 1,
                     equilibrium = equilibrium, n_runs = 10)

# Graph

plot_sis <- sis_graph(test_sis)

# Get summary numbers

summary_sis <- sis_summary(test_sis, n_runs = 10)

direct_effect <- sis_direct(summary_sis, n_runs = 10)
indirect_effect <- sis_indirect(summary_sis, n_runs = 10)
total_effect <- sis_total(summary_sis, n_runs = 10)
overall_effect <- sis_overall(summary_sis, n_runs = 10)
