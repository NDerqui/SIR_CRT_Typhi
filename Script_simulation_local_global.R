
## Cluster Randomised Trial for a Typhoid vaccine in Vellore ##
      ## MRes Biomedical Research - EECID - Project 2 ##


## SIR model to simulate TCV vaccine allocation to clusters ##

                  ## Model main function ##

## Run SIR open model without vaccination to reach equilibrium ##
  ## Then, allocate vaccine to some 1/2 clusters and follow  ##




# SET UP ------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(sandwich)
library(lmtest)
library(forcats)
library(DescTools)
library(svMisc)
library(wesanderson)
library(ggpubr)




# PARAMETERS --------------------------------------------------------------


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

death <- 0.007        # 7 per 1,000


#### Transmission ####

# (infections/time): beta = R0/Duration of infectiousness

R0 <- 2            # Basic reproduction number: to explore
dur_inf <- 7       # Duration of infectiousness (days)

# Percentage of local transmission (in the same cluster),
# the rest of global transmission (external clusters)

per_local <- 0.95


#### Vaccination ####

p_vax <- 0.5       # Proportion of vaccinated population in vaccine clusters
p_clusvax <- 0.5   # Proportion of clusters assigned to vaccine
vax_eff <- 0.7     # Direct vaccine efficacy (infection)


#### Detected infections ####

p_sym <- 0.50       # Probability of being symptomatic - CHECK
p_test <- 0.50      # Probability of seeking a test - CHECK
p_positive <- 0.60  # Probability of test being positive


#### Population ####

# Total population in the study and number of clusters

N <- 200000               # Total population in the study
C <- 100                  # Number of clusters
var <- 0.01               # Variability in cluster populations: 1% or 20%

random_cluster <- 1

# If 1 the first 50 clusters are vaccinated, but randomly distributed in map
# If 0 the vaccine clusters are allocated like a chessboard



# FUNCTIONS -------------------------------------------------------------


main <- function(N, C, var, random_cluster = 1, # Population and cluster characteristics
                 incidence, birth, death,       # Incidence, birth and death rate
                 R0, dur_inf, per_local,        # Infection parameters & % of local trans
                 p_vax, p_clusvax, vax_eff,     # Vaccination parameters
                 p_sym, p_test, p_positive,     # Observed infections
                 time_step = 1, years1, years2, # Time: equilibrium sim y, vaccine sim y
                 n_runs) {   
  
  
  ## Cluster data
  
  # Vector of clusters
  
  cluster_no <- seq(1:C)
  
  # Population in each cluster vector
  
  n <- as.integer(rnorm(n = C, mean = N/C, sd = var*N/C))
  cluster_n <- abs(n)
  
  # Cluster (and vaccine allocation) distribution
  
  if (random_cluster == 1) {
    
    # Cluster map
    # Random location of clusters
    
    if (sqrt(C) %% 1 == 0) {
      cluster_map <- matrix(sample(cluster_no), ncol = sqrt(C), nrow = sqrt(C),
                            dimnames = list(seq(1:sqrt(C)), seq(1:sqrt(C))))
    } else {
      cluster_map <- matrix(c(sample(cluster_no), rep(NA, times = (ceiling(sqrt(C))*(floor(sqrt(C))+1) - C))),
                            ncol = ceiling(sqrt(C)), nrow = (floor(sqrt(C)) + 1),
                            dimnames = list(seq(1:(floor(sqrt(C))+1)), seq(1:ceiling(sqrt(C)))))
    }
    
    # Vaccination status of clusters
    # Vaccination of the first half of clusters (randomly located)
    
    V <- as.integer(C*p_clusvax)  
    
    if (C %% 2 == 0) {
      cluster_vstatus <- c(rep(1, times = V), rep(0, times = V))   
    } else {
      cluster_vstatus <- c(rep(1, times = V+1), rep(0, times = V)) 
    }
    
  } else {
    
    # Cluster map
    # Location of clusters is not random, follows the order
    
    if (sqrt(C) %% 1 == 0) {
      cluster_map <- matrix(cluster_no, ncol = sqrt(C), nrow = sqrt(C),
                            dimnames = list(seq(1:sqrt(C)), seq(1:sqrt(C))))
    } else {
      cluster_map <- matrix(c(cluster_no, rep(NA, times = (ceiling(sqrt(C))*(floor(sqrt(C))+1) - C))),
                            ncol = ceiling(sqrt(C)), nrow = (floor(sqrt(C)) + 1),
                            dimnames = list(seq(1:(floor(sqrt(C))+1)), seq(1:ceiling(sqrt(C)))))
    }
    
    # Vaccination status of clusters
    # Vaccination like a chessboard
     
    if (nrow(cluster_map) %% 2 == 0) {
      cluster_vstatus <- rep(c(rep(c(0,1), times = nrow(cluster_map)/2),
                               rep(c(1,0), times = nrow(cluster_map)/2)),
                             times = ncol(cluster_map)/2)
      cluster_vstatus <- cluster_vstatus[1:length(cluster_no)]
    } else {
      if (C %% 2 == 0) {
        cluster_vstatus <- c(rep(c(0, 1), times = C/2)) # Flag for vax clusters
      } else {
        cluster_vstatus <- c(rep(c(0, 1), times = C/2), 0) # Flag for vax clusters
      }
    }
    
  }
  
  # Cluster distance matrix
  
  cluster_dis <- matrix(1, nrow = C, ncol = C,
                        dimnames = list(seq(1:C), seq(1:C)))
  
  for (i in 1:C) {                          # Each cluster distance
    for (j in 1:C) {                        # with each of the others
      
      for (k in 1:nrow(cluster_map)) {                # Look for location of first cluster
        for (l in 1:ncol(cluster_map)) {
          if (i == cluster_map[k, l] & !is.na(cluster_map[k, l])) {
            
            for (m in 1:nrow(cluster_map)) {          # Look for location of second cluster
              for (n in 1:ncol(cluster_map)) {
                if (j == cluster_map[m, n] & !is.na(cluster_map[m, n])) {
                  
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
  
  # Cluster data for reference
  
  cluster_data <- data.frame(cluster = cluster_no,
                             vaccine = cluster_vstatus,
                             pop = cluster_n,
                             cluster_dis)
  colnames(cluster_data) <- c("cluster", "vaccine", "pop",
                              paste0("dis_", cluster_no))
  
  # Clean
  
  rm(i, j, k, l, m, n, vertical, horizontal, distance)
  
  
  ## Considerations
  
  # Calculated parameters
  
  R0 <- rnorm(n = C, mean = R0, sd = 0.1*R0) # Varying R0
  beta <- R0/dur_inf                         # Transmission rate
  mu <- p_sym*p_test*p_positive              # Prob of detecting I
  
  
  
  ## For several simulations  
  
  sir_output <- data.frame(cluster = 0, vaccine = 0, time_seq = 0,
                           susceptible = 0, vaccinated = 0, infected = 0, observed = 0,
                           inc_sus_inf = 0, inc_vax_inf = 0, run = 0)
  
  
  for (n in 1:n_runs) {
    
    
    ## Step 1: reach equilibrium
    
    # Total time
    
    time_seq1 <- seq(from = 1, to = 365*years1, by = time_step)
    
    # Build the array
    
    names_row1 <- paste0("t_", time_seq1)
    
    names_column1 <- c("cluster", "vax_status", "time_seq", 
                       "no_N", "no_S", "no_I", "no_R",
                       "hazard_inf", "inc_S_event", "inc_SI", "inc_SD",
                       "inc_I_event2", "inc_IR", "inc_ID",
                       "inc_RD", "inc_NB")
    
    names_matrix1 <- paste0("cluster_", cluster_no)
    
    sir_first <- array(NA, dim = c(length(time_seq1), 16, C),
                       dimnames = list(names_row1, names_column1, names_matrix1))
    
    # Assign initial values
    
    sir_first[, 3, ] <- time_seq1
    
    for (i in 1:C) {
      sir_first[, 1, i] = cluster_no[i]
      sir_first[, 2, i] = cluster_vstatus[i]
      sir_first[1, 4, i] = cluster_n[i]
      sir_first[1, 6, i] = abs(as.integer(incidence*sir_first[1, 4, i]*dur_inf))# Initial I depends on incidence rate x duration infectiousness
      sir_first[1, 5, i] = sir_first[1, 4 ,i] - sir_first[1, 6, i]              # Initial S depends on N - I
      sir_first[1, 7, i] = 0
    }
    
    # Put the model to run
    
    for (i in 2:length(time_seq1)) {
      
      for (j in 1:C) {
        
        # Hazard of infection
        
        sir_first[i, 8, j] = beta[j]*(per_local*sir_first[i-1, 6, j]/sir_first[i-1, 4, j] + (1 - per_local)*sum(sir_first[i-1, 6, -j])/sum(sir_first[i-1, 4, -j]))
        
        # From S to event (I or death, allow competing hazards) (Prob of I: (1 - exp(-(haz_inf)*time_step)))
        
        sir_first[i, 9, j] = abs(as.integer(rbinom(n = 1, size = sir_first[i-1, 5, j], prob = (1 - exp(-(sir_first[i, 8, j] + death)*time_step)))))
        
        # From event to I (prob is hazard_I/hazard_I+hazard_death)
        
        sir_first[i, 10, j] = abs(as.integer(rbinom(n = 1, size = sir_first[i, 9, j], prob = sir_first[i, 8, j]/(sir_first[i, 8, j] + death))))
        
        # From event to D (prob is hazard_death/hazard_I+hazard_death)
        
        sir_first[i, 11, j] = abs(as.integer(rbinom(n = 1, size = sir_first[i, 9, j], prob = death/(sir_first[i, 8, j] + death))))
        
        # From I to event2 (R or death, allow competing hazards) (Prob of recov: (1 - exp(-(1/dur_inf)*time_step)))
        
        sir_first[i, 12, j] = abs(as.integer(rbinom(n = 1, size = sir_first[i-1, 6, j], prob = (1 - exp(-((1/dur_inf) + death)*time_step)))))
        
        # From event2 to R (prob is hazard_R/hazard_R+hazard_death)
        
        sir_first[i, 13, j] = abs(as.integer(rbinom(n = 1, size = sir_first[i, 12, j], prob = (1/dur_inf)/((1/dur_inf) + death))))
        
        # From event2 to D (prob is hazard_death/hazard_R+hazard_death)
        
        sir_first[i, 14, j] = abs(as.integer(rbinom(n = 1, size = sir_first[i, 12, j], prob = death/((1/dur_inf) + death))))
        
        # Deaths from R (Prob of death the same across: (1 - exp(-death*time_step)))
        
        sir_first[i, 15, j] = abs(as.integer(rbinom(n = 1, size = sir_first[i-1, 7, j], prob = (1 - exp(-death*time_step)))))
        
        # Births from N (Prob of birth the same across: (1 - exp(-birth*time_step)))
        
        sir_first[i, 16, j] = abs(as.integer(rbinom(n = 1, size = sir_first[i-1, 4, j], prob = (1 - exp(-birth*time_step))))) 
        
        # Model equations
        
        sir_first[i, 5, j] = abs(as.integer(sir_first[i-1, 5, j] - sir_first[i, 10, j] - sir_first[i, 11, j] + sir_first[i, 16, j]))
        
        sir_first[i, 6, j] = abs(as.integer(sir_first[i-1, 6, j] + sir_first[i, 10, j] - sir_first[i, 13, j] - sir_first[i, 14, j]))
        
        sir_first[i, 7, j] = abs(as.integer(sir_first[i-1, 7, j] + sir_first[i, 13, j] - sir_first[i, 15, j]))
        
        sir_first[i, 4, j] = sir_first[i, 5, j] + sir_first[i, 6, j] + sir_first[i, 7, j]
        
      }
    }
    
    #  Store the results: input the last row for the vaccine simulation
    
    equilibrium_result <- as.data.frame(matrix(0, nrow = C, ncol = 16))
    colnames(equilibrium_result) <- names_column1
    for (i in 1:C) {
      equilibrium_result[i, ] <- sir_first[length(time_seq1), , i]         
    }
    
    
    ## Step 2: vaccine introduction
    
    
    # Total time
    
    time_seq2 <- seq(from = 1, to = 365*years2, by = time_step)
    
    # Build the array
    
    names_row2 <- paste0("t_", time_seq2)
    
    names_column2 <- c("cluster", "vaccine", "time_seq",
                       "no_N", "no_S", "no_V", "no_I", "no_R",
                       "haz_inf", "inc_S_event", "inc_SI", "inc_SD",
                       "inc_V_event2", "inc_VI", "inc_VD",
                       "inc_I_event3", "inc_IR", "inc_ID",
                       "inc_RD", "inc_NB")
    
    names_matrix2 <- paste0("cluster_", cluster_no)
    
    sir <- array(NA, dim = c(length(time_seq2), 20, C),
                 dimnames = list(names_row2, names_column2, names_matrix2))
    
    # Assign initial values
    
    sir[, 3, ] <- time_seq2
    
    for (i in 1:C) {
      sir[, 1, i] = cluster_no[i]
      sir[, 2, i] = cluster_vstatus[i]
      sir[1, 4, i] = equilibrium_result[i, 4] # Total N from simulation
      sir[1, 7, i] = equilibrium_result[i, 6] # I from the simulation
      sir[1, 8, i] = equilibrium_result[i, 7] # R from the simulation
    }
    
    for (i in 1:C) {
      if (sir[1, 2, i] == 1) {                                                                   # In vaccinated clusters
        sir[1, 5, i] = equilibrium_result[i, 5] - abs(as.integer(equilibrium_result[i, 5]*p_vax)) # S=S*(1-COVERAGE)
        sir[1, 6, i] = abs(as.integer(equilibrium_result[i, 5]*p_vax))                            # V=S*COVERAGE
        
      } else {                                   # In non-vax clusters
        sir[1, 5, i] = equilibrium_result[i, 5]   # Susceptible are all S from simulation
        sir[1, 6, i] = as.integer(0)              # None vaccinated
      }
    }
    
    # Put the model to run
    
    for (i in 2:length(time_seq2)) {
      
      for (j in 1:C) {
        
        # Hazard of infection
        
        sir[i, 9, j] = beta[j]*(per_local*sir[i-1, 7, j]/sir[i-1, 4, j] + (1 - per_local)*sum(sir[i-1, 7, -j])/sum(sir[i-1, 4, -j])) 
        
        # From S to event (I or death, allow competing hazards) (Prob of I: (1 - exp(-(haz_inf)*time_step)))
        
        sir[i, 10, j] = abs(as.integer(rbinom(n = 1, size = sir[i-1, 5, j], prob = (1 - exp(-(sir[i, 9, j] + death)*time_step)))))
        
        # From event to I (prob is hazard_I/hazard_I+hazard_death)
        
        sir[i, 11, j] = abs(as.integer(rbinom(n = 1, size = sir[i, 10, j], prob = sir[i, 9, j]/(sir[i, 9, j] + death))))
        
        # From event to D (prob is hazard_death/hazard_I+hazard_death)
        
        sir[i, 12, j] = abs(as.integer(rbinom(n = 1, size = sir[i, 10, j], prob = death/(sir[i, 9, j] + death))))
        
        # From V to event2 (I or death, allow competing hazards) (Prob of I: (1 - exp(-(haz_inf*(1 - vax_eff))*time_step)))
        
        sir[i, 13, j] = abs(as.integer(rbinom(n = 1, size = sir[i-1, 6, j], prob = (1 - exp(-(sir[i, 9, j]*(1 - vax_eff) + death)*time_step)))))
        
        # From event to I (prob is hazard_I*(1 - vax_eff)/hazard_I*(1 - vax_eff)+hazard_death)
        
        sir[i, 14, j] = abs(as.integer(rbinom(n = 1, size = sir[i, 13, j], prob = sir[i, 9, j]*(1 - vax_eff)/(sir[i, 9, j]*(1 - vax_eff) + death))))
        
        # From event to D (prob is hazard_death/hazard_I*(1 - vax_eff)+hazard_death)
        
        sir[i, 15, j] = abs(as.integer(rbinom(n = 1, size = sir[i, 13, j], prob = death/(sir[i, 9, j]*(1 - vax_eff) + death))))
        
        # From I to event3 (R or death, allow competing hazards) (Prob of recov: (1 - exp(-(1/dur_inf)*time_step)))
        
        sir[i, 16, j] = abs(as.integer(rbinom(n = 1, size = sir[i-1, 7, j], prob = (1 - exp(-((1/dur_inf) + death)*time_step)))))
        
        # From event2 to R (prob is hazard_R/hazard_R+hazard_death)
        
        sir[i, 17, j] = abs(as.integer(rbinom(n = 1, size = sir[i, 16, j], prob = (1/dur_inf)/((1/dur_inf) + death))))
        
        # From event2 to D (prob is hazard_death/hazard_R+hazard_death)
        
        sir[i, 18, j] = abs(as.integer(rbinom(n = 1, size = sir[i, 16, j], prob = death/((1/dur_inf) + death))))
        
        # Deaths from R (Prob of death the same across: (1 - exp(-death*time_step)))
        
        sir[i, 19, j] = abs(as.integer(rbinom(n = 1, size = sir[i-1, 8, j], prob = (1 - exp(-death*time_step)))))
        
        # Births from N (Prob of birth the same across: (1 - exp(-birth*time_step)))
        
        sir[i, 20, j] = abs(as.integer(rbinom(n = 1, size = sir[i-1, 4, j], prob = (1 - exp(-birth*time_step))))) 
        
        # Model equations
        if (sir[i, 2, j] == 1) {   # In vaccine clusters, births are divided into S and V
          
          sir[i, 5, j] = abs(as.integer(sir[i-1, 5, j] - sir[i, 11, j] - sir[i, 12, j] + abs(as.integer(sir[i, 20, j]*(1-p_vax)))))
          
          sir[i, 6, j] = abs(as.integer(sir[i-1, 6, j] - sir[i, 14, j] - sir[i, 15, j] + (sir[i, 20, j] - abs(as.integer(sir[i, 20, j]*(1-p_vax))))))
        
        } else {                   # In non-vaccine clusters, all births go to S
          
          sir[i, 5, j] = abs(as.integer(sir[i-1, 5, j] - sir[i, 11, j] - sir[i, 12, j] + sir[i, 20, j]))
          
          sir[i, 6, j] = sir[i-1, 6, j]
        }
        
        sir[i, 7, j] = abs(as.integer(sir[i-1, 7, j] + sir[i, 11, j] + sir[i, 14, j] - sir[i, 17, j] - sir[i, 18, j]))
        
        sir[i, 8, j] = abs(as.integer(sir[i-1, 8, j] + sir[i, 17, j] - sir[i, 19, j]))
        
        sir[i, 4, j] = sir[i, 5, j] + sir[i, 6, j] + sir[i, 7, j] + sir[i, 8, j]
        
      }
    }
    
    # Store the results
    
    ## Total S in each cluster
    sir_res_susceptible <- as.data.frame(matrix(0, nrow = length(time_seq2), ncol = C))
    colnames(sir_res_susceptible) <- names_matrix2
    for (i in 1:C) {
      sir_res_susceptible[,i] <- sir[,5,i]         
    }
    
    ## Total V in each cluster
    sir_res_vaccinated <- as.data.frame(matrix(0, nrow = length(time_seq2), ncol = C))
    colnames(sir_res_vaccinated) <- names_matrix2
    for (i in 1:C) {
      sir_res_vaccinated[,i] <- sir[,6,i]         
    }
    
    ## Infected
    sir_res_infected <- as.data.frame(matrix(0, nrow = length(time_seq2), ncol = C))
    colnames(sir_res_infected) <- names_matrix2
    for (i in 1:C) {
      sir_res_infected[,i] <- sir[,7,i]         
    }
    
    ## Detected infections
    sir_res_observed <- as.data.frame(matrix(0, nrow = length(time_seq2), ncol = C))
    colnames(sir_res_observed) <- names_matrix2
    for (i in 1:C) {
      sir_res_observed[,i] <- round(rbinom(n = 1, size = sir[,7,i], prob = mu), digits = 0)
    }
    
    ## Detected incidence of infection from S
    sir_res_inc_SI <- as.data.frame(matrix(0, nrow = length(time_seq2), ncol = C))
    colnames(sir_res_inc_SI) <- names_matrix2
    for (i in 1:C) {
      sir_res_inc_SI[,i] <- round(rbinom(n = 1, size = sir[,11,i], prob = mu), digits = 0)
    }
    
    ## Detected incidence of infection from V
    sir_res_inc_VI <- as.data.frame(matrix(0, nrow = length(time_seq2), ncol = C))
    colnames(sir_res_inc_VI) <- names_matrix2
    for (i in 1:C) {
      sir_res_inc_VI[,i] <- round(rbinom(n = 1, size = sir[,13,i], prob = mu), digits = 0)
    }
    
    ## Pivot the results and merge all together
    sir_res_susceptible <- sir_res_susceptible %>%
      mutate(time_seq = time_seq2) %>%
      pivot_longer(-time_seq, names_to = "cluster", values_to = "susceptible") %>%
      mutate(cluster = readr::parse_number(cluster)) %>%
      arrange(cluster) %>%
      merge(y = cluster_data[, c("cluster", "vaccine")], by = "cluster", all = TRUE) %>%
      select(c("cluster", "vaccine", "time_seq", "susceptible"))
    sir_res_vaccinated <- sir_res_vaccinated %>%
      mutate(time_seq = time_seq2) %>%
      pivot_longer(-time_seq, names_to = "cluster", values_to = "vaccinated") %>%
      mutate(cluster = readr::parse_number(cluster)) %>%
      arrange(cluster) %>%
      merge(y = cluster_data[, c("cluster", "vaccine")], by = "cluster", all = TRUE) %>%
      select(c("cluster", "vaccine", "time_seq", "vaccinated"))
    sir_res_infected <- sir_res_infected %>%
      mutate(time_seq = time_seq2) %>%
      pivot_longer(-time_seq, names_to = "cluster", values_to = "infected") %>%
      mutate(cluster = readr::parse_number(cluster)) %>%
      arrange(cluster) %>%
      merge(y = cluster_data[, c("cluster", "vaccine")], by = "cluster", all = TRUE) %>%
      select(c("cluster", "vaccine", "time_seq", "infected"))
    sir_res_observed <- sir_res_observed %>%
      mutate(time_seq = time_seq2) %>%
      pivot_longer(-time_seq, names_to = "cluster", values_to = "observed") %>%
      mutate(cluster = readr::parse_number(cluster)) %>%
      arrange(cluster) %>%
      merge(y = cluster_data[, c("cluster", "vaccine")], by = "cluster", all = TRUE) %>%
      select(c("cluster", "vaccine", "time_seq", "observed"))
    sir_res_inc_SI <- sir_res_inc_SI %>%
      mutate(time_seq = time_seq2) %>%
      pivot_longer(-time_seq, names_to = "cluster", values_to = "inc_sus_inf") %>%
      mutate(cluster = readr::parse_number(cluster)) %>%
      arrange(cluster) %>%
      merge(y = cluster_data[, c("cluster", "vaccine")], by = "cluster", all = TRUE) %>%
      select(c("cluster", "vaccine", "time_seq", "inc_sus_inf"))
    sir_res_inc_VI <- sir_res_inc_VI %>%
      mutate(time_seq = time_seq2) %>%
      pivot_longer(-time_seq, names_to = "cluster", values_to = "inc_vax_inf") %>%
      mutate(cluster = readr::parse_number(cluster)) %>%
      arrange(cluster) %>%
      merge(y = cluster_data[, c("cluster", "vaccine")], by = "cluster", all = TRUE) %>%
      select(c("cluster", "vaccine", "time_seq", "inc_vax_inf"))
    
    sir_result <- merge(sir_res_susceptible, sir_res_vaccinated,
                        by = c("cluster", "vaccine", "time_seq"), all = TRUE)
    sir_result <- merge(sir_result, sir_res_infected,
                        by = c("cluster", "vaccine", "time_seq"), all = TRUE)
    sir_result <- merge(sir_result, sir_res_observed,
                        by = c("cluster", "vaccine", "time_seq"), all = TRUE)
    sir_result <- merge(sir_result, sir_res_inc_SI,
                        by = c("cluster", "vaccine", "time_seq"), all = TRUE)
    sir_result <- merge(sir_result, sir_res_inc_VI,
                        by = c("cluster", "vaccine", "time_seq"), all = TRUE)
    sir_result <- arrange(sir_result, cluster, time_seq)
    
    
    ## Combine the results of several runs
    
    sir_result <- sir_result %>%
      mutate(run = n)
    
    sir_output <- merge(sir_output, sir_result,
                        by = c("cluster", "vaccine", "time_seq", "susceptible", "vaccinated",
                               "infected", "observed", "inc_sus_inf", "inc_vax_inf", "run"),
                        all = TRUE)
    
    sir_output <- sir_output %>%
      select("run", "cluster", "vaccine", "time_seq", "susceptible", "vaccinated",
             "infected", "observed", "inc_sus_inf", "inc_vax_inf") %>%
      arrange(run, cluster, time_seq)
    
    
    ## Progress
    
    progress(value = n, max.value = n_runs,
             progress.bar = TRUE)
    Sys.sleep(0.01)
    
    if (n == n_runs) {cat("Simulations done!\n")}
    
  }
  
  # Clean
  
  rm(i, j, n)
  
  sir_output <- sir_output[2:nrow(sir_output),]
  
  
  
  ## Results from all the runs
  
  # Summarize the results
  
  ## Total incidence of infection in the 2 years
  incidence_plot <- sir_output %>%
    # Eliminate simulations (run-cluster level) in which there was a problem
    mutate(eliminate = case_when(is.na(infected) ~ 1,
                                 infected >= 0 ~ 0)) %>%
    group_by(run, cluster) %>%
    mutate(eliminate = max(eliminate)) %>%
    filter(eliminate == 0) %>%
    ungroup() %>%
    # Sum all the incident infections for each vax arm by run
    group_by(run, vaccine) %>%
    mutate(sum_SI = sum(inc_sus_inf, na.rm = TRUE)) %>%
    mutate(sum_VI = sum(inc_vax_inf, na.rm = TRUE)) %>%
    mutate(total_obs_inc_by_vax = sum_SI + sum_VI) %>%
    # Get only the observed
    mutate(total_obs_inc_by_vax = total_obs_inc_by_vax*mu) %>%
    ungroup() %>%
    # Sum all the incident infections for both arms by run
    group_by(run) %>%
    mutate(sum_SI = sum(inc_sus_inf, na.rm = TRUE)) %>%
    mutate(sum_VI = sum(inc_vax_inf, na.rm = TRUE)) %>%
    mutate(total_obs_inc_all = sum_SI + sum_VI) %>%
    # Get only the observed
    mutate(total_obs_inc_all = total_obs_inc_all*mu) %>%
    ungroup() %>%
    # Delete unnecessary cols and rows
    group_by(run, vaccine) %>%
    filter(row_number() == 1) %>%
    select(run, vaccine, total_obs_inc_by_vax, total_obs_inc_all) %>%
    ungroup() %>%
    # Incidence per year
    mutate(total_obs_inc_by_vax = total_obs_inc_by_vax/years2) %>%
    mutate(total_obs_inc_all = total_obs_inc_all/years2) %>%
    # Pivot the results
    pivot_wider(id_cols = c(run, total_obs_inc_all),
                names_from = vaccine, values_from = total_obs_inc_by_vax)
  colnames(incidence_plot) <- c("run", "Total_inc_year", "Vax_inc_year", "NoVax_inc_year")
  # Incidence per 1,000 people
  incidence_plot <- incidence_plot %>%
    mutate(Total_inc_year_1000 = Total_inc_year*1000/N) %>%
    mutate(Vax_inc_year_1000 = Vax_inc_year*1000/N) %>%
    mutate(NoVax_inc_year_1000 = NoVax_inc_year*1000/N) 
  
  summary <- incidence_plot %>%
    mutate(Total_inc_year = mean(Total_inc_year, na.rm = TRUE)) %>%
    mutate(Total_inc_year_1000 = mean(Total_inc_year_1000, na.rm = TRUE)) %>%
    mutate(Vax_inc_year = mean(Vax_inc_year, na.rm = TRUE)) %>%
    mutate(Vax_inc_year_1000 = mean(Vax_inc_year_1000, na.rm = TRUE)) %>%
    mutate(NoVax_inc_year = mean(NoVax_inc_year, na.rm = TRUE)) %>%
    mutate(NoVax_inc_year_1000 = mean(NoVax_inc_year_1000, na.rm = TRUE)) %>%
    filter(row_number() == 1) %>%
    select(-run)
  
  ## Data for the Poisson regression  
  for_poisson <- sir_output %>%
    # Eliminate simulations (run-cluster level) in which there was a problem
    mutate(eliminate = case_when(is.na(infected) ~ 1,
                                 infected >= 0 ~ 0)) %>%
    group_by(run, cluster) %>%
    mutate(eliminate = max(eliminate)) %>%
    filter(eliminate == 0) %>%
    ungroup() %>%
    # Aim: : estimate incidence per person-year
    # Calculate sum of all susceptible (or vaccinated) in year
    group_by(run, vaccine) %>%
    mutate(sus_risk = susceptible) %>%
    mutate(vax_risk = vaccinated) %>%
    mutate(sus_risk = sum(sus_risk, na.rm = TRUE)*time_step/365) %>%
    mutate(vax_risk = sum(vax_risk, na.rm = TRUE)*time_step/365) %>%
    # Sum all incidence of infection: per vax status and total
    mutate(sum_SI = sum(inc_sus_inf, na.rm = TRUE)) %>%
    mutate(sum_VI = sum(inc_vax_inf, na.rm = TRUE)) %>%
    mutate(sum_all_inc = sum_SI + sum_VI) %>%
    # Get only the observed
    mutate(sum_SI = rbinom(n = 1, size = sum_SI, prob = mu)) %>%
    mutate(sum_VI = rbinom(n = 1, size = sum_VI, prob = mu)) %>%
    mutate(sum_all_inc = rbinom(n = 1, size = sum_all_inc, prob = mu)) %>%
    # Only one obs per run and cluster
    filter(row_number() == 1) %>%
    select(-time_seq, -cluster, -eliminate, -susceptible, -vaccinated,
           -infected, -observed, -inc_sus_inf, -inc_vax_inf) %>%
    ungroup()
  
  ## Calculating ICC and DesignEffect
  icc <- sir_output %>%
    mutate(sum_all_inc = inc_sus_inf + inc_vax_inf, na.rm = TRUE) %>%
    select(run, cluster, sum_all_inc) %>%
    group_by(run, cluster) %>%
    mutate(mean_clus = mean(sum_all_inc, na.rm = TRUE)) %>%
    mutate(group_mean_centered = sum_all_inc - mean_clus) %>%
    ungroup() %>%
    group_by(run) %>%
    mutate(within_var = var(group_mean_centered, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(run, cluster) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    group_by(run) %>%
    mutate(between_var = var(mean_clus, na.rm = TRUE)) %>%
    filter(row_number() == 1) %>%
    select(run, within_var, between_var) %>%
    mutate(icc = between_var/(between_var + within_var)) %>%
    mutate(des_eff = 1 + (mean(cluster_n) - 1)*icc) %>%
    ungroup() %>%
    as.data.frame()
  
  # Plot
  
  plot_a_data <- incidence_plot %>%
    select(run, Vax_inc_year, NoVax_inc_year, Total_inc_year) %>%
    pivot_longer(-run)
  plot_a <- ggplot(data = plot_a_data) +
    geom_boxplot(mapping = aes(x = factor(name, level = c("Vax_inc_year", "NoVax_inc_year", "Total_inc_year")),
                               y = value, color = name, fill = name), alpha = 0.2) +
    geom_jitter(mapping = aes(x = factor(name, level = c("Vax_inc_year", "NoVax_inc_year", "Total_inc_year")),
                              y = value, color = name, fill = name), size = 1.5) +
    scale_color_manual(name = NULL,
                       breaks = c("Vax_inc_year", "NoVax_inc_year", "Total_inc_year"),
                       values = c(wes_palette("GrandBudapest1", n = 3)),
                       labels = c("Incidence in vaccine clusters",
                                  "Incidence in non-vaccine clusters",
                                  "Incidence in all clusters")) +
    scale_fill_manual(name = NULL,
                       breaks = c("Vax_inc_year", "NoVax_inc_year", "Total_inc_year"),
                       values = c(wes_palette("GrandBudapest1", n = 3)),
                       labels = c("Incidence in vaccine clusters",
                                  "Incidence in non-vaccine clusters",
                                  "Incidence in all clusters")) +
    scale_x_discrete(name = NULL,
                     breaks = c("Vax_inc_year", "NoVax_inc_year", "Total_inc_year"),
                     labels = c("Vaccine clusters", "Non-vaccine clusters", "All clusters")) +
    scale_y_continuous(name = "Incidence per year") +
    labs(title = "A") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = rel(1.2), face="bold"),
      axis.title.x = element_text(size = rel(1.1), face="bold"),
      axis.title.y = element_text(size = rel(1.1), face="bold"),
      axis.text = element_text(size=rel(1), face="bold"),
      legend.position = "right",
      legend.text = element_text(size=rel(1)))
  plot_b_data <- incidence_plot %>%
    select(run, Vax_inc_year_1000, NoVax_inc_year_1000, Total_inc_year_1000) %>%
    pivot_longer(-run)
  plot_b <- ggplot(data = plot_b_data) +
    geom_boxplot(mapping = aes(x = factor(name, level = c("Vax_inc_year_1000", "NoVax_inc_year_1000", "Total_inc_year_1000")),
                               y = value, color = name, fill = name), alpha = 0.2) +
    geom_jitter(mapping = aes(x = factor(name, level = c("Vax_inc_year_1000", "NoVax_inc_year_1000", "Total_inc_year_1000")),
                              y = value, color = name, fill = name), size = 1.5) +
    scale_color_manual(name = NULL,
                       breaks = c("Vax_inc_year_1000", "NoVax_inc_year_1000", "Total_inc_year_1000"),
                       values = c(wes_palette("GrandBudapest1", n = 3)),
                       labels = c("Incidence in vaccine clusters",
                                  "Incidence in non-vaccine clusters",
                                  "Incidence in all clusters")) +
    scale_fill_manual(name = NULL,
                      breaks = c("Vax_inc_year_1000", "NoVax_inc_year_1000", "Total_inc_year_1000"),
                      values = c(wes_palette("GrandBudapest1", n = 3)),
                      labels = c("Incidence in vaccine clusters",
                                 "Incidence in non-vaccine clusters",
                                 "Incidence in all clusters")) +
    scale_x_discrete(name = NULL,
                     breaks = c("Vax_inc_year_1000", "NoVax_inc_year_1000", "Total_inc_year_1000"),
                     labels = c("Vaccine clusters", "Non-vaccine clusters", "All clusters")) +
    scale_y_continuous(name = "Incidence per 1,000 & per year") +
    labs(title = "B") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = rel(1.2), face="bold"),
      axis.title.x = element_text(size = rel(1.1), face="bold"),
      axis.title.y = element_text(size = rel(1.1), face="bold"),
      axis.text = element_text(size=rel(1), face="bold"),
      legend.position = "right",
      legend.text = element_text(size=rel(1)))
  incidence_plot <- ggarrange(plot_a, plot_b,
            common.legend = TRUE, legend = "bottom")
  
  
  
  ## Poisson regressions
  
  
  # Direct effect: vax vs unvax in vaccine cluster
  
  sir_direct <- for_poisson %>%
    # Only measured in vaccine clusters
    filter(vaccine == 1) %>%
    select(run, vaccine, sus_risk, vax_risk, sum_SI, sum_VI) %>%
    # Pivot to get the incidence SI vs VI only
    pivot_longer(-c(run, vaccine, sus_risk, vax_risk),
                 names_to = "vax_incluster", values_to = "incidence") %>%
    # Change names
    mutate(vax_incluster = case_when(vax_incluster == "sum_VI" ~ 1,
                                     vax_incluster == "sum_SI" ~ 0)) %>%
    # Get the total
    mutate(total = case_when(vax_incluster == 1 ~ vax_risk,
                             vax_incluster == 0 ~ sus_risk)) %>%
    # Clean
    select(run, vax_incluster, incidence, total)
  
  # Output vector to store the results of the test
  
  output_direct <- matrix(0, ncol = 6, nrow = n_runs)
  
  for (i in 1:n_runs) {
    
    model <- glm(formula = as.formula((incidence/total) ~ vax_incluster),
                 family = "poisson"(link = "log"),
                 data = filter(sir_direct, run == i))
    
    # Do the robust standard errors
    
    robust <- coeftest(model, vcov. = sandwich)
    
    # Get the coefficients of the model
    
    x <- exp(robust[2, 1])
    y <- robust[2, 2:4]
    z <- exp(confint(robust))
    
    # Store them
    
    output_direct[i, 1]   <- x
    output_direct[i, 2:4] <- y
    output_direct[i, 5:6] <- z[2,]
  }
  
  # Clean the result vector
  
  rownames(output_direct) <- paste0("run_", seq(1:n_runs))
  colnames(output_direct) <- c("exp(Estimate)", "SE", "z Val", "Pr(>|z|)",
                              "2.5%", "97.5%")
  

  # Indirect effect: unvax in vaccine clusters vs non-vaccine clusters
    
  sir_indirect <- for_poisson %>%
    # Focus therefore only on incidence S to I
    select(run, vaccine, sum_SI, sus_risk)
    
  # Output vector to store the results of the test
    
  output_indirect <- matrix(0, ncol = 6, nrow = n_runs)
    
  for (i in 1:n_runs) {
      
  model <- glm(formula = as.formula((sum_SI/sus_risk) ~ vaccine),
               family = "poisson"(link = "log"),
               data = filter(sir_indirect, run == i))
      
    # Do the robust standard errors
      
    robust <- robust <- coeftest(model, vcov. = sandwich)
      
    # Get the coefficients of the model
      
    x <- exp(robust[2, 1])
    y <- robust[2, 2:4]
    z <- exp(confint(robust))
      
    # Store them
      
    output_indirect[i, 1]   <- x
    output_indirect[i, 2:4] <- y
    output_indirect[i, 5:6] <- z[2,]
  }
    
  # Clean the result vector
    
  rownames(output_indirect) <- paste0("run_", seq(1:n_runs))
  colnames(output_indirect) <- c("exp(Estimate)", "SE", "z Val", "Pr(>|z|)",
                          "2.5%", "97.5%")

  
  # Total effect: vax in vaccine cluster vs unvax in non-vaccine
    
  sir_total <- for_poisson %>%
    # Select
    select(run, vaccine, sus_risk, vax_risk, sum_VI, sum_all_inc) %>%
    # Get one var with VI only from vaccine clusters,
    # and all incidence from non-vaccine
    mutate(incidence = case_when(vaccine == 1 ~ sum_VI,
                                 vaccine == 0 ~ sum_all_inc)) %>%
    # Get total
    mutate(total = case_when(vaccine == 1 ~ vax_risk,
                             vaccine == 0 ~ sus_risk)) %>%
    # Clean
    select(run, vaccine, incidence, total)
    
  # Output vector to store the results of the test
    
  output_total <- matrix(0, ncol = 6, nrow = n_runs)
    
  for (i in 1:n_runs) {
      
    model <- glm(formula = as.formula((incidence/total) ~ vaccine),
                 family = "poisson"(link = "log"),
                 data = filter(sir_total, run == i))
      
    # Do the robust standard errors
      
    robust <- robust <- coeftest(model, vcov. = sandwich)
      
    # Get the coefficients of the model
      
    x <- exp(robust[2, 1])
    y <- robust[2, 2:4]
    z <- exp(confint(robust))
      
    # Store them
      
    output_total[i, 1]   <- x
    output_total[i, 2:4] <- y
    output_total[i, 5:6] <- z[2,]
  }
    
  # Clean the result vector
    
  rownames(output_total) <- paste0("run_", seq(1:n_runs))
  colnames(output_total) <- c("exp(Estimate)", "SE", "z Val", "Pr(>|z|)",
                                "2.5%", "97.5%")
     
  
 # Overall effect: all in vaccine cluster vs all in non-vaccine
    
 sir_overall <- for_poisson %>%
   # We are interested in overall incidence
   select(run, vaccine, sus_risk, vax_risk, sum_all_inc) %>%
   # Get total
   mutate(total = sus_risk + vax_risk) %>%
   # Clean
   select(run, vaccine, sum_all_inc, total)
    
  # Output vector to store the results of the test
    
  output_overall <- matrix(0, ncol = 6, nrow = n_runs)
    
    for (i in 1:n_runs) {
      
      model <- glm(formula = as.formula((sum_all_inc/total) ~ vaccine),
                   family = "poisson"(link = "log"), 
                   data = filter(sir_overall, run == i))
      
      # Do the robust standard errors
      
      robust <- robust <- coeftest(model, vcov. = sandwich)
      
      # Get the coefficients of the model
      
      x <- exp(robust[2, 1])
      y <- robust[2, 2:4]
      z <- exp(confint(robust))
      
      # Store them
      
      output_overall[i, 1]   <- x
      output_overall[i, 2:4] <- y
      output_overall[i, 5:6] <- z[2,]
   }
    
  # Clean the result vector
    
  rownames(output_overall) <- paste0("run_", seq(1:n_runs))
  colnames(output_overall) <- c("exp(Estimate)", "SE", "z Val", "Pr(>|z|)",
                          "2.5%", "97.5%")
    

  # Summary of all effects
  
  poisson_summary <- matrix(0, nrow = 4, ncol = 4)
  
  # Fill with the RR estimates
  poisson_summary[1, 1] <- mean(output_direct[,1], na.rm = TRUE)
  poisson_summary[1, 2] <- quantile(output_direct[,1], probs = c(0.025, 0.975), na.rm = TRUE)[1]
  poisson_summary[1, 3] <- quantile(output_direct[,1], probs = c(0.025, 0.975), na.rm = TRUE)[2]
  poisson_summary[2, 1] <- mean(output_indirect[,1], na.rm = TRUE)
  poisson_summary[2, 2] <- quantile(output_indirect[,1], probs = c(0.025, 0.975), na.rm = TRUE)[1]
  poisson_summary[2, 3] <- quantile(output_indirect[,1], probs = c(0.025, 0.975), na.rm = TRUE)[2]
  poisson_summary[3, 1] <- mean(output_total[,1], na.rm = TRUE)
  poisson_summary[3, 2] <- quantile(output_total[,1], probs = c(0.025, 0.975), na.rm = TRUE)[1]
  poisson_summary[3, 3] <- quantile(output_total[,1], probs = c(0.025, 0.975), na.rm = TRUE)[2]
  poisson_summary[4, 1] <- mean(output_overall[,1], na.rm = TRUE)
  poisson_summary[4, 2] <- quantile(output_overall[,1], probs = c(0.025, 0.975), na.rm = TRUE)[1]
  poisson_summary[4, 3] <- quantile(output_overall[,1], probs = c(0.025, 0.975), na.rm = TRUE)[2]
  
  # Get the power: prop of sims in which the CI do not cross one
  
  counts_dir <- n_runs
  for (i in 1:n_runs) {
    if(output_direct[i,5] < 1 & output_direct[i,6] >= 1) {
      counts_dir <- counts_dir - 1
    }
  }
  counts_ind <- n_runs
  for (i in 1:n_runs) {
    if(output_indirect[i,5] < 1 & output_indirect[i,6] >= 1) {
      counts_ind <- counts_ind - 1
    }
  }
  counts_tot <- n_runs
  for (i in 1:n_runs) {
    if(output_total[i,5] < 1 & output_total[i,6] >= 1) {
      counts_tot <- counts_tot - 1
    }
  }
  counts_ove <- n_runs
  for (i in 1:n_runs) {
    if(output_overall[i,5] < 1 & output_overall[i,6] >= 1) {
      counts_ove <- counts_ove - 1
    }
  }
  
  poisson_summary[1, 4] <- counts_dir/n_runs*100
  poisson_summary[2, 4] <- counts_ind/n_runs*100
  poisson_summary[3, 4] <- counts_tot/n_runs*100
  poisson_summary[4, 4] <- counts_ove/n_runs*100
  
  rownames(poisson_summary) <- c("Direct", "Indirect", "Total", "Overall")
  colnames(poisson_summary) <- c("Mean Effect", "Lower CI from Mean", "Upper CI from Mean", "Power %")
  
  poisson_summary <- round(poisson_summary, digits = 4)
  
  
  
  ## ICC-DE summary
  
  other_summary <- matrix(0, nrow = 3, ncol = 4)
  
  other_summary[1, 1] <- mean(icc[,4], na.rm = TRUE)
  other_summary[1, 2] <- quantile(icc[,4], probs = c(0.025, 0.975), na.rm = TRUE)[1]
  other_summary[1, 3] <- quantile(icc[,4], probs = c(0.025, 0.975), na.rm = TRUE)[2]
  other_summary[2, 1] <- mean(icc[,5], na.rm = TRUE)
  other_summary[2, 2] <- quantile(icc[,5], probs = c(0.025, 0.975), na.rm = TRUE)[1]
  other_summary[2, 3] <- quantile(icc[,5], probs = c(0.025, 0.975), na.rm = TRUE)[2]
  other_summary[3, 1] <- mean(R0, na.rm = TRUE)
  other_summary[3, 2] <- quantile(R0, probs = c(0.025, 0.975), na.rm = TRUE)[1]
  other_summary[3, 3] <- quantile(R0, probs = c(0.025, 0.975), na.rm = TRUE)[2]
  
  rownames(other_summary) <- c("ICC", "DEsign Effect", "R0")
  colnames(other_summary) <- c("Mean Effect", "Lower CI from Mean", "Upper CI from Mean", "Power %")
  
  
  
  ## Returned objects
  
  name_simulation <- paste0("C=", C, " var=", var, " PerLocal=", per_local, " N=", N,
                            " VE=", vax_eff, " Cover=", p_vax, " VaxArm=", p_clusvax)
  
  other_characteristics <- paste0("Incidence=", incidence, " BirthRate=", birth, " DeathRate=", death,
                                  " Mean_R0=", mean(R0), " DurationInf=", dur_inf, 
                                  " ProbSymp=", p_sym, " ProbTest=", p_test, " ProbPos=", p_positive,
                                  " YearsEq=", years1, " YearsVax=", years2)
  
  list <- list(name_simulation = name_simulation,
               other_characteristics = other_characteristics,
               cluster_pop = cluster_n,
               incidence_plot = incidence_plot,
               summary_infections = summary,
               output_direct = output_direct,
               output_indirect = output_indirect,
               output_overall = output_overall,
               output_total = output_total,
               poisson_summary = poisson_summary,
               icc_design = icc,
               other_summary = other_summary,
               reference_data = list(cluster_map = cluster_map,
                                     cluster_data = cluster_data,
                                     equilibrium_result = equilibrium_result,
                                     sir_output = sir_output)
  )

  list
}




# SIMULATIONS -------------------------------------------------------------


#### Run simulations ####

# With these characteristics

N
C
var
vax_eff
p_vax
p_clusvax

# Run

run <- main(N = 200000, C = 10, var = 0.2,
            incidence = incidence, birth = birth, death = death,
            R0 = R0, dur_inf = dur_inf, per_local = 0.5,
            p_vax = 0.5, p_clusvax = p_clusvax, vax_eff = 0.8,
            p_sym = p_sym, p_test = p_test, p_positive = p_positive,
            years1 = 3, years2 = 2, n_runs = 10)


#### Save results #### 

library(here)
library(openxlsx)

dir.create(here("Results/", Sys.Date()),recursive = TRUE)

# Give name to simulation

name <- run[[1]]

dir.create(here(paste0("Results/", Sys.Date(), "/", name)),recursive = TRUE)

# Save

write.table(run[[2]], file = paste0("Results/", Sys.Date(), "/", name,"/characteristics.txt"))

png(paste0("Results/", Sys.Date(), "/", name,"/Cluster_Pop_Hist.png"),
    width = 9, height = 5, units = 'in', res = 600)
hist(run[[3]], main = "Histogram of clusters' population",
     xlab = "Clusters' population", ylab = "Frequency of clusters")
dev.off()

png(paste0("Results/", Sys.Date(), "/", name,"/Incidence_plot.png"),
    width = 9.2, height = 4, units = 'in', res = 600)
run[[4]]
dev.off()

write.xlsx(as.data.frame(run[[5]]), rowNames = TRUE,
           paste0("Results/", Sys.Date(), "/", name,"/Total_Incidence.xlsx"))

write.xlsx(as.data.frame(run[[6]]), rowNames = TRUE,
           paste0("Results/", Sys.Date(), "/", name,"/RR_Direct_Effect.xlsx"))
write.xlsx(as.data.frame(run[[7]]), rowNames = TRUE,
           paste0("Results/", Sys.Date(), "/", name,"/RR_Indirect_Effect.xlsx"))
write.xlsx(as.data.frame(run[[8]]), rowNames = TRUE,
           paste0("Results/", Sys.Date(), "/", name,"/RR_Overall_Effect.xlsx"))
write.xlsx(as.data.frame(run[[9]]), rowNames = TRUE,
           paste0("Results/", Sys.Date(), "/", name,"/RR_Total_Effect.xlsx"))
write.xlsx(as.data.frame(run[[10]]), rowNames = TRUE,
           paste0("Results/", Sys.Date(), "/", name,"/Summary_RR_Effects.xlsx"))

write.xlsx(as.data.frame(run[[11]]), rowNames = TRUE,
           paste0("Results/", Sys.Date(), "/", name,"/IntraCC_DesignEffect.xlsx"))
write.xlsx(as.data.frame(run[[12]]), rowNames = TRUE,
           paste0("Results/", Sys.Date(), "/", name,"/Other_sum.xlsx"))

saveRDS(run[[13]], file = paste0("Results/", Sys.Date(), "/", name,"/All_data_reference.Rds"))


#### Do in a loop ####

# These are the conditions in which I am exploring R0:
# 3 y to eq; 2 y with vax; 30 runs

# R0 general analyses to be continued with 1.7, 1.5, 1.2

pop_list <- c(200000, 100000, 50000, 20000, 10000)
var_list <- c(0.01, 0.2)
cover_list <- c(0.5, 0.7, 0.9)

for (i in 1:length(pop_list)) {
  
  for (j in 1:length(var_list)) {
    
    for (k in 1:length(cover_list)) {
      
      run <- main(N = pop_list[i], C = 80, var = var_list[j],
                  incidence = incidence, birth = birth, death = death,
                  R0 = 2, dur_inf = dur_inf, per_local = 0.5,
                  p_vax = cover_list[k], p_clusvax = p_clusvax, vax_eff = 0.7,
                  p_sym = p_sym, p_test = p_test, p_positive = p_positive,
                  years1 = 3, years2 = 2, n_runs = 30)
      
      # Directory
      
      dir.create(here("Results/", Sys.Date()),recursive = TRUE)
      
      # Give name to simulation
      
      name <- run[[1]]
      
      dir.create(here(paste0("Results/", Sys.Date(), "/", name)),recursive = TRUE)
      
      # Save
      
      write.table(run[[2]], file = paste0("Results/", Sys.Date(), "/", name,"/characteristics.txt"))
      
      png(paste0("Results/", Sys.Date(), "/", name,"/Cluster_Pop_Hist.png"),
          width = 9, height = 5, units = 'in', res = 600)
      hist(run[[3]], main = "Histogram of clusters' population",
           xlab = "Clusters' population", ylab = "Frequency of clusters")
      dev.off()
      
      png(paste0("Results/", Sys.Date(), "/", name,"/Infections_after_vax.png"),
          width = 14, height = 9, units = 'in', res = 600)
      print(run[[4]])
      dev.off()
      
      write.xlsx(as.data.frame(run[[5]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/Total_Infections.xlsx"))
      
      write.xlsx(as.data.frame(run[[6]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/RR_Direct_Effect.xlsx"))
      write.xlsx(as.data.frame(run[[7]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/RR_Indirect_Effect.xlsx"))
      write.xlsx(as.data.frame(run[[8]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/RR_Overall_Effect.xlsx"))
      write.xlsx(as.data.frame(run[[9]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/RR_Total_Effect.xlsx"))
      write.xlsx(as.data.frame(run[[10]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/Summary_RR_Effects.xlsx"))
      
      write.xlsx(as.data.frame(run[[11]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/IntraCC_DesignEffect.xlsx"))
      write.xlsx(as.data.frame(run[[12]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/Other_sum.xlsx"))
      
      saveRDS(run[[13]], file = paste0("Results/", Sys.Date(), "/", name,"/All_data_reference.Rds"))
      
    }
  }
}


# Loop to explore the trial's conditions

# First, explore R0 with N = 72,000 (70 clusters)

r0_list <- c(0.2, 0.3, 0.4, 0.5, 0.7, 1, 1.2, 1.5, 2)
var_list <- c(0.01, 0.2)
cover_list <- c(0.5, 0.7, 0.9)

for (i in 1:length(r0_list)) {
  
  for (j in 1:length(var_list)) {
    
    for (k in 1:length(cover_list)) {
      
      run <- main(N = 72000, C = 100, var = var_list[j],
                  incidence = incidence, birth = birth, death = death,
                  R0 = r0_list[i], dur_inf = dur_inf, per_local = 0.5,
                  p_vax = cover_list[k], p_clusvax = p_clusvax, vax_eff = 0.7,
                  p_sym = p_sym, p_test = p_test, p_positive = p_positive,
                  years1 = 3, years2 = 2, n_runs = 30)
      
      # Directory
      
      dir.create(here("Results/", Sys.Date()),recursive = TRUE)
      
      # Give name to simulation
      
      name <- paste0("R=", r0_list[i], " ", run[[1]])
      
      dir.create(here(paste0("Results/", Sys.Date(), "/", name)),recursive = TRUE)
      
      # Save
      
      write.table(run[[2]], file = paste0("Results/", Sys.Date(), "/", name,"/characteristics.txt"))
      
      png(paste0("Results/", Sys.Date(), "/", name,"/Cluster_Pop_Hist.png"),
          width = 9, height = 5, units = 'in', res = 600)
      hist(run[[3]], main = "Histogram of clusters' population",
           xlab = "Clusters' population", ylab = "Frequency of clusters")
      dev.off()
      
      png(paste0("Results/", Sys.Date(), "/", name,"/Infections_after_vax.png"),
          width = 14, height = 9, units = 'in', res = 600)
      print(run[[4]])
      dev.off()
      
      write.xlsx(as.data.frame(run[[5]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/Total_Infections.xlsx"))
      
      write.xlsx(as.data.frame(run[[6]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/RR_Direct_Effect.xlsx"))
      write.xlsx(as.data.frame(run[[7]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/RR_Indirect_Effect.xlsx"))
      write.xlsx(as.data.frame(run[[8]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/RR_Overall_Effect.xlsx"))
      write.xlsx(as.data.frame(run[[9]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/RR_Total_Effect.xlsx"))
      write.xlsx(as.data.frame(run[[10]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/Summary_RR_Effects.xlsx"))
      
      write.xlsx(as.data.frame(run[[11]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/IntraCC_DesignEffect.xlsx"))
      write.xlsx(as.data.frame(run[[12]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/Other_sum.xlsx"))
      
      saveRDS(run[[13]], file = paste0("Results/", Sys.Date(), "/", name,"/All_data_reference.Rds"))
      
    }
  }
}

# Then, explore C with N = 72,000 with R0 = 0.2

cluster_list <- c(36, 52, 68, 84, 120)
var_list <- c(0.01, 0.2)
cover_list <- c(0.7, 0.75, 0.80, 0.85)

for (i in 1:length(cluster_list)) {
  
  for (j in 1:length(var_list)) {
    
    for (k in 1:length(cover_list)) {
      
      run <- main(N = 72000, C = cluster_list[i], var = var_list[j],
                  incidence = incidence, birth = birth, death = death,
                  R0 = 0.2, dur_inf = dur_inf, per_local = 0.5,
                  p_vax = cover_list[k], p_clusvax = p_clusvax, vax_eff = 0.7,
                  p_sym = p_sym, p_test = p_test, p_positive = p_positive,
                  years1 = 3, years2 = 2, n_runs = 30)
      
      # Directory
      
      dir.create(here("Results/", Sys.Date()),recursive = TRUE)
      
      # Give name to simulation
      
      name <- run[[1]]
      
      dir.create(here(paste0("Results/", Sys.Date(), "/", name)),recursive = TRUE)
      
      # Save
      
      write.table(run[[2]], file = paste0("Results/", Sys.Date(), "/", name,"/characteristics.txt"))
      
      png(paste0("Results/", Sys.Date(), "/", name,"/Cluster_Pop_Hist.png"),
          width = 9, height = 5, units = 'in', res = 600)
      hist(run[[3]], main = "Histogram of clusters' population",
           xlab = "Clusters' population", ylab = "Frequency of clusters")
      dev.off()
      
      png(paste0("Results/", Sys.Date(), "/", name,"/Infections_after_vax.png"),
          width = 14, height = 9, units = 'in', res = 600)
      print(run[[4]])
      dev.off()
      
      write.xlsx(as.data.frame(run[[5]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/Total_Infections.xlsx"))
      
      write.xlsx(as.data.frame(run[[6]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/RR_Direct_Effect.xlsx"))
      write.xlsx(as.data.frame(run[[7]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/RR_Indirect_Effect.xlsx"))
      write.xlsx(as.data.frame(run[[8]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/RR_Overall_Effect.xlsx"))
      write.xlsx(as.data.frame(run[[9]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/RR_Total_Effect.xlsx"))
      write.xlsx(as.data.frame(run[[10]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/Summary_RR_Effects.xlsx"))
      
      write.xlsx(as.data.frame(run[[11]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/IntraCC_DesignEffect.xlsx"))
      write.xlsx(as.data.frame(run[[12]]), rowNames = TRUE,
                 paste0("Results/", Sys.Date(), "/", name,"/Other_sum.xlsx"))
      
      saveRDS(run[[13]], file = paste0("Results/", Sys.Date(), "/", name,"/All_data_reference.Rds"))
      
    }
  }
}