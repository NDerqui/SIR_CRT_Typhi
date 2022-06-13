
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
                        incidence, birth, death,
                        R0, dur_inf, imp_rate,
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
  
  ## Hazards 
  
  haz_inf <- rep(0, times = length(time_seq))

  ##  Probabilities: 1 - exp(hazard)
  
  prob_event <- rep(0, times = length(time_seq))
  prob_x_to_I <- rep(0, times = length(time_seq))
  prob_x_to_D <- rep(0, times = length(time_seq))
  prob_recov <- rep(0, times = length(time_seq))
  prob_death <- rep(0, times = length(time_seq))
  prob_birth <- rep(0, times = length(time_seq))
  
  ## State variables: incidence/recovery
  
  inc_Sx <- as.integer(rep(0, times = length(time_seq)))
  inc_SI <- as.integer(rep(0, times = length(time_seq)))
  inc_SD <- as.integer(rep(0, times = length(time_seq)))
  inc_IR <- as.integer(rep(0, times = length(time_seq)))
  inc_ID <- as.integer(rep(0, times = length(time_seq)))
  inc_RD <- as.integer(rep(0, times = length(time_seq)))
  inc_NB <- as.integer(rep(0, times = length(time_seq)))
  
  
  # Build the array
  
  ## Basic structure
  
  names_row <- paste0("t_", time_seq)
  
  names_column <- c("cluster", "vax_status", "time_seq", # 1 2 3
                    "no_N", "no_S", "no_I", "no_R",      # 4 5 6 7
                    "haz_inf", "prob_event", "inc_Sx",   # 8 9 10
                    "prob_x_to_I", "inc_SI",             # 11 12
                    "prob_x_to_D", "inc_SD",             # 13 14
                    "prob_recov", "inc_IR",              # 15 16
                    "prob_death", "inc_ID", "inc_RD",    # 17 18 19
                    "prob_birth", "inc_NB")              # 20 21    
  
  names_matrix <- paste0("cluster_", cluster_no)
  
  sir <- array(c(cluster, v_cluster, time_seq,
                 no_N, no_S, no_I, no_R,
                 haz_inf, prob_event, inc_Sx,
                 prob_x_to_I, inc_SI,
                 prob_x_to_D, inc_SD,
                 prob_recov, inc_IR,
                 prob_death, inc_ID, inc_RD,
                 prob_birth, inc_NB),
               dim = c(length(time_seq), length(names_column), C),
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
      
      # From S to event(x)
      
      sir[i, 8, j] = beta*sir[i-1, 6, j]/sir[i-1, 4, j]     # Initial FOI: from cluster
      for (k in 1:C) {                                      # Loop to add external FOI
        if (k != j) {
          sir[i, 8, j] = sir[i, 8, j] + (imp_rate/cluster_dis[k,j])*beta*sir[i-1, 6, k]/sir[i-1, 4 ,k] 
        }
      }
      sir[i, 9, j] = (1 - exp(-(sir[i, 8, j] + death)*time_step))
      sir[i, 10, j] = round(rbinom(n = 1, size = sir[i-1, 5, j], prob = sir[i, 9, j]), digits = 0)
      
      # From event to I
      
      sir[i, 11, j] = sir[i, 8, j]/(sir[i, 8, j] + death)  
      sir[i, 12, j] = round(rbinom(n = 1, size = sir[i-1, 10, j], prob = sir[i, 11, j]), digits = 0) 
      
      # From event to D
      
      sir[i, 13, j] = death/(sir[i, 8, j] + death)  
      sir[i, 14, j] = round(rbinom(n = 1, size = sir[i-1, 10, j], prob = sir[i, 13, j]), digits = 0) 
      
      # From I to R
      
      sir[i, 15, j] = (1 - exp(-(1/dur_inf)*time_step))  
      sir[i, 16, j] = round(rbinom(n = 1, size = sir[i-1, 6, j], prob = sir[i, 15, j]), digits = 0) 
      
      # Deaths
      
      sir[i, 17, j] = (1 - exp(-death*time_step))  
      sir[i, 18, j] = round(rbinom(n = 1, size = sir[i-1, 6, j], prob = sir[i, 17, j]), digits = 0) 
      sir[i, 19, j] = round(rbinom(n = 1, size = sir[i-1, 7, j], prob = sir[i, 17, j]), digits = 0) 
      
      # Births
      
      sir[i, 20, j] = (1 - exp(-birth*time_step))  
      sir[i, 21, j] = round(rbinom(n = 1, size = sir[i-1, 4, j], prob = sir[i, 17, j]), digits = 0) 
      
      # Model equations
      sir[i, 5, j] = sir[i-1, 5, j] - sir[i, 10, j] + sir[i, 21, j]
      sir[i, 6, j] = sir[i-1, 6, j] + sir[i, 12, j] - sir[i, 16, j] - sir[i, 18, j]
      sir[i, 7, j] = sir[i-1, 7, j] + sir[i, 16, j] - sir[i, 19, j]
      sir[i, 4, j] = sir[i, 5, j] + sir[i, 6, j] + sir[i, 7, j]
      
    }
  }
  
  
  # Store the results
  
  ## Graph to check
  
  sir_result_graph <- as.data.frame(matrix(0, nrow = length(time_seq), ncol = C))
  colnames(sir_result_graph) <- names_matrix
  for (i in 1:C) {
    sir_result_graph[,i] <- sir[,6,i]
  }

  sir_result_graph <- sir_result_graph %>%
    mutate(time_seq = time_seq) %>%
    pivot_longer(-time_seq, names_to = "cluster", values_to = "infected") %>%
    mutate(num_cluster = readr::parse_number(cluster)) %>%
    arrange(num_cluster)
  
  plot <- ggplot(data = sir_result_graph) +
    geom_line(mapping = aes(x = time_seq, y = infected, group = num_cluster),
              color = "firebrick") +
    theme_classic() +
    labs(title = "Incidence over time",
         x = paste0("Time over ", years, " years (days)"),
         y = "Number of infections") +
    theme(
      plot.title = element_text(size = rel(1.2), face="bold", hjust = 0.5),
      axis.title.x = element_text(size = rel(1.1), face="bold"),
      axis.title.y = element_text(size = rel(1.1), face="bold"),
      axis.text = element_text(size=rel(1)))
  
  ## Input for the vaccine simulation
  
  sir_result <- as.data.frame(matrix(0, nrow = C, ncol = 21))
  colnames(sir_result) <- names_column
  for (i in 1:C) {
    sir_result[i,] <- sir[length(time_seq),,i]         
  }
  
  ## Returned object
  
  sir_list <- list(plot, sir_result)
  
}



# FUNCTIONS: CRT ----------------------------------------------------------


#### Main functions ####

sir_model <- function(N, C, cluster_no, cluster_n, cluster_vstatus, cluster_dis,
                      R0, dur_inf, imp_rate,
                      p_vax, p_clusvax, vax_eff,
                      p_sym, p_test, p_positive,
                      birth, death,
                      time_step, years, equilibrium_result) {

  
  time_seq <- seq(from = 1, to = 365*years, by = time_step) # Total time   
    
  
  # Calculated parameters
  
  beta <- R0/dur_inf                  # Force of infection
  mu <- p_sym*p_test*p_positive       # Prob of detecting I
  
  
  # Cluster vectors

  cluster_no <- cluster_no            # Vector of clusters
  cluster_n <- cluster_n              # Vector of cluster populations
  cluster_vstatus <- cluster_vstatus  # Flag for vax clusters
  
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
  
  ## SIR model compartments
  
  no_N <- as.integer(rep(0, times = length(time_seq)))     # Total n per cluster
  no_S <- as.integer(rep(0, times = length(time_seq)))     # Susceptible 
  no_V <- as.integer(rep(0, times = length(time_seq)))     # Vaccinated  
  no_I <- as.integer(rep(0, times = length(time_seq)))     # Infected
  no_R <- as.integer(rep(0, times = length(time_seq)))     # Recovered
  
  ## Hazards 
  
  haz_inf <- rep(0, times = length(time_seq))
  
  ##  Probabilities: 1 - exp(hazard)
  
  prob_infec <- rep(0, times = length(time_seq))
  prob_vaxin <- rep(0, times = length(time_seq))
  prob_recov <- rep(0, times = length(time_seq))
  prob_death <- rep(0, times = length(time_seq))
  prob_birth <- rep(0, times = length(time_seq))
  
  ## State variables: incidence/recovery
  
  inc_SI <- as.integer(rep(0, times = length(time_seq)))
  inc_VI <- as.integer(rep(0, times = length(time_seq)))
  inc_IR <- as.integer(rep(0, times = length(time_seq)))
  inc_SD <- as.integer(rep(0, times = length(time_seq)))
  inc_VD <- as.integer(rep(0, times = length(time_seq)))
  inc_ID <- as.integer(rep(0, times = length(time_seq)))
  inc_RD <- as.integer(rep(0, times = length(time_seq)))
  inc_NB <- as.integer(rep(0, times = length(time_seq)))
  
  
  # Build the array
  
  ## Basic structure
  
  names_row <- paste0("t_", time_seq)
  names_column <- c("cluster", "vaccine", "time_seq",
                    "no_N", "no_S", "no_V", "no_I", "no_R",
                    "haz_inf", "prob_infec", "inc_SI", "prob_vaxin", "inc_VI",
                    "prob_recov", "inc_IR",
                    "prob_death", "inc_SD", "inc_VD", "inc_ID", "inc_RD",
                    "prob_birth", "inc_NB")
  names_matrix <- paste0("cluster_", cluster_no)
  sir <- array(c(cluster, v_cluster, time_seq, no_N, no_S, no_V, no_I, no_R,
                 haz_inf, prob_infec, inc_SI, prob_vaxin, inc_VI, prob_recov, inc_IR,
                 prob_death, inc_SD, inc_VD, inc_ID, inc_RD, prob_birth, inc_NB),
               dim = c(length(time_seq), 22, C),
               dimnames = list(names_row, names_column, names_matrix))
  
  ## Assign initial values
  
  # Cluster number, vaccine status and no_N
  for (i in 1:C) {
    sir[, 1, i] = cluster_no[i]
    sir[, 2, i] = cluster_vstatus[i]
    sir[, 4, i] = equilibrium_result[i, 4]
  }
  # I an R
  for (i in 1:C) {
    sir[, 7, i] = equilibrium_result[i, 6]
    sir[, 8, i] = equilibrium_result[i, 7]
  }
  # S and V population - depends on vaccination cluster status
  for (i in 1:C) {
    if (sir[1, 2, i] == 1) {                                                # In vaccinated clusters
      sir[, 5, i] = round(equilibrium_result[i, 5]*(1 - p_vax), digits = 0) # Susceptible (non-vax)
      sir[, 6, i] = round(equilibrium_result[i, 5]*p_vax, digits = 0)       # Vaccinated
      
    } else {                                   # In non-vax clusters
      sir[, 5, i] = equilibrium_result[i, 5]   # Susceptible (non-vax)
      sir[, 6, i] = as.integer(0)              # Vaccinated
    }
  }
  

  # Put the model to run
  
  for (j in 1:C) {
    
    for (i in 2:length(time_seq)) {
      
      # From S to I
      
      sir[i, 9, j] = beta*sir[i-1, 7, j]/sir[i-1, 4 ,j]
      for (k in 1:C) { # Loop to add external FOI
       if (k != j) {
          sir[i, 9, j] = sir[i, 9, j] + (imp_rate/cluster_dis[k,j])*beta*sir[i-1, 7, k]/sir[i-1, 4 ,k] 
        }
      }
      sir[i, 10, j] = (1 - exp(-sir[i, 9, j]*time_step))
      sir[i, 11, j] = round(rbinom(n = 1, size = sir[i-1, 5, j], prob = sir[i, 10, j]), digits = 0)
      
      # From V to I
      
      sir[i, 12, j] = (1 - exp(-sir[i, 9, j]*(1 - vax_eff)*time_step))
      sir[i, 13, j] = round(rbinom(n = 1, size = sir[i-1, 6, j], prob = sir[i, 12, j]), digits = 0)
      
      # From I to R
      
      sir[i, 14, j] = (1 - exp(-(1/dur_inf)*time_step))  
      sir[i, 15, j] = round(rbinom(n = 1, size = sir[i-1, 7, j], prob = sir[i, 14, j]), digits = 0)  
      
      # Deaths
      
      sir[i, 16, j] = (1 - exp(-death*time_step))  
      sir[i, 17, j] = round(rbinom(n = 1, size = sir[i-1, 5, j], prob = sir[i, 16, j]), digits = 0) 
      sir[i, 18, j] = round(rbinom(n = 1, size = sir[i-1, 6, j], prob = sir[i, 16, j]), digits = 0) 
      sir[i, 19, j] = round(rbinom(n = 1, size = sir[i-1, 7, j], prob = sir[i, 16, j]), digits = 0) 
      sir[i, 20, j] = round(rbinom(n = 1, size = sir[i-1, 8, j], prob = sir[i, 16, j]), digits = 0) 
      
      # Births
      
      sir[i, 21, j] = (1 - exp(-birth*time_step))  
      sir[i, 22, j] = round(rbinom(n = 1, size = sir[i-1, 4, j], prob = sir[i, 21, j]), digits = 0) 
      
      # Model equations
      sir[i, 5, j] = sir[i-1, 5, j] - sir[i, 11, j] - sir[i, 17, j] + round(sir[i, 22, j]*(1-p_vax), digits = 0)
      sir[i, 6, j] = sir[i-1, 6, j] - sir[i, 13, j] - sir[i, 18, j] + round(sir[i, 22, j]*(p_vax), digits = 0)
      sir[i, 7, j] = sir[i-1, 7, j] + sir[i, 11, j] + sir[i, 13, j] - sir[i, 15, j] - sir[i, 19, j]
      sir[i, 8, j] = sir[i-1, 8, j] + sir[i, 15, j] - sir[i, 20, j]
      sir[i, 4, j] = sir[i, 5, j] + sir[i, 6, j] + sir[i, 7, j] + sir[i, 8, j]
      
    }
  }
  
  
  # Store the results
  
  ## Infected
  
  sir_res_infected <- as.data.frame(matrix(0, nrow = length(time_seq), ncol = C))
  colnames(sir_res_infected) <- names_matrix
  
  for (i in 1:C) {
    sir_res_infected[,i] <- sir[,7,i]         
  }
  
  ## Detected infections
  
  sir_res_observed <- as.data.frame(matrix(0, nrow = length(time_seq), ncol = C))
  colnames(sir_res_observed) <- names_matrix
  
  for (i in 1:C) {
    sir_res_observed[,i] <- round(sir[,7,i]*mu, digits = 0)
  }
  
  ## Incidence of infection from S
  
  sir_res_inc_SI <- as.data.frame(matrix(0, nrow = length(time_seq), ncol = C))
  colnames(sir_res_inc_SI) <- names_matrix
  
  for (i in 1:C) {
    sir_res_inc_SI[,i] <- sir[,11,i]
  }
  
  ## Incidence of infection from V
  
  sir_res_inc_VI <- as.data.frame(matrix(0, nrow = length(time_seq), ncol = C))
  colnames(sir_res_inc_VI) <- names_matrix
  
  for (i in 1:C) {
    sir_res_inc_VI[,i] <- sir[,13,i]
  }
  
  ## Pivot the results and merge all together
  
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
  
  sir_res_inc_SI <- sir_res_inc_SI %>%
    mutate(time_seq = time_seq) %>%
    pivot_longer(-time_seq, names_to = "cluster", values_to = "inc_sus_inf") %>%
    mutate(cluster = readr::parse_number(cluster)) %>%
    arrange(cluster) %>%
    merge(y = cluster_data[, c("cluster", "vaccine")], by = "cluster", all = TRUE) %>%
    select(c("cluster", "vaccine", "time_seq", "inc_sus_inf"))
  
  sir_res_inc_VI <- sir_res_inc_VI %>%
    mutate(time_seq = time_seq) %>%
    pivot_longer(-time_seq, names_to = "cluster", values_to = "inc_vax_inf") %>%
    mutate(cluster = readr::parse_number(cluster)) %>%
    arrange(cluster) %>%
    merge(y = cluster_data[, c("cluster", "vaccine")], by = "cluster", all = TRUE) %>%
    select(c("cluster", "vaccine", "time_seq", "inc_vax_inf"))
  
  sir_result <- merge(sir_res_infected, sir_res_observed,
                      by = c("cluster", "vaccine", "time_seq"), all = TRUE)
  sir_result <- merge(sir_result, sir_res_inc_SI,
                      by = c("cluster", "vaccine", "time_seq"), all = TRUE)
  sir_result <- merge(sir_result, sir_res_inc_VI,
                      by = c("cluster", "vaccine", "time_seq"), all = TRUE)
  sir_result <- arrange(sir_result, cluster, time_seq)
  
  
  # Returned object
  
  sir_result
}


# Function to run the SIR several times
# Append the results of each run to each other

sir_many <- function(..., n_runs) {
  
  # Generate output vector
  
  sir_result <- sir_model(...)
  sir_result <- sir_result %>%
    mutate(run = 1)
  sir_output <- sir_result
  
  for (i in 2:n_runs) {
    
    # Get result of one run
    sir_result <- sir_model(...)
    sir_result <- sir_result %>%
      mutate(run = i)
    
    # Merge to previous
    sir_output <- merge(sir_output, sir_result,
                        by = c("cluster", "vaccine", "time_seq",
                               "infected", "observed", "inc_sus_inf", "inc_vax_inf",
                               "run"), all = TRUE)
  }
  
  # Clean the result
  
  sir_many_output <- sir_output %>%
    select("run", "cluster", "vaccine", "time_seq",
           "infected", "observed", "inc_sus_inf", "inc_vax_inf") %>%
    arrange(run, cluster, time_seq)
  
  sir_many_output
}


#### Results functions ####


# Graph: mainly to check if run was okay

sir_graph <- function(sir_many_result) {
  
  ggplot() +
    geom_line(data = filter(sir_many_result, vaccine == 1),
              mapping = aes(x = time_seq, y = infected, group = run,
                            color = "Inf_Vax")) +
    geom_line(data = filter(sir_many_result, vaccine == 0),
              mapping = aes(x = time_seq, y = infected, group = run,
                            color = "Inf_No")) +
    geom_line(data = filter(sir_many_result, vaccine == 1),
              mapping = aes(x = time_seq, y = observed, group = run,
                            color = "Obs_Vax")) +
    geom_line(data = filter(sir_many_result, vaccine == 0),
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
    labs(title = "Incidence over time (SIR)",
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

sir_summary <- function(sir_many_result, n_runs) {
  
  summary <- sir_many_result %>%
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

sir_direct <- function(sir_summary_result, n_runs) {
  
  # Direct effect: vax vs unvax in vaccine cluster
  
  sir_summary_result <- sir_summary_result %>%
    
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
                 family = "poisson", data = filter(sir_summary_result, run == i))
    
    # Do the robust standard errors
    
    robust <- robust <- coeftest(model, vcov. = sandwich)
    
    # Get the coefficients of the model
    
    x <- exp(robust[2, 1])
    y <- robust[2, 2:4]
    z <- exp(confint(robust))
    
    # Store them
    
    output[i, 1]   <- x
    output[i, 2:4] <- y
    output[i, 5:6] <- z[2,]
  }
  
  # Clean the result vector
  
  rownames(output) <- paste0("run_", seq(1:10))
  colnames(output) <- c("exp(Estimate)", "SE", "z Val", "Pr(>|z|)",
                        "2.5%", "97.5%")
  
  output
}


# Poisson regression : indirect effects

sir_indirect <- function(sir_summary_result, n_runs) {
  
  # Indirect effect: unvax in vaccine clusters vs non-vaccine clusters
  
  sir_summary_result <- sir_summary_result %>%
    
    # Focus therefore only on incidence S to I
    select(run, cluster, vaccine, sum_SI) 
  
  # Output vector to store the results of the test
  
  output <- matrix(0, ncol = 6, nrow = n_runs)
  
  for (i in 1:n_runs) {
    
    model <- glm(formula = sum_SI ~ vaccine,
                 family = "poisson", data = filter(sir_summary_result, run == i))
    
    # Do the robust standard errors
    
    robust <- robust <- coeftest(model, vcov. = sandwich)
    
    # Get the coefficients of the model
    
    x <- exp(robust[2, 1])
    y <- robust[2, 2:4]
    z <- exp(confint(robust))
    
    # Store them
    
    output[i, 1]   <- x
    output[i, 2:4] <- y
    output[i, 5:6] <- z[2,]
  }
  
  # Clean the result vector
  
  rownames(output) <- paste0("run_", seq(1:10))
  colnames(output) <- c("exp(Estimate)", "SE", "z Val", "Pr(>|z|)",
                        "2.5%", "97.5%")
  
  output
}


# Poisson regression : total effects

sir_total <- function(sir_summary_result, n_runs) {
  
  # Total effect: vax in vaccine cluster vs unvax in non-vaccine
  
  sir_summary_result <- sir_summary_result %>%
    select(run, cluster, vaccine, sum_SI, sum_VI, sum_all_inc) %>%
    
    # Get one var with VI only from vaccine clusters,
    # and all incidence from non-vaccine
    mutate(incidence = case_when(vaccine == 1 ~ sum_VI,
                                 vaccine == 0 ~ sum_all_inc))
  
  # Output vector to store the results of the test
  
  output <- matrix(0, ncol = 6, nrow = n_runs)
  
  for (i in 1:n_runs) {
    
    model <- glm(formula = incidence ~ vaccine,
                 family = "poisson", data = filter(sir_summary_result, run == i))
    
    # Do the robust standard errors
    
    robust <- robust <- coeftest(model, vcov. = sandwich)
    
    # Get the coefficients of the model
    
    x <- exp(robust[2, 1])
    y <- robust[2, 2:4]
    z <- exp(confint(robust))
    
    # Store them
    
    output[i, 1]   <- x
    output[i, 2:4] <- y
    output[i, 5:6] <- z[2,]
  }
  
  # Clean the result vector
  
  rownames(output) <- paste0("run_", seq(1:10))
  colnames(output) <- c("exp(Estimate)", "SE", "z Val", "Pr(>|z|)",
                        "2.5%", "97.5%")
  
  output
}


# Poisson regression : overall effects

sir_overall <- function(sir_summary_result, n_runs) {
  
  # Overall effect: all in vaccine cluster vs all in non-vaccine
  
  sir_summary_result <- sir_summary_result %>%
    
    # We are interested in overall incidence
    select(run, cluster, vaccine, sum_all_inc)
  
  # Output vector to store the results of the test
  
  output <- matrix(0, ncol = 6, nrow = n_runs)
  
  for (i in 1:n_runs) {
    
    model <- glm(formula = sum_all_inc ~ vaccine,
                 family = "poisson", data = filter(sir_summary_result, run == i))
    
    # Do the robust standard errors
    
    robust <- robust <- coeftest(model, vcov. = sandwich)
    
    # Get the coefficients of the model
    
    x <- exp(robust[2, 1])
    y <- robust[2, 2:4]
    z <- exp(confint(robust))
    
    # Store them
    
    output[i, 1]   <- x
    output[i, 2:4] <- y
    output[i, 5:6] <- z[2,]
  }
  
  # Clean the result vector
  
  rownames(output) <- paste0("run_", seq(1:10))
  colnames(output) <- c("exp(Estimate)", "SE", "z Val", "Pr(>|z|)",
                        "2.5%", "97.5%")
  
  output
}



# CRT: SIMULATIONS --------------------------------------------------------


#### Save cluster details #### 

library(here)
library(openxlsx)

dir.create(here("Results"),recursive = TRUE)
dir.create(here("Results/2step_simulation_nocompeting"),recursive = TRUE)

# Save the cluster data

write.xlsx(as.data.frame(cluster_map), rowNames = TRUE,
           "Results/2step_simulation_nocompeting/Cluster_map.xlsx")

png("Results/2step_simulation_nocompeting/Cluster_populations.png",
    width = 10, height = 6, units = 'in', res = 600)
hist(cluster_n, main = "Histogram of clusters' population",
     xlab = "Clusters' population", ylab = "Frequency of clusters")
dev.off()

# With these characteristics

R0
p_vax <- 0.95
vax_eff


#### Run simulations ####

equilibrium <- equilibrium(N = N, C = C, cluster_no = cluster_no, cluster_n = cluster_n,
                           cluster_vstatus = cluster_vstatus, cluster_dis = cluster_dis,
                           incidence = incidence, birth = birth, death = death,
                           R0 = R0, dur_inf = dur_inf, imp_rate = imp_rate,
                           time_step = time_step, years = 5)

# Check results of equilibrium

equilibrium[[1]]

result_equilibrium <- equilibrium[[2]]

# Run the vaccine model

test <- sir_many(N = N, C = C, cluster_no = cluster_no, cluster_n = cluster_n,
                cluster_vstatus = cluster_vstatus, cluster_dis = cluster_dis,
                R0 = R0, dur_inf = dur_inf, imp_rate = imp_rate,
                p_vax = p_vax, p_clusvax = p_clusvax, vax_eff = vax_eff,
                p_sym = p_sym, p_test = p_test, p_positive = p_positive,
                birth = birth, death = death,
                time_step = time_step, years = 1,
                equilibrium_result = result_equilibrium, n_runs = 10)

# Graph

plot_sir <- sir_graph(test)

# Get summary numbers

summary_sir <- sir_summary(test, n_runs = 10)

# Poisson regression

direct_effect <- sir_direct(summary_sir, n_runs = 10)
indirect_effect <- sir_indirect(summary_sir, n_runs = 10)
total_effect <- sir_total(summary_sir, n_runs = 10)
overall_effect <- sir_overall(summary_sir, n_runs = 10)


#### Save results ####

# Give name to simulation

characteristics <- "R0=2 Cover=0.95 VE=0.8"

dir.create(here(paste0("Results/2step_simulation_nocompeting/", characteristics)),recursive = TRUE)

# Save

png(paste0("Results/2step_simulation_nocompeting/", characteristics,"/SIR_equilibrium.png"),
    width = 20, height = 12, units = 'in', res = 600)
equilibrium[[1]]
dev.off()

png(paste0("Results/2step_simulation_nocompeting/", characteristics,"/SIR_after_vax.png"),
    width = 20, height = 12, units = 'in', res = 600)
plot_sir
dev.off()

write.xlsx(as.data.frame(direct_effect), rowNames = TRUE,
           paste0("Results/2step_simulation_nocompeting/", characteristics,"/Direct_Effect.xlsx"))
write.xlsx(as.data.frame(indirect_effect), rowNames = TRUE,
           paste0("Results/2step_simulation_nocompeting/", characteristics,"/Indirect_Effect.xlsx"))
write.xlsx(as.data.frame(total_effect), rowNames = TRUE,
           paste0("Results/2step_simulation_nocompeting/", characteristics,"/Total_Effect.xlsx"))
write.xlsx(as.data.frame(overall_effect), rowNames = TRUE,
           paste0("Results/2step_simulation_nocompeting/", characteristics,"/Overall_Effect.xlsx"))
