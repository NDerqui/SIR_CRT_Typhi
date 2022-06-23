
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
vax_eff <- 0.7     # Vaccine effectiveness (infection)


#### Detected infections ####

p_sym <- 0.55       # Probability of being symptomatic - CHECK
p_test <- 0.50      # Probability of seeking a test - CHECK
p_positive <- 0.60  # Probability of test being positive


#### Time frame ####

time_step <- 1      # Time step change (days)


#### Population ####

# Total population in the study and number of clusters

N <- 200000                       # Total population in the study
C <- 100                  # Number of clusters

n <- round(rnorm(n = C, mean = N/C, sd = 100), digits = 0)  # Pop in each cluster
cluster_n <- abs(n)                                         # Vector of cluster populations


#### Clusters ####

random_cluster <- 0

# If 1 the first 50 clusters are vaccinated, but randomly distributed in map
# If 0 the vaccine clusters are allocated like a chessboard

# Cluster (and vaccine allocation) distribution

if (random_cluster == 1) {
  
  # Vector of clusters
  
  cluster_no <- seq(1:C)
  
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
  # Vaccination of the first half of clusters (randmly located)
  
  V <- C*p_clusvax                    
  
  if (C %% 2 == 0) {
    cluster_vstatus <- c(rep(1, times = V), rep(0, times = V))   
  } else {
    cluster_vstatus <- c(rep(1, times = V+1), rep(0, times = V)) 
  }
  
} else {
  
  # Vector of clusters
  
  cluster_no <- seq(1:C)
  
  # Cluster map
  # Location of clusters is not random, follows an order
  
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
  
  V <- C*p_clusvax                      
  
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
diag(cluster_dis) <- 0




# FUNCTION -------------------------------------------------------------


main <- function(N, C, cluster_no, cluster_dis, # Cluster characteristics
                 cluster_n, cluster_vstatus, 
                 incidence, birth, death,       # Incidence, birth and death rate
                 R0, dur_inf, imp_rate,         # Infection parameters
                 p_vax, p_clusvax, vax_eff,     # Vaccination parameters
                 p_sym, p_test, p_positive,     # Observed infections
                 time_step, years1, years2,     # Time: years for equi, years for vax
                 n_runs) {   
  
  
  ## Considerations

  # Calculated parameters
  
  beta <- R0/dur_inf                  # Force of infection
  mu <- p_sym*p_test*p_positive       # Prob of detecting I
  
  
  
  ## Cluster data for reference
  
  cluster_data <- data.frame(cluster = cluster_no,
                             vaccine = cluster_vstatus,
                             pop = cluster_n,
                             cluster_dis)
  colnames(cluster_data) <- c("cluster", "vaccine", "pop",
                              paste0("dis_", cluster_no))
  
  # Histogram of populations
  
  cluster_pop <- hist(cluster_n, main = "Histogram of clusters' population",
                      xlab = "Clusters' population", ylab = "Frequency of clusters")
  
  # Cluster map
  
  cluster_pop_map <- ggplot(data = data.frame(number = cluster_no, pop = cluster_n, vax = cluster_vstatus)) +
    geom_point(aes(x = pop, y = vax, color = as.factor(vax)), size = rel(20)) +
    geom_label(aes(x = pop, y = vax, label = pop, color = as.factor(vax)), size = rel(5)) +
    xlim(c(0.6*mean(cluster_n), 1.4*mean(cluster_n))) +
    ylim(c(-10, 10)) +
    theme_void() +
    labs(x = "Cluster population",
         y = NULL) +
    scale_color_manual(name = "Cluster vaccination status",
                       labels=c("Non-vaccinated", "Vaccinated"),
                       values = c("firebrick", "limegreen")) +
    theme(
      axis.title.x = element_text(size = rel(1.1), face="bold"),
      axis.title.y = element_text(size = rel(1.1), face="bold"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = rel(1), face="bold"),
      legend.text = element_text(size=rel(1)),
      strip.text.x = element_text(size = rel(2))) +
    facet_wrap(~ fct_relevel(as.character(number), as.character(cluster_map)),
                 ncol = nrow(cluster_map)) 
  cluster_pop_map
  
  
  ## For several simulations
  
  sir_output <- data.frame(cluster = 0, vaccine = 0, time_seq = 0,
                           infected = 0, observed = 0,
                           inc_sus_inf = 0, inc_vax_inf = 0, run = 0)
  
  for (n in 1:n_runs) {
    
    
    ## Step 1: reach equilibrium
    
    # Total time
    
    time_seq1 <- seq(from = 1, to = 365*years1, by = time_step)
    
    # Build the array
    
    names_row1 <- paste0("t_", time_seq1)
    
    names_column1 <- c("cluster", "vax_status", "time_seq", 
                       "no_N", "no_S", "no_I", "no_R",
                       "haz_inf","prob_infec", "inc_SI",
                       "prob_recov", "inc_IR",
                       "prob_death", "inc_SD", "inc_ID", "inc_RD",
                       "prob_birth", "inc_NB")
    
    names_matrix1 <- paste0("cluster_", cluster_no)
    
    sir_first <- array(0, dim = c(length(time_seq1), 18, C),
                       dimnames = list(names_row1, names_column1, names_matrix1))
    
    # Assign initial values
    
    sir_first[, 3, ] <- time_seq1
    
    for (i in 1:C) {
      sir_first[, 1, i] = cluster_no[i]
      sir_first[, 2, i] = cluster_vstatus[i]
      sir_first[, 4, i] = cluster_n[i]
      sir_first[, 6, i] = round(incidence*sir_first[, 4, i], digits = 0)           # Initial I depends on incidence rate
      sir_first[, 5, i] = round((sir_first[, 4 ,i]-sir_first[, 6, i]), digits = 0) # Initial S depends on N - I
      sir_first[, 7, i] = 0
    }
    
    # Put the model to run
    
    for (j in 1:C) {
      
      for (i in 2:length(time_seq1)) {
        
        # From S to I
        
        sir_first[i, 8, j] = beta*sir_first[i-1, 6, j]/sir_first[i-1, 4, j]     # Initial FOI: from cluster
        for (k in 1:C) {                                                        # Loop to add external FOI
          if (k != j) {
            sir_first[i, 8, j] = sir_first[i, 8, j] + (imp_rate/cluster_dis[k,j])*beta*sir_first[i-1, 6, k]/sir_first[i-1, 4 ,k] 
          }
        }
        sir_first[i, 9, j] = (1 - exp(-sir_first[i, 8, j]*time_step))
        sir_first[i, 10, j] = round(rbinom(n = 1, size = sir_first[i-1, 5, j], prob = sir_first[i, 9, j]), digits = 0)
        
        # From I to R
        
        sir_first[i, 11, j] = (1 - exp(-(1/dur_inf)*time_step))  
        sir_first[i, 12, j] = round(rbinom(n = 1, size = sir_first[i-1, 6, j], prob = sir_first[i, 11, j]), digits = 0) 
        
        # Deaths
        
        sir_first[i, 13, j] = (1 - exp(-death*time_step))  
        sir_first[i, 14, j] = round(rbinom(n = 1, size = sir_first[i-1, 5, j], prob = sir_first[i, 13, j]), digits = 0) 
        sir_first[i, 15, j] = round(rbinom(n = 1, size = sir_first[i-1, 6, j], prob = sir_first[i, 13, j]), digits = 0) 
        sir_first[i, 16, j] = round(rbinom(n = 1, size = sir_first[i-1, 7, j], prob = sir_first[i, 13, j]), digits = 0) 
        
        # Births
        
        sir_first[i, 17, j] = (1 - exp(-birth*time_step))  
        sir_first[i, 18, j] = round(rbinom(n = 1, size = sir_first[i-1, 4, j], prob = sir_first[i, 17, j]), digits = 0) 
        
        # Model equations
        sir_first[i, 5, j] = sir_first[i-1, 5, j] - sir_first[i, 10, j] - sir_first[i, 14, j] + sir_first[i, 18, j]
        sir_first[i, 6, j] = sir_first[i-1, 6, j] + sir_first[i, 10, j] - sir_first[i, 12, j] - sir_first[i, 15, j]
        sir_first[i, 7, j] = sir_first[i-1, 7, j] + sir_first[i, 12, j] - sir_first[i, 16, j]
        sir_first[i, 4, j] = sir_first[i, 5, j] + sir_first[i, 6, j] + sir_first[i, 7, j]
        
      }
    }
    
    #  Store the results: input the last row for the vaccine simulation
    
    equilibrium_result <- as.data.frame(matrix(0, nrow = C, ncol = 18))
    
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
                       "haz_inf", "prob_infec", "inc_SI", "prob_vaxin", "inc_VI",
                       "prob_recov", "inc_IR",
                       "prob_death", "inc_SD", "inc_VD", "inc_ID", "inc_RD",
                       "prob_birth", "inc_NB")
    
    names_matrix2 <- paste0("cluster_", cluster_no)
    
    sir <- array(0, dim = c(length(time_seq2), 22, C),
                 dimnames = list(names_row2, names_column2, names_matrix2))
    
    # Assign initial values
    
    sir[, 3, ] <- time_seq2
    
    for (i in 1:C) {
      sir[, 1, i] = cluster_no[i]
      sir[, 2, i] = cluster_vstatus[i]
      sir[, 4, i] = equilibrium_result[i, 4] # Total N from simulation
      sir[, 7, i] = equilibrium_result[i, 6] # I from the simulation
      sir[, 8, i] = equilibrium_result[i, 7] # R from the simulation
    }
    
    for (i in 1:C) {
      if (sir[1, 2, i] == 1) {                                                # In vaccinated clusters
        sir[, 5, i] = round(equilibrium_result[i, 5]*(1 - p_vax), digits = 0) # S=S*(1-COVERAGE)
        sir[, 6, i] = round(equilibrium_result[i, 5]*p_vax, digits = 0)       # V=S*COVERAGE
        
      } else {                                   # In non-vax clusters
        sir[, 5, i] = equilibrium_result[i, 5]   # Susceptible are all S from simulation
        sir[, 6, i] = as.integer(0)              # None vaccinated
      }
    }
    
    # Put the model to run
    
    for (j in 1:C) {
      
      for (i in 2:length(time_seq2)) {
        
        # From S to I
        
        sir[i, 9, j] = beta*sir[i-1, 7, j]/sir[i-1, 4 ,j]
        for (k in 1:C) { 
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
    sir_res_infected <- as.data.frame(matrix(0, nrow = length(time_seq2), ncol = C))
    colnames(sir_res_infected) <- names_matrix2
    for (i in 1:C) {
      sir_res_infected[,i] <- sir[,7,i]         
    }
    
    ## Detected infections
    sir_res_observed <- as.data.frame(matrix(0, nrow = length(time_seq2), ncol = C))
    colnames(sir_res_observed) <- names_matrix2
    for (i in 1:C) {
      sir_res_observed[,i] <- round(sir[,7,i]*mu, digits = 0)
    }
    
    ## Incidence of infection from S
    sir_res_inc_SI <- as.data.frame(matrix(0, nrow = length(time_seq2), ncol = C))
    colnames(sir_res_inc_SI) <- names_matrix2
    for (i in 1:C) {
      sir_res_inc_SI[,i] <- sir[,11,i]
    }
    
    ## Incidence of infection from V
    sir_res_inc_VI <- as.data.frame(matrix(0, nrow = length(time_seq2), ncol = C))
    colnames(sir_res_inc_VI) <- names_matrix2
    for (i in 1:C) {
      sir_res_inc_VI[,i] <- sir[,13,i]
    }
    
    ## Pivot the results and merge all together
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
    
    sir_result <- merge(sir_res_infected, sir_res_observed,
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
                        by = c("cluster", "vaccine", "time_seq", "infected",
                               "observed", "inc_sus_inf", "inc_vax_inf", "run"),
                        all = TRUE)
    
    sir_output <- sir_output %>%
      select("run", "cluster", "vaccine", "time_seq",
             "infected", "observed", "inc_sus_inf", "inc_vax_inf") %>%
      arrange(run, cluster, time_seq)
    
  }
  
  sir_output <- sir_output[2:nrow(sir_output),]
  
  
  
  ## Results from all the runs
  
  # Plot
  
  sir_vax_plot <- ggplot(data = sir_output) +
    geom_line(data = filter(sir_output, vaccine == 1),
              mapping = aes(x = time_seq, y = infected, group = run,
                            color = "Inf_Vax")) +
    geom_line(data = filter(sir_output, vaccine == 0),
              mapping = aes(x = time_seq, y = infected, group = run,
                            color = "Inf_No")) +
    geom_line(data = filter(sir_output, vaccine == 1),
              mapping = aes(x = time_seq, y = observed, group = run,
                            color = "Obs_Vax")) +
    geom_line(data = filter(sir_output, vaccine == 0),
              mapping = aes(x = time_seq, y = observed, group = run,
                            color = "Obs_No")) +
    scale_color_manual(name = NULL,
                       breaks = c("Inf_Vax", "Inf_No", "Obs_Vax", "Obs_No"),
                       values = c("Inf_Vax" = "limegreen",
                                  "Inf_No" = "firebrick",
                                  "Obs_Vax" = "greenyellow",
                                  "Obs_No" = "red"),
                       labels = c("Infections in vaccine clusters",
                                  "Infections in non-vaccine clusters",
                                  "Detected infections in vaccine clusters",
                                  "Detected infections in non-vaccine cluster")) +
    theme_classic() +
    labs(title = paste0("Incidence over time in n = ", C, " clusters"),
         x = "Time (days)",
         y = "Number of infections/detected infections") +
    theme(
      plot.title = element_text(size = rel(1.2), face="bold", hjust = 0.5),
      axis.title.x = element_text(size = rel(1.1), face="bold"),
      axis.title.y = element_text(size = rel(1.1), face="bold"),
      axis.text = element_text(size=rel(1)),
      legend.position = "bottom",
      legend.text = element_text(size=rel(1))) +
    facet_wrap(~fct_relevel(as.character(cluster), as.character(cluster_map)),
               ncol = sqrt(C), nrow = sqrt(C))
  
  # Summarize the results
  
  summary <- sir_output %>%
    group_by(run, cluster) %>%
    mutate(sum_SI = sum(inc_sus_inf)) %>%
    mutate(sum_VI = sum(inc_vax_inf)) %>%
    mutate(sum_all_inc = sum_SI + sum_VI) %>%
    filter(row_number() == 1) %>%
    select(-time_seq, -infected, -observed, -inc_sus_inf, -inc_vax_inf) %>%
    ungroup()
  
  
  
  ## Poisson regressions
  
  
  # Direct effect: vax vs unvax in vaccine cluster
  
  sir_direct <- summary %>%
    
    # Only measured in vaccine clusters
    filter(vaccine == 1) %>%
    select(run, cluster, vaccine, sum_SI, sum_VI) %>%
    
    # Pivot to get the incidence SI vs VI only
    pivot_longer(-c(run, cluster, vaccine),
                 names_to = "vax_incluster", values_to = "incidence")
  
  # Output vector to store the results of the test
  
  output_direct <- matrix(0, ncol = 6, nrow = n_runs)
  
  for (i in 1:n_runs) {
    
    model <- glm(formula = incidence ~ vax_incluster,
                 family = "poisson", data = filter(sir_direct, run == i))
    
    # Do the robust standard errors
    
    robust <- robust <- coeftest(model, vcov. = sandwich)
    
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
    
  sir_indirect <- summary %>%
      
    # Focus therefore only on incidence S to I
    select(run, cluster, vaccine, sum_SI) 
    
  # Output vector to store the results of the test
    
  output_indirect <- matrix(0, ncol = 6, nrow = n_runs)
    
  for (i in 1:n_runs) {
      
  model <- glm(formula = sum_SI ~ vaccine,
               family = "poisson", data = filter(sir_indirect, run == i))
      
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
    
  sir_total <- summary %>%
      select(run, cluster, vaccine, sum_SI, sum_VI, sum_all_inc) %>%
      
      # Get one var with VI only from vaccine clusters,
      # and all incidence from non-vaccine
      mutate(incidence = case_when(vaccine == 1 ~ sum_VI,
                                   vaccine == 0 ~ sum_all_inc))
    
  # Output vector to store the results of the test
    
  output_total <- matrix(0, ncol = 6, nrow = n_runs)
    
  for (i in 1:n_runs) {
      
    model <- glm(formula = incidence ~ vaccine,
                 family = "poisson", data = filter(sir_total, run == i))
      
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
    
 sir_overall <- summary %>%
      
    # We are interested in overall incidence
    select(run, cluster, vaccine, sum_all_inc)
    
  # Output vector to store the results of the test
    
  output_overall <- matrix(0, ncol = 6, nrow = n_runs)
    
    for (i in 1:n_runs) {
      
      model <- glm(formula = sum_all_inc ~ vaccine,
                   family = "poisson", data = filter(sir_overall, run == i))
      
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
  
  poisson_summary[1, 1:3] <- MeanCI(output_direct[,1])
  poisson_summary[2, 1:3] <- MeanCI(output_indirect[,1])
  poisson_summary[3, 1:3] <- MeanCI(output_total[,1])
  poisson_summary[4, 1:3] <- MeanCI(output_overall[,1])
  
  poisson_summary[1, 4] <- MeanCI(output_direct[,4])[1]
  poisson_summary[2, 4] <- MeanCI(output_indirect[,4])[1]
  poisson_summary[3, 4] <- MeanCI(output_total[,4])[1]
  poisson_summary[4, 4] <- MeanCI(output_overall[,4])[1]
  
  rownames(poisson_summary) <- c("Direct", "Indirect", "Total", "Overall")
  colnames(poisson_summary) <- c("Mean Effect", "Lower CI", "Upper CI", "P-value")
  
  poisson_summary <- round(poisson_summary, digits = 4)
  
  
  
  ## Returned objects
  
  name_simulation <- paste0("N=", N, " C=", C, " VE=", vax_eff, " Cover=", p_vax, " VaxArm=", p_clusvax)
  
  other_characteristics <- paste0("Incidence=", incidence, " BirthRate=", birth, " DeathRate=", death, " R0=", R0, " DurationInf=", dur_inf, " YearsEq=", years1, " YearsVax=", years2)
  
  list <- list(name_simulation,
               other_characteristics,
               cluster_data,
               cluster_pop,
               cluster_pop_map,
               sir_vax_plot,
               output_direct,
               output_indirect,
               output_overall,
               output_total,
               poisson_summary
  )

  list
}






  
  
  






# CRT: SIMULATIONS --------------------------------------------------------


#### Save cluster details #### 

library(here)
library(openxlsx)
library(forcats)

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

png("Results/2step_simulation_nocompeting/Cluster_map.png",
    width = 10, height = 12, units = 'in', res = 600)
ggplot(data = data.frame(number = cluster_no, pop = cluster_n, vax = cluster_vstatus)) +
  geom_point(aes(x = pop, y = vax, color = as.factor(vax)), size = rel(30)) +
  geom_label(aes(x = pop, y = vax, label = pop, color = as.factor(vax)), size = rel(8)) +
  xlim(c(0.8*mean(cluster_n), 1.2*mean(cluster_n))) +
  ylim(c(-2, 3)) +
  theme_classic() +
  labs(x = "Cluster population",
       y = NULL) +
  scale_color_manual(name = "Cluster vaccination status",
                     labels=c("Non-vaccinated", "Vaccinated"),
                     values = c("firebrick", "limegreen")) +
  theme(
    axis.title.x = element_text(size = rel(1.1), face="bold"),
    axis.title.y = element_text(size = rel(1.1), face="bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = rel(1), face="bold"),
    legend.text = element_text(size=rel(1))) +
  facet_wrap(~fct_relevel(as.character(number), as.character(cluster_map)),
             ncol = sqrt(C), nrow = sqrt(C))
dev.off()


#### Run simulations ####

# With these characteristics

R0
p_vax 
vax_eff
N
C

# Equilibrium

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
                time_step = time_step, years = 2,
                equilibrium_result = result_equilibrium, n_runs = 10)

# Graph

plot_sir <- sir_graph(test, cluster_map)

# Get summary numbers

summary_sir <- sir_summary(test, n_runs = 10)

# Poisson regression

direct_effect <- sir_direct(summary_sir, n_runs = 10)
indirect_effect <- sir_indirect(summary_sir, n_runs = 10)
total_effect <- sir_total(summary_sir, n_runs = 10)
overall_effect <- sir_overall(summary_sir, n_runs = 10)


#### Save results ####

# Give name to simulation

characteristics <- "N=200k C=100 sd=100 R0=2 Cover=0.8 VE=0.8"

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
