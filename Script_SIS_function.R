
## Cluster Randomised Trial for a Typhoid vaccine in Vellore ##
      ## MRes Biomedical Research - EECID - Project 2 ##

                      ## SIS Model ##



# SET UP ------------------------------------------------------------------

rm(list = ls())

library(tidyverse)



# PARAMETERS --------------------------------------------------------------


#### Basic ####

# Total population and number of clusters

N <- 150000                       # Total population in the study
C <- 80                          # Number of clusters

# The population in each cluster should be the same across runs

n <- round(rnorm(n = C, mean = N/C, sd = 100), digits = 0)  # Pop in each cluster
cluster_n <- abs(n)                                         # Vector of cluster populations


#### Clusters ####

# Distance matrix across clusters (should be the same across runs)

cluster_dis <- matrix(1, nrow = C, ncol = C,                   # Matrix with distance between each pair of clusters
                      dimnames = list(seq(1:C), seq(1:C)))

cluster_dis[lower.tri(cluster_dis, diag = FALSE)] <- runif(    # Filled with random numbers in a symmetrical way
  n = (C^2 - C)/2, min = 1, max = 100)
cluster_dis[upper.tri(cluster_dis, diag = FALSE)] <- t(cluster_dis)[upper.tri(cluster_dis)]
cluster_dis <- round(cluster_dis, digits = 1)


#### Prevalence ####

# Starting seed of epidemic

number_infected <- as.integer(10) # Number of initially infected in the cluster


#### Force of infection ####

# (infections/time): beta = R0/Duration of infectiousness

R0 <- 2            # Basic reproduction number
dur_inf <- 3       # Duration of infectiousness (days)

# Importation rate: from external clusters to a given one

imp_rate <- 0.25


#### Vaccination ####

# coverage and effect

p_vax <- 0.5       # Proportion of vaccinated population in vaccine clusters
p_clusvax <- 0.5   # Proportion of clusters assigned to vaccine
vax_eff <- 0.8     # Vaccine effectiveness (infection)


#### Detected infections ####

p_sym <- 0.55       # Probability of being symptomatic - CHECK
p_test <- 0.50      # Probability of seeking a test - CHECK
p_positive <- 0.80  # Probability of test being positive


#### Time frame ####

time_step <- 1      # Time step change (days)
years <- 1          # Total duration of simulation (years)



# FUNCTION ----------------------------------------------------------------


# From SIR model develop SIS allowing R people to go back to S or V

sis_model <- function(N, C, cluster_n, cluster_dis,
                      number_infected, R0, dur_inf, imp_rate,
                      p_vax, p_clusvax, vax_eff, p_sym, p_test, p_positive,
                      time_step, years) {

  
  time_seq <- seq(from = 1, to = 365*years, by = time_step) # Total time   
    
  
  # Calculated parameters
  
  beta <- R0/dur_inf                                    # Force of infection
  
  mu <- p_sym*p_test*p_positive                         # Prob of detecting I
  
  
  # Cluster vectors

  cluster_no <- seq(1:C)                                      # Vector of clusters
  
  cluster_n <- cluster_n                                      # Vector of cluster populations
  
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
  
  # State variables: incidence/recovery
  
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
    sis[, 7, i] = number_infected
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
      sis[i, 5, j] = sis[i-1, 5, j] - sis[i, 14, j] + round(sis[i-1, 8, j]*(1 - p_vax), digits = 0)
      sis[i, 6, j] = sis[i-1, 6, j] - sis[i, 15, j] + round(sis[i-1, 8, j]*p_vax, digits = 0)
      
      sis[i, 7, j] = sis[i-1, 7, j] + sis[i, 14, j] + sis[i, 15, j] - sis[i, 16, j]
      
      sis[i, 8, j] = sis[i-1, 8, j] + sis[i, 16, j] - sis[i-1, 8, j]
      
      sis[i, 4, j] = sis[i, 5, j] + sis[i, 6, j] + sis[i, 7, j] + sis[i, 8, j]
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
  
  sis_result <- merge(sis_res_infected, sis_res_observed,
                      by = c("cluster", "vaccine", "time_seq"), all = TRUE)
  sis_result <- arrange(sis_result, cluster, time_seq)
  

  sis_result
}



# SEVERAL RUNS ------------------------------------------------------------


# Function to run the sis several times

# One run generates a matrix of time_points x clusters rows
# and five cols (cluster, vaccine, time_seq, infected and observed)

# To store the results of run we would need an array
# [number of rows, number of cols, num of runs]

sis_many <- function(..., n_runs) {
  
  sis_result <- sis_model(...)
  
  sis_result <- sis_result %>%
    mutate(run = 1)
  
  sis_output <- sis_result
  
  for (i in 2:n_runs) {
    
    sis_result <- sis_model(...)
    
    sis_result <- sis_result %>%
      mutate(run = i)
    
    sis_output <- merge(sis_output, sis_result,
                        by = c("cluster", "vaccine", "time_seq",
                               "infected", "observed", "run"), all = TRUE)
    
  }
  
  sis_output <- sis_output %>%
    select("run", "cluster", "vaccine", "time_seq", "infected", "observed") %>%
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
    mutate(mena_obs = mean(observed)) %>%
    mutate(sum_inf = sum(infected)) %>%
    mutate(sum_obs = sum(observed)) %>%
    filter(row_number() == 1) %>%
    select(-time_seq, -infected, -observed) %>%
    ungroup()
  
}


# Poisson regression

sis_stats <- function(sis_summary_result, n_runs) {
  
  output <- matrix(0, ncol = 6, nrow = n_runs)
  
  for (i in 1:n_runs) {
    
    model <- glm(formula = sum_obs ~ vaccine,
                 family = "poisson", data = filter(sis_summary_result, run == i))
    
    x <- exp(summary(model)$coef)
    
    y <- exp(confint(model))
    
    output[i, 1:4] <- x[2,]
    
    output[i, 5:6] <- y[2,]
  }
  
  rownames(output) <- paste0("run_", seq(1:10))
  
  colnames(output) <- c("exp(Estimate)", "SE", "z Val", "Pr(>|z|)",
                        "2.5%", "97.5%")
  
  output
}


# SAMPLE RUNS -------------------------------------------------------------


# How many simulations?

n_runs <- 10

# Run!

test_sis_80 <- sis_many(N = N, C = C, cluster_n = cluster_n, cluster_dis = cluster_dis,
                 number_infected = number_infected,
                 R0 = R0, dur_inf = dur_inf, imp_rate = imp_rate,
                 p_vax = p_vax, p_clusvax = p_clusvax, vax_eff = vax_eff,
                 p_sym = p_sym, p_test = p_test, p_positive = p_positive,
                 time_step = time_step, years = years,
                 n_runs = n_runs)

# Graph

plot_sis_80 <- sis_graph(test_sis_80)

# Get summary stats

summary_sis_80 <- sis_summary(test_sis_80, n_runs)

stats_sis_80 <- sis_stats(summary_sis_80, n_runs)
