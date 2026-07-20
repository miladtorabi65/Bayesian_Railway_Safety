# ============================================================
# FIT BAYESIAN MODELS TO THE REAL DATA
# M1: Poisson
# M2: Negative Binomial
# M3: Hierarchical Negative Binomial
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

library(rjags)
library(coda)


# ============================================================
# 2. Reproducibility
# ============================================================

set.seed(2026)


# ============================================================
# 3. Paths
# ============================================================

DATA_PATH <- paste0(
  "../data/processed/",
  "railway_country_year_2010_2024.csv"
)

MODEL_DIR <- "../models"

RESULT_DIR <- "../results/mcmc"

dir.create(
  RESULT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 4. MCMC settings
# ============================================================

N_CHAINS <- 4

N_ADAPT <- 5000

N_BURNIN <- 10000

N_SAMPLES <- 50000

THIN <- 1


cat(
  "MCMC configuration\n"
)

cat(
  "Chains:",
  N_CHAINS,
  "\n"
)

cat(
  "Adaptation:",
  N_ADAPT,
  "\n"
)

cat(
  "Burn-in:",
  N_BURNIN,
  "\n"
)

cat(
  "Samples per chain:",
  N_SAMPLES,
  "\n"
)

cat(
  "Thinning:",
  THIN,
  "\n"
)


# ============================================================
# 5. Load the frozen analysis dataset
# ============================================================

data <- read.csv(
  DATA_PATH
)


# ============================================================
# 6. Sort the observations
# ============================================================

data <- data[
  order(
    data$country_id,
    data$year
  ),
]


# ============================================================
# 7. Validation
# ============================================================

stopifnot(
  nrow(data) == 405
)

stopifnot(
  length(
    unique(
      data$country_id
    )
  ) == 27
)

stopifnot(
  length(
    unique(
      data$year
    )
  ) == 15
)

stopifnot(
  all(
    data$train_km > 0
  )
)

stopifnot(
  all(
    data$accidents >= 0
  )
)

stopifnot(
  !anyNA(data)
)


cat(
  "\nDataset validation passed.\n"
)

cat(
  "Observations:",
  nrow(data),
  "\n"
)

cat(
  "Countries:",
  length(
    unique(
      data$country_id
    )
  ),
  "\n"
)


# ============================================================
# 8. JAGS data
# ============================================================

N_country <- length(
  unique(
    data$country_id
  )
)


jags_data <- list(
  
  N = nrow(
    data
  ),
  
  N_country = N_country,
  
  accidents = as.integer(
    data$accidents
  ),
  
  train_km = as.numeric(
    data$train_km
  ),
  
  year_centered = as.numeric(
    data$year_centered
  ),
  
  country_id = as.integer(
    data$country_id
  )
)



# ============================================================
# 9. Country lookup table
# ============================================================

country_lookup <- unique(
  
  data[
    ,
    c(
      "country_id",
      "country",
      "country_name"
    )
  ]
)


country_lookup <- country_lookup[
  
  order(
    country_lookup$country_id
  ),
  
]


write.csv(
  
  country_lookup,
  
  file.path(
    RESULT_DIR,
    "country_lookup.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 10. Initial-value function
# ============================================================

make_initial_values <- function(
    model_name,
    chain_seed
) {
  
  common_values <- list(
    
    alpha = rnorm(
      1,
      mean = log(0.5),
      sd = 0.25
    ),
    
    beta = rnorm(
      1,
      mean = 0,
      sd = 0.03
    ),
    
    .RNG.name = (
      "base::Mersenne-Twister"
    ),
    
    .RNG.seed = (
      chain_seed
    )
  )
  
  
  # ==========================================================
  # Model 1
  # ==========================================================
  
  if (
    model_name
    == "M1_Poisson"
  ) {
    
    return(
      common_values
    )
  }
  
  
  # ==========================================================
  # Model 2
  # ==========================================================
  
  if (
    model_name
    == "M2_Negative_Binomial"
  ) {
    
    common_values$log_r <- rnorm(
      
      1,
      
      mean = log(10),
      
      sd = 0.4
    )
    
    
    return(
      common_values
    )
  }
  
  
  # ==========================================================
  # Model 3
  # ==========================================================
  
  if (
    model_name
    ==
    "M3_Hierarchical_Negative_Binomial"
  ) {
    
    common_values$log_r <- rnorm(
      1,
      mean = log(10),
      sd = 0.4
    )
    
    
    common_values$sigma_country <- runif(
      1,
      min = 0.2,
      max = 1.0
    )
    
    
    common_values$u_raw <- rnorm(
      N_country,
      mean = 0,
      sd = 0.5
    )
    
    
    return(
      common_values
    )
  }
  
  
  stop(
    paste(
      "Unknown model:",
      model_name
    )
  )
}




# ============================================================
# 11. Model-fitting function
# ============================================================

fit_jags_model <- function(
    
  model_name,
  
  model_file,
  
  parameters_to_monitor,
  
  seed_offset = 0
  
) {
  
  
  cat(
    "\n",
    paste(
      rep(
        "=",
        60
      ),
      collapse = ""
    ),
    "\n"
  )
  
  
  cat(
    "Fitting:",
    model_name,
    "\n"
  )
  
  
  cat(
    paste(
      rep(
        "=",
        60
      ),
      collapse = ""
    ),
    "\n"
  )
  
  
  # ==========================================================
  # Create output folder
  # ==========================================================
  
  model_result_dir <- file.path(
    
    RESULT_DIR,
    
    model_name
  )
  
  
  dir.create(
    
    model_result_dir,
    
    recursive = TRUE,
    
    showWarnings = FALSE
  )
  
  
  # ==========================================================
  # Initial values
  # ==========================================================
  
  chain_seeds <- c(
    
    1101,
    
    2202,
    
    3303,
    
    4404
    
  ) + seed_offset
  
  
  initial_values <- lapply(
    
    chain_seeds,
    
    function(
    seed
    ) {
      
      make_initial_values(
        
        model_name = model_name,
        
        chain_seed = seed
        
      )
    }
  )
  
  
  # ==========================================================
  # Compile model
  # ==========================================================
  
  cat(
    "Compiling model...\n"
  )
  
  
  jags_model <- jags.model(
    
    file = model_file,
    
    data = jags_data,
    
    inits = initial_values,
    
    n.chains = N_CHAINS,
    
    n.adapt = N_ADAPT
  )
  
  
  # ==========================================================
  # Burn-in
  # ==========================================================
  
  cat(
    "Running burn-in...\n"
  )
  
  
  update(
    
    jags_model,
    
    n.iter = N_BURNIN
  )
  
  
  # ==========================================================
  # Posterior sampling
  # ==========================================================
  
  cat(
    "Drawing posterior samples...\n"
  )
  
  
  sampling_time <- system.time(
    
    samples <- coda.samples(
      
      model = jags_model,
      
      variable.names = (
        parameters_to_monitor
      ),
      
      n.iter = N_SAMPLES,
      
      thin = THIN
    )
  )
  
  
  # ==========================================================
  # Save posterior samples
  # ==========================================================
  
  saveRDS(
    
    samples,
    
    file.path(
      
      model_result_dir,
      
      "posterior_samples.rds"
    )
  )
  
  
  # ==========================================================
  # Save textual MCMC summary
  # ==========================================================
  
  capture.output(
    
    summary(
      samples
    ),
    
    file = file.path(
      
      model_result_dir,
      
      "posterior_summary.txt"
    )
  )
  
  
  # ==========================================================
  # Save runtime
  # ==========================================================
  
  runtime_table <- data.frame(
    
    model = model_name,
    
    elapsed_seconds = (
      
      sampling_time[
        "elapsed"
      ]
    ),
    
    chains = N_CHAINS,
    
    adaptation_iterations = (
      N_ADAPT
    ),
    
    burnin_iterations = (
      N_BURNIN
    ),
    
    sampling_iterations = (
      N_SAMPLES
    ),
    
    thinning = THIN
  )
  
  
  write.csv(
    
    runtime_table,
    
    file.path(
      
      model_result_dir,
      
      "runtime_and_settings.csv"
    ),
    
    row.names = FALSE
  )
  
  
  cat(
    model_name,
    "completed.\n"
  )
  
  
  return(
    
    list(
      
      samples = samples,
      
      model = jags_model,
      
      runtime = sampling_time
      
    )
  )
}


# ============================================================
# 12. Model 1 parameters
# ============================================================

parameters_M1 <- c(
  
  "alpha",
  
  "beta",
  
  "baseline_rate",
  
  "annual_rate_ratio",
  
  "annual_percent_change"
)


# ============================================================
# 13. Model 2 parameters
# ============================================================

parameters_M2 <- c(
  
  "alpha",
  
  "beta",
  
  "baseline_rate",
  
  "annual_rate_ratio",
  
  "annual_percent_change",
  
  "log_r",
  
  "r"
)


# ============================================================
# 14. Model 3 parameters
# ============================================================

parameters_M3 <- c(
  
  "alpha",
  
  "beta",
  
  "baseline_rate",
  
  "annual_rate_ratio",
  
  "annual_percent_change",
  
  "log_r",
  
  "r",
  
  "sigma_country",
  
  "u",
  
  "country_rate_ratio"
)


# ============================================================
# 15. Fit Model 1 — Poisson
# ============================================================

fit_M1 <- fit_jags_model(
  
  model_name = (
    "M1_Poisson"
  ),
  
  model_file = file.path(
    
    MODEL_DIR,
    
    "model1_poisson.jags"
  ),
  
  parameters_to_monitor = (
    parameters_M1
  ),
  
  seed_offset = 0
)


samples_M1 <- (
  fit_M1$samples
)



# ============================================================
# 16. Fit Model 2 — Negative Binomial
# ============================================================

fit_M2 <- fit_jags_model(
  
  model_name = (
    "M2_Negative_Binomial"
  ),
  
  model_file = file.path(
    
    MODEL_DIR,
    
    "model2_negative_binomial.jags"
  ),
  
  parameters_to_monitor = (
    parameters_M2
  ),
  
  seed_offset = 10000
)


samples_M2 <- (
  fit_M2$samples
)

# ============================================================
# 17. Fit Model 3 — Hierarchical Negative Binomial
# ============================================================

fit_M3 <- fit_jags_model(
  
  model_name = (
    "M3_Hierarchical_Negative_Binomial"
  ),
  
  model_file = file.path(
    
    MODEL_DIR,
    
    "model3_hierarchical_negative_binomial.jags"
  ),
  
  parameters_to_monitor = (
    parameters_M3
  ),
  
  seed_offset = 20000
)


samples_M3 <- (
  fit_M3$samples
)


# ============================================================
# 18. Basic posterior summaries
# ============================================================

cat(
  "\n\nMODEL 1 SUMMARY\n"
)

print(
  summary(
    samples_M1
  )
)


cat(
  "\n\nMODEL 2 SUMMARY\n"
)

print(
  summary(
    samples_M2
  )
)


cat(
  "\n\nMODEL 3 SUMMARY\n"
)

print(
  summary(
    samples_M3
  )
)


# ============================================================
# 19. Save session information
# ============================================================

capture.output(
  
  sessionInfo(),
  
  file = file.path(
    
    RESULT_DIR,
    
    "R_session_info.txt"
  )
)


exists("fit_M3")
exists("fit_M3") && !is.null(fit_M3$model)

parameters_M3_full <- c(
  "alpha",
  "beta",
  "baseline_rate",
  "annual_rate_ratio",
  "annual_percent_change",
  "log_r",
  "r",
  "sigma_country",
  "u",
  "country_rate_ratio"
)

samples_M3_full <- coda.samples(
  
  model = fit_M3$model,
  
  variable.names = parameters_M3_full,
  
  n.iter = 20000,
  
  thin = 1
)




colnames(
  as.matrix(
    samples_M3_full
  )
)




sum(
  startsWith(
    colnames(
      as.matrix(
        samples_M3_full
      )
    ),
    "u["
  )
)


saveRDS(
  
  samples_M3_full,
  
  file.path(
    RESULT_DIR,
    "M3_Hierarchical_Negative_Binomial",
    "posterior_samples_full.rds"
  )
)









