# ============================================================
# PARAMETER-RECOVERY SMOKE TEST
# Model 3: Hierarchical Negative Binomial
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

library(rjags)
library(coda)
library(ggplot2)


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

MODEL_PATH <- paste0(
  "../models/",
  "model3_hierarchical_negative_binomial.jags"
)

RESULT_DIR <- paste0(
  "../results/",
  "parameter_recovery_smoke_test"
)

FIGURE_DIR <- paste0(
  "../figures/",
  "parameter_recovery_smoke_test"
)

dir.create(
  RESULT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 4. Load the frozen analysis dataset
# ============================================================

data <- read.csv(
  DATA_PATH
)

stopifnot(
  nrow(data) == 405
)

stopifnot(
  length(unique(data$country_id)) == 27
)

stopifnot(
  all(data$train_km > 0)
)

data <- data[
  order(
    data$country_id,
    data$year
  ),
]


# ============================================================
# 5. Known true parameter values
# ============================================================

alpha_true <- log(0.5)

beta_true <- log(0.97)

r_true <- 10

sigma_country_true <- 0.5


# Derived true quantities

baseline_rate_true <- exp(
  alpha_true
)

annual_rate_ratio_true <- exp(
  beta_true
)

annual_percent_change_true <- (
  100 *
    (
      annual_rate_ratio_true
      - 1
    )
)


cat(
  "True baseline rate:",
  baseline_rate_true,
  "\n"
)

cat(
  "True annual rate ratio:",
  annual_rate_ratio_true,
  "\n"
)

cat(
  "True annual percent change:",
  annual_percent_change_true,
  "%\n"
)

cat(
  "True NB dispersion r:",
  r_true,
  "\n"
)

cat(
  "True country SD:",
  sigma_country_true,
  "\n"
)


# ============================================================
# 6. Simulate true country effects
# ============================================================

N_country <- length(
  unique(
    data$country_id
  )
)


z_true <- rnorm(
  N_country,
  mean = 0,
  sd = 1
)


u_true <- (
  sigma_country_true
  * z_true
)


country_rate_ratio_true <- exp(
  u_true
)


# ============================================================
# 7. Calculate the true expected counts
# ============================================================

log_mu_true <- (
  log(
    data$train_km
  )
  +
    alpha_true
  +
    beta_true
  *
    data$year_centered
  +
    u_true[
      data$country_id
    ]
)


mu_true <- exp(
  log_mu_true
)


# ============================================================
# 8. Simulate negative-binomial accident counts
# ============================================================

# Same parameterization as JAGS:
#
# p = r / (r + mu)
#
# E(Y) = mu
#
# Var(Y) = mu + mu^2 / r


p_true <- (
  r_true
  /
    (
      r_true
      +
        mu_true
    )
)


accidents_simulated <- rnbinom(
  n = nrow(data),
  size = r_true,
  prob = p_true
)


# Quick checks

cat(
  "\nSimulated total accidents:",
  sum(accidents_simulated),
  "\n"
)

cat(
  "Simulated mean:",
  mean(accidents_simulated),
  "\n"
)

cat(
  "Simulated variance:",
  var(accidents_simulated),
  "\n"
)

cat(
  "Simulated zero counts:",
  sum(
    accidents_simulated == 0
  ),
  "\n"
)


# ============================================================
# 9. Save the simulated dataset
# ============================================================

simulated_data <- data

simulated_data$accidents_observed <- (
  simulated_data$accidents
)

simulated_data$accidents <- (
  accidents_simulated
)

simulated_data$mu_true <- (
  mu_true
)


write.csv(
  simulated_data,
  file.path(
    RESULT_DIR,
    "simulated_smoke_test_dataset.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. Prepare data for JAGS
# ============================================================

jags_data <- list(
  
  N = nrow(
    simulated_data
  ),
  
  N_country = N_country,
  
  accidents = as.integer(
    simulated_data$accidents
  ),
  
  train_km = as.numeric(
    simulated_data$train_km
  ),
  
  year_centered = as.numeric(
    simulated_data$year_centered
  ),
  
  country_id = as.integer(
    simulated_data$country_id
  )
)


# ============================================================
# 11. Initial values for four chains
# ============================================================

make_initial_values <- function(
    chain_seed
) {
  
  list(
    
    alpha = rnorm(
      1,
      mean = log(0.5),
      sd = 0.3
    ),
    
    beta = rnorm(
      1,
      mean = 0,
      sd = 0.03
    ),
    
    log_r = rnorm(
      1,
      mean = log(10),
      sd = 0.3
    ),
    
    sigma_country = runif(
      1,
      min = 0.1,
      max = 1
    ),
    
    z = rnorm(
      N_country
    ),
    
    .RNG.name = (
      "base::Mersenne-Twister"
    ),
    
    .RNG.seed = chain_seed
  )
}


initial_values <- list(
  
  make_initial_values(
    1001
  ),
  
  make_initial_values(
    2002
  ),
  
  make_initial_values(
    3003
  ),
  
  make_initial_values(
    4004
  )
)


# ============================================================
# 12. Compile the JAGS model
# ============================================================

jags_model <- jags.model(
  
  file = MODEL_PATH,
  
  data = jags_data,
  
  inits = initial_values,
  
  n.chains = 4,
  
  n.adapt = 2000
)


# ============================================================
# 13. Burn-in
# ============================================================

update(
  jags_model,
  n.iter = 5000
)


# ============================================================
# 14. Draw posterior samples
# ============================================================

parameters_to_monitor <- c(
  
  "alpha",
  
  "beta",
  
  "baseline_rate",
  
  "annual_rate_ratio",
  
  "annual_percent_change",
  
  "r",
  
  "sigma_country",
  
  "u"
)


samples <- coda.samples(
  
  model = jags_model,
  
  variable.names = (
    parameters_to_monitor
  ),
  
  n.iter = 10000,
  
  thin = 1
)


# ============================================================
# 15. Basic MCMC summary
# ============================================================

print(
  summary(
    samples
  )
)



# ============================================================
# 16. Scalar parameters for diagnostics
# ============================================================

scalar_parameters <- c(
  
  "alpha",
  
  "beta",
  
  "r",
  
  "sigma_country",
  
  "annual_rate_ratio"
)


scalar_samples <- samples[
  ,
  scalar_parameters
]


# ============================================================
# 17. Gelman-Rubin diagnostic
# ============================================================

rhat_results <- gelman.diag(
  
  scalar_samples,
  
  autoburnin = FALSE,
  
  multivariate = FALSE
)


print(
  rhat_results
)


# ============================================================
# 18. Effective sample size
# ============================================================

ess_results <- effectiveSize(
  scalar_samples
)


print(
  ess_results
)


# ============================================================
# 19. Trace plots
# ============================================================

png(
  file.path(
    FIGURE_DIR,
    "smoke_test_traceplots.png"
  ),
  width = 1600,
  height = 1200,
  res = 150
)


traceplot(
  scalar_samples
)


dev.off()



# ============================================================
# 20. Convert posterior samples to matrix
# ============================================================

posterior_matrix <- as.matrix(
  samples
)


# ============================================================
# 21. True values
# ============================================================

true_values <- c(
  
  alpha = alpha_true,
  
  beta = beta_true,
  
  baseline_rate = (
    baseline_rate_true
  ),
  
  annual_rate_ratio = (
    annual_rate_ratio_true
  ),
  
  annual_percent_change = (
    annual_percent_change_true
  ),
  
  r = r_true,
  
  sigma_country = (
    sigma_country_true
  )
)


# ============================================================
# 22. Recovery summary
# ============================================================

recovery_rows <- list()


for (
  parameter
  in names(
    true_values
  )
) {
  
  draws <- posterior_matrix[
    ,
    parameter
  ]
  
  
  lower <- quantile(
    draws,
    0.025
  )
  
  
  upper <- quantile(
    draws,
    0.975
  )
  
  
  recovery_rows[[
    parameter
  ]] <- data.frame(
    
    parameter = parameter,
    
    true_value = (
      true_values[
        parameter
      ]
    ),
    
    posterior_mean = (
      mean(draws)
    ),
    
    posterior_median = (
      median(draws)
    ),
    
    lower_95 = lower,
    
    upper_95 = upper,
    
    true_inside_95_interval = (
      
      true_values[
        parameter
      ]
      >= lower
      
      &
        
        true_values[
          parameter
        ]
      <= upper
    )
  )
}


recovery_summary <- do.call(
  
  rbind,
  
  recovery_rows
)


rownames(
  recovery_summary
) <- NULL


print(
  recovery_summary
)


write.csv(
  
  recovery_summary,
  
  file.path(
    RESULT_DIR,
    "parameter_recovery_summary.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 23. Extract posterior country effects
# ============================================================

u_columns <- grep(
  "^u\\[",
  colnames(
    posterior_matrix
  )
)


u_posterior <- posterior_matrix[
  ,
  u_columns,
  drop = FALSE
]


# ============================================================
# 24. Posterior summaries for country effects
# ============================================================

u_posterior_mean <- apply(
  
  u_posterior,
  
  2,
  
  mean
)


u_lower <- apply(
  
  u_posterior,
  
  2,
  
  quantile,
  
  probs = 0.025
)


u_upper <- apply(
  
  u_posterior,
  
  2,
  
  quantile,
  
  probs = 0.975
)


country_lookup <- unique(
  
  data[
    ,
    c(
      "country_id",
      "country_name"
    )
  ]
)


country_lookup <- country_lookup[
  
  order(
    country_lookup$country_id
  ),
  
]


country_recovery <- data.frame(
  
  country_id = (
    1:N_country
  ),
  
  country_name = (
    country_lookup$country_name
  ),
  
  true_u = u_true,
  
  posterior_mean = (
    u_posterior_mean
  ),
  
  lower_95 = u_lower,
  
  upper_95 = u_upper
)


country_recovery[
  "true_inside_95_interval"
] <- (
  
  country_recovery$true_u
  >=
    country_recovery$lower_95
  
  &
    
    country_recovery$true_u
  <=
    country_recovery$upper_95
)


print(
  country_recovery
)


write.csv(
  
  country_recovery,
  
  file.path(
    RESULT_DIR,
    "country_effect_recovery.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 25. Country-effect recovery plot
# ============================================================

recovery_plot <- ggplot(
  
  country_recovery,
  
  aes(
    x = true_u,
    y = posterior_mean
  )
  
) +
  
  geom_point(
    size = 3
  ) +
  
  geom_abline(
    
    slope = 1,
    
    intercept = 0,
    
    linetype = "dashed"
  ) +
  
  labs(
    
    title = "Recovery of Simulated Country Effects",
    
    x = (
      "True country effect"
    ),
    
    y = (
      "Posterior mean"
    )
  ) +
  
  theme_minimal()


print(
  recovery_plot
)


ggsave(
  
  filename = file.path(
    FIGURE_DIR,
    "country_effect_recovery.png"
  ),
  
  plot = recovery_plot,
  
  width = 8,
  
  height = 6,
  
  dpi = 300
)
