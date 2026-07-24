# ============================================================
# SENSITIVITY ANALYSIS
# Robustness check: refit M3 after excluding the year 2020
# ============================================================

library(rjags)
library(coda)
library(ggplot2)

DATA_PATH <- "../data/processed/railway_country_year_2010_2024.csv"
MODEL_FILE <- "../models/model3_hierarchical_negative_binomial.jags"
MAIN_SUMMARY_PATH <- "../tables/posterior_inference/posterior_main_summary.csv"

RESULT_DIR <- "../results/sensitivity_analysis"
TABLE_DIR <- "../tables/sensitivity_analysis"
FIGURE_DIR <- "../figures/sensitivity_analysis"

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(2026)

# Shorter than the final real-data fit because this is a focused
# robustness check. Convergence is still checked explicitly.
N_CHAINS <- 4
N_ADAPT <- 3000
N_BURNIN <- 5000
N_SAMPLES <- 10000
THIN <- 1

# ============================================================
# 1. Load data and exclude 2020
# ============================================================

data_full <- read.csv(DATA_PATH)

data_full <- data_full[
  order(
    data_full$country_id,
    data_full$year
  ),
]

data_sensitivity <- data_full[
  data_full$year != 2020,
]

data_sensitivity <- data_sensitivity[
  order(
    data_sensitivity$country_id,
    data_sensitivity$year
  ),
]

stopifnot(
  nrow(data_sensitivity) == 378,
  length(unique(data_sensitivity$country_id)) == 27,
  length(unique(data_sensitivity$year)) == 14,
  !anyNA(data_sensitivity),
  all(data_sensitivity$train_km > 0)
)

N_country <- length(
  unique(
    data_sensitivity$country_id
  )
)

stopifnot(
  identical(
    sort(
      unique(
        as.integer(
          data_sensitivity$country_id
        )
      )
    ),
    seq_len(N_country)
  )
)

cat("Sensitivity dataset validation passed.\n")
cat("Rows:", nrow(data_sensitivity), "\n")
cat("Countries:", N_country, "\n")

# ============================================================
# 2. JAGS data
# ============================================================

jags_data <- list(
  N = nrow(data_sensitivity),
  N_country = N_country,
  accidents = as.integer(data_sensitivity$accidents),
  train_km = as.numeric(data_sensitivity$train_km),
  year_centered = as.numeric(data_sensitivity$year_centered),
  country_id = as.integer(data_sensitivity$country_id)
)

# ============================================================
# 3. Initial values
# ============================================================

make_initial_values <- function(chain_seed) {

  list(
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

    log_r = rnorm(
      1,
      mean = log(10),
      sd = 0.4
    ),

    sigma_country = runif(
      1,
      min = 0.2,
      max = 1.0
    ),

    u_raw = rnorm(
      N_country,
      mean = 0,
      sd = 0.5
    ),

    .RNG.name = "base::Mersenne-Twister",
    .RNG.seed = chain_seed
  )
}

chain_seeds <- c(
  5101,
  6202,
  7303,
  8404
)

initial_values <- lapply(
  chain_seeds,
  make_initial_values
)

# ============================================================
# 4. Parameters to monitor
# ============================================================

parameters_to_monitor <- c(
  "alpha",
  "beta",
  "baseline_rate",
  "annual_rate_ratio",
  "annual_percent_change",
  "r",
  "sigma_country"
)

# ============================================================
# 5. Fit sensitivity model
# ============================================================

cat("\nCompiling sensitivity model...\n")

jags_model <- jags.model(
  file = MODEL_FILE,
  data = jags_data,
  inits = initial_values,
  n.chains = N_CHAINS,
  n.adapt = N_ADAPT
)

cat("Running burn-in...\n")

update(
  jags_model,
  n.iter = N_BURNIN
)

cat("Drawing sensitivity posterior samples...\n")

sampling_time <- system.time(

  samples_sensitivity <- coda.samples(
    model = jags_model,
    variable.names = parameters_to_monitor,
    n.iter = N_SAMPLES,
    thin = THIN
  )
)

saveRDS(
  samples_sensitivity,
  file.path(
    RESULT_DIR,
    "m3_excluding_2020_posterior_samples.rds"
  )
)

# ============================================================
# 6. MCMC diagnostics
# ============================================================

diagnostic_parameters <- c(
  "alpha",
  "beta",
  "r",
  "sigma_country"
)

diagnostic_samples <- samples_sensitivity[
  ,
  diagnostic_parameters
]

gelman_results <- gelman.diag(
  diagnostic_samples,
  autoburnin = FALSE,
  multivariate = FALSE
)

rhat_table <- gelman_results$psrf

ess_results <- effectiveSize(
  diagnostic_samples
)

diagnostics <- data.frame(
  parameter = diagnostic_parameters,

  rhat = rhat_table[
    diagnostic_parameters,
    "Point est."
  ],

  rhat_upper = rhat_table[
    diagnostic_parameters,
    "Upper C.I."
  ],

  ess = ess_results[
    diagnostic_parameters
  ]
)

diagnostics$status <- ifelse(
  diagnostics$rhat < 1.01
  &
  diagnostics$ess >= 1000,
  "Good",
  "Review"
)

write.csv(
  diagnostics,
  file.path(
    TABLE_DIR,
    "m3_excluding_2020_mcmc_diagnostics.csv"
  ),
  row.names = FALSE
)

cat("\nSensitivity-analysis MCMC diagnostics:\n")
print(diagnostics, row.names = FALSE)

# ============================================================
# 7. Posterior summary helper
# ============================================================

summarize_draws <- function(
    draws,
    parameter_name
) {

  data.frame(
    parameter = parameter_name,
    mean = mean(draws),
    median = median(draws),
    sd = sd(draws),
    lower_95 = unname(
      quantile(
        draws,
        0.025
      )
    ),
    upper_95 = unname(
      quantile(
        draws,
        0.975
      )
    )
  )
}

# ============================================================
# 8. Posterior summaries
# ============================================================

posterior <- as.matrix(
  samples_sensitivity
)

beta_draws <- posterior[
  ,
  "beta"
]

annual_rate_ratio_draws <- exp(
  beta_draws
)

annual_percent_change_draws <- (
  100
  *
  (
    annual_rate_ratio_draws
    -
    1
  )
)

r_draws <- posterior[
  ,
  "r"
]

sigma_country_draws <- posterior[
  ,
  "sigma_country"
]

sensitivity_summary <- rbind(

  summarize_draws(
    beta_draws,
    "beta"
  ),

  summarize_draws(
    annual_rate_ratio_draws,
    "annual_rate_ratio"
  ),

  summarize_draws(
    annual_percent_change_draws,
    "annual_percent_change"
  ),

  summarize_draws(
    r_draws,
    "r"
  ),

  summarize_draws(
    sigma_country_draws,
    "sigma_country"
  )
)

write.csv(
  sensitivity_summary,
  file.path(
    TABLE_DIR,
    "m3_excluding_2020_posterior_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 9. Posterior probability of decline
# ============================================================

sensitivity_probability <- data.frame(
  analysis = "M3_excluding_2020",

  probability_beta_below_zero = mean(
    beta_draws < 0
  ),

  probability_annual_decline = mean(
    annual_rate_ratio_draws < 1
  )
)

write.csv(
  sensitivity_probability,
  file.path(
    TABLE_DIR,
    "m3_excluding_2020_probability_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 10. Compare with the main M3 analysis
# ============================================================

main_summary <- read.csv(
  MAIN_SUMMARY_PATH
)

main_M3 <- main_summary[
  main_summary$model
  ==
  "M3_Hierarchical_Negative_Binomial",
]

parameters_for_comparison <- c(
  "annual_percent_change",
  "r",
  "sigma_country"
)

main_comparison <- main_M3[
  main_M3$parameter
  %in%
  parameters_for_comparison,
  c(
    "parameter",
    "median",
    "lower_95",
    "upper_95"
  )
]

main_comparison$analysis <- "Main analysis"

sensitivity_comparison <- sensitivity_summary[
  sensitivity_summary$parameter
  %in%
  parameters_for_comparison,
  c(
    "parameter",
    "median",
    "lower_95",
    "upper_95"
  )
]

sensitivity_comparison$analysis <- "Excluding 2020"

comparison_table <- rbind(
  main_comparison,
  sensitivity_comparison
)

comparison_table <- comparison_table[
  ,
  c(
    "analysis",
    "parameter",
    "median",
    "lower_95",
    "upper_95"
  )
]

write.csv(
  comparison_table,
  file.path(
    TABLE_DIR,
    "main_vs_excluding_2020_comparison.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 11. Comparison figure
# ============================================================

plot_data <- comparison_table[
  comparison_table$parameter
  ==
  "annual_percent_change",
]

sensitivity_plot <- ggplot(
  plot_data,
  aes(
    x = analysis,
    y = median
  )
) +

  geom_errorbar(
    aes(
      ymin = lower_95,
      ymax = upper_95
    ),
    width = 0.15
  ) +

  geom_point(
    size = 3
  ) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +

  labs(
    title = "Sensitivity of the Estimated Annual Railway Accident Trend",
    subtitle = "Main hierarchical model compared with the same model excluding 2020",
    x = NULL,
    y = "Annual percentage change (%)"
  ) +

  theme_minimal(
    base_size = 12
  )

print(
  sensitivity_plot
)

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "annual_trend_sensitivity_excluding_2020.png"
  ),
  plot = sensitivity_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# 12. Print final results
# ============================================================

cat("\n========================================\n")
cat("SENSITIVITY ANALYSIS: EXCLUDING 2020\n")
cat("========================================\n\n")

cat("Posterior summary:\n\n")
print(
  sensitivity_summary,
  row.names = FALSE
)

cat("\nPosterior probability of decline:\n\n")
print(
  sensitivity_probability,
  row.names = FALSE
)

cat("\nMain analysis versus excluding 2020:\n\n")
print(
  comparison_table,
  row.names = FALSE
)

cat(
  "\nSampling runtime in seconds:",
  sampling_time["elapsed"],
  "\n"
)

# ============================================================
# 13. Save runtime and session information
# ============================================================

runtime_table <- data.frame(
  analysis = "M3_excluding_2020",
  elapsed_seconds = sampling_time["elapsed"],
  chains = N_CHAINS,
  adaptation_iterations = N_ADAPT,
  burnin_iterations = N_BURNIN,
  sampling_iterations = N_SAMPLES,
  thinning = THIN
)

write.csv(
  runtime_table,
  file.path(
    RESULT_DIR,
    "runtime_and_settings.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    RESULT_DIR,
    "R_session_info.txt"
  )
)

cat(
  "\nSensitivity analysis completed successfully.\n"
)
