# ============================================================
# PARAMETER-RECOVERY STUDY
# Preferred model: M3 Hierarchical Negative Binomial
#
# I save progress after every attempted replication.
#
# Recommended workflow:
#   1. Keep N_TARGET_SUCCESS <- 2 for a short pilot run.
#   2. Inspect runtime and convergence.
#   3. Change N_TARGET_SUCCESS <- 20 and rerun this same file.
#      Existing successful replications will be kept and the
#      script will continue from where it stopped.
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

library(rjags)
library(coda)
library(ggplot2)


# ============================================================
# 2. Paths
# ============================================================

DATA_PATH <- paste0(
  "../data/processed/",
  "railway_country_year_2010_2024.csv"
)

MODEL_FILE <- paste0(
  "../models/",
  "model3_hierarchical_negative_binomial.jags"
)

RESULT_DIR <- "../results/parameter_recovery"

TABLE_DIR <- "../tables/parameter_recovery"

FIGURE_DIR <- "../figures/parameter_recovery"


dir.create(
  RESULT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. Recovery-study settings
# ============================================================

N_TARGET_SUCCESS <- 10


# Maximum total attempts allowed.
#
# When the final target is changed to 20, a value of 30 allows
# some replications to fail convergence without stopping the
# whole study immediately.

MAX_ATTEMPTS <- 30


# Moderate MCMC settings for repeated recovery.
#
# These are intentionally shorter than the final real-data fit.
# The purpose here is repeated simulation-based validation,
# not final real-data inference.

N_CHAINS <- 3

N_ADAPT <- 1000

N_BURNIN <- 2000

N_SAMPLES <- 3000

THIN <- 1


# Convergence criteria for accepting a replication

R_HAT_THRESHOLD <- 1.05

ESS_THRESHOLD <- 300


# ============================================================
# 4. Known true parameter values
# ============================================================

# These values are fixed for the entire repeated study.
#
# They represent a realistic hierarchical negative-binomial
# data-generating process with:
#
# - baseline rate around 0.5 accidents per million train-km;
# - about a 3% annual decline;
# - substantial country heterogeneity;
# - moderate residual negative-binomial overdispersion.

ALPHA_TRUE <- log(
  0.5
)

BETA_TRUE <- log(
  0.97
)

R_TRUE <- 30

SIGMA_COUNTRY_TRUE <- 0.8


TRUE_PARAMETER_VALUES <- c(

  alpha = ALPHA_TRUE,

  beta = BETA_TRUE,

  r = R_TRUE,

  sigma_country = SIGMA_COUNTRY_TRUE
)


# ============================================================
# 5. Reproducibility
# ============================================================

MASTER_SEED <- 2026

set.seed(
  MASTER_SEED
)


# ============================================================
# 6. Load and validate the frozen analysis design
# ============================================================

data <- read.csv(
  DATA_PATH
)


data <- data[
  order(
    data$country_id,
    data$year
  ),
]


stopifnot(
  nrow(
    data
  ) == 405
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
  !anyNA(
    data
  )
)

stopifnot(
  all(
    data$train_km > 0
  )
)


N <- nrow(
  data
)

N_COUNTRY <- length(
  unique(
    data$country_id
  )
)


country_ids <- sort(
  unique(
    as.integer(
      data$country_id
    )
  )
)


stopifnot(
  identical(
    country_ids,
    seq_len(
      N_COUNTRY
    )
  )
)


train_km <- as.numeric(
  data$train_km
)

year_centered <- as.numeric(
  data$year_centered
)

country_id <- as.integer(
  data$country_id
)


cat(
  "Frozen simulation design loaded successfully.\n"
)

cat(
  "Observations:",
  N,
  "\n"
)

cat(
  "Countries:",
  N_COUNTRY,
  "\n"
)


# ============================================================
# 7. Progress-file paths
# ============================================================

ATTEMPT_LOG_PATH <- file.path(
  RESULT_DIR,
  "recovery_attempt_log.csv"
)

GLOBAL_RESULTS_PATH <- file.path(
  RESULT_DIR,
  "recovery_global_parameter_results.csv"
)

COUNTRY_RESULTS_PATH <- file.path(
  RESULT_DIR,
  "recovery_country_effect_results.csv"
)

COUNTRY_METRICS_PATH <- file.path(
  RESULT_DIR,
  "recovery_country_effect_metrics_by_replication.csv"
)


# ============================================================
# 8. Helper: append and save a data frame
# ============================================================

append_and_save <- function(
    existing_data,
    new_data,
    path
) {

  if (
    is.null(
      existing_data
    )
  ) {

    combined_data <- new_data

  } else {

    combined_data <- rbind(
      existing_data,
      new_data
    )
  }


  write.csv(

    combined_data,

    path,

    row.names = FALSE
  )


  return(
    combined_data
  )
}


# ============================================================
# 9. Load previous progress if it exists
# ============================================================

attempt_log <- NULL

global_results <- NULL

country_results <- NULL

country_metrics <- NULL


if (
  file.exists(
    ATTEMPT_LOG_PATH
  )
) {

  attempt_log <- read.csv(
    ATTEMPT_LOG_PATH
  )
}


if (
  file.exists(
    GLOBAL_RESULTS_PATH
  )
) {

  global_results <- read.csv(
    GLOBAL_RESULTS_PATH
  )
}


if (
  file.exists(
    COUNTRY_RESULTS_PATH
  )
) {

  country_results <- read.csv(
    COUNTRY_RESULTS_PATH
  )
}


if (
  file.exists(
    COUNTRY_METRICS_PATH
  )
) {

  country_metrics <- read.csv(
    COUNTRY_METRICS_PATH
  )
}


if (
  is.null(
    attempt_log
  )
) {

  current_successes <- 0

  next_attempt_id <- 1

} else {

  current_successes <- sum(
    attempt_log$success
  )

  next_attempt_id <- (
    max(
      attempt_log$attempt_id
    )
    +
    1
  )
}


cat(
  "\nExisting successful replications:",
  current_successes,
  "\n"
)

cat(
  "Current target:",
  N_TARGET_SUCCESS,
  "\n"
)


# ============================================================
# 10. Generate true country effects
#
# ============================================================

generate_true_country_effects <- function(
    sigma_country,
    n_country
) {

  u_raw <- rnorm(

    n_country,

    mean = 0,

    sd = sigma_country
  )


  scale_adjust <- sqrt(

    n_country

    /

    (
      n_country
      -
      1
    )
  )


  u <- (

    scale_adjust

    *

    (
      u_raw
      -
      mean(
        u_raw
      )
    )
  )


  return(
    u
  )
}


# ============================================================
# 11. Generate chain-specific initial values
# ============================================================

make_recovery_initial_values <- function(
    chain_seed
) {

  list(

    alpha = rnorm(
      1,
      mean = log(
        0.5
      ),
      sd = 0.25
    ),

    beta = rnorm(
      1,
      mean = 0,
      sd = 0.03
    ),

    log_r = rnorm(
      1,
      mean = log(
        10
      ),
      sd = 0.4
    ),

    sigma_country = runif(
      1,
      min = 0.2,
      max = 1.2
    ),

    u_raw = rnorm(
      N_COUNTRY,
      mean = 0,
      sd = 0.5
    ),

    .RNG.name = "base::Mersenne-Twister",

    .RNG.seed = chain_seed
  )
}


# ============================================================
# 12. Helper: order u[1], ..., u[27]
# ============================================================

find_ordered_u_parameters <- function(
    parameter_names
) {

  u_parameters <- parameter_names[
    startsWith(
      parameter_names,
      "u["
    )
  ]


  u_numbers <- as.integer(

    sub(

      "u\\[([0-9]+)\\]",

      "\\1",

      u_parameters
    )
  )


  u_parameters[
    order(
      u_numbers
    )
  ]
}


# ============================================================
# 13. Fit one simulated replication
# ============================================================

run_one_recovery_replication <- function(
    attempt_id
) {

  replication_seed <- (
    MASTER_SEED
    +
    10000
    +
    attempt_id
  )


  set.seed(
    replication_seed
  )


  cat(
    "\n============================================================\n"
  )

  cat(
    "Starting recovery attempt:",
    attempt_id,
    "\n"
  )

  cat(
    "============================================================\n"
  )


  # ----------------------------------------------------------
  # A. Generate known true country effects
  # ----------------------------------------------------------

  u_true <- generate_true_country_effects(

    sigma_country = SIGMA_COUNTRY_TRUE,

    n_country = N_COUNTRY
  )


  # ----------------------------------------------------------
  # B. Generate synthetic accident counts
  # ----------------------------------------------------------

  mu_true <- (

    train_km

    *

    exp(

      ALPHA_TRUE

      +

      BETA_TRUE
      *
      year_centered

      +

      u_true[
        country_id
      ]
    )
  )


  accidents_simulated <- rnbinom(

    N,

    size = R_TRUE,

    mu = mu_true
  )


  # ----------------------------------------------------------
  # C. JAGS data
  # ----------------------------------------------------------

  jags_data <- list(

    N = N,

    N_country = N_COUNTRY,

    accidents = as.integer(
      accidents_simulated
    ),

    train_km = train_km,

    year_centered = year_centered,

    country_id = country_id
  )


  # ----------------------------------------------------------
  # D. Initial values
  # ----------------------------------------------------------

  chain_seeds <- (

    c(
      1101,
      2202,
      3303
    )

    +

    (
      attempt_id
      *
      100
    )
  )


  initial_values <- lapply(

    chain_seeds,

    make_recovery_initial_values
  )


  # ----------------------------------------------------------
  # E. Fit model safely
  # ----------------------------------------------------------

  fit_result <- tryCatch(

    {

      jags_model <- jags.model(

        file = MODEL_FILE,

        data = jags_data,

        inits = initial_values,

        n.chains = N_CHAINS,

        n.adapt = N_ADAPT,

        quiet = TRUE
      )


      update(

        jags_model,

        n.iter = N_BURNIN,

        progress.bar = "none"
      )


      samples <- coda.samples(

        model = jags_model,

        variable.names = c(
          "alpha",
          "beta",
          "r",
          "sigma_country",
          "u"
        ),

        n.iter = N_SAMPLES,

        thin = THIN,

        progress.bar = "none"
      )


      list(

        success = TRUE,

        samples = samples,

        error_message = NA_character_
      )
    },

    error = function(
      error
    ) {

      list(

        success = FALSE,

        samples = NULL,

        error_message = conditionMessage(
          error
        )
      )
    }
  )


  # ----------------------------------------------------------
  # F. Return immediately when JAGS itself failed
  # ----------------------------------------------------------

  if (
    !fit_result$success
  ) {

    return(

      list(

        success = FALSE,

        convergence_passed = FALSE,

        max_rhat = NA_real_,

        min_ess = NA_real_,

        error_message = fit_result$error_message,

        global_results = NULL,

        country_results = NULL,

        country_metrics = NULL
      )
    )
  }


  samples <- fit_result$samples


  # ----------------------------------------------------------
  # G. Convergence diagnostics for the main global parameters
  # ----------------------------------------------------------

  main_parameters <- c(
    "alpha",
    "beta",
    "r",
    "sigma_country"
  )


  selected_samples <- samples[
    ,
    main_parameters
  ]


  diagnostic_result <- tryCatch(

    {

      gelman_result <- gelman.diag(

        selected_samples,

        autoburnin = FALSE,

        multivariate = FALSE
      )


      rhat_values <- gelman_result$psrf[
        ,
        "Point est."
      ]


      ess_values <- effectiveSize(
        selected_samples
      )


      list(

        success = TRUE,

        max_rhat = max(
          rhat_values,
          na.rm = TRUE
        ),

        min_ess = min(
          ess_values,
          na.rm = TRUE
        )
      )
    },

    error = function(
      error
    ) {

      list(

        success = FALSE,

        max_rhat = NA_real_,

        min_ess = NA_real_
      )
    }
  )


  convergence_passed <- (

    diagnostic_result$success

    &&

    is.finite(
      diagnostic_result$max_rhat
    )

    &&

    is.finite(
      diagnostic_result$min_ess
    )

    &&

    diagnostic_result$max_rhat
    <=
    R_HAT_THRESHOLD

    &&

    diagnostic_result$min_ess
    >=
    ESS_THRESHOLD
  )


  # ----------------------------------------------------------
  # H. Convert posterior to matrix
  # ----------------------------------------------------------

  posterior <- as.matrix(
    samples
  )


  # ----------------------------------------------------------
  # I. Summarize global parameter recovery
  # ----------------------------------------------------------

  global_rows <- list()


  for (
    parameter
    in main_parameters
  ) {

    draws <- posterior[
      ,
      parameter
    ]


    true_value <- TRUE_PARAMETER_VALUES[
      parameter
    ]


    lower_95 <- unname(

      quantile(
        draws,
        0.025
      )
    )


    upper_95 <- unname(

      quantile(
        draws,
        0.975
      )
    )


    global_rows[[
        parameter
      ]] <- data.frame(

      attempt_id = attempt_id,

      parameter = parameter,

      true_value = as.numeric(
        true_value
      ),

      posterior_mean = mean(
        draws
      ),

      posterior_median = median(
        draws
      ),

      lower_95 = lower_95,

      upper_95 = upper_95,

      interval_width = (
        upper_95
        -
        lower_95
      ),

      covered = (

        true_value
        >=
        lower_95

        &&

        true_value
        <=
        upper_95
      ),

      convergence_passed = convergence_passed
    )
  }


  global_summary <- do.call(

    rbind,

    global_rows
  )


  rownames(
    global_summary
  ) <- NULL


  # ----------------------------------------------------------
  # J. Summarize country-effect recovery
  # ----------------------------------------------------------

  u_parameters <- find_ordered_u_parameters(

    colnames(
      posterior
    )
  )


  stopifnot(
    length(
      u_parameters
    ) == N_COUNTRY
  )


  country_rows <- list()


  for (
    j
    in seq_len(
      N_COUNTRY
    )
  ) {

    draws <- posterior[
      ,
      u_parameters[
        j
      ]
    ]


    lower_95 <- unname(

      quantile(
        draws,
        0.025
      )
    )


    upper_95 <- unname(

      quantile(
        draws,
        0.975
      )
    )


    country_rows[[
        j
      ]] <- data.frame(

      attempt_id = attempt_id,

      country_id = j,

      true_u = u_true[
        j
      ],

      posterior_mean_u = mean(
        draws
      ),

      lower_95 = lower_95,

      upper_95 = upper_95,

      interval_width = (
        upper_95
        -
        lower_95
      ),

      covered = (

        u_true[
          j
        ]
        >=
        lower_95

        &&

        u_true[
          j
        ]
        <=
        upper_95
      ),

      convergence_passed = convergence_passed
    )
  }


  country_summary <- do.call(

    rbind,

    country_rows
  )


  country_effect_metrics <- data.frame(

    attempt_id = attempt_id,

    correlation = cor(

      country_summary$true_u,

      country_summary$posterior_mean_u
    ),

    rmse = sqrt(

      mean(

        (

          country_summary$posterior_mean_u

          -

          country_summary$true_u

        )^2
      )
    ),

    coverage = mean(
      country_summary$covered
    ),

    average_interval_width = mean(
      country_summary$interval_width
    ),

    convergence_passed = convergence_passed
  )


  return(

    list(

      success = TRUE,

      convergence_passed = convergence_passed,

      max_rhat = diagnostic_result$max_rhat,

      min_ess = diagnostic_result$min_ess,

      error_message = NA_character_,

      global_results = global_summary,

      country_results = country_summary,

      country_metrics = country_effect_metrics
    )
  )
}


# ============================================================
# 14. Run repeated recovery attempts
# ============================================================

attempt_id <- next_attempt_id


while (

  current_successes
  <
  N_TARGET_SUCCESS

  &&

  attempt_id
  <=
  MAX_ATTEMPTS

) {

  attempt_start_time <- Sys.time()


  result <- run_one_recovery_replication(

    attempt_id = attempt_id
  )


  attempt_end_time <- Sys.time()


  elapsed_seconds <- as.numeric(

    difftime(

      attempt_end_time,

      attempt_start_time,

      units = "secs"
    )
  )


  replication_success <- (

    result$success

    &&

    result$convergence_passed
  )


  attempt_row <- data.frame(

    attempt_id = attempt_id,

    success = replication_success,

    model_fit_completed = result$success,

    convergence_passed = result$convergence_passed,

    max_rhat = result$max_rhat,

    min_ess = result$min_ess,

    elapsed_seconds = elapsed_seconds,

    error_message = ifelse(

      is.na(
        result$error_message
      ),

      "",

      result$error_message
    )
  )


  attempt_log <- append_and_save(

    attempt_log,

    attempt_row,

    ATTEMPT_LOG_PATH
  )


  if (
    result$success
  ) {

    global_results <- append_and_save(

      global_results,

      result$global_results,

      GLOBAL_RESULTS_PATH
    )


    country_results <- append_and_save(

      country_results,

      result$country_results,

      COUNTRY_RESULTS_PATH
    )


    country_metrics <- append_and_save(

      country_metrics,

      result$country_metrics,

      COUNTRY_METRICS_PATH
    )
  }


  if (
    replication_success
  ) {

    current_successes <- (
      current_successes
      +
      1
    )
  }


  cat(
    "\nAttempt:",
    attempt_id,
    "\n"
  )

  cat(
    "Accepted successful replications:",
    current_successes,
    "/",
    N_TARGET_SUCCESS,
    "\n"
  )

  cat(
    "Maximum R-hat:",
    result$max_rhat,
    "\n"
  )

  cat(
    "Minimum ESS:",
    result$min_ess,
    "\n"
  )

  cat(
    "Elapsed seconds:",
    round(
      elapsed_seconds,
      1
    ),
    "\n"
  )


  attempt_id <- (
    attempt_id
    +
    1
  )
}


# ============================================================
# 15. Check whether the target was reached
# ============================================================

if (
  current_successes
  <
  N_TARGET_SUCCESS
) {

  warning(

    paste(

      "The recovery study stopped with",

      current_successes,

      "successful replications out of the target of",

      N_TARGET_SUCCESS,

      "after reaching MAX_ATTEMPTS."
    )
  )
}


# ============================================================
# 16. Reload all saved progress before final summaries
# ============================================================

attempt_log <- read.csv(
  ATTEMPT_LOG_PATH
)

global_results <- read.csv(
  GLOBAL_RESULTS_PATH
)

country_results <- read.csv(
  COUNTRY_RESULTS_PATH
)

country_metrics <- read.csv(
  COUNTRY_METRICS_PATH
)


# ============================================================
# 17. Keep only accepted converged replications
# ============================================================

successful_attempt_ids <- attempt_log$attempt_id[
  attempt_log$success
]


successful_global_results <- global_results[
  global_results$attempt_id
  %in%
  successful_attempt_ids,
]


successful_country_results <- country_results[
  country_results$attempt_id
  %in%
  successful_attempt_ids,
]


successful_country_metrics <- country_metrics[
  country_metrics$attempt_id
  %in%
  successful_attempt_ids,
]


# ============================================================
# 18. Summarize global parameter recovery
# ============================================================

recovery_summary_rows <- list()


for (
  parameter
  in names(
    TRUE_PARAMETER_VALUES
  )
) {

  parameter_data <- successful_global_results[

    successful_global_results$parameter
    ==
    parameter,

  ]


  true_value <- TRUE_PARAMETER_VALUES[
    parameter
  ]


  mean_estimate <- mean(
    parameter_data$posterior_mean
  )


  recovery_summary_rows[[
      parameter
    ]] <- data.frame(

    parameter = parameter,

    true_value = as.numeric(
      true_value
    ),

    mean_posterior_estimate = mean_estimate,

    bias = (
      mean_estimate
      -
      true_value
    ),

    rmse = sqrt(

      mean(

        (

          parameter_data$posterior_mean

          -

          true_value

        )^2
      )
    ),

    coverage_95 = mean(
      parameter_data$covered
    ),

    average_interval_width = mean(
      parameter_data$interval_width
    ),

    successful_replications = nrow(
      parameter_data
    )
  )
}


parameter_recovery_summary <- do.call(

  rbind,

  recovery_summary_rows
)


rownames(
  parameter_recovery_summary
) <- NULL


# ============================================================
# 19. Overall convergence-failure rate
# ============================================================

convergence_failure_rate <- (

  1

  -

  mean(
    attempt_log$success
  )
)


recovery_study_status <- data.frame(

  attempted_replications = nrow(
    attempt_log
  ),

  successful_replications = sum(
    attempt_log$success
  ),

  convergence_failure_rate = convergence_failure_rate,

  mean_runtime_seconds = mean(
    attempt_log$elapsed_seconds
  ),

  median_runtime_seconds = median(
    attempt_log$elapsed_seconds
  )
)


# ============================================================
# 20. Overall country-effect recovery summary
# ============================================================

country_effect_recovery_summary <- data.frame(

  mean_correlation = mean(
    successful_country_metrics$correlation
  ),

  median_correlation = median(
    successful_country_metrics$correlation
  ),

  mean_rmse = mean(
    successful_country_metrics$rmse
  ),

  mean_coverage_95 = mean(
    successful_country_metrics$coverage
  ),

  average_interval_width = mean(
    successful_country_metrics$average_interval_width
  ),

  successful_replications = nrow(
    successful_country_metrics
  )
)


# ============================================================
# 21. Save final summary tables
# ============================================================

write.csv(

  parameter_recovery_summary,

  file.path(
    TABLE_DIR,
    "parameter_recovery_summary.csv"
  ),

  row.names = FALSE
)


write.csv(

  recovery_study_status,

  file.path(
    TABLE_DIR,
    "parameter_recovery_study_status.csv"
  ),

  row.names = FALSE
)


write.csv(

  country_effect_recovery_summary,

  file.path(
    TABLE_DIR,
    "country_effect_recovery_summary.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 22. Recovery interval figure
# ============================================================

recovery_plot_data <- successful_global_results


recovery_plot_data$attempt_id <- factor(
  recovery_plot_data$attempt_id
)


parameter_recovery_plot <- ggplot(

  recovery_plot_data,

  aes(
    x = attempt_id,
    y = posterior_mean
  )

) +

  geom_errorbar(

    aes(
      ymin = lower_95,
      ymax = upper_95
    ),

    width = 0.2
  ) +

  geom_point(
    size = 1.8
  ) +

  geom_hline(

    aes(
      yintercept = true_value
    ),

    linetype = "dashed"
  ) +

  facet_wrap(

    ~ parameter,

    scales = "free_y"
  ) +

  labs(

    title = "Repeated Parameter-Recovery Study",

    subtitle = (
      "Posterior means and 95% credible intervals; dashed lines show the known true values"
    ),

    x = "Successful simulation replication",

    y = "Parameter value"
  ) +

  theme_minimal(
    base_size = 11
  )


print(
  parameter_recovery_plot
)


ggsave(

  filename = file.path(
    FIGURE_DIR,
    "parameter_recovery_intervals.png"
  ),

  plot = parameter_recovery_plot,

  width = 10,

  height = 8,

  dpi = 300
)


# ============================================================
# 23. Print final results
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "FORMAL PARAMETER-RECOVERY SUMMARY\n"
)

cat(
  "============================================================\n\n"
)


print(
  parameter_recovery_summary,
  row.names = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "COUNTRY-EFFECT RECOVERY SUMMARY\n"
)

cat(
  "============================================================\n\n"
)


print(
  country_effect_recovery_summary,
  row.names = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "RECOVERY-STUDY STATUS\n"
)

cat(
  "============================================================\n\n"
)


print(
  recovery_study_status,
  row.names = FALSE
)


# ============================================================
# 24. Save session information
# ============================================================

capture.output(

  sessionInfo(),

  file = file.path(
    RESULT_DIR,
    "R_session_info.txt"
  )
)

