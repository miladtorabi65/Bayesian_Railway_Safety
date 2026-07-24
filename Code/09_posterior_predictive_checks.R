# ============================================================
# POSTERIOR PREDICTIVE CHECKS
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

library(ggplot2)


# ============================================================
# 2. Paths
# ============================================================

DATA_PATH <- paste0(
  "../data/processed/",
  "railway_country_year_2010_2024.csv"
)

RESULT_DIR <- "../results/mcmc"

TABLE_DIR <- "../tables/posterior_predictive_checks"

FIGURE_DIR <- "../figures/posterior_predictive_checks"


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
# 3. Reproducibility
# ============================================================

set.seed(2026)


# ============================================================
# 4. Load and validate data
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
  !anyNA(data)
)

stopifnot(
  all(
    data$train_km > 0
  )
)


y <- as.integer(
  data$accidents
)

exposure <- as.numeric(
  data$train_km
)

year_centered <- as.numeric(
  data$year_centered
)

country_id <- as.integer(
  data$country_id
)

years <- sort(
  unique(
    data$year
  )
)


# ============================================================
# 5. Load posterior samples
# ============================================================

posterior_M1 <- as.matrix(
  readRDS(
    file.path(
      RESULT_DIR,
      "M1_Poisson",
      "posterior_samples.rds"
    )
  )
)


posterior_M2 <- as.matrix(
  readRDS(
    file.path(
      RESULT_DIR,
      "M2_Negative_Binomial",
      "posterior_samples.rds"
    )
  )
)


posterior_M3 <- as.matrix(
  readRDS(
    file.path(
      RESULT_DIR,
      "M3_Hierarchical_Negative_Binomial",
      "posterior_samples_full.rds"
    )
  )
)


cat(
  "Posterior samples loaded successfully.\n"
)


# ============================================================
# 6. Number of posterior predictive simulations
# ============================================================

N_PPC <- 2000


# ============================================================
# 7. Find and order M3 country effects
# ============================================================

u_parameters <- colnames(
  posterior_M3
)


u_parameters <- u_parameters[
  startsWith(
    u_parameters,
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


u_parameters <- u_parameters[
  order(
    u_numbers
  )
]


stopifnot(
  length(
    u_parameters
  ) == 27
)


# ============================================================
# 8. Posterior predictive simulation function
# ============================================================

simulate_posterior_predictive <- function(
    posterior,
    model_name,
    n_simulations,
    exposure,
    year_centered,
    country_id,
    u_parameters = NULL
) {

  selected_draws <- sample(
    seq_len(
      nrow(
        posterior
      )
    ),
    size = n_simulations,
    replace = FALSE
  )


  N <- length(
    exposure
  )


  y_rep <- matrix(
    NA_integer_,
    nrow = n_simulations,
    ncol = N
  )


  for (
    s
    in seq_len(
      n_simulations
    )
  ) {

    draw <- selected_draws[
      s
    ]


    alpha <- posterior[
      draw,
      "alpha"
    ]


    beta <- posterior[
      draw,
      "beta"
    ]


    eta <- (
      alpha
      +
      beta
      *
      year_centered
    )


    # Add country effects only for M3
    if (
      !is.null(
        u_parameters
      )
    ) {

      u <- posterior[
        draw,
        u_parameters
      ]


      eta <- (
        eta
        +
        u[
          country_id
        ]
      )
    }


    mu <- (
      exposure
      *
      exp(
        eta
      )
    )


    if (
      model_name
      ==
      "M1_Poisson"
    ) {

      y_rep[
        s,
      ] <- rpois(
        N,
        lambda = mu
      )

    } else {

      r <- posterior[
        draw,
        "r"
      ]


      y_rep[
        s,
      ] <- rnbinom(
        N,
        size = r,
        mu = mu
      )
    }
  }


  return(
    y_rep
  )
}


# ============================================================
# 9. Generate posterior predictive datasets
# ============================================================

cat(
  "Simulating M1 posterior predictive datasets...\n"
)

y_rep_M1 <- simulate_posterior_predictive(

  posterior = posterior_M1,

  model_name = "M1_Poisson",

  n_simulations = N_PPC,

  exposure = exposure,

  year_centered = year_centered,

  country_id = country_id
)


cat(
  "Simulating M2 posterior predictive datasets...\n"
)

y_rep_M2 <- simulate_posterior_predictive(

  posterior = posterior_M2,

  model_name = "M2_Negative_Binomial",

  n_simulations = N_PPC,

  exposure = exposure,

  year_centered = year_centered,

  country_id = country_id
)


cat(
  "Simulating M3 posterior predictive datasets...\n"
)

y_rep_M3 <- simulate_posterior_predictive(

  posterior = posterior_M3,

  model_name = "M3_Hierarchical_Negative_Binomial",

  n_simulations = N_PPC,

  exposure = exposure,

  year_centered = year_centered,

  country_id = country_id,

  u_parameters = u_parameters
)


# ============================================================
# 10. Calculate dataset-level statistics
# ============================================================

calculate_ppc_statistics <- function(
    y_values,
    exposure,
    country_id
) {

  country_accidents <- tapply(
    y_values,
    country_id,
    sum
  )


  country_exposure <- tapply(
    exposure,
    country_id,
    sum
  )


  country_rates <- (
    country_accidents
    /
    country_exposure
  )


  c(

    total = sum(
      y_values
    ),

    mean = mean(
      y_values
    ),

    variance = var(
      y_values
    ),

    zeros = sum(
      y_values == 0
    ),

    maximum = max(
      y_values
    ),

    country_rate_variance = var(
      country_rates
    )
  )
}


observed_statistics <- calculate_ppc_statistics(

  y_values = y,

  exposure = exposure,

  country_id = country_id
)


cat(
  "\nObserved statistics:\n"
)

print(
  observed_statistics
)


# ============================================================
# 11. Summarize PPC statistics
# ============================================================

summarize_ppc <- function(
    y_rep,
    model_name
) {

  replicated_statistics <- t(
    apply(
      y_rep,
      1,
      calculate_ppc_statistics,
      exposure = exposure,
      country_id = country_id
    )
  )


  results <- data.frame(

    model = model_name,

    statistic = colnames(
      replicated_statistics
    ),

    observed = as.numeric(
      observed_statistics
    ),

    predictive_median = apply(
      replicated_statistics,
      2,
      median
    ),

    lower_95 = apply(
      replicated_statistics,
      2,
      quantile,
      probs = 0.025
    ),

    upper_95 = apply(
      replicated_statistics,
      2,
      quantile,
      probs = 0.975
    )
  )


  results$observed_inside_95 <- (

    results$observed
    >=
    results$lower_95

    &

    results$observed
    <=
    results$upper_95
  )


  return(
    list(

      summary = results,

      raw = replicated_statistics
    )
  )
}


ppc_M1 <- summarize_ppc(
  y_rep_M1,
  "M1_Poisson"
)


ppc_M2 <- summarize_ppc(
  y_rep_M2,
  "M2_Negative_Binomial"
)


ppc_M3 <- summarize_ppc(
  y_rep_M3,
  "M3_Hierarchical_Negative_Binomial"
)


ppc_summary <- rbind(

  ppc_M1$summary,

  ppc_M2$summary,

  ppc_M3$summary
)


cat(
  "\nPosterior predictive summary:\n"
)

print(
  ppc_summary,
  row.names = FALSE
)


write.csv(

  ppc_summary,

  file.path(
    TABLE_DIR,
    "posterior_predictive_summary.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 12. Yearly aggregate-rate PPC for M3
# ============================================================

calculate_yearly_rates <- function(
    y_values
) {

  yearly_accidents <- tapply(
    y_values,
    data$year,
    sum
  )


  yearly_exposure <- tapply(
    data$train_km,
    data$year,
    sum
  )


  yearly_accidents / yearly_exposure
}


observed_yearly_rates <- calculate_yearly_rates(
  y
)


m3_yearly_rates <- t(
  apply(
    y_rep_M3,
    1,
    calculate_yearly_rates
  )
)


yearly_ppc_data <- data.frame(

  year = years,

  observed = as.numeric(
    observed_yearly_rates
  ),

  predictive_median = apply(
    m3_yearly_rates,
    2,
    median
  ),

  lower_95 = apply(
    m3_yearly_rates,
    2,
    quantile,
    probs = 0.025
  ),

  upper_95 = apply(
    m3_yearly_rates,
    2,
    quantile,
    probs = 0.975
  )
)


write.csv(

  yearly_ppc_data,

  file.path(
    TABLE_DIR,
    "m3_yearly_rate_ppc.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 13. Yearly PPC figure
# ============================================================

yearly_ppc_plot <- ggplot(

  yearly_ppc_data,

  aes(
    x = year
  )

) +

  geom_ribbon(

    aes(
      ymin = lower_95,
      ymax = upper_95
    ),

    alpha = 0.25
  ) +

  geom_line(

    aes(
      y = predictive_median
    ),

    linewidth = 1
  ) +

  geom_point(

    aes(
      y = observed
    ),

    size = 2
  ) +

  geom_line(

    aes(
      y = observed
    ),

    linetype = "dashed"
  ) +

  labs(

    title = (
      "Posterior Predictive Check of the Aggregate Yearly Accident Rate"
    ),

    subtitle = (
      "Observed rates compared with the M3 posterior predictive distribution"
    ),

    x = "Year",

    y = "Significant accidents per million train-km"
  ) +

  theme_minimal(
    base_size = 12
  )


print(
  yearly_ppc_plot
)


ggsave(

  filename = file.path(
    FIGURE_DIR,
    "m3_yearly_rate_posterior_predictive_check.png"
  ),

  plot = yearly_ppc_plot,

  width = 10,

  height = 6,

  dpi = 300
)


# ============================================================
# 14. Final status
# ============================================================

cat(
  "\nPosterior predictive checks completed successfully.\n"
)

cat(
  "Saved table:\n"
)

cat(
  file.path(
    TABLE_DIR,
    "posterior_predictive_summary.csv"
  ),
  "\n"
)

cat(
  file.path(
    TABLE_DIR,
    "m3_yearly_rate_ppc.csv"
  ),
  "\n"
)

cat(
  "Saved figure:\n"
)

cat(
  file.path(
    FIGURE_DIR,
    "m3_yearly_rate_posterior_predictive_check.png"
  ),
  "\n"
)
