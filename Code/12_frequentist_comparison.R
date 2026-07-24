# ============================================================
# FREQUENTIST COMPARISON
#
# Frequentist counterparts:
#   F1: Poisson GLM
#   F2: Negative Binomial regression
#   F3: Negative Binomial mixed-effects model
#
# These models use the same outcome, exposure offset, time
# covariate, and country grouping as the Bayesian models.
# ============================================================


# ============================================================
# 1. Required packages
# ============================================================

required_packages <- c(
  "MASS",
  "glmmTMB",
  "ggplot2"
)


missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]


if (
  length(
    missing_packages
  ) > 0
) {

  stop(
    paste0(
      "Missing required package(s): ",
      paste(
        missing_packages,
        collapse = ", "
      ),
      ". Install them first with: install.packages(c(",
      paste0(
        '"',
        missing_packages,
        '"',
        collapse = ", "
      ),
      "))"
    )
  )
}


library(MASS)
library(glmmTMB)
library(ggplot2)


# ============================================================
# 2. Paths
# ============================================================

DATA_PATH <- paste0(
  "../data/processed/",
  "railway_country_year_2010_2024.csv"
)

BAYESIAN_SUMMARY_PATH <- paste0(
  "../tables/posterior_inference/",
  "posterior_main_summary.csv"
)

TABLE_DIR <- "../tables/frequentist_comparison"

FIGURE_DIR <- "../figures/frequentist_comparison"

RESULT_DIR <- "../results/frequentist_comparison"


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

dir.create(
  RESULT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. Load and validate data
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

stopifnot(
  all(
    data$accidents >= 0
  )
)


# The grouping variable must be a factor for the mixed model.

data$country_id <- factor(
  data$country_id
)


cat(
  "Dataset validation passed.\n"
)


# ============================================================
# 4. Fit F1 — Poisson GLM
# ============================================================

cat(
  "\nFitting frequentist Poisson GLM...\n"
)


fit_F1 <- glm(

  accidents
  ~
  year_centered
  +
  offset(
    log(
      train_km
    )
  ),

  family = poisson(
    link = "log"
  ),

  data = data
)


# ============================================================
# 5. Fit F2 — Negative Binomial regression
# ============================================================

cat(
  "Fitting frequentist Negative Binomial regression...\n"
)


fit_F2 <- MASS::glm.nb(

  accidents
  ~
  year_centered
  +
  offset(
    log(
      train_km
    )
  ),

  data = data,

  link = log
)


# ============================================================
# 6. Fit F3 — Hierarchical Negative Binomial model
#
# nbinom2 uses:
#
# Var(Y) = mu + mu^2 / dispersion
#
# which is directly comparable to the Bayesian M3 parameter r.
# ============================================================

cat(
  "Fitting frequentist Negative Binomial mixed model...\n"
)


fit_F3 <- glmmTMB::glmmTMB(

  accidents
  ~
  year_centered
  +
  offset(
    log(
      train_km
    )
  )
  +
  (
    1
    |
    country_id
  ),

  family = glmmTMB::nbinom2(
    link = "log"
  ),

  data = data
)


cat(
  "\nAll frequentist models fitted.\n"
)


# ============================================================
# 7. Check mixed-model optimization
# ============================================================

mixed_model_convergence <- data.frame(

  optimizer_convergence_code = fit_F3$fit$convergence,

  positive_definite_hessian = fit_F3$sdr$pdHess
)


write.csv(

  mixed_model_convergence,

  file.path(
    TABLE_DIR,
    "frequentist_mixed_model_convergence.csv"
  ),

  row.names = FALSE
)


cat(
  "\nMixed-model convergence check:\n"
)


print(
  mixed_model_convergence,
  row.names = FALSE
)


# ============================================================
# 8. Helper: summarize the temporal coefficient
# ============================================================

summarize_time_effect <- function(
    estimate,
    standard_error,
    model_name
) {

  lower_beta <- (
    estimate
    -
    1.96
    *
    standard_error
  )


  upper_beta <- (
    estimate
    +
    1.96
    *
    standard_error
  )


  annual_rate_ratio <- exp(
    estimate
  )


  lower_rate_ratio <- exp(
    lower_beta
  )


  upper_rate_ratio <- exp(
    upper_beta
  )


  annual_percent_change <- (
    100
    *
    (
      annual_rate_ratio
      -
      1
    )
  )


  lower_percent_change <- (
    100
    *
    (
      lower_rate_ratio
      -
      1
    )
  )


  upper_percent_change <- (
    100
    *
    (
      upper_rate_ratio
      -
      1
    )
  )


  data.frame(

    model = model_name,

    beta = estimate,

    beta_se = standard_error,

    beta_lower_95 = lower_beta,

    beta_upper_95 = upper_beta,

    annual_rate_ratio = annual_rate_ratio,

    annual_rate_ratio_lower_95 = lower_rate_ratio,

    annual_rate_ratio_upper_95 = upper_rate_ratio,

    annual_percent_change = annual_percent_change,

    annual_percent_change_lower_95 = lower_percent_change,

    annual_percent_change_upper_95 = upper_percent_change
  )
}


# ============================================================
# 9. Extract the time effect from F1
# ============================================================

F1_coefficient_table <- summary(
  fit_F1
)$coefficients


F1_time_summary <- summarize_time_effect(

  estimate = F1_coefficient_table[
    "year_centered",
    "Estimate"
  ],

  standard_error = F1_coefficient_table[
    "year_centered",
    "Std. Error"
  ],

  model_name = "F1_Poisson_GLM"
)


# ============================================================
# 10. Extract the time effect from F2
# ============================================================

F2_coefficient_table <- summary(
  fit_F2
)$coefficients


F2_time_summary <- summarize_time_effect(

  estimate = F2_coefficient_table[
    "year_centered",
    "Estimate"
  ],

  standard_error = F2_coefficient_table[
    "year_centered",
    "Std. Error"
  ],

  model_name = "F2_Negative_Binomial"
)


# ============================================================
# 11. Extract the time effect from F3
# ============================================================

F3_coefficient_table <- summary(
  fit_F3
)$coefficients$cond


F3_time_summary <- summarize_time_effect(

  estimate = F3_coefficient_table[
    "year_centered",
    "Estimate"
  ],

  standard_error = F3_coefficient_table[
    "year_centered",
    "Std. Error"
  ],

  model_name = "F3_Hierarchical_Negative_Binomial"
)


# ============================================================
# 12. Combine frequentist temporal results
# ============================================================

frequentist_time_summary <- rbind(

  F1_time_summary,

  F2_time_summary,

  F3_time_summary
)


write.csv(

  frequentist_time_summary,

  file.path(
    TABLE_DIR,
    "frequentist_time_trend_summary.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 13. Extract dispersion parameters
#
# For MASS::glm.nb:
#   Var(Y) = mu + mu^2 / theta
#
# For glmmTMB::nbinom2:
#   Var(Y) = mu + mu^2 / dispersion
# ============================================================

dispersion_summary <- data.frame(

  model = c(
    "F2_Negative_Binomial",
    "F3_Hierarchical_Negative_Binomial"
  ),

  parameter = c(
    "theta",
    "dispersion"
  ),

  estimate = c(

    fit_F2$theta,

    sigma(
      fit_F3
    )
  )
)


write.csv(

  dispersion_summary,

  file.path(
    TABLE_DIR,
    "frequentist_dispersion_summary.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 14. Extract the country random-effect SD from F3
# ============================================================

country_variance_matrix <- VarCorr(
  fit_F3
)[[
  c(
    "cond",
    "country_id"
  )
]]


country_sd <- as.numeric(
  attr(
    country_variance_matrix,
    "stddev"
  )[
    1
  ]
)


country_heterogeneity_summary <- data.frame(

  model = "F3_Hierarchical_Negative_Binomial",

  parameter = "sigma_country",

  estimate = country_sd,

  one_sd_rate_multiplier = exp(
    country_sd
  )
)


write.csv(

  country_heterogeneity_summary,

  file.path(
    TABLE_DIR,
    "frequentist_country_heterogeneity_summary.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 15. Frequentist model fit statistics
#
# These are descriptive comparisons among the frequentist
# models. They are not direct replacements for the Bayesian
# DIC analysis.
# ============================================================

fit_statistics <- data.frame(

  model = c(
    "F1_Poisson_GLM",
    "F2_Negative_Binomial",
    "F3_Hierarchical_Negative_Binomial"
  ),

  log_likelihood = c(

    as.numeric(
      logLik(
        fit_F1
      )
    ),

    as.numeric(
      logLik(
        fit_F2
      )
    ),

    as.numeric(
      logLik(
        fit_F3
      )
    )
  ),

  AIC = c(

    AIC(
      fit_F1
    ),

    AIC(
      fit_F2
    ),

    AIC(
      fit_F3
    )
  ),

  BIC = c(

    BIC(
      fit_F1
    ),

    BIC(
      fit_F2
    ),

    BIC(
      fit_F3
    )
  )
)


fit_statistics$delta_AIC <- (

  fit_statistics$AIC

  -

  min(
    fit_statistics$AIC
  )
)


fit_statistics$delta_BIC <- (

  fit_statistics$BIC

  -

  min(
    fit_statistics$BIC
  )
)


fit_statistics <- fit_statistics[
  order(
    fit_statistics$AIC
  ),
]


rownames(
  fit_statistics
) <- NULL


write.csv(

  fit_statistics,

  file.path(
    TABLE_DIR,
    "frequentist_model_fit_statistics.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 16. Load Bayesian posterior summaries
# ============================================================

bayesian_summary <- read.csv(

  BAYESIAN_SUMMARY_PATH
)


bayesian_annual_change <- bayesian_summary[

  bayesian_summary$parameter
  ==
  "annual_percent_change",

]


bayesian_comparison <- data.frame(

  approach = "Bayesian",

  model = bayesian_annual_change$model,

  estimate = bayesian_annual_change$median,

  lower_95 = bayesian_annual_change$lower_95,

  upper_95 = bayesian_annual_change$upper_95
)


frequentist_comparison <- data.frame(

  approach = "Frequentist",

  model = c(
    "M1_Poisson",
    "M2_Negative_Binomial",
    "M3_Hierarchical_Negative_Binomial"
  ),

  estimate = frequentist_time_summary$annual_percent_change,

  lower_95 = frequentist_time_summary$annual_percent_change_lower_95,

  upper_95 = frequentist_time_summary$annual_percent_change_upper_95
)


bayesian_vs_frequentist <- rbind(

  bayesian_comparison,

  frequentist_comparison
)


model_order <- c(

  "M1_Poisson",

  "M2_Negative_Binomial",

  "M3_Hierarchical_Negative_Binomial"
)


bayesian_vs_frequentist$model <- factor(

  bayesian_vs_frequentist$model,

  levels = model_order
)


write.csv(

  bayesian_vs_frequentist,

  file.path(
    TABLE_DIR,
    "bayesian_vs_frequentist_annual_trend.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 17. Compare key M2 and M3 structural parameters
# ============================================================

bayesian_r <- bayesian_summary[

  bayesian_summary$parameter
  ==
  "r",

]


bayesian_sigma_country <- bayesian_summary[

  bayesian_summary$model
  ==
  "M3_Hierarchical_Negative_Binomial"

  &

  bayesian_summary$parameter
  ==
  "sigma_country",

]


structural_comparison <- data.frame(

  quantity = c(
    "M2_dispersion",
    "M3_dispersion",
    "M3_sigma_country"
  ),

  bayesian_estimate = c(

    bayesian_r$median[
      bayesian_r$model
      ==
      "M2_Negative_Binomial"
    ],

    bayesian_r$median[
      bayesian_r$model
      ==
      "M3_Hierarchical_Negative_Binomial"
    ],

    bayesian_sigma_country$median
  ),

  frequentist_estimate = c(

    fit_F2$theta,

    sigma(
      fit_F3
    ),

    country_sd
  )
)


write.csv(

  structural_comparison,

  file.path(
    TABLE_DIR,
    "bayesian_vs_frequentist_structural_parameters.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 18. Annual-trend comparison figure
# ============================================================

comparison_plot <- ggplot(

  bayesian_vs_frequentist,

  aes(
    x = model,
    y = estimate,
    shape = approach,
    group = approach
  )

) +

  geom_errorbar(

    aes(
      ymin = lower_95,
      ymax = upper_95
    ),

    position = position_dodge(
      width = 0.35
    ),

    width = 0.12
  ) +

  geom_point(

    position = position_dodge(
      width = 0.35
    ),

    size = 3
  ) +

  geom_hline(

    yintercept = 0,

    linetype = "dashed"
  ) +

  labs(

    title = "Bayesian and Frequentist Estimates of the Annual Railway Accident Trend",

    subtitle = "Points show annual percentage changes with 95% uncertainty intervals",

    x = NULL,

    y = "Annual percentage change (%)",

    shape = "Approach"
  ) +

  theme_minimal(
    base_size = 12
  )


print(
  comparison_plot
)


ggsave(

  filename = file.path(
    FIGURE_DIR,
    "bayesian_vs_frequentist_annual_trend.png"
  ),

  plot = comparison_plot,

  width = 10,

  height = 6,

  dpi = 300
)


# ============================================================
# 19. Save textual model summaries
# ============================================================

capture.output(

  summary(
    fit_F1
  ),

  file = file.path(
    RESULT_DIR,
    "F1_Poisson_GLM_summary.txt"
  )
)


capture.output(

  summary(
    fit_F2
  ),

  file = file.path(
    RESULT_DIR,
    "F2_Negative_Binomial_summary.txt"
  )
)


capture.output(

  summary(
    fit_F3
  ),

  file = file.path(
    RESULT_DIR,
    "F3_Hierarchical_Negative_Binomial_summary.txt"
  )
)


capture.output(

  sessionInfo(),

  file = file.path(
    RESULT_DIR,
    "R_session_info.txt"
  )
)


# ============================================================
# 20. Print key results
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "FREQUENTIST TIME-TREND RESULTS\n"
)

cat(
  "============================================================\n\n"
)


print(
  frequentist_time_summary,
  row.names = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "FREQUENTIST MODEL FIT STATISTICS\n"
)

cat(
  "============================================================\n\n"
)


print(
  fit_statistics,
  row.names = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "BAYESIAN VS FREQUENTIST STRUCTURAL PARAMETERS\n"
)

cat(
  "============================================================\n\n"
)


print(
  structural_comparison,
  row.names = FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "MIXED-MODEL CONVERGENCE\n"
)

cat(
  "============================================================\n\n"
)


print(
  mixed_model_convergence,
  row.names = FALSE
)


cat(
  "\nFrequentist comparison completed successfully.\n"
)
