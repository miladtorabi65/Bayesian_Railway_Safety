# ============================================================
# MAIN BAYESIAN POSTERIOR INFERENCE
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

library(coda)
library(ggplot2)


# ============================================================
# 2. Paths
# ============================================================

RESULT_DIR <- "../results/mcmc"

TABLE_DIR <- "../tables/posterior_inference"

FIGURE_DIR <- "../figures/posterior_inference"


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
# 3. Load final posterior samples
# ============================================================

samples_M1 <- readRDS(
  file.path(
    RESULT_DIR,
    "M1_Poisson",
    "posterior_samples.rds"
  )
)


samples_M2 <- readRDS(
  file.path(
    RESULT_DIR,
    "M2_Negative_Binomial",
    "posterior_samples.rds"
  )
)


samples_M3 <- readRDS(
  file.path(
    RESULT_DIR,
    "M3_Hierarchical_Negative_Binomial",
    "posterior_samples_full.rds"
  )
)


country_lookup <- read.csv(
  file.path(
    RESULT_DIR,
    "country_lookup.csv"
  )
)


cat(
  "Posterior samples loaded successfully.\n"
)


# ============================================================
# 4. Posterior summary function
# ============================================================

summarize_parameter <- function(
    draws,
    model_name,
    parameter_name
) {
  
  data.frame(
    
    model = model_name,
    
    parameter = parameter_name,
    
    mean = mean(
      draws
    ),
    
    median = median(
      draws
    ),
    
    sd = sd(
      draws
    ),
    
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
# 5. Main model-summary function
# ============================================================

create_model_summary <- function(
    samples,
    model_name
) {
  
  posterior <- as.matrix(
    samples
  )
  
  
  alpha <- posterior[
    ,
    "alpha"
  ]
  
  
  beta <- posterior[
    ,
    "beta"
  ]
  
  
  # Baseline rate in the centered year: 2017
  baseline_rate <- exp(
    alpha
  )
  
  
  # Annual multiplicative change
  annual_rate_ratio <- exp(
    beta
  )
  
  
  # Annual percentage change
  annual_percent_change <- (
    100 *
      (
        annual_rate_ratio -
          1
      )
  )
  
  
  # ----------------------------------------------------------
  # Overall change from 2010 to 2024
  #
  # 2024 - 2010 = 14 years
  # ----------------------------------------------------------
  
  rate_ratio_2010_2024 <- exp(
    14 *
      beta
  )
  
  
  percent_change_2010_2024 <- (
    100 *
      (
        rate_ratio_2010_2024 -
          1
      )
  )
  
  
  # ----------------------------------------------------------
  # Average/reference rate in 2010 and 2024
  #
  # year_centered:
  #
  # 2010 = -7
  # 2017 =  0
  # 2024 = +7
  # ----------------------------------------------------------
  
  rate_2010 <- exp(
    alpha +
      beta *
      (-7)
  )
  
  
  rate_2024 <- exp(
    alpha +
      beta *
      7
  )
  
  
  results <- list(
    
    summarize_parameter(
      alpha,
      model_name,
      "alpha"
    ),
    
    summarize_parameter(
      beta,
      model_name,
      "beta"
    ),
    
    summarize_parameter(
      baseline_rate,
      model_name,
      "baseline_rate_2017"
    ),
    
    summarize_parameter(
      annual_rate_ratio,
      model_name,
      "annual_rate_ratio"
    ),
    
    summarize_parameter(
      annual_percent_change,
      model_name,
      "annual_percent_change"
    ),
    
    summarize_parameter(
      rate_ratio_2010_2024,
      model_name,
      "rate_ratio_2010_2024"
    ),
    
    summarize_parameter(
      percent_change_2010_2024,
      model_name,
      "percent_change_2010_2024"
    ),
    
    summarize_parameter(
      rate_2010,
      model_name,
      "rate_2010"
    ),
    
    summarize_parameter(
      rate_2024,
      model_name,
      "rate_2024"
    )
  )
  
  
  # Negative-binomial dispersion parameter
  if (
    "r" %in%
    colnames(
      posterior
    )
  ) {
    
    results[[
      length(
        results
      ) + 1
    ]] <- summarize_parameter(
      
      posterior[
        ,
        "r"
      ],
      
      model_name,
      
      "r"
    )
  }
  
  
  # Hierarchical country heterogeneity
  if (
    "sigma_country" %in%
    colnames(
      posterior
    )
  ) {
    
    sigma_country <- posterior[
      ,
      "sigma_country"
    ]
    
    
    results[[
      length(
        results
      ) + 1
    ]] <- summarize_parameter(
      
      sigma_country,
      
      model_name,
      
      "sigma_country"
    )
    
    
    # Multiplicative interpretation of one SD
    results[[
      length(
        results
      ) + 1
    ]] <- summarize_parameter(
      
      exp(
        sigma_country
      ),
      
      model_name,
      
      "country_sd_multiplier"
    )
  }
  
  
  do.call(
    rbind,
    results
  )
}





# ============================================================
# 6. Main posterior summaries
# ============================================================

posterior_summary_M1 <- create_model_summary(
  
  samples_M1,
  
  "M1_Poisson"
)


posterior_summary_M2 <- create_model_summary(
  
  samples_M2,
  
  "M2_Negative_Binomial"
)


posterior_summary_M3 <- create_model_summary(
  
  samples_M3,
  
  "M3_Hierarchical_Negative_Binomial"
)


posterior_summary_all <- rbind(
  
  posterior_summary_M1,
  
  posterior_summary_M2,
  
  posterior_summary_M3
)


print(
  posterior_summary_all
)


write.csv(
  
  posterior_summary_all,
  
  file.path(
    TABLE_DIR,
    "posterior_main_summary.csv"
  ),
  
  row.names = FALSE
)




# ============================================================
# 7. Posterior probability of a declining trend
# ============================================================

calculate_decline_probability <- function(
    samples,
    model_name
) {
  
  posterior <- as.matrix(
    samples
  )
  
  
  beta <- posterior[
    ,
    "beta"
  ]
  
  
  data.frame(
    
    model = model_name,
    
    probability_beta_below_zero = mean(
      beta
      <
        0
    ),
    
    probability_annual_rate_ratio_below_one = mean(
      exp(
        beta
      )
      <
        1
    ),
    
    probability_2010_2024_decline = mean(
      exp(
        14
        *
          beta
      )
      <
        1
    )
  )
}


probability_summary <- rbind(
  
  calculate_decline_probability(
    samples_M1,
    "M1_Poisson"
  ),
  
  calculate_decline_probability(
    samples_M2,
    "M2_Negative_Binomial"
  ),
  
  calculate_decline_probability(
    samples_M3,
    "M3_Hierarchical_Negative_Binomial"
  )
)


print(
  probability_summary
)

write.csv(
  
  probability_summary,
  
  file.path(
    TABLE_DIR,
    "posterior_probability_summary.csv"
  ),
  
  row.names = FALSE
)



# ============================================================
# 8. Country-specific posterior inference
# ============================================================

posterior_M3 <- as.matrix(
  samples_M3
)


parameter_names <- colnames(
  posterior_M3
)


u_parameters <- parameter_names[
  startsWith(
    parameter_names,
    "u["
  )
]


rr_parameters <- parameter_names[
  startsWith(
    parameter_names,
    "country_rate_ratio["
  )
]


stopifnot(
  length(
    u_parameters
  ) == 27
)


stopifnot(
  length(
    rr_parameters
  ) == 27
)





# ============================================================
# 9. Summarize country effects
# ============================================================

country_results <- list()


for (
  j
  in seq_len(
    27
  )
) {
  
  u_draws <- posterior_M3[
    ,
    u_parameters[
      j
    ]
  ]
  
  
  rr_draws <- posterior_M3[
    ,
    rr_parameters[
      j
    ]
  ]
  
  
  country_results[[
      j
    ]] <- data.frame(
    
    country_id = j,
    
    u_mean = mean(
      u_draws
    ),
    
    u_median = median(
      u_draws
    ),
    
    u_lower_95 = unname(
      quantile(
        u_draws,
        0.025
      )
    ),
    
    u_upper_95 = unname(
      quantile(
        u_draws,
        0.975
      )
    ),
    
    rate_ratio_mean = mean(
      rr_draws
    ),
    
    rate_ratio_median = median(
      rr_draws
    ),
    
    rate_ratio_lower_95 = unname(
      quantile(
        rr_draws,
        0.025
      )
    ),
    
    rate_ratio_upper_95 = unname(
      quantile(
        rr_draws,
        0.975
      )
    ),
    
    probability_above_average = mean(
      rr_draws
      >
        1
    )
  )
}


country_posterior_summary <- do.call(
  rbind,
  country_results
)


country_posterior_summary <- merge(
  
  country_lookup,
  
  country_posterior_summary,
  
  by = "country_id",
  
  all.x = TRUE,
  
  sort = FALSE
)


country_posterior_summary <- (
  country_posterior_summary[
    order(
      -country_posterior_summary$rate_ratio_median
    ),
  ]
)




print(
  country_posterior_summary
)



write.csv(
  
  country_posterior_summary,
  
  file.path(
    TABLE_DIR,
    "m3_country_posterior_summary.csv"
  ),
  
  row.names = FALSE
)



# ============================================================
# 10. Country relative-rate forest plot
# ============================================================

country_plot_data <- country_posterior_summary


country_plot_data$country_name <- factor(
  
  country_plot_data$country_name,
  
  levels = rev(
    country_plot_data$country_name
  )
)


country_effect_plot <- ggplot(
  
  country_plot_data,
  
  aes(
    x = rate_ratio_median,
    y = country_name
  )
  
) +
  
  geom_errorbarh(
    
    aes(
      xmin = rate_ratio_lower_95,
      xmax = rate_ratio_upper_95
    ),
    
    height = 0.2
  ) +
  
  geom_point(
    size = 2
  ) +
  
  geom_vline(
    
    xintercept = 1,
    
    linetype = "dashed"
  ) +
  
  labs(
    
    title = (
      "Posterior Country-Specific Railway Accident Rate Ratios"
    ),
    
    subtitle = (
      "Relative to the overall baseline after accounting for exposure and time"
    ),
    
    x = "Posterior rate ratio",
    
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 12
  )


print(
  country_effect_plot
)




ggsave(
  
  filename = file.path(
    FIGURE_DIR,
    "m3_country_relative_rates.png"
  ),
  
  plot = country_effect_plot,
  
  width = 9,
  
  height = 10,
  
  dpi = 300
)











# ============================================================
# 11. Annual trend comparison
# ============================================================

annual_change_table <- posterior_summary_all[
  posterior_summary_all$parameter
  ==
    "annual_percent_change",
]


annual_change_plot <- ggplot(
  
  annual_change_table,
  
  aes(
    x = model,
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
    
    title = (
      "Posterior Annual Change in Significant Railway Accident Rates"
    ),
    
    x = NULL,
    
    y = "Annual percentage change (%)"
  ) +
  
  theme_minimal(
    base_size = 12
  )


print(
  annual_change_plot
)



ggsave(
  
  filename = file.path(
    FIGURE_DIR,
    "annual_percent_change_by_model.png"
  ),
  
  plot = annual_change_plot,
  
  width = 9,
  
  height = 6,
  
  dpi = 300
)



# ============================================================
# 12. Key results for interpretation
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "POSTERIOR PROBABILITY OF DECLINE\n"
)

cat(
  "========================================\n"
)


print(
  probability_summary
)


cat(
  "\n========================================\n"
)

cat(
  "M3 MAIN POSTERIOR RESULTS\n"
)

cat(
  "========================================\n"
)


print(
  
  posterior_summary_M3[
    posterior_summary_M3$parameter
    %in%
      c(
        "baseline_rate_2017",
        "annual_rate_ratio",
        "annual_percent_change",
        "percent_change_2010_2024",
        "r",
        "sigma_country",
        "country_sd_multiplier"
      ),
  ]
)


cat(
  "\n========================================\n"
)

cat(
  "COUNTRIES WITH HIGHEST POSTERIOR RATES\n"
)

cat(
  "========================================\n"
)


print(
  
  head(
    country_posterior_summary,
    10
  )
)


cat(
  "\n========================================\n"
)

cat(
  "COUNTRIES WITH LOWEST POSTERIOR RATES\n"
)

cat(
  "========================================\n"
)


print(
  
  tail(
    country_posterior_summary,
    10
  )
)