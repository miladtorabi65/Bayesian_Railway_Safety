# ============================================================
# MCMC DIAGNOSTICS
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

library(coda)


# ============================================================
# 2. Paths
# ============================================================

RESULT_DIR <- "../results/mcmc"

FIGURE_DIR <- "../figures/mcmc_diagnostics"

TABLE_DIR <- "../tables/mcmc_diagnostics"


dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 3. Load saved posterior samples
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


cat("Posterior samples loaded successfully.\n")


# ============================================================
# 4. Important parameters
# ============================================================

parameters_M1 <- c(
  "alpha",
  "beta"
)


parameters_M2 <- c(
  "alpha",
  "beta",
  "r"
)


parameters_M3 <- c(
  "alpha",
  "beta",
  "r",
  "sigma_country"
)



# ============================================================
# 5. Numerical diagnostic function
# ============================================================

calculate_diagnostics <- function(
    samples,
    parameters,
    model_name
) {
  
  selected_samples <- samples[
    ,
    parameters
  ]
  
  
  # ==========================================================
  # R-hat / Gelman-Rubin
  # ==========================================================
  
  gelman_results <- gelman.diag(
    selected_samples,
    autoburnin = FALSE,
    multivariate = FALSE
  )
  
  
  rhat_table <- (
    gelman_results$psrf
  )
  
  
  # ==========================================================
  # Effective sample size
  # ==========================================================
  
  ess_results <- effectiveSize(
    selected_samples
  )
  
  
  # ==========================================================
  # Posterior standard deviation
  # ==========================================================
  
  posterior_matrix <- as.matrix(
    selected_samples
  )
  
  
  posterior_sd <- apply(
    posterior_matrix,
    2,
    sd
  )
  
  
  # ==========================================================
  # Approximate Monte Carlo standard error
  #
  # MCSE ≈ posterior SD / sqrt(ESS)
  # ==========================================================
  
  mcse <- (
    posterior_sd
    /
      sqrt(
        ess_results
      )
  )
  
  
  # ==========================================================
  # Create final table
  # ==========================================================
  
  diagnostics <- data.frame(
    
    model = model_name,
    
    parameter = parameters,
    
    rhat = (
      rhat_table[
        parameters,
        "Point est."
      ]
    ),
    
    rhat_upper = (
      rhat_table[
        parameters,
        "Upper C.I."
      ]
    ),
    
    ess = (
      ess_results[
        parameters
      ]
    ),
    
    posterior_sd = (
      posterior_sd[
        parameters
      ]
    ),
    
    mcse = (
      mcse[
        parameters
      ]
    )
  )
  
  
  # MCSE relative to posterior SD
  diagnostics[
    "mcse_to_sd_percent"
  ] <- (
    
    100
    *
      diagnostics$mcse
    /
      diagnostics$posterior_sd
    
  )
  
  
  rownames(
    diagnostics
  ) <- NULL
  
  
  return(
    diagnostics
  )
}



# ============================================================
# 6. Calculate diagnostics
# ============================================================

diagnostics_M1 <- calculate_diagnostics(
  
  samples = samples_M1,
  
  parameters = parameters_M1,
  
  model_name = "M1_Poisson"
)


diagnostics_M2 <- calculate_diagnostics(
  
  samples = samples_M2,
  
  parameters = parameters_M2,
  
  model_name = "M2_Negative_Binomial"
)


diagnostics_M3 <- calculate_diagnostics(
  
  samples = samples_M3,
  
  parameters = parameters_M3,
  
  model_name = (
    "M3_Hierarchical_Negative_Binomial"
  )
)



all_diagnostics <- rbind(
  
  diagnostics_M1,
  
  diagnostics_M2,
  
  diagnostics_M3
)


print(
  all_diagnostics
)



write.csv(
  
  all_diagnostics,
  
  file.path(
    TABLE_DIR,
    "mcmc_diagnostics_summary.csv"
  ),
  
  row.names = FALSE
)




# ============================================================
# 7. Diagnostic status
# ============================================================

all_diagnostics$status <- ifelse(
  
  all_diagnostics$rhat < 1.01
  &
    all_diagnostics$ess >= 1000
  &
    all_diagnostics$mcse_to_sd_percent < 5,
  
  "Good",
  
  "Review"
)


print(
  all_diagnostics
)


write.csv(
  
  all_diagnostics,
  
  file.path(
    TABLE_DIR,
    "mcmc_diagnostics_summary.csv"
  ),
  
  row.names = FALSE
)




# ============================================================
# 8. Function for separate trace plots
# ============================================================

save_trace_plots <- function(
    samples,
    parameters,
    model_name
) {
  
  model_dir <- file.path(
    FIGURE_DIR,
    model_name,
    "traceplots"
  )
  
  
  dir.create(
    model_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  for (
    parameter
    in parameters
  ) {
    
    png(
      
      filename = file.path(
        model_dir,
        paste0(
          "trace_",
          parameter,
          ".png"
        )
      ),
      
      width = 1600,
      
      height = 900,
      
      res = 150
    )
    
    
    traceplot(
      samples[
        ,
        parameter
      ],
      
      main = paste(
        "Trace plot:",
        model_name,
        "-",
        parameter
      )
    )
    
    
    dev.off()
  }
}



save_trace_plots(
  samples_M1,
  parameters_M1,
  "M1_Poisson"
)


save_trace_plots(
  samples_M2,
  parameters_M2,
  "M2_Negative_Binomial"
)


save_trace_plots(
  samples_M3,
  parameters_M3,
  "M3_Hierarchical_Negative_Binomial"
)




# ============================================================
# 9. Posterior density plots
# ============================================================

save_density_plots <- function(
    samples,
    parameters,
    model_name
) {
  
  model_dir <- file.path(
    FIGURE_DIR,
    model_name,
    "density"
  )
  
  
  dir.create(
    model_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  for (
    parameter
    in parameters
  ) {
    
    png(
      
      filename = file.path(
        model_dir,
        paste0(
          "density_",
          parameter,
          ".png"
        )
      ),
      
      width = 1400,
      
      height = 900,
      
      res = 150
    )
    
    
    densplot(
      
      samples[
        ,
        parameter
      ],
      
      main = paste(
        "Posterior density:",
        model_name,
        "-",
        parameter
      )
    )
    
    
    dev.off()
  }
}



save_density_plots(
  samples_M1,
  parameters_M1,
  "M1_Poisson"
)


save_density_plots(
  samples_M2,
  parameters_M2,
  "M2_Negative_Binomial"
)


save_density_plots(
  samples_M3,
  parameters_M3,
  "M3_Hierarchical_Negative_Binomial"
)


# ============================================================
# 10. Autocorrelation plots
# ============================================================

save_autocorrelation_plots <- function(
    samples,
    parameters,
    model_name
) {
  
  model_dir <- file.path(
    FIGURE_DIR,
    model_name,
    "autocorrelation"
  )
  
  
  dir.create(
    model_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  for (
    parameter
    in parameters
  ) {
    
    png(
      
      filename = file.path(
        model_dir,
        paste0(
          "autocorrelation_",
          parameter,
          ".png"
        )
      ),
      
      width = 1600,
      
      height = 1200,
      
      res = 150
    )
    
    
    autocorr.plot(
      
      samples[
        ,
        parameter
      ],
      
      main = paste(
        "Autocorrelation:",
        model_name,
        "-",
        parameter
      )
    )
    
    
    dev.off()
  }
}




save_autocorrelation_plots(
  samples_M1,
  parameters_M1,
  "M1_Poisson"
)


save_autocorrelation_plots(
  samples_M2,
  parameters_M2,
  "M2_Negative_Binomial"
)


save_autocorrelation_plots(
  samples_M3,
  parameters_M3,
  "M3_Hierarchical_Negative_Binomial"
)


# ============================================================
# 11. Country-effect diagnostics
# ============================================================

all_parameter_names <- colnames(
  as.matrix(
    samples_M3
  )
)

print(
  all_parameter_names
)

u_parameters <- all_parameter_names[
  startsWith(
    all_parameter_names,
    "u["
  )
]


length(
  u_parameters
)




country_effect_diagnostics <- (
  calculate_diagnostics(
    
    samples = samples_M3,
    
    parameters = u_parameters,
    
    model_name = "M3_country_effects"
    
  )
)

print(
  country_effect_diagnostics
)


write.csv(
  
  country_effect_diagnostics,
  
  file.path(
    TABLE_DIR,
    "country_effect_mcmc_diagnostics.csv"
  ),
  
  row.names = FALSE
)




# ============================================================
# 12. Worst country-effect diagnostics
# ============================================================

cat(
  "\nWorst R-hat values:\n"
)


print(
  
  head(
    
    country_effect_diagnostics[
      
      order(
        -country_effect_diagnostics$rhat
      ),
      
    ],
    
    10
  )
)


cat(
  "\nLowest effective sample sizes:\n"
)


print(
  
  head(
    
    country_effect_diagnostics[
      
      order(
        country_effect_diagnostics$ess
      ),
      
    ],
    
    10
  )
)

