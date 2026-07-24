# ============================================================
# MODEL COMPARISON USING DIC
# ============================================================


# ============================================================
# 1. Paths
# ============================================================

DATA_PATH <- paste0(
  "../data/processed/",
  "railway_country_year_2010_2024.csv"
)

RESULT_DIR <- "../results/mcmc"

TABLE_DIR <- "../tables/model_comparison"


dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Load data
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


posterior_M1 <- as.matrix(
  samples_M1
)

posterior_M2 <- as.matrix(
  samples_M2
)

posterior_M3 <- as.matrix(
  samples_M3
)


cat(
  "Posterior samples loaded successfully.\n"
)


# ============================================================
# 4. Poisson mean deviance
# ============================================================

calculate_dbar_poisson <- function(
    posterior,
    y,
    exposure,
    year_centered,
    chunk_size = 5000
) {
  
  n_draws <- nrow(
    posterior
  )
  
  N <- length(
    y
  )
  
  total_deviance <- 0
  
  
  starts <- seq(
    1,
    n_draws,
    by = chunk_size
  )
  
  
  for (
    start
    in starts
  ) {
    
    end <- min(
      start + chunk_size - 1,
      n_draws
    )
    
    
    index <- start:end
    
    
    k <- length(
      index
    )
    
    
    alpha <- posterior[
      index,
      "alpha"
    ]
    
    
    beta <- posterior[
      index,
      "beta"
    ]
    
    
    eta <- (
      matrix(
        alpha,
        nrow = k,
        ncol = N
      )
      +
        beta
      %o%
        year_centered
    )
    
    
    mu <- (
      exp(
        eta
      )
      *
        matrix(
          exposure,
          nrow = k,
          ncol = N,
          byrow = TRUE
        )
    )
    
    
    y_matrix <- matrix(
      y,
      nrow = k,
      ncol = N,
      byrow = TRUE
    )
    
    
    log_likelihood <- dpois(
      y_matrix,
      lambda = mu,
      log = TRUE
    )
    
    
    total_deviance <- (
      total_deviance
      +
        sum(
          -2
          *
            rowSums(
              log_likelihood
            )
        )
    )
  }
  
  
  total_deviance/n_draws
}


# ============================================================
# 5. Negative-binomial mean deviance
# ============================================================

calculate_dbar_nb <- function(
    posterior,
    y,
    exposure,
    year_centered,
    country_id = NULL,
    u_parameters = NULL,
    chunk_size = 5000
) {
  
  n_draws <- nrow(
    posterior
  )
  
  N <- length(
    y
  )
  
  total_deviance <- 0
  
  
  starts <- seq(
    1,
    n_draws,
    by = chunk_size
  )
  
  
  for (
    start
    in starts
  ) {
    
    end <- min(
      start + chunk_size - 1,
      n_draws
    )
    
    
    index <- start:end
    
    
    k <- length(
      index
    )
    
    
    alpha <- posterior[
      index,
      "alpha"
    ]
    
    
    beta <- posterior[
      index,
      "beta"
    ]
    
    
    r <- posterior[
      index,
      "r"
    ]
    
    
    eta <- (
      matrix(
        alpha,
        nrow = k,
        ncol = N
      )
      +
        beta
      %o%
        year_centered
    )
    
    
    # Add country effects only for M3
    if (
      !is.null(
        u_parameters
      )
    ) {
      
      u_draws <- posterior[
        index,
        u_parameters,
        drop = FALSE
      ]
      
      
      eta <- (
        eta
        +
          u_draws[
            ,
            country_id,
            drop = FALSE
          ]
      )
    }
    
    
    mu <- (
      exp(
        eta
      )
      *
        matrix(
          exposure,
          nrow = k,
          ncol = N,
          byrow = TRUE
        )
    )
    
    
    y_matrix <- matrix(
      y,
      nrow = k,
      ncol = N,
      byrow = TRUE
    )
    
    
    r_matrix <- matrix(
      r,
      nrow = k,
      ncol = N
    )
    
    
    log_likelihood <- dnbinom(
      y_matrix,
      size = r_matrix,
      mu = mu,
      log = TRUE
    )
    
    
    total_deviance <- (
      total_deviance
      +
        sum(
          -2
          *
            rowSums(
              log_likelihood
            )
        )
    )
  }
  
  
  total_deviance/n_draws
}


# ============================================================
# 6. Deviance at posterior mean — M1
# ============================================================

calculate_dhat_M1 <- function(
    posterior
) {
  
  alpha_mean <- mean(
    posterior[
      ,
      "alpha"
    ]
  )
  
  
  beta_mean <- mean(
    posterior[
      ,
      "beta"
    ]
  )
  
  
  mu <- (
    exposure
    *
      exp(
        alpha_mean
        +
          beta_mean
        *
          year_centered
      )
  )
  
  
  -2 * sum(
      dpois(
        y,
        lambda = mu,
        log = TRUE
      )
    )
}


# ============================================================
# 7. Deviance at posterior mean — M2
# ============================================================

calculate_dhat_M2 <- function(
    posterior
) {
  
  alpha_mean <- mean(
    posterior[
      ,
      "alpha"
    ]
  )
  
  
  beta_mean <- mean(
    posterior[
      ,
      "beta"
    ]
  )
  
  
  r_mean <- mean(
    posterior[
      ,
      "r"
    ]
  )
  
  
  mu <- (
    exposure
    *
      exp(
        alpha_mean
        +
          beta_mean
        *
          year_centered
      )
  )
  
  
  -2 * sum(
      dnbinom(
        y,
        size = r_mean,
        mu = mu,
        log = TRUE
      )
    )
}


# ============================================================
# 8. Find and order M3 country effects
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
# 9. Deviance at posterior mean — M3
# ============================================================

calculate_dhat_M3 <- function(
    posterior,
    u_parameters
) {
  
  alpha_mean <- mean(
    posterior[
      ,
      "alpha"
    ]
  )
  
  
  beta_mean <- mean(
    posterior[
      ,
      "beta"
    ]
  )
  
  
  r_mean <- mean(
    posterior[
      ,
      "r"
    ]
  )
  
  
  u_mean <- colMeans(
    posterior[
      ,
      u_parameters,
      drop = FALSE
    ]
  )
  
  
  mu <- (
    exposure
    *
      exp(
        alpha_mean
        +
          beta_mean
        *
          year_centered
        +
          u_mean[
            country_id
          ]
      )
  )
  
  
  -2 * sum(
      dnbinom(
        y,
        size = r_mean,
        mu = mu,
        log = TRUE
      )
    )
}


# ============================================================
# 10. Calculate DIC for M1
# ============================================================

cat(
  "\nCalculating DIC for M1...\n"
)


Dbar_M1 <- calculate_dbar_poisson(
  posterior_M1,
  y,
  exposure,
  year_centered
)


Dhat_M1 <- calculate_dhat_M1(
  posterior_M1
)


pD_M1 <- (
  Dbar_M1
  -
    Dhat_M1
)


DIC_M1 <- (
  Dbar_M1
  +
    pD_M1
)


# ============================================================
# 11. Calculate DIC for M2
# ============================================================

cat(
  "\nCalculating DIC for M2...\n"
)


Dbar_M2 <- calculate_dbar_nb(
  posterior = posterior_M2,
  y = y,
  exposure = exposure,
  year_centered = year_centered
)


Dhat_M2 <- calculate_dhat_M2(
  posterior_M2
)


pD_M2 <- (
  Dbar_M2
  -
    Dhat_M2
)


DIC_M2 <- (
  Dbar_M2
  +
    pD_M2
)


# ============================================================
# 12. Calculate DIC for M3
# ============================================================

cat(
  "\nCalculating DIC for M3...\n"
)


Dbar_M3 <- calculate_dbar_nb(
  posterior = posterior_M3,
  y = y,
  exposure = exposure,
  year_centered = year_centered,
  country_id = country_id,
  u_parameters = u_parameters
)


Dhat_M3 <- calculate_dhat_M3(
  posterior_M3,
  u_parameters
)


pD_M3 <- (
  Dbar_M3
  -
    Dhat_M3
)


DIC_M3 <- (
  Dbar_M3
  +
    pD_M3
)


# ============================================================
# 13. Build comparison table
# ============================================================

dic_comparison <- data.frame(
  
  model = c(
    "M1_Poisson",
    "M2_Negative_Binomial",
    "M3_Hierarchical_Negative_Binomial"
  ),
  
  Dbar = c(
    Dbar_M1,
    Dbar_M2,
    Dbar_M3
  ),
  
  Dhat = c(
    Dhat_M1,
    Dhat_M2,
    Dhat_M3
  ),
  
  pD = c(
    pD_M1,
    pD_M2,
    pD_M3
  ),
  
  DIC = c(
    DIC_M1,
    DIC_M2,
    DIC_M3
  )
)


# Lower DIC is better

dic_comparison$delta_DIC <- (
  dic_comparison$DIC
  -
    min(
      dic_comparison$DIC
    )
)


dic_comparison <- dic_comparison[
  order(
    dic_comparison$DIC
  ),
]


rownames(
  dic_comparison
) <- NULL


# ============================================================
# 14. Print and save results
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "DIC MODEL COMPARISON\n"
)

cat(
  "========================================\n\n"
)


print(
  dic_comparison,
  row.names = FALSE
)


write.csv(
  
  dic_comparison,
  
  file.path(
    TABLE_DIR,
    "dic_model_comparison.csv"
  ),
  
  row.names = FALSE
)


cat(
  "\nDIC analysis completed successfully.\n"
)