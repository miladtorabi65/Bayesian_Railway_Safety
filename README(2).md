# Bayesian Railway Safety Analysis

## Overview

This project studies significant railway accident counts across European countries from 2010 to 2024.

The main goal is to estimate how accident rates changed over time after accounting for railway traffic exposure, and to determine whether persistent differences remained between countries.

Three Bayesian count models were compared:

- Poisson regression
- Negative Binomial regression
- Hierarchical Negative Binomial regression with country random effects

The hierarchical Negative Binomial model provided the best overall description of the data.

---

## Research question

How have significant railway accident rates changed over time across European countries, and how much persistent between-country variation remains after accounting for train-km exposure?

---

## Data

The analysis uses the European Union Agency for Railways Common Safety Indicators dataset.

The final dataset contains:

- 27 countries
- 15 years, from 2010 to 2024
- 405 country-year observations
- 0 missing values
- 0 duplicated country-year rows
- 7 zero-accident observations

The main variables are:

| Variable | Description |
|---|---|
| `accidents` | Total number of significant railway accidents |
| `train_km` | Railway traffic exposure |
| `year` | Calendar year |
| `year_centered` | Year centered at 2017 |
| `country_id` | Numeric identifier used for country effects |

The processed dataset is stored in:

```text
data/processed/railway_country_year_2010_2024.csv
```

The Channel Tunnel was excluded because it is not a country. The United Kingdom was excluded because complete train-km exposure was not available for the full study period.

---

## Why train-km is included as an offset

Raw accident counts are not directly comparable across countries because railway systems differ greatly in size.

A country with more train traffic has more opportunities for accidents, even when its underlying accident rate is similar to that of a smaller railway system.

For this reason, the models use train-km as an exposure offset.

The shared mean structure is:

```text
log(mu_it) = log(E_it) + alpha + beta * x_t
```

For the hierarchical model:

```text
log(mu_it) = log(E_it) + alpha + beta * x_t + u_i
```

```math
\log(\mu_{it})
=
\log(E_{it})
+
\alpha
+
\beta x_t
+
u_i.
```

```math
\mu_{it}
=
E_{it}
\exp(\alpha + \beta x_t + u_i)
```

where:

- `mu_it` is the expected accident count for country `i` in year `t`
- `E_it` is train-km exposure
- `x_t` is year centered at 2017
- `alpha` is the baseline log accident rate
- `beta` is the common yearly trend
- `u_i` is the country-specific random effect

---

## Exploratory findings

The exploratory analysis showed three important patterns:

1. accident rates declined over time;
2. accident counts were much more variable than a simple Poisson model would suggest;
3. countries showed persistent differences in exposure-adjusted accident rates.

Main descriptive statistics:

| Statistic | Value |
|---|---:|
| Total accidents | 27,121 |
| Mean count | 66.97 |
| Median count | 36 |
| Variance | 6,579.34 |
| Maximum count | 488 |
| Zero observations | 7 |
| Variance-to-mean ratio | 98.25 |

The aggregate accident rate declined from about 0.617 accidents per million train-km in 2010 to about 0.374 in 2024.

This descriptive pattern motivated the time-trend model, while the large variation and country differences motivated the Negative Binomial and hierarchical extensions.

---

## Models

### M1: Poisson regression

```text
Y_it ~ Poisson(mu_it)

log(mu_it) = log(E_it) + alpha + beta * x_t
```

```math
\log(\mu_{it})
=
\log(E_{it})
+
\alpha
+
\beta x_t
```

The Poisson model assumes:

```text
Var(Y_it | mu_it) = mu_it
```

This was used as the baseline model.

---

### M2: Negative Binomial regression

```text
Y_it ~ NegativeBinomial(mu_it, r)

log(mu_it) = log(E_it) + alpha + beta * x_t
```

Its variance is:

```text
Var(Y_it) = mu_it + (mu_it^2 / r)
```

The parameter `r` controls overdispersion.

- smaller `r` means stronger overdispersion;
- larger `r` makes the model closer to Poisson.

---

### M3: Hierarchical Negative Binomial regression

```text
Y_it ~ NegativeBinomial(mu_it, r)

log(mu_it) = log(E_it) + alpha + beta * x_t + u_i

u_i ~ Normal(0, sigma_country^2)
```

The country effect `u_i` is shared across all years for the same country.

It captures persistent differences that remain after accounting for exposure and the common time trend.

The final model uses a sum-to-zero constraint:

```text
sum(u_i) = 0
```

This separates the global intercept from the country effects and improves MCMC mixing.

---

## Priors

The final priors were:

| Parameter | Prior |
|---|---|
| `alpha` | Normal(log(0.5), 0.5^2) |
| `beta` | Normal(0, 0.05^2) |
| `log_r` | Normal(log(10), 1^2) |
| `sigma_country` | Half-Normal(0, 0.5) |

The original priors for `alpha` and `sigma_country` were wider. Prior predictive checks showed that they produced unrealistically extreme simulated datasets, so they were revised before fitting the final models.

Further details are available in:

```text
models/prior_specification.md
```

---

## Prior predictive checks

Prior predictive simulation was used to check whether the priors generated plausible accident datasets before fitting the observed outcomes.

For each model, 2,000 prior predictive datasets were generated using the real exposure, year, and country structure.

The checks included:

- total accidents
- mean count
- variance
- number of zeros
- maximum count
- country-rate variance
- yearly aggregate rates

The final priors were broad enough to include realistic datasets without being dominated by implausibly extreme outcomes.

---

## Parameter recovery

Parameter recovery was used to test whether the model implementation could recover known values from simulated data.

The formal study used:

- 10 independently simulated datasets
- the real 405-row panel structure
- known values of `alpha`, `beta`, `r`, and `sigma_country`
- convergence checks for every fitted model

All 10 replications passed the convergence criteria.

Main recovery results:

| Parameter | True value | Mean estimate | Bias | RMSE | 95% coverage |
|---|---:|---:|---:|---:|---:|
| `alpha` | -0.6931 | -0.6892 | 0.0039 | 0.0168 | 90% |
| `beta` | -0.03046 | -0.03004 | 0.00042 | 0.00291 | 100% |
| `r` | 30.0 | 30.66 | 0.66 | 2.53 | 100% |
| `sigma_country` | 0.80 | 0.770 | -0.030 | 0.094 | 90% |

Country-effect recovery was also strong:

- mean correlation between true and estimated country effects: about 0.995
- mean RMSE: about 0.074
- empirical 95% coverage: about 95.6%

These simulations showed that the implementation could recover the main parameters and country effects under the selected data-generating conditions.

---

## MCMC fitting and diagnostics

The Bayesian models were fitted in JAGS through R using `rjags` and `coda`.

The final analysis used:

- 4 chains
- 5,000 adaptation iterations
- 10,000 burn-in iterations
- up to 50,000 retained iterations per chain
- no thinning

Convergence was assessed using:

- trace plots
- posterior density plots
- autocorrelation plots
- Gelman-Rubin diagnostics
- effective sample sizes

The final diagnostics were strong:

- R-hat values were essentially 1
- effective sample sizes were large
- the chains mixed well
- country-effect diagnostics were also satisfactory

The complete final M3 posterior is stored in:

```text
results/mcmc/M3_Hierarchical_Negative_Binomial/posterior_samples_full.rds
```

---

## Main results

### Time trend

For the preferred hierarchical model:

```text
annual rate ratio = exp(beta) ≈ 0.967
```

This corresponds to an estimated annual decline of about:

```text
3.26% per year
```

The 95% credible interval was approximately:

```text
2.73% to 3.78% decline per year
```

The posterior probability of a negative time trend was effectively 1.

Across 2010 to 2024, the estimated reduction in the underlying accident rate was about:

```text
37.1%
```

---

### Country heterogeneity

The posterior median of the country-effect standard deviation was:

```text
sigma_country ≈ 0.819
```

On the multiplicative scale:

```text
exp(0.819) ≈ 2.27
```

This means that a one-standard-deviation difference in country effect corresponds to roughly a 2.27-fold difference in the underlying exposure-adjusted accident rate.

These country effects are residual model estimates. They should not be interpreted as causal effects or as definitive national safety rankings.

---

### Dispersion

The estimated Negative Binomial dispersion changed substantially after country effects were added:

| Model | Estimated `r` |
|---|---:|
| M2 | about 1.62 |
| M3 | about 32.84 |

This suggests that much of the variation treated as unstructured overdispersion in M2 was explained by persistent country differences in M3.

---

## Model comparison

DIC results:

| Model | DIC | Delta DIC | Rank |
|---|---:|---:|---:|
| M3 Hierarchical Negative Binomial | 2,981.30 | 0 | 1 |
| M2 Negative Binomial | 3,948.39 | 967.09 | 2 |
| M1 Poisson | 17,956.81 | 14,975.51 | 3 |

Lower DIC indicates a better balance between fit and effective model complexity.

M3 had the lowest DIC by a large margin.

---

## Posterior predictive checks

Posterior predictive checks compared replicated datasets from each fitted model with the observed data.

The main checks were:

- total accidents
- mean count
- variance
- number of zeros
- maximum count
- country-rate variance

Summary:

| Model | Checks reproduced within the final 95% predictive criterion |
|---|---:|
| M1 Poisson | 2 of 6 |
| M2 Negative Binomial | 1 of 6 |
| M3 Hierarchical Negative Binomial | 6 of 6 |

For country-rate variance:

- observed value: about 0.384
- M3 predictive median: about 0.377
- observed value was inside the M3 predictive interval

All 15 observed yearly aggregate rates were also inside the M3 posterior predictive intervals.

---

## Sensitivity analysis

The preferred M3 model was refitted after excluding 2020.

| Quantity | Main M3 | Excluding 2020 |
|---|---:|---:|
| Annual change | -3.26% | -3.14% |
| `r` | 32.85 | 33.46 |
| `sigma_country` | 0.819 | 0.825 |

The main conclusions remained almost unchanged.

This suggests that the estimated decline was not driven by the unusual 2020 observation.

---

## Frequentist comparison

Frequentist counterparts of the three Bayesian models were also fitted.

| Model | Bayesian annual change | Frequentist annual change |
|---|---:|---:|
| Poisson | -3.59% | -3.59% |
| Negative Binomial | -4.10% | -4.24% |
| Hierarchical Negative Binomial | -3.26% | -3.27% |

For the preferred hierarchical model:

| Quantity | Bayesian | Frequentist |
|---|---:|---:|
| Dispersion | 32.85 | 33.63 |
| Country-effect SD | 0.819 | about 0.84 |

The close agreement shows that the main conclusion was not strongly dependent on the inferential framework.

---

## Final model choice

The preferred model was M3, the hierarchical Negative Binomial model.

It was selected because it:

- converged well;
- recovered known parameters in simulation;
- had the lowest DIC;
- reproduced all six main posterior predictive statistics;
- captured persistent country differences;
- remained stable after excluding 2020;
- agreed closely with the corresponding frequentist model.

---

## Main conclusion

After accounting for train-km exposure, significant railway accident rates declined by about 3.3% per year across the studied European countries between 2010 and 2024.

This corresponds to an overall decline of about 37% over the full study period.

At the same time, substantial persistent differences remained between countries.

The hierarchical Negative Binomial model provided the most convincing description of both the common time trend and the remaining country-level heterogeneity.

---

## Limitations

The main limitations are:

- the model assumes one common linear time trend for all countries;
- country-specific time slopes were not included;
- serial dependence between consecutive years was not modeled explicitly;
- country effects are associative rather than causal;
- time-varying explanatory variables such as infrastructure investment, regulation, weather, and network composition were not included;
- reporting practices may differ between countries;
- countries without complete exposure data were excluded from the balanced panel.

---

## Repository structure

```text
Bayesian_Railway_Safety/
├── README.md
├── repository_tree.txt
├── Code/
│   ├── 01_preprocessing.py
│   ├── 02_exploratory_analysis.ipynb
│   ├── 03_prior_predictive_checks.ipynb
│   ├── 04_parameter_recovery_smoke_test.R
│   ├── 05_fit_bayesian_models.R
│   ├── 06_mcmc_diagnostics.R
│   ├── 07_posterior_inference.R
│   ├── 08_model_comparison_dic.R
│   ├── 09_posterior_predictive_checks.R
│   ├── 10_formal_parameter_recovery.R
│   ├── 11_sensitivity_analysis_exclude_2020.R
│   └── 12_frequentist_comparison.R
├── data/
│   ├── raw/
│   └── processed/
│       └── railway_country_year_2010_2024.csv
├── figures/
│   ├── eda/
│   ├── frequentist_comparison/
│   ├── mcmc_diagnostics/
│   ├── parameter_recovery/
│   ├── parameter_recovery_smoke_test/
│   ├── posterior_inference/
│   ├── posterior_predictive_checks/
│   ├── prior_predictive/
│   └── sensitivity_analysis/
├── models/
│   ├── model1_poisson.jags
│   ├── model2_negative_binomial.jags
│   ├── model3_hierarchical_negative_binomial.jags
│   ├── prior_specification.md
│   └── README.md
├── report/
├── results/
│   ├── frequentist_comparison/
│   ├── mcmc/
│   ├── parameter_recovery/
│   ├── parameter_recovery_smoke_test/
│   ├── prior_predictive/
│   └── sensitivity_analysis/
└── tables/
    ├── eda/
    ├── frequentist_comparison/
    ├── mcmc_diagnostics/
    ├── model_comparison/
    ├── parameter_recovery/
    ├── posterior_inference/
    ├── posterior_predictive_checks/
    ├── prior_predictive/
    └── sensitivity_analysis/
```

---

## Analysis order

Run the project in the following order:

```text
01_preprocessing.py
02_exploratory_analysis.ipynb
03_prior_predictive_checks.ipynb
04_parameter_recovery_smoke_test.R
05_fit_bayesian_models.R
06_mcmc_diagnostics.R
07_posterior_inference.R
08_model_comparison_dic.R
09_posterior_predictive_checks.R
10_formal_parameter_recovery.R
11_sensitivity_analysis_exclude_2020.R
12_frequentist_comparison.R
```

Each script writes its results to the corresponding `results/`, `tables/`, and `figures/` directories.

---

## Reproducing the analysis

### 1. Place the raw ERA files in `data/raw/`

The preprocessing script expects the original Common Safety Indicators workbook in the raw-data directory.

### 2. Install Python packages

```bash
pip install pandas numpy matplotlib jupyter openpyxl
```

### 3. Install R packages

```r
install.packages(
  c(
    "rjags",
    "coda",
    "ggplot2",
    "MASS",
    "glmmTMB"
  )
)
```

JAGS must also be installed separately.

### 4. Set the working directory

The scripts use relative paths and should be run from:

```text
Bayesian_Railway_Safety/Code/
```

### 5. Run the files in the documented order

The final M3 fit and formal recovery study are the most time-consuming parts of the analysis.

Existing posterior samples and output tables can be used to review the results without rerunning the complete workflow.

---

## Key output files

### Processed dataset

```text
data/processed/railway_country_year_2010_2024.csv
```

### Final M3 posterior

```text
results/mcmc/M3_Hierarchical_Negative_Binomial/posterior_samples_full.rds
```

### Main posterior summaries

```text
tables/posterior_inference/posterior_main_summary.csv
tables/posterior_inference/posterior_probability_summary.csv
tables/posterior_inference/m3_country_posterior_summary.csv
```

### MCMC diagnostics

```text
tables/mcmc_diagnostics/mcmc_diagnostics_summary.csv
tables/mcmc_diagnostics/country_effect_mcmc_diagnostics.csv
```

### Model comparison

```text
tables/model_comparison/dic_model_comparison.csv
```

### Posterior predictive checks

```text
tables/posterior_predictive_checks/posterior_predictive_summary.csv
tables/posterior_predictive_checks/m3_yearly_rate_ppc.csv
```

### Parameter recovery

```text
tables/parameter_recovery/parameter_recovery_summary.csv
tables/parameter_recovery/country_effect_recovery_summary.csv
tables/parameter_recovery/parameter_recovery_study_status.csv
```

### Sensitivity analysis

```text
tables/sensitivity_analysis/main_vs_excluding_2020_comparison.csv
```

### Frequentist comparison

```text
tables/frequentist_comparison/bayesian_vs_frequentist_annual_trend.csv
tables/frequentist_comparison/bayesian_vs_frequentist_structural_parameters.csv
tables/frequentist_comparison/frequentist_model_fit_statistics.csv
```

---

## Software

The project uses:

- Python
- Jupyter Notebook
- R
- JAGS
- `rjags`
- `coda`
- `ggplot2`
- `MASS`
- `glmmTMB`

---

## Author

**Milad Torabi**  
FSL II Final Project  
Sapienza University of Rome

Professor: Luca Tardella

---

## Data acknowledgment

The analysis uses data from the European Union Agency for Railways Common Safety Indicators programme.
