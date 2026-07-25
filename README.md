# Bayesian Railway Safety Analysis

## Overview

This project investigates long-term changes in significant railway accident rates across European countries using fully Bayesian count models.

The analysis focuses on two questions:

1. **How did significant railway accident rates change from 2010 to 2024 after accounting for railway traffic exposure?**
2. **How much persistent between-country variation remained after accounting for exposure and the common time trend?**

Three Bayesian models of increasing complexity were fitted and compared:

- **M1 — Poisson regression**
- **M2 — Negative Binomial regression**
- **M3 — Hierarchical Negative Binomial regression with country random effects**

The preferred model was **M3**, because it provided the best balance of fit, predictive performance, convergence, interpretability, and robustness.

---

## Research Question

> How have significant railway accident rates changed over time across European countries, and how much persistent between-country variation remains after accounting for railway traffic exposure?

---

## Data

### Source

The project uses the **European Union Agency for Railways Common Safety Indicators dataset**.

The analysis uses:

- **Outcome:** `N00` — total significant railway accidents
- **Exposure:** `R01` — total train-km
- **Observation unit:** one country in one year
- **Period:** 2010–2024
- **Countries:** 27
- **Final balanced panel:** 405 country-year observations

The processed analysis dataset is stored at:

```text
data/processed/railway_country_year_2010_2024.csv
```

### Final validation checks

The locked analysis dataset contains:

- 405 rows
- 27 countries
- 15 years
- 0 missing values
- 0 duplicated country-year observations
- 7 zero-accident observations

The Channel Tunnel was excluded because it is not a country. The United Kingdom was excluded because complete train-km exposure was not available for the full selected period.

### Main processed variables

| Variable | Description |
|---|---|
| `year` | Calendar year |
| `country` | Country code |
| `country_name` | Country name |
| `country_id` | Integer identifier used for hierarchical effects |
| `accidents` | Total significant railway accidents |
| `train_km` | Railway traffic exposure |
| `year_centered` | `year - 2017` |
| `accidents_per_million_train_km` | Descriptive exposure-adjusted rate |

Centering year at 2017 makes the model intercept interpretable as the baseline log accident rate near the middle of the study period.

---

## Why an Exposure Offset Is Required

Raw accident counts are not directly comparable across countries because railway systems differ substantially in traffic volume.

The expected count is therefore modeled as:

```math
\mu_{it}
=
E_{it}
\exp(\alpha + \beta x_t + u_i)
```

where $E_{it}$ is train-km exposure.

Equivalently, on the log scale:

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

The term $\log(E_{it})$ is an **offset** with coefficient fixed at one. This makes the model describe accident **rates** while retaining an appropriate count-data likelihood.

---

## Exploratory Findings

The exploratory analysis identified three main features:

1. **A declining aggregate accident rate over time**
2. **Strong overdispersion in accident counts**
3. **Persistent differences between countries**

### Main descriptive statistics

| Statistic | Observed value |
|---|---:|
| Total accidents | 27,121 |
| Mean count | 66.97 |
| Median count | 36 |
| Variance | 6,579.34 |
| Maximum count | 488 |
| Zero-count observations | 7 |
| Variance-to-mean ratio | 98.25 |

The aggregate exposure-adjusted rate declined from approximately:

- **0.617 accidents per million train-km in 2010**
- to **0.374 accidents per million train-km in 2024**

This corresponds to an unadjusted descriptive reduction of about 39%.

Observed country-level rates ranged from approximately 0.12 to 2.25 accidents per million train-km. These descriptive differences motivated the hierarchical model but should not be interpreted as causal or definitive safety rankings.

---

## Bayesian Models

### M1 — Poisson Regression

```math
Y_{it}
\sim
\operatorname{Poisson}(\mu_{it})
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

with:

```math
\operatorname{Var}(Y_{it}\mid\mu_{it})
=
\mu_{it}
```

M1 is the simplest exposure-adjusted count model. Its main limitation is the assumption that the conditional variance equals the conditional mean.

### M2 — Negative Binomial Regression

```math
Y_{it}
\sim
\operatorname{NB}(\mu_{it}, r)
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

with:

```math
\operatorname{Var}(Y_{it})
=
\mu_{it}
+
\frac{\mu_{it}^{2}}{r}
```

The parameter $r$ controls extra-Poisson variation:

- smaller $r$: stronger overdispersion
- larger $r$: distribution closer to Poisson

M2 allows extra variation but does not explicitly model persistent country structure.

### M3 — Hierarchical Negative Binomial Regression

```math
Y_{it}
\sim
\operatorname{NB}(\mu_{it}, r)
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
u_i
```

```math
u_i
\sim
\operatorname{Normal}(0,\sigma_{\text{country}}^2)
```

The country effect $u_i$ is shared by all observations from country $i$. It represents a persistent deviation from the common exposure-adjusted time trend.

For interpretation:

```math
\exp(u_i)
```

is the country-specific relative rate.

The final implementation constrains the country effects to sum to zero:

```math
\sum_i u_i = 0
```

This separates the global intercept from the country effects and improves MCMC mixing and identifiability.

---

## Priors

The final priors were selected after prior predictive simulation.

| Parameter | Final prior |
|---|---|
| `alpha` | Normal(log(0.5), 0.5²) |
| `beta` | Normal(0, 0.05²) |
| `log_r` | Normal(log(10), 1²) |
| `sigma_country` | Half-Normal(0, 0.5) |

The initial priors for $\alpha$ and $\sigma_{\text{country}}$ were broader. Prior predictive checks showed that they generated implausibly extreme hierarchical datasets, so they were tightened before fitting the observed accident outcomes.

Detailed prior documentation is available in:

```text
models/prior_specification.md
```

---

## Prior Predictive Checks

Two thousand prior predictive datasets were generated for each model using the real exposure, time, and country structure.

The checks included:

- total accidents
- mean count
- variance
- number of zeros
- maximum count
- country-rate variance
- yearly aggregate rates

The final priors generated broadly plausible railway accident datasets.

An important structural finding was that the observed country-rate variance could be generated by M3, but not adequately by M1 or M2. This supported the inclusion of country random effects before fitting the observed accident outcomes.

---

## Bayesian Computation

The Bayesian models were fitted in **JAGS** through **R** using `rjags` and `coda`.

The final analysis used:

- multiple independent chains
- adaptation
- burn-in
- retained posterior samples without thinning
- trace plots
- density plots
- autocorrelation plots
- Gelman-Rubin diagnostics
- effective sample sizes

Exact settings and runtimes are stored in:

```text
results/mcmc/*/runtime_and_settings.csv
```

The complete final M3 posterior is:

```text
results/mcmc/M3_Hierarchical_Negative_Binomial/posterior_samples_full.rds
```

This is the canonical M3 posterior used for final diagnostics, inference, DIC, and posterior predictive checks.

---

## MCMC Diagnostics

The final chains showed good mixing and stable posterior sampling.

The main diagnostics indicated:

- $\hat R$ values essentially equal to 1
- large effective sample sizes
- no visible non-stationarity in the main trace plots
- well-overlapping chain-specific posterior densities
- acceptable autocorrelation behavior

Diagnostics for the 27 country effects were also evaluated.

Outputs are stored in:

```text
tables/mcmc_diagnostics/
figures/mcmc_diagnostics/
```

---

## Parameter Recovery

Two forms of parameter recovery were performed.

### Smoke test

A single simulated dataset with known parameters was used to verify the complete model implementation and country-effect extraction.

### Formal repeated recovery

The final formal study used:

- 10 independently simulated datasets
- the real 405-row exposure, year, and country design
- known values of $\alpha$, $\beta$, $r$, and $\sigma_{\text{country}}$
- convergence checks for every replication

All 10 replications passed the predefined convergence criteria.

### Global parameter recovery

| Parameter | True value | Mean recovered estimate | Bias | RMSE | 95% coverage |
|---|---:|---:|---:|---:|---:|
| `alpha` | -0.6931 | -0.6892 | 0.0039 | 0.0168 | 90% |
| `beta` | -0.03046 | -0.03004 | 0.00042 | 0.00291 | 100% |
| `r` | 30.0 | 30.66 | 0.66 | 2.53 | 100% |
| `sigma_country` | 0.80 | 0.770 | -0.030 | 0.094 | 90% |

### Country-effect recovery

- Mean correlation between true and estimated country effects: **0.995**
- Mean RMSE: **0.074**
- Empirical 95% interval coverage: **95.6%**
- Convergence failure rate: **0%**

These results validate the implementation under realistic simulated conditions. They do not prove that the real-data model is true.

---

## Main Bayesian Results

### Common temporal trend

Under the preferred M3 model:

```math
\exp(\beta)
\approx
0.967.
```

This corresponds to an estimated annual rate change of:

```math
-3.26\%
```

with a 95% credible interval of approximately:

```math
[-3.78\%, -2.73\%].
```

The posterior probability of a declining trend was effectively 1 in the retained posterior sample.

Across 2010–2024, the estimated reduction in the underlying accident rate was approximately:

```math
37.1\%
```

with a 95% credible interval of approximately 32.2% to 41.7%.

### Persistent country heterogeneity

For M3:

```math
\sigma_{\text{country}}
\approx
0.819.
```

A one-standard-deviation difference on the country-effect scale corresponds to:

```math
\exp(0.819)
\approx
2.27.
```

Therefore, a one-standard-deviation difference in country effect corresponds to roughly a 2.27-fold multiplicative difference in the underlying exposure-adjusted accident rate.

Country effects should be interpreted as persistent residual associations after accounting for exposure and the common time trend. They are not causal effects or definitive national safety rankings.

### Dispersion before and after country effects

| Model | Estimated `r` |
|---|---:|
| M2 | approximately 1.62 |
| M3 | approximately 32.84 |

Because larger $r$ implies less residual overdispersion, this suggests that much of the variation treated as unstructured overdispersion in M2 was explained by persistent country differences in M3.

---

## Model Comparison

### DIC

| Model | DIC | ΔDIC | Rank |
|---|---:|---:|---:|
| M3 Hierarchical Negative Binomial | 2,981.30 | 0 | 1 |
| M2 Negative Binomial | 3,948.39 | 967.09 | 2 |
| M1 Poisson | 17,956.81 | 14,975.51 | 3 |

Lower DIC indicates a better balance between fit and effective complexity.

DIC strongly favored M3, but final model selection was not based on DIC alone. It also considered convergence, parameter recovery, posterior predictive performance, robustness, and scientific interpretation.

---

## Posterior Predictive Checks

Replicated datasets were generated from the fitted posterior distributions and compared with the observed data.

The checks included:

- total accidents
- mean count
- variance
- zeros
- maximum count
- country-rate variance
- yearly aggregate rate trajectory

### Summary

| Model | Observed features inside 95% predictive intervals |
|---|---:|
| M1 Poisson | 2 of 6 |
| M2 Negative Binomial | 1 of 6 under the final criterion |
| M3 Hierarchical Negative Binomial | 6 of 6 |

For country-rate variance:

- observed: approximately 0.384
- M3 predictive median: approximately 0.377
- observed value inside M3 95% predictive interval

All 15 observed yearly aggregate rates were also inside the corresponding M3 posterior predictive intervals.

This was the strongest evidence that M3 reproduced the important structure of the observed dataset.

---

## Sensitivity Analysis

The preferred M3 model was refitted after excluding 2020.

| Quantity | Main M3 | Excluding 2020 |
|---|---:|---:|
| Annual change | -3.26% | -3.14% |
| `r` | 32.85 | 33.46 |
| `sigma_country` | 0.819 | 0.825 |

The sensitivity model converged successfully, and the main conclusions were essentially unchanged.

Therefore, the estimated decline was not driven by the unusual 2020 observation.

---

## Frequentist Comparison

Frequentist counterparts were fitted using:

- Poisson GLM
- Negative Binomial regression
- Negative Binomial mixed-effects regression

The frequentist hierarchical model converged successfully.

### Annual trend comparison

| Model | Bayesian annual change | Frequentist annual change |
|---|---:|---:|
| Poisson | -3.59% | -3.59% |
| Negative Binomial | -4.10% | -4.24% |
| Hierarchical Negative Binomial | -3.26% | -3.27% |

### Structural comparison for the preferred model

| Quantity | Bayesian | Frequentist |
|---|---:|---:|
| M3 dispersion | 32.85 | 33.63 |
| Country-effect SD | 0.819 | approximately 0.84 |

AIC and BIC also preferred the hierarchical Negative Binomial model.

The close agreement suggests that the principal substantive conclusions were not strongly dependent on the Bayesian prior specification or inferential framework.

---

## Preferred Model

The final preferred model is:

```math
\boxed{\text{M3 — Hierarchical Negative Binomial}}
```

It was selected because it:

- had excellent MCMC convergence
- successfully recovered known parameters in simulation
- had the lowest DIC
- reproduced all major observed-data features in posterior predictive checks
- captured persistent country heterogeneity
- remained stable after excluding 2020
- agreed closely with the corresponding frequentist model

---

## Main Conclusion

After accounting for train-km exposure:

- significant railway accident rates declined by approximately **3.3% per year**
- the underlying rate declined by approximately **37% from 2010 to 2024**
- substantial persistent differences remained between countries
- the preferred hierarchical model estimated a country-effect standard deviation of approximately **0.82**

The analysis therefore supports both a broad improvement in railway safety over time and meaningful remaining cross-country heterogeneity.

---

## Limitations

The main limitations are:

- one common log-linear time trend was assumed for all countries
- country-specific slopes were not included
- serial dependence across consecutive years was not explicitly modeled
- country effects are associative rather than causal
- no time-varying explanatory covariates were included for infrastructure, regulation, weather, network composition, or reporting practices
- differences in reporting systems may remain
- the balanced panel excludes countries without complete exposure data

These limitations should be considered when interpreting the country-specific posterior effects.

---

## Repository Structure

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

## Analysis Workflow

Run the analysis in the following order:

```text
01_preprocessing.py
        ↓
02_exploratory_analysis.ipynb
        ↓
03_prior_predictive_checks.ipynb
        ↓
04_parameter_recovery_smoke_test.R
        ↓
05_fit_bayesian_models.R
        ↓
06_mcmc_diagnostics.R
        ↓
07_posterior_inference.R
        ↓
08_model_comparison_dic.R
        ↓
09_posterior_predictive_checks.R
        ↓
10_formal_parameter_recovery.R
        ↓
11_sensitivity_analysis_exclude_2020.R
        ↓
12_frequentist_comparison.R
```

The preprocessing script creates the locked processed dataset. The remaining scripts and notebooks read from that file and write outputs to the corresponding `results/`, `tables/`, and `figures/` directories.

---

## Reproducing the Analysis

### 1. Clone or download the project

```bash
git clone <repository-url>
cd Bayesian_Railway_Safety
```

### 2. Place the raw ERA files in `data/raw/`

The preprocessing script expects the original Common Safety Indicators workbook and supporting documentation in the raw-data directory.

### 3. Install Python dependencies

```bash
pip install pandas numpy matplotlib jupyter openpyxl
```

Additional notebook-specific packages may be listed in the notebook import cells.

### 4. Install R packages

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

JAGS must also be installed separately and available to `rjags`.

### 5. Run from the `Code/` directory

The scripts use relative paths. Set the working directory to:

```text
Bayesian_Railway_Safety/Code/
```

Then run the files in the documented analysis order.

### 6. Large and time-consuming steps

The following steps can take substantial time:

- final M3 sampling
- full MCMC diagnostics
- repeated parameter recovery

Existing posterior samples and generated tables can be used for review without rerunning all models.

---

## Key Output Files

### Processed dataset

```text
data/processed/railway_country_year_2010_2024.csv
```

### Final M3 posterior

```text
results/mcmc/M3_Hierarchical_Negative_Binomial/posterior_samples_full.rds
```

### Main posterior results

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

### DIC comparison

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

The analysis uses:

- Python
- Jupyter Notebook
- R
- JAGS
- `rjags`
- `coda`
- `ggplot2`
- `MASS`
- `glmmTMB`

Exact R session information is stored in the corresponding results directories.

---

## Interpretation Warning

Country-specific relative rates describe persistent residual differences after adjustment for exposure and time.

They should **not** be interpreted as:

- causal country effects
- definitive national safety rankings
- proof that one railway system is intrinsically safer or less safe than another

Unmeasured structural, regulatory, operational, reporting, and environmental factors may contribute to the estimated country effects.

---

## Author

**Milad Torabi**  
FSL II — Final Project  
Sapienza University of Rome

Professor: Luca Tardella

---

## Acknowledgment

The analysis uses data from the European Union Agency for Railways Common Safety Indicators programme.
