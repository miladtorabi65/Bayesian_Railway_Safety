# Prior Specification

The candidate Bayesian models use weakly informative priors chosen
on scales that are interpretable for the railway accident data.

## Shared regression parameters

### Baseline log rate

alpha ~ Normal(log(0.5), 1^2)

This implies a broad prior for the baseline accident rate in 2017.
Approximately 95% of the prior mass for exp(alpha) lies between
0.07 and 3.55 accidents per million train-km.

### Annual time trend

beta ~ Normal(0, 0.05^2)

The prior is centered on no temporal trend.

Approximately 95% of the implied annual rate ratios exp(beta)
lie between 0.91 and 1.10.

## Negative-binomial dispersion

For Models 2 and 3:

log(r) ~ Normal(log(10), 1^2)

r = exp(log(r))

This guarantees a positive dispersion parameter while allowing a
broad range of extra-Poisson variability.

The negative-binomial variance is:

Var(Y) = mu + mu^2 / r

Smaller r implies greater overdispersion.

## Country-level heterogeneity

For Model 3:

sigma_country ~ Half-Normal(0, 1)

The country effects are represented using the non-centred
parameterization:

z_i ~ Normal(0, 1)

u_i = sigma_country * z_i

The quantity exp(u_i) represents the country-specific relative
accident rate.

## Prior validation

These priors are provisional until prior predictive checks are
performed.

The next step will simulate hypothetical railway accident datasets
from the priors using the actual 2010–2024 year, country and train-km
structure.

Priors will be revised only if their prior predictive implications
are clearly implausible.


revised values:
alpha ~ dnorm(log(0.5), 4)

beta ~ dnorm(0, 400)

log_r ~ dnorm(log(10), 1)
r <- exp(log_r)

sigma_country ~ dnorm(0, 4) T(0,)