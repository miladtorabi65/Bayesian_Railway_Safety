# Bayesian Model Specifications

This folder contains the three candidate Bayesian models used in the
European railway safety project.

## Model 1

Poisson regression with:

- significant accident count as the outcome;
- log(train-km) as an exposure offset;
- centred year as the temporal predictor.

## Model 2

Negative binomial regression with the same mean structure as Model 1,
plus an overdispersion parameter.

The parameterization is:

E(Y) = mu

Var(Y) = mu + mu^2 / r

For JAGS:

p = r / (r + mu)

and:

Y ~ dnegbin(p, r)

## Model 3

Hierarchical negative binomial regression with:

- the same temporal trend;
- the same exposure adjustment;
- negative binomial overdispersion;
- a country-specific random intercept.

Country effects follow:

u_i ~ Normal(0, sigma_country^2)

The quantity exp(u_i) represents the country's model-based relative
accident rate compared with the average country.

## Status

The likelihood structures were frozen in Step 4.

The prior distributions will be selected and added in Step 5 after
prior predictive analysis.