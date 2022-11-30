[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)[![Downloads Total](https://cranlogs.r-pkg.org/badges/grand-total/latentFactoR?color=brightgreen)](https://cran.r-project.org/package=latentFactoR) [![Downloads per month](http://cranlogs.r-pkg.org/badges/latentFactoR?color=brightgreen)](https://cran.r-project.org/package=latentFactoR) 

### latentFactoR: Data Simulation based on Latent Factors

# How to Use

The main goal of the {latentFactoR} package in R is to provide a straightforward, modular approach to simulate data from latent factor models. The main function `simulate_factors` serves as the general data simulation function. From the output obtained from `simulate_factors`, other conditions can be added to the data in a piecemeal fashion. Local dependence, for example, can be added to the output using the `add_local_dependence` function. Afterward, population error can be added to the local dependence output using the `add_population_error`. This step-wise procedure aims to build conditions that are useful for testing data under optimal conditons (data output from `simulate_factors`) against more realistic conditions (e.g., substantial cross-loadings, local dependence, population error, wording effects).

## Simulating Latent Factor Data

Examples are provided in the documentation of each function. Below are some of these examples to demonstrate the modular nature of the package. First, generating some data:

```
# Two factor model
simulated <- simulate_factors(
  factors = 2, # number of factors
  variables = 6, # number of variables per factor
  loadings = 0.50 # average loading (+/- 0.10)
  cross_loadings = 0, # magnitude of cross-loadings (value = SD of normal distribution)
  correlations = 0.30, # correlations between factors
  sample_size = 500, # number of cases
  variable_categories = 5 # number of categories in the variables
  # (Inf = continuous; 5 = polytomous; 2 = dichtomous)
)
```

The above example generates a simple structure two factor model with six variables per factor and moderate loadings and correlations between factors. The goal of the package, however, is to provide maximum flexibility and generalization. The `variables` argument, for example, will accept a single value that will be used for all factor or values the same length as the number of factors (specifying the number of items that should belong to each respective factor). Similarly, the `loadings` argument will accept a single value (applied to all factors), values the same length as the number of factors (applied to each respective factor), or loading matrix (pre-specified).

In addition to full specification of parameters, there are arguments for most parameters that allow a `range` to be specified so that the parameter is generated randomly within a certain range. `variables_range`, for instance, will randomly generate the number of variables from a uniform distribution for each factor (e.g., `c(3, 8)` will randomly generate between 3 and 8 variables on each factor).

## Adding Conditions

There are several conditions that can be added on top of the latent factor data generated from `simulate_factors`. Conditions currently available in the package in include adding local dependence (`add_local_dependence`), population error (`add_population_error`), and wording effects (`add_wording_effects`).

EXAMPLES COMING SOON... (see documentation for the respective functions for now)



