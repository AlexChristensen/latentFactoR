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

### Local Dependence

Local dependence is the extent to which variables are correlated after accounting for latent variables in a model. Adding local dependence to simualted data can be achieved using the following code:

```
# Add local dependence
two_factor_LD <- add_local_dependence(
  lf_object = simulated, # object from `simulate_factors`
  proportion_LD = 0.25, # proportion of variables to be locally dependent on *each* factor
  add_residuals = 0.20, # magnitude of the residuals to add between locally dependent variables
  allow_multiple = FALSE # whether a variable can be locally dependent with more than one other variable
)
```

### Population Error
Population error is when there is misfit between the population and sample model. This misfit is often characterized by small local dependencies between variables that are not enough to change the factor loadings but result in model misspecification (e.g., over- or underfactoring). Adding population error to simulated data can be achieved using the following code:

```
# Add small population error using Cudeck method
two_factor_Cudeck <- add_population_error(
  lf_object = simulated, # object from `simulate_factors`
  cfa_method = "minres", # minimum residual method for determining model fit
  fit = "rmsr", # measure used to determine misfit
  misfit = "close", # amount of misfit (can be "close", "acceptable", or a numeric value between 0 and 1)
  error_method = "cudeck" # method to add misfit ("cudeck" and "yuan" are available; "cudeck" is recommend)
)
```

### Wording Effects
Wording effects occur for a number of reasons: participants generally agree with statements (`"acquiescence"`), have difficulty discerning level due to negating words (`"difficulty"`), don't care and respond randomly (`"random_careless"`), or select the same response for everything (`"straight_line"`). These different wording effects can wreak havoc on the accurate detection of dimensions, precision of parameter estimates, and (mis)specification of models. Adding wording effects to simulated data can be achieved using the following code:

```
# Add wording effects using acquiescence method
two_factor_acquiescence <- add_wording_effects(
  lf_object = simulated, # object from `simulate_factors`
  proportion_negative = 0.50, # proportion of negatively worded variables on each factor
  proportion_biased_cases = 0.10, # proportion of total cases that are biased
  proportion_biased_variables = 1, # proportion of variables that show bias (out of possible variables to have bias)
  proportion_biased_person = 1, # person-specific parameter of how much bias they show (proportion of variables to show bias)
  method = "acquiescence" # method to introduce bias
)
```
