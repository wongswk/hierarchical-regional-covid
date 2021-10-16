# hierarchical-regional-covid
Hierarchical Bayesian model for assessing regional COVID-19 physical distancing measures

This repository contains the code to reproduce the results of our paper on regional COVID-19 modelling in British Columbia: https://arxiv.org/abs/2104.10878

Our hierarchical Bayesian multi-region model extends the modeling framework presented in the paper
```
Anderson SC, Edwards AM, Yerlanov M, Mulberry N, Stockdale JE, et al. (2020) Quantifying the impact of COVID-19 control measures using a Bayesian model of physical distancing. PLOS Computational Biology 16(12): e1008274. https://doi.org/10.1371/journal.pcbi.1008274
```
Portions of our code are adapted from these authors' original Github at https://github.com/carolinecolijn/distancing-impact-covid19 as permitted under the GPL-3.0 License.

### Data 

The input data files for each region are stored inside `data-generated`, e.g., `data-generated/daily-cases-Fraser.csv`, which contain the case counts for each day.
Data for the five British Columbia regions were obtained from
```
BC Centre for Disease Control (2021). BC COVID-19 Data. http://www.bccdc.ca/health-info/diseases-conditions/covid-19/data
```

### Usage

Analyses are carried out in R, using the `rstan` package to facilitate the fitting of the main model.

The following R commands call the various scripts provided to do model fitting. The fitted model R objects are saved inside `data-generated`.
```
# Main model: Fit hierarchical multiregion model for the BC data 
source("analysis/run-multiregion.R")

# Provincial-wide model: Fit a single model to the combined BC data 
source("analysis/run-provincewide.R")

# Simulation study: Fit hierarchical multiregion model to simulated data with known parameters
source("analysis/run-multiregion-simulated.R")

# Main model (prediction): Fit hierarchical multiregion model for the BC data and generate predictions for November
source("analysis/run-multiregion-predNov.R")

# Provincial-wide model (prediction): Fit a single model to the combined BC data and generate predictions for November
source("analysis/run-multiregion-predNov.R")

# Individual models: Fit a model to each BC region separately
source("analysis/run-individual-region.R")

```

After the models are fitted, plots and estimates can be generated via the following scripts. The figures are saved inside `figs-ms`.
```
# Make figures for multiregion fit with comparison to provincial-wide fit
source("analysis/multiregion-makeplots.R")

# Calculate R0 values implied by the model parameters
source("analysis/calc-R0.R")

# Make figures for November prediction using multiregion and provincial-wide fits
source("analysis/multiregion-predictionNov-makeplots.R")

# Make figures for models fitted to individual regions with comparison to provincial-wide fit
source("analysis/individual-region-makeplots.R")
```

The parameter values for the simulation study can be adjusted by modifying the values in `analysis/simulation-params.R`.



