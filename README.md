# HB_HABs
Hierarchical Bayesian model for forecasting harmful algal blooms in Scotland

These models aim to predict HABs and their associated toxins. There are currently four models (Hierarchical Bayesian cumulative ordinal, Hierarchical Bayesian dual logistic, Random Forest, and Random Forest split) along with an ensemble model. Each uses seasonally-varying effects, including autoregressive (2) terms, environmental conditions, particle tracking simulations of offshore influx, with site-level random effects. 

The ultimate goal is to produce an operational model that generates forecasts automatically based on the new monitoring counts and the WeStCOMS model outputs as they become available.

This involves the following steps:

## Set up
### code/01_initial_setup.R
  1. Create particle tracking directories, properties file
  2. Run particle tracking simulations
  3. Compile particle tracking output
  4. Extract WeStCOMS variables
  5. Extract Copernicus variables
  6. Extract algal densities
  7. Extract toxin concentrations
  8. Load traffic light thresholds
  9. Compile datasets for each species
### code/02-1_fit_initial_models.R
  1. Load compiled datasets
  2. Define predictors
  3. for(year in years) { fit(data[-year]); predict(data[year])}
### (code/02-2_variable_selection.R)
  1. Load compiled datasets
  2. Define sets of predictors
  3. for(preds in predSets) { for(year in years) { ... } }
### code/03_calc_model_weights.R
  1. Load predictions
  2. Calculate overall predictive likelihood
  3. Calculate per-site predictive likelihood
  4. Calculate per-month predictive likelihood

## Weekly updates
### code/04_weekly_main.R
  - runs each of the following scripts in succession
### code/04-1_get_new_data.R
  1. Set dates
  2. Download WeStCOMS (days -7:5)
  3. Download Copernicus (days -7:5)
  4. Download new algal densities
  5. Download new toxin concentrations
### code/04-2_simulate_particles.R
  1. Create particle tracking directories, properties file
  2. Copy previous particle locations file
  3. Run particle tracking simulations with hot start
### code/04-3_compile_data.R
  1. Compile particle tracking output
  2. Extract and compile WeStCOMS, Copernicus
  3. Append new days to previous dataset, replacing with new data as needed
  4. Compile datasets for each species
### code/04-4_fit_new_models.R
  1. Load compiled datasets
  2. Define predictors
  3. Fit models with observed data
  4. Predict probabilities for the next week
### code/04-5_export_predictions.R
  1. Load predictions
  2. Calculate ensemble predictions
  3. Clean and summarise
  4. Export to csv, map, etc 


# Models
We compare several Hierarchical Bayesian and Machine Learning models, as well as different methods for generating ensemble predictions.

  - **HB~ord~**: Cumulative ordinal with traffic light warning categories, post-hoc aggregated to bloom states based on regulatory thresholds  
  - **HB~2lgs~**: Dual logistic model predicting bloom state, split to submodels by Bloom~t-1~  
  - **ML~RF~**: Random forest for bloom state  
  - **ML~2RF~**: Random forest for bloom state, split to submodels by Bloom~t-1~  
  - **EN~avg~**: Ensemble using a simple average of predictions from all models
  - **EN~mo~**: Ensemble using LL ratios as model weights, averaging by month  
  - **EN~xy~**: Ensemble using LL ratios as model weights, averaging by site  
  - **EN~xy-mo**: Ensemble using LL ratios as model weights, averaging by site-month  
  - **EN-d~avg~**: Dynamically updated version of **EN~avg~**  
  - **EN-d~mo~**: Dynamically updated version of **EN~mo~**  
  - **EN-d~xy~**: Dynamically updated version of **EN~xy~**  
  - **EN-d~xy-mo~**: Dynamically updated version of **EN~xy-mo~**  
  
Static ensembles models weight based on the initial model fits only, maintaining those weights for all predictions. Dynamic ensemble models recalculate weights with each new set of predictions.  

I don't know if it's possible to generate a pseudo-posterior distribution for RF models, given the mean prediction and an error... If so, ensembles could be fully Bayesian, including uncertainty. 