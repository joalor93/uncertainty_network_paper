README - Simulated Example
==========================

This folder contains a complete simulated example illustrating the graph-selection and model-averaging procedures used in the spatiotemporal analysis.

FILES
-----

selecting_trees_ram_friendly.R
This script generates the simulated dataset and performs the candidate-tree selection step.

Main tasks:
- Loads the true spatial tree and a predefined partition of the graph.
- Simulates a spatiotemporal Gaussian response using the specified covariance structure and regression parameters.
- Splits the data into an estimation period (first 20 time points) and a forecasting period (last 7 time points).
- Fits the spatiotemporal DAGAR model in Stan using the true graph to obtain posterior parameter draws.
- Generates a collection of candidate spatial trees.
- Evaluates each candidate tree using WAIC through the Kalman-based likelihood calculation.
- Retains the best candidate trees according to WAIC and also stores the worst tree as a benchmark.
- Saves the simulated response, covariates, selected candidate trees, WAIC values, and forecasting information in sim_001.rds.

The script is designed to be memory efficient by keeping only the posterior quantities required for graph evaluation.

stacking_and_best_tree.R
This script uses the output from selecting_trees_ram_friendly.R to fit and compare the selected candidate graph models.

Main tasks:
- Loads sim_001.rds and the candidate trees selected in the previous step.
- Fits the spatiotemporal model separately for the selected candidate graphs using Stan.
- Obtains posterior predictive distributions for the 7 forecasting time points.
- Computes predictive log-likelihoods for each fitted candidate graph.
- Evaluates the candidate models using PSIS-LOO and ELPD.
- Computes stacking and pseudo-BMA model weights.
- Constructs model-averaged posterior predictive distributions.
- Identifies the best candidate graph according to predictive performance.
- Fits additional benchmark models using the worst graph and the true graph.
- Compares stacking, pseudo-BMA, the best candidate graph, the worst graph, and the true graph using predictive accuracy and uncertainty measures.
- Reports RMSE, ELPD, predictive interval length, and coverage indicators.

WORKFLOW
--------

The scripts should be run in the following order:

1. selecting_trees_ram_friendly.R
   Produces sim_001.rds containing the simulated data and the selected candidate trees.

2. stacking_and_best_tree.R
   Uses sim_001.rds to fit the candidate graph models, compute model-averaging weights, and evaluate predictive performance.

MAIN SUPPORTING FILES
---------------------

The scripts also require the following supporting files:

- sourcegraphs_spatiotemporal_beta.R
  Auxiliary functions for graph manipulation, covariance construction, WAIC computation, and Kalman-based predictive likelihood evaluation.

- truetree_2.rds
  True spatial tree used to generate the simulated data.

- partition45.rds
  Graph partition used to generate candidate trees.

- dagar_beta_spatiotemporal_with_intercept.stan
  Stan model used during the candidate-tree selection stage.

- dagar_beta_spatiotemporal_only.stan
  Stan model used for fitting the selected candidate graphs and generating forecasts.

OUTPUT
------

The main intermediate output is:

- sim_001.rds
  Contains the simulated training/test data, covariates, selected candidate trees, their WAIC values, the worst candidate tree, and related information.

The final script produces predictive summaries comparing:
- Stacking
- Pseudo-BMA
- Best selected tree
- Worst tree
- True tree

using RMSE, ELPD, predictive interval length, and coverage.
