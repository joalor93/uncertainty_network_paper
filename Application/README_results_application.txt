README - results_application.R
==============================

File
----
results_application.R

Overview
--------
This script performs the final predictive comparison and spatial visualization for the application study under uncertainty in the spatial graph structure.

The analysis combines results from several fitted spatiotemporal models, including:
- the reported/original spatial graph,
- a randomly selected spatial graph,
- a model without spatial dependence,
- and a set of selected candidate spatial graphs.

Main tasks
----------
1. Reads the fitted candidate models and combines their posterior predictive draws and pointwise log-likelihood values.
2. Computes PSIS-LOO and ELPD for each candidate graph.
3. Estimates model-averaging weights using:
   - Bayesian stacking,
   - pseudo-BMA.
4. Constructs model-averaged posterior predictive distributions.
5. Identifies the best candidate graph according to ELPD.
6. Compares the competing approaches using:
   - ELPD,
   - RMSE,
   - empirical coverage of 95% predictive intervals,
   - average predictive interval length.
7. Performs pairwise PSIS-LOO comparisons between the best candidate graph and the alternative approaches.
8. Merges the predictive results with the spatial shapefile.
9. Produces maps of:
   - predictions from the best candidate graph,
   - stacking predictions,
   - observed values,
   - spatial prediction errors.
10. Saves the spatial figures in EPS format.

Candidate-model files
---------------------
The script expects ten fitted candidate-model files named:

modelo_grafo_2_selec_2_1.rds
...
modelo_grafo_2_selec_2_10.rds

Each file must contain at least:
- pred_dist: posterior predictive draws,
- log_lik: pointwise log-likelihood values.

Benchmark-model files
---------------------
modelo_grafo_2_intercepto_selec_reported_2.rds
Model fitted using the reported/original spatial graph.

modelo_grafo_2_intercepto_selec_random_2.rds
Model fitted using a randomly selected spatial graph.

modelo_grafo_2_intercepto_selec_ind_temp_2.rds
Model without spatial dependence, retaining the temporal structure.

Additional input files
----------------------
datapred.rds
Contains the observed values and prediction-related variables for the forecasting period.

chile_shp_union.shp
Spatial polygons used to generate the application maps. The associated .dbf, .shx, and .prj files should be kept in the same directory.

Required R packages
-------------------
loo
data.table
dplyr
igraph
ggplot2
sf

Main outputs
------------
The script produces:
- a summary table containing ELPD, RMSE, coverage, and predictive interval length for the competing approaches,
- pairwise PSIS-LOO comparisons,
- faceted spatial maps,
- individual date-specific maps saved in EPS format.

Compared approaches
-------------------
The final comparison includes:
- Reported graph
- Random graph
- Best candidate graph
- Stacking
- Pseudo-BMA
- No-spatial-dependence model

Purpose
-------
The script is intended to evaluate how uncertainty in the spatial network affects out-of-sample prediction and to compare graph selection with model averaging in the application study.
