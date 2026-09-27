README - results_application.R
==============================

File: results_application.R

Description
-----------
This script summarizes the predictive analysis for the application study.

It loads the fitted models corresponding to different spatial graph assumptions, including the reported graph, a random graph, a model without spatial dependence, and a collection of candidate graph models used for model averaging.

The script performs the following main tasks:

1. Computes PSIS-LOO and ELPD values for the candidate spatial graph models.
2. Calculates stacking and pseudo-BMA model weights.
3. Builds model-averaged posterior predictive distributions.
4. Identifies the best candidate graph according to ELPD.
5. Compares predictive performance across the reported graph, random graph, best candidate graph, stacking, pseudo-BMA, and the non-spatial model.
6. Evaluates forecasting performance using:
   - ELPD
   - RMSE
   - 95% predictive interval coverage
   - Mean predictive interval length
7. Produces pairwise PSIS-LOO comparisons between the best candidate graph and the alternative approaches.
8. Generates spatial maps for:
   - Best-graph predictions
   - Stacking predictions
   - Observed values
   - Prediction errors
9. Saves the resulting figures in EPS format.

Main input files
----------------
- modelo_grafo_2_intercepto_selec_random_2.rds
- modelo_grafo_2_intercepto_selec_ind_temp_2.rds
- modelo_grafo_2_intercepto_selec_reported_2.rds
- modelo_grafo_2_for_averaging.rds
- datapred.rds
- chile_shp_union.shp and its associated shapefile components

Required R packages
-------------------
loo
data.table
dplyr
igraph
ggplot2
sf

Output
------
The script produces predictive performance summaries and spatial figures comparing the different graph-based modeling strategies and model-averaging approaches.
