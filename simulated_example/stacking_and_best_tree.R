############################################################
# Required packages
############################################################

library(RcppArmadillo)
library(mvtnorm)
library(parallel)
library(foreach)
library(doParallel)
library(doSNOW)
library(igraph)
library(invgamma)
library(modeest)
library(rbenchmark)
library(Matrix)
library(TreeSampleR)
library(loo)
library(rstan)


############################################################
# Working directory and model specification
############################################################

address = "C:/Users/alejo/Dropbox/postdoc20242/random_graphs_new/scripts_random_graphs_prof_mauricio_ordenados/new_simulations/Codigo github biom/simulated_example"
setwd(address)

# Auxiliary functions used throughout the spatiotemporal analysis
source("sourcegraphs_spatiotemporal_beta.R")

# True spatial tree used to generate the simulated data
tree <- readRDS("truetree_2.rds")

# Compile the Stan model only once, before parallel computation
stan_code_obs <- readLines("dagar_beta_spatiotemporal_only.stan")
model_obs <- stan_model(model_code = stan_code_obs)


############################################################
# Spatial and temporal structure
############################################################

# Adjacency matrix associated with the true spatial graph
adj_matrix <- as_adjacency_matrix(tree, sparse = FALSE)
adj_matrix = as.matrix(adj_matrix)

N_spa = nrow(adj_matrix)

# Total number of temporal observations
N_time = 27

# Keep only one triangular part of the adjacency matrix.
# This avoids counting each undirected edge twice when the
# graph structure is passed to the DAGAR model.
adjm1 = adj_matrix
adjm1[upper.tri(adjm1)] = 0

# First 20 time points are used for model estimation,
# whereas the remaining 7 observations are used for forecasting.
N_time_obs = 20
pred = N_time - N_time_obs

nobs = N_spa * N_time_obs


############################################################
# MCMC and model settings
############################################################

iter = 10000
burn = 1000
thin = 20

thetaini = c(0.1, 0.5, 0.2)

region <- 1:N_spa
time   <- 1:N_time

# Spatial and temporal indicators used in the regression component
s_ind <- as.integer(region > 6)
d_ind <- as.integer(time > 10)

time_obs = 1:N_time_obs

# Prevent nested parallelization from BLAS/OpenMP while several
# Stan models are already being fitted in parallel.
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)


############################################################
# Load one simulated dataset
############################################################

aux1 = readRDS("sim_001.rds")
indipred = aux1$indipred


############################################################
# Prior specification
############################################################

mualpha <- 0
sigmaalpha <- 100

mutheta = rep(0, 4)
sigmatheta = rep(100, 4)

# Quantities that will be retained from the Stan fit.
# In particular, y_fore contains posterior predictive draws,
# while the remaining quantities are required to reconstruct
# the predictive log-likelihood.
pars_keep <- c(
  "y_fore",
  "log_lik_real",
  "theta",
  "std_dev_w_2",
  "rho",
  "phi",
  "tau2"
)


############################################################
# Covariates for the forecasting period
############################################################

x_fore = aux1$xtot

# Alternative representation of forecasting covariates,
# with one vector for each time point.
x_fore_stan <- lapply(
  seq_len(nrow(x_fore)),
  function(t) as.vector(x_fore[t, ])
)

# Seven future time points are predicted
L_fore = 7
time_fore <- 21:27


############################################################
# Fit candidate graph models in parallel
############################################################

# Ten candidate spatial graphs are fitted independently.
# Since each Stan fit is computationally expensive, the models
# are distributed across different CPU cores.
n.cores <- 10
my.cluster <- makeCluster(n.cores)

itersim = 10

registerDoSNOW(my.cluster)


foreach(
  i = 1:itersim,
  .multicombine = TRUE,
  .packages = c(
    "parallel",
    "doParallel",
    "Matrix",
    "rstan"
  )
) %dopar% {
  
  ##########################################################
  # Load simulated response and candidate graph
  ##########################################################
  
  aux = readRDS("sim_001.rds")
  
  yobs  = aux$youtcome$yobs
  yreal = aux$youtcome$yreal
  
  # Observed response arranged as time x space
  ymatobs = matrix(yobs, N_time_obs, N_spa)
  
  # True future values used only for evaluating forecasting performance
  ymatreal = matrix(yreal, pred, N_spa)
  
  xmatobs = aux$xobs
  
  # Each parallel iteration corresponds to a different
  # candidate spatial tree.
  adj_matrix_rand = aux$candidatetrees[[i]]
  
  
  ##########################################################
  # Convert candidate graph into the format required by Stan
  ##########################################################
  
  # Retain only one direction of each undirected edge
  adj_matrix_rand[upper.tri(adj_matrix_rand)] = 0
  
  N_edges <- sum(adj_matrix_rand)
  
  # Vector containing the neighbors of each spatial unit
  nei <- unlist(
    lapply(seq_len(nrow(adj_matrix_rand)), function(j) {
      which(adj_matrix_rand[j, ] == 1)
    }),
    use.names = FALSE
  )
  
  # Edge indices and number of neighbors for each area
  adj.ends <- adj_index(adj_matrix_rand)
  N_nei <- rowSums(adj_matrix_rand)
  
  
  ##########################################################
  # Data passed to the Stan model
  ##########################################################
  
  data_model <- list(
    N = N_spa,
    N_time = N_time_obs,
    s_ind = s_ind,
    time = time_obs,
    
    N_nei = N_nei,
    N_edges = N_edges,
    nei = nei,
    adjacency_ends = adj.ends,
    
    y = ymatobs,
    x_st = xmatobs,
    
    mutheta = mutheta,
    sigmatheta = sigmatheta,
    mualpha = mualpha,
    sigalpha = sigmaalpha,
    
    # Forecasting is performed jointly with model fitting
    do_forecast = 1,
    L_fore = 7,
    x_st_fore = x_fore,
    time_fore = time_fore,
    
    # Future observations are supplied for predictive calculations
    y_real = ymatreal
  )
  
  
  ##########################################################
  # Bayesian model fitting
  ##########################################################
  
  fit <- sampling(
    model_obs,
    data = data_model,
    iter = 10000,
    warmup = 3000,
    chains = 1,
    thin = 10,
    save_warmup = FALSE,
    pars = pars_keep,
    include = TRUE,
    refresh = 0
  )
  
  
  ##########################################################
  # Extract posterior predictive draws and model parameters
  ##########################################################
  
  # Posterior draws for the seven forecasting periods
  yfore_chosen <- as.matrix(
    fit,
    pars = "y_fore"
  )
  
  # Posterior draws required to reconstruct the
  # spatiotemporal covariance structure.
  theta_draws <- as.matrix(
    fit,
    pars = c(
      "theta",
      "std_dev_w_2",
      "rho",
      "phi",
      "tau2"
    )
  )
  
  
  ##########################################################
  # Graph information used in predictive likelihood
  ##########################################################
  
  graph_rand = list(
    N_nei = N_nei,
    nei = nei,
    adjacency_ends = adj.ends
  )
  
  
  ##########################################################
  # Compute pointwise predictive log-likelihood
  ##########################################################
  
  # The Kalman-based function evaluates the predictive
  # log-likelihood for posterior draws of the parameters.
  #
  # Only 500 posterior draws are used here to reduce the
  # computational cost of the LOO calculations.
  ll_500 <- kalman_loglik_from_draws_spatiotemporal_beta(
    theta_draws = theta_draws,
    y            = ymatobs,
    x_st         = xmatobs,
    s_ind        = s_ind,
    d_ind        = time_obs,
    graph        = graph_rand,
    n_iter       = 500
  )
  
  
  ##########################################################
  # Save results for this candidate graph
  ##########################################################
  
  outcome_2 = list(
    pred_dist = yfore_chosen,
    log_lik = ll_500
  )
  
  saveRDS(
    outcome_2,
    paste0(
      "chains_simulation_stan_selected_",
      i,
      "_",
      ".rds"
    )
  )
  
  return(outcome_2)
}


# Release the parallel workers after all candidate models
# have been fitted.
stopCluster(my.cluster)
gc()


####################################################################
# Combine results from the candidate spatial graphs
####################################################################

# Read the posterior predictive distributions and log-likelihood
# matrices obtained from the ten candidate graphs.
listaux = list()

for (i in 1:10) {
  
  outcome3 = readRDS(
    paste0(
      "chains_simulation_stan_selected_",
      i,
      "_.rds"
    )
  )
  
  listaux[[i]] = list(
    pred_dist = outcome3$pred_dist,
    log_lik = outcome3$log_lik
  )
}

stackpred = listaux


####################################################################
# Compute model-specific PSIS-LOO quantities
####################################################################

prop_concordantes = 0
itersim = 10

h = stackpred

aux = readRDS("sim_001.rds")

yobs  = aux$youtcome$yobs
yreal = aux$youtcome$yreal

yobsmat = matrix(
  yobs,
  N_time_obs,
  N_spa
)

ymatreal = matrix(
  yreal,
  pred,
  N_spa
)

xmatobs = aux$xobs

# Stores the loo objects associated with each candidate graph
arg_weights = list()

mse = 0
elpdind = 0
se = 0


# Evaluate each candidate model separately using PSIS-LOO.
# The resulting loo objects will subsequently be used to obtain
# stacking and pseudo-BMA model weights.
for (k in 1:length(h)) {
  
  arg_weights[[k]] = loo(
    h[[k]]$log_lik
  )
  
  elpdind[k] =
    arg_weights[[k]]$estimates[
      "elpd_loo",
      "Estimate"
    ]
}


####################################################################
# Ranking of candidate graph structures
####################################################################

# Sort candidate graphs according to their ELPD.
# Larger ELPD values correspond to better predictive performance.
elpdsort = sort(
  elpdind,
  decreasing = TRUE
)

# Rank based on the original candidate ordering
ranking_1 <- rank(
  -elpdind,
  ties.method = "average"
)

# Rank after sorting the ELPD values
ranking_2 <- rank(
  -elpdsort,
  ties.method = "average"
)

# Kendall's tau is used to measure agreement between rankings.
tau_kendall <- cor(
  ranking_1,
  ranking_2,
  method = "kendall",
  use = "complete.obs"
)

# When there are no ties, Kendall's tau can be converted into
# the proportion of concordant pairs using this expression.
prop_concordantes <- (1 + tau_kendall) / 2


# Identify the candidate graphs with the three and five
# highest individual ELPD values.
pos3 = which(
  elpdind %in% elpdsort[1:3]
)

pos5 = which(
  elpdind %in% elpdsort[1:5]
)


####################################################################
# Model averaging weights
####################################################################

# Bayesian stacking chooses weights by maximizing the
# leave-one-out predictive performance of the model combination.
model_weights_1 = loo_model_weights(
  arg_weights
)

# Pseudo-BMA provides an alternative set of model weights
# based on the models' LOO predictive performance.
model_weights_2 = loo_model_weights(
  arg_weights,
  method = "pseudobma"
)

gc()

weights = as.numeric(model_weights_1)
weights_2 = as.numeric(model_weights_2)

gc()


####################################################################
# Construct model-averaged predictive distributions
####################################################################

ypreddistfin = 0
ypreddistfin2 = 0
ypreddistfin3 = 0

elpd1aux = 0
elpd3aux = 0

elpdind = 0


for (s in 1:length(weights)) {
  
  # Stacking predictive distribution
  ypreddistfin =
    ypreddistfin +
    weights[s] * h[[s]]$pred_dist
  
  # Pseudo-BMA predictive distribution
  ypreddistfin3 =
    ypreddistfin3 +
    weights_2[s] * h[[s]]$pred_dist
  
  # Corresponding weighted log-likelihood matrices
  elpd1aux =
    elpd1aux +
    weights[s] * h[[s]]$log_lik
  
  elpd3aux =
    elpd3aux +
    weights_2[s] * h[[s]]$log_lik
}


####################################################################
# Predictive performance of model averaging
####################################################################

elpd1aux2 = loo(elpd1aux)
elpd3aux2 = loo(elpd3aux)

elpd1 =
  elpd1aux2$estimates[
    "elpd_loo",
    "Estimate"
  ]

elpd3 =
  elpd3aux2$estimates[
    "elpd_loo",
    "Estimate"
  ]


# Posterior medians of the model-averaged predictive distributions
medianypred = apply(
  ypreddistfin,
  2,
  median
)

medianypred3 = apply(
  ypreddistfin3,
  2,
  median
)


# Out-of-sample mean squared prediction error
msefindist =
  mean(
    (medianypred - yreal)^2
  )

msefindist3 =
  mean(
    (medianypred3 - yreal)^2
  )


####################################################################
# Fit the model associated with the worst graph
####################################################################

# The worst graph is stored separately in the simulated dataset.
# It provides a benchmark showing how much predictive performance
# deteriorates when an unfavorable spatial structure is assumed.
adj_matrix_worse = aux$worsttree

adj_matrix_worse[
  upper.tri(adj_matrix_worse)
] = 0

N_edges_worse <- sum(
  adj_matrix_worse
)

nei_worse <- neighbors(
  adj_matrix_worse
)

adj.ends_worse <- adj_index(
  adj_matrix_worse
)

N_nei_worse <- rowSums(
  adj_matrix_worse
)


# Data for the Stan model under the worst spatial graph
data_model_worse <- list(
  
  N = N_spa,
  N_time = N_time_obs,
  s_ind = s_ind,
  time = time_obs,
  
  N_nei = N_nei_worse,
  N_edges = N_edges_worse,
  nei = nei_worse,
  adjacency_ends = adj.ends_worse,
  
  y = yobsmat,
  x_st = xmatobs,
  
  mutheta = mutheta,
  sigmatheta = sigmatheta,
  
  mualpha = mualpha,
  sigalpha = sigmaalpha,
  
  do_forecast = 1,
  
  L_fore = 7,
  x_st_fore = x_fore,
  time_fore = time_fore,
  
  y_real = ymatreal
)


fit_worse <- sampling(
  
  model_obs,
  data = data_model_worse,
  
  iter = 10000,
  warmup = 3000,
  chains = 1,
  
  thin = 10,
  
  save_warmup = FALSE,
  
  pars = pars_keep,
  include = TRUE,
  
  refresh = 0
)


# Posterior predictive draws under the worst graph
yfore_worse <- as.matrix(
  fit_worse,
  pars = "y_fore"
)


theta_draws_worse <- as.matrix(
  fit_worse,
  pars = c(
    "theta",
    "std_dev_w_2",
    "rho",
    "phi",
    "tau2"
  )
)


graph_worse = list(
  N_nei = N_nei_worse,
  nei = nei_worse,
  adjacency_ends = adj.ends_worse
)


# Predictive log-likelihood under the worst graph
ll_500_worse <-
  kalman_loglik_from_draws_spatiotemporal_beta(
    
    theta_draws = theta_draws_worse,
    
    y = yobsmat,
    x_st = xmatobs,
    
    s_ind = s_ind,
    d_ind = time_obs,
    
    graph = graph_worse,
    
    n_iter = 500
  )


# Posterior median prediction under the worst graph
predmcmcR2worse = apply(
  yfore_worse,
  2,
  quantile,
  probs = 0.5
)


# Add the worst graph to the collection of fitted models.
# Candidate graphs occupy positions 1,...,10, so position 11
# is reserved for this additional benchmark model.
chainrand = 11

h[[chainrand]] = list(
  pred = predmcmcR2worse,
  pred_dist = yfore_worse,
  log_lik = ll_500_worse
)


####################################################################
# Performance of the worst graph
####################################################################

ypredonedist =
  h[[chainrand]]$pred_dist

mseone =
  mean(
    (
      h[[chainrand]]$pred -
        yreal
    )^2
  )

ypredonemean = apply(
  ypredonedist,
  2,
  median
)

mseonedist =
  mean(
    (
      ypredonemean -
        yreal
    )^2
  )

elpdoneaux = loo(
  h[[chainrand]]$log_lik
)

elpdone =
  elpdoneaux$estimates[
    "elpd_loo",
    "Estimate"
  ]


####################################################################
# Performance of the best candidate graph
####################################################################

# Best graph among the original candidate structures according
# to its individual PSIS-LOO ELPD.
chainbest =
  which(
    elpdind == max(elpdind)
  )

ypredbestdist =
  h[[chainbest]]$pred_dist

msebest =
  mean(
    (
      h[[chainbest]]$pred -
        yreal
    )^2
  )


ypredbestmean = apply(
  ypredbestdist,
  2,
  median
)


msebest =
  mean(
    (
      ypredbestmean -
        yreal
    )^2
  )


elpdbestaux =
  loo(
    h[[chainbest]]$log_lik
  )

elpdbest =
  elpdbestaux$estimates[
    "elpd_loo",
    "Estimate"
  ]


####################################################################
# Predictive intervals
####################################################################

# 95% posterior predictive intervals for stacking
ypreddistint = apply(
  ypreddistfin,
  2,
  quantile,
  probs = c(0.025, 0.975)
)

# 95% posterior predictive intervals for the worst graph
ypredoneint = apply(
  ypredonedist,
  2,
  quantile,
  probs = c(0.025, 0.975)
)

# 95% posterior predictive intervals for the best candidate graph
ypredbestint = apply(
  ypredbestdist,
  2,
  quantile,
  probs = c(0.025, 0.975)
)

# 95% posterior predictive intervals for pseudo-BMA
ypreddistint3 = apply(
  ypreddistfin3,
  2,
  quantile,
  probs = c(0.025, 0.975)
)


####################################################################
# Coverage probabilities
####################################################################

# Indicator equal to one when the true future observation falls
# inside the corresponding 95% predictive interval.

intervalindstack =
  0 + (
    yreal >= ypreddistint[1, ] &
      yreal <= ypreddistint[2, ]
  )

intervalindone =
  0 + (
    yreal >= ypredoneint[1, ] &
      yreal <= ypredoneint[2, ]
  )

intervalindbest =
  0 + (
    yreal >= ypredbestint[1, ] &
      yreal <= ypredbestint[2, ]
  )

intervalindbma =
  0 + (
    yreal >= ypreddistint3[1, ] &
      yreal <= ypreddistint3[2, ]
  )


####################################################################
# Predictive interval lengths
####################################################################

# Median interval length is used to evaluate the sharpness of
# each forecasting approach.
lengthstack =
  median(
    ypreddistint[2, ] -
      ypreddistint[1, ]
  )

lengthbma =
  median(
    ypreddistint3[2, ] -
      ypreddistint3[1, ]
  )

lengthone =
  median(
    ypredoneint[2, ] -
      ypredoneint[1, ]
  )

lengthbest =
  median(
    ypredbestint[2, ] -
      ypredbestint[1, ]
  )


####################################################################
# Fit the model using the true spatial graph
####################################################################

# The true graph provides an oracle benchmark. It represents
# the predictive performance that can be obtained when the
# actual spatial dependence structure is known.
adj_matrix_true = adjm1

adj_matrix_true[
  upper.tri(adj_matrix_true)
] = 0

N_edges_true <-
  sum(adj_matrix_true)

nei_true <-
  neighbors(adj_matrix_true)

adj.ends_true <-
  adj_index(adj_matrix_true)

N_nei_true <-
  rowSums(adj_matrix_true)


# Data passed to Stan under the true graph
data_model_true <- list(
  
  N = N_spa,
  N_time = N_time_obs,
  
  s_ind = s_ind,
  time = time_obs,
  
  N_nei = N_nei_true,
  N_edges = N_edges_true,
  nei = nei_true,
  
  adjacency_ends = adj.ends_true,
  
  y = yobsmat,
  x_st = xmatobs,
  
  mutheta = mutheta,
  sigmatheta = sigmatheta,
  
  mualpha = mualpha,
  sigalpha = sigmaalpha,
  
  do_forecast = 1,
  
  L_fore = 7,
  x_st_fore = x_fore,
  time_fore = time_fore,
  
  y_real = ymatreal
)


fit_true <- sampling(
  
  model_obs,
  data = data_model_true,
  
  iter = 10000,
  warmup = 3000,
  chains = 1,
  
  thin = 10,
  
  save_warmup = FALSE,
  
  pars = pars_keep,
  include = TRUE,
  
  refresh = 0
)


####################################################################
# Predictions and likelihood under the true graph
####################################################################

yfore_true <- as.matrix(
  fit_true,
  pars = "y_fore"
)


theta_draws_true <- as.matrix(
  fit_true,
  pars = c(
    "theta",
    "std_dev_w_2",
    "rho",
    "phi",
    "tau2"
  )
)


graph_true = list(
  N_nei = N_nei_true,
  nei = nei_true,
  adjacency_ends = adj.ends_true
)


ll_500_true <-
  kalman_loglik_from_draws_spatiotemporal_beta(
    
    theta_draws = theta_draws_true,
    
    y = yobsmat,
    x_st = xmatobs,
    
    s_ind = s_ind,
    d_ind = time_obs,
    
    graph = graph_true,
    
    n_iter = 500
  )


####################################################################
# Oracle predictive performance
####################################################################

# Posterior predictive 2.5%, 50%, and 97.5% quantiles
predmcmcR2 = apply(
  yfore_true,
  2,
  quantile,
  probs = c(
    0.025,
    0.5,
    0.975
  )
)


# Coverage indicator for the true-graph model
intervalindtrue =
  0 + (
    yreal >= predmcmcR2[1, ] &
      yreal <= predmcmcR2[3, ]
  )


# Prediction error based on posterior median forecasts
msetrue =
  mean(
    (
      predmcmcR2[2, ] -
        yreal
    )^2
  )


# Median width of the 95% predictive intervals
lengthtrue =
  median(
    predmcmcR2[3, ] -
      predmcmcR2[1, ]
  )


# PSIS-LOO ELPD under the true graph
elpdtrueaux =
  loo(ll_500_true)

elpdtrue =
  elpdtrueaux$estimates[
    "elpd_loo",
    "Estimate"
  ]


####################################################################
# Collect performance measures
####################################################################

# The final vector summarizes predictive accuracy, predictive
# performance according to ELPD, and interval sharpness for:
#
#   - Stacking
#   - Pseudo-BMA
#   - True graph
#   - Best candidate graph
#   - Worst graph
#
# RMSE is reported instead of MSE by taking the square root.

resinte = c(
  
  sqrt(
    c(
      msefindist,
      msefindist3,
      msetrue,
      msebest,
      mseonedist
    )
  ),
  
  elpd1,
  elpd3,
  elpdone,
  elpdbest,
  elpdtrue,
  
  lengthstack,
  lengthbma,
  lengthone,
  lengthbest,
  lengthtrue
)


names(resinte) = c(
  
  "stackdist",
  "bmadist",
  "truemse",
  "msebest",
  "mseone",
  
  "elpdstack",
  "elpdbma",
  "elpdone",
  "elpdbest",
  "elpdtrue",
  
  "lengthstack",
  "lengthbma",
  "lengthone",
  "lengthbest",
  "lengthtrue"
)

resinte
####################################################################
# Coverage indicators
####################################################################

intervalindstack
intervalindone
intervalindbest
intervalindbma