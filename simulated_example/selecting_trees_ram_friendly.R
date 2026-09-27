library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)
library(Matrix)
library(loo)
library(rstan)
library(igraph)


address = "C:/Users/alejo/Dropbox/postdoc20242/random_graphs_new/scripts_random_graphs_prof_mauricio_ordenados/new_simulations/Codigo github biom/simulated_example"

setwd(address)

### e SOURCE FUNCTIONS 
source("sourcegraphs_spatiotemporal_beta.R")
### STAN CODE 
stan_code_obs <- readLines("dagar_beta_spatiotemporal_with_intercept.stan")
model_obs <- stan_model(model_code = stan_code_obs)


rstan_options(auto_write = TRUE)
options(mc.cores = 1)

#### True tree and simulated partitions ##############
tree <- readRDS("truetree_2.rds")
subtrees2 <- readRDS("partition45.rds")

# adjacency matrix formed from the subtrees 
W_block <- make_block_adjacency(subtrees2)
W_block_matrix <- as.matrix(W_block$W_original_order)

#### true adjacency matrix
adj_matrix <- as.matrix(igraph::as_adjacency_matrix(tree, sparse = FALSE))


### Number of spatial and temporal observations
N_spa <- nrow(adj_matrix)
N_time <- 27 ### total temporal observations
N_time_obs <- 20 ### used for estimation 
pred <- N_time - N_time_obs ### used as test dataset 


adjm1 <- adj_matrix
adjm1[upper.tri(adjm1)] <- 0


##### covariates and parameters ##########
sigma2 <- 1.5
theta_true <- c(0.5, 0.5, 0.7, 0.3)
rho_true <- 0.8
phi_true <- 0.8
tau2_true <- 0.5
psi <- tau2_true / sigma2

N <- N_spa * N_time
xvar <-  rnorm(N)
xvarmat = matrix(xvar,N_time,N_spa)
region <- 1:N_spa
time  = 1:N_time
s_ind <- as.integer(region > 6)
time_obs = 1:N_time_obs



#### simulating the mean of the process
mumat = matrix(0,N_time,N_spa)
for (t in 1:N_time) {
  for (i in 1:N_spa) {
mumat[t,i] <- (theta_true[1] + theta_true[2]*s_ind[i] + theta_true[3]*time[t] + theta_true[4]*s_ind[i]*time[t])*xvarmat[t,i]
  }
}

mu = as.vector(mumat)
  
nneigh <- rowSums(adjm1)

#### covariate matrix of the process 
covmat_cpp <- sigma2 * spatimecovar_2(
  lag = N_time, adjmatinf = adjm1, rho = rho_true, psi = psi, phi = phi_true
)

mualpha <- 0
sigmaalpha <- 100

mutheta = rep(0,4)
sigmatheta = rep(100,4)


#### parameters to keep (RAM friendly), in this step predictions are not necessary 
pars_keep <- c("theta","std_dev_w_2","rho","phi","tau2")



### the final list contains ##### 

## youtcome: The simulated outcome Y
## indipred: the prediction indicator
## xobs: the observed covariates matrix X 
## candidatetrees : the selected candidate trees under the waic criterion
## candidates_waic: the waic criterion for those trees 
## worsttree: the worst tree over all the trees 
## waic_val: waic for all the graphs


  #simulated responde
  ytot <- as.vector(mvtnorm::rmvnorm(1, mean = mu, sigma = covmat_cpp))
  y1mat <- matrix(ytot, N_time, N_spa)
  y1matreal = y1mat[(N_time - pred + 1):N_time, ]
  yvecreal = as.vector(y1matreal)
  y1mat[(N_time - pred + 1):N_time, ] <- NA
  
  
  yvec <- as.vector(y1mat)
  miss <- is.na(yvec)

  # prediction indicator
  indipred <- as.integer(miss)     # 1 if predicted, 0 if observed
  
  yreal <- yvecreal[miss] #### test response
  ## observed x and Y (training set)
  yobs <- yvec[!miss]
  # matriz observada T_obs x N_spa (para WAIC)
  ymatobs <- matrix(yobs, N_time_obs, N_spa)
  
  
  xtot <- xvarmat 
  xreal = xtot[miss] #### test set of covariates 
  xobs = xtot[!miss]### training set of covariates 
  x1mat <- matrix(xtot, N_time, N_spa)
  xtotreal = x1mat[(N_time - pred + 1):N_time, ]
  x1mat[(N_time - pred + 1):N_time, ] <- NA
  xobsmat = matrix(xobs,N_time_obs,N_spa)
  #adjmatknown <- W_block_matrix
  #adjmatknown[upper.tri(adjmatknown)] <- 0
  youtcome <- list(yobs = yobs, yreal = as.vector(y1matreal))
  
  youtcome <- list(yobs = yobs, yreal = as.vector(y1matreal),xobs = xobs,xreal =xreal)
  

  N_edges <- sum(adjm1)
  nei <- neighbors(adjm1)
  adj.ends <- adj_index(adjm1)
  N_nei <- rowSums(adjm1)

  x_fore <- array(
    numeric(0),
    dim = c(0, N_spa)
  )
  
  time_fore = integer(0)
  ##### income to fit the model with the known list of subtrees in stan ####
  data_model <- list(
    N = N_spa, N_time = N_time_obs,s_ind = s_ind, time = time_obs,
    N_nei = N_nei, N_edges = N_edges, nei = nei, adjacency_ends = adj.ends,
    y = ymatobs, x = xobsmat, mutheta = mutheta, sigmatheta = sigmatheta,
    mualpha = mualpha, sigalpha = sigmaalpha,
    do_forecast = 0L,
    L_fore = 0L,
    x_fore = x_fore,
    time_fore = time_fore
  )


  # Fit Stan
  fit <- sampling(
    model_obs, data = data_model,
    iter = 10000, warmup = 1000, chains = 1,
    thin = 10,
    save_warmup = FALSE,
    pars = pars_keep, include = TRUE,
    refresh = 0
  )

 # Extracting the parameters in a RAM friendly way 
  theta_cols <- paste0("theta[", 1:4, "]")
  draws <- as.matrix(fit, pars = c(theta_cols, "rho", "phi", "std_dev_w_2", "tau2"))

  theta_draws <- cbind(
    draws[, theta_cols, drop = FALSE],
    sigma = sqrt(draws[, "std_dev_w_2"]),
    rho   = draws[, "rho"],
    phi   = draws[, "phi"],
    tau2  = draws[, "tau2"]
  )

  rm(fit, draws); gc()

  # Getting the candidate trees, we keep the best  Kkeep
  Kkeep <- 10 ## the best candidate trees 
  Q = 50 ### set of candidate trees 
  top_waic <- numeric(0)
  top_adj  <- list()
  worst_waic <- -Inf
  worst_adj  <- NULL

  waic_val <- numeric(Q)

  #### computing the waic for 1000 trees
  for (i in 1:Q) {
    treerand <- generating_tree_subtrees(subtrees2)
    treerand <- igraph::simplify(
      treerand,
      remove.multiple = TRUE,
      remove.loops = TRUE
    )

    Adj <- as.matrix(igraph::as_adjacency_matrix(treerand, sparse = FALSE))
    Adj <- Adj[order(as.numeric(rownames(Adj))), order(as.numeric(colnames(Adj)))]
    Adj[upper.tri(Adj)] <- 0

    N_edges_r <- sum(Adj)
    
    #Adj_lower <- Adj
    #Adj_lower[upper.tri(Adj_lower, diag = TRUE)] <- 0
    
    N_nei_r <- rowSums(Adj)
    
    nei_r <- unlist(
      lapply(seq_len(nrow(Adj)), function(j) {
        which(Adj[j, ] == 1)
      }),
      use.names = FALSE
    )
    
    
    
    adj.ends_r <- adj_index(Adj)
    N_nei_r <- rowSums(Adj)

    graph_r <- list(N_nei = N_nei_r, nei = nei_r, adjacency_ends = adj.ends_r)

    grafo = list(N_nei = N_nei_r,nei = nei_r, adjacency_ends = adj.ends_r)
    a = waic_one_graph_kalman_spatiotemporal_beta(
    y = ymatobs, x_st = xobsmat, s_ind = s_ind, d_ind = time_obs, graph= grafo, theta_draws,
    S_use = NULL)
    

    w <- a$estimates["waic", "Estimate"]
    waic_val[i] <- w

    upd <- update_topk(top_waic, top_adj, w, Adj, K = Kkeep)
    top_waic <- upd$top_waic
    top_adj  <- upd$top_adj

    if (w > worst_waic) {
      worst_waic <- w
      worst_adj  <- Adj
    }

    rm(treerand, Adj, graph_r, a, N_edges_r, nei_r, adj.ends_r, N_nei_r)
    if (i %% 50 == 0) gc()
    print(i)
  }

  ### gettin the positions of the best WAIC
  ord <- order(top_waic)

  ### creating the output
  sim_out <- list(
    youtcome = youtcome,
    indipred = indipred,
    xobs = xobsmat,
    xtot = xtotreal, 
    candidatetrees  = top_adj[ord],
    candidates_waic = top_waic[ord],
    worsttree       = worst_adj,
    waic_val        = waic_val
  )
  
##### waic for all the trees and for the best ten selected

  sim_out$candidates_waic
  sim_out$waic_val
  
  path <- file.path(sprintf("sim_%03d.rds", 1))
  ### saving the output (this one i am gonna use it to feed algorithms such as stacking and the best tree procedure)
  saveRDS(sim_out, path)
