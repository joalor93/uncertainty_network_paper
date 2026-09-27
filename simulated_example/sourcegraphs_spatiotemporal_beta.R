dagar_Sigma1 <- function(
    N,
    N_nei,
    nei,
    adjacency_ends,
    rho
) {

  N_nei <- as.numeric(N_nei)
  nei <- as.integer(nei)
  adjacency_ends <- as.integer(adjacency_ends)

  if (length(N_nei) != N) {
    stop("'N_nei' debe tener longitud N.")
  }

  if (length(adjacency_ends) != N) {
    stop("'adjacency_ends' debe tener longitud N.")
  }

  if (abs(rho) >= 1) {
    stop("'rho' debe satisfacer abs(rho) < 1.")
  }

  # Varianza condicional de cada nodo
  vec_var <- (1 - rho^2) /
    (1 + (N_nei - 1) * rho^2)

  # Coeficiente asignado a cada vecino anterior
  b <- rho /
    (1 + (N_nei - 1) * rho^2)

  # Matriz triangular A = I - B
  A <- diag(1, N)

  if (N >= 2L) {
    for (i in 2:N) {

      start_i <- adjacency_ends[i - 1L] + 1L
      end_i   <- adjacency_ends[i]

      if (start_i <= end_i) {
        parents_i <- nei[start_i:end_i]

        A[i, parents_i] <- A[i, parents_i] - b[i]
      }
    }
  }

  # Covarianza DAGAR:
  # Sigma = A^{-1} F^{-1} A^{-T},
  # donde F^{-1} = diag(vec_var)
  A_inv <- solve(A)

  Sigma1 <- A_inv %*%
    diag(vec_var, nrow = N, ncol = N) %*%
    t(A_inv)

  # Corrección numérica de simetría
  Sigma1 <- 0.5 * (Sigma1 + t(Sigma1))

  Sigma1
}

inclattice=function(m){
  n=m^2
  Minc=matrix(0,n,n)
  for(i in 1:(m-1))	for(j in 1:(m-1)) Minc[(i-1)*m+j,(i-1)*m+j+1]=Minc[(i-1)*m+j,i*m+j]=1
  for(i in 1:(m-1)) Minc[(i-1)*m+m,i*m+m]=1
  for(j in 1:(m-1)) Minc[(m-1)*m+j,(m-1)*m+j+1]=1
  Minc+t(Minc)
}



node_as_subtree <- function(g, v) {
  # Ensure names exist (helps with plotting/consistency)
  if (is.null(V(g)$name)) V(g)$name <- as.character(seq_len(vcount(g)))
  
  # Resolve to a vertex id
  v_id <- if (is.numeric(v)) {
    as.integer(v)
  } else {
    match(as.character(v), V(g)$name)
  }
  
  if (is.na(v_id) || v_id < 1 || v_id > vcount(g)) {
    stop("Vertex not found. Check the id/name you passed.")
  }
  
  # Induce a graph on that single vertex
  induced_subgraph(g, vids = v_id)
}


# Return nodes in g1 that are NOT in g2 (by name)
node_diff <- function(g1, g2) {
  n1 <- V(g1)$name
  n2 <- V(g2)$name
  if (is.null(n1) || is.null(n2)) stop("Give your vertices names first: V(g)$name <- as.character(1:vcount(g))")
  setdiff(n1, n2)
}




generating_tree_subtrees = function(subtrees){

  ############# getting a tree from the subtreees ########
  w = sample(1:length(subtrees),length(subtrees))

  # OPTION 2: Or use a function to pick nodes, e.g., highest degree node
  vertex = rbind()
  for(i in 2:length(subtrees)){
    # set.seed(rnorm(1))
    node_a <- names(sample(degree(subtrees[[w[i-1]]]),1))
    node_b <- names(sample(degree(subtrees[[w[i]]]),1))
    vertex = rbind(vertex,c(node_a,node_b))

  }




  # Merge all subtrees into one disconnected graph
  forest <- subtrees[[1]]
  for (i in 2:length(subtrees)) {
    forest <- forest %u% subtrees[[i]]
  }


  forest <- add_edges(forest, vertex)


  return(forest)

}


rowsums_mine <- function(a) {
  rowSums(a)
}


matrixb2 <- function(adjmatinf, rho) {
  
  nneigh <- rowsums_mine(adjmatinf)
  
  b <- rho / (1 + nneigh * rho^2)
  
  B <- adjmatinf
  
  for (i in seq_len(nrow(B))) {
    
    indices <- which(B[i, ] == 1)
    
    if (length(indices) > 0) {
      B[i, indices] <- b[i]
    }
  }
  
  return(B)
}


matrixFinv <- function(nneigh, rho) {
  
  tau <- (1 + (nneigh - 1) * rho^2) /
    (1 - rho^2)
  
  diag(1 / tau)
}


varcovspadagar <- function(adjmatinf, rho) {
  
  nneigh <- rowsums_mine(adjmatinf)
  
  matB <- matrixb2(
    adjmatinf = adjmatinf,
    rho = rho
  )
  
  matFinv <- matrixFinv(
    nneigh = nneigh,
    rho = rho
  )
  
  n <- nrow(matFinv)
  
  I <- diag(n)
  
  I_minus_B <- I - matB
  
  # En C++:
  # inv_I_minus_B = inv(I_minus_B)
  inv_I_minus_B <- solve(I_minus_B)
  
  Sigma <- inv_I_minus_B %*%
    matFinv %*%
    t(inv_I_minus_B)
  
  # Estandarización para convertir Sigma
  # en una matriz de correlación
  diag_Sigma <- diag(Sigma)
  
  inv_sqrt_diag <- 1 / sqrt(diag_Sigma)
  
  D <- outer(
    inv_sqrt_diag,
    inv_sqrt_diag
  )
  
  Sigma <- Sigma * D
  
  return(Sigma)
}


ar1_correlation_matrix_2 <- function(phi, lag) {
  
  outer(
    1:lag,
    1:lag,
    FUN = function(i, j) {
      phi^abs(i - j)
    }
  )
}



spatimecovar_2 <- function(lag,
                           adjmatinf,
                           rho,
                           psi,
                           phi) {
  
  # Correlación espacial DAGAR
  A <- varcovspadagar(
    adjmatinf = adjmatinf,
    rho = rho
  )
  
  # Correlación temporal AR(1)
  B <- ar1_correlation_matrix_2(
    phi = phi,
    lag = lag
  )
  
  # Producto de Kronecker
  # IMPORTANTE:
  # C++ usa arma::kron(A, B)
  timespa <- kronecker(A, B)
  
  # Nugget
  I <- diag(nrow(timespa))
  
  timespa + psi * I
}









generating_tree_subtrees_prob = function(subtrees){

  numvertex = unlist(lapply(subtrees,vcount))
  probtreeind = numvertex/sum(numvertex)

  ############# getting a tree from the subtreees ########
  w = sample(1:length(subtrees),length(subtrees),prob = probtreeind)

  # OPTION 2: Or use a function to pick nodes, e.g., highest degree node
  vertex = rbind()
  for(i in 2:length(subtrees)){
    # set.seed(rnorm(1))
    aux1 = degree(subtrees[[w[i-1]]])
    aux2 = degree(subtrees[[w[i]]])
    probtree1 = aux1/sum(aux1)
    probtree2 = aux2/sum(aux2)
    node_a <- names(sample(aux1,1,probtree1,replace = T))
    node_b <- names(sample(aux2,1,probtree2,replace = T))
    vertex = rbind(vertex,c(node_a,node_b))

  }




  # Merge all subtrees into one disconnected graph
  forest <- subtrees[[1]]
  for (i in 2:length(subtrees)) {
    forest <- forest %u% subtrees[[i]]
  }


  forest <- add_edges(forest, vertex)


  return(forest)

}

library(igraph)

partition_tree <- function(tree, k, min_size = 1, seed = NULL) {
  
  # Validaciones
  if (!inherits(tree, "igraph")) {
    stop("'tree' debe ser un objeto igraph.")
  }
  
  if (is_directed(tree)) {
    stop("'tree' debe ser no dirigido.")
  }
  
  if (!is_tree(tree)) {
    stop("'tree' debe ser un árbol.")
  }
  
  n <- vcount(tree)
  
  if (length(k) != 1L || k != as.integer(k)) {
    stop("'k' debe ser un entero.")
  }
  
  if (length(min_size) != 1L || min_size != as.integer(min_size)) {
    stop("'min_size' debe ser un entero.")
  }
  
  k <- as.integer(k)
  min_size <- as.integer(min_size)
  
  if (k < 1L || k > n) {
    stop("'k' debe satisfacer 1 <= k <= número de nodos.")
  }
  
  if (min_size < 1L) {
    stop("'min_size' debe ser al menos 1.")
  }
  
  if (k * min_size > n) {
    stop(
      "No es posible formar ", k,
      " subárboles con al menos ", min_size,
      " nodos cada uno."
    )
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Caso sin cortes
  if (k == 1L) {
    return(list(subtree_1 = tree))
  }
  
  # Todas las combinaciones de k - 1 aristas
  cut_combinations <- combn(
    seq_len(ecount(tree)),
    k - 1L,
    simplify = FALSE
  )
  
  # Conservar solamente las particiones que cumplen min_size
  valid_partitions <- lapply(cut_combinations, function(edge_ids) {
    
    forest <- delete_edges(tree, edge_ids)
    
    comp <- components(forest)
    
    if (
      comp$no == k &&
      all(comp$csize >= min_size)
    ) {
      return(forest)
    }
    
    NULL
  })
  
  valid_partitions <- Filter(Negate(is.null), valid_partitions)
  
  if (length(valid_partitions) == 0L) {
    stop(
      "No existe una partición válida para k = ", k,
      " y min_size = ", min_size, "."
    )
  }
  
  # Elegir aleatoriamente una partición válida
  forest <- valid_partitions[[
    sample.int(length(valid_partitions), size = 1L)
  ]]
  
  comp <- components(forest)
  
  # Retornar solamente la lista de subárboles
  subtrees <- lapply(seq_len(comp$no), function(j) {
    
    nodes_j <- V(forest)[comp$membership == j]
    
    induced_subgraph(
      forest,
      vids = nodes_j
    )
  })
  
  names(subtrees) <- paste0("subtree_", seq_along(subtrees))
  
  subtrees
}



# ============================================================================
# Funciones para el modelo Stan dagar_beta_spatiotemporal_only
#
# Media:
#   mu[t,i] = beta_st[t,i] * x_st[t,i]
#   beta_st[t,i] = theta0 + theta_s*s_ind[i] + theta_t*d_ind[t]
#                  + theta_st*s_ind[i]*d_ind[t]
#
# Orden esperado de theta_draws:
#   theta[1], theta[2], theta[3], theta[4], rho, phi, sigma, tau2
# donde sigma = sqrt(std_dev_w_2).
# ============================================================================

# Construye la matriz T x N de la media fija beta_st[t,i] * x_st[t,i].
.beta_spatiotemporal_mean <- function(x_st, s_ind, d_ind, theta) {
  x_st <- as.matrix(x_st)
  Tt <- nrow(x_st)
  N  <- ncol(x_st)

  if (length(theta) != 4L)
    stop("'theta' debe tener longitud 4: theta0, theta_s, theta_t, theta_st.")
  if (length(s_ind) != N)
    stop("'s_ind' debe tener longitud ncol(x_st).")
  if (length(d_ind) != Tt)
    stop("'d_ind' debe tener longitud nrow(x_st).")

  s_ind <- as.numeric(s_ind)
  d_ind <- as.numeric(d_ind)

  beta_st <- outer(rep(1, Tt), rep(1, N)) * theta[1] +
    outer(rep(1, Tt), s_ind) * theta[2] +
    outer(d_ind, rep(1, N)) * theta[3] +
    outer(d_ind, s_ind) * theta[4]

  beta_st * x_st
}

# Construye Sigma espacial DAGAR desde el objeto graph usado en las funciones
# originales. Se define aquí para que las nuevas funciones sean autocontenidas.
.dagar_sigma_from_graph <- function(graph, rho, N = NULL) {
  required <- c("N_nei", "nei", "adjacency_ends")
  if (!all(required %in% names(graph))) {
    stop("'graph' debe contener N_nei, nei y adjacency_ends.")
  }

  if (is.null(N)) N <- length(graph$N_nei)
  N_nei <- as.numeric(graph$N_nei)
  nei <- as.integer(graph$nei)
  adjacency_ends <- as.integer(graph$adjacency_ends)

  if (length(N_nei) != N || length(adjacency_ends) != N)
    stop("Las dimensiones de graph no coinciden con N.")

  vec_var <- (1 - rho^2) / (1 + (N_nei - 1) * rho^2)
  b <- rho / (1 + (N_nei - 1) * rho^2)

  A <- diag(1, N)
  if (N >= 2L) {
    for (i in 2:N) {
      start <- adjacency_ends[i - 1L] + 1L
      stop_i <- adjacency_ends[i]
      if (start <= stop_i) {
        parents <- nei[start:stop_i]
        A[i, parents] <- A[i, parents] - b[i]
      }
    }
  }

  Ainv <- solve(A)
  Ainv %*% diag(vec_var, N, N) %*% t(Ainv)
}

# Log-densidad normal multivariada a partir del Cholesky superior de Sigma.
.mvn_logpdf_chol_stbeta <- function(x, mean, cholS) {
  z <- backsolve(cholS, x - mean, transpose = TRUE)
  -0.5 * (length(x) * log(2 * pi) +
            2 * sum(log(diag(cholS))) +
            sum(z^2))
}

# ----------------------------------------------------------------------------
# 1) Análoga a kalman_loglik_by_time
#
# Devuelve log p(y_t | y_1,...,y_{t-1}, parametros) para cada tiempo.
# Integra el campo latente w mediante el filtro de Kalman.
# ----------------------------------------------------------------------------
kalman_loglik_by_time_spatiotemporal_beta <- function(
    y, x_st, s_ind, d_ind, theta,
    rho, phi, sigma2, tau2,
    graph,
    N = ncol(y),
    Tt = nrow(y)) {
  
  y <- as.matrix(y)
  x_st <- as.matrix(x_st)
  
  if (!identical(dim(y), dim(x_st))) {
    stop("'y' y 'x_st' deben tener las mismas dimensiones T x N.")
  }
  
  if (nrow(y) != Tt || ncol(y) != N) {
    stop("N o Tt no coinciden con las dimensiones de 'y'.")
  }
  
  if (!all(c("N_nei", "nei", "adjacency_ends") %in% names(graph))) {
    stop(
      "'graph' debe ser una lista con los componentes: ",
      "'N_nei', 'nei' y 'adjacency_ends'."
    )
  }
  
  if (length(graph$N_nei) != N) {
    stop("La longitud de 'graph$N_nei' debe ser igual a N.")
  }
  
  if (abs(phi) >= 1) {
    stop("'phi' debe estar estrictamente entre -1 y 1.")
  }
  
  if (sigma2 <= 0 || tau2 <= 0) {
    stop("'sigma2' y 'tau2' deben ser positivos.")
  }
  
  mu_fixed <- .beta_spatiotemporal_mean(
    x_st  = x_st,
    s_ind = s_ind,
    d_ind = d_ind,
    theta = theta
  )
  
  Sigma1 <- dagar_Sigma1(
    N              = N,
    N_nei          = graph$N_nei,
    nei            = graph$nei,
    adjacency_ends = graph$adjacency_ends,
    rho            = rho
  )
  
  P0 <- sigma2 * Sigma1
  Q  <- sigma2 * (1 - phi^2) * Sigma1
  R  <- tau2 * diag(N)
  
  loglik_t <- numeric(Tt)
  
  m_prev <- numeric(N)
  C_prev <- P0
  
  for (t in seq_len(Tt)) {
    
    ytilde <- as.numeric(
      y[t, ] - mu_fixed[t, ]
    )
    
    if (t == 1L) {
      a <- numeric(N)
      P <- P0
    } else {
      a <- phi * m_prev
      P <- phi^2 * C_prev + Q
    }
    
    S_t <- P + R
    cholS <- chol(S_t)
    
    loglik_t[t] <- .mvn_logpdf_chol_stbeta(
      ytilde,
      a,
      cholS
    )
    
    innovation <- ytilde - a
    
    Sinv_innovation <- backsolve(
      cholS,
      forwardsolve(
        t(cholS),
        innovation
      )
    )
    
    m_prev <- a + as.numeric(
      P %*% Sinv_innovation
    )
    
    Sinv_P <- backsolve(
      cholS,
      forwardsolve(
        t(cholS),
        P
      )
    )
    
    C_prev <- P - P %*% Sinv_P
    C_prev <- 0.5 * (C_prev + t(C_prev))
  }
  
  loglik_t
}


kalman_loglik_from_draws_spatiotemporal_beta <- function(
    theta_draws,
    y,
    x_st,
    s_ind,
    d_ind,
    graph,
    n_iter = nrow(theta_draws)) {
  
  theta_draws <- as.matrix(theta_draws)
  y <- as.matrix(y)
  x_st <- as.matrix(x_st)
  
  if (ncol(theta_draws) < 8L) {
    stop(
      "'theta_draws' debe tener al menos 8 columnas en el orden: ",
      "theta[1], theta[2], theta[3], theta[4], ",
      "std_dev_w_2, rho, phi, tau2."
    )
  }
  
  if (!identical(dim(y), dim(x_st))) {
    stop("'y' y 'x_st' deben tener las mismas dimensiones T x N.")
  }
  
  N  <- ncol(y)
  Tt <- nrow(y)
  
  if (length(s_ind) != N) {
    stop("'s_ind' debe tener longitud N.")
  }
  
  if (length(d_ind) != Tt) {
    stop("'d_ind' debe tener longitud Tt.")
  }
  
  n_iter <- as.integer(n_iter)
  
  if (length(n_iter) != 1L || is.na(n_iter) || n_iter < 1L) {
    stop("'n_iter' debe ser un entero positivo.")
  }
  
  n_iter <- min(n_iter, nrow(theta_draws))
  
  ll_mat <- matrix(
    NA_real_,
    nrow = n_iter,
    ncol = Tt
  )
  
  for (s in seq_len(n_iter)) {
    
    par <- theta_draws[s, ]
    
    theta  <- par[1:4]
    sigma2 <- par[5]  # std_dev_w_2
    rho    <- par[6]
    phi    <- par[7]
    tau2   <- par[8]
    
    ll_mat[s, ] <- kalman_loglik_by_time_spatiotemporal_beta(
      y      = y,
      x_st   = x_st,
      s_ind  = s_ind,
      d_ind  = d_ind,
      theta  = theta,
      rho    = rho,
      phi    = phi,
      sigma2 = sigma2,
      tau2   = tau2,
      graph  = graph,
      N      = N,
      Tt     = Tt
    )
  }
  
  rownames(ll_mat) <- paste0("draw_", seq_len(n_iter))
  colnames(ll_mat) <- paste0("time_", seq_len(Tt))
  
  ll_mat
}
# ----------------------------------------------------------------------------
# 3) Análoga a loglikes_kron_from_draws_R
#
# Devuelve una matriz draws x (T*N) de log p(y_j | y_-j, parametros),
# con las observaciones vectorizadas por columnas:
#   tiempo 1,...,T del nodo 1; luego tiempo 1,...,T del nodo 2; etc.
#
# Esta es la densidad condicional gaussiana punto a punto calculada desde
# la precisión de la distribución marginal
#   y ~ N(mu, sigma^2 (Sigma_space kron Sigma_time) + tau2 I).
# ----------------------------------------------------------------------------
loglikes_kron_from_draws_spatiotemporal_beta_R <- function(
    theta_draws, x_st, y, s_ind, d_ind,
    graph = NULL, adjmatinf = NULL,
    n_iter = NULL) {

  theta_draws <- as.matrix(theta_draws)
  x_st <- as.matrix(x_st)

  if (is.matrix(y)) {
    y_mat <- y
  } else {
    if (length(y) != length(x_st))
      stop("Si 'y' es vector, su longitud debe ser igual a length(x_st).")
    y_mat <- matrix(y, nrow = nrow(x_st), ncol = ncol(x_st), byrow = FALSE)
  }

  if (!identical(dim(y_mat), dim(x_st)))
    stop("'y' y 'x_st' deben tener las mismas dimensiones T x N.")
  if (ncol(theta_draws) < 8L)
    stop("'theta_draws' debe tener 8 columnas: theta0, theta_s, theta_t, theta_st, rho, phi, sigma, tau2.")

  n_t <- nrow(x_st)
  n_s <- ncol(x_st)
  Ntot <- n_t * n_s
  M <- nrow(theta_draws)
  use_M <- if (is.null(n_iter)) M else min(M, as.integer(n_iter))

  # Compatibilidad con la versión antigua que recibía adjmatinf.
  if (is.null(graph)) {
    if (is.null(adjmatinf))
      stop("Debe proporcionar 'graph' o 'adjmatinf'.")

    adjmatinf <- as.matrix(adjmatinf)
    if (!all(dim(adjmatinf) == c(n_s, n_s)))
      stop("'adjmatinf' debe ser una matriz N x N.")

    # La convención es que los padres del nodo i están en columnas j < i.
    adj_lower <- adjmatinf
    adj_lower[upper.tri(adj_lower, diag = TRUE)] <- 0
    N_nei <- rowSums(adj_lower)
    nei <- unlist(lapply(seq_len(n_s), function(i) which(adj_lower[i, ] == 1)),
                  use.names = FALSE)
    graph <- list(
      N_nei = as.numeric(N_nei),
      nei = as.integer(nei),
      adjacency_ends = as.integer(cumsum(N_nei))
    )
  }

  y_vec <- as.vector(y_mat)
  out <- matrix(NA_real_, nrow = use_M, ncol = Ntot)

  for (m in seq_len(use_M)) {
    theta <- as.numeric(theta_draws[m, 1:4])
    rho   <- as.numeric(theta_draws[m, 5])
    phi   <- as.numeric(theta_draws[m, 6])
    sigma <- as.numeric(theta_draws[m, 7])
    tau2  <- as.numeric(theta_draws[m, 8])

    if (abs(phi) >= 1 || sigma <= 0 || tau2 <= 0) {
      stop("Draw inválido: se requiere |phi| < 1, sigma > 0 y tau2 > 0.")
    }

    mu_mat <- .beta_spatiotemporal_mean(x_st, s_ind, d_ind, theta)
    residual <- y_vec - as.vector(mu_mat)

    Sigma_space <- .dagar_sigma_from_graph(graph, rho, n_s)
    idx_t <- seq_len(n_t)
    Sigma_time <- phi^abs(outer(idx_t, idx_t, "-"))

    eS <- eigen(Sigma_space, symmetric = TRUE)
    eT <- eigen(Sigma_time, symmetric = TRUE)

    lamS <- pmax(eS$values, .Machine$double.eps)
    lamT <- pmax(eT$values, .Machine$double.eps)
    U <- eS$vectors
    V <- eT$vectors

    # Por el orden columna-major: Cov(vec(Y)) = Sigma_space kron Sigma_time.
    inv_eigenvalues <- 1 / (sigma^2 * outer(lamT, lamS, "*") + tau2)

    Rmat <- matrix(residual, nrow = n_t, ncol = n_s, byrow = FALSE)
    transformed <- (t(V) %*% Rmat) %*% U
    transformed <- transformed * inv_eigenvalues
    qmat <- (V %*% transformed) %*% t(U)
    q <- as.vector(qmat)

    # Diagonal de Q sin construir la matriz de precisión completa.
    Qdiag_mat <- ((V^2) %*% inv_eigenvalues) %*% t(U^2)
    Qii <- as.vector(Qdiag_mat)

    cond_var <- 1 / Qii
    cond_mean <- y_vec - q * cond_var

    out[m, ] <- dnorm(y_vec, mean = cond_mean,
                      sd = sqrt(cond_var), log = TRUE)
  }

  rownames(out) <- paste0("draw_", seq_len(use_M))
  colnames(out) <- paste0(
    "t", rep(seq_len(n_t), times = n_s),
    "_s", rep(seq_len(n_s), each = n_t)
  )
  out
}

update_topk <- function(
    top_waic,
    top_adj,
    new_waic,
    new_adj,
    K = 10
) {

  if (length(top_waic) < K) {

    # Todavía hay espacio entre los K mejores
    top_waic <- c(top_waic, new_waic)
    top_adj  <- c(top_adj, list(new_adj))

  } else {

    # Como un WAIC menor es mejor, localiza el peor
    # de los modelos actualmente seleccionados
    worst_in_top <- which.max(top_waic)

    # Sustituye el peor si el nuevo modelo tiene menor WAIC
    if (new_waic < top_waic[worst_in_top]) {
      top_waic[worst_in_top] <- new_waic
      top_adj[[worst_in_top]] <- new_adj
    }
  }

  list(
    top_waic = top_waic,
    top_adj = top_adj
  )
}

# ----------------------------------------------------------------------------
# 4) Versión directa de kalman_loglik_by_time para
#    dagar_beta_spatiotemporal_only.stan
#
# Calcula, para cada tiempo t,
#   log p(y_t | y_1,...,y_{t-1}, theta, rho, phi, sigma, tau2, M_m),
# donde y_t es el vector espacial completo de dimensión N y
#   mu_ti = [theta0 + theta_s s_i + theta_t d_t
#            + theta_st s_i d_t] x_st[t,i].
#
# Esta función es un alias explícito y documentado de la versión adaptada
# anterior, conservando una interfaz cercana a la función original.
# ----------------------------------------------------------------------------
kalman_loglik_by_time_dagar_beta_spatiotemporal <- function(
    y, x_st, s_ind, d_ind, theta,
    rho, phi, sigma, tau2,
    graph, N = ncol(y), Tt = nrow(y)) {

  kalman_loglik_by_time_spatiotemporal_beta(
    y = y,
    x_st = x_st,
    s_ind = s_ind,
    d_ind = d_ind,
    theta = theta,
    rho = rho,
    phi = phi,
    sigma = sigma,
    tau2 = tau2,
    graph = graph,
    N = N,
    Tt = Tt
  )
}

waic_one_graph_kalman_spatiotemporal_beta <- function(
    y,
    x_st,
    s_ind,
    d_ind,
    graph,
    theta_draws,
    S_use = NULL
) {
  
  theta_draws <- as.matrix(theta_draws)
  y <- as.matrix(y)
  x_st <- as.matrix(x_st)
  
  # ------------------------------------------------------------
  # Validaciones
  # ------------------------------------------------------------
  
  if (ncol(theta_draws) < 8L) {
    stop(
      "'theta_draws' debe tener al menos 8 columnas en el orden: ",
      "theta[1], theta[2], theta[3], theta[4], ",
      "std_dev_w_2, rho, phi, tau2."
    )
  }
  
  if (!identical(dim(y), dim(x_st))) {
    stop(
      "'y' y 'x_st' deben tener las mismas dimensiones T x N."
    )
  }
  
  Tt <- nrow(y)
  N  <- ncol(y)
  S  <- nrow(theta_draws)
  
  if (length(s_ind) != N) {
    stop("'s_ind' debe tener longitud N = ncol(y).")
  }
  
  if (length(d_ind) != Tt) {
    stop("'d_ind' debe tener longitud Tt = nrow(y).")
  }
  
  if (!all(
    c("N_nei", "nei", "adjacency_ends") %in% names(graph)
  )) {
    stop(
      "'graph' debe contener los elementos ",
      "'N_nei', 'nei' y 'adjacency_ends'."
    )
  }
  
  if (is.null(S_use)) {
    S_use <- S
  }
  
  S_use <- as.integer(S_use)
  
  if (
    length(S_use) != 1L ||
    is.na(S_use) ||
    S_use < 1L
  ) {
    stop("'S_use' debe ser un entero positivo.")
  }
  
  S_use <- min(S_use, S)
  
  # ------------------------------------------------------------
  # Log-verosimilitud por draw y por tiempo
  #
  # Filas: draws posteriores
  # Columnas: tiempos t = 1,...,T
  # ------------------------------------------------------------
  
  ll_mat <- kalman_loglik_from_draws_spatiotemporal_beta(
    theta_draws = theta_draws,
    y            = y,
    x_st         = x_st,
    s_ind        = s_ind,
    d_ind        = d_ind,
    graph        = graph,
    n_iter       = S_use
  )
  
  # ------------------------------------------------------------
  # WAIC
  # ------------------------------------------------------------
  
  waic_result <- loo::waic(ll_mat)
  
  # Se añade la matriz porque puede ser útil para diagnóstico
  waic_result$log_lik <- ll_mat
  
  waic_result
}





make_block_adjacency <- function(subtrees2) {

  library(igraph)
  library(Matrix)

  if (!is.list(subtrees2)) {
    stop("El objeto debe ser una lista de grafos igraph.")
  }

  if (!all(vapply(subtrees2, inherits, logical(1), what = "igraph"))) {
    stop("Todos los elementos deben ser grafos igraph.")
  }

  if (any(vapply(subtrees2, function(g) is.null(V(g)$name), logical(1)))) {
    stop("Todos los grafos deben tener nombres de nodos.")
  }

  node_names <- as.character(
    sort(
      as.integer(
        unique(
          unlist(lapply(subtrees2, function(g) V(g)$name))
        )
      )
    )
  )

  W <- Matrix(
    0,
    nrow = length(node_names),
    ncol = length(node_names),
    sparse = TRUE,
    dimnames = list(node_names, node_names)
  )

  for (g in subtrees2) {

    nodes_g <- V(g)$name

    Wg <- as_adjacency_matrix(
      g,
      sparse = TRUE,
      attr = NULL
    )

    W[nodes_g, nodes_g] <- Wg[nodes_g, nodes_g, drop = FALSE]
  }

  node_order <- unlist(
    lapply(subtrees2, function(g) {
      as.character(sort(as.integer(V(g)$name)))
    })
  )

  W_block <- W[node_order, node_order, drop = FALSE]

  list(
    W_original_order = W,
    W_block_order = W_block,
    node_order = node_order,
    permutation = data.frame(
      new_position = seq_along(node_order),
      original_node = node_order
    )
  )
}


### produces number of directed neighbors and a vector where each the adjacency of each node ends in the nei vector
adj_nei=function(Minc){
  Minc.low=Minc
  Minc.low[upper.tri(Minc)]=0
  return(list(rowSums(Minc.low),Minc = Minc.low))
}


adj_index=function(Minc){
  Minc.low=Minc
  Minc.low[upper.tri(Minc)]=0
  #list(N_nei=rowSums(Minc.low),adj_index=cumsum(rowSums(Minc.low)))
 return (cumsum(rowSums(Minc.low)))
}

neighbors=function(Minc){
  n=nrow(Minc)
  unlist(lapply(2:n,function(i) which(Minc[i,1:(i-1)]==1)))
}


mvn_logpdf_chol <- function(y, mu, cholSigma) {

  n <- length(y)

  z <- y - mu

  # chol() en R devuelve U tal que:
  # t(U) %*% U = Sigma

  tmp <- forwardsolve(
    t(cholSigma),
    z
  )

  quad <- sum(tmp^2)

  logdet <- 2 * sum(
    log(diag(cholSigma))
  )

  -0.5 * (
    n * log(2 * pi) +
      logdet +
      quad
  )
}

kalman_loglik_by_time <- function(
    y, Xt_list, beta, rho, phi, sigma, tau2,
    graph, N, Tt) {

  Sigma1 <- dagar_Sigma1(
    N,
    graph$N_nei,
    graph$nei,
    graph$adjacency_ends,
    rho
  )

  P0 <- sigma^2 * Sigma1
  Q  <- sigma^2 * (1 - phi^2) * Sigma1
  R  <- tau2 * diag(N)

  loglik_t <- numeric(Tt)

  m_prev <- rep(0.0, N)
  C_prev <- P0

  for (t in 1:Tt) {

    Xt <- Xt_list[[t]]

    ytilde <- as.numeric(
      y[t, ] - Xt %*% beta
    )

    if (t == 1) {

      a <- rep(0.0, N)
      P <- P0

    } else {

      a <- phi * m_prev
      P <- phi^2 * C_prev + Q
    }

    S <- P + R

    cholS <- chol(S)

    loglik_t[t] <- mvn_logpdf_chol(
      ytilde,
      a,
      cholS
    )

    v <- ytilde - a

    Sinv_v <- backsolve(
      cholS,
      forwardsolve(
        t(cholS),
        v
      )
    )

    m_new <- a + as.numeric(
      P %*% Sinv_v
    )

    Sinv_P <- backsolve(
      cholS,
      forwardsolve(
        t(cholS),
        P
      )
    )

    C_new <- P - P %*% Sinv_P

    # útil para evitar pequeñas asimetrías numéricas
    C_new <- 0.5 * (C_new + t(C_new))

    m_prev <- m_new
    C_prev <- C_new
  }

  loglik_t
}

kalman_loglik_from_draws <- function(
    theta_draws,
    y,
    Xt_list,
    graph,
    n_iter = nrow(theta_draws)) {
  
  theta_draws <- as.matrix(theta_draws)
  y <- as.matrix(y)
  
  if (ncol(theta_draws) < 5L) {
    stop(
      "'theta_draws' debe tener al menos 5 columnas en el orden: ",
      "beta[1], ..., beta[p], rho, phi, sigma, tau2."
    )
  }
  
  N  <- ncol(y)
  Tt <- nrow(y)
  
  if (!is.list(Xt_list)) {
    stop("'Xt_list' debe ser una lista.")
  }
  
  if (length(Xt_list) != Tt) {
    stop("'Xt_list' debe tener longitud Tt.")
  }
  
  # Número de covariables/betas inferido desde Xt_list
  p <- ncol(Xt_list[[1]])
  
  if (is.null(p)) {
    stop("Cada elemento de 'Xt_list' debe ser una matriz N x p.")
  }
  
  # Verificar dimensiones de todas las matrices X_t
  for (t in seq_len(Tt)) {
    
    Xt <- as.matrix(Xt_list[[t]])
    
    if (!identical(dim(Xt), c(N, p))) {
      stop(
        "Cada Xt_list[[t]] debe tener dimensión N x p. ",
        "Problema en t = ", t, "."
      )
    }
  }
  
  # Se necesitan p betas + rho + phi + sigma + tau2
  if (ncol(theta_draws) < p + 4L) {
    stop(
      "'theta_draws' debe tener al menos p + 4 columnas en el orden: ",
      "beta[1], ..., beta[p], rho, phi, sigma, tau2."
    )
  }
  
  n_iter <- as.integer(n_iter)
  
  if (length(n_iter) != 1L || is.na(n_iter) || n_iter < 1L) {
    stop("'n_iter' debe ser un entero positivo.")
  }
  
  n_iter <- min(n_iter, nrow(theta_draws))
  
  ll_mat <- matrix(
    NA_real_,
    nrow = n_iter,
    ncol = Tt
  )
  
  for (s in seq_len(n_iter)) {
    
    par <- theta_draws[s, ]
    
    beta  <- par[seq_len(p)]
    rho   <- par[p + 1L]
    phi   <- par[p + 2L]
    sigma <- par[p + 3L]
    tau2  <- par[p + 4L]
    
    ll_mat[s, ] <- kalman_loglik_by_time(
      y       = y,
      Xt_list = Xt_list,
      beta    = beta,
      rho     = rho,
      phi     = phi,
      sigma   = sigma,
      tau2    = tau2,
      graph   = graph,
      N       = N,
      Tt      = Tt
    )
  }
  
  rownames(ll_mat) <- paste0("draw_", seq_len(n_iter))
  colnames(ll_mat) <- paste0("time_", seq_len(Tt))
  
  ll_mat
}

