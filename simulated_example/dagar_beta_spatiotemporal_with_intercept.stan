functions {
  vector rowsums_mine(matrix a) {
    int N = rows(a);
    vector[N] row_sums;
    for (i in 1:N) row_sums[i] = sum(a[i]);
    return row_sums;
  }

  // These are kept for compatibility with your file, but NOT used in the scalable conditional prior below.
  matrix matrixb2(matrix adjmatinf, real rho) {
    int N = rows(adjmatinf);
    vector[N] nneigh;
    matrix[N, N] result;
    for (i in 1:N) nneigh[i] = sum(adjmatinf[i]);
    vector[N] b = rho ./ (1 + nneigh .* square(rho));
    result = adjmatinf;
    for (i in 1:N)
      for (j in 1:N)
        if (result[i, j] == 1) result[i, j] = b[i];
    return result;
  }

  matrix matrixFinv(vector nneigh, real rho) {
    int N = num_elements(nneigh);
    vector[N] tau;
    for (i in 1:N)
      tau[i] = (1 + (nneigh[i] - 1) * square(rho)) / (1 - square(rho));
    return diag_matrix(1 ./ tau);
  }

  matrix matrixF(vector nneigh, real rho) {
    int N = num_elements(nneigh);
    vector[N] tau;
    for (i in 1:N)
      tau[i] = (1 + (nneigh[i] - 1) * square(rho)) / (1 - square(rho));
    return diag_matrix(tau);
  }

  matrix matrixFinvsqrt(vector nneigh, real rho) {
    int N = num_elements(nneigh);
    vector[N] tau;
    for (i in 1:N)
      tau[i] = (1 + (nneigh[i] - 1) * square(rho)) / (1 - square(rho));
    return diag_matrix(1 ./ sqrt(tau));
  }

  matrix varcovspadagar(matrix adjmatinf, real rho) {
    int N = rows(adjmatinf);
    vector[N] nneigh = rowsums_mine(adjmatinf);
    matrix[N, N] matB    = matrixb2(adjmatinf, rho);
    matrix[N, N] matFinv = matrixFinv(nneigh, rho);
    matrix[N, N] I = diag_matrix(rep_vector(1.0, N));
    matrix[N, N] I_minus_B = I - matB;
    matrix[N, N] inv_I_minus_B = inverse(I_minus_B);
    return inv_I_minus_B * matFinv * inv_I_minus_B';
  }
}

data {
  int<lower=1> N;                  // spatial nodes (must be in DAG order)
  int<lower=1> N_time;             // observed time points

  vector[N] N_nei;                 // number of parents per node
  int<lower=0> N_edges;
  array[N_edges] int<lower=1, upper=N> nei;
  array[N] int<lower=0, upper=N_edges> adjacency_ends;

  int<lower=0, upper=1> do_forecast;
  int<lower=0> L_fore;

  matrix[N_time, N] y;             // y[t, i]
  matrix[N_time, N] x;             // x[t, i]

  // Spatial and temporal indicators defining beta_{ti}
  array[N] int<lower=0, upper=1> s_ind;
  array[N_time] int <lower=0> time;

  // Future covariate and temporal indicator
  array[L_fore] vector[N] x_fore;
  array[L_fore] int<lower=0, upper=1> time_fore;

  // Prior hyperparameters
  real mualpha;
  real<lower=0> sigalpha;
  vector[4] mutheta;
  vector<lower=0>[4] sigmatheta;
}

parameters {
  real<lower=0, upper=1> rho;      // spatial DAGAR parameter
  real<lower=0, upper=1> phiaux;   // auxiliary parameter for AR(1)
  real<lower=0> std_dev_w_2;       // sigma^2, latent variance
  real<lower=0> tau2;              // nugget variance

  real alpha;
  vector[4] theta;                 // theta0, theta_s, theta_t, theta_st
  matrix[N_time, N] w;             // latent spatio-temporal field
}

transformed parameters {
  matrix[N_time, N] beta_st;
  matrix[N_time, N] mufix;

  // map [0,1] -> [-1,1]
  real<lower=-1, upper=1> phi = -1 + 2 * phiaux;
  real<lower=0> std_dev_w = sqrt(std_dev_w_2);

  for (t in 1:N_time) {
    for (i in 1:N) {
      beta_st[t, i] = theta[1]
                      + theta[2] * s_ind[i]
                      + theta[3] * time[t]
                      + theta[4] * s_ind[i] * time[t];

      // mu_{ti} = alpha + beta_{ti} x_{ti}
      mufix[t, i] = alpha + beta_st[t, i] * x[t, i];
    }
  }
}

model {
  vector[N] b;
  vector[N] vec_var;

  // ---- priors ----
  alpha       ~ normal(mualpha, sigalpha);
  theta       ~ normal(mutheta, sigmatheta);
  // rho      ~ beta(3, 2);
  // phiaux   ~ beta(3, 2);
  // tau2     ~ inv_gamma(2, 2);
  target += -log(std_dev_w_2);
  // ---- DAGAR weights from parent counts ----
  // vec_var[i] = (1-rho^2) / (1 + (N_nei[i]-1)*rho^2)
  // b[i]       = rho       / (1 + (N_nei[i]-1)*rho^2)
  vec_var = (1 - square(rho)) ./ (1 + (N_nei - rep_vector(1.0, N)) * square(rho));
  b       =  rho              ./ (1 + (N_nei - rep_vector(1.0, N)) * square(rho));

  // ============================================================
  // Latent prior: separable DAGAR (space) ⊗ AR(1) (time)
  //
  // w[1, ]            ~ N(0, sigma^2 * Sigma1)  in DAGAR conditional form
  // eta[t, ] = w[t,]-phi*w[t-1,] ~ N(0, sigma^2*(1-phi^2)*Sigma1) in DAGAR conditional form
  // ============================================================

  // ---- t = 1: DAGAR conditional for w[1, ] ----
  // Root node (must have zero parents in your DAG order)
  w[1, 1] ~ normal(0, std_dev_w);

  for (i in 2:N) {
    int start = adjacency_ends[i - 1] + 1;
    int stop  = adjacency_ends[i];
    real acc  = 0;

    if (start <= stop) {
      for (j in start:stop)
        acc += w[1, nei[j]];
    }

    w[1, i] ~ normal(b[i] * acc, std_dev_w * sqrt(vec_var[i]));
  }

  // ---- t >= 2: DAGAR conditional for innovations eta[t, ] ----
  for (t in 2:N_time) {

    // root innovation
    (w[t, 1] - phi * w[t - 1, 1])
      ~ normal(0, std_dev_w * sqrt(1 - square(phi)));

    for (i in 2:N) {
      int start = adjacency_ends[i - 1] + 1;
      int stop  = adjacency_ends[i];
      real acc  = 0;

      if (start <= stop) {
        for (j in start:stop) {
          // innovation at parent node
          acc += (w[t, nei[j]] - phi * w[t - 1, nei[j]]);
        }
      }

      (w[t, i] - phi * w[t - 1, i])
        ~ normal(b[i] * acc,
                 std_dev_w * sqrt((1 - square(phi)) * vec_var[i]));
    }
  }

  // ---- Observation model ----
  // y_{ti} ~ N(alpha + beta_{ti} x_{ti} + w_{ti}, tau2)
  for (i in 1:N)
    y[, i] ~ normal(to_vector(mufix[, i]) + to_vector(w[, i]), sqrt(tau2));
}

generated quantities {
  matrix[L_fore, N] w_fore;
  matrix[L_fore, N] eta_fore;
  matrix[L_fore, N] y_fore;

  if (do_forecast == 1) {

    // Recalcular b y vec_var (no se guardaron del model block)
    vector[N] b;
    vector[N] vec_var;

    vec_var = (1 - square(rho)) ./ (1 + (N_nei - rep_vector(1.0, N)) * square(rho));
    b       =  rho              ./ (1 + (N_nei - rep_vector(1.0, N)) * square(rho));

    // último estado observado del campo latente
    vector[N] w_prev = to_vector(w[N_time, ]);

    for (ell in 1:L_fore) {
      vector[N] eta;
      vector[N] w_new;

      // ============================================================
      // 1) Simular eta ~ N(0, std_dev_w^2*(1-phi^2)*Sigma_DAGAR)
      //    usando la forma condicional DAGAR (recursiva en el DAG)
      // ============================================================

      // raíz del DAG
      eta[1] = normal_rng(0, std_dev_w * sqrt(1 - square(phi)));

      // nodos restantes (condicional a "padres" en el DAG)
      for (i in 2:N) {
        int start = adjacency_ends[i - 1] + 1;
        int stop  = adjacency_ends[i];
        real acc  = 0;

        if (start <= stop) {
          for (j in start:stop)
            acc += eta[nei[j]];
        }

        eta[i] = normal_rng(b[i] * acc,
                            std_dev_w * sqrt((1 - square(phi)) * vec_var[i]));
      }

      eta_fore[ell, ] = to_row_vector(eta);

      // ============================================================
      // 2) Propagar AR(1): w_{T+ell} = phi * w_{T+ell-1} + eta_{ell}
      // ============================================================
      w_new = phi * w_prev + eta;
      w_fore[ell, ] = to_row_vector(w_new);

      // ============================================================
      // 3) Simular y_fore dado w_new y X_fore
      // ============================================================
      {
        real sigma = sqrt(tau2);
        vector[N] beta_fore;
        vector[N] mu_fore;

        for (i in 1:N) {
          beta_fore[i] = theta[1]
                           + theta[2] * s_ind[i]
                           + theta[3] * time_fore[ell]
                           + theta[4] * s_ind[i] * time_fore[ell];

          mu_fore[i] = alpha + beta_fore[i] * x_fore[ell][i] + w_new[i];
          y_fore[ell, i] = normal_rng(mu_fore[i], sigma);
        }
      }

      // actualizar para el siguiente paso
      w_prev = w_new;
    }

  } else {
    w_fore   = rep_matrix(0.0, L_fore, N);
    eta_fore = rep_matrix(0.0, L_fore, N);
    y_fore   = rep_matrix(0.0, L_fore, N);
  }
}
