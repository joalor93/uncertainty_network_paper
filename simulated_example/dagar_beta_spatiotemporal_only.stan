functions {
  vector rowsums_mine(matrix a) {
    int N = rows(a);
    vector[N] row_sums;
    for (i in 1:N) row_sums[i] = sum(a[i]);
    return row_sums;
  }

  // Kept for compatibility (NOT used in scalable conditional prior below)
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
  int<lower=1> N;                 // spatial nodes (MUST be in DAG order)
  int<lower=1> N_time;            // time points
  vector[N] N_nei;                // number of parents per node (in DAG order)
  int<lower=0> N_edges;
  int<lower=1, upper=N> nei[N_edges];
  int<lower=0, upper=N_edges> adjacency_ends[N]; // cumulative ends per i

  int<lower=0, upper=1> do_forecast;

  // Permitir L_fore = 0 cuando do_forecast = 0
  int<lower=0> L_fore;

  matrix[L_fore, N] y_real;
  matrix[N_time, N] y;            // y[t, i]

  // Covariable whose coefficient varies by spatial group and temporal period
  matrix[N_time, N] x_st;         // x_st[t, i]
  array[L_fore] vector[N] x_st_fore;

  // Binary indicators used in beta_st[t,i]
  array[N] int<lower=0, upper=1> s_ind;             // spatial group
  array[N_time] int<lower=0> time;        // temporal period (training)
  array[L_fore] int<lower=0> time_fore;   // temporal period (forecast)

  // Priors for theta = (theta_0, theta_s, theta_t, theta_st)
  vector[4] mutheta;
  vector<lower=0>[4] sigmatheta;
}

parameters {
  real<lower=0, upper=1> rho;          // spatial DAGAR parameter
  real<lower=0, upper=1> phiaux;       // auxiliary for AR(1)
  real<lower=0> std_dev_w_2;           // sigma^2 (latent scale squared)
  real<lower=0> tau2;                  // nugget variance

  vector[4] theta;                 // theta_0, theta_s, theta_t, theta_st
  matrix[N_time, N] w;                 // latent field
}

transformed parameters {
  matrix[N_time, N] beta_st;
  matrix[N_time, N] mu;

  // beta_st[t,i] = theta_0 + theta_s s_i + theta_t d_t
  //                + theta_st s_i d_t
  for (t in 1:N_time) {
    for (i in 1:N) {
      beta_st[t, i] = theta[1]
                      + theta[2] * s_ind[i]
                      + theta[3] * time[t]
                      + theta[4] * s_ind[i] * time[t];

      mu[t, i] = beta_st[t, i] * x_st[t, i] + w[t, i];
    }
  }

  // map [0,1] -> [-1,1]
  real<lower=-1, upper=1> phi = -1 + 2 * phiaux;

  // sigma
  real<lower=0> std_dev_w = sqrt(std_dev_w_2);
}

model {
  vector[N] b;
  vector[N] vec_var;

  // ---- priors ----
  target += -log(std_dev_w_2);
  theta ~ normal(mutheta, sigmatheta);

  // ---- DAGAR weights from parent counts ----
  vec_var = (1 - square(rho)) ./ (1 + (N_nei - rep_vector(1.0, N)) * square(rho));
  b       =  rho              ./ (1 + (N_nei - rep_vector(1.0, N)) * square(rho));

  // ============================================================
  // Latent prior: separable DAGAR (space) ⊗ AR(1) (time)
  // ============================================================

  // ---- t = 1: DAGAR conditional for w[1, ] ----
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
          acc += (w[t, nei[j]] - phi * w[t - 1, nei[j]]);
        }
      }

      (w[t, i] - phi * w[t - 1, i])
        ~ normal(b[i] * acc,
                 std_dev_w * sqrt((1 - square(phi)) * vec_var[i]));
    }
  }

  // ---- Observation model: nugget tau2 * I ----
  for (i in 1:N)
    y[, i] ~ normal(to_vector(mu[, i]), sqrt(tau2));
}

generated quantities {
  // -------------------------------
  // (A) log-likelihood per obs (train)
  // -------------------------------
  matrix[N_time, N] log_lik;

  {
    real sigma = sqrt(tau2);

    for (i in 1:N) {
      vector[N_time] mu_i = to_vector(mu[, i]);

      for (t in 1:N_time) {
        log_lik[t, i] = normal_lpdf(y[t, i] | mu_i[t], sigma);
      }
    }
  }

  // -------------------------------
  // (B) Forecast + test log predictive density per draw
  // -------------------------------
  matrix[L_fore, N] w_fore;
  matrix[L_fore, N] eta_fore;
  matrix[L_fore, N] y_fore;

  // log p(y_real | theta, w_fore) por celda (ell,i)
  matrix[L_fore, N] log_lik_real;

  // suma total en este draw (opcional)
  real log_lik_real_sum;

  // init (por si do_forecast=0 o L_fore=0)
  log_lik_real_sum = 0.0;

  if (do_forecast == 1 && L_fore > 0) {

    vector[N] b;
    vector[N] vec_var;

    vec_var = (1 - square(rho)) ./ (1 + (N_nei - rep_vector(1.0, N)) * square(rho));
    b       =  rho              ./ (1 + (N_nei - rep_vector(1.0, N)) * square(rho));

    vector[N] w_prev = to_vector(w[N_time, ]);
    real sigma = sqrt(tau2);

    for (ell in 1:L_fore) {
      vector[N] eta;
      vector[N] w_new;

      // 1) simulate innovations eta (DAGAR conditional)
      eta[1] = normal_rng(0, std_dev_w * sqrt(1 - square(phi)));

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

      // 2) propagate AR(1)
      w_new = phi * w_prev + eta;
      w_fore[ell, ] = to_row_vector(w_new);

      // 3) predictive mean for this horizon
      {
        vector[N] beta_st_fore;
        vector[N] mu_real;

        for (i in 1:N) {
          beta_st_fore[i] = theta[1]
                            + theta[2] * s_ind[i]
                            + theta[3] * time_fore[ell]
                            + theta[4] * s_ind[i] * time_fore[ell];
        }

        mu_real = beta_st_fore .* x_st_fore[ell] + w_new;

        // simulate posterior predictive y_fore (optional)
        for (i in 1:N) {
          y_fore[ell, i] = normal_rng(mu_real[i], sigma);

          // log predictive density of the REAL test observation
          log_lik_real[ell, i] =
            normal_lpdf(y_real[ell, i] | mu_real[i], sigma);

          log_lik_real_sum += log_lik_real[ell, i];
        }
      }

      w_prev = w_new;
    }

  } else {
    // make outputs defined even when not forecasting
    w_fore        = rep_matrix(0.0, L_fore, N);
    eta_fore      = rep_matrix(0.0, L_fore, N);
    y_fore        = rep_matrix(0.0, L_fore, N);
    log_lik_real  = rep_matrix(0.0, L_fore, N);
  }
}
