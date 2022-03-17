

data {
  int<lower=0> N;
  int<lower=0> nCov;
  matrix[N,nCov] X;
  int<lower=0> y[N];
}

parameters {
  vector[nCov] beta;
  vector[N] ln_lambda;
}

model {
  beta ~ normal(0, 1);
  ln_lambda ~ normal(3, 2);
  
  y ~ poisson_log_glm(X, ln_lambda, beta);
}

generated quantities {
  int<lower=0> y_pred[N];
  vector[N] log_lik;
  
  for(i in 1:N) {
    y_pred[i] = poisson_log_rng(ln_lambda[i] + X[i,] * beta);
    log_lik[i] = poisson_log_lpmf(y[i] | ln_lambda[i] + X[i,] * beta);
  }
}

