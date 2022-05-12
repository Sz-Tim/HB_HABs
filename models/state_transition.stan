//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> y;
  vector<lower=0, upper=1>[N] y_lag;
  matrix[N,2] x;
}

transformed data {
  vector[N] ones = x[,1];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[2] b_01;
  vector[2] b_11;
}

transformed parameters {
  vector<lower=0, upper=1>[N] p_01;
  vector<lower=0, upper=1>[N] p_11;
  
  p_01 = inv_logit(x * b_01);
  p_11 = inv_logit(x * b_11);
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  // b_start ~ normal(0, 1);
  // b_end ~ normal(0, 1);
  
  y ~ bernoulli((ones-y_lag).*p_01 + y_lag .* p_11);
  
  // for(i in 1:N) {
  //   if(y_lag[i] == 0) {
  //     y[i] ~ bernoulli(p_01[i]);
  //   } else {
  //     y[i] ~ bernoulli(p_11[i]);
  //   }
  // }
}

generated quantities {
  array[N] int<lower=0, upper=1> y_pred;
  for(i in 1:N) {
    if(y_lag[i] == 0) {
      y_pred[i] = bernoulli_rng(p_01[i]);
    } else {
      y_pred[i] = bernoulli_rng(p_11[i]);
    }
  }
}

