  /* hurdle frechet log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   mu: mean parameter of the frechet distribution 
   *   nu: scale parameter of the frechet distribution
   *   hu: hurdle probability, pr(y=0)
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
  real hurdle_frechet_lpdf(real y, real mu, real nu, real hu) {
    if (y == 0) { 
      return bernoulli_lpmf(1 | hu); 
    } else { 
      return bernoulli_lpmf(0 | hu) +  
             frechet_lpdf(y | mu, nu); 
    }
  }

  /* hurdle frechet log-PDF of a single response
   * logit parameterization of the hurdle part
   * Args: 
   *   y: the response value 
   *   mu: mean parameter of the frechet distribution 
   *   nu: scale parameter of the frechet distribution
   *   hu: linear predictor for the hurdle part, pr(y=0)
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
  real hurdle_frechet_logit_lpdf(real y, real mu, real nu, real hu) { 
    if (y == 0) { 
      return bernoulli_logit_lpmf(1 | hu); 
    } else { 
      return bernoulli_logit_lpmf(0 | hu) +  
             frechet_lpdf(y | mu, nu); 
    } 
  } 

  // hurdle frechet log-CCDF and log-CDF functions 
  real hurdle_frechet_lccdf(real y, real mu, real nu, real hu) { 
    return bernoulli_lpmf(0 | hu) + frechet_lccdf(y | mu, nu); 
  }
  real hurdle_frechet_lcdf(real y, real mu, real nu, real hu) { 
    return log1m_exp(hurdle_frechet_lccdf(y | mu, nu, hu));
  }

  // hurdle frechet random number generator
  real hurdle_frechet_rng(real mu, real nu, real hu) {
    return (1-bernoulli_rng(hu)) * frechet_rng(mu, nu);
  }
  