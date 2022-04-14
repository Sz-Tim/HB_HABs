# HABReports Bayesian modelling
# Hurdle Frechet functions
# Tim Szewczyk





# Hurdle Frechet family function
hurdle_frechet <- function(name="hurdle_frechet", 
                           dpars=c("mu", "nu", "hu"),
                           links=c("log", "logm1", "probit"), 
                           type="real", 
                           lb=c(NA, 1, 0), 
                           ub=c(NA, NA, 1)) {
  brms::custom_family(
    name, dpars=dpars,
    links=links, type=type, 
    lb=lb, ub=ub
  )
}



# Hurdle Frechet log likelihood
log_lik_hurdle_frechet <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  nu <- brms::get_dpar(prep, "nu", i = i)
  hu <- brms::get_dpar(prep, "hu", i = i)
  y <- prep$data$Y[i]
  hurdle_frechet_lpdf(y, mu, nu, hu)
}



# Hurdle Frechet posterior predictions
posterior_predict_hurdle_frechet <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  nu <- brms::get_dpar(prep, "nu", i = i)
  hu <- brms::get_dpar(prep, "hu", i = i)
  hurdle_frechet_rng(mu, nu, hu)
}



# Hurdle Frechet posterior expected posterior predictions
posterior_epred_hurdle_frechet <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  hu <- brms::get_dpar(prep, "hu")
  mu * (1-hu)
}
