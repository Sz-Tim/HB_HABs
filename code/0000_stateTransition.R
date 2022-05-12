## -----------------------------------------------------------------------------
## project
## scriptDescription
## Tim Szewczyk
## -----------------------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(bayesplot)


# Simulate data
N <- 2e2
x <- cbind(1, sort(runif(N, -2, 2)))
b_01 <- c(boot::logit(0.1), -3)
b_11 <- c(boot::logit(0.3), -1.5)
p_01 <- boot::inv.logit(x %*% b_01)
p_11 <- boot::inv.logit(x %*% b_11)
y <- rep(0, N+1)
for(i in 1:N) {
  if(y[i] == 0) {
    y[i+1] <- rbinom(1, 1, p_01[i])
  } else {
    y[i+1] <- rbinom(1, 1, p_11[i])
  }
}


plot(y, type="l", ylim=c(0,1))
plot(x[,2], y[2:(N+1)])
plot(x[,2], y[2:(N+1)]-y[1:N])
plot(x[,2], p_01, ylim=c(0,1))
plot(x[,2], p_11, ylim=c(0,1))

data.ls <- list(N=N,
                y=y[-1],
                y_lag=y[-(N+1)],
                x=x)


mod <- cmdstan_model("models/state_transition.stan")

fit <- mod$sample(data=data.ls, refresh=1000, parallel_chains=4)

tibble(true=c(b_01, b_11)) %>%
  bind_cols(fit$summary(c(paste0("b_01[", 1:2, "]"), 
                          paste0("b_11[", 1:2, "]"))))

fit$summary("y_pred") %>%
  mutate(true=data.ls$y,
         y_lag=data.ls$y_lag) %>%
  ggplot(aes(mean, true)) + geom_point(alpha=0.2) + xlim(0, 1) +
  stat_smooth(method="glm", method.args=list(family="binomial"), fullrange=T)

fit$summary("y_pred") %>%
  mutate(true=data.ls$y,
         y_lag=data.ls$y_lag) %>%
  filter(y_lag==0 & true==1) %>%
  ggplot(aes(mean)) + geom_density()

fit$summary("p_01") %>%
  mutate(true=p_01) %>%
  ggplot(aes(mean, true)) + geom_point(alpha=0.2) +
  xlim(0,1) + ylim(0,1)

fit$summary("p_11") %>%
  mutate(true=p_11) %>%
  ggplot(aes(mean, true)) + geom_point(alpha=0.2) +
  xlim(0,1) + ylim(0,1)




sim.df <- tibble(x=x[,2], 
                 y=data.ls$y) %>%
  mutate(y_lag=data.ls$y_lag)
noBloom.df <- sim.df %>% filter(y_lag==0)
bloom.df <- sim.df %>% filter(y_lag==1)
noBloom.brm <- brms::brm(y ~ x, data=noBloom.df, family=brms::bernoulli())
bloom.brm <- brms::brm(y ~ x, data=bloom.df, family=brms::bernoulli())


