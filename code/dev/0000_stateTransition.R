## -----------------------------------------------------------------------------
## project
## scriptDescription
## Tim Szewczyk
## -----------------------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(bayesplot); theme_set(theme_classic())


# Simulate data
N <- 5e2
x <- cbind(1, 
           cos(1:N * 2 * pi / N),
           sin(1:N * 2 * pi / N),
           rnorm(N))
colnames(x) <- c("Intercept", paste0("X", 2:ncol(x)))
b_01 <- c(boot::logit(0.05), rnorm(ncol(x)-1))
b_11 <- c(boot::logit(0.3), rnorm(ncol(x)-1))
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
plot(x[,2], p_01, ylim=c(0,1))
plot(x[,3], p_01, ylim=c(0,1))
plot(x[,4], p_01, ylim=c(0,1))
plot(x[,2], p_11, ylim=c(0,1))
plot(x[,3], p_11, ylim=c(0,1))
plot(x[,4], p_11, ylim=c(0,1))

data.ls <- list(N=N,
                K=ncol(x)-1,
                y=y[-1],
                y_lag=y[-(N+1)],
                x=x)


mod <- cmdstan_model("models/state_transition.stan")

fit <- mod$sample(data=data.ls, refresh=1000, parallel_chains=4)

tibble(true=c(b_01, b_11)) %>%
  bind_cols(fit$summary(c(paste0("b_01[", seq_along(b_01), "]"), 
                          paste0("b_11[", seq_along(b_11), "]"))))

fit$summary("y_pred") %>%
  mutate(true=data.ls$y,
         y_lag=data.ls$y_lag) %>%
  ggplot(aes(mean, true)) + geom_point(alpha=0.2) + xlim(0, 1) +
  stat_smooth(method="glm", method.args=list(family="binomial"), fullrange=T)

fit$summary("y_pred") %>%
  mutate(true=data.ls$y,
         y_lag=data.ls$y_lag,
         bloom_t=true==1,
         bloom_tm1=y_lag==1) %>%
  ggplot(aes(mean, colour=bloom_t)) + geom_density() + facet_wrap(~bloom_tm1)

fit$summary("p_01") %>%
  mutate(true=p_01) %>%
  ggplot(aes(mean, true)) + geom_point(alpha=0.2) + geom_abline() + 
  xlim(0,1) + ylim(0,1)

fit$summary("p_11") %>%
  mutate(true=p_11) %>%
  ggplot(aes(mean, true)) + geom_point(alpha=0.2) + geom_abline() + 
  xlim(0,1) + ylim(0,1)




sim.df <- as_tibble(x[,-1]) %>%
  mutate(y=data.ls$y,
         y_lag=data.ls$y_lag)
noBloom.df <- sim.df %>% filter(y_lag==0)
bloom.df <- sim.df %>% filter(y_lag==1)
brm.form <- paste("y ~", paste(colnames(x)[-1], collapse=" + "))
noBloom.brm <- brms::brm(brm.form , 
                         data=noBloom.df, family=brms::bernoulli())
bloom.brm <- brms::brm(brm.form, 
                       data=bloom.df, family=brms::bernoulli())


