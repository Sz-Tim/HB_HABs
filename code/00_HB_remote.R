# HABReports Bayesian modelling
# Bayesian HAB model exploration
# Tim Szewczyk


# This script is a translation of parts of docs/02_HB_exploration.Rmd for running remotely



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "LaplacesDemon", "brms", "bayesplot")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
walk(dir("code", "00[0-9]_fn", full.names=T), source)

# minch2:    2013-06-20 to 2019-07-02
# WeStCOMS2: 2019-04-01 to 2022-01-26
if(.Platform$OS.type=="unix") {
  sep <- "/"
  westcoms.dir <- c("/media/archiver/common/sa01da-work/minch2/Archive/",
                    "/media/archiver/common/sa01da-work/WeStCOMS2/Archive/")
  mesh.f <- c("/home/sa04ts/FVCOM_meshes/WeStCOMS_mesh.gpkg",
              "/home/sa04ts/FVCOM_meshes/WeStCOMS2_mesh.gpkg")
} else {
  sep <- "\\"
  westcoms.dir <- c("D:\\hydroOut\\minch2\\Archive\\",
                    "D:\\hydroOut\\WestCOMS2\\Archive\\")
  mesh.f <- c("..\\..\\01_FVCOM\\data\\WeStCOMS_Mesh.gpkg",
              "..\\..\\01_FVCOM\\data\\WeStCOMS2_Mesh.gpkg")
}





# sampling_df -------------------------------------------------------------

sampling.df <- read_csv(glue("data{sep}obs_df.csv")) %>%
  # mutate(target_sp=pseudo_nitzschia_sp,
  #        target_cat=case_when(target_sp==0 ~ 1,
  #                             between(target_sp, 0, 40e3) ~ 2,
  #                             between(target_sp, 40e3, 50e3) ~ 3,
  #                             between(target_sp, 50e3, 100e3) ~ 4,
  #                             between(target_sp, 100e3, 200e3) ~ 5,
  #                             between(target_sp, 200e3, 500e3) ~ 6,
  #                             between(target_sp, 500e3, 1e6) ~ 7,
  #                             target_sp > 1e6 ~ 8)) %>%
  mutate(target_sp=karenia_mikimotoi,
         target_catFull=case_when(target_sp==0 ~ 0,
                                  between(target_sp, 0, 2.5e3) ~ 1,
                                  between(target_sp, 2.5e3, 5e3) ~ 2,
                                  between(target_sp, 5e3, 20e3) ~ 3,
                                  between(target_sp, 20e3, 100e3) ~ 4,
                                  between(target_sp, 100e3, 500e3) ~ 5,
                                  between(target_sp, 500e3, 1e6) ~ 6,
                                  target_sp > 1e6 ~ 7)) %>%
  mutate(target_sp=karenia_mikimotoi,
         target_cat=case_when(target_sp==0 ~ 0,
                              between(target_sp, 0, 2.5e3) ~ 1,
                              between(target_sp, 2.5e3, 5e3) ~ 2,
                              between(target_sp, 5e3, Inf) ~ 3)) %>%
  mutate(yday=yday(date_collected),
         ydayCos=cos(2*pi*yday/365),
         ydaySin=sin(2*pi*yday/365),
         wind_lag_0=sqrt(uwind_speed_lag_0^2 + vwind_speed_lag_0^2),
         wind_lag_1=sqrt(uwind_speed_lag_1^2 + vwind_speed_lag_1^2),
         wind_lag_2=sqrt(uwind_speed_lag_2^2 + vwind_speed_lag_2^2),
         wind_lag_3=sqrt(uwind_speed_lag_3^2 + vwind_speed_lag_3^2),
         wind_lag_4=sqrt(uwind_speed_lag_4^2 + vwind_speed_lag_4^2),
         wind_lag_5=sqrt(uwind_speed_lag_5^2 + vwind_speed_lag_5^2),
         windDir_lag_0=atan2(vwind_speed_lag_0, uwind_speed_lag_0),
         windDir_lag_1=atan2(vwind_speed_lag_1, uwind_speed_lag_1),
         windDir_lag_2=atan2(vwind_speed_lag_2, uwind_speed_lag_2),
         windDir_lag_3=atan2(vwind_speed_lag_3, uwind_speed_lag_3),
         windDir_lag_4=atan2(vwind_speed_lag_4, uwind_speed_lag_4),
         windDir_lag_5=atan2(vwind_speed_lag_5, uwind_speed_lag_5),
         water_lag_0=sqrt(u_lag_0^2 + v_lag_0^2),
         water_lag_1=sqrt(u_lag_1^2 + v_lag_1^2),
         water_lag_2=sqrt(u_lag_2^2 + v_lag_2^2),
         water_lag_3=sqrt(u_lag_3^2 + v_lag_3^2),
         water_lag_4=sqrt(u_lag_4^2 + v_lag_4^2),
         water_lag_5=sqrt(u_lag_5^2 + v_lag_5^2),
         waterDir_lag_0=atan2(v_lag_0, u_lag_0),
         waterDir_lag_1=atan2(v_lag_1, u_lag_1),
         waterDir_lag_2=atan2(v_lag_2, u_lag_2),
         waterDir_lag_3=atan2(v_lag_3, u_lag_3),
         waterDir_lag_4=atan2(v_lag_4, u_lag_4),
         waterDir_lag_5=atan2(v_lag_5, u_lag_5),
         short_wave_wk=short_wave_lag_1+short_wave_lag_2+short_wave_lag_3+short_wave_lag_4+short_wave_lag_5,
         temp_wk=temp_lag_1+temp_lag_2+temp_lag_3+temp_lag_4+temp_lag_5,
         wind_wk=wind_lag_1+wind_lag_2+wind_lag_3+wind_lag_4+wind_lag_5,
         windDir_wk=windDir_lag_1+windDir_lag_2+windDir_lag_3+windDir_lag_4+windDir_lag_5,
         water_wk=water_lag_1+water_lag_2+water_lag_3+water_lag_4+water_lag_5,
         waterDir_wk=waterDir_lag_1+waterDir_lag_2+waterDir_lag_3+waterDir_lag_4+waterDir_lag_5) %>%
  mutate(across(starts_with(c("time", "depth", "tides", "short_wave", "zeta", "temp", "wind", "water")),
                CenterScale,
                .names="{.col}_sc")) %>%
  arrange(site.id, date_collected) %>%
  group_by(site.id) %>% 
  mutate(target_sp_lag1=lag(log(target_sp+1), 1),
         target_sp_lag2=lag(log(target_sp+1), 2),
         target_sp_lag3=lag(log(target_sp+1), 3),
         target_cat_lag1=lag(target_cat, 1),
         target_cat_lag2=lag(target_cat, 2),
         target_cat_lag3=lag(target_cat, 3),
         target_catFull_lag1=lag(target_catFull, 1),
         target_catFull_lag2=lag(target_catFull, 2),
         target_catFull_lag3=lag(target_catFull, 3),
         target_sp_binary=target_cat > 0,
         target_sp_binary_lag1=target_cat_lag1 > 0,
         date_lag1=as.numeric(date_collected-lag(date_collected, 1)),
         date_lag2=as.numeric(date_collected-lag(date_collected, 2)),
         date_lag3=as.numeric(date_collected-lag(date_collected, 3)),
         target_sp_lagWt1=(1+target_sp_lag1)/date_lag1,
         target_sp_lagWt2=(1+target_sp_lag2)/date_lag2,
         target_sp_lagWt3=(1+target_sp_lag3)/date_lag3) %>%
  mutate(nObs=n(),
         obs_id=row_number()) %>%
  ungroup %>%
  filter(complete.cases(.)) %>%
  mutate(wk=floor(as.numeric(date - min(date))/7),
         target_cat_f=factor(target_cat, levels=0:3, ordered=T),
         target_cat_f_lag1=factor(target_cat_lag1, levels=0:3, ordered=T),
         target_cat_f_lag2=factor(target_cat_lag2, levels=0:3, ordered=T),
         target_cat_shifted=target_cat+1,
         target_cat_shifted_f=factor(target_cat_shifted, levels=1:4, ordered=T),
         target_sp_log=log(target_sp+1),
         target_sp=round(target_sp)) %>%
  select(site.id, lon, lat, date, year, obs_id, nObs, wk,
         starts_with("target_"), starts_with("date_"), starts_with("yday"),
         ends_with("_wk_sc"))







# hurdle frechet ----------------------------------------------------------

week_i <- 153
train.df <- sampling.df %>% filter(year %in% c(2013, 2016, 2019, 2021))
  # filter(wk < week_i)
test.df <- sampling.df %>% filter(year %in% c(2013, 2016, 2019, 2021))
  # filter(wk == week_i)
train.df_no0 <- sampling.df %>%  filter(year %in% c(2013, 2016, 2019, 2021)) %>%
  # filter(wk < week_i) %>%
  filter(target_cat > 0)
test.df_no0 <- sampling.df %>% filter(year %in% c(2013, 2016, 2019, 2021)) %>%
  # filter(wk == week_i) %>%
  filter(target_cat > 0)


predictors_main <- "ydayCos + ydaySin +
                 
                 short_wave_wk_sc:ydayCos +
                 temp_wk_sc:ydayCos +
                 
                 target_sp_lagWt1:ydayCos + 
                 target_sp_lagWt2:ydayCos + 
                 target_cat_f_lag1:ydayCos +
                 target_cat_f_lag2:ydayCos +

                 target_cat_f_lag1:target_cat_f_lag2 +
                 
                 target_sp_lagWt1:temp_wk_sc +
                 target_sp_lagWt1:short_wave_wk_sc +
                 
                 (1|site.id)"
predictors_hu <- "ydayCos + ydaySin +
                  ydayCos:target_cat_f_lag1 +
                  ydaySin:target_cat_f_lag1"

form_hu_frechet <- bf(glue("target_sp_log ~ {predictors_main}"),
                      glue("hu ~ {predictors_hu}"))


# Promising: hurdle Frechet. Cut and paste from hurdle_lognormal and frechet.

hurdle_frechet <- custom_family(
  "hurdle_frechet", dpars=c("mu", "nu", "hu"),
  links=c("log", "logm1", "probit"), type="real", 
  lb=c(NA, 1, 0), ub=c(NA, NA, 1)
)


stanvars <- stanvar(scode=readr::read_file("models\\hurdle_frechet_fn.stan"), block="functions")

out3.hu_frechet <- brm(form_hu_frechet, data=train.df, 
                      family=hurdle_frechet, stanvars=stanvars,
                      chains=4, cores=4, init="0",
                      control=list(adapt_delta=0.9, max_treedepth=20),
                      prior=c(prior(normal(2, 1), "Intercept", dpar="hu"),
                              prior(horseshoe(1), "b")),
                      save_model=glue("temp{sep}00_mod_hu_frechet.stan"))
saveRDS(out3.hu_frechet, "temp/out_hu_frechet.rds")
expose_functions(out3.hu_frechet, vectorize=T)

log_lik_hurdle_frechet <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  nu <- brms::get_dpar(prep, "nu", i = i)
  hu <- brms::get_dpar(prep, "hu", i = i)
  y <- prep$data$Y[i]
  hurdle_frechet_lpdf(y, mu, nu, hu)
}

posterior_predict_hurdle_frechet <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  nu <- brms::get_dpar(prep, "nu", i = i)
  hu <- brms::get_dpar(prep, "hu", i = i)
  hurdle_frechet_rng(mu, nu, hu)
}

posterior_epred_hurdle_frechet <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  hu <- brms::get_dpar(prep, "hu")
  mu * (1-hu)
}

pp_check(out3.hu_frechet, ndraws=100) + 
  scale_x_continuous(trans="log1p", limits=c(0, max(train.df$target_sp_log)))
conditional_effects(out3.hu_frechet)

mcmc_areas(out3.hu_frechet, regex_pars="b_")
mcmc_trace(out.hu_frechet, regex_pars="b_")


out.pred <- posterior_epred(out3.hu_frechet, newdata=test.df, allow_new_levels=T)
out_summary <- summary(out.hu_frechet)

pred.df <- test.df %>%
  mutate(pred_mn=colMeans(out.pred), 
         pred_se=apply(out.pred, 2, function(x) sd(x)/sqrt(length(x))),
         pred_q10=apply(out.pred, 2, function(x) quantile(x, probs=0.1)),
         pred_q50=apply(out.pred, 2, function(x) quantile(x, probs=0.5)),
         pred_q90=apply(out.pred, 2, function(x) quantile(x, probs=0.9)),
         pred_P=colMeans(out.pred > 0.5),
         epred_mn=colMeans(exp(out.pred)-1), 
         epred_se=apply(exp(out.pred)-1, 2, function(x) sd(x)/sqrt(length(x))),
         epred_q10=apply(exp(out.pred)-1, 2, function(x) quantile(x, probs=0.1)),
         epred_q50=apply(exp(out.pred)-1, 2, function(x) quantile(x, probs=0.5)),
         epred_q90=apply(exp(out.pred)-1, 2, function(x) quantile(x, probs=0.9)),
         epred_P=colMeans((exp(out.pred)-1) > 1))
ggplot(pred.df, aes(epred_q50, target_sp_log)) + geom_point() + stat_smooth(method="lm") + geom_abline()
ggplot(pred.df, aes(epred_P, as.numeric(target_sp_binary))) + geom_point(alpha=0.5) +
  stat_smooth(method="glm", se=F, method.args=list(family=binomial), size=0.5)
ggplot(pred.df, aes(target_cat_f, pred_mn)) + geom_boxplot()
ggplot(pred.df, aes(target_cat_f, epred_mn)) + geom_boxplot()

ggplot(pred.df, aes(pred_mn, as.numeric(target_sp_binary),
                    xmin=pred_q10, xmax=pred_q90)) + 
  geom_point(aes(colour=target_cat_shifted_f)) + geom_linerange(aes(colour=target_cat_shifted_f)) + 
  stat_smooth(method="glm", se=F, method.args=list(family=binomial), size=0.5) +
  scale_colour_manual("Obs. cat.", values=c("1"="green3", "2"="gold1", "3"="orange", "4"="red")) +
  labs(y="Observed binary category", x="Predicted mean + 80% CI (prob)")


out.pred <- posterior_epred(out.hu_frechet)

pred.df <- train.df %>%
  mutate(pred_mn=colMeans(out.pred), 
         pred_se=apply(out.pred, 2, function(x) sd(x)/sqrt(length(x))),
         pred_q10=apply(out.pred, 2, function(x) quantile(x, probs=0.1)),
         pred_q50=apply(out.pred, 2, function(x) quantile(x, probs=0.5)),
         pred_q90=apply(out.pred, 2, function(x) quantile(x, probs=0.9)),
         pred_P=colMeans(out.pred > 0.5),
         epred_mn=colMeans(exp(out.pred)-1), 
         epred_se=apply(exp(out.pred)-1, 2, function(x) sd(x)/sqrt(length(x))),
         epred_q10=apply(exp(out.pred)-1, 2, function(x) quantile(x, probs=0.1)),
         epred_q50=apply(exp(out.pred)-1, 2, function(x) quantile(x, probs=0.5)),
         epred_q90=apply(exp(out.pred)-1, 2, function(x) quantile(x, probs=0.9)),
         epred_P=colMeans((exp(out.pred)-1) > 1))
ggplot(pred.df, aes(epred_q50, target_sp_log)) + geom_abline() +
  geom_point(aes(colour=target_cat_shifted_f)) + stat_smooth(method="lm") +
  scale_colour_manual("Obs. cat.", values=c("1"="green3", "2"="gold1", "3"="orange", "4"="red"))
ggplot(pred.df, aes(epred_P, as.numeric(target_sp_binary))) + geom_point(alpha=0.5) +
  stat_smooth(method="glm", se=F, method.args=list(family=binomial), size=0.5)
ggplot(pred.df, aes(target_cat_f, pred_q50)) + geom_boxplot()
ggplot(pred.df, aes(target_cat_f, epred_P)) + geom_boxplot()

ggplot(pred.df, aes(pred_q50, as.numeric(target_sp_binary),
                    xmin=pred_q10, xmax=pred_q90)) + 
  geom_point(aes(colour=target_cat_shifted_f)) + geom_linerange(aes(colour=target_cat_shifted_f)) + 
  stat_smooth(method="glm", se=F, method.args=list(family=binomial), size=0.5) +
  scale_colour_manual("Obs. cat.", values=c("1"="green3", "2"="gold1", "3"="orange", "4"="red")) +
  labs(y="Observed binary category", x="Predicted mean + 80% CI (prob)")

ggplot(pred.df, aes(date, epred_q50-target_sp_log, colour=target_cat_shifted_f)) + 
  geom_point() +
  scale_colour_manual("Obs. cat.", values=c("1"="green3", "2"="gold1", "3"="orange", "4"="red")) +
  facet_wrap(~target_cat_shifted_f)




# model exploration -------------------------------------------------------

week_i <- 153
train.df <- sampling.df %>% 
  filter(wk < week_i)
test.df <- sampling.df %>%
  filter(wk == week_i)
train.df_no0 <- sampling.df %>% 
  filter(wk < week_i) %>%
  filter(target_cat > 0)
test.df_no0 <- sampling.df %>%
  filter(wk == week_i) %>%
  filter(target_cat > 0)


predictors <- "ydayCos + ydaySin +
                 
                 short_wave_wk_sc:ydayCos +
                 temp_wk_sc:ydayCos +
                 windDir_wk_sc:ydayCos +
                 waterDir_wk_sc:ydayCos +
                 
                 target_sp_lag1:target_sp_lag2:target_sp_lag3:ydayCos +
                 
                 target_sp_lagWt1:ydayCos + 
                 target_sp_lagWt2:ydayCos + 
                 target_sp_lagWt3:ydayCos +
                 
                 target_sp_lagWt1:temp_wk_sc +
                 target_sp_lagWt1:short_wave_wk_sc +
                 target_sp_lagWt1:wind_wk_sc +
                 
                 (1|site.id)"

form_ord_1to3 <- bf(glue("target_cat | thres(2) ~ {predictors}"))
form_ord_1to4 <- bf(glue("target_cat | thres(3) ~ {predictors}"))
form_dens_log <- bf(glue("target_sp_log ~ {predictors}"))
form_count_zi_0to7 <- bf(glue("target_catFull | trunc(0,7) ~ {predictors}"),
                         zi ~ ydayCos*ydaySin*target_sp_lag1)
form_count_hu_0to7 <- bf(glue("target_catFull | trunc(0,7) ~ {predictors}"),
                         hu ~ ydayCos*ydaySin*target_sp_lag1)
form_bern <- bf(glue("target_sp_binary ~ {predictors}"))

# Ordinal models: GYOR No 0's
mods.ord <- c("logit", "probit", "probit_approx", "cauchit") %>%
  setNames(glue("ord_{.}"))
out.ord_no0 <- imap(mods.ord, 
                    ~brm(form_ord_1to3, family=cumulative(.x),
                         train.df_no0, cores=4, inits="0",
                         file=glue("temp{sep}00_mod_{.y}_no0"),
                         save_model=glue("temp{sep}00_mod_{.y}.stan")))
loo_compare(map(out.ord_no0, loo))



# Ordinal models: GYOR 0's
mods.ord <- c("logit", "probit", "probit_approx", "cauchit") %>%
  setNames(glue("ord_{.}"))
out.ord <- imap(mods.ord, 
                ~brm(form_ord_1to4, family=cumulative(.x),
                     train.df, cores=4, inits="0",
                     file=glue("temp{sep}00_mod_{.y}"),
                     save_model=glue("temp{sep}00_mod_{.y}.stan")))
loo_compare(map(out.ord, loo))



# Density models: no 0's
mods.dens <- c(weibull=weibull, frechet=frechet, sk_norm=skew_normal)
out.dens_no0 <- imap(mods.dens, 
                     ~brm(form_dens_log, family=.x,
                          train.df_no0, cores=4, inits="0",
                          file=glue("temp{sep}00_mod_{.y}_no0"),
                          save_model=glue("temp{sep}00_mod_{.y}.stan")))
loo_compare(map(out.dens_no0, loo))



# Density models: 0's
mods.dens <- c(h_gam=hurdle_gamma, h_lnorm=hurdle_lognormal)
out.dens <- imap(mods.dens, 
                 ~brm(form_dens_log, family=.x,
                      train.df, cores=4, inits="0",
                      file=glue("temp{sep}00_mod_{.y}"),
                      save_model=glue("temp{sep}00_mod_{.y}.stan")))
loo_compare(map(out.dens, loo))



# Binary models
mods.bern <- c(bern_logit="logit", bern_probit="probit")
out.bern <- imap(mods.bern, 
                 ~brm(form_bern, family=bernoulli(.x),
                      train.df, cores=4, inits="0",
                      file=glue("temp{sep}00_mod_{.y}"),
                      save_model=glue("temp{sep}00_mod_{.y}.stan")))
loo_compare(map(out.bern, loo))



# Count models: all cat 0's
mod.zi <- c(zi_pois=zero_inflated_poisson, zi_nb=zero_inflated_negbinomial)
out.zi <- imap(mod.zi, 
               ~brm(form_count_zi_0to7, family=.x,
                    train.df, cores=4, inits="0",
                    file=glue("temp{sep}00_mod_{.y}"),
                    save_model=glue("temp{sep}00_mod_{.y}.stan")))
mod.hu <- c(h_pois=hurdle_poisson, hu_nb=hurdle_negbinomial)
out.hu <- imap(mod.hu, 
               ~brm(form_count_hu_0to7, family=.x,
                    train.df, cores=4, inits="0",
                    file=glue("temp{sep}00_mod_{.y}"),
                    save_model=glue("temp{sep}00_mod_{.y}.stan")))
loo_compare(c(map(out.zi, loo), map(out.hu, loo)))





# model runs --------------------------------------------------------------

weeks <- sort(unique(sampling.df$wk))[-(1:110)]
out <- vector("list", length(weeks))

for(i in 1:length(weeks)) {
  week_i <- weeks[i]
  
  if(i==1) {
    train.df <- sampling.df %>% 
      filter(wk < week_i)
  } else {
    train.df <- sampling.df %>%
      filter(wk == week_i-1)
  }
  test.df <- sampling.df %>%
    filter(wk == week_i) 
  date_range <- paste(as.character(range(test.df$date)), collapse=" to ")
  cat("Starting", i, "= week", as.character(week_i), "=", date_range, "\n")
  cat("  train:", nrow(train.df), "rows with", n_distinct(train.df$site.id), "sites",
      "    ", table(train.df$target_cat_f), "\n")
  cat("  test:", nrow(test.df), "rows with", n_distinct(test.df$site.id), "sites",
      "    ", table(test.df$target_cat_f), "\n")
  
  
  if(i == 1) {
    prior_i <- set_prior("normal(0,2)", class="b")
  } else {
    post_im1 <- readRDS(glue("temp{sep}out_{weeks[i-1]}_ord.rds"))
    b_im1 <- as_draws_df(post_im1, variable="b_", regex=T) %>%
      pivot_longer(cols=starts_with("b_"), names_to="param", values_to="val") %>%
      group_by(param) %>%
      summarise(mn=mean(val), sd=sd(val)) %>%
      mutate(coef=str_remove(param, "b_"))
    eff_im1 <- b_im1 %>% filter(!grepl("Intercept", param))
    int_im1 <- b_im1 %>% filter(grepl("Intercept", param)) %>%
      mutate(coef=str_remove(str_remove(str_remove(coef, "Intercept"), "\\["), "]"))
    sd_im1 <- summary(post_im1)$random$site.id
    prior.ls <- vector("list", nrow(b_im1)+1)
    for(j in 1:nrow(int_im1)) {
      prior.ls[[j]] <- prior_string(glue("normal({int_im1$mn[j]},{int_im1$sd[j]})"),
                                    class="Intercept", coef=int_im1$coef[j])
    }
    for(j in 1:nrow(eff_im1)) {
      prior.ls[[j+nrow(int_im1)]] <- prior_string(glue("normal({eff_im1$mn[j]},{eff_im1$sd[j]})"),
                                                  class="b", coef=eff_im1$coef[j])
    }
    prior.ls[[nrow(b_im1)+1]] <- prior_string(glue("normal({sd_im1$Estimate},{sqrt(sd_im1$Est.Error)})"), 
                                              class="sd")
    prior_i <- do.call(rbind, prior.ls)
  }
  
  if(i==1) {
    out[[i]] <- brm(target_cat | thres(3) ~
                      
                      ydayCos + ydaySin +
                      
                      short_wave_wk_sc:ydayCos +
                      temp_wk_sc:ydayCos +
                      windDir_wk_sc:ydayCos +
                      waterDir_wk_sc:ydayCos +
                      
                      target_sp_lag1:target_sp_lag2:target_sp_lag3:ydayCos +
                      
                      target_sp_lagWt1:ydayCos + 
                      target_sp_lagWt2:ydayCos + 
                      target_sp_lagWt3:ydayCos +
                      
                      # target_sp_lag1:date_lag1 + 
                      # target_sp_lag2:date_lag2 +
                      # target_sp_lag3:date_lag3 +
                      # target_sp_lag1:date_lag1:ydayCos + 
                      # target_sp_lag2:date_lag2:ydayCos + 
                      # target_sp_lag3:date_lag3:ydayCos +
                      
                      target_sp_lagWt1:temp_wk_sc +
                      target_sp_lagWt1:short_wave_wk_sc +
                      target_sp_lagWt1:wind_wk_sc +
                      
                      (1|site.id),
                    data=train.df, cores=4,  family=cumulative(probit),
                    prior=prior_i, file=glue("temp{sep}out_{week_i}_ord"),
                    save_model=glue("temp{sep}mod_{week_i}_ord.stan")) 
  } else {
    out[[i]] <- update(out[[i-1]], newdata=train.df, cores=4, 
                       prior=prior_i, file=glue("temp{sep}out_{week_i}_ord"))
  }
  
  cat("  Fitted model \n")
  
  out.logitpred <- posterior_linpred(out[[i]], transform=F, newdata=test.df, allow_new_levels=T)
  out.pred <- posterior_linpred(out[[i]], transform=T, newdata=test.df, allow_new_levels=T)
  out_summary <- summary(out[[i]])
  # cutoffs <- c(-Inf, 
  #              out_summary$fixed[grep("Intercept", rownames(out_summary$fixed)),1],
  #              Inf)
  
  cat("  Made predictions \n")
  
  pred.df <- test.df %>%
    mutate(pred_mn=colMeans(out.pred), 
           pred_se=apply(out.pred, 2, function(x) sd(x)/sqrt(length(x))),
           pred_q10=apply(out.pred, 2, function(x) quantile(x, probs=0.1)),
           pred_q90=apply(out.pred, 2, function(x) quantile(x, probs=0.9))) %>%
    mutate(lpred_mn=colMeans(out.logitpred), 
           lpred_se=apply(out.logitpred, 2, function(x) sd(x)/sqrt(length(x))),
           lpred_q10=apply(out.logitpred, 2, function(x) quantile(x, probs=0.1)),
           lpred_q90=apply(out.logitpred, 2, function(x) quantile(x, probs=0.9))) %>%
    # bind_cols(imap_dfc(setNames(2:length(cutoffs), glue("pr{1:(length(cutoffs)-1)}")), 
    #                    ~apply(out.pred, 2, 
    #                           function(x) sum(x>cutoffs[.x-1] & x<cutoffs[.x])/nrow(out.pred)))) %>%
    mutate(target_cat_f=factor(target_cat, levels=1:4))
  
  p <- ggplot(pred.df, aes(log(target_sp+1), pred_mn,
                           ymin=pred_q10, ymax=pred_q90)) + 
    geom_point(aes(colour=target_cat_f)) + geom_errorbar(aes(colour=target_cat_f), width=0.1) + 
    stat_smooth(method="glm", se=F, method.args=list(family=binomial), size=0.5) +
    scale_colour_manual("Obs. cat.", values=c("1"="green3", "2"="gold1", "3"="orange", "4"="red")) +
    ylim(0, 1) +
    labs(title=glue("Week {week_i}: {date_range}"), 
         x="observed log(N+1)", y="Predicted mean + 80% CI (prob)")
  ggsave(glue("temp{sep}plot_lpred_{week_i}.png"), p, height=4, width=5, dpi=200)
  
  p <- ggplot(pred.df, aes(pred_mn, as.numeric(target_sp_binary),
                      xmin=pred_q10, xmax=pred_q90)) + 
    geom_point(aes(colour=target_cat_f)) + geom_linerange(aes(colour=target_cat_f)) + 
    stat_smooth(method="glm", se=F, method.args=list(family=binomial), size=0.5) +
    scale_colour_manual("Obs. cat.", values=c("1"="green3", "2"="gold1", "3"="orange", "4"="red")) +
    ylim(0, 1) + xlim(0, 1) +
    labs(title=glue("Week {week_i}: {date_range}"),
         y="Observed binary category", x="Predicted mean + 80% CI (prob)")
  ggsave(glue("temp{sep}plot_lpred2_{week_i}.png"), p, height=4, width=5, dpi=200)
  
  p <- ggplot(pred.df, aes(log(target_sp+1), lpred_mn,
                           ymin=lpred_q10, ymax=lpred_q90)) + 
    geom_point(aes(colour=target_cat_f)) + geom_errorbar(aes(colour=target_cat_f), width=0.1) + 
    stat_smooth(se=F, method="lm", formula=y~x, size=0.5) +
    scale_colour_manual("Obs. cat.", values=c("1"="green3", "2"="gold1", "3"="orange", "4"="red")) +
    labs(title=glue("Week {week_i}: {date_range}"), 
         x="observed log(N+1)", y="Predicted mean + 80% CI (logit)")
  ggsave(glue("temp{sep}plot_pred_{week_i}.png"), p, height=4, width=5, dpi=200)
  
  p <- ggplot(pred.df, aes(lpred_mn, log(target_sp+1), 
                           xmin=lpred_q10, xmax=lpred_q90)) + 
    geom_point(aes(colour=target_cat_f)) + geom_linerange(aes(colour=target_cat_f)) + 
    stat_smooth(se=F, method="lm", formula=y~x, size=0.5) +
    scale_colour_manual("Obs. cat.", values=c("1"="green3", "2"="gold1", "3"="orange", "4"="red")) +
    labs(title=glue("Week {week_i}: {date_range}"), 
         y="observed log(N+1)", x="Predicted mean + 80% CI (logit)")
  ggsave(glue("temp{sep}plot_pred2_{week_i}.png"), p, height=4, width=5, dpi=200)
  
  cat("  Saving output \n")
  saveRDS(out_summary$fixed, glue("temp{sep}fixed_eff_{week_i}.rds"))
  saveRDS(pred.df, glue("temp{sep}pred_{week_i}_ord.rds"))
}






# post-hoc ----------------------------------------------------------------

coef.df <- dir("temp", "fixed.*rds") %>%
  map_dfr(~readRDS(glue("temp{sep}{.x}")) %>%
            rownames_to_column("param") %>%
            mutate(week=as.numeric(str_sub(.x, 11, -5))))

coef.df <- read_csv(glue("temp{sep}fixed_eff_aggregated.csv"))
ggplot(coef.df, aes(week, y=Estimate, ymin=`l-95% CI`, ymax=`u-95% CI`)) +
  geom_hline(yintercept=0) +
  geom_ribbon(alpha=0.25) +
  geom_line() + facet_wrap(~param, scales="free_y")


pred.df <- dir("temp", "pred.*rds") %>%
  map_dfr(~readRDS(glue("temp{sep}{.x}")) %>%
            rownames_to_column("param") %>%
            mutate(week=as.numeric(str_sub(.x, 11, -5))))
pred.df <- read_csv(glue("temp{sep}pred_df_aggregated.csv"))
ggplot(pred.df, aes(lpred_mn, y=target_cat)) + geom_point()
ggplot(pred.df, aes(lpred_mn, y=log(target_sp+1))) + geom_point()
ggplot(pred.df, aes(factor(target_cat), lpred_mn)) + geom_boxplot() 
ggplot(pred.df, aes(date, lpred_mn-log(target_sp+1), colour=target_sp_lag1)) + 
  geom_point() + facet_wrap(~target_cat) + scale_colour_viridis_c()
ggplot(pred.df, aes(date, lpred_mn, colour=target_cat_f)) + geom_point() +
  scale_colour_viridis_c()


# Events for K. mikimotoi only occurred in 2013, 2016, 2019, and 2021
# Updating priors weekly ends up with intercept wandering because of the lack
# of 3's and 4's occurring. Maybe set the prior based on only the above years?
sampling.df %>% group_by(year, target_cat) %>% summarise(N=n()) %>% arrange(desc(target_cat), year)

pred.df %>% filter(year %in% c(2013, 2016, 2019, 2021)) %>%
  ggplot(aes(lpred_mn, y=target_cat)) + geom_point()
pred.df %>% filter(year %in% c(2013, 2016, 2019, 2021)) %>%
  ggplot(aes(lpred_mn, y=log(target_sp+1))) + geom_point()
pred.df %>% filter(year %in% c(2013, 2016, 2019, 2021)) %>%
  ggplot(aes(factor(target_cat), lpred_mn)) + geom_boxplot() 
pred.df %>% filter(year %in% c(2013, 2016, 2019, 2021)) %>%
  ggplot(aes(date, lpred_mn-log(target_sp+1))) + geom_point() + facet_wrap(~target_cat)






# gbm ---------------------------------------------------------------------

# library(gbm)
# week_i <- 158
# train.df <- sampling.df %>%
#   filter(wk < week_i)
# test.df <- sampling.df %>%
#   filter(wk == week_i) %>%
#   filter(site.id %in% train.df$site.id)
# 
# train.gbm <- train.df %>% ungroup %>% 
#   select(-obs_id, -nObs, -target_sp, -target_cat, -date_collected, -yday, -wk, -date)
# mod.gbm <- gbm(target_sp_binary~., 
#                "bernoulli", 
#                data=train.gbm, 
#                n.trees=1e3, 
#                cv.folds=10,
#                interaction.depth=5)
# print(mod.gbm)
# summary(mod.gbm) %>% filter(rel.inf > 0)
# 
# test.gbm <- test.df %>% ungroup %>%
#   select(-obs_id, -nObs, -target_sp, -target_cat, -date_collected, -yday, -wk, -date)
# test.x <- test.gbm %>% select(-target_sp_binary)
# test.y <- test.gbm %>% select(target_sp_binary)
# 
# pred.y = predict.gbm(mod.gbm, test.x)
# x.ax = 1:length(pred.y)
# plot(x.ax, as.numeric(test.y$target_sp_binary), col="blue", pch=20, cex=.9, ylim=c(0,1))
# points(x.ax, as.numeric(pred.y>0), col="red", pch=1, cex=1.25) 
# 
# 
# 
# train.gbm <- train.df %>% ungroup %>% mutate(target_sp=log(target_sp+1)) %>%
#   select(-obs_id, -nObs, -target_sp_binary, -target_cat, -date_collected, -yday, -wk, -date)
# mod.gbm <- gbm(target_sp~., 
#                "tdist", 
#                data=train.gbm, 
#                n.trees=1e3, 
#                cv.folds=10,
#                interaction.depth=5)
# print(mod.gbm)
# summary(mod.gbm) %>% filter(rel.inf > 0)
# 
# test.gbm <- test.df %>% ungroup %>%
#   select(-obs_id, -nObs, -target_sp_binary, -target_cat, -date_collected, -yday, -wk, -date)
# test.x <- test.gbm %>% select(-target_sp)
# test.y <- test.gbm %>% select(target_sp) %>% mutate(target_sp=log(target_sp+1))
# 
# pred.y = predict.gbm(mod.gbm, test.x)
# x.ax = 1:length(pred.y)
# plot(x.ax, test.y$target_sp, col="blue", pch=20, cex=.9)
# points(x.ax, pred.y, col="red", pch=1, cex=1.25) 
# plot(pred.y, test.y$target_sp); abline(a=0, b=1)
# plot(test.y$target_sp, pred.y-test.y$target_sp); abline(h=0)





# nnet --------------------------------------------------------------------

# library(nnet)
# 
# train.nnet <- train.df %>% ungroup %>% mutate(target_cat=factor(target_cat)) %>%
#   select(-obs_id, -nObs, -target_sp, -target_cat, -date_collected, -yday, -wk, -date)
# test.nnet <- test.df %>% ungroup %>% mutate(target_cat=factor(target_cat)) %>%
#   select(-obs_id, -nObs, -target_sp, -target_cat, -date_collected, -yday, -wk, -date)
# 
# mod.nnet <- nnet(target_sp_binary~., train.nnet, size=3, entropy=T, weights=train.nnet$target_sp_binary+1)
# pred.nnet <- predict(mod.nnet, newdata=test.nnet)
# plot(test.nnet$target_cat, pred.nnet[,1])








# scratch -----------------------------------------------------------------

sampling.df %>%
  ggplot(aes(target_cat_f, target_sp_lag1/(date_lag1))) + geom_boxplot()

sampling.df %>%
  ggplot(aes(target_sp_lag1/(date_lag1), log(target_sp+1), colour=target_cat_f)) + 
  geom_point(alpha=0.5)

sampling.df %>%
  ggplot(aes(target_sp_lagWt1, log(target_sp+1), colour=target_sp_lagWt2)) +
  geom_point() + scale_colour_viridis_c()

sampling.df %>% filter(target_cat>1) %>%
  ggplot(aes(target_cat_f, fill=target_sp_lagWt1==0)) + geom_bar(position="dodge")

sampling.df %>% filter(target_cat==4) %>% select(starts_with("target_sp_lag")) %>%
  pivot_longer(1:3, names_to="lag", values_to="val") %>%
  mutate(lagNum=as.numeric(str_sub(lag, -1, -1))) %>%
  ggplot(aes(val, fill=lag)) + geom_histogram() + facet_grid(.~lag) +
  ggtitle("Category 4 occurrences")


sampling.df %>% filter(target_cat==4) %>% select(starts_with("target_cat_lag"), year) %>%
  pivot_longer(1:3, names_to="lag", values_to="cat") %>%
  mutate(lagNum=as.numeric(str_sub(lag, -1, -1))) %>%
  ggplot(aes(factor(cat), fill=lag)) + geom_bar(position="dodge") + facet_grid(lag~year) +
  ggtitle("Category 4 occurrences")


sampling.df %>% filter(target_cat==3) %>% select(starts_with("target_sp_lag")) %>%
  pivot_longer(1:3, names_to="lag", values_to="val") %>%
  mutate(lagNum=as.numeric(str_sub(lag, -1, -1))) %>%
  ggplot(aes(val, fill=lag)) + geom_histogram() + facet_grid(.~lag) +
  ggtitle("Category 3 occurrences")


sampling.df %>% filter(target_cat==3) %>% select(starts_with("target_cat_lag"), year) %>%
  pivot_longer(1:3, names_to="lag", values_to="cat") %>%
  mutate(lagNum=as.numeric(str_sub(lag, -1, -1))) %>%
  ggplot(aes(factor(cat), fill=lag)) + geom_bar(position="dodge") + facet_grid(lag~year) +
  ggtitle("Category 3 occurrences")
