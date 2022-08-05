# HABReports Bayesian modelling
# Bayesian HAB model exploration
# Tim Szewczyk


# This script is a translation of parts of docs/02_HB_exploration.Rmd for running remotely



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "ncdf4", "sf", "LaplacesDemon",
          "cmdstanr", "bayesplot", "posterior", "jsonlite")
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


fsa.df <- fromJSON(glue("data{sep}copy_fsa.txt")) %>% 
  as_tibble %>% 
  select(-geom) %>%
  filter(easting < 4e10) %>% # entry error: Camb - Mid Yell Voe - 2013-04-16 
  filter(easting > 0 & northing > 0) %>%
  filter(!is.na(date_collected)) %>%
  filter(karenia_mikimotoi >= 0) %>% # -99 in karenia...?
  mutate(datetime_collected=as_datetime(date_collected),
         date_collected=date(datetime_collected),
         year=year(date_collected),
         month=month(date_collected),
         hour=pmin(20, pmax(5, hour(datetime_collected))),
         minutes=minute(datetime_collected),
         grid=if_else(date_collected < "2019-05-01", "minch2", "WeStCOMS2")) %>%
  filter(datetime_collected >= "2013-07-20") %>% 
  group_by(sin, area, site) %>% 
  mutate(lon=mean(easting), lat=mean(northing)) %>%
  ungroup %>% 
  mutate(site.id=as.numeric(factor(paste(sin, area, site)))) 







# munge datasets ----------------------------------------------------------

par.ls <- setSimParams()
# generate sampling details and lnlambdas
sampling.df <- fsa.df %>% 
  mutate(date=date_collected,
         dateChar=str_remove_all(date, "-"),
         time=hour + minutes/60,
         mode=sample(par.ls$modes, n(), T, c(0.7, 0.2, 0.1)),
         depth=sample(par.ls$depths, n(), T),
         tides=sample(par.ls$tides, n(), T, c(0.2, 0.05, 0.5, 0.05, 0.3)),
         lnlambda=rnorm(n(), 7, 2),
         lnlambdap1=log(exp(lnlambda)+1))

# extract element IDs, trinodes for each observation
mesh.sf <- map(mesh.f, st_read)
mesh <- map(mesh.f, loadMesh)
site.xy <- fsa.df %>% select("grid", "site.id", "lon", "lat") %>%
  group_by(grid, site.id) %>% slice_head(n=1) %>%
  group_by(grid) %>%
  group_split() %>%
  map2(.x=., .y=mesh.sf, 
       ~st_as_sf(.x, coords=c("lon", "lat"), remove=F, crs=27700) %>% 
         st_join(., .y %>% 
                   select(-area) %>% 
                   rename(depth.elem=depth,
                          site.elem=i)) %>% 
         filter(!is.na(site.elem)) %>%
         st_drop_geometry) %>%
  bind_rows
sampling.df <- sampling.df %>%
  right_join(., site.xy, by=c("grid", "site.id", "lon", "lat")) %>%
  group_by(grid) %>% group_split
site_trinode <- map2(.x=mesh, .y=sampling.df,
                     ~ncvar_get(.x, "trinodes")[.y$site.elem,])
walk(mesh, nc_close)

# load relevant data
out.df <- vector("list", 2)
for(grid in 1:2) {
  sampling_dates <- unique(sampling.df[[grid]]$dateChar)
  hydroVars.ls <- vector("list", length(sampling_dates))
  pb <- txtProgressBar(max=length(sampling_dates))
  for(i in 1:length(sampling_dates)) {
    date_i <- sampling_dates[i]
    rows_i <-  which(sampling.df[[grid]]$dateChar == date_i)
    hydroVars.ls[[i]] <- loadHydroVars(date_i, 
                                       site_trinode[[grid]][rows_i,,drop=F], 
                                       sampling.df[[grid]]$hour[rows_i],
                                       sampling.df[[grid]]$depth[rows_i],
                                       westcoms.dir[[grid]], 
                                       sep, 
                                       vars=c("temp", "short_wave", "zeta", 
                                              "uwind_speed", "vwind_speed", 
                                              "u", "v"), 
                                       lag=5) %>%
      mutate(rows=rows_i)
    setTxtProgressBar(pb, i)
  }
  out.df[[grid]] <- sampling.df[[grid]] %>%
    bind_cols(do.call(rbind, hydroVars.ls) %>% arrange(rows))
}


write_csv(do.call('rbind', out.df), glue("data{sep}obs_df.csv"))




# prepare model -----------------------------------------------------------

# test_obs identifies which observation to use as a test, with everything prior
# used to train. test_id = max(obs_id) - test_obs

# Need to know which thresholds to use. Different for each species on HABReports,
# and it isn't clear if they're arbitrary or related to something real.
# This same methodology could be used for toxins too.
test_obs <- 150
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
         target_cat=case_when(target_sp==0 ~ 1,
                              between(target_sp, 0, 2.5e3) ~ 2,
                              between(target_sp, 2.5e3, 5e3) ~ 3,
                              between(target_sp, 5e3, 20e3) ~ 4,
                              between(target_sp, 20e3, 100e3) ~ 5,
                              between(target_sp, 100e3, 500e3) ~ 6,
                              between(target_sp, 500e3, 1e6) ~ 7,
                              target_sp > 1e6 ~ 8)) %>%
  # mutate(target_sp=karenia_mikimotoi,
  #        target_cat=case_when(target_sp==0 ~ 1,
  #                             between(target_sp, 0, 2.5e3) ~ 2,
  #                             between(target_sp, 2.5e3, 5e3) ~ 3,
  #                             between(target_sp, 5e3, Inf) ~ 4)) %>%
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
         date_lag1=as.numeric(date_collected-lag(date_collected, 1)),
         date_lag2=as.numeric(date_collected-lag(date_collected, 2)),
         date_lag3=as.numeric(date_collected-lag(date_collected, 3)),
         target_sp_lagWt1=target_sp_lag1/date_lag1,
         target_sp_lagWt2=target_sp_lag2/date_lag2,
         target_sp_lagWt3=target_sp_lag3/date_lag3) %>%
  mutate(nObs=n(),
         obs_id=row_number()) %>%
  ungroup %>%
  filter(complete.cases(.))

train.df <- sampling.df %>%
  group_by(site.id) %>%
  filter(obs_id < (max(obs_id)-test_obs)) %>%
  ungroup
test.df <- sampling.df %>%
  group_by(site.id) %>%
  filter(obs_id == (max(obs_id)-test_obs)) %>%
  ungroup %>%
  filter(site.id %in% train.df$site.id)

table(train.df$target_cat)
table(test.df$target_cat)

# Do some analysis here: what are the preceding observations before cat > 2?
sampling.df %>% 
  mutate(cat_0_x=lag(target_cat, 1)==1 & 
           lag(target_cat, 2)==1 & 
           lag(target_cat, 3)==1 & 
           lag(target_cat, 4)==1 & 
           lag(target_cat, 5)==1) %>% 
  ggplot(aes(as.factor(cat_0_x))) + geom_bar() + facet_wrap(~target_cat, scales="free_y")





# ordinal model -----------------------------------------------------------


library(brms)

out <- brm(target_cat ~ 
             
             ydayCos + ydaySin +
             
             short_wave_wk_sc:ydayCos + temp_wk_sc:ydayCos +
             windDir_wk_sc:ydayCos + waterDir_wk_sc:ydayCos +
             
             # target_sp_lagWt1 + target_sp_lagWt2 + target_sp_lagWt3 +
             target_sp_lag1:target_sp_lag2:target_sp_lag3:ydayCos +
             
             target_sp_lagWt1:ydayCos + target_sp_lagWt2:ydayCos + target_sp_lagWt3:ydayCos +
             
             target_sp_lag1:date_lag1 + target_sp_lag2:date_lag2 + target_sp_lag3:date_lag3 +
             target_sp_lag1:date_lag1:ydayCos + target_sp_lag2:date_lag2:ydayCos + target_sp_lag3:date_lag3:ydayCos +
             
             
             target_sp_lagWt1:temp_wk_sc +
             target_sp_lagWt1:short_wave_wk_sc + target_sp_lagWt1:wind_wk_sc +
           
               (1|site.id),
           data=train.df, family=cumulative(), cores=4,
           prior=set_prior("normal(0,2)", class="b"))
saveRDS(out, "temp\\brms_mod.rds")
pp_check(out, type="bars")
summary(out)
conditional_effects(out, categorical=T)
conditional_effects(out)
mcmc_areas(out, regex_pars="b_*")

#out.pred <- posterior_predict(out, newdata=test.df) 
out.pred <- posterior_linpred(out, transform=F, newdata=test.df, allow_new_levels=F)
cutoffs <- summary(out)$fixed[1:3,1]
pred.df <- test.df %>%
  mutate(pred_mn=colMeans(out.pred), 
         pred_se=apply(out.pred, 2, function(x) sd(x)/sqrt(length(x))),
         p1=apply(out.pred, 2, function(x) sum(x<cutoffs[1])/nrow(out.pred)),
         p2=apply(out.pred, 2, function(x) sum(x>cutoffs[1])/nrow(out.pred)),
         p3=apply(out.pred, 2, function(x) sum(x>cutoffs[2])/nrow(out.pred)),
         p4=apply(out.pred, 2, function(x) sum(x>cutoffs[3])/nrow(out.pred)))
# p1=apply(out.pred, 2, function(x) sum(x<cutoffs[1])/nrow(out.pred)),
# p2=apply(out.pred, 2, function(x) sum(x>cutoffs[1] & x<cutoffs[2])/nrow(out.pred)),
# p3=apply(out.pred, 2, function(x) sum(x>cutoffs[2] & x<cutoffs[3])/nrow(out.pred)),
# p4=apply(out.pred, 2, function(x) sum(x>cutoffs[3])/nrow(out.pred)))
ggplot(pred.df, aes(log(target_sp+1), pred_mn, ymin=pred_mn-2*pred_se, ymax=pred_mn+2*pred_se)) + 
  geom_point() + geom_linerange() + stat_smooth(se=F, method="lm")
ggplot(pred.df, aes(target_cat, pred_mn, ymin=pred_mn-2*pred_se, ymax=pred_mn+2*pred_se)) + 
  geom_point() + geom_linerange() + stat_smooth(se=F, method="lm")
ggplot(pred.df, aes(target_cat, p1)) + 
  geom_point() + stat_smooth(se=F) + ylim(0, 1) + xlim(1, max(test.df$target_cat))
ggplot(pred.df, aes(target_cat, p2)) + 
  geom_point() + stat_smooth(se=F) + ylim(0, 1) + xlim(1, max(test.df$target_cat))
ggplot(pred.df, aes(target_cat, p3)) + 
  geom_point() + stat_smooth(se=F) + ylim(0, 1) + xlim(1, max(test.df$target_cat))
ggplot(pred.df, aes(target_cat, p4)) + 
  geom_point() + stat_smooth(se=F) + ylim(0, 1) + xlim(1, max(test.df$target_cat))


plot(colMeans(posterior_linpred(out)), log(1+train.df$target_sp),
     xlab="Fitted mean", ylab="log(N[obs])")


plot(colMeans(posterior_linpred(out, transform=T)), train.df$target_cat,
     xlab="Fitted mean (0-1)", ylab="Observed category")
plot(factor(train.df$target_cat), colMeans(posterior_linpred(out, transform=T)),
     xlab="Observed category", ylab="Fitted mean (0-1)")

out.fit <- posterior_linpred(out, transform=F)
plot(apply(out.fit, 2, function(x) sum(x>cutoffs[3])/nrow(out.fit)), 
     train.df$target_cat,
     xlab="Pr(>2nd)", ylab="Observed category")
fit.df <- train.df %>%
  mutate(pred_mn=colMeans(out.fit), 
         pred_se=apply(out.fit, 2, function(x) sd(x)/sqrt(length(x))),
         p1=apply(out.fit, 2, function(x) sum(x<cutoffs[1])/nrow(out.fit)),
         p2=apply(out.fit, 2, function(x) sum(x>cutoffs[1])/nrow(out.fit)),
         p3=apply(out.fit, 2, function(x) sum(x>cutoffs[2])/nrow(out.fit)),
         p4=apply(out.fit, 2, function(x) sum(x>cutoffs[3])/nrow(out.fit)))
ggplot(fit.df, aes(log(target_sp+1), pred_mn, ymin=pred_mn-2*pred_se, ymax=pred_mn+2*pred_se)) + 
  geom_point() + geom_linerange() + stat_smooth(se=F)
ggplot(fit.df, aes(target_cat, pred_mn, ymin=pred_mn-2*pred_se, ymax=pred_mn+2*pred_se)) + 
  geom_point() + geom_linerange() + stat_smooth(se=F)
ggplot(fit.df, aes(target_cat, p1)) + 
  geom_point() + stat_smooth(se=F) + ylim(0, 1) + xlim(1, max(test.df$target_cat))
ggplot(fit.df, aes(target_cat, p2)) + 
  geom_point() + stat_smooth(se=F) + ylim(0, 1) + xlim(1, max(test.df$target_cat))
ggplot(fit.df, aes(target_cat, p3)) + 
  geom_point() + stat_smooth(se=F) + ylim(0, 1) + xlim(1, max(test.df$target_cat))
ggplot(fit.df, aes(target_cat, p4)) + 
  geom_point() + stat_smooth(se=F) + ylim(0, 1) + xlim(1, max(test.df$target_cat))





dates <- sort(unique(sampling.df$date_collected))[325:380]

for(i in 1:length(dates)) {
  date_i <- dates[i]
  train.df <- sampling.df %>%
    filter(date <= date_i) %>%
    group_by(site.id) %>%
    slice_head(n=-1) %>%
    ungroup
  test.df <- sampling.df %>%
    filter(date <= date_i) %>%
    group_by(site.id) %>%
    slice_tail(n=1) %>%
    ungroup %>%
    filter(site.id %in% train.df$site.id)
  
  if(i == 1) {
    prior_i <- set_prior("normal(0,2)", class="b")
  } else {
    prior_im1 <- readRDS(glue("temp{sep}prior_{dates[i-1]}.rds"))
    post_im1 <- readRDS(glue("temp{sep}out_{dates[i-1]}.rds"))
    b_im1 <- as_draws_df(post_im1, variable="b_", regex=T) %>%
      pivot_longer(cols=starts_with("b_"), names_to="param", values_to="val") %>%
      group_by(param) %>%
      summarise(mn=mean(val), sd=sd(val)) %>%
      mutate(coef=str_remove(param, "b_"))
    eff_im1 <- b_im1 %>% filter(!grepl("Intercept", param))
    int_im1 <- b_im1 %>% filter(grepl("Intercept", param)) %>%
      mutate(coef=str_remove(str_remove(coef, "Intercept\\["), "]"))
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
  
  out <- brm(target_cat ~ 
               
               ydayCos + ydaySin +
               
               # short_wave_wk_sc:ydayCos + 
               # temp_wk_sc:ydayCos +
               # windDir_wk_sc:ydayCos +
               # waterDir_wk_sc:ydayCos +
               
               target_sp_lag1:target_sp_lag2:target_sp_lag3:ydayCos +
               
               target_sp_lagWt1:ydayCos + 
               target_sp_lagWt2:ydayCos + 
               target_sp_lagWt3:ydayCos +
               
               # target_sp_lag1:date_lag1 + 
               # target_sp_lag2:date_lag2 +
               # target_sp_lag3:date_lag3 +
               target_sp_lag1:date_lag1:ydayCos + 
               target_sp_lag2:date_lag2:ydayCos + 
               target_sp_lag3:date_lag3:ydayCos +
               
               # target_sp_lagWt1:temp_wk_sc +
               # target_sp_lagWt1:short_wave_wk_sc + 
               # target_sp_lagWt1:wind_wk_sc +
               
               (1|site.id),
             data=train.df, family=cumulative(), cores=4,
             prior=prior_i)
  out.pred <- posterior_linpred(out, transform=F, newdata=test.df, allow_new_levels=F)
  cutoffs <- c(-Inf, 
               summary(out)$fixed[grep("Intercept", rownames(summary(out)$fixed)),1],
               Inf)
  
  pred.df <- test.df %>%
    mutate(pred_mn=colMeans(out.pred), 
           pred_se=apply(out.pred, 2, function(x) sd(x)/sqrt(length(x)))) %>%
    bind_cols(imap_dfc(setNames(2:length(cutoffs), glue("pr{1:(length(cutoffs)-1)}")), 
                       ~apply(out.pred, 2, 
                              function(x) sum(x>cutoffs[.x-1] & x<cutoffs[.x])/nrow(out.pred))))
  saveRDS(out, glue("temp{sep}out_{date_i}.rds"))
  saveRDS(prior_summary(out), glue("temp{sep}prior_{date_i}.rds"))
  saveRDS(summary(out)$fixed, glue("temp{sep}fixed_eff_{date_i}.rds"))
  saveRDS(pred.df, glue("temp{sep}pred_{date_i}.rds"))
}



























# sampling model ---------------------------------------------------------------

samp.mx <- model.matrix(~ 0 + time_sc + I(time_sc^2) +
                          short_wave_sc + zeta_sc + temp_sc,
                        data=sampling.df)

mod <- cmdstan_model(glue("models{sep}00_sampling_simple.stan"))
sampling_hydro.data <- list(
  N=nrow(sampling.df),
  X=samp.mx,
  nCov=ncol(samp.mx),
  y=sampling.df$target_sp,
  prior_ln_lambda=c(mean(log(sampling.df$target_sp)),
                    sd(log(sampling.df$target_sp)))
)

fit <- mod$sample(data=sampling_hydro.data,
                  chains=4, parallel_chains=4,
                  max_treedepth=20,
                  refresh=500, iter_sampling=2000)
fit$cmdstan_diagnose()
out.sum <- fit$summary()

mcmc_areas(fit$draws("beta")) + scale_y_discrete(labels=colnames(samp.mx))
ppc_dens_overlay(y=log(sampling_hydro.data$y+1),
                 yrep=log(fit$draws("y_pred", format="matrix")[1:100,]+1))

out.df <- sampling.df %>%
  mutate(lnlambda_mean=filter(out.sum, grepl("ln_lambda", variable))$mean,
         lnlambda_q5=filter(out.sum, grepl("ln_lambda", variable))$q5,
         lnlambda_q95=filter(out.sum, grepl("ln_lambda", variable))$q95)
ggplot(out.df, aes(log(target_sp), lnlambda_mean)) + geom_point() +
  geom_linerange(aes(ymin=lnlambda_q5, ymax=lnlambda_q95)) +
  geom_abline()
ggplot(out.df, aes(target_sp, exp(lnlambda_mean))) + geom_point() +
  geom_linerange(aes(ymin=exp(lnlambda_q5), ymax=exp(lnlambda_q95))) +
  geom_abline()

ggplot(out.df, aes(temp_sc, log(target_sp)-lnlambda_mean)) +
  geom_point()
ggplot(out.df, aes(time_sc, log(target_sp)-lnlambda_mean)) +
  geom_point()
ggplot(out.df, aes(short_wave_sc, log(target_sp)-lnlambda_mean)) +
  geom_point()
ggplot(out.df, aes(zeta_sc, log(target_sp)-lnlambda_mean)) +
  geom_point()



# Is this doing *anything* besides shrinking values toward the mean?





