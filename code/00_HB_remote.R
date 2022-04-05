# HABReports Bayesian modelling
# Bayesian HAB model exploration
# Tim Szewczyk


# This script is a translation of parts of docs/02_HB_exploration.Rmd for running remotely



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "ncdf4", "sf", "LaplacesDemon",
          "brms", "bayesplot", "jsonlite")
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
  # mutate(target_sp=karenia_mikimotoi,
  #        target_cat=case_when(target_sp==0 ~ 1,
  #                             between(target_sp, 0, 2.5e3) ~ 2,
  #                             between(target_sp, 2.5e3, 5e3) ~ 3,
  #                             between(target_sp, 5e3, 20e3) ~ 4,
  #                             between(target_sp, 20e3, 100e3) ~ 5,
  #                             between(target_sp, 100e3, 500e3) ~ 6,
  #                             between(target_sp, 500e3, 1e6) ~ 7,
  #                             target_sp > 1e6 ~ 8)) %>%
  mutate(target_sp=karenia_mikimotoi,
         target_cat=case_when(target_sp==0 ~ 1,
                              between(target_sp, 0, 2.5e3) ~ 2,
                              between(target_sp, 2.5e3, 5e3) ~ 3,
                              between(target_sp, 5e3, Inf) ~ 4)) %>%
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






dates <- sort(unique(sampling.df$date_collected))[350:370]

for(i in 2:length(dates)) {
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
  cat("Starting", i, "--", as.character(date_i), "\n")
  cat("  train:", nrow(train.df), "rows with", n_distinct(train.df$site.id), "sites\n")
  cat("  test:", nrow(test.df), "rows with", n_distinct(test.df$site.id), "sites\n")
  
  if(i == 1) {
    prior_i <- set_prior("normal(0,2)", class="b")
  } else {
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
  cat("  Set priors\n")
  prior_i
  
  Sys.sleep(2)
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
  
  Sys.sleep(2)
  cat("  Fitted model \n")
  out.pred <- posterior_linpred(out, transform=F, newdata=test.df, allow_new_levels=T)
  # out_summary <- summary(out)
  # cutoffs <- c(-Inf, 
  #              out_summary$fixed[grep("Intercept", rownames(out_summary$fixed)),1],
  #              Inf)
  
  cat("  Made predictions \n")
  Sys.sleep(2)
  pred.df <- test.df %>%
    mutate(pred_mn=colMeans(out.pred), 
           pred_se=apply(out.pred, 2, function(x) sd(x)/sqrt(length(x))),
           pred_q10=apply(out.pred, 2, function(x) quantile(x, probs=0.1)),
           pred_q90=apply(out.pred, 2, function(x) quantile(x, probs=0.9))) %>%
    # bind_cols(imap_dfc(setNames(2:length(cutoffs), glue("pr{1:(length(cutoffs)-1)}")), 
    #                    ~apply(out.pred, 2, 
    #                           function(x) sum(x>cutoffs[.x-1] & x<cutoffs[.x])/nrow(out.pred)))) %>%
    mutate(target_cat_f=factor(target_cat, levels=1:4))
  
  p <- ggplot(pred.df, aes(log(target_sp+1), pred_mn,
                           ymin=pred_q10, ymax=pred_q90)) + 
    geom_point(aes(colour=target_cat_f)) + geom_errorbar(aes(colour=target_cat_f), width=0.1) + 
    stat_smooth(se=F, method="lm") +
    scale_colour_manual("Obs. cat.", values=c("1"="green3", "2"="gold1", "3"="orange", "4"="red")) +
    labs(title=date_i, x="observed log(N+1)", y="Predicted mean + 80% CI")
  ggsave(glue("temp{sep}plot_pred_{date_i}.png"), p, height=4, width=5, dpi=200)
  
  cat("  Saving output \n")
  saveRDS(out, glue("temp{sep}out_{date_i}.rds"))
  # saveRDS(out_summary$fixed, glue("temp{sep}fixed_eff_{date_i}.rds"))
  saveRDS(pred.df, glue("temp{sep}pred_{date_i}.rds"))
}



