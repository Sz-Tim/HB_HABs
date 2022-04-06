# HABReports Bayesian modelling
# Bayesian HAB model exploration
# Tim Szewczyk


# This script is a translation of parts of docs/02_HB_exploration.Rmd for running remotely



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "LaplacesDemon", "brms")
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
  filter(complete.cases(.)) %>%
  select(site.id, lon, lat, date, year, obs_id, nObs, 
         starts_with("target_"), starts_with("date_"), starts_with("yday"),
         ends_with("_wk_sc")) %>%
  mutate(target_sp_binary=target_cat > 3,
         wk=floor(as.numeric(date - min(date))/7))





# model runs --------------------------------------------------------------

weeks <- sort(unique(sampling.df$wk))[-c(1:110)]
out <- vector("list", length(weeks))

for(i in 1:length(weeks)) {
  week_i <- weeks[i]
  train.df <- sampling.df %>%
    filter(wk < week_i)
  test.df <- sampling.df %>%
    filter(wk == week_i) %>%
    filter(site.id %in% train.df$site.id)
  date_range <- paste(as.character(range(test.df$date)), collapse=" to ")
  cat("Starting", i, "= week", as.character(week_i), "=", date_range, "\n")
  cat("  train:", nrow(train.df), "rows with", n_distinct(train.df$site.id), "sites\n")
  cat("  test:", nrow(test.df), "rows with", n_distinct(test.df$site.id), "sites\n")
  
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
  cat("  Set priors\n")
  
  out[[i]] <- brm(target_cat ~
               
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
                  data=train.df, cores=4,  family=cumulative(),
                  prior=prior_i, file=glue("temp{sep}out_{week_i}_ord"),
                  save_model=glue("temp{sep}mod_{week_i}_ord.stan"))
  
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
# 
# coef.df <- dir("temp", "fixed") %>% 
#   map_dfr(~readRDS(paste0("temp\\", .x)) %>% 
#             rownames_to_column("param") %>%
#             mutate(week=as.numeric(str_sub(.x, 11, -5))))
# 
# ggplot(coef.df, aes(week, y=Estimate, ymin=`l-95% CI`, ymax=`u-95% CI`)) + 
#   geom_hline(yintercept=0) + 
#   geom_ribbon(alpha=0.25) + 
#   geom_line() + facet_wrap(~param)









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

