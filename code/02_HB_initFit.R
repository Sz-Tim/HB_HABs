# HABReports Bayesian modelling
# Bayesian HAB model initial fits
# Tim Szewczyk


# This script fits models with data from 2013-2017 using horseshoe priors



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "LaplacesDemon", "brms", 
          "randomForest", "caret")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
walk(dir("code", "000_fn", full.names=T), source)

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





# covariates --------------------------------------------------------------

species <- c("alexandrium_sp", "dinophysis_sp", "karenia_mikimotoi",
             "prorocentrum_lima", "pseudo_nitzschia_sp")

sampling.df <- read_csv(glue("data{sep}sampling_local.csv"))

thresh.df <- read_csv(glue("data{sep}hab_tf_thresholds.csv")) %>%
  filter(!is.na(tl)) %>%
  group_by(hab_parameter, tl) %>%
  slice_head(n=1) %>%
  ungroup %>%
  mutate(bloom=as.numeric(!is.na(alert)))

hydro.df <- dir("data", "hydro_", full.names=T) %>% 
  map(~read_csv(.x, show_col_types=F)) %>%
  reduce(., full_join) %>% 
  pivot_longer(c(contains("_L"), contains("_R")), names_to="param", values_to="value") %>%
  mutate(parType=str_sub(param, 1, -4),
         res=str_sub(param, -2, -2),
         lag=str_sub(param, -1, -1)) %>%
  select(-param) %>%
  pivot_wider(names_from="parType", values_from="value") %>%
  mutate(water=sqrt(u^2 + v^2 + ww^2),
         wind=sqrt(uwind_speed^2 + vwind_speed^2),
         waterDir=atan2(v, u),
         windDir=atan2(vwind_speed, uwind_speed),
         timespan=if_else(lag=="0", "0", "wk")) %>% 
  select(-u, -v, -ww, -vwind_speed, -uwind_speed) %>%
  group_by(obs.id, res, timespan) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm=T))) %>%
  ungroup %>%
  pivot_wider(names_from=c(res, timespan), values_from=4:ncol(.),
              names_glue="{.value}_{res}_{timespan}")
connect.df <- read_csv(glue("data{sep}connectivity_5e3parts.csv")) %>%
  group_by(site.id) %>%
  mutate(influx_wk=zoo::rollsum(influx, 7, align='right', fill="extend")) %>%
  mutate(across(starts_with("influx"), ~log(.x+1)))
cprn.df <- read_csv("data/copernicus_westcoms1.csv") %>%
  mutate(site.id=site.id_v1) %>%
  select(site.id, date, attn_wk, chl_wk, dino_wk, o2_wk, ph_wk, po4_wk)
# >0.85 cor(res,time) = temp, salinity
# >0.85 cor(res) = short_wave, km, precip, wind, windDir
# >0.85 cor(time) = water
# waterDir is more variable across all 4
covars <- c("temp_L_wk", "salinity_L_wk", "short_wave_L_wk", "km_L_wk",
            "precip_L_wk", 
            "wind_L_wk", "windDir_L_wk",
            "water_L_wk", "waterDir_L_wk", 
            "water_R_wk", "waterDir_R_wk",
            "influx_wk", 
            "attn_wk", "chl_wk", "dino_wk", "o2_wk", "ph_wk", "po4_wk")

predictors_int <- c(
  "ydayCos", "ydaySin", 
  
  "tempLwk:ydayCos:ydaySin",
  "salinityLwk:ydayCos:ydaySin",
  "shortwaveLwk:ydayCos:ydaySin",
  "kmLwk:ydayCos:ydaySin",
  "precipLwk:ydayCos:ydaySin",
  "waterLwk:waterDirLwk:ydayCos:ydaySin",
  "waterRwk:waterDirRwk:ydayCos:ydaySin",
  "windLwk:windDirLwk:ydayCos:ydaySin",
  
  "fetch:ydayCos:ydaySin", 
  "influxwk:ydayCos:ydaySin",
  
  "attnwk:ydayCos:ydaySin",
  "chlwk:ydayCos:ydaySin",
  "dinowk:ydayCos:ydaySin",
  "o2wk:ydayCos:ydaySin",
  "phwk:ydayCos:ydaySin",
  "po4wk:ydayCos:ydaySin",
  
  "mo(NcatF1):ydayCos:ydaySin",
  "mo(NcatF2):ydayCos:ydaySin",
  
  "NlnWt1:ydayCos:ydaySin",
  "NlnWt2:ydayCos:ydaySin", 
  
  "mo(NcatF1):mo(NcatF2)",
  
  "windLwk:windDirLwk:fetch",
  "waterLwk:waterDirLwk:fetch",
  "waterRwk:waterDirRwk:fetch",
  
  "(1|siteid)"
)

predictors_s <- c(
  "tempLwk", "salinityLwk", "shortwaveLwk",
  "windVel", "waterVelL", "waterVelR",
  "influxwk", "fetch",
  "attnwk", "chlwk", "dinowk", "o2wk", "phwk", "po4wk",
  "Nbloom1", "Nbloom2"
)
predictors_s <- c(
  "tempLwk", "phwk", "Nbloom1", "Nbloom2"
)

s_b <- glue("b{predictors_s}") 
s_flist <- map(s_b, ~as.formula(paste0(.x, "~s(ydayCos,ydaySin) + (1|siteid)")))
s_form_ord <- bf(
  glue("NcatNum | thres(3) ~ bydayC*ydayCos + bydayS*ydaySin + bydaySC*ydaySC +",
       "{paste(s_b, predictors_s, sep='*', collapse='+')}"),
  bydayC ~ 1 + (1|siteid),
  bydayS ~ 1 + (1|siteid),
  bydaySC ~ 1 + (1|siteid),
  flist=s_flist,
  nl=T)
s_form_bern <- bf(
  glue("Nbloom ~ bydayC*ydayCos + bydayS*ydaySin + bydaySC*ydaySC +",
       "{paste(s_b, predictors_s, sep='*', collapse='+')}"),
  bydayC ~ 1 + (1|siteid),
  bydayS ~ 1 + (1|siteid),
  bydaySC ~ 1 + (1|siteid),
  flist=s_flist,
  nl=T)

form_ordinal <- bf(glue("NcatNum | thres(3) ~ {paste(predictors_int, collapse=' + ')}"))
form_bern <- bf(glue("Nbloom ~ {paste(predictors_int, collapse=' + ')}"))
form_bern_noCatF <- bf(glue("Nbloom ~ {paste(grep('catF_1', predictors_int, value=T, invert=T), collapse=' + ')}"))

priors <- c(prior(horseshoe(3, par_ratio=0.2), class="b"),
            prior(normal(0, 1), class="Intercept"))
s_priors <- map(s_b, 
                ~prior_string(prior="normal(0, 1)", nlpar=.x)) %>% 
  do.call('c', .) %>%
  c(., prior(normal(0, 1), nlpar="bydayC"), 
    prior(normal(0, 1), nlpar="bydayS"),
    prior(normal(0, 1), nlpar="bydaySC"))

# Model details
ctrl <- list(adapt_delta=0.95, max_treedepth=20)
chains <- 4
iter <- 100
warmup <- iter/2
refresh <- 100





# initial fit -------------------------------------------------------------

out.bern01 <- out.bern11 <- out.sbern01 <- out.sbern11 <- out.ord <- vector("list", length(species))

for(sp in 1:length(species)) {
  target <- species[sp]
  target.tf <- thresh.df %>% filter(hab_parameter==target)
  
  target.df <- sampling.df %>%
    filter(grid==1) %>% 
    rename(N=!!target) %>%
    select(obs.id, site.id, date, hour, grid, lon, lat, fetch, bearing, N) %>%
    mutate(yday=yday(date),
           ydayCos=cos(2*pi*yday/365),
           ydaySin=sin(2*pi*yday/365),
           year=year(date),
           bearing=bearing*pi/180,
           N=round(N),
           N.ln=log(N+1),
           N.PA=as.numeric(N>0)) %>%
    rowwise() %>%
    mutate(N.cat=target.tf$tl[max(which(N >= target.tf$min_ge))]) %>%
    ungroup %>%
    mutate(N.catF=factor(N.cat, levels=unique(target.tf$tl), ordered=T),
           N.catNum=as.numeric(N.catF),
           N.bloom=target.tf$bloom[match(N.cat, target.tf$tl)]) %>%
    arrange(site.id, date) %>%
    group_by(site.id) %>%
    multijetlag(N.ln, N.PA, N.cat, N.catF, N.bloom, date, n=2) %>%
    ungroup %>%
    mutate(across(starts_with("date_"), ~as.numeric(date-.x)),
           # I don't love this since small if N.ln_x is small OR date_x is large
           N.lnWt_1=N.ln_1/date_1,
           N.lnWt_2=N.ln_2/date_2) %>%
    full_join(hydro.df) %>%
    full_join(connect.df) %>%
    full_join(cprn.df) %>%
    mutate(across(contains("Dir_"), ~cos(.x-bearing))) %>%
    mutate(across(one_of(grep("Dir", covars, invert=T, value=T)), CenterScale)) %>%
    arrange(site.id, date) %>%
    filter(complete.cases(.)) %>%
    select(site.id, lon, lat, date, year, obs.id, fetch, bearing,
           starts_with("N"), starts_with("date_"), starts_with("yday"),
           one_of(covars)) %>%
    rename_with(~str_remove_all(.x, "\\.|_")) %>%
    mutate(ydaySC=ydaySin*ydayCos,
           windVel=windLwk*windDirLwk,
           waterVelL=waterLwk*waterDirLwk,
           waterVelR=waterRwk*waterDirRwk)
  
  write_csv(target.df, glue("out{sep}test_full{sep}dataset_{target}.csv"))
  
  train.df <- target.df %>% filter(year <= 2017)
  test.df <- target.df %>% filter(year > 2017)
  
  out.ord[[sp]] <- brm(form_ordinal, data=train.df,
                       family=cumulative("probit"), prior=priors, 
                       iter=iter, warmup=warmup, refresh=refresh, init=0,
                       control=ctrl, chains=chains, cores=chains,
                       file=glue("out{sep}test_full{sep}ord_{target}"))
  out.sord[[sp]] <- brm(s_form_ord, data=train.df, 
                        family=cumulative("probit"), prior=s_priors, 
                        iter=iter, warmup=warmup, refresh=refresh, init=0,
                        control=ctrl, chains=chains, cores=chains,
                        file=glue("out{sep}test_full{sep}sord_{target}"))
  if(n_distinct(filter(train.df, N.bloom_1==0)$N.catF_1)==1) {
    form_01 <- form_bern_noCatF
  } else {
    form_01 <- form_bern
  }
  out.bern01[[sp]] <- brm(form_01, data=train.df %>% filter(Nbloom1==0),
                          family=bernoulli("probit"), prior=priors, 
                          iter=iter, warmup=warmup, refresh=refresh, init=0,
                          control=ctrl, chains=chains, cores=chains,
                          file=glue("out{sep}test_full{sep}bern01_{target}"))
  out.bern11[[sp]] <- brm(form_bern, data=train.df %>% filter(Nbloom1==1),
                          family=bernoulli("probit"), prior=priors, 
                          iter=iter, warmup=warmup, refresh=refresh, init=0,
                          control=ctrl, chains=chains, cores=chains,
                          file=glue("out{sep}test_full{sep}bern11_{target}"))
  out.sbern01[[sp]] <- brm(s_form_bern, data=train.df %>% filter(Nbloom1==0), 
                           family=bernoulli("probit"), prior=s_priors, 
                           iter=iter, warmup=warmup, refresh=refresh, init=0,
                           control=ctrl, chains=chains, cores=chains,
                           file=glue("out{sep}test_full{sep}sbern01_{target}"))
  out.sbern11[[sp]] <- brm(s_form_bern, data=train.df %>% filter(Nbloom1==1), 
                           family=bernoulli("probit"), prior=s_priors, 
                           iter=iter, warmup=warmup, refresh=refresh, init=0,
                           control=ctrl, chains=chains, cores=chains,
                           file=glue("out{sep}test_full{sep}sbern11_{target}"))
  
  # RF
  train.rf <- train.df %>%
    select(Nbloom, Nbloom1,
           lon, lat, 
           ydaySin, ydayCos, fetch,
           Nln1, Ncat1, 
           Nln2, Ncat2, 
           NlnWt1, NlnWt2, 
           one_of(covars)) %>%
    mutate(Nbloom=factor(Nbloom)) %>%
    as.data.frame()
  test.rf <- test.df %>%
    select(Nbloom, Nbloom1,
           lon, lat, 
           ydaySin, ydayCos, fetch,
           Nln1, Ncat1, 
           Nln2, Ncat2, 
           NlnWt1, NlnWt2, 
           one_of(covars)) %>%
    mutate(Nbloom=factor(Nbloom)) %>%
    as.data.frame()
  train.rf01 <- train.rf %>% filter(Nbloom1==0)
  train.rf11 <- train.rf %>% filter(Nbloom1==1)
  test.rf01 <- test.rf %>% filter(Nbloom1==0)
  test.rf11 <- test.rf %>% filter(Nbloom1==1)
  rf <- randomForest(Nbloom ~ ., data=train.rf, 
                     xtest=test.rf[,-1], ytest=test.rf[,1], 
                     proximity=T, keep.forest=T)
  p.rf <- predict(rf, test.rf, type="prob")
  
  rf.01 <- randomForest(Nbloom ~ ., 
                        data=train.rf01, 
                        xtest=test.rf01[,-1], 
                        ytest=test.rf01[,1], 
                        proximity=T, keep.forest=T)
  p.rf.01 <- predict(rf.01, test.rf01, type="prob")
  
  rf.11 <- randomForest(Nbloom ~ ., 
                        data=train.rf11, 
                        xtest=test.rf11[,-1], 
                        ytest=test.rf11[,1], 
                        proximity=T, keep.forest=T)
  p.rf.11 <- predict(rf.11, test.rf11, type="prob")
  
  pred.df <- test.df
  pred.ord <- posterior_epred(out.ord[[sp]], 
                              newdata=test.df, 
                              allow_new_levels=T)
  cat.iter <- apply(pred.ord, 1:2, which.max)
  pred.df <- pred.df %>%
    mutate(ord_mnpr0=c(colMeans(pred.ord[,,1,drop=F])),
           ord_mnpr1=c(colMeans(pred.ord[,,2,drop=F])),
           ord_mnpr2=c(colMeans(pred.ord[,,3,drop=F])),
           ord_mnpr3=c(colMeans(pred.ord[,,4,drop=F])),
           ord_prmax0=colMeans(cat.iter==1),
           ord_prmax1=colMeans(cat.iter==2),
           ord_prmax2=colMeans(cat.iter==3),
           ord_prmax3=colMeans(cat.iter==4),
           ord_mnCat=colMeans(cat.iter),
           ord_wtmnCat=ord_mnpr0 + 2*ord_mnpr1 + 3*ord_mnpr2 + 4*ord_mnpr3,
           rf_mnpr=p.rf[,2])
  pred.sord <- posterior_epred(out.sord[[sp]], 
                              newdata=test.df, 
                              allow_new_levels=T)
  cat.iter <- apply(pred.sord, 1:2, which.max)
  pred.df <- pred.df %>%
    mutate(sord_mnpr0=c(colMeans(pred.sord[,,1,drop=F])),
           sord_mnpr1=c(colMeans(pred.sord[,,2,drop=F])),
           sord_mnpr2=c(colMeans(pred.sord[,,3,drop=F])),
           sord_mnpr3=c(colMeans(pred.sord[,,4,drop=F])),
           sord_prmax0=colMeans(cat.iter==1),
           sord_prmax1=colMeans(cat.iter==2),
           sord_prmax2=colMeans(cat.iter==3),
           sord_prmax3=colMeans(cat.iter==4),
           sord_mnCat=colMeans(cat.iter),
           sord_wtmnCat=sord_mnpr0 + 2*sord_mnpr1 + 3*sord_mnpr2 + 4*sord_mnpr3)
  pred.bern01 <- posterior_epred(out.bern01[[sp]], 
                                 newdata=test.df %>% filter(Nbloom1==0) %>% droplevels, 
                                 allow_new_levels=T) 
  pred.sbern01 <- posterior_epred(out.sbern01[[sp]], 
                                 newdata=test.df %>% filter(Nbloom1==0) %>% droplevels, 
                                 allow_new_levels=T) 
  pred.01 <- test.df %>% filter(Nbloom1==0) %>% select(obsid) %>%
    mutate(bern_mnpr=colMeans(pred.bern01),
           sbern_mnpr=colMeans(pred.sbern01))
  pred.bern11 <- posterior_epred(out.bern11[[sp]], 
                                 newdata=test.df %>% filter(Nbloom1==1) %>% droplevels, 
                                 allow_new_levels=T) 
  pred.sbern11 <- posterior_epred(out.sbern11[[sp]], 
                                 newdata=test.df %>% filter(Nbloom1==1) %>% droplevels, 
                                 allow_new_levels=T) 
  pred.11 <- test.df %>% filter(N.bloom_1==1) %>% select(obs.id) %>%
    mutate(bern_mnpr=colMeans(pred.bern11),
           sbern_mnpr=colMeans(pred.sbern11))
  pred.split <- bind_rows(pred.01, pred.11) %>%
    mutate(rf_split_mnpr=c(p.rf.01[,2], p.rf.11[,2]))
  pred.df <- full_join(pred.df, pred.split)
  
  write_csv(pred.df, glue("out{sep}test_full{sep}pred_{target}.csv"))
  
  cat("Finished", target, "\n")
}



# 
# 
# # NOT RUN -----------------------------------------------------------------
# 

# # library(lisa)
# cm <- readRDS("../../fishCircadianTemps/figs/cmr_cmaps.RDS")
# # mod.col <- rev(lisa$LeonardodaVinci)
# mod.col <- c(cm$waterlily[length(cm$waterlily) * c(0.2, 0.35, 0.65, 0.8)], "black")
# tl.col <- c("green3", "gold1", "orange", "red")
# sams.cols <- c("alexandrium_sp"="#2A5883",
#                "dinophysis_sp"="#46BFDE",
#                "karenia_mikimotoi"="#A9DAE0",
#                "prorocentrum_lima"="#E77334",
#                "pseudo_nitzschia_sp"="#D21E4C")
# 
# data.df <- dir("out", "^dataset_", full.names=T) %>%
#   imap_dfr(~read_csv(.x, show_col_types=F) %>%
#              mutate(species=str_remove(str_remove(.x, "out/dataset_"), ".csv")))
# pred.f <- dir(glue("out{sep}initFit"), "^pred_")
# pred.df <- map_dfr(pred.f, ~read_csv(glue("out{sep}initFit{sep}{.x}"), show_col_types=F) %>%
#                      select(obs.id, N.bloom, N.catF, N.catNum, N.bloom_1,
#                             contains("_mnpr")) %>%
#                      mutate(species=str_remove(str_remove(.x, "pred_"), ".csv"),
#                             bloomThresh=max((!N.bloom)*N.catNum))) %>%
#   pivot_longer(starts_with("ord"), names_to="ord_cat", values_to="ord_pr") %>%
#   mutate(ord_cat=as.numeric(str_remove(ord_cat, "ord_mnpr"))) %>%
#   filter(ord_cat >= bloomThresh) %>%
#   group_by(species, obs.id) %>%
#   summarise(ord_pr=sum(ord_pr),
#             across(everything(), ~first(.x))) %>%
#   mutate(avg_pr=(bern_mnpr + ord_pr + rf_mnpr + rf_split_mnpr)/4)
# 
# pred.df %>%
#   pivot_longer(ends_with("pr"), names_to="model", values_to="pr_bloom") %>%
#   mutate(model=factor(model, 
#                       levels=c("ord_pr", "bern_mnpr", "rf_mnpr", "rf_split_mnpr", "avg_pr"),
#                       labels=c("Ordinal", "Dual logistic", "RF", "Dual RF", "Average"))) %>%
#   ggplot(aes(pr_bloom, N.bloom, colour=model)) + xlim(0,1) +
#   geom_jitter(alpha=0.25, width=0, height=0.05) +
#   stat_smooth(method="glm", method.args=list(family="binomial"), fullrange=T, se=F) +
#   facet_wrap(~species) +
#   labs(x="Predicted pr(bloom)", y="Bloom") +
#   scale_colour_manual("Model", values=mod.col) +
#   theme_classic() +
#   theme(legend.position=c(0.85, 0.2))
# ggsave(glue("figs{sep}performance{sep}pr_v_obs_bloom.png"),
#        width=8, height=5, dpi=300)
# 
# pred.df %>%
#   pivot_longer(ends_with("pr"), names_to="model", values_to="pr_bloom") %>%
#   mutate(model=factor(model, 
#                       levels=c("ord_pr", "bern_mnpr", "rf_mnpr", "rf_split_mnpr", "avg_pr"),
#                       labels=c("Ordinal", "Dual logistic", "RF", "Dual RF", "Average"))) %>%
#   left_join(., sampling.df %>% select(obs.id, site.id, date)) %>%
#   mutate(month=rollback(date, roll_to_first=T)) %>%
#   ggplot(aes(date, pr_bloom, colour=as.factor(N.bloom), shape=as.factor(N.bloom))) +
#   geom_hline(yintercept=0, colour="grey", size=0.5) +
#   geom_point(alpha=0.5, size=0.9) + 
#   scale_colour_manual("Bloom", values=c("grey50", "red3"), 
#                       guide=guide_legend(override.aes=list(size=1, alpha=1))) +
#   scale_shape_manual("Bloom", values=c(1, 19)) +
#   facet_grid(species~model) +
#   theme_bw() +
#   theme(axis.title.x=element_blank(), panel.grid=element_blank()) +
#   scale_y_continuous(breaks=c(0, 0.5, 1), limits=c(0, 1)) +
#   labs(y="Pr(bloom)", title="All observations")
# ggsave(glue("figs{sep}performance{sep}pr_v_obs_bloom_dates.png"),
#        width=10, height=7.25, dpi=300)
# 
# pred.df %>%
#   filter(N.bloom_1==0) %>%
#   pivot_longer(ends_with("pr"), names_to="model", values_to="pr_bloom") %>%
#   mutate(model=factor(model, 
#                       levels=c("ord_pr", "bern_mnpr", "rf_mnpr", "rf_split_mnpr", "avg_pr"),
#                       labels=c("Ordinal", "Dual logistic", "RF", "Dual RF", "Average"))) %>%
#   left_join(., sampling.df %>% select(obs.id, site.id, date)) %>%
#   mutate(month=rollback(date, roll_to_first=T)) %>%
#   ggplot(aes(date, pr_bloom, colour=as.factor(N.bloom), shape=as.factor(N.bloom))) +
#   geom_hline(yintercept=0, colour="grey", size=0.5) +
#   geom_point(alpha=0.5, size=0.9) + 
#   scale_colour_manual("Bloom", values=c("grey50", "red3"), 
#                       guide=guide_legend(override.aes=list(size=1, alpha=1))) +
#   scale_shape_manual("Bloom", values=c(1, 19)) +
#   facet_grid(species~model) +
#   theme_bw() +
#   theme(axis.title.x=element_blank(), panel.grid=element_blank()) +
#   scale_y_continuous(breaks=c(0, 0.5, 1), limits=c(0, 1)) +
#   labs(y="Pr(bloom)", title="Bloom initiation")
# ggsave(glue("figs{sep}performance{sep}pr_v_obs_bloom_dates_init.png"),
#        width=10, height=7.25, dpi=300)
# 
# pred.df %>%
#   pivot_longer(ends_with("pr"), names_to="model", values_to="pr_bloom") %>%
#   mutate(model=factor(model, 
#                       levels=c("ord_pr", "bern_mnpr", "rf_mnpr", "rf_split_mnpr", "avg_pr"),
#                       labels=c("Ordinal", "Dual logistic", "RF", "Dual RF", "Average"))) %>%
#   left_join(., sampling.df %>% select(obs.id, site.id, date)) %>%
#   mutate(month=rollback(date, roll_to_first=T)) %>%
#   group_by(species, date, N.bloom, model) %>%
#   summarise(pr_mn=mean(pr_bloom), 
#             pr_lo=quantile(pr_bloom, 0.05),
#             pr_hi=quantile(pr_bloom, 0.95)) %>%
#   ggplot(aes(date, pr_mn, colour=as.factor(N.bloom), shape=as.factor(N.bloom))) +
#   geom_hline(yintercept=0, colour="grey", size=0.5) +
#   geom_point(alpha=0.5, size=0.9) + 
#   scale_colour_manual("Bloom", values=c("grey50", "red3"), 
#                       guide=guide_legend(override.aes=list(size=1, alpha=1))) +
#   scale_shape_manual("Bloom", values=c(1, 19)) +
#   facet_grid(species~model) +
#   theme_bw() +
#   theme(axis.title.x=element_blank(), panel.grid=element_blank()) +
#   scale_y_continuous(breaks=c(0, 0.5, 1), limits=c(0, 1)) +
#   labs(y="Pr(bloom)", title="All observations")
# ggsave(glue("figs{sep}performance{sep}pr_v_obs_bloom_dates_dayMeans.png"),
#        width=10, height=7.25, dpi=300)
# 
# 
# pred.df %>% 
#   filter(N.bloom_1==0) %>%
#   pivot_longer(ends_with("pr"), names_to="model", values_to="pr_bloom") %>%
#   mutate(pred_round=((pr_bloom*100) %/% 5 * 5)/100) %>%
#   mutate(model=factor(model, 
#                       levels=c("ord_pr", "bern_mnpr", "rf_mnpr", "rf_split_mnpr", "avg_pr"),
#                       labels=c("Ordinal", "Dual logistic", "RF", "Dual RF", "Average"))) %>%
#   group_by(species, model, pred_round) %>%
#   summarise(obs_bloom=mean(N.bloom), N=n()) %>%
#   ggplot(aes(pred_round, obs_bloom, colour=model)) + 
#   geom_abline(linetype=2) + 
#   geom_point(alpha=0.5, aes(size=N)) + 
#   stat_smooth(method="lm", se=F) + 
#   xlim(0,1) + ylim(0,1) +
#   scale_colour_manual("Model", values=mod.col) +
#   scale_size_continuous(breaks=c(1, 100, 500)) +
#   facet_wrap(~species) +
#   theme_classic() + 
#   theme(legend.position="bottom") + 
#   labs(title="Bloom initiation", x="Mean prediction", y="Proportion observed blooms")
# ggsave(glue("figs{sep}performance{sep}pred_v_obs_bloom_init.png"),
#        width=8, height=6, dpi=300)
# 
# pred.df %>%
#   filter(N.bloom_1==0) %>%
#   pivot_longer(ends_with("pr"), names_to="model", values_to="pr_bloom") %>%
#   mutate(model=factor(model, 
#                       levels=c("ord_pr", "bern_mnpr", "rf_mnpr", "rf_split_mnpr", "avg_pr"),
#                       labels=c("Ordinal", "Dual logistic", "RF", "Dual RF", "Average"))) %>%
#   left_join(., sampling.df %>% select(obs.id, site.id, date)) %>%
#   mutate(month=rollback(date, roll_to_first=T)) %>%
#   group_by(species, date, N.bloom, model) %>%
#   summarise(pr_mn=mean(pr_bloom), 
#             pr_lo=quantile(pr_bloom, 0.05),
#             pr_hi=quantile(pr_bloom, 0.95)) %>%
#   ggplot(aes(date, pr_mn, colour=as.factor(N.bloom), shape=as.factor(N.bloom))) +
#   geom_hline(yintercept=0, colour="grey", size=0.5) +
#   geom_point(alpha=0.5, size=0.9) + 
#   scale_colour_manual("Bloom", values=c("grey50", "red3"), 
#                       guide=guide_legend(override.aes=list(size=1, alpha=1))) +
#   scale_shape_manual("Bloom", values=c(1, 19)) +
#   facet_grid(species~model) +
#   theme_bw() +
#   theme(axis.title.x=element_blank(), panel.grid=element_blank()) +
#   scale_y_continuous(breaks=c(0, 0.5, 1), limits=c(0, 1)) +
#   labs(y="Pr(bloom)", title="Bloom initiation")
# ggsave(glue("figs{sep}performance{sep}pr_v_obs_bloom_dates_dayMeans.png"),
#        width=10, height=7.25, dpi=300)
# 
# pred.df %>% 
#   pivot_longer(ends_with("pr"), names_to="model", values_to="pr_bloom") %>%
#   mutate(pred_round=((pr_bloom*100) %/% 5 * 5)/100) %>%
#   mutate(model=factor(model, 
#                       levels=c("ord_pr", "bern_mnpr", "rf_mnpr", "rf_split_mnpr", "avg_pr"),
#                       labels=c("Ordinal", "Dual logistic", "RF", "Dual RF", "Average"))) %>%
#   group_by(species, model, pred_round) %>%
#   summarise(obs_bloom=mean(N.bloom), N=n()) %>%
#   ggplot(aes(pred_round, obs_bloom, colour=model)) + 
#   geom_abline(linetype=2) + 
#   geom_point(alpha=0.5, aes(size=N)) + 
#   stat_smooth(method="lm", se=F) + 
#   xlim(0,1) + ylim(0,1) +
#   scale_colour_manual("Model", values=mod.col) +
#   scale_size_continuous(breaks=c(1, 100, 500)) +
#   facet_wrap(~species) +
#   theme_classic() + 
#   theme(legend.position="bottom") + 
#   labs(title="All predictions", x="Mean prediction", y="Proportion observed blooms")
# ggsave(glue("figs{sep}performance{sep}pred_v_obs.png"),
#        width=8, height=6, dpi=300)
# 
# pred.df %>%
#   pivot_longer(ends_with("pr"), names_to="model", values_to="pr_bloom") %>%
#   mutate(model=factor(model, 
#                       levels=c("ord_pr", "bern_mnpr", "rf_mnpr", "rf_split_mnpr", "avg_pr"),
#                       labels=c("Ordinal", "Dual logistic", "RF", "Dual RF", "Average"))) %>%
#   ggplot(aes(model, pr_bloom, fill=factor(N.bloom))) +
#   geom_boxplot() + facet_grid(N.bloom_1~species)
# 
# ggplot(pred.df, aes(as.factor(N.bloom_1), bern_mnpr, fill=N.catF)) +
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   scale_fill_manual(values=tl.col) +
#   ggtitle("Dual bernoulli")
# ggplot(pred.df, aes(as.factor(N.bloom_1), ord_pr, fill=N.catF)) +
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   scale_fill_manual(values=tl.col) +
#   ggtitle("Ordinal")
# ggplot(pred.df, aes(as.factor(N.bloom_1), rf_mnpr, fill=N.catF)) +
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   scale_fill_manual(values=tl.col) +
#   ggtitle("RF")
# ggplot(pred.df, aes(as.factor(N.bloom_1), rf_split_mnpr, fill=N.catF)) +
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   scale_fill_manual(values=tl.col) +
#   ggtitle("RF split")
# ggplot(pred.df, aes(as.factor(N.bloom_1), avg_pr, fill=N.catF)) +
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   scale_fill_manual(values=tl.col) +
#   ggtitle("Avg")
# 
# 
# ggplot(pred.df, aes(as.factor(N.bloom_1), bern_mnpr, fill=factor(N.bloom))) +
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   ggtitle("Dual bernoulli")
# ggplot(pred.df, aes(as.factor(N.bloom_1), ord_pr, fill=factor(N.bloom))) +
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   ggtitle("Ordinal")
# ggplot(pred.df, aes(as.factor(N.bloom_1), rf_mnpr, fill=factor(N.bloom))) +
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   ggtitle("RF")
# ggplot(pred.df, aes(as.factor(N.bloom_1), rf_mnpr, fill=factor(N.bloom))) +
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   ggtitle("RF split")
# ggplot(pred.df, aes(as.factor(N.bloom_1), avg_pr, fill=factor(N.bloom))) +
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   ggtitle("Avg")
# 
# 
# binary.df <- pred.df %>% #filter(N.bloom_1==0) %>%
#   select(obs.id, species, ord_pr, bern_mnpr, rf_mnpr, rf_split_mnpr, avg_pr, N.bloom) %>%
#   group_by(species) %>%
#   group_split()
# 
# roc.ls <- list(
#   ord=map(binary.df, ~pROC::roc(.x$N.bloom ~ .x$ord_pr)),
#   bern=map(binary.df, ~pROC::roc(.x$N.bloom ~ .x$bern_mnpr)),
#   rf=map(binary.df, ~pROC::roc(.x$N.bloom ~ .x$rf_mnpr)),
#   rf_s=map(binary.df, ~pROC::roc(.x$N.bloom ~ .x$rf_split_mnpr)),
#   avg=map(binary.df, ~pROC::roc(.x$N.bloom ~ .x$avg_pr))
# )
# 
# for(i in 1:5) {
#   par(mfrow=c(1,1)) 
#   png(glue("figs{sep}performance{sep}01_ROC_{species[i]}.png"), width=5, height=5, res=300, units="in")
#   plot(NA, NA, xlim=c(0, 1), ylim=c(0, 1), 
#        xlab="1 - Specificity", ylab="Sensitivity", axes=F, main=species[i])
#   axis(side=1, at=c(0, 0.5, 1))
#   axis(side=2, at=c(0, 0.5, 1))
#   abline(a=0, b=1, col="grey30", lwd=0.5)
#   map(1:5, 
#       ~lines(1-roc.ls[[.x]][[i]]$specificities, roc.ls[[.x]][[i]]$sensitivities,
#              col=mod.col[.x], lwd=2.5))
#   legend("bottomright", 
#          paste(c("Ordinal", "Dual logistic", "RF", "Dual RF", "Average"),
#                map_dbl(roc.ls, ~round(.x[[i]]$auc, 2)),  sep=": "),
#          col=mod.col, lty=1, lwd=2, bty="n", cex=0.95)
#   dev.off()
# }
# 
# auc.df <- tibble(
#   bern=map(roc.ls[["bern"]], ~.x$auc) %>% do.call(c, .) %>% round(3),
#   ord=map(roc.ls[["ord"]], ~.x$auc) %>% do.call(c, .) %>% round(3),
#   rf=map(roc.ls[["rf"]], ~.x$auc) %>% do.call(c, .) %>% round(3),
#   rf_s=map(roc.ls[["rf_s"]], ~.x$auc) %>% do.call(c, .) %>% round(3),
#   avg=map(roc.ls[["avg"]], ~.x$auc) %>% do.call(c, .) %>% round(3)
# ) %>%
#   mutate(species=species) %>%
#   pivot_longer(1:5, names_to="model", values_to="AUC")
# write_csv(auc.df, glue("figs{sep}performance{sep}auc_all.csv"))
# 
# par(mfrow=c(1,1))
# 
# 
# 
# 
# 
# auc.f <- dir(glue("figs{sep}performance{sep}"), "auc_.*csv", full.names=T)
# auc.df <- map_dfr(auc.f, ~read_csv(.x, show_col_types=F) %>%
#                     mutate(type=str_split_fixed(.x, "auc_", 2)[,2])) %>%
#   mutate(data_subset=if_else(grepl("01", type), "01", "all"),
#          covars=if_else(grepl("dayOnly", type), "dayOnly", "full"))
# ggplot(auc.df, aes(model, AUC, colour=covars)) +
#   geom_point() +
#   facet_grid(species~data_subset)
# 
# auc.df %>% select(-type) %>%
#   pivot_wider(names_from="covars", values_from="AUC") %>%
#   mutate(AUC_diff=full-dayOnly) %>%
#   ggplot(aes(model, AUC_diff, colour=data_subset)) + 
#   geom_hline(yintercept=0) +
#   geom_point() + 
#   facet_wrap(~species) + theme_classic()
# 
# 
# 
# 
# 
# 
# out.b01 <- map(dir(glue("out{sep}initFit{sep}"), "bern01.*rds", full.names=T),
#                readRDS)
# p.ls <- map(out.b01, ~plot(conditional_effects(.x, surface=T), stype="raster", plot=F))
# walk2(p.ls, species, ~ggpubr::ggarrange(plotlist=.x, ncol=4, nrow=5) %>%
#         ggsave(filename=glue("figs{sep}temp{sep}bern01_{.y}.png"), plot=., width=16, height=16))
# 
# 
# out.b11 <- map(dir(glue("out{sep}initFit{sep}"), "bern11.*rds", full.names=T),
#                readRDS)
# p.ls <- map(out.b11, ~plot(conditional_effects(.x, surface=T), stype="raster", plot=F))
# walk2(p.ls, species, ~ggpubr::ggarrange(plotlist=.x, ncol=4, nrow=5) %>%
#         ggsave(filename=glue("figs{sep}temp{sep}bern11_{.y}.png"), plot=., width=16, height=16))
# 
# 
# out.ord <- map(dir(glue("out{sep}initFit{sep}"), "ord.*rds", full.names=T),
#                readRDS)
# p.ls <- map(out.ord, ~plot(conditional_effects(.x, surface=T), stype="raster", plot=F))
# walk2(p.ls, species, ~ggpubr::ggarrange(plotlist=.x, ncol=4, nrow=5) %>%
#         ggsave(filename=glue("figs{sep}temp{sep}ord_{.y}.png"), plot=., width=16, height=16))
# 
# 
# 
# 
# 
# 
# 
# 
# obs.df <- map_dfr(species, ~read_csv(glue("out{sep}dataset_{.x}.csv")) %>% 
#                     mutate(species=.x))
# 
# obs.df %>% mutate(week=week(date)) %>%
#   group_by(week, species) %>%
#   summarise(mn=mean(N.bloom),
#             q_lo=Hmisc::binconf(sum(N.bloom), n(), alpha=0.05)[1,2],
#             q_hi=Hmisc::binconf(sum(N.bloom), n(), alpha=0.05)[1,3]) %>%
#   ggplot(aes(week)) + 
#   geom_point(aes(y=mn)) +
#   geom_linerange(aes(ymin=q_lo, ymax=q_hi)) +
#   facet_wrap(~species, nrow=1) + theme_classic() +
#   scale_x_continuous(breaks=seq(1,50,length.out=5), labels=c("Jan", "Apr", "Jul", "Oct", "Dec")) +
#   ylim(0,1) +
#   labs(x="Week", y="Proportion of sites with blooms") +
#   theme(strip.text=element_text(size=12))
# ggsave(glue("figs{sep}temp{sep}obs_bloom_proportion.png"), width=12, height=3, dpi=300)
# 
# 
# obs.df %>% mutate(week=week(date)) %>%
#   filter(species %in% c("dinophysis_sp", "karenia_mikimotoi", "pseudo_nitzschia_sp")) %>%
#   group_by(week, species) %>%
#   summarise(mn=mean(N.bloom),
#             q_lo=Hmisc::binconf(sum(N.bloom), n(), alpha=0.05)[1,2],
#             q_hi=Hmisc::binconf(sum(N.bloom), n(), alpha=0.05)[1,3]) %>%
#   ggplot(aes(week)) + 
#   geom_point(aes(y=mn)) +
#   geom_linerange(aes(ymin=q_lo, ymax=q_hi)) +
#   facet_wrap(~species, nrow=1) + theme_classic() +
#   scale_x_continuous(breaks=seq(1,50,length.out=5), labels=c("Jan", "Apr", "Jul", "Oct", "Dec")) +
#   ylim(0,1) +
#   labs(x="Week", y="Proportion of sites with blooms", title="Observed bloom frequency") +
#   theme(strip.text=element_text(size=10)) 
# ggsave(glue("figs{sep}temp{sep}obs_bloom_proportion_narrow.png"), width=5, height=3.5, dpi=300)
