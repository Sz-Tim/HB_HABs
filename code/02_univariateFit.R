# HABReports Bayesian modelling
# Bayesian HAB model initial fits
# Tim Szewczyk


# This script fits models with data from 2013-2017 using horseshoe priors



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "brms", 
          "randomForest", "caret", "e1071", "xgboost")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
walk(dir("code", "*00_fn", full.names=T), source)

# Model details
ctrl <- list(adapt_delta=0.95, max_treedepth=20)
chains <- 4
iter <- 2000
warmup <- iter/2
refresh <- 50

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

# tribble(~full, ~abbr, ~clean,
#         "alexandrium_sp", "Alsp", "Alexandrium sp.",
#         "dinophysis_sp", "Disp", "Dinophysis sp.",
#         "karenia_mikimotoi", "Kami", "Karenia mikimotoi",
#         "prorocentrum_lima", "Prli", "Prorocentrum lima",
#         "pseudo_nitzschia_sp", "Pssp", "Pseudo-nitzschia sp.") %>%
#   write_csv("data/sp_i.csv")
sp.i <- read_csv("data/sp_i.csv")




# covariates --------------------------------------------------------------

species <- sp.i$full

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
         waterDir=atan2(v, u), # moving N = 0; moving S = pi
         windDir=atan2(vwind_speed, uwind_speed), # blowing N = 0; blowing S = pi
         timespan=if_else(lag=="0", "0", "wk")) %>% 
  select(-u, -v, -ww, -vwind_speed, -uwind_speed) %>%
  group_by(obs.id, res, timespan) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm=T))) %>%
  ungroup %>%
  pivot_wider(names_from=c(res, timespan), values_from=4:ncol(.),
              names_glue="{.value}_{res}_{timespan}")
# TODO: this is only for WC 1 -- update for > 2019; pretty sure I have this (more or less) on salmon HAB_operational/ptrack/out/init_*
connect.df <- read_csv(glue("data{sep}connectivity_5e3parts.csv")) %>%
  group_by(site.id) %>%
  mutate(influx_wk=zoo::rollsum(influx, 7, align='right', fill="extend")) %>%
  mutate(across(starts_with("influx"), ~log(.x+1)))
cprn.df <- read_csv("data/cprn.csv") %>%
  select(site.id, date, attn_wk, chl_wk, dino_wk, o2_wk, ph_wk, po4_wk)
site.100k <-  read_csv("data/fsa_site_pairwise_distances.csv") %>% 
  bind_rows(tibble(origin=unique(.$origin), 
                   destination=unique(.$origin), 
                   distance=0)) %>%
  filter(distance < 100e3) %>% 
  select(-distance) %>% 
  group_by(origin) %>% 
  nest(data=destination) %>%
  mutate(destination_c=c(data[[1]])) %>% 
  select(-data) %>%
  ungroup

# >0.85 cor(res,time) = temp, salinity
# >0.85 cor(res) = short_wave, km, precip, wind, windDir
# >0.85 cor(time) = water
# waterDir is more variable across all 4
covars.all <- c("temp_L_wk", "salinity_L_wk", "short_wave_L_wk", "km_L_wk",
            "precip_L_wk", "tempStrat20m_L_wk", "tempStrat20m_R_wk",
            "wind_L_wk", "windDir_L_wk",
            "water_L_wk", "waterDir_L_wk", 
            "water_R_wk", "waterDir_R_wk",
            "fetch", #"influxwk",  
            "attn_wk", "chl_wk", "dino_wk", "o2_wk", "ph_wk", "po4_wk")
covar_int.all <- c(
  "ydayCos", "ydaySin",
  paste0(c("tempLwk", "salinityLwk", "shortwaveLwk", "kmLwk", "precipLwk",
           "tempStrat20mLwk", "tempStrat20mRwk",
           "windVel", "waterVelL", "waterVelR", "windLwk", "waterLwk", "waterRwk",
           "fetch", #"influxwk", 
           "attnwk", "chlwk", "dinowk", "o2wk", "phwk", "po4wk",
           # "mo(NcatF1)", "mo(NcatF2)",
           "Nbloom1", "Nbloom2",
           "NlnWt1", "NlnWt2",
           "NlnRAvg1", "NlnRAvg2"), 
         ":ydayCos:ydaySin"),
  # "mo(NcatF1):mo(NcatF2)",
  paste0(c("windVel", "waterVelL", "waterVelR", "windLwk", "waterLwk", "waterRwk"), 
         ":fetch")
)
covar_s.all <- c(
  "tempLwk", "salinityLwk", "shortwaveLwk", "kmLwk", "precipLwk",
  "tempStrat20mLwk", "tempStrat20mRwk",
  "windVel", "waterVelL", "waterVelR", "windLwk", "waterLwk", "waterRwk",
  "fetch", #"influxwk", 
  "attnwk", "chlwk", "dinowk", "o2wk", "phwk", "po4wk",
  "Nbloom1", "Nbloom2",
  "NlnWt1", "NlnWt2",
  "NlnRAvg1", "NlnRAvg2"
)


covariate_sets <- list(
  # null="NA",
  # date="^yday",
  # autoreg="^N|mo",
  # external="fetch|influx|water|wind",
  # local="temp|salinity|shortwave|precip|km",
  # cprn="attn|chl|dino|o2|ph|po4",
  all="."
)



for(i in length(covariate_sets)) {
  i.name <- names(covariate_sets)[i]
  i.covs <- covariate_sets[[i]]
  covars <- covars.all[grepl(i.covs, str_remove_all(covars.all, "\\.|_"))]
  covar_int <- grep(i.covs, covar_int.all, value=T)
  covar_s <- grep(i.covs, covar_s.all, value=T)
  covar_date <- NULL
  
  # Smoothers
  s_b <- glue("b{covar_s}") 
  s_p <- glue("p{covar_s}")
  
  # Interaction priors
  if(i.name=="null") { 
    priors <- c(prior(normal(0, 1), class="Intercept"),
                prior(normal(0, 0.1), class="sd"))
  } else {
    priors <- c(prior(horseshoe(3, par_ratio=0.2), class="b"),
                prior(normal(0, 1), class="Intercept"),
                prior(normal(0, 0.1), class="sd"))
  }
  
  # Indicator variable: p * b * X
  flist.P <- c(
    map(s_b, ~as.formula(paste0(.x, "~ s(ydayCos,ydaySin) + (1|siteid)"))),
    map(s_p, ~as.formula(paste0(.x, "~ 1 + (1|siteid)"))),
    map(1, ~"bIntercept ~ 1 + s(ydayCos,ydaySin) + (1|siteid)")
  )
  priors.P <- c(
    prior_string("normal(0,1)", class="b", nlpar="bIntercept"),
    prior_string("normal(0,.5)", class="sd", nlpar="bIntercept", lb=0),
    prior_string("double_exponential(0,0.1)", class="sds", nlpar="bIntercept", lb=0),
    map(covar_s, 
        ~c(prior_string("beta(0.2,1)", nlpar=paste0("p", .x), lb=0, ub=1),
           prior_string("normal(0,1)", class="b", nlpar=paste0("b", .x)),
           prior_string("normal(0,.5)", class="sd", nlpar=paste0("b", .x), lb=0),
           prior_string("normal(0,.5)", class="sd", nlpar=paste0("p", .x), lb=0),
           prior_string("double_exponential(0,0.1)", class="sds", nlpar=paste0("b", .x), lb=0))) %>%
      do.call('c', .)
  )
  
  
  
  
  
  
  # initial fit -------------------------------------------------------------
  
  for(sp in 1) {
    target <- species[sp]
    target.tf <- thresh.df %>% filter(hab_parameter==target)
    
    target.df <- sampling.df %>%
      rename(N=!!target) %>%
      select(obs.id, site.id, date, hour, grid, lon, lat, fetch, bearing, N) %>%
      mutate(yday=yday(date),
             ydayCos=cos(2*pi*yday/365),
             ydaySin=sin(2*pi*yday/365),
             year=year(date),
             bearing=bearing*pi/180, # Nearest shore: N=pi, S=0
             N=round(N),
             N.ln=log(N+1)) %>%
      rowwise() %>%
      mutate(N.cat=target.tf$tl[max(which(N >= target.tf$min_ge))]) %>%
      ungroup %>%
      mutate(N.catF=factor(N.cat, levels=unique(target.tf$tl), ordered=T),
             N.catNum=as.numeric(N.catF),
             N.bloom=target.tf$bloom[match(N.cat, target.tf$tl)]) %>%
      arrange(site.id, date) %>%
      group_by(site.id) %>%
      multijetlag(N.ln, N.cat, N.catF, N.bloom, date, n=2) %>%
      ungroup %>%
      mutate(across(starts_with("date_"), ~as.numeric(date-.x)),
             # I don't love this since small if N.ln_x is small OR date_x is large
             N.lnWt_1=N.ln_1/date_1,
             N.lnWt_2=N.ln_2/date_2) %>%
      full_join(hydro.df) %>%
      full_join(connect.df) %>%
      full_join(cprn.df) %>%
      mutate(across(contains("Dir_"), ~cos(.x-bearing))) %>% # -1 = toward shore, 1 = away from shore
      mutate(across(one_of(grep("Dir", covars, invert=T, value=T)), LaplacesDemon::CenterScale)) %>%
      arrange(site.id, date) %>%
      rename_with(~str_remove_all(.x, "\\.|_")) %>%
      mutate(ydaySC=ydaySin*ydayCos,
             windVel=windLwk*windDirLwk,
             waterVelL=waterLwk*waterDirLwk,
             waterVelR=waterRwk*waterDirRwk) %>%
      select(siteid, lon, lat, date, year, obsid,
             starts_with("N"), starts_with("date_"), starts_with("yday"),
             one_of(covars, covar_int, covar_s)) %>% #, covar_date)) %>%
      filter(complete.cases(.)) %>%
      mutate(covarSet=i.name,
             NlnRAvg1=NA,
             NlnRAvg2=NA,
             lon_sc=LaplacesDemon::CenterScale(lon),
             lat_sc=LaplacesDemon::CenterScale(lat)) %>%
      filter(!siteid %in% c(70, 74, 75, 80, 88))
    # TODO: I don't remember why these sites are excluded
    
    if("NlnRAvg1" %in% covar_s) {
      for(j in 1:nrow(target.df)) {
        site_j <- target.df$siteid[j]
        date_j <- target.df$date[j]
        target_wk <- target.df %>% select(siteid, date, Nln1, Nln2) %>%
          filter(siteid %in% site.100k$destination_c[site.100k$origin==site_j][[1]]) %>%
          filter(date <= date_j & date > date_j-7) 
        target.df$NlnRAvg1[j] <- mean(target_wk$Nln1)
        target.df$NlnRAvg2[j] <- mean(target_wk$Nln2)
        if(j %% 100 == 0) {cat(j, "of", nrow(target.df), "\n")}
      } 
    }
    
    write_csv(target.df, glue("out{sep}full{sep}dataset_{i.name}_{target}.csv"))
    
    train.df <- target.df %>% filter(year <= 2019)
    test.df <- target.df %>% filter(year > 2019)
    bloomThresh <- max((!target.df$Nbloom)*target.df$NcatNum) # 1:4, maximum considered 'No bloom'
    
    # Full fit for predictions
    
    # Formulas with interactions: errors if missing NcatF levels
    form_ord <- makeFormula(train.df, covar_int, "NcatNum | thres(3)")
    form_01 <- makeFormula(filter(train.df, Nbloom1==0), covar_int, "Nbloom")
    form_11 <- makeFormula(filter(train.df, Nbloom1==1), covar_int, "Nbloom")
    
    # Smoother formulas
    form_ordP <- makeFormula(train.df, covar_s, "NcatNum | thres(3)", covar_date,
                             flist.P, list(b=s_b, p=s_p))
    form_bernP <- makeFormula(train.df, covar_s, "Nbloom", covar_date,
                              flist.P, list(b=s_b, p=s_p))
    
    out.ord <- brm(form_ord, data=train.df,
                   family=cumulative("probit"), prior=priors, 
                   iter=iter, warmup=warmup, refresh=refresh, init=0,
                   control=ctrl, chains=chains, cores=chains,
                   file=glue("out{sep}full{sep}ord_{i.name}_{target}"))
    out.ordP <- brm(form_ordP, data=train.df, 
                    family=cumulative("probit"), prior=priors.P, 
                    iter=iter, warmup=warmup, refresh=refresh, init=0,
                    control=ctrl, chains=chains, cores=chains,
                    file=glue("out{sep}full{sep}ordP_{i.name}_{target}"))
    out.bern01 <- brm(form_01, data=train.df %>% filter(Nbloom1==0),
                      family=bernoulli("probit"), prior=priors, 
                      iter=iter, warmup=warmup, refresh=refresh, init=0,
                      control=ctrl, chains=chains, cores=chains,
                      file=glue("out{sep}full{sep}bern01_{i.name}_{target}"))
    out.bernP01 <- brm(form_bernP, data=train.df %>% filter(Nbloom1==0), 
                       family=bernoulli("probit"), prior=priors.P, 
                       iter=iter, warmup=warmup, refresh=refresh, init=0,
                       control=ctrl, chains=chains, cores=chains,
                       file=glue("out{sep}full{sep}bernP01_{i.name}_{target}"))
    out.bern11 <- brm(form_11, data=train.df %>% filter(Nbloom1==1),
                      family=bernoulli("probit"), prior=priors, 
                      iter=iter, warmup=warmup, refresh=refresh, init=0,
                      control=ctrl, chains=chains, cores=chains,
                      file=glue("out{sep}full{sep}bern11_{i.name}_{target}"))
    out.bernP11 <- brm(form_bernP, data=train.df %>% filter(Nbloom1==1), 
                       family=bernoulli("probit"), prior=priors.P, 
                       iter=iter, warmup=warmup, refresh=refresh, init=0,
                       control=ctrl, chains=chains, cores=chains,
                       file=glue("out{sep}full{sep}bernP11_{i.name}_{target}"))
    
    # Machine Learning: Random Forest, Support Vector Machine, XGBoost
    ML_vars <- c("Nbloom", "Nbloom1", "lon_sc", "lat_sc", covar_date, covar_s)
    train.ML <- train.df %>% select(one_of(ML_vars)) %>% 
      mutate(Nbloom=factor(Nbloom)) %>% as.data.frame()
    train.ML01 <- train.df %>% filter(Nbloom1==0) %>% select(one_of(ML_vars)) %>% 
      mutate(Nbloom=factor(Nbloom)) %>% as.data.frame()
    train.ML11 <- train.df %>% filter(Nbloom1==1) %>% select(one_of(ML_vars)) %>% 
      mutate(Nbloom=factor(Nbloom)) %>% as.data.frame()
    test.ML <- test.df %>% select(one_of(ML_vars)) %>% 
      mutate(Nbloom=factor(Nbloom)) %>% as.data.frame()
    test.ML01 <- test.df %>% filter(Nbloom1==0) %>% select(one_of(ML_vars)) %>% 
      mutate(Nbloom=factor(Nbloom)) %>% as.data.frame()
    test.ML11 <- test.df %>% filter(Nbloom1==1) %>% select(one_of(ML_vars)) %>% 
      mutate(Nbloom=factor(Nbloom)) %>% as.data.frame()
    
    rf <- tuneRF(x=train.ML[,-1], y=train.ML[,1], doBest=T, trace=F, plot=F)
    rf.01 <- tuneRF(x=train.ML01[,-1], y=train.ML01[,1], doBest=T, trace=F, plot=F)
    rf.11 <- tuneRF(x=train.ML11[,-1], y=train.ML11[,1], doBest=T, trace=F, plot=F)
    saveRDS(rf, glue("out{sep}full{sep}rf_{i.name}_{target}.rds"))
    saveRDS(rf.01, glue("out{sep}full{sep}rf01_{i.name}_{target}.rds"))
    saveRDS(rf.11, glue("out{sep}full{sep}rf11_{i.name}_{target}.rds"))
    
    svm_rng <- list(epsilon=seq(0,0.5,0.05), cost=2^(seq(3,7,0.2)))
    svm_ <- tune(svm, train.ML[,-1], train.ML[,1], probability=T, ranges=svm_rng)
    svm_01 <- tune(svm, train.ML01[,-1], train.ML01[,1], probability=T, ranges=svm_rng)
    svm_11 <- tune(svm, train.ML11[,-1], train.ML11[,1], probability=T, ranges=svm_rng)
    saveRDS(svm_, glue("out{sep}full{sep}svm_{i.name}_{target}.rds"))
    saveRDS(svm_01, glue("out{sep}full{sep}svm01_{i.name}_{target}.rds"))
    saveRDS(svm_11, glue("out{sep}full{sep}svm11_{i.name}_{target}.rds"))
    
    xg <- xgboost(data=as.matrix(train.ML[,-1]), label=as.numeric(train.ML[,1])-1, 
                  max.depth=2, eta=1, nthread=6, nrounds=10, objective="binary:logistic")
    xg.01 <- xgboost(data=as.matrix(train.ML01[,-1]), label=as.numeric(train.ML01[,1])-1, 
                     max.depth=2, eta=1, nthread=6, nrounds=10, objective="binary:logistic")
    xg.11 <- xgboost(data=as.matrix(train.ML11[,-1]), label=as.numeric(train.ML11[,1])-1, 
                     max.depth=2, eta=1, nthread=6, nrounds=10, objective="binary:logistic")
    saveRDS(xg, glue("out{sep}full{sep}xgb_{i.name}_{target}.rds"))
    saveRDS(xg.01, glue("out{sep}full{sep}xgb01_{i.name}_{target}.rds"))
    saveRDS(xg.11, glue("out{sep}full{sep}xgb11_{i.name}_{target}.rds"))
    
    # Fitted
    fit.ord <- posterior_epred(out.ord)
    fit.ordP <- posterior_epred(out.ordP)
    fit.bern01 <- posterior_epred(out.bern01) 
    fit.bernP01 <- posterior_epred(out.bernP01) 
    fit.bern11 <- posterior_epred(out.bern11) 
    fit.bernP11 <- posterior_epred(out.bernP11) 
    
    fit.df <- full_join(
      train.df %>%
        mutate(ord_mnpr=calc_ord_mnpr(fit.ord, bloomThresh),
               ordP_mnpr=calc_ord_mnpr(fit.ordP, bloomThresh),
               rf_mnpr=rf$votes[,2],
               svm_mnpr=attr(predict(svm_$best.model, newdata=train.ML[,-1], probability=T),
                             "probabilities")[,2],
               xgb_mnpr=predict(xg, as.matrix(train.ML[,-1]))),
      tibble(obsid=c(filter(train.df, Nbloom1==0)$obsid, filter(train.df, Nbloom1==1)$obsid),
             bern_mnpr=c(colMeans(fit.bern01), colMeans(fit.bern11)),
             bernP_mnpr=c(colMeans(fit.bernP01), colMeans(fit.bernP11)),
             rf_split_mnpr=c(rf.01$votes[,2], rf.11$votes[,2]),
             svm_split_mnpr=c(attr(predict(svm_01$best.model, newdata=train.ML01[,-1], probability=T),
                                   "probabilities")[,2],
                              attr(predict(svm_11$best.model, newdata=train.ML11[,-1], probability=T),
                                   "probabilities")[,2]),
             xgb_split_mnpr=c(predict(xg.01, as.matrix(train.ML01[,-1])),
                              predict(xg.11, as.matrix(train.ML11[,-1]))))
    ) %>%
      mutate(covarSet=i.name)
    write_csv(fit.df, glue("out{sep}full{sep}fit_{i.name}_{target}.csv"))
    
    # OOS predictions
    test.df01 <- test.df %>% filter(Nbloom1==0) %>% droplevels
    test.df11 <- test.df %>% filter(Nbloom1==1) %>% droplevels
    pred.ord <- posterior_epred(out.ord, newdata=test.df, allow_new_levels=T)
    pred.ordP <- posterior_epred(out.ordP, newdata=test.df, allow_new_levels=T)
    pred.bern01 <- colMeans(posterior_epred(out.bern01, newdata=test.df01, allow_new_levels=T))
    pred.bernP01 <- colMeans(posterior_epred(out.bernP01, newdata=test.df01, allow_new_levels=T))
    pred.rf <- predict(rf, newdata=test.ML, type="prob")[,2]
    pred.rf_01 <- predict(rf.01, newdata=test.ML01, type="prob")[,2]
    pred.svm <- attr(predict(svm_$best.model, newdata=test.ML[,-1], probability=T),
                     "probabilities")[,2]
    pred.svm_01 <- attr(predict(svm_01$best.model, newdata=test.ML01[,-1], probability=T),
                        "probabilities")[,2]
    pred.xgb <- predict(xg, as.matrix(test.ML[,-1]))
    pred.xgb_01 <- predict(xg.01, as.matrix(test.ML01[,-1]))
    if(nrow(test.df11) > 0) {
      pred.bern11 <- colMeans(posterior_epred(out.bern11, newdata=test.df11, allow_new_levels=T))
      pred.bernP11 <- colMeans(posterior_epred(out.bernP11, newdata=test.df11, allow_new_levels=T))
      pred.rf_11 <- predict(rf.11, newdata=test.ML11, type="prob")[,2]
      pred.svm_11 <- attr(predict(svm_11$best.model, newdata=test.ML11[,-1], probability=T),
                          "probabilities")[,2]
      pred.xgb_11 <- predict(xg.11, as.matrix(test.ML11[,-1]))
    } else {
      pred.bern11 <- numeric(0)
      pred.bernP11 <- numeric(0)
      pred.rf_11 <- numeric(0)
      pred.svm_11 <- numeric(0)
      pred.xgb_11 <- numeric(0)
    }
    
    pred.df <- full_join(
      test.df %>%
        mutate(ord_mnpr=calc_ord_mnpr(pred.ord, bloomThresh),
               ordP_mnpr=calc_ord_mnpr(pred.ordP, bloomThresh),
               rf_mnpr=pred.rf,
               svm_mnpr=pred.svm,
               xgb_mnpr=pred.xgb),
      tibble(obsid=c(filter(test.df, Nbloom1==0)$obsid, filter(test.df, Nbloom1==1)$obsid),
             bern_mnpr=c(pred.bern01, pred.bern11),
             bernP_mnpr=c(pred.bernP01, pred.bernP11),
             rf_split_mnpr=c(pred.rf_01, pred.rf_11),
             svm_split_mnpr=c(pred.svm_01, pred.svm_11),
             xgb_split_mnpr=c(pred.xgb_01, pred.xgb_11))
    ) %>%
      mutate(covarSet=i.name)
    
    write_csv(pred.df, glue("out{sep}full{sep}pred_{i.name}_{target}.csv"))
    
    # Cross-validation by year
    yrCV <- unique(train.df$year)
    cv_pred <- map(yrCV, ~NULL)
    for(k in 1:length(yrCV)) {
      yr <- yrCV[k]
      cv_train.df <- train.df %>% filter(year != yr)
      cv_test.df <- train.df %>% filter(year == yr)
      
      # Formulas with interactions: errors if missing NcatF levels
      form_ord <- makeFormula(cv_train.df, covar_int, "NcatNum | thres(3)")
      form_01 <- makeFormula(filter(cv_train.df, Nbloom1==0), covar_int, "Nbloom")
      form_11 <- makeFormula(filter(cv_train.df, Nbloom1==1), covar_int, "Nbloom")
      
      # Smoother formulas
      form_ordP <- makeFormula(cv_train.df, covar_s, "NcatNum | thres(3)", covar_date,
                               flist.P, list(b=s_b, p=s_p))
      form_bernP <- makeFormula(cv_train.df, covar_s, "Nbloom", covar_date,
                                flist.P, list(b=s_b, p=s_p))
      
      cv.ord <- brm(form_ord, data=cv_train.df,
                    family=cumulative("probit"), prior=priors, 
                    iter=iter, warmup=warmup, refresh=refresh, init=0,
                    control=ctrl, chains=chains, cores=chains,
                    file=glue("out{sep}full{sep}ord_CV{k}_{i.name}_{target}"))
      cv.ordP <- brm(form_ordP, data=cv_train.df, 
                     family=cumulative("probit"), prior=priors.P, 
                     iter=iter, warmup=warmup, refresh=refresh, init=0,
                     control=ctrl, chains=chains, cores=chains,
                     file=glue("out{sep}full{sep}ordP_CV{k}_{i.name}_{target}"))
      cv.bern01 <- brm(form_01, data=cv_train.df %>% filter(Nbloom1==0),
                       family=bernoulli("probit"), prior=priors, 
                       iter=iter, warmup=warmup, refresh=refresh, init=0,
                       control=ctrl, chains=chains, cores=chains,
                       file=glue("out{sep}full{sep}bern01_CV{k}_{i.name}_{target}"))
      cv.bernP01 <- brm(form_bernP, data=cv_train.df %>% filter(Nbloom1==0), 
                        family=bernoulli("probit"), prior=priors.P, 
                        iter=iter, warmup=warmup, refresh=refresh, init=0,
                        control=ctrl, chains=chains, cores=chains,
                        file=glue("out{sep}full{sep}bernP01_CV{k}_{i.name}_{target}"))
      cv.bern11 <- brm(form_11, data=cv_train.df %>% filter(Nbloom1==1),
                       family=bernoulli("probit"), prior=priors, 
                       iter=iter, warmup=warmup, refresh=refresh, init=0,
                       control=ctrl, chains=chains, cores=chains,
                       file=glue("out{sep}full{sep}bern11_CV{k}_{i.name}_{target}"))
      cv.bernP11 <- brm(form_bernP, data=cv_train.df %>% filter(Nbloom1==1), 
                        family=bernoulli("probit"), prior=priors.P, 
                        iter=iter, warmup=warmup, refresh=refresh, init=0,
                        control=ctrl, chains=chains, cores=chains,
                        file=glue("out{sep}full{sep}bernP11_CV{k}_{i.name}_{target}"))
      
      # Machine Learning: Random Forest, Support Vector Machine
      ML_vars <- c("Nbloom", "Nbloom1", "lon_sc", "lat_sc", covar_date, covar_s)
      train.ML <- cv_train.df %>% select(one_of(ML_vars)) %>% 
        mutate(Nbloom=factor(Nbloom)) %>% as.data.frame()
      train.ML01 <- cv_train.df %>% filter(Nbloom1==0) %>% select(one_of(ML_vars)) %>% 
        mutate(Nbloom=factor(Nbloom)) %>% as.data.frame()
      train.ML11 <- cv_train.df %>% filter(Nbloom1==1) %>% select(one_of(ML_vars)) %>% 
        mutate(Nbloom=factor(Nbloom)) %>% as.data.frame()
      test.ML <- cv_test.df %>% select(one_of(ML_vars)) %>% 
        mutate(Nbloom=factor(Nbloom)) %>%  as.data.frame()
      test.ML01 <- cv_test.df %>% filter(Nbloom1==0) %>% select(one_of(ML_vars)) %>% 
        mutate(Nbloom=factor(Nbloom)) %>% as.data.frame()
      test.ML11 <- cv_test.df %>% filter(Nbloom1==1) %>% select(one_of(ML_vars)) %>% 
        mutate(Nbloom=factor(Nbloom)) %>% as.data.frame()
      
      rf <- tuneRF(x=train.ML[,-1], y=train.ML[,1], doBest=T, trace=F, plot=F, ntree=1000)
      rf.01 <- tuneRF(x=train.ML01[,-1], y=train.ML01[,1], doBest=T, trace=F, plot=F, ntree=1000)
      rf.11 <- tuneRF(x=train.ML11[,-1], y=train.ML11[,1], doBest=T, trace=F, plot=F, ntree=1000)
      
      svm_rng <- list(epsilon=seq(0,0.5,0.05), cost=2^(seq(3,7,0.2)))
      svm_ <- tune(svm, train.ML[,-1], train.ML[,1], probability=T, ranges=svm_rng)
      svm_01 <- tune(svm, train.ML01[,-1], train.ML01[,1], probability=T, ranges=svm_rng)
      svm_11 <- tune(svm, train.ML11[,-1], train.ML11[,1], probability=T, ranges=svm_rng)
      
      xg <- xgboost(data=as.matrix(train.ML[,-1]), label=as.numeric(train.ML[,1])-1, 
                    max.depth=2, eta=1, nthread=6, nrounds=10, objective="binary:logistic")
      xg.01 <- xgboost(data=as.matrix(train.ML01[,-1]), label=as.numeric(train.ML01[,1])-1, 
                       max.depth=2, eta=1, nthread=6, nrounds=10, objective="binary:logistic")
      xg.11 <- xgboost(data=as.matrix(train.ML11[,-1]), label=as.numeric(train.ML11[,1])-1, 
                       max.depth=2, eta=1, nthread=6, nrounds=10, objective="binary:logistic")
      
      
      # Cross-validation predictions
      cv_test.df01 <- cv_test.df %>% filter(Nbloom1==0) %>% droplevels
      cv_test.df11 <- cv_test.df %>% filter(Nbloom1==1) %>% droplevels
      pred.ord <- posterior_epred(cv.ord, newdata=cv_test.df, allow_new_levels=T)
      pred.ordP <- posterior_epred(cv.ordP, newdata=cv_test.df, allow_new_levels=T)
      pred.bern01 <- colMeans(posterior_epred(cv.bern01, newdata=cv_test.df01, allow_new_levels=T))
      pred.bernP01 <- colMeans(posterior_epred(cv.bernP01, newdata=cv_test.df01, allow_new_levels=T))
      pred.rf <- predict(rf, newdata=test.ML, type="prob")[,2]
      pred.rf_01 <- predict(rf.01, newdata=test.ML01, type="prob")[,2]
      pred.svm <- attr(predict(svm_$best.model, newdata=test.ML[,-1], probability=T),
                       "probabilities")[,2]
      pred.svm_01 <- attr(predict(svm_01$best.model, newdata=test.ML01[,-1], probability=T),
                          "probabilities")[,2]
      pred.xgb <- predict(xg, as.matrix(test.ML[,-1]))
      pred.xgb_01 <- predict(xg.01, as.matrix(test.ML01[,-1]))
      if(nrow(cv_test.df11) > 0) {
        pred.bern11 <- colMeans(posterior_epred(cv.bern11, newdata=cv_test.df11, allow_new_levels=T))
        pred.bernP11 <- colMeans(posterior_epred(cv.bernP11, newdata=cv_test.df11, allow_new_levels=T))
        pred.rf_11 <- predict(rf.11, newdata=test.ML11, type="prob")[,2]
        pred.svm_11 <- attr(predict(svm_11$best.model, newdata=test.ML11[,-1], probability=T),
                            "probabilities")[,2]
        pred.xgb_11 <- predict(xg.11, as.matrix(test.ML11[,-1]))
      } else {
        pred.bern11 <- numeric(0)
        pred.bernP11 <- numeric(0)
        pred.rf_11 <- numeric(0)
        pred.svm_11 <- numeric(0)
        pred.xgb_11 <- numeric(0)
      }
      
      cv_pred[[k]] <- full_join(
        cv_test.df %>%
          mutate(ord_mnpr=calc_ord_mnpr(pred.ord, bloomThresh),
                 ordP_mnpr=calc_ord_mnpr(pred.ordP, bloomThresh),
                 rf_mnpr=pred.rf,
                 svm_mnpr=pred.svm,
                 xgb_mnpr=pred.xgb),
        tibble(obsid=c(filter(cv_test.df, Nbloom1==0)$obsid, filter(cv_test.df, Nbloom1==1)$obsid),
               bern_mnpr=c(pred.bern01, pred.bern11),
               bernP_mnpr=c(pred.bernP01, pred.bernP11),
               rf_split_mnpr=c(pred.rf_01, pred.rf_11),
               svm_split_mnpr=c(pred.svm_01, pred.svm_11),
               xgb_split_mnpr=c(pred.xgb_01, pred.xgb_11))
      ) %>%
        mutate(covarSet=i.name)
    }
    cv_pred %>% do.call('rbind', .) %>%
      write_csv(glue("out{sep}full{sep}CV_{i.name}_{target}.csv"))
    
    cat("Finished", target, "\n")
  }
  
}




