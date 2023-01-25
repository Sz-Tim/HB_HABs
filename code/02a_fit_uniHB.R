# HABReports Bayesian modelling
# Bayesian HAB model initial fits
# Tim Szewczyk


# This script fits models with data from 2013-2021 using horseshoe priors



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "brms")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
walk(dir("code", "*00_fn", full.names=T), source)

remake_dataset <- F

# Model details
ctrl <- list(adapt_delta=0.95, max_treedepth=20)
chains <- 4
iter <- 2000
warmup <- iter/2
refresh <- 50
prior_strength <- "full" # 1-4 or 'full'
sp <- 1

# minch2:    2013-06-20 to 2019-07-02
# WeStCOMS2: 2019-04-01 to 2022-01-26
sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")

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

if(remake_dataset) {
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
                names_glue="{.value}_{res}_{timespan}") %>%
    select(obs.id, ends_with("wk")) %>%
    left_join(., sampling.df %>% select(obs.id, site.id, date)) %>%
    filter(complete.cases(.)) %>%
    mutate(day_cumul=as.numeric(date - ymd("2010-01-01"))) %>%
    group_by(site.id) %>%
    mutate(across(ends_with("wk"),
                  ~detrend_loess(day_cumul, .x, span="adapt"), 
                  .names="{.col}_dt")) %>%
    ungroup
  connect.df <- read_csv(glue("data{sep}influx_20220801.csv")) %>%
    mutate(day_cumul=as.numeric(date - ymd("2010-01-01"))) %>%
    group_by(site.id) %>%
    mutate(influx_wk=zoo::rollsum(influx, 7, align='right', fill="extend")) %>%
    mutate(across(starts_with("influx"), ~log(.x+1))) %>%
    ungroup
  cprn.df <- read_csv("data/cprn_20211231.csv") %>%
    select(site.id, date, attn_wk, chl_wk, dino_wk, o2_wk, ph_wk, po4_wk) %>%
    filter(complete.cases(.)) %>%
    mutate(across(ends_with("wk"), log1p)) %>%
    mutate(day_cumul=as.numeric(date - ymd("2010-01-01"))) %>%
    group_by(site.id) %>%
    mutate(across(ends_with("wk"),
                  ~detrend_loess(day_cumul, .x, span=0.1), 
                  .names="{.col}_dt")) %>%
    ungroup
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
}


# >0.85 cor(res,time) = temp, salinity
# >0.85 cor(res) = short_wave, km, precip, wind, windDir
# >0.85 cor(time) = water
# waterDir is more variable across all 4
covars.all <- c(
  # WeStCOMS weekly averages
  "temp_L_wk", "salinity_L_wk", "short_wave_L_wk", "km_L_wk",
  "precip_L_wk", "tempStrat20m_L_wk", "tempStrat20m_R_wk",
  "wind_L_wk", "windDir_L_wk",
  "water_L_wk", "waterDir_L_wk", 
  "water_R_wk", "waterDir_R_wk",
  # external forcings          
  "fetch", "influx_wk",
  # Copernicus
  "attn_wk", "chl_wk", "dino_wk", "o2_wk", "ph_wk", "po4_wk",
  # detrended WeStCOMS, Copernicus (cor < 0.8 with original)          
  "temp_L_wk_dt", "salinity_L_wk_dt", "short_wave_L_wk_dt", 
  "precip_L_wk_dt", "tempStrat20m_L_wk_dt", "tempStrat20m_R_wk_dt",
  "wind_L_wk_dt", "water_L_wk_dt", "water_R_wk_dt", 
  "chl_wk_dt", "o2_wk_dt", "ph_wk_dt", "po4_wk_dt")
covar_int.all <- c(
  "ydayCos", "ydaySin",
  "tempLwk", "salinityLwk", "shortwaveLwk", "kmLwk", "precipLwk",
  "tempStrat20mLwk", "tempStrat20mRwk",
  "windVel", "waterVelL", "waterVelR", "windLwk", "waterLwk", "waterRwk",
  "fetch", "influxwk", 
  "attnwk", "chlwk", "dinowk", "o2wk", "phwk", "po4wk",
  "tempLwkdt", "salinityLwkdt", "shortwaveLwkdt", 
  "precipLwkdt", "tempStrat20mLwkdt", "tempStrat20mRwkdt",
  "windLwkdt", "waterLwkdt", "waterRwkdt", 
  "chlwkdt", "o2wkdt", "phwkdt", "po4wkdt",
  "NlnWt1", "NlnWt2",
  "NlnRAvg1", "NlnRAvg2"
)
covar_s.all <- c(
  "tempLwk", "salinityLwk", "shortwaveLwk", "kmLwk", "precipLwk",
  "tempStrat20mLwk", "tempStrat20mRwk",
  "windVel", "waterVelL", "waterVelR", "windLwk", "waterLwk", "waterRwk",
  "fetch", "influxwk", 
  "attnwk", "chlwk", "dinowk", "o2wk", "phwk", "po4wk",
  "tempLwkdt", "salinityLwkdt", "shortwaveLwkdt", 
  "precipLwkdt", "tempStrat20mLwkdt", "tempStrat20mRwkdt",
  "windLwkdt", "waterLwkdt", "waterRwkdt", 
  "chlwkdt", "o2wkdt", "phwkdt", "po4wkdt",
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
  

# priors ------------------------------------------------------------------
  
  # Smoothers
  s_b <- glue("b{covar_s}") 
  s_p <- glue("p{covar_s}")
  
  # Prior strengths
  pri.var <- switch(prior_strength,
                    "1"=list(hs1=0.5, hs2=0.6, b=0.75, de=0.3, i="1-loose"),
                    "2"=list(hs1=1, hs2=0.4, b=0.5, de=0.2, i="2-medium"),
                    "3"=list(hs1=3, hs2=0.2, b=0.2, de=0.1, i="3-tight"),
                    "4"=list(hs1=5, hs2=0.1, b=0.1, de=0.05, i="4-tighter"),
                    "full"=list(hs1=3, hs2=0.2, b=0.2, de=0.1, i="full")) 
  
  # Interaction priors
  if(i.name=="null") { 
    priors <- c(prior(normal(0, 1), class="Intercept"),
                prior(normal(0, 0.1), class="sd"))
  } else {
    priors <- c(prior_string(glue("horseshoe({pri.var$hs1}, par_ratio={pri.var$hs2})"), class="b"),
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
    prior_string(glue("double_exponential(0,{pri.var$de})"), class="sds", nlpar="bIntercept", lb=0),
    map(covar_s, 
        ~c(prior_string(glue("beta({pri.var$b},1)"), nlpar=paste0("p", .x), lb=0, ub=1),
           prior_string("normal(0,1)", class="b", nlpar=paste0("b", .x)),
           prior_string("normal(0,.5)", class="sd", nlpar=paste0("b", .x), lb=0),
           prior_string("normal(0,.5)", class="sd", nlpar=paste0("p", .x), lb=0),
           prior_string(glue("double_exponential(0,{pri.var$de})"), class="sds", nlpar=paste0("b", .x), lb=0))) %>%
      do.call('c', .)
  )
  
  
  
  
  
  

# dataset -----------------------------------------------------------------

  target <- species[sp]
  f.prefix <- glue("out{sep}{pri.var$i}{sep}")
  f.suffix <- glue("{i.name}_{target}")
  
  if(remake_dataset) {
    target.tf <- thresh.df %>% filter(hab_parameter==target)
    
    target.df <- sampling.df %>%
      rename(N=!!target) %>%
      select(obs.id, site.id, date, hour, grid, lon, lat, fetch, openBearing, N) %>%
      mutate(yday=yday(date),
             ydayCos=cos(2*pi*yday/365),
             ydaySin=sin(2*pi*yday/365),
             year=year(date),
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
      mutate(across(contains("Dir_"), ~cos(.x-openBearing))) %>% # 1 = toward open ocean, -1 = in from open ocean
      mutate(ydaySC=ydaySin*ydayCos,
             windVel=wind_L_wk*windDir_L_wk,
             waterVelL=water_L_wk*waterDir_L_wk,
             waterVelR=water_R_wk*waterDir_R_wk) %>%
      mutate(across(one_of(grep("Dir", covars, invert=T, value=T)), LaplacesDemon::CenterScale)) %>%
      rename_with(~str_remove_all(.x, "\\.|_")) %>%
      arrange(siteid, date) %>%
      select(siteid, lon, lat, date, year, obsid,
             starts_with("N"), starts_with("date_"), starts_with("yday"),
             one_of(covars, covar_int, covar_s)) %>%
      filter(complete.cases(.)) %>%
      mutate(covarSet=i.name,
             NlnRAvg1=NA,
             NlnRAvg2=NA,
             lon_sc=LaplacesDemon::CenterScale(lon),
             lat_sc=LaplacesDemon::CenterScale(lat),
             species=target) 
    
    
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
    
    write_csv(target.df, glue("{f.prefix}dataset_{i.name}_{target}.csv"))
  } else {
    target.df <- read_csv(glue("{f.prefix}dataset_{i.name}_{target}.csv"))
  }
  
  
  train.df <- target.df %>% filter(year <= 2019)
  test.df <- target.df %>% filter(year > 2019)
  bloomThresh <- max((!target.df$Nbloom)*target.df$NcatNum) # 1:4, maximum considered 'No bloom'
  
  
  

# full fits ---------------------------------------------------------------
  
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
                 file=glue("{f.prefix}ord_{f.suffix}"))
  out.ordP <- brm(form_ordP, data=train.df, 
                  family=cumulative("probit"), prior=priors.P, 
                  iter=iter, warmup=warmup, refresh=refresh, init=0,
                  control=ctrl, chains=chains, cores=chains,
                  file=glue("{f.prefix}ordP_{f.suffix}"))
  out.bern01 <- brm(form_01, data=train.df %>% filter(Nbloom1==0),
                    family=bernoulli("probit"), prior=priors, 
                    iter=iter, warmup=warmup, refresh=refresh, init=0,
                    control=ctrl, chains=chains, cores=chains,
                    file=glue("{f.prefix}bern01_{f.suffix}"))
  out.bernP01 <- brm(form_bernP, data=train.df %>% filter(Nbloom1==0), 
                     family=bernoulli("probit"), prior=priors.P, 
                     iter=iter, warmup=warmup, refresh=refresh, init=0,
                     control=ctrl, chains=chains, cores=chains,
                     file=glue("{f.prefix}bernP01_{f.suffix}"))
  out.bern11 <- brm(form_11, data=train.df %>% filter(Nbloom1==1),
                    family=bernoulli("probit"), prior=priors, 
                    iter=iter, warmup=warmup, refresh=refresh, init=0,
                    control=ctrl, chains=chains, cores=chains,
                    file=glue("{f.prefix}bern11_{f.suffix}"))
  out.bernP11 <- brm(form_bernP, data=train.df %>% filter(Nbloom1==1), 
                     family=bernoulli("probit"), prior=priors.P, 
                     iter=iter, warmup=warmup, refresh=refresh, init=0,
                     control=ctrl, chains=chains, cores=chains,
                     file=glue("{f.prefix}bernP11_{f.suffix}"))
  
  
  # Fitted
  fit.ord <- posterior_epred(out.ord)
  fit.ordP <- posterior_epred(out.ordP)
  fit.bern01 <- posterior_epred(out.bern01) 
  fit.bernP01 <- posterior_epred(out.bernP01) 
  fit.bern11 <- posterior_epred(out.bern11) 
  fit.bernP11 <- posterior_epred(out.bernP11) 
  
  fit.df <- full_join(
    train.df %>%
      mutate(ord_mnpr=calc_ord_mnpr(fit.ord, bloomThresh, summaryStat="median"),
             ordq90_mnpr=calc_ord_mnpr(fit.ord, bloomThresh, summaryStat="q90"),
             ordP_mnpr=calc_ord_mnpr(fit.ordP, bloomThresh, summaryStat="median"),
             ordPq90_mnpr=calc_ord_mnpr(fit.ordP, bloomThresh, summaryStat="q90")),
    tibble(obsid=c(filter(train.df, Nbloom1==0)$obsid, filter(train.df, Nbloom1==1)$obsid),
           bern_mnpr=c(apply(fit.bern01, 2, median), 
                       apply(fit.bern11, 2, median)),
           bernq90_mnpr=c(apply(fit.bern01, 2, quantile, probs=0.9), 
                          apply(fit.bern11, 2, quantile, probs=0.9)),
           bernP_mnpr=c(apply(fit.bernP01, 2, median), 
                        apply(fit.bernP11, 2, median)),
           bernPq90_mnpr=c(apply(fit.bernP01, 2, quantile, probs=0.9), 
                           apply(fit.bernP11, 2, quantile, probs=0.9)))
  ) %>%
    mutate(covarSet=i.name,
           species=target)
  write_csv(fit.df, glue("{f.prefix}fit_HBuv_{f.suffix}.csv"))
  
  # OOS predictions
  test.df01 <- test.df %>% filter(Nbloom1==0) %>% droplevels
  test.df11 <- test.df %>% filter(Nbloom1==1) %>% droplevels
  pred.ord <- posterior_epred(out.ord, newdata=test.df, allow_new_levels=T)
  pred.ordP <- posterior_epred(out.ordP, newdata=test.df, allow_new_levels=T)
  pred.bern01 <- posterior_epred(out.bern01, newdata=test.df01, allow_new_levels=T)
  pred.bernP01 <- posterior_epred(out.bernP01, newdata=test.df01, allow_new_levels=T)
  if(nrow(test.df11) > 0) {
    pred.bern11 <- posterior_epred(out.bern11, newdata=test.df11, allow_new_levels=T)
    pred.bernP11 <- posterior_epred(out.bernP11, newdata=test.df11, allow_new_levels=T)
    pred.bern11.md <- apply(pred.bern11, 2, median)
    pred.bernP11.md <- apply(pred.bernP11, 2, median)
    pred.bern11.q90 <- apply(pred.bern11, 2, quantile, probs=0.9)
    pred.bernP11.q90 <- apply(pred.bernP11, 2, quantile, probs=0.9)
  } else {
    pred.bern11.md <- pred.bern11.q90 <- numeric(0)
    pred.bernP11.md <- pred.bernP11.q90 <- numeric(0)
  }
  
  pred.df <- full_join(
    test.df %>%
      mutate(ord_mnpr=calc_ord_mnpr(pred.ord, bloomThresh, summaryStat="median"),
             ordq90_mnpr=calc_ord_mnpr(pred.ord, bloomThresh, summaryStat="q90"),
             ordP_mnpr=calc_ord_mnpr(pred.ordP, bloomThresh, summaryStat="median"),
             ordPq90_mnpr=calc_ord_mnpr(pred.ordP, bloomThresh, summaryStat="q90")),
    tibble(obsid=c(filter(test.df, Nbloom1==0)$obsid, filter(test.df, Nbloom1==1)$obsid),
           bern_mnpr=c(apply(pred.bern01, 2, median), pred.bern11.md),
           bernq90_mnpr=c(apply(pred.bern01, 2, quantile, probs=0.9), pred.bern11.q90),
           bernP_mnpr=c(apply(pred.bernP01, 2, median), pred.bernP11.md),
           bernPq90_mnpr=c(apply(pred.bernP01, 2, quantile, probs=0.9), pred.bernP11.q90))
  ) %>%
    mutate(covarSet=i.name,
           species=target)
  
  write_csv(pred.df, glue("{f.prefix}pred_HBuv_{f.suffix}.csv"))
  
  

# cross-validation --------------------------------------------------------

  # Cross-validation by year
  yrCV <- unique(filter(train.df, year>2014)$year)
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
                  file=glue("{f.prefix}ord_CV{k}_{f.suffix}"))
    cv.ordP <- brm(form_ordP, data=cv_train.df, 
                   family=cumulative("probit"), prior=priors.P, 
                   iter=iter, warmup=warmup, refresh=refresh, init=0,
                   control=ctrl, chains=chains, cores=chains,
                   file=glue("{f.prefix}ordP_CV{k}_{f.suffix}"))
    cv.bern01 <- brm(form_01, data=cv_train.df %>% filter(Nbloom1==0),
                     family=bernoulli("probit"), prior=priors, 
                     iter=iter, warmup=warmup, refresh=refresh, init=0,
                     control=ctrl, chains=chains, cores=chains,
                     file=glue("{f.prefix}bern01_CV{k}_{f.suffix}"))
    cv.bernP01 <- brm(form_bernP, data=cv_train.df %>% filter(Nbloom1==0), 
                      family=bernoulli("probit"), prior=priors.P, 
                      iter=iter, warmup=warmup, refresh=refresh, init=0,
                      control=ctrl, chains=chains, cores=chains,
                      file=glue("{f.prefix}bernP01_CV{k}_{f.suffix}"))
    cv.bern11 <- brm(form_11, data=cv_train.df %>% filter(Nbloom1==1),
                     family=bernoulli("probit"), prior=priors, 
                     iter=iter, warmup=warmup, refresh=refresh, init=0,
                     control=ctrl, chains=chains, cores=chains,
                     file=glue("{f.prefix}bern11_CV{k}_{f.suffix}"))
    cv.bernP11 <- brm(form_bernP, data=cv_train.df %>% filter(Nbloom1==1), 
                      family=bernoulli("probit"), prior=priors.P, 
                      iter=iter, warmup=warmup, refresh=refresh, init=0,
                      control=ctrl, chains=chains, cores=chains,
                      file=glue("{f.prefix}bernP11_CV{k}_{f.suffix}"))
    
    # Cross-validation predictions
    cv_test.df01 <- cv_test.df %>% filter(Nbloom1==0) %>% droplevels
    cv_test.df11 <- cv_test.df %>% filter(Nbloom1==1) %>% droplevels
    pred.ord <- posterior_epred(cv.ord, newdata=cv_test.df, allow_new_levels=T)
    pred.ordP <- posterior_epred(cv.ordP, newdata=cv_test.df, allow_new_levels=T)
    pred.bern01 <- posterior_epred(cv.bern01, newdata=cv_test.df01, allow_new_levels=T)
    pred.bernP01 <- posterior_epred(cv.bernP01, newdata=cv_test.df01, allow_new_levels=T)
    if(nrow(cv_test.df11) > 0) {
      pred.bern11 <- posterior_epred(cv.bern11, newdata=cv_test.df11, allow_new_levels=T)
      pred.bernP11 <- posterior_epred(cv.bernP11, newdata=cv_test.df11, allow_new_levels=T)
      pred.bern11.md <- apply(pred.bern11, 2, median)
      pred.bernP11.md <- apply(pred.bernP11, 2, median)
      pred.bern11.q90 <- apply(pred.bern11, 2, quantile, probs=0.9)
      pred.bernP11.q90 <- apply(pred.bernP11, 2, quantile, probs=0.9)
    } else {
      pred.bern11.md <- pred.bern11.q90 <- numeric(0)
      pred.bernP11.md <- pred.bernP11.q90 <- numeric(0)
    }
    
    cv_pred[[k]] <- full_join(
      cv_test.df %>%
        mutate(ord_mnpr=calc_ord_mnpr(pred.ord, bloomThresh, summaryStat="median"),
               ordq90_mnpr=calc_ord_mnpr(pred.ord, bloomThresh, summaryStat="q90"),
               ordP_mnpr=calc_ord_mnpr(pred.ordP, bloomThresh, summaryStat="median"),
               ordPq90_mnpr=calc_ord_mnpr(pred.ordP, bloomThresh, summaryStat="q90")),
      tibble(obsid=c(filter(test.df, Nbloom1==0)$obsid, filter(test.df, Nbloom1==1)$obsid),
             bern_mnpr=c(apply(pred.bern01, 2, median), pred.bern11.md),
             bernq90_mnpr=c(apply(pred.bern01, 2, quantile, probs=0.9), pred.bern11.q90),
             bernP_mnpr=c(apply(pred.bernP01, 2, median), pred.bernP11.md),
             bernPq90_mnpr=c(apply(pred.bernP01, 2, quantile, probs=0.9), pred.bernP11.q90))
    ) %>%
      mutate(covarSet=i.name,
             species=target)
  }
  cv_pred %>% do.call('rbind', .) %>%
    write_csv(glue("{f.prefix}CV_HBuv_{f.suffix}.csv"))
  
  cat("Finished", target, "\n")
}





