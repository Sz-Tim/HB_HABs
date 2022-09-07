# HABReports Bayesian modelling
# Bayesian HAB model initial fits
# Tim Szewczyk


# This script fits models with data from 2013-2017 using horseshoe priors



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "brms", 
          "randomForest", "caret")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
walk(dir("code", "*00_fn", full.names=T), source)

# Model details
ctrl <- list(adapt_delta=0.95, max_treedepth=20)
chains <- 4
iter <- 2000
warmup <- iter/2
refresh <- 500

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
         waterDir=atan2(v, u), # moving N = 0; moving S = pi
         windDir=atan2(vwind_speed, uwind_speed), # blowing N = 0; blowing S = pi
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
            "precip_L_wk", 
            "wind_L_wk", "windDir_L_wk",
            "water_L_wk", "waterDir_L_wk", 
            "water_R_wk", "waterDir_R_wk",
            "influx_wk", 
            "attn_wk", "chl_wk", "dino_wk", "o2_wk", "ph_wk", "po4_wk")
covar_int.all <- c(
  "ydayCos", "ydaySin",
  paste0(c("tempLwk", "salinityLwk", "shortwaveLwk", "kmLwk", "precipLwk",
           "windVel", "waterVelL", "waterVelR",
           "influxwk", "fetch",
           "attnwk", "chlwk", "dinowk", "o2wk", "phwk", "po4wk",
           # "mo(NcatF1)", "mo(NcatF2)",
           "Nbloom1", "Nbloom2",
           "NlnWt1", "NlnWt2",
           "NlnRAvg1", "NlnRAvg2"), 
         ":ydayCos:ydaySin"),
  # "mo(NcatF1):mo(NcatF2)",
  paste0(c("windVel", "waterVelL", "waterVelR"), 
         ":fetch")
)
covar_s.all <- c(
  "tempLwk", "salinityLwk", "shortwaveLwk", "kmLwk", "precipLwk",
  "windVel", "waterVelL", "waterVelR",
  "influxwk", "fetch",
  "attnwk", "chlwk", "dinowk", "o2wk", "phwk", "po4wk",
  "Nbloom1", "Nbloom2",
  "NlnWt1", "NlnWt2",
  "NlnRAvg1", "NlnRAvg2"
)


covariate_sets <- list(
  null="NA",
  date="^yday",
  autoreg="^N|mo",
  external="fetch|influx|water|wind",
  local="temp|salinity|shortwave|precip|km",
  cprn="attn|chl|dino|o2|ph|po4",
  all="."
)

#TODO: External, local with WeStCOMS2 = all sites


for(i in 3:length(covariate_sets)) {
  i.name <- names(covariate_sets)[i]
  i.covs <- covariate_sets[[i]]
  covars <- covars.all[grepl(i.covs, str_remove_all(covars.all, "\\.|_"))]
  covar_int <- grep(i.covs, covar_int.all, value=T)
  covar_s <- grep(i.covs, covar_s.all, value=T)
  covar_date <- grep(i.covs, c("ydayCos", "ydaySin"), value=T)
  
  # Smoothers
  s_b <- glue("b{covar_s}") 
  s_yday <- glue("b{covar_date}")
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
  
  # Laplace priors only
  flist.LP <- c(
    map(s_yday, ~as.formula(paste0(.x, "~1+(1|siteid)"))),
    map(s_b, ~as.formula(paste0(.x, "~s(ydayCos,ydaySin) + (1+s(ydayCos,ydaySin)|siteid)"))),
    map(1, ~"bIntercept~1 + (1|siteid)")
  )
  priors.LP <- c(
    prior_string("normal(0,1)", class="b", nlpar="bIntercept"),
    map(covar_s, 
        ~c(prior_string("double_exponential(0,0.1)", class="b", nlpar=paste0("b", .x)),
           prior_string("normal(0,.1)", class="sd", nlpar=paste0("b", .x), lb=0),
           prior_string("double_exponential(0,0.1)", class="sds", nlpar=paste0("b", .x), lb=0))) %>%
      do.call('c', .), 
    map(covar_date, ~prior_string("double_exponential(0,0.1)", class="b", nlpar=paste0("b", .x))) %>% 
      do.call('c', .)
  )
  
  # Indicator variable: p * b * X
  flist.P <- c(
    map(s_yday, ~as.formula(paste0(.x, "~1+(1|siteid)"))),
    map(s_b, ~as.formula(paste0(.x, "~s(ydayCos,ydaySin) + (1+s(ydayCos,ydaySin)|siteid)"))),
    map(s_p, ~as.formula(paste0(.x, "~ 1"))),
    map(1, ~"bIntercept~1 + (1|siteid)")
  )
  priors.P <- c(
    prior_string("normal(0,1)", class="b", nlpar="bIntercept"),
    map(covar_s, 
        ~c(prior_string("beta(0.1,1)", nlpar=paste0("p", .x), lb=0, ub=1),
           prior_string("normal(0,1)", class="b", nlpar=paste0("b", .x)),
           prior_string("normal(0,.1)", class="sd", nlpar=paste0("b", .x), lb=0),
           prior_string("double_exponential(0,0.1)", class="sds", nlpar=paste0("b", .x), lb=0))) %>%
      do.call('c', .),
    map(covar_date, ~prior_string("double_exponential(0,0.1)", class="b", nlpar=paste0("b", .x))) %>% 
      do.call('c', .)
  )
  
  
  
  
  
  
  # initial fit -------------------------------------------------------------
  
  for(sp in 1:length(species)) {
    target <- species[sp]
    target.tf <- thresh.df %>% filter(hab_parameter==target)
    
    target.df <- sampling.df %>%
      # filter(grid==1) %>%
      rename(N=!!target) %>%
      select(obs.id, site.id, date, hour, grid, lon, lat, fetch, bearing, N) %>%
      mutate(yday=yday(date),
             ydayCos=cos(2*pi*yday/365),
             ydaySin=sin(2*pi*yday/365),
             year=year(date),
             bearing=bearing*pi/180, # Nearest shore: N=pi, S=0
             # TODO: replace 'bearing' with longest direction via fetch layer
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
             one_of(covars, covar_int, covar_s, covar_date)) %>%
      filter(complete.cases(.)) %>%
      mutate(covarSet=i.name,
             NlnRAvg1=NA,
             NlnRAvg2=NA) %>%
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
    
    write_csv(target.df, glue("out{sep}test_full{sep}dataset_{i.name}_{target}.csv"))
    
    train.df <- target.df %>% filter(year <= 2017)
    test.df <- target.df %>% filter(year > 2017)
    bloomThresh <- max((!target.df$Nbloom)*target.df$NcatNum) # 1:4, maximum considered 'No bloom'
    
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
      form_ordLP <- makeFormula(cv_train.df, covar_s, "NcatNum | thres(3)", covar_date,
                                flist.LP, list(b=s_b, yday=s_yday, p=NULL))
      form_bernLP <- makeFormula(cv_train.df, covar_s, "Nbloom", covar_date,
                                 flist.LP, list(b=s_b, yday=s_yday, p=NULL))
      form_ordP <- makeFormula(cv_train.df, covar_s, "NcatNum | thres(3)", covar_date,
                               flist.P, list(b=s_b, yday=s_yday, p=s_p))
      form_bernP <- makeFormula(cv_train.df, covar_s, "Nbloom", covar_date,
                                flist.P, list(b=s_b, yday=s_yday, p=s_p))
      
      cv.ord <- brm(form_ord, data=cv_train.df,
                    family=cumulative("probit"), prior=priors, 
                    iter=iter, warmup=warmup, refresh=refresh, init=0,
                    control=ctrl, chains=chains, cores=chains,
                    file=glue("out{sep}test_full{sep}ord_CV{k}_{i.name}_{target}"))
      cv.ordP <- brm(form_ordP, data=cv_train.df, 
                     family=cumulative("probit"), prior=priors.P, 
                     iter=iter, warmup=warmup, refresh=refresh, init=0,
                     control=ctrl, chains=chains, cores=chains,
                     file=glue("out{sep}test_full{sep}ordP_CV{k}_{i.name}_{target}"))
      cv.ordLP <- brm(form_ordLP, data=cv_train.df, 
                      family=cumulative("probit"), prior=priors.LP, 
                      iter=iter, warmup=warmup, refresh=refresh, init=0,
                      control=ctrl, chains=chains, cores=chains,
                      file=glue("out{sep}test_full{sep}ordLP_CV{k}_{i.name}_{target}"))
      
      cv.bern01 <- brm(form_01, data=cv_train.df %>% filter(Nbloom1==0),
                       family=bernoulli("probit"), prior=priors, 
                       iter=iter, warmup=warmup, refresh=refresh, init=0,
                       control=ctrl, chains=chains, cores=chains,
                       file=glue("out{sep}test_full{sep}bern01_CV{k}_{i.name}_{target}"))
      cv.bernP01 <- brm(form_bernP, data=cv_train.df %>% filter(Nbloom1==0), 
                        family=bernoulli("probit"), prior=priors.P, 
                        iter=iter, warmup=warmup, refresh=refresh, init=0,
                        control=ctrl, chains=chains, cores=chains,
                        file=glue("out{sep}test_full{sep}bernP01_CV{k}_{i.name}_{target}"))
      cv.bernLP01 <- brm(form_bernLP, data=cv_train.df %>% filter(Nbloom1==0), 
                         family=bernoulli("probit"), prior=priors.LP, 
                         iter=iter, warmup=warmup, refresh=refresh, init=0,
                         control=ctrl, chains=chains, cores=chains,
                         file=glue("out{sep}test_full{sep}bernLP01_CV{k}_{i.name}_{target}"))
      cv.bern11 <- brm(form_11, data=cv_train.df %>% filter(Nbloom1==1),
                       family=bernoulli("probit"), prior=priors, 
                       iter=iter, warmup=warmup, refresh=refresh, init=0,
                       control=ctrl, chains=chains, cores=chains,
                       file=glue("out{sep}test_full{sep}bern11_CV{k}_{i.name}_{target}"))
      cv.bernP11 <- brm(form_bernP, data=cv_train.df %>% filter(Nbloom1==1), 
                        family=bernoulli("probit"), prior=priors.P, 
                        iter=iter, warmup=warmup, refresh=refresh, init=0,
                        control=ctrl, chains=chains, cores=chains,
                        file=glue("out{sep}test_full{sep}bernP11_CV{k}_{i.name}_{target}"))
      cv.bernLP11 <- brm(form_bernLP, data=cv_train.df %>% filter(Nbloom1==1), 
                         family=bernoulli("probit"), prior=priors.LP, 
                         iter=iter, warmup=warmup, refresh=refresh, init=0,
                         control=ctrl, chains=chains, cores=chains,
                         file=glue("out{sep}test_full{sep}bernLP11_CV{k}_{i.name}_{target}"))
      
      # RF
      rf_vars <- c("Nbloom", "Nbloom1", "lon", "lat", covar_date, covar_s)
      train.rf <- cv_train.df %>% select(one_of(rf_vars)) %>% 
        mutate(Nbloom=factor(Nbloom), lon=c(scale(lon)), lat=c(scale(lat))) %>%
        as.data.frame()
      train.rf01 <- cv_train.df %>% filter(Nbloom1==0) %>% select(one_of(rf_vars)) %>% 
        mutate(Nbloom=factor(Nbloom), lon=c(scale(lon)), lat=c(scale(lat))) %>%
        as.data.frame()
      train.rf11 <- cv_train.df %>% filter(Nbloom1==1) %>% select(one_of(rf_vars)) %>% 
        mutate(Nbloom=factor(Nbloom), lon=c(scale(lon)), lat=c(scale(lat))) %>%
        as.data.frame()
      test.rf <- cv_test.df %>% select(one_of(rf_vars)) %>% 
        mutate(Nbloom=factor(Nbloom), lon=c(scale(lon)), lat=c(scale(lat))) %>%
        as.data.frame()
      test.rf01 <- cv_test.df %>% filter(Nbloom1==0) %>% select(one_of(rf_vars)) %>% 
        mutate(Nbloom=factor(Nbloom), lon=c(scale(lon)), lat=c(scale(lat))) %>%
        as.data.frame()
      test.rf11 <- cv_test.df %>% filter(Nbloom1==1) %>% select(one_of(rf_vars)) %>% 
        mutate(Nbloom=factor(Nbloom), lon=c(scale(lon)), lat=c(scale(lat))) %>%
        as.data.frame()
      
      rf <- tuneRF(x=train.rf[,-1], y=train.rf[,1], doBest=T, trace=F, plot=F, ntree=1000)
      rf.01 <- tuneRF(x=train.rf01[,-1], y=train.rf01[,1], doBest=T, trace=F, plot=F, ntree=1000)
      rf.11 <- tuneRF(x=train.rf11[,-1], y=train.rf11[,1], doBest=T, trace=F, plot=F, ntree=1000)
      
      # Cross-validation predictions
      cv_test.df01 <- cv_test.df %>% filter(Nbloom1==0) %>% droplevels
      cv_test.df11 <- cv_test.df %>% filter(Nbloom1==1) %>% droplevels
      pred.ord <- posterior_epred(cv.ord, newdata=cv_test.df, allow_new_levels=T)
      pred.ordP <- posterior_epred(cv.ordP, newdata=cv_test.df, allow_new_levels=T)
      pred.ordLP <- posterior_epred(cv.ordLP, newdata=cv_test.df, allow_new_levels=T)
      pred.bern01 <- posterior_epred(cv.bern01, newdata=cv_test.df01, allow_new_levels=T) 
      pred.bernP01 <- posterior_epred(cv.bernP01, newdata=cv_test.df01, allow_new_levels=T) 
      pred.bernLP01 <- posterior_epred(cv.bernLP01, newdata=cv_test.df01, allow_new_levels=T) 
      pred.bern11 <- posterior_epred(cv.bern11, newdata=cv_test.df11, allow_new_levels=T) 
      pred.bernP11 <- posterior_epred(cv.bernP11, newdata=cv_test.df11, allow_new_levels=T) 
      pred.bernLP11 <- posterior_epred(cv.bernLP11, newdata=cv_test.df11, allow_new_levels=T) 
      
      cv_pred[[k]] <- full_join(
        cv_test.df %>%
          mutate(ord_mnpr=calc_ord_mnpr(pred.ord, bloomThresh),
                 ordP_mnpr=calc_ord_mnpr(pred.ordP, bloomThresh),
                 ordLP_mnpr=calc_ord_mnpr(pred.ordLP, bloomThresh),
                 rf_mnpr=predict(rf, newdata=test.rf, type="prob")[,2]),
        tibble(obsid=c(filter(cv_test.df, Nbloom1==0)$obsid, filter(cv_test.df, Nbloom1==1)$obsid),
               bern_mnpr=c(colMeans(pred.bern01), colMeans(pred.bern11)),
               bernP_mnpr=c(colMeans(pred.bernP01), colMeans(pred.bernP11)),
               bernLP_mnpr=c(colMeans(pred.bernLP01), colMeans(pred.bernLP11)),
               rf_split_mnpr=c(predict(rf.01, newdata=test.rf01, type="prob")[,2], 
                               predict(rf.11, newdata=test.rf11, type="prob")[,2]))
      ) %>%
        mutate(covarSet=i.name)
    }
    cv_pred %>% do.call('rbind', .) %>%
      write_csv(glue("out{sep}test_full{sep}CV_{i.name}_{target}.csv"))
    
    
    # Full fit for predictions
    
    # Formulas with interactions: errors if missing NcatF levels
    form_ord <- makeFormula(train.df, covar_int, "NcatNum | thres(3)")
    form_01 <- makeFormula(filter(train.df, Nbloom1==0), covar_int, "Nbloom")
    form_11 <- makeFormula(filter(train.df, Nbloom1==1), covar_int, "Nbloom")
    
    # Smoother formulas
    form_ordLP <- makeFormula(train.df, covar_s, "NcatNum | thres(3)", covar_date,
                              flist.LP, list(b=s_b, yday=s_yday, p=NULL))
    form_bernLP <- makeFormula(train.df, covar_s, "Nbloom", covar_date,
                               flist.LP, list(b=s_b, yday=s_yday, p=NULL))
    form_ordP <- makeFormula(train.df, covar_s, "NcatNum | thres(3)", covar_date,
                             flist.P, list(b=s_b, yday=s_yday, p=s_p))
    form_bernP <- makeFormula(train.df, covar_s, "Nbloom", covar_date,
                              flist.P, list(b=s_b, yday=s_yday, p=s_p))
    
    out.ord <- brm(form_ord, data=train.df,
                   family=cumulative("probit"), prior=priors, 
                   iter=iter, warmup=warmup, refresh=refresh, init=0,
                   control=ctrl, chains=chains, cores=chains,
                   file=glue("out{sep}test_full{sep}ord_{i.name}_{target}"))
    out.ordP <- brm(form_ordP, data=train.df, 
                    family=cumulative("probit"), prior=priors.P, 
                    iter=iter, warmup=warmup, refresh=refresh, init=0,
                    control=ctrl, chains=chains, cores=chains,
                    file=glue("out{sep}test_full{sep}ordP_{i.name}_{target}"))
    out.ordLP <- brm(form_ordLP, data=train.df, 
                     family=cumulative("probit"), prior=priors.LP, 
                     iter=iter, warmup=warmup, refresh=refresh, init=0,
                     control=ctrl, chains=chains, cores=chains,
                     file=glue("out{sep}test_full{sep}ordLP_{i.name}_{target}"))
    out.bern01 <- brm(form_01, data=train.df %>% filter(Nbloom1==0),
                      family=bernoulli("probit"), prior=priors, 
                      iter=iter, warmup=warmup, refresh=refresh, init=0,
                      control=ctrl, chains=chains, cores=chains,
                      file=glue("out{sep}test_full{sep}bern01_{i.name}_{target}"))
    out.bernP01 <- brm(form_bernP, data=train.df %>% filter(Nbloom1==0), 
                       family=bernoulli("probit"), prior=priors.P, 
                       iter=iter, warmup=warmup, refresh=refresh, init=0,
                       control=ctrl, chains=chains, cores=chains,
                       file=glue("out{sep}test_full{sep}bernP01_{i.name}_{target}"))
    out.bernLP01 <- brm(form_bernLP, data=train.df %>% filter(Nbloom1==0), 
                        family=bernoulli("probit"), prior=priors.LP, 
                        iter=iter, warmup=warmup, refresh=refresh, init=0,
                        control=ctrl, chains=chains, cores=chains,
                        file=glue("out{sep}test_full{sep}bernLP01_{i.name}_{target}"))
    out.bern11 <- brm(form_11, data=train.df %>% filter(Nbloom1==1),
                      family=bernoulli("probit"), prior=priors, 
                      iter=iter, warmup=warmup, refresh=refresh, init=0,
                      control=ctrl, chains=chains, cores=chains,
                      file=glue("out{sep}test_full{sep}bern11_{i.name}_{target}"))
    out.bernP11 <- brm(form_bernP, data=train.df %>% filter(Nbloom1==1), 
                       family=bernoulli("probit"), prior=priors.P, 
                       iter=iter, warmup=warmup, refresh=refresh, init=0,
                       control=ctrl, chains=chains, cores=chains,
                       file=glue("out{sep}test_full{sep}bernP11_{i.name}_{target}"))
    out.bernLP11 <- brm(form_bernLP, data=train.df %>% filter(Nbloom1==1), 
                        family=bernoulli("probit"), prior=priors.LP, 
                        iter=iter, warmup=warmup, refresh=refresh, init=0,
                        control=ctrl, chains=chains, cores=chains,
                        file=glue("out{sep}test_full{sep}bernLP11_{i.name}_{target}"))
    
    # RF
    rf_vars <- c("Nbloom", "Nbloom1", "lon", "lat", covar_date, covar_s)
    train.rf <- train.df %>% select(one_of(rf_vars)) %>% 
      mutate(Nbloom=factor(Nbloom), lon=c(scale(lon)), lat=c(scale(lat))) %>%
      as.data.frame()
    train.rf01 <- train.df %>% filter(Nbloom1==0) %>% select(one_of(rf_vars)) %>% 
      mutate(Nbloom=factor(Nbloom), lon=c(scale(lon)), lat=c(scale(lat))) %>%
      as.data.frame()
    train.rf11 <- train.df %>% filter(Nbloom1==1) %>% select(one_of(rf_vars)) %>% 
      mutate(Nbloom=factor(Nbloom), lon=c(scale(lon)), lat=c(scale(lat))) %>%
      as.data.frame()
    test.rf <- test.df %>% select(one_of(rf_vars)) %>% 
      mutate(Nbloom=factor(Nbloom), lon=c(scale(lon)), lat=c(scale(lat))) %>%
      as.data.frame()
    test.rf01 <- test.df %>% filter(Nbloom1==0) %>% select(one_of(rf_vars)) %>% 
      mutate(Nbloom=factor(Nbloom), lon=c(scale(lon)), lat=c(scale(lat))) %>%
      as.data.frame()
    test.rf11 <- test.df %>% filter(Nbloom1==1) %>% select(one_of(rf_vars)) %>% 
      mutate(Nbloom=factor(Nbloom), lon=c(scale(lon)), lat=c(scale(lat))) %>%
      as.data.frame()
    
    rf <- tuneRF(x=train.rf[,-1], y=train.rf[,1], doBest=T, trace=F, plot=F)
    rf.01 <- tuneRF(x=train.rf01[,-1], y=train.rf01[,1], doBest=T, trace=F, plot=F)
    rf.11 <- tuneRF(x=train.rf11[,-1], y=train.rf11[,1], doBest=T, trace=F, plot=F)
    saveRDS(rf, glue("out{sep}test_full{sep}rf_{i.name}_{target}.rds"))
    saveRDS(rf.01, glue("out{sep}test_full{sep}rf01_{i.name}_{target}.rds"))
    saveRDS(rf.11, glue("out{sep}test_full{sep}rf11_{i.name}_{target}.rds"))
    
    # Fitted
    fit.ord <- posterior_epred(out.ord)
    fit.ordP <- posterior_epred(out.ordP)
    fit.ordLP <- posterior_epred(out.ordLP)
    fit.bern01 <- posterior_epred(out.bern01) 
    fit.bernP01 <- posterior_epred(out.bernP01) 
    fit.bernLP01 <- posterior_epred(out.bernLP01) 
    fit.bern11 <- posterior_epred(out.bern11) 
    fit.bernP11 <- posterior_epred(out.bernP11) 
    fit.bernLP11 <- posterior_epred(out.bernLP11) 
    
    fit.df <- full_join(
      train.df %>%
        mutate(ord_mnpr=calc_ord_mnpr(fit.ord, bloomThresh),
               ordP_mnpr=calc_ord_mnpr(fit.ordP, bloomThresh),
               ordLP_mnpr=calc_ord_mnpr(fit.ordLP, bloomThresh),
               rf_mnpr=rf$votes[,2]),
      tibble(obsid=c(filter(train.df, Nbloom1==0)$obsid, filter(train.df, Nbloom1==1)$obsid),
             bern_mnpr=c(colMeans(fit.bern01), colMeans(fit.bern11)),
             bernP_mnpr=c(colMeans(fit.bernP01), colMeans(fit.bernP11)),
             bernLP_mnpr=c(colMeans(fit.bernLP01), colMeans(fit.bernLP11)),
             rf_split_mnpr=c(rf.01$votes[,2], rf.11$votes[,2]))
    ) %>%
      mutate(covarSet=i.name)
    write_csv(fit.df, glue("out{sep}test_full{sep}fit_{i.name}_{target}.csv"))
    
    # OOS predictions
    test.df01 <- test.df %>% filter(Nbloom1==0) %>% droplevels
    test.df11 <- test.df %>% filter(Nbloom1==1) %>% droplevels
    pred.ord <- posterior_epred(out.ord, newdata=test.df, allow_new_levels=T)
    pred.ordP <- posterior_epred(out.ordP, newdata=test.df, allow_new_levels=T)
    pred.ordLP <- posterior_epred(out.ordLP, newdata=test.df, allow_new_levels=T)
    pred.bern01 <- posterior_epred(out.bern01, newdata=test.df01, allow_new_levels=T) 
    pred.bernP01 <- posterior_epred(out.bernP01, newdata=test.df01, allow_new_levels=T) 
    pred.bernLP01 <- posterior_epred(out.bernLP01, newdata=test.df01, allow_new_levels=T) 
    pred.bern11 <- posterior_epred(out.bern11, newdata=test.df11, allow_new_levels=T) 
    pred.bernP11 <- posterior_epred(out.bernP11, newdata=test.df11, allow_new_levels=T) 
    pred.bernLP11 <- posterior_epred(out.bernLP11, newdata=test.df11, allow_new_levels=T) 
    
    pred.df <- full_join(
      test.df %>%
        mutate(ord_mnpr=calc_ord_mnpr(pred.ord, bloomThresh),
               ordP_mnpr=calc_ord_mnpr(pred.ordP, bloomThresh),
               ordLP_mnpr=calc_ord_mnpr(pred.ordLP, bloomThresh),
               rf_mnpr=predict(rf, newdata=test.rf, type="prob")[,2]),
      tibble(obsid=c(filter(test.df, Nbloom1==0)$obsid, filter(test.df, Nbloom1==1)$obsid),
             bern_mnpr=c(colMeans(pred.bern01), colMeans(pred.bern11)),
             bernP_mnpr=c(colMeans(pred.bernP01), colMeans(pred.bernP11)),
             bernLP_mnpr=c(colMeans(pred.bernLP01), colMeans(pred.bernLP11)),
             rf_split_mnpr=c(predict(rf.01, newdata=test.rf01, type="prob")[,2], 
                             predict(rf.11, newdata=test.rf11, type="prob")[,2]))
    ) %>%
      mutate(covarSet=i.name)
    
    write_csv(pred.df, glue("out{sep}test_full{sep}pred_{i.name}_{target}.csv"))
    
    cat("Finished", target, "\n")
  }
  
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

# for(sp in 1:length(species)) {
#   target <- species[sp]
#   target.tf <- thresh.df %>% filter(hab_parameter==target)
#   
#   target.df <- sampling.df %>%
#     rename(N=!!target) %>%
#     select(obs.id, site.id, date, hour, grid, lon, lat, fetch, bearing, N) %>%
#     mutate(yday=yday(date),
#            ydayCos=cos(2*pi*yday/365),
#            ydaySin=sin(2*pi*yday/365),
#            year=year(date),
#            bearing=bearing*pi/180,
#            N=round(N),
#            N.ln=log(N+1)) %>%
#     rowwise() %>%
#     mutate(N.cat=target.tf$tl[max(which(N >= target.tf$min_ge))]) %>%
#     ungroup %>%
#     mutate(N.bloom=target.tf$bloom[match(N.cat, target.tf$tl)]) %>%
#     arrange(site.id, date) %>%
#     group_by(site.id) %>%
#     multijetlag(N.ln, N.cat, N.bloom, date, n=2) %>%
#     ungroup %>%
#     arrange(site.id, date) %>%
#     rename_with(~str_remove_all(.x, "\\.|_"))
#   
#   write_csv(target.df, glue("out{sep}dataset_{target}.csv"))
# }
# obs.df <- map_dfr(species, ~read_csv(glue("out{sep}dataset_{.x}.csv")) %>%
#                     mutate(species=.x))
# 
# obs.df %>% mutate(week=week(date)) %>%
#   group_by(week, species) %>%
#   summarise(mn=mean(Nbloom),
#             q_lo=Hmisc::binconf(sum(Nbloom), n(), alpha=0.05)[1,2],
#             q_hi=Hmisc::binconf(sum(Nbloom), n(), alpha=0.05)[1,3],
#             N=n()) %>%
#   ungroup %>%
#   mutate(propN=N/max(N)) %>%
#   ggplot(aes(week)) +
#   geom_point(aes(y=mn, alpha=N)) +
#   geom_linerange(aes(ymin=q_lo, ymax=q_hi, alpha=N)) +
#   facet_wrap(~species, nrow=1) + theme_classic() +
#   scale_x_continuous(breaks=seq(1,50,length.out=5), labels=c("Jan", "Apr", "Jul", "Oct", "Dec")) +
#   scale_alpha_continuous(range=c(0.25, 1)) +
#   ylim(0,1) +
#   labs(x="Week", y="Proportion of sites with blooms") +
#   theme(strip.text=element_text(size=12))
# ggsave(glue("figs{sep}obs_bloom_proportion.png"), width=12, height=3, dpi=300)
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
