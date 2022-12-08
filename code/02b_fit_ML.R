# HABReports Bayesian modelling
# Bayesian HAB model initial fits
# Tim Szewczyk


# This script fits models with data from 2013-2017 using horseshoe priors



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue",
          "randomForest", "caret", "e1071", "xgboost")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
walk(dir("code", "*00_fn", full.names=T), source)
options(mc.cores=10)
cores <- 2

# Model details

# minch2:    2013-06-20 to 2019-07-02
# WeStCOMS2: 2019-04-01 to 2022-01-26
sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")

sp.i <- read_csv("data/sp_i.csv")

sp <- 5



# covariates --------------------------------------------------------------

species <- sp.i$full

covars.all <- c(
  # WeStCOMS weekly averages
  "temp_L_wk", "salinity_L_wk", "short_wave_L_wk", "km_L_wk",
  "precip_L_wk", "tempStrat20m_L_wk", "tempStrat20m_R_wk",
  "wind_L_wk", "windDir_L_wk",
  "water_L_wk", "waterDir_L_wk", "waterVelL",
  "water_R_wk", "waterDir_R_wk", "waterVelR",
  # external forcings          
  "fetch", "influx_wk",
  # Copernicus
  "attn_wk", "chl_wk", "dino_wk", "o2_wk", "ph_wk", "po4_wk",
  # detrended WeStCOMS, Copernicus (cor < 0.8 with original)          
  "temp_L_wk_dt", "salinity_L_wk_dt", "short_wave_L_wk_dt", 
  "precip_L_wk_dt", "tempStrat20m_L_wk_dt", "tempStrat20m_R_wk_dt",
  "wind_L_wk_dt", "water_L_wk_dt", "water_R_wk_dt", 
  "chl_wk_dt", "o2_wk_dt", "ph_wk_dt", "po4_wk_dt")

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
  covars <- covars.all[grepl(i.covs, str_remove_all(covars.all, "\\.|_"))] %>%
    str_remove_all("\\.|_")
  
  # initial fit -------------------------------------------------------------
  
  target <- species[sp]
  f.prefix <- glue("out{sep}1-loose{sep}")
  f.suffix <- glue("{i.name}_{target}")
  
  target.df <- read_csv(glue("{f.prefix}dataset_{i.name}_{target}.csv"))
  
  train.df <- target.df %>% filter(year <= 2019)
  test.df <- target.df %>% filter(year > 2019)
  bloomThresh <- max((!target.df$Nbloom)*target.df$NcatNum) # 1:4, maximum considered 'No bloom'
  
  # Full fit for predictions
  # Machine Learning: Random Forest, Support Vector Machine, XGBoost
  ML_vars <- c("Nbloom", "Nbloom1", "lon_sc", "lat_sc", covars)
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
  saveRDS(rf, glue("{f.prefix}rf_{f.suffix}.rds"))
  saveRDS(rf.01, glue("{f.prefix}rf01_{f.suffix}.rds"))
  saveRDS(rf.11, glue("{f.prefix}rf11_{f.suffix}.rds"))
  cat("Finished rf:", target, "\n")
  
  #svm_rng <- list(epsilon=seq(0,0.5,0.05), cost=2^(seq(3,7,0.2)))
  svm_rng <- list(cost=10^(-2:2), gamma=2^(-4:2))
  svm_ <- tune(svm, train.ML[,-1], train.ML[,1], probability=T, ranges=svm_rng)
  svm_01 <- tune(svm, train.ML01[,-1], train.ML01[,1], probability=T, ranges=svm_rng)
  svm_11 <- tune(svm, train.ML11[,-1], train.ML11[,1], probability=T, ranges=svm_rng)
  saveRDS(svm_, glue("{f.prefix}svm_{f.suffix}.rds"))
  saveRDS(svm_01, glue("{f.prefix}svm01_{f.suffix}.rds"))
  saveRDS(svm_11, glue("{f.prefix}svm11_{f.suffix}.rds"))
  cat("Finished svm:", target, "\n")
  
  xg <- best_xgb(train.ML, nthread=cores)
  xg.01 <- best_xgb(train.ML01, nthread=cores)
  xg.11 <- best_xgb(train.ML11, nthread=cores)
  saveRDS(xg, glue("{f.prefix}xgb_{f.suffix}.rds"))
  saveRDS(xg.01, glue("{f.prefix}xgb01_{f.suffix}.rds"))
  saveRDS(xg.11, glue("{f.prefix}xgb11_{f.suffix}.rds"))
  cat("Finished xgb:", target, "\n")
  
  # Fitted
  fit.df <- full_join(
    train.df %>%
      mutate(rf_mnpr=rf$votes[,2],
             svm_mnpr=attr(predict(svm_$best.model, newdata=train.ML[,-1], probability=T),
                           "probabilities")[,2],
             xgb_mnpr=predict(xg, as.matrix(train.ML[,-1]))),
    tibble(obsid=c(filter(train.df, Nbloom1==0)$obsid, filter(train.df, Nbloom1==1)$obsid),
           rf_split_mnpr=c(rf.01$votes[,2], rf.11$votes[,2]),
           svm_split_mnpr=c(attr(predict(svm_01$best.model, newdata=train.ML01[,-1], probability=T),
                                 "probabilities")[,2],
                            attr(predict(svm_11$best.model, newdata=train.ML11[,-1], probability=T),
                                 "probabilities")[,2]),
           xgb_split_mnpr=c(predict(xg.01, as.matrix(train.ML01[,-1])),
                            predict(xg.11, as.matrix(train.ML11[,-1]))))
  ) %>%
    mutate(covarSet=i.name,
           species=target)
  write_csv(fit.df, glue("{f.prefix}fit_ML_{f.suffix}.csv"))
  
  # OOS predictions
  test.df01 <- test.df %>% filter(Nbloom1==0) %>% droplevels
  test.df11 <- test.df %>% filter(Nbloom1==1) %>% droplevels
  pred.rf <- predict(rf, newdata=test.ML, type="prob")[,2]
  pred.rf_01 <- predict(rf.01, newdata=test.ML01, type="prob")[,2]
  pred.svm <- attr(predict(svm_$best.model, newdata=test.ML[,-1], probability=T),
                   "probabilities")[,2]
  pred.svm_01 <- attr(predict(svm_01$best.model, newdata=test.ML01[,-1], probability=T),
                      "probabilities")[,2]
  pred.xgb <- predict(xg, as.matrix(test.ML[,-1]))
  pred.xgb_01 <- predict(xg.01, as.matrix(test.ML01[,-1]))
  if(nrow(test.df11) > 0) {
    pred.rf_11 <- predict(rf.11, newdata=test.ML11, type="prob")[,2]
    pred.svm_11 <- attr(predict(svm_11$best.model, newdata=test.ML11[,-1], probability=T),
                        "probabilities")[,2]
    pred.xgb_11 <- predict(xg.11, as.matrix(test.ML11[,-1]))
  } else {
    pred.rf_11 <- numeric(0)
    pred.svm_11 <- numeric(0)
    pred.xgb_11 <- numeric(0)
  }
  
  pred.df <- full_join(
    test.df %>%
      mutate(rf_mnpr=pred.rf,
             svm_mnpr=pred.svm,
             xgb_mnpr=pred.xgb),
    tibble(obsid=c(filter(test.df, Nbloom1==0)$obsid, filter(test.df, Nbloom1==1)$obsid),
           rf_split_mnpr=c(pred.rf_01, pred.rf_11),
           svm_split_mnpr=c(pred.svm_01, pred.svm_11),
           xgb_split_mnpr=c(pred.xgb_01, pred.xgb_11))
  ) %>%
    mutate(covarSet=i.name,
           species=target)
  
  write_csv(pred.df, glue("{f.prefix}pred_ML_{f.suffix}.csv"))
  
  
  # Cross-validation by year
  yrCV <- unique(train.df$year)
  cv_pred <- map(yrCV, ~NULL)
  for(k in 1:length(yrCV)) {
    cat("Starting cross-validation:", k, "of", length(yrCV), "\n")
    yr <- yrCV[k]
    cv_train.df <- train.df %>% filter(year != yr)
    cv_test.df <- train.df %>% filter(year == yr)
    
    # Machine Learning: Random Forest, Support Vector Machine
    ML_vars <- c("Nbloom", "Nbloom1", "lon_sc", "lat_sc", covars)
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
    
    svm_rng <- list(cost=10^(-2:2), gamma=2^(-4:2))
    svm_ <- tune(svm, train.ML[,-1], train.ML[,1], probability=T, ranges=svm_rng, kernel="radial")
    svm_01 <- tune(svm, train.ML01[,-1], train.ML01[,1], probability=T, ranges=svm_rng, kernel="radial")
    svm_11 <- tune(svm, train.ML11[,-1], train.ML11[,1], probability=T, ranges=svm_rng, kernel="radial")
    
    xg <- best_xgb(train.ML, nthread=cores)
    xg.01 <- best_xgb(train.ML01, nthread=cores)
    xg.11 <- best_xgb(train.ML11, nthread=cores)
    
    
    # Cross-validation predictions
    pred.rf <- predict(rf, newdata=test.ML, type="prob")[,2]
    pred.rf_01 <- predict(rf.01, newdata=test.ML01, type="prob")[,2]
    pred.svm <- attr(predict(svm_$best.model, newdata=test.ML[,-1], probability=T),
                     "probabilities")[,2]
    pred.svm_01 <- attr(predict(svm_01$best.model, newdata=test.ML01[,-1], probability=T),
                        "probabilities")[,2]
    pred.xgb <- predict(xg, as.matrix(test.ML[,-1]))
    pred.xgb_01 <- predict(xg.01, as.matrix(test.ML01[,-1]))
    if(nrow(test.ML11) > 0) {
      pred.rf_11 <- predict(rf.11, newdata=test.ML11, type="prob")[,2]
      pred.svm_11 <- attr(predict(svm_11$best.model, newdata=test.ML11[,-1], probability=T),
                          "probabilities")[,2]
      pred.xgb_11 <- predict(xg.11, as.matrix(test.ML11[,-1]))
    } else {
      pred.rf_11 <- numeric(0)
      pred.svm_11 <- numeric(0)
      pred.xgb_11 <- numeric(0)
    }
    
    cv_pred[[k]] <- full_join(
      cv_test.df %>%
        mutate(rf_mnpr=pred.rf,
               svm_mnpr=pred.svm,
               xgb_mnpr=pred.xgb),
      tibble(obsid=c(filter(cv_test.df, Nbloom1==0)$obsid, filter(cv_test.df, Nbloom1==1)$obsid),
             rf_split_mnpr=c(pred.rf_01, pred.rf_11),
             svm_split_mnpr=c(pred.svm_01, pred.svm_11),
             xgb_split_mnpr=c(pred.xgb_01, pred.xgb_11))
    ) %>%
      mutate(covarSet=i.name,
             species=target)
  }
  cv_pred %>% do.call('rbind', .) %>%
    write_csv(glue("{f.prefix}CV_ML_{f.suffix}.csv"))
  
  
  cat("Finished", target, "\n")
}





