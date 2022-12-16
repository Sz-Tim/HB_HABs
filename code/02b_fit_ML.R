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

sp <- 1



# covariates --------------------------------------------------------------

species <- sp.i$full

covars.all <- c(
  "ydaySin", "ydayCos",
  # WeStCOMS weekly averages
  "temp_L_wk", "salinity_L_wk", "short_wave_L_wk", "km_L_wk",
  "precip_L_wk", "tempStrat20m_L_wk", "tempStrat20m_R_wk",
  "wind_L_wk", 
  "water_L_wk", "waterVelL",
  "water_R_wk", "waterVelR",
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


i <- 1
# for(i in length(covariate_sets)) {
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
    mutate(Nbloom=factor(Nbloom, levels=0:1, labels=c("X0", "X1"))) %>% as.data.frame()
  train.ML01 <- train.df %>% filter(Nbloom1==0) %>% select(one_of(ML_vars)) %>% 
    mutate(Nbloom=factor(Nbloom, levels=0:1, labels=c("X0", "X1"))) %>% as.data.frame()
  train.ML11 <- train.df %>% filter(Nbloom1==1) %>% select(one_of(ML_vars)) %>% 
    mutate(Nbloom=factor(Nbloom, levels=0:1, labels=c("X0", "X1"))) %>% as.data.frame()
  test.ML <- test.df %>% select(one_of(ML_vars)) %>% 
    mutate(Nbloom=factor(Nbloom, levels=0:1, labels=c("X0", "X1"))) %>% as.data.frame()
  test.ML01 <- test.df %>% filter(Nbloom1==0) %>% select(one_of(ML_vars)) %>% 
    mutate(Nbloom=factor(Nbloom, levels=0:1, labels=c("X0", "X1"))) %>% as.data.frame()
  test.ML11 <- test.df %>% filter(Nbloom1==1) %>% select(one_of(ML_vars)) %>% 
    mutate(Nbloom=factor(Nbloom, levels=0:1, labels=c("X0", "X1"))) %>% as.data.frame()

  folds <- createFoldsByYear(train.df)
  folds01 <- createFoldsByYear(train.df %>% filter(Nbloom1==0))
  folds11 <- createFoldsByYear(train.df %>% filter(Nbloom1==1))
  
  ctrl <- trainControl("repeatedcv", repeats=3, classProbs=T, 
                       index=folds$i.in, indexOut=folds$i.out)
  ctrl01 <- trainControl("repeatedcv", repeats=3, classProbs=T, 
                         index=folds01$i.in, indexOut=folds01$i.out)
  ctrl11 <- trainControl("repeatedcv", repeats=3, classProbs=T, 
                         index=folds11$i.in, indexOut=folds11$i.out)
  rf.tuneGrid <- expand.grid(.mtry=1:(ncol(train.ML)/2))
  rf_ <- train(Nbloom~., data=train.ML, method="rf", 
               trControl=ctrl, tuneGrid=rf.tuneGrid, ntree=2000)
  rf_01 <- train(Nbloom~., data=train.ML01 %>% select(-Nbloom1), method="rf", 
               trControl=ctrl01, tuneGrid=rf.tuneGrid, ntree=2000)
  rf_11 <- train(Nbloom~., data=train.ML11 %>% select(-Nbloom1), method="rf", 
               trControl=ctrl11, tuneGrid=rf.tuneGrid, ntree=2000)
  saveRDS(rf_, glue("{f.prefix}rf_{f.suffix}.rds"))
  saveRDS(rf_01, glue("{f.prefix}rf01_{f.suffix}.rds"))
  saveRDS(rf_11, glue("{f.prefix}rf11_{f.suffix}.rds"))
  cat("Finished rf:", target, "\n")
  
  
  svm.tuneGrid <- expand.grid(C=2^(seq(-5,5,length.out=100)))
  svm_ <- train(Nbloom~., data=train.ML, 
                method="svmRadialCost", trControl=ctrl, tuneGrid=svm.tuneGrid)
  svm_01 <- train(Nbloom~., data=train.ML01 %>% select(-Nbloom1), 
                  method="svmRadialCost", trControl=ctrl01, tuneGrid=svm.tuneGrid)
  svm_11 <- train(Nbloom~., data=train.ML11 %>% select(-Nbloom1), 
                  method="svmRadialCost", trControl=ctrl11, tuneGrid=svm.tuneGrid)
  saveRDS(svm_, glue("{f.prefix}svm_{f.suffix}.rds"))
  saveRDS(svm_01, glue("{f.prefix}svm01_{f.suffix}.rds"))
  saveRDS(svm_11, glue("{f.prefix}svm11_{f.suffix}.rds"))
  cat("Finished svm:", target, "\n")
  

  xgb.tuneGrid <- expand.grid(min_child_weight=1, subsample=0.8, colsample_bytree=0.8,
                              eta=2^(seq(-10,0,length.out=10)),
                              max_depth=2:5,
                              nrounds=seq(50, 200, length.out=4),
                              gamma=2^(seq(-10,2,length.out=10)))
  xgb_ <- train(Nbloom~., data=train.ML,
                method="xgbTree", trControl=ctrl, tuneGrid=xgb.tuneGrid, verbosity=0)
  xgb_01 <- train(Nbloom~., data=train.ML01,
                  method="xgbTree", trControl=ctrl01, tuneGrid=xgb.tuneGrid, verbosity=0)
  xgb_11 <- train(Nbloom~., data=train.ML11,
                  method="xgbTree", trControl=ctrl11, tuneGrid=xgb.tuneGrid, verbosity=0)
  saveRDS(xgb_, glue("{f.prefix}xgb_{f.suffix}.rds"))
  saveRDS(xgb_01, glue("{f.prefix}xgb01_{f.suffix}.rds"))
  saveRDS(xgb_11, glue("{f.prefix}xgb11_{f.suffix}.rds"))
  cat("Finished xgb:", target, "\n")
  
  # Fitted
  fit.df <- full_join(
    train.df %>%
      mutate(rf_mnpr=predict(rf_, train.ML, type="prob")[,2],
             svm_mnpr=predict(svm_, train.ML, type="prob")[,2],
             xgb_mnpr=predict(xgb_, train.ML, type="prob")[,2]),
    tibble(obsid=c(filter(train.df, Nbloom1==0)$obsid, filter(train.df, Nbloom1==1)$obsid),
           rf_split_mnpr=c(predict(rf_01, train.ML01, type="prob")[,2],
                           predict(rf_11, train.ML11, type="prob")[,2]),
           svm_split_mnpr=c(predict(svm_01, train.ML01, type="prob")[,2],
                            predict(svm_11, train.ML11, type="prob")[,2]),
           xgb_split_mnpr=c(predict(xgb_01, train.ML01, type="prob")[,2],
                            predict(xgb_11, train.ML11, type="prob")[,2]))
  ) %>%
    mutate(covarSet=i.name,
           species=target)
  write_csv(fit.df, glue("{f.prefix}fit_ML_{f.suffix}.csv"))
  
  # OOS predictions
  pred.rf <- predict(rf_, test.ML, type="prob")[,2]
  pred.rf_01 <- predict(rf_01, test.ML01, type="prob")[,2]
  pred.svm <- predict(svm_, test.ML, type="prob")[,2]
  pred.svm_01 <- predict(svm_01, test.ML01, type="prob")[,2]
  pred.xgb <- predict(xgb_, test.ML, type="prob")[,2]
  pred.xgb_01 <- predict(xgb_01, test.ML01, type="prob")[,2]
  if(nrow(test.ML11) > 0) {
    pred.rf_11 <- predict(rf_11, test.ML11, type="prob")[,2]
    pred.svm_11 <- predict(svm_11, test.ML11, type="prob")[,2]
    pred.xgb_11 <- predict(xgb_11, test.ML11, type="prob")[,2]
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
  yrCV <- unique(filter(train.df, year>2014)$year)
  cv_pred <- map(yrCV, ~NULL)
  for(k in 1:length(yrCV)) {
    cat("Starting cross-validation:", k, "of", length(yrCV), "\n")
    yr <- yrCV[k]
    cv_train.df <- train.df %>% filter(year != yr)
    cv_test.df <- train.df %>% filter(year == yr)
    
    # Machine Learning: Random Forest, Support Vector Machine
    ML_vars <- c("Nbloom", "Nbloom1", "lon_sc", "lat_sc", covars)
    train.ML <- cv_train.df %>% select(one_of(ML_vars)) %>% 
      mutate(Nbloom=factor(Nbloom, levels=0:1, labels=c("X0", "X1"))) %>% as.data.frame()
    train.ML01 <- cv_train.df %>% filter(Nbloom1==0) %>% select(one_of(ML_vars)) %>% 
      mutate(Nbloom=factor(Nbloom, levels=0:1, labels=c("X0", "X1"))) %>% as.data.frame()
    train.ML11 <- cv_train.df %>% filter(Nbloom1==1) %>% select(one_of(ML_vars)) %>% 
      mutate(Nbloom=factor(Nbloom, levels=0:1, labels=c("X0", "X1"))) %>% as.data.frame()
    test.ML <- cv_test.df %>% select(one_of(ML_vars)) %>% 
      mutate(Nbloom=factor(Nbloom, levels=0:1, labels=c("X0", "X1"))) %>%  as.data.frame()
    test.ML01 <- cv_test.df %>% filter(Nbloom1==0) %>% select(one_of(ML_vars)) %>% 
      mutate(Nbloom=factor(Nbloom, levels=0:1, labels=c("X0", "X1"))) %>% as.data.frame()
    test.ML11 <- cv_test.df %>% filter(Nbloom1==1) %>% select(one_of(ML_vars)) %>% 
      mutate(Nbloom=factor(Nbloom, levels=0:1, labels=c("X0", "X1"))) %>% as.data.frame()
    
    rf_ <- train(Nbloom~., data=train.ML, method="rf", 
                 trControl=rf.ctrl, tuneGrid=rf.tuneGrid, ntree=2000)
    rf_01 <- train(Nbloom~., data=train.ML01 %>% select(-Nbloom1), method="rf", 
                   trControl=rf.ctrl, tuneGrid=rf.tuneGrid, ntree=2000)
    rf_11 <- train(Nbloom~., data=train.ML11 %>% select(-Nbloom1), method="rf", 
                   trControl=rf.ctrl, tuneGrid=rf.tuneGrid, ntree=2000)
    
    svm_ <- train(Nbloom~., data=train.ML, 
                  method="svmRadialCost", trControl=svm.ctrl, tuneGrid=svm.tuneGrid)
    svm_01 <- train(Nbloom~., data=train.ML01 %>% select(-Nbloom1), 
                    method="svmRadialCost", trControl=svm.ctrl, tuneGrid=svm.tuneGrid)
    svm_11 <- train(Nbloom~., data=train.ML11 %>% select(-Nbloom1), 
                    method="svmRadialCost", trControl=svm.ctrl, tuneGrid=svm.tuneGrid)
    
    xgb_ <- train(Nbloom~., data=train.ML,
                  method="xgbTree", trControl=xgb.ctrl, tuneGrid=xgb.tuneGrid, verbosity=0)
    xgb_01 <- train(Nbloom~., data=train.ML01,
                    method="xgbTree", trControl=xgb.ctrl, tuneGrid=xgb.tuneGrid, verbosity=0)
    xgb_11 <- train(Nbloom~., data=train.ML11,
                    method="xgbTree", trControl=xgb.ctrl, tuneGrid=xgb.tuneGrid, verbosity=0)
    
    
    # Cross-validation predictions
    pred.rf <- predict(rf_, test.ML, type="prob")[,2]
    pred.rf_01 <- predict(rf_01, test.ML01, type="prob")[,2]
    pred.svm <- predict(svm_, test.ML, type="prob")[,2]
    pred.svm_01 <- predict(svm_01, test.ML01, type="prob")[,2]
    pred.xgb <- predict(xgb_, test.ML, type="prob")[,2]
    pred.xgb_01 <- predict(xgb_01, test.ML01, type="prob")[,2]
    if(nrow(test.ML11) > 0) {
      pred.rf_11 <- predict(rf_11, test.ML11, type="prob")[,2]
      pred.svm_11 <- predict(svm_11, test.ML11, type="prob")[,2]
      pred.xgb_11 <- predict(xgb_11, test.ML11, type="prob")[,2]
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
# }





