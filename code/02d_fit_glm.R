# HABReports Bayesian modelling
# Bayesian HAB model initial fits
# Tim Szewczyk


# This script fits models with data from 2013-2017 using horseshoe priors



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "glmnet", "caret")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
walk(dir("code", "*00_fn", full.names=T), source)

# Model details

# minch2:    2013-06-20 to 2019-07-02
# WeStCOMS2: 2019-04-01 to 2022-01-26
sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")

sp.i <- read_csv("data/sp_i.csv")




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
  covar_int <- grep(i.covs, covar_int.all, value=T)
  
  # initial fit -------------------------------------------------------------
  
  for(sp in 1:5) {
    target <- species[sp]
    f.prefix <- glue("out{sep}1-loose{sep}")
    f.suffix <- glue("{i.name}_{target}")
    
    target.df <- read_csv(glue("{f.prefix}dataset_{i.name}_{target}.csv"))
    
    train.df <- target.df %>% filter(year <= 2019)
    test.df <- target.df %>% filter(year > 2019)
    bloomThresh <- max((!target.df$Nbloom)*target.df$NcatNum) # 1:4, maximum considered 'No bloom'
    
    # Full fit for predictions
    train.df_ <- train.df %>% 
      mutate(lonlat=lon_sc*lat_sc,
             Nbloom=factor(Nbloom, labels=c("X0", "X1"))) %>%
      select(Nbloom, all_of(covar_int), lonlat, Nbloom1)
    train.df01 <- train.df_ %>% filter(Nbloom1==0)
    train.df11 <- train.df_ %>% filter(Nbloom1==1)
    
    elast.ctrl <- trainControl("repeatedcv", number=6, repeats=3, classProbs=T)
    elast.grid <- expand.grid(alpha=seq(0, 1, length.out=26),
                              lambda=2^(seq(-15,-1,length.out=25)))
    elast_ <- train(Nbloom~., data=train.df_, 
                    method="glmnet", trControl=elast.ctrl, tuneGrid=elast.grid)
    elast01 <- train(Nbloom~., data=train.df01, 
                     method="glmnet", trControl=elast.ctrl, tuneGrid=elast.grid)
    elast11 <- train(Nbloom~., data=train.df11, 
                     method="glmnet", trControl=elast.ctrl, tuneGrid=elast.grid)
    
    saveRDS(elast_, glue("{f.prefix}glmElast_{f.suffix}.rds"))
    saveRDS(elast01, glue("{f.prefix}glmElast01_{f.suffix}.rds"))
    saveRDS(elast11, glue("{f.prefix}glmElast11_{f.suffix}.rds"))
    
    # Fitted
    fit.df <- full_join(
      train.df %>%
        mutate(glmElast_mnpr=predict(elast_, train.df_, type="prob")[,2]),
      tibble(obsid=c(filter(train.df, Nbloom1==0)$obsid, filter(train.df, Nbloom1==1)$obsid),
             glmElast_split_mnpr=c(
               predict(elast01, train.df01, type="prob")[,2],
               predict(elast11, train.df11, type="prob")[,2]))
    ) %>%
      mutate(covarSet=i.name,
             species=target)
    write_csv(fit.df, glue("{f.prefix}fit_glm_{f.suffix}.csv"))
    
    # OOS predictions
    test.df_ <- test.df %>% 
      mutate(lonlat=lon_sc*lat_sc,
             Nbloom=factor(Nbloom, labels=c("X0", "X1"))) %>%
      select(Nbloom, all_of(covar_int), lonlat, Nbloom1) 
    test.df01 <- test.df_ %>% filter(Nbloom1==0)
    test.df11 <- test.df_ %>% filter(Nbloom1==1)
    pred.elast <- predict(elast_, test.df_, type="prob")[,2]
    pred.elast_01 <- predict(elast01, test.df01, type="prob")[,2]
    if(nrow(test.df11) > 0) {
      pred.elast_11 <- predict(elast11, test.df11, type="prob")[,2]
    } else {
      pred.elast_11 <- numeric(0)
    }
    
    pred.df <- full_join(
      test.df %>%
        mutate(glmElast_mnpr=pred.elast),
      tibble(obsid=c(filter(test.df, Nbloom1==0)$obsid, filter(test.df, Nbloom1==1)$obsid),
             glmElast_split_mnpr=c(pred.elast_01, pred.elast_11))
    ) %>%
      mutate(covarSet=i.name,
             species=target)
    
    write_csv(pred.df, glue("{f.prefix}pred_glm_{f.suffix}.csv"))
    
    
    # Cross-validation by year
    yrCV <- unique(filter(train.df, year>2014)$year)
    cv_pred <- map(yrCV, ~NULL)
    for(k in 1:length(yrCV)) {
      yr <- yrCV[k]
      cv_train.df <- train.df %>% filter(year != yr)
      cv_test.df <- train.df %>% filter(year == yr)
      
      train.df_ <- cv_train.df %>% 
        mutate(lonlat=lon_sc*lat_sc,
               Nbloom=factor(Nbloom, labels=c("X0", "X1"))) %>%
        select(Nbloom, all_of(covar_int), lonlat, Nbloom1)
      train.df01 <- train.df_ %>% filter(Nbloom1==0)
      train.df11 <- train.df_ %>% filter(Nbloom1==1)
      
      elast_ <- train(Nbloom~., data=train.df_, 
                      method="glmnet", trControl=elast.ctrl, tuneGrid=elast.grid)
      elast01 <- train(Nbloom~., data=train.df01, 
                       method="glmnet", trControl=elast.ctrl, tuneGrid=elast.grid)
      elast11 <- train(Nbloom~., data=train.df11, 
                       method="glmnet", trControl=elast.ctrl, tuneGrid=elast.grid)
      
      # Cross-validation predictions
      test.df_ <- cv_test.df %>% 
        mutate(lonlat=lon_sc*lat_sc,
               Nbloom=factor(Nbloom, labels=c("X0", "X1"))) %>%
        select(Nbloom, all_of(covar_int), lonlat, Nbloom1) 
      test.df01 <- test.df_ %>% filter(Nbloom1==0)
      test.df11 <- test.df_ %>% filter(Nbloom1==1)
      pred.elast <- predict(elast_, test.df_, type="prob")[,2]
      pred.elast_01 <- predict(elast01, test.df01, type="prob")[,2]
      if(nrow(test.df11) > 0) {
        pred.elast_11 <- predict(elast11, test.df11, type="prob")[,2]
      } else {
        pred.elast_11 <- numeric(0)
      }
      
      cv_pred[[k]] <- full_join(
        cv_test.df %>%
          mutate(glmElast_mnpr=pred.elast),
        tibble(obsid=c(filter(cv_test.df, Nbloom1==0)$obsid, filter(cv_test.df, Nbloom1==1)$obsid),
               glmElast_split_mnpr=c(pred.elast_01, pred.elast_11))
      ) %>%
        mutate(covarSet=i.name,
               species=target)
    }
    cv_pred %>% do.call('rbind', .) %>%
      write_csv(glue("{f.prefix}CV_glm_{f.suffix}.csv"))
    
    
    cat("Finished", target, "\n")
  }
  
}




