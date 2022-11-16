# HABReports Bayesian modelling
# Bayesian HAB model initial fits
# Tim Szewczyk


# This script fits models with data from 2013-2017 using horseshoe priors



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "lme4")
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
           "NlnWt1", "NlnWt2",
           "NlnRAvg1", "NlnRAvg2"), 
         ":ydayCos:ydaySin")
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
    f.prefix <- glue("out{sep}full{sep}")
    f.suffix <- glue("{i.name}_{target}")
    
    target.df <- read_csv(glue("{f.prefix}dataset_{i.name}_{target}.csv"))
    
    train.df <- target.df %>% filter(year <= 2019)
    test.df <- target.df %>% filter(year > 2019)
    bloomThresh <- max((!target.df$Nbloom)*target.df$NcatNum) # 1:4, maximum considered 'No bloom'
    
    
    
    
    # Full fit for predictions
    form <- glue("Nbloom ~ {paste0(covar_int, collapse=' + ')} + lon_sc*lat_sc")
    out <- glm(form, data=train.df, family="binomial")
    out.01 <- glm(form, data=train.df %>% filter(Nbloom1==0), family="binomial")
    out.11 <- glm(form, data=train.df %>% filter(Nbloom1==1), family="binomial")
    
    saveRDS(out, glue("{f.prefix}glm_{f.suffix}.rds"))
    saveRDS(out.01, glue("{f.prefix}glm01_{f.suffix}.rds"))
    saveRDS(out.11, glue("{f.prefix}glm11_{f.suffix}.rds"))
    
    # Fitted
    fit.df <- full_join(
      train.df %>%
        mutate(glm_mnpr=fitted(out)),
      tibble(obsid=c(filter(train.df, Nbloom1==0)$obsid, filter(train.df, Nbloom1==1)$obsid),
             glm_split_mnpr=c(fitted(out.01), fitted(out.11)))
    ) %>%
      mutate(covarSet=i.name,
             species=target)
    write_csv(fit.df, glue("{f.prefix}fit_glm_{f.suffix}.csv"))
    
    # OOS predictions
    test.df01 <- test.df %>% filter(Nbloom1==0) %>% droplevels
    test.df11 <- test.df %>% filter(Nbloom1==1) %>% droplevels
    pred.glm <- predict(out, newdata=test.df, type="response", allow.new.levels=T)
    pred.glm_01 <- predict(out.01, newdata=test.df01, type="response", allow.new.levels=T)
    if(nrow(test.df11) > 0) {
      pred.glm_11 <- predict(out.11, newdata=test.df11, type="response", allow.new.levels=T)
    } else {
      pred.glm_11 <- numeric(0)
    }
    
    pred.df <- full_join(
      test.df %>%
        mutate(glm_mnpr=pred.glm),
      tibble(obsid=c(filter(test.df, Nbloom1==0)$obsid, filter(test.df, Nbloom1==1)$obsid),
             glm_split_mnpr=c(pred.glm_01, pred.glm_11))
    ) %>%
      mutate(covarSet=i.name,
             species=target)
    
    write_csv(pred.df, glue("{f.prefix}pred_glm_{f.suffix}.csv"))
    
    
    # Cross-validation by year
    yrCV <- unique(train.df$year)
    cv_pred <- map(yrCV, ~NULL)
    for(k in 1:length(yrCV)) {
      yr <- yrCV[k]
      cv_train.df <- train.df %>% filter(year != yr)
      cv_test.df <- train.df %>% filter(year == yr)
      
      form <- glue("Nbloom ~ {paste0(covar_int, collapse=' + ')} + lon_sc*lat_sc")
      out <- glm(form, data=cv_train.df, family="binomial")
      out.01 <- glm(form, data=cv_train.df %>% filter(Nbloom1==0), family="binomial")
      out.11 <- glm(form, data=cv_train.df %>% filter(Nbloom1==1), family="binomial")
      
      # Cross-validation predictions
      cv_test.df01 <- cv_test.df %>% filter(Nbloom1==0) %>% droplevels
      cv_test.df11 <- cv_test.df %>% filter(Nbloom1==1) %>% droplevels
      pred.glm <- predict(out, newdata=cv_test.df, type="response", allow.new.levels=T)
      pred.glm_01 <- predict(out.01, newdata=cv_test.df01, type="response", allow.new.levels=T)
      if(nrow(cv_test.df11) > 0) {
        pred.glm_11 <- predict(out.11, newdata=cv_test.df11, type="response", allow.new.levels=T)
      } else {
        pred.glm_11 <- numeric(0)
      }
      
      cv_pred[[k]] <- full_join(
        cv_test.df %>%
          mutate(glm_mnpr=pred.glm),
        tibble(obsid=c(filter(cv_test.df, Nbloom1==0)$obsid, filter(cv_test.df, Nbloom1==1)$obsid),
               glm_split_mnpr=c(pred.glm_01, pred.glm_11))
      ) %>%
        mutate(covarSet=i.name,
               species=target)
    }
    cv_pred %>% do.call('rbind', .) %>%
      write_csv(glue("{f.prefix}CV_glm_{f.suffix}.csv"))
    
    
    cat("Finished", target, "\n")
  }
  
}




