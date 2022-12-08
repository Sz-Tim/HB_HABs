# HABReports Bayesian modelling
# Bayesian HAB model initial fits
# Tim Szewczyk


# This script fits models with data from 2013-2017 using horseshoe priors



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "brms")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
walk(dir("code", "*00_fn", full.names=T), source)

# Model details
ctrl <- list(adapt_delta=0.95, max_treedepth=20)
chains <- 4
iter <- 2000
warmup <- iter/2
refresh <- 20
prior_strength <- c(1,2,3)[1]

# minch2:    2013-06-20 to 2019-07-02
# WeStCOMS2: 2019-04-01 to 2022-01-26
sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")

sp.i <- read_csv("data/sp_i.csv")




# covariates --------------------------------------------------------------

species <- sp.i$abbr

data.df <- dir("out/full", "dataset.*csv", full.names=T) %>% 
  map_dfr(~read_csv(.x, show_col_type=F) %>%
            mutate(month=lubridate::month(date),
                   week=lubridate::week(date)) %>%
            mutate(species=str_sub(str_split_fixed(.x, "_", 3)[,3], 1, -5))) %>%
  mutate(species=sp.i$abbr[match(species, sp.i$full)]) %>%
  pivot_wider(names_from="species", values_from=starts_with("N")) %>%
  rename_with(.cols=starts_with("N"), .fn=~str_remove_all(.x, "_"))

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
           # paste0("Nbloom1", species), paste0("Nbloom2", species),
           paste0("NlnWt1", species), paste0("NlnWt2", species),
           paste0("NlnRAvg2", species), paste0("NlnRAvg2", species)), 
         ":ydayCos:ydaySin"),
  paste0(c("windVel", "waterVelL", "waterVelR", "windLwk", "waterLwk", "waterRwk"), 
         ":fetch")
)
covar_s.all <- c(
  "tempLwk", "salinityLwk", "shortwaveLwk", "kmLwk", "precipLwk",
  "tempStrat20mLwk", "tempStrat20mRwk",
  "windVel", "waterVelL", "waterVelR", "windLwk", "waterLwk", "waterRwk",
  "fetch",# "influxwk", 
  "attnwk", "chlwk", "dinowk", "o2wk", "phwk", "po4wk",
  # paste0("Nbloom1", species), paste0("Nbloom2", species),
  paste0("NlnWt1", species), paste0("NlnWt2", species),
  paste0("NlnRAvg1", species), paste0("NlnRAvg2", species)
)


covariate_sets <- list(
  # null="NA",
  # date="^yday",
  # autoreg="^N|mo",
  # external="fetch|influx|water|wind",
  # local="temp|salinity|shortwave|precip|km",
  # cprn="attn|chl|dino|o2|ph|po4",
  test="^NlnRAvg1|temp|fetch",
  all="."
)



# for(i in length(covariate_sets)) {
i <- 1
i.name <- names(covariate_sets)[i]
i.covs <- covariate_sets[[i]]
covars <- covars.all[grepl(i.covs, str_remove_all(covars.all, "\\.|_"))]
covar_int <- grep(i.covs, covar_int.all, value=T)
covar_s <- grep(i.covs, covar_s.all, value=T)
covar_date <- NULL


# Smoothers
s_b <- glue("b{covar_s}") 
s_p <- glue("p{covar_s}")
mv_i <- expand_grid(sp=species, 
                    covar=covar_s) %>%
  mutate(ordThres=glue("NcatNum{sp} | thres(3)"),
         ord=glue("NcatNum{sp}"),
         bern=glue("Nbloom{sp}"),
         nl_b=glue("b{covar}"),
         nl_p=glue("p{covar}"))

# Prior strengths
pri.var <- switch(prior_strength,
                  "1"=list(hs1=1, hs2=0.4, b=0.5, de=0.5, i="1-loose"),
                  "2"=list(hs1=3, hs2=0.2, b=0.2, de=0.1, i="2-medium"),
                  "3"=list(hs1=5, hs2=0.1, b=0.05, de=0.02, i="3-tight"))

# Interaction priors
priors.ord <- map(
  unique(mv_i$ord),
  ~c(prior_string(glue("horseshoe({pri.var$hs1}, par_ratio={pri.var$hs2})"), class="b", resp=.x),
     prior_string("normal(0,1)", class="Intercept", resp=.x),
     prior_string("normal(0,0.5)", class="sd", lb=0, resp=.x))) %>% 
  do.call('c',. )
priors.bern <- map(
  unique(mv_i$bern),
  ~c(prior_string(glue("horseshoe({pri.var$hs1}, par_ratio={pri.var$hs2})"), class="b", resp=.x),
     prior_string("normal(0,1)", class="Intercept", resp=.x),
     prior_string("normal(0,0.5)", class="sd", lb=0, resp=.x))) %>% 
  do.call('c',. )

# Indicator variable: p * b * X
flist.P <- c(
  map(s_b, ~as.formula(paste0(.x, "~ s(ydayCos,ydaySin) + (1|siteid)"))),
  map(s_p, ~as.formula(paste0(.x, "~ 1 + (1|siteid)"))),
  map(1, ~"bIntercept ~ 1 + s(ydayCos,ydaySin) + (1|siteid)")
)
priors.ordP <- c(
  map(unique(mv_i$ord), 
      ~c(prior_string("normal(0,1)", class="b", nlpar="bIntercept", resp=.x),
         prior_string("normal(0,.5)", class="sd", nlpar="bIntercept", lb=0, resp=.x),
         prior_string(glue("double_exponential(0,{pri.var$de})"), class="sds", nlpar="bIntercept", lb=0, resp=.x))) %>%
    do.call('c', .),
  map(1:nrow(mv_i),
      ~c(prior_string(glue("beta({pri.var$b},1)"), nlpar=mv_i$nl_p[.x], resp=mv_i$ord[.x], lb=0, ub=1),
         prior_string("normal(0,1)", class="b", nlpar=mv_i$nl_b[.x], resp=mv_i$ord[.x]),
         prior_string("normal(0,.5)", class="sd", nlpar=mv_i$nl_p[.x], resp=mv_i$ord[.x], lb=0),
         prior_string("normal(0,.5)", class="sd", nlpar=mv_i$nl_b[.x], resp=mv_i$ord[.x], lb=0),
         prior_string(glue("double_exponential(0,{pri.var$de})"), class="sds", nlpar=mv_i$nl_b[.x], resp=mv_i$ord[.x], lb=0))) %>%
    do.call('c', .)
)
priors.bernP <- c(
  map(unique(mv_i$bern),
      ~c(prior_string("normal(0,1)", class="b", nlpar="bIntercept", resp=.x),
         prior_string("normal(0,.5)", class="sd", nlpar="bIntercept", lb=0, resp=.x),
         prior_string(glue("double_exponential(0,{pri.var$de})"), class="sds", nlpar="bIntercept", lb=0, resp=.x))) %>%
    do.call('c', .),
  map(1:nrow(mv_i),
      ~c(prior_string(glue("beta({pri.var$b},1)"), nlpar=mv_i$nl_p[.x], resp=mv_i$bern[.x], lb=0, ub=1),
         prior_string("normal(0,1)", class="b", nlpar=mv_i$nl_b[.x], resp=mv_i$bern[.x]),
         prior_string("normal(0,.5)", class="sd", nlpar=mv_i$nl_p[.x], resp=mv_i$bern[.x], lb=0),
         prior_string("normal(0,.5)", class="sd", nlpar=mv_i$nl_b[.x], resp=mv_i$bern[.x], lb=0),
         prior_string(glue("double_exponential(0,{pri.var$de})"), class="sds", nlpar=mv_i$nl_b[.x], resp=mv_i$bern[.x], lb=0))) %>%
    do.call('c', .)
)



# initial fit -------------------------------------------------------------

f.prefix <- glue("out{sep}{pri.var$i}{sep}")
f.suffix <- glue("{i.name}")

train.df <- data.df %>% filter(year <= 2019)
test.df <- data.df %>% filter(year > 2019)

form_bern <- map(unique(mv_i$bern), 
                 ~makeFormula(train.df, covar_int, .x) + bernoulli("probit"))
form_ord <- map(unique(mv_i$ordThres), 
                ~makeFormula(train.df, covar_int, .x) + cumulative("probit")) 
form_bernP <- map(unique(mv_i$bern), 
                  ~makeFormula(train.df, covar_s, .x, NULL, 
                               flist.P, list(b=s_b, yday=NULL, p=s_p)) + 
                    bernoulli("probit"))
form_ordP <- map(unique(mv_i$ordThres), 
                 ~makeFormula(train.df, covar_s, .x, NULL, 
                              flist.P, list(b=s_b, yday=NULL, p=s_p)) + 
                   cumulative("probit"))  

out.bern <- brm(reduce(form_bern, `+`) + set_rescor(FALSE), 
                data=train.df, prior=priors.bern, 
                iter=iter, warmup=warmup, refresh=refresh, init=0,
                control=ctrl, chains=chains, cores=chains,
                file=glue("{f.prefix}mv_bern_{f.suffix}"))
out.ord <- brm(reduce(form_ord, `+`) + set_rescor(FALSE), 
               data=train.df, prior=priors.ord, 
               iter=iter, warmup=warmup, refresh=refresh, init=0,
               control=ctrl, chains=chains, cores=chains,
               file=glue("{f.prefix}mv_ord_{f.suffix}"))  
out.bernP <- brm(reduce(form_bernP, `+`) + set_rescor(FALSE), 
                 data=train.df, prior=priors.bernP, 
                 iter=iter, warmup=warmup, refresh=refresh, init=0,
                 control=ctrl, chains=chains, cores=chains,
                 file=glue("{f.prefix}mv_bernP_{f.suffix}"))
out.ordP <- brm(reduce(form_ordP, `+`) + set_rescor(FALSE), 
                data=train.df, prior=priors.ordP, 
                iter=iter, warmup=warmup, refresh=refresh, init=0,
                control=ctrl, chains=chains, cores=chains,
                file=glue("{f.prefix}mv_ordP_{f.suffix}"))


# Evaluation set up
data.sp <- dir("out/full", "dataset.*csv", full.names=T) %>% 
  map(read_csv)
train.sp <- map(data.sp, ~filter(.x, year <= 2019))
test.sp <- map(data.sp, ~filter(.x, year > 2019))
bloomThresh <- map(data.sp, ~max((!.x$Nbloom)*.x$NcatNum)) 


# Fitted
fits.ord <- map(unique(mv_i$ord), ~posterior_epred(out.ord, resp=.x))
fits.ordP <- map(unique(mv_i$ord), ~posterior_epred(out.ordP, resp=.x))
fits.bern <- map(unique(mv_i$bern), ~posterior_epred(out.bern, resp=.x))
fits.bernP <- map(unique(mv_i$bern), ~posterior_epred(out.bernP, resp=.x))

walk(1:nrow(sp.i), 
    ~train.sp[[.x]] %>%
      mutate(MVord_mnpr=calc_ord_mnpr(fits.ord[[.x]], bloomThresh[[.x]]),
             MVordP_mnpr=calc_ord_mnpr(fits.ordP[[.x]], bloomThresh[[.x]]),
             MVbern_mnpr=colMeans(fits.bern[[.x]]),
             MVbernP_mnpr=colMeans(fits.bernP[[.x]]),
             covarSet=i.name) %>%
      write_csv(glue("{f.prefix}fit_HBmv_all_{sp.i$full[.x]}.csv"))) 


# OOS prediction
preds.ord <- map(unique(mv_i$ord), ~posterior_epred(out.ord, resp=.x, newdata=test.df, allow_new_levels=T))
preds.ordP <- map(unique(mv_i$ord), ~posterior_epred(out.ordP, resp=.x, newdata=test.df, allow_new_levels=T))
preds.bern <- map(unique(mv_i$bern), ~posterior_epred(out.bern, resp=.x, newdata=test.df, allow_new_levels=T))
preds.bernP <- map(unique(mv_i$bern), ~posterior_epred(out.bernP, resp=.x, newdata=test.df, allow_new_levels=T))

walk(1:nrow(sp.i), 
     ~test.sp[[.x]] %>%
       mutate(MVord_mnpr=calc_ord_mnpr(preds.ord[[.x]], bloomThresh[[.x]]),
              MVordP_mnpr=calc_ord_mnpr(preds.ordP[[.x]], bloomThresh[[.x]]),
              MVbern_mnpr=colMeans(preds.bern[[.x]]),
              MVbernP_mnpr=colMeans(preds.bernP[[.x]]),
              covarSet=i.name) %>%
       write_csv(glue("{f.prefix}pred_HBmv_all_{sp.i$full[.x]}.csv"))) 



# cross validation --------------------------------------------------------
yrCV <- unique(train.df$year)
cv_pred <- map(yrCV, ~NULL)

for(k in 1:length(yrCV)) {
  yr <- yrCV[k]
  cv_train.df <- train.df %>% filter(year != yr)
  cv_test.df <- train.df %>% filter(year == yr)
  
  form_bern <- map(unique(mv_i$bern), 
                   ~makeFormula(cv_train.df, covar_int, .x) + bernoulli("probit"))
  form_ord <- map(unique(mv_i$ordThres), 
                  ~makeFormula(cv_train.df, covar_int, .x) + cumulative("probit")) 
  form_bernP <- map(unique(mv_i$bern), 
                    ~makeFormula(cv_train.df, covar_s, .x, NULL, 
                                 flist.P, list(b=s_b, yday=NULL, p=s_p)) + 
                      bernoulli("probit"))
  form_ordP <- map(unique(mv_i$ordThres), 
                   ~makeFormula(cv_train.df, covar_s, .x, NULL, 
                                flist.P, list(b=s_b, yday=NULL, p=s_p)) + 
                     cumulative("probit"))  
  
  cv.bern <- brm(reduce(form_bern, `+`) + set_rescor(FALSE), 
                  data=cv_train.df, prior=priors.bern, 
                  iter=iter, warmup=warmup, refresh=refresh, init=0,
                  control=ctrl, chains=chains, cores=chains,
                  file=glue("{f.prefix}mv_bern_CV{k}_{f.suffix}"))
  cv.ord <- brm(reduce(form_ord, `+`) + set_rescor(FALSE), 
                 data=cv_train.df, prior=priors.ord, 
                 iter=iter, warmup=warmup, refresh=refresh, init=0,
                 control=ctrl, chains=chains, cores=chains,
                 file=glue("{f.prefix}mv_ord_CV{k}_{f.suffix}"))  
  cv.bernP <- brm(reduce(form_bernP, `+`) + set_rescor(FALSE), 
                   data=cv_train.df, prior=priors.bernP, 
                   iter=iter, warmup=warmup, refresh=refresh, init=0,
                   control=ctrl, chains=chains, cores=chains,
                   file=glue("{f.prefix}mv_bernP_CV{k}_{f.suffix}"))
  cv.ordP <- brm(reduce(form_ordP, `+`) + set_rescor(FALSE), 
                  data=cv_train.df, prior=priors.ordP, 
                  iter=iter, warmup=warmup, refresh=refresh, init=0,
                  control=ctrl, chains=chains, cores=chains,
                  file=glue("{f.prefix}mv_ordP_CV{k}_{f.suffix}"))
  
  
  
  # Evaluation set up
  data.sp <- dir("out/full", "dataset.*csv", full.names=T) %>% 
    map(read_csv)
  cv_test.sp <- map(data.sp, ~filter(.x, year==yr))
  bloomThresh <- map(data.sp, ~max((!.x$Nbloom)*.x$NcatNum)) 
  
  preds.ord <- map(unique(mv_i$ord), ~posterior_epred(cv.ord, resp=.x, newdata=cv_test.df, allow_new_levels=T))
  preds.ordP <- map(unique(mv_i$ord), ~posterior_epred(cv.ordP, resp=.x, newdata=cv_test.df, allow_new_levels=T))
  preds.bern <- map(unique(mv_i$bern), ~posterior_epred(cv.bern, resp=.x, newdata=cv_test.df, allow_new_levels=T))
  preds.bernP <- map(unique(mv_i$bern), ~posterior_epred(cv.bernP, resp=.x, newdata=cv_test.df, allow_new_levels=T))
  
  cv_pred[[k]] <- map(
    1:nrow(sp.i), 
    ~cv_test.sp[[.x]] %>%
      mutate(MVord_mnpr=calc_ord_mnpr(preds.ord[[.x]], bloomThresh[[.x]]),
             MVordP_mnpr=calc_ord_mnpr(preds.ordP[[.x]], bloomThresh[[.x]]),
             MVbern_mnpr=colMeans(preds.bern[[.x]]),
             MVbernP_mnpr=colMeans(preds.bernP[[.x]]),
             covarSet=i.name)
    )
}

saveRDS(cv_pred, glue("{f.prefix}CV_PRED_HBmv_list.rds"))
for(i in 1:nrow(sp.i)) {
  map(cv_pred, ~.x[[i]]) %>% 
    do.call('rbind', .) %>%
    write_csv(glue("{f.prefix}CV_HBmv_all_{sp.i$full[i]}.csv"))
}
  
  









