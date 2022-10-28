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
rebalance_pr0 <- c(NA, 0.4)[1]
ctrl <- list(adapt_delta=0.95, max_treedepth=20)
chains <- 2
iter <- 100
warmup <- iter/2
refresh <- 1

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

sp.i <- read_csv("data/sp_i.csv")




# covariates --------------------------------------------------------------

species <- sp.i$abbr

data.df <- dir("out/full", "dataset.*csv", full.names=T) %>% 
  grep("rebal", ., value=T, invert=ifelse(is.na(rebalance_pr0), T, F)) %>%
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
           paste0("Nbloom1", species), paste0("Nbloom2", species),
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
  paste0("Nbloom1", species), paste0("Nbloom2", species),
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

f.prefix <- glue("out{sep}full{sep}")
f.suffix <- glue("{i.name}",
                 "{ifelse(is.na(rebalance_pr0),'',paste0('_rebal',rebalance_pr0*100))}")

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

# Interaction priors
priors.ord <- map(
  unique(mv_i$ord),
  ~c(prior_string("horseshoe(3,par_ratio=0.2)", class="b", resp=.x),
     prior_string("normal(0,1)", class="Intercept", resp=.x),
     prior_string("normal(0,0.5)", class="sd", lb=0, resp=.x))) %>% 
  do.call('c',. )
priors.bern <- map(
  unique(mv_i$bern),
  ~c(prior_string("horseshoe(3,par_ratio=0.2)", class="b", resp=.x),
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
         prior_string("double_exponential(0,0.1)", class="sds", nlpar="bIntercept", lb=0, resp=.x))) %>%
    do.call('c', .),
  map(1:nrow(mv_i),
      ~c(prior_string("beta(0.2,1)", nlpar=mv_i$nl_p[.x], resp=mv_i$ord[.x], lb=0, ub=1),
         prior_string("normal(0,1)", class="b", nlpar=mv_i$nl_b[.x], resp=mv_i$ord[.x]),
         prior_string("normal(0,.5)", class="sd", nlpar=mv_i$nl_p[.x], resp=mv_i$ord[.x], lb=0),
         prior_string("normal(0,.5)", class="sd", nlpar=mv_i$nl_b[.x], resp=mv_i$ord[.x], lb=0),
         prior_string("double_exponential(0,0.1)", class="sds", nlpar=mv_i$nl_b[.x], resp=mv_i$ord[.x], lb=0))) %>%
    do.call('c', .)
)
priors.bernP <- c(
  map(unique(mv_i$bern),
      ~c(prior_string("normal(0,1)", class="b", nlpar="bIntercept", resp=.x),
         prior_string("normal(0,.5)", class="sd", nlpar="bIntercept", lb=0, resp=.x),
         prior_string("double_exponential(0,0.1)", class="sds", nlpar="bIntercept", lb=0, resp=.x))) %>%
    do.call('c', .),
  map(1:nrow(mv_i),
      ~c(prior_string("beta(0.2,1)", nlpar=mv_i$nl_p[.x], resp=mv_i$bern[.x], lb=0, ub=1),
         prior_string("normal(0,1)", class="b", nlpar=mv_i$nl_b[.x], resp=mv_i$bern[.x]),
         prior_string("normal(0,.5)", class="sd", nlpar=mv_i$nl_p[.x], resp=mv_i$bern[.x], lb=0),
         prior_string("normal(0,.5)", class="sd", nlpar=mv_i$nl_b[.x], resp=mv_i$bern[.x], lb=0),
         prior_string("double_exponential(0,0.1)", class="sds", nlpar=mv_i$nl_b[.x], resp=mv_i$bern[.x], lb=0))) %>%
    do.call('c', .)
)



# initial fit -------------------------------------------------------------

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


# loop through species to load target.df, get fits and predictions with resp=sp

# initial fit -------------------------------------------------------------
# 
# for(sp in 1) {
#   bloomThresh <- max((!target.df$Nbloom)*target.df$NcatNum) # 1:4, maximum considered 'No bloom'
#   
#   # Full fit for predictions
#   
#   
#   # Fitted
#   fit.ord <- posterior_epred(out.ord)
#   fit.ordP <- posterior_epred(out.ordP)
#   fit.bern01 <- posterior_epred(out.bern01) 
#   fit.bernP01 <- posterior_epred(out.bernP01) 
#   fit.bern11 <- posterior_epred(out.bern11) 
#   fit.bernP11 <- posterior_epred(out.bernP11) 
#   
#   fit.df <- full_join(
#     train.df %>%
#       mutate(ord_mnpr=calc_ord_mnpr(fit.ord, bloomThresh),
#              ordP_mnpr=calc_ord_mnpr(fit.ordP, bloomThresh)),
#     tibble(obsid=c(filter(train.df, Nbloom1==0)$obsid, filter(train.df, Nbloom1==1)$obsid),
#            bern_mnpr=c(colMeans(fit.bern01), colMeans(fit.bern11)),
#            bernP_mnpr=c(colMeans(fit.bernP01), colMeans(fit.bernP11)))
#   ) %>%
#     mutate(covarSet=i.name)
#   write_csv(fit.df, glue("{f.prefix}fit_HBuv_{f.suffix}.csv"))
#   
#   # OOS predictions
#   pred.ord <- posterior_epred(out.ord, newdata=test.df, allow_new_levels=T)
#   pred.ordP <- posterior_epred(out.ordP, newdata=test.df, allow_new_levels=T)
#   pred.bern01 <- colMeans(posterior_epred(out.bern01, newdata=test.df01, allow_new_levels=T))
#   pred.bernP01 <- colMeans(posterior_epred(out.bernP01, newdata=test.df01, allow_new_levels=T))
#   if(nrow(test.df11) > 0) {
#     pred.bern11 <- colMeans(posterior_epred(out.bern11, newdata=test.df11, allow_new_levels=T))
#     pred.bernP11 <- colMeans(posterior_epred(out.bernP11, newdata=test.df11, allow_new_levels=T))
#   } else {
#     pred.bern11 <- numeric(0)
#     pred.bernP11 <- numeric(0)
#   }
#   
#   pred.df <- full_join(
#     test.df %>%
#       mutate(ord_mnpr=calc_ord_mnpr(pred.ord, bloomThresh),
#              ordP_mnpr=calc_ord_mnpr(pred.ordP, bloomThresh)),
#     tibble(obsid=c(filter(test.df, Nbloom1==0)$obsid, filter(test.df, Nbloom1==1)$obsid),
#            bern_mnpr=c(pred.bern01, pred.bern11),
#            bernP_mnpr=c(pred.bernP01, pred.bernP11))
#   ) %>%
#     mutate(covarSet=i.name)
#   
#   write_csv(pred.df, glue("{f.prefix}pred_HBuv_{f.suffix}.csv"))
#   
#   # Cross-validation by year
#   yrCV <- unique(train.df$year)
#   cv_pred <- map(yrCV, ~NULL)
#   for(k in 1:length(yrCV)) {
#     yr <- yrCV[k]
#     cv_train.df <- train.df %>% filter(year != yr)
#     cv_test.df <- train.df %>% filter(year == yr)
#     
#     # Formulas with interactions: errors if missing NcatF levels
#     form_ord <- makeFormula(cv_train.df, covar_int, "NcatNum | thres(3)")
#     form_01 <- makeFormula(filter(cv_train.df, Nbloom1==0), covar_int, "Nbloom")
#     form_11 <- makeFormula(filter(cv_train.df, Nbloom1==1), covar_int, "Nbloom")
#     
#     # Smoother formulas
#     form_ordP <- makeFormula(cv_train.df, covar_s, "NcatNum | thres(3)", covar_date,
#                              flist.P, list(b=s_b, p=s_p))
#     form_bernP <- makeFormula(cv_train.df, covar_s, "Nbloom", covar_date,
#                               flist.P, list(b=s_b, p=s_p))
#     
#     cv.ord <- brm(form_ord, data=cv_train.df,
#                   family=cumulative("probit"), prior=priors, 
#                   iter=iter, warmup=warmup, refresh=refresh, init=0,
#                   control=ctrl, chains=chains, cores=chains,
#                   file=glue("{f.prefix}ord_CV{k}_{f.suffix}"))
#     cv.ordP <- brm(form_ordP, data=cv_train.df, 
#                    family=cumulative("probit"), prior=priors.P, 
#                    iter=iter, warmup=warmup, refresh=refresh, init=0,
#                    control=ctrl, chains=chains, cores=chains,
#                    file=glue("{f.prefix}ordP_CV{k}_{f.suffix}"))
#     cv.bern01 <- brm(form_01, data=cv_train.df %>% filter(Nbloom1==0),
#                      family=bernoulli("probit"), prior=priors, 
#                      iter=iter, warmup=warmup, refresh=refresh, init=0,
#                      control=ctrl, chains=chains, cores=chains,
#                      file=glue("{f.prefix}bern01_CV{k}_{f.suffix}"))
#     cv.bernP01 <- brm(form_bernP, data=cv_train.df %>% filter(Nbloom1==0), 
#                       family=bernoulli("probit"), prior=priors.P, 
#                       iter=iter, warmup=warmup, refresh=refresh, init=0,
#                       control=ctrl, chains=chains, cores=chains,
#                       file=glue("{f.prefix}bernP01_CV{k}_{f.suffix}"))
#     cv.bern11 <- brm(form_11, data=cv_train.df %>% filter(Nbloom1==1),
#                      family=bernoulli("probit"), prior=priors, 
#                      iter=iter, warmup=warmup, refresh=refresh, init=0,
#                      control=ctrl, chains=chains, cores=chains,
#                      file=glue("{f.prefix}bern11_CV{k}_{f.suffix}"))
#     cv.bernP11 <- brm(form_bernP, data=cv_train.df %>% filter(Nbloom1==1), 
#                       family=bernoulli("probit"), prior=priors.P, 
#                       iter=iter, warmup=warmup, refresh=refresh, init=0,
#                       control=ctrl, chains=chains, cores=chains,
#                       file=glue("{f.prefix}bernP11_CV{k}_{f.suffix}"))
#     
#     # Cross-validation predictions
#     cv_test.df01 <- cv_test.df %>% filter(Nbloom1==0) %>% droplevels
#     cv_test.df11 <- cv_test.df %>% filter(Nbloom1==1) %>% droplevels
#     pred.ord <- posterior_epred(cv.ord, newdata=cv_test.df, allow_new_levels=T)
#     pred.ordP <- posterior_epred(cv.ordP, newdata=cv_test.df, allow_new_levels=T)
#     pred.bern01 <- colMeans(posterior_epred(cv.bern01, newdata=cv_test.df01, allow_new_levels=T))
#     pred.bernP01 <- colMeans(posterior_epred(cv.bernP01, newdata=cv_test.df01, allow_new_levels=T))
#     if(nrow(cv_test.df11) > 0) {
#       pred.bern11 <- colMeans(posterior_epred(cv.bern11, newdata=cv_test.df11, allow_new_levels=T))
#       pred.bernP11 <- colMeans(posterior_epred(cv.bernP11, newdata=cv_test.df11, allow_new_levels=T))
#     } else {
#       pred.bern11 <- numeric(0)
#       pred.bernP11 <- numeric(0)
#     }
#     
#     cv_pred[[k]] <- full_join(
#       cv_test.df %>%
#         mutate(ord_mnpr=calc_ord_mnpr(pred.ord, bloomThresh),
#                ordP_mnpr=calc_ord_mnpr(pred.ordP, bloomThresh)),
#       tibble(obsid=c(filter(cv_test.df, Nbloom1==0)$obsid, filter(cv_test.df, Nbloom1==1)$obsid),
#              bern_mnpr=c(pred.bern01, pred.bern11),
#              bernP_mnpr=c(pred.bernP01, pred.bernP11))
#     ) %>%
#       mutate(covarSet=i.name)
#   }
#   cv_pred %>% do.call('rbind', .) %>%
#     write_csv(glue("{f.prefix}CV_HBuv_{f.suffix}.csv"))
#   
#   cat("Finished", target, "\n")
# }
# 
# 
# 
# 
