# HABReports Bayesian modelling
# Bayesian HAB model weekly forecast 
# Tim Szewczyk


# This script fits models with weekly data from 2018-2019 to assess forecasting performance



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "LaplacesDemon", "brms", "WeStCOMS")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
walk(dir("code", "000_fn", full.names=T), source)

huf <- F
ord <- T

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
  ungroup

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
# >0.85 cor(res,time) = temp, salinity
# >0.85 cor(res) = short_wave, km, precip, wind, windDir
# >0.85 cor(time) = water
# waterDir is more variable across all 4
covars <- c("temp_L_wk", "salinity_L_wk",
            "short_wave_L_0", "short_wave_L_wk",
            "km_L_0", "km_L_wk",
            "precip_L_0", "precip_L_wk",
            "wind_L_0", "wind_L_wk",
            "windDir_L_0", "windDir_L_wk",
            "water_L_wk", "water_R_wk",
            "waterDir_L_0", "waterDir_R_0", "waterDir_L_wk", "waterDir_R_wk")

predictors_main <- "
  ydayCos + ydaySin + 
                 
  short_wave_L_wk:ydayCos +
  temp_L_wk:ydayCos +
  precip_L_wk:ydayCos +
  salinity_L_wk:ydayCos +
  water_L_wk:ydayCos +
  water_R_wk:ydayCos +
  wind_L_wk:ydayCos +
  windDir_L_wk:ydayCos +
  waterDir_L_wk:ydayCos +
  waterDir_R_wk:ydayCos +
  fetch +
                 
  N.PA_1:ydayCos +
  N.PA_2:ydayCos +

  N.lnWt_1:ydayCos + 
  N.lnWt_2:ydayCos + 

  N.PA_1:N.PA_2 +

  wind_L_wk:fetch +
  water_L_wk:fetch +
  windDir_L_wk:fetch +
  waterDir_L_wk:fetch +
  waterDir_R_wk:fetch +
                 
  (1|site.id)
"
predictors_hu <- "
  ydayCos + ydaySin +
  N.PA_1:ydayCos +
  N.PA_1:ydaySin +
  N.PA_1:N.PA_2 +
  fetch +
  wind_L_wk:fetch +
  water_L_wk:fetch
"

form_huf <- bf(glue("N.ln ~ {predictors_main}"),
                      glue("hu ~ {predictors_hu}"))
form_ordinal <- bf(glue("N.catNum | thres(3) ~ {predictors_main}"))






# weekly updates ----------------------------------------------------------

sampling.df <- sampling.df %>% 
  mutate(wk=floor(as.numeric(date - min(date))/7),
         month=month(date))
weeks <- sort(unique(filter(sampling.df, year(date)>2017 & between(month, 3, 11))$wk))

for(sp in 1:length(species)) {
  target <- species[sp]
  target.tf <- thresh.df %>% filter(hab_parameter==target)
  out.ord <- out.huf <- vector("list", length(weeks))
  for(i in 99:(length(weeks)-1)) {
    week_fit <- weeks[i]
    week_pred <- weeks[i+1]
    
    target.df <- sampling.df %>% 
      rename(N=!!target) %>%
      select(obs.id, site.id, date, wk, hour, grid, lon, lat, fetch, bearing, N) %>%
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
      mutate(N.catF=factor(N.cat, levels=unique(target.tf$tl), ordered=T),
             N.catNum=as.numeric(N.catF)) %>%
      arrange(site.id, date) %>%
      group_by(site.id) %>%
      multijetlag(N.ln, N.PA, N.cat, N.catF, date, n=2) %>%
      ungroup %>%
      mutate(across(starts_with("date_"), ~as.numeric(date-.x)),
             # I don't love this since small if N.ln_x is small OR date_x is large
             N.lnWt_1=N.ln_1/date_1,
             N.lnWt_2=N.ln_2/date_2) %>%
      left_join(hydro.df) %>%
      mutate(across(contains("Dir_"), ~cos(.x-bearing))) %>%
      mutate(across(one_of(grep("Dir", covars, invert=T, value=T)), CenterScale)) %>%
      arrange(site.id, date) %>%
      filter(complete.cases(.)) %>%
      select(site.id, lon, lat, date, year, wk, obs.id, fetch, bearing,
             starts_with("N"), starts_with("date_"), starts_with("yday"),
             one_of(covars)) 
    
    
    if(i == 1) {
      train.df <- target.df %>% filter(wk <= week_fit)
      if(ord) {
        out_ord_im1 <- readRDS(glue("out{sep}initFit{sep}ordinal_{target}.rds")) 
      }
      if(huf) {
        out_huf_im1 <- readRDS(glue("out{sep}initFit{sep}huf_{target}.rds")) 
      }
    } else {
      train.df <- target.df %>% filter(wk == week_fit)
      if(ord) {
        out_ord_im1 <- readRDS(glue("out{sep}weekFit{sep}ordinal_{target}_{weeks[i-1]}.rds"))  
      }
      if(huf) {
        out_huf_im1 <- readRDS(glue("out{sep}weekFit{sep}huf_{target}_{weeks[i-1]}.rds"))  
      }
    }
    test.df <- target.df %>% filter(wk == week_pred)
    
    # ordinal priors
    if(ord) {
      b_ord_im1 <- as_draws_df(out_ord_im1, variable="b_", regex=T) %>%
        pivot_longer(cols=starts_with("b_"), names_to="param", values_to="val") %>%
        group_by(param) %>%
        summarise(mn=mean(val), sd=sd(val)) %>%
        mutate(coef=str_remove(param, "b_"))
      eff_ord_im1 <- b_ord_im1 %>% filter(!grepl("Intercept", param))
      int_ord_im1 <- b_ord_im1 %>% filter(grepl("Intercept", param)) %>%
        mutate(coef=str_remove(str_remove(str_remove(coef, "Intercept"), "\\["), "]"))
      sd_im1 <- summary(out_ord_im1)$random$site.id
      prior_ord.ls <- vector("list", nrow(b_ord_im1)+1)
      for(j in 1:nrow(int_ord_im1)) {
        prior_ord.ls[[j]] <- prior_string(glue("normal({int_ord_im1$mn[j]},{int_ord_im1$sd[j]})"),
                                          class="Intercept", coef=int_ord_im1$coef[j])
      }
      for(j in 1:nrow(eff_ord_im1)) {
        prior_ord.ls[[j+nrow(int_ord_im1)]] <- prior_string(glue("normal({eff_ord_im1$mn[j]},{eff_ord_im1$sd[j]})"),
                                                            class="b", coef=eff_ord_im1$coef[j])
      }
      prior_ord.ls[[nrow(b_ord_im1)+1]] <- prior_string(glue("normal({sd_im1$Estimate},{sqrt(sd_im1$Est.Error)})"),
                                                        class="sd")
      prior_ord <- do.call(rbind, prior_ord.ls)
    }
    
    # hurdle frechet priors
    if(huf) {
      b_huf_im1 <- as_draws_df(out_huf_im1, variable="b_", regex=T) %>%
        pivot_longer(cols=starts_with("b_"), names_to="param", values_to="val") %>%
        group_by(param) %>%
        summarise(mn=mean(val), sd=sd(val)) %>%
        mutate(coef=str_remove(param, "b_"))
      fr_int_huf_im1 <- b_huf_im1 %>% filter(grepl("Intercept", param) & !grepl("_hu_", param)) %>%
        mutate(coef="")
      hu_int_huf_im1 <- b_huf_im1 %>% filter(grepl("Intercept", param) & grepl("_hu_", param)) %>%
        mutate(coef="")
      fr_eff_huf_im1 <- b_huf_im1 %>% filter(!grepl("Intercept", param) & !grepl("_hu_", param))
      hu_eff_huf_im1 <- b_huf_im1 %>% filter(!grepl("Intercept", param) & grepl("_hu_", param)) %>%
        mutate(coef=str_remove(coef, "hu_"))
      sd_im1 <- summary(out_huf_im1)$random$site.id
      nu_im1 <- summary(out_huf_im1)$spec_pars
      prior_huf.ls <- vector("list", nrow(b_huf_im1))
      prior_huf.ls[[1]] <- prior_string(glue("normal({fr_int_huf_im1$mn[1]},{fr_int_huf_im1$sd[1]})"),
                                        class="Intercept", coef="")
      prior_huf.ls[[2]] <- prior_string(glue("normal({hu_int_huf_im1$mn[1]},{hu_int_huf_im1$sd[1]})"),
                                        class="Intercept", coef="", dpar="hu")
      k <- 3
      for(j in 1:nrow(fr_eff_huf_im1)) {
        prior_huf.ls[[k]] <- prior_string(glue("normal({fr_eff_huf_im1$mn[j]},{fr_eff_huf_im1$sd[j]})"),
                                          class="b", coef=fr_eff_huf_im1$coef[j])
        k <- k+1
      }
      for(j in 1:nrow(hu_eff_huf_im1)) {
        prior_huf.ls[[k]] <- prior_string(glue("normal({hu_eff_huf_im1$mn[j]},{hu_eff_huf_im1$sd[j]})"),
                                          class="b", coef=hu_eff_huf_im1$coef[j], dpar="hu")
        k <- k+1
      }
      prior_huf.ls[[k]] <- prior_string(glue("normal({sd_im1$Estimate},{sqrt(sd_im1$Est.Error)})"),
                                        class="sd")
      prior_huf.ls[[k+1]] <- prior_string(glue("normal({nu_im1$Estimate},{sqrt(nu_im1$Est.Error)})"),
                                          class="nu")
      
      prior_huf <- do.call(rbind, prior_huf.ls)
    }
    
    
    date_range <- paste(as.character(range(train.df$date)), collapse=" to ")
    cat("Starting", i, "= week", as.character(week_fit), "=", date_range, "\n")
    cat("  train:", nrow(train.df), "rows with", n_distinct(train.df$site.id), "sites",
        "    ", table(train.df$N.catF), "\n")
    if(i == 1) {
      if(ord) {
        out.ord[[i]] <- brm(form_ordinal, train.df, cumulative("probit"),
                            prior=prior_ord, chains=4, cores=4, 
                            iter=3000, warmup=2000, refresh=500, inits="0",
                            control=list(adapt_delta=0.95, max_treedepth=20),
                            file=glue("out{sep}weekFit{sep}ordinal_{target}_{week_fit}"))
      }
      if(huf) {
        stanvars <- stanvar(scode=readr::read_file(glue("models{sep}hurdle_frechet_fn.stan")),
                            block="functions")
        out.huf[[i]] <- brm(form_huf, train.df, hurdle_frechet,
                            prior=prior_huf, chains=4, cores=4, 
                            iter=3000, warmup=2000, refresh=500, inits="0",
                            control=list(adapt_delta=0.95, max_treedepth=20),
                            stanvars=stanvars,
                            file=glue("out{sep}weekFit{sep}huf_{target}_{week_fit}")) 
      }
    } else {
      if(ord) {
        out.ord[[i]] <- update(out_ord_im1, newdata=train.df, cores=4, 
                               prior=prior_ord, file=glue("out{sep}weekFit{sep}ordinal_{target}_{week_fit}")) 
      }
      if(huf) {
        out.huf[[i]] <- update(out_huf_im1, newdata=train.df, cores=4,
                               prior=prior_huf, file=glue("out{sep}weekFit{sep}huf_{target}_{week_fit}")) 
      }
    }
    
    pred.df <- test.df %>% select(obs.id)
    if(ord) {
      pred.ord <- posterior_epred(out.ord[[i]], newdata=test.df, allow_new_levels=T)
      cat.iter <- apply(pred.ord, 1:2, which.max)
      pred.df <- pred.df %>%
        mutate(ord_mnpr0=colMeans(pred.ord[,,1,drop=F]),
               ord_mnpr1=colMeans(pred.ord[,,2,drop=F]),
               ord_mnpr2=colMeans(pred.ord[,,3,drop=F]),
               ord_mnpr3=colMeans(pred.ord[,,4,drop=F]),
               ord_prmax0=colMeans(cat.iter==1),
               ord_prmax1=colMeans(cat.iter==2),
               ord_prmax2=colMeans(cat.iter==3),
               ord_prmax3=colMeans(cat.iter==4),
               ord_mnCat=colMeans(cat.iter),
               ord_wtmnCat=ord_mnpr0 + 2*ord_mnpr1 + 3*ord_mnpr2 + 4*ord_mnpr3)
    }
    if(huf) {
      pred.hu.huf <- posterior_epred(out.huf[[i]], newdata=test.df, allow_new_levels=T, dpar="hu")
      pred.mu.huf <- posterior_epred(out.huf[[i]], newdata=test.df, allow_new_levels=T, dpar="mu")
      pred.huf <- (1-rbinom(prod(dim(pred.hu.huf)), 1, pred.hu.huf)) * pred.mu.huf
      pred.huf.exp <- (1-rbinom(prod(dim(pred.hu.huf)), 1, pred.hu.huf)) * exp(pred.mu.huf)
      cat.iter <- apply(pred.huf.exp, 1:2, function(x) (1:4)[max(which(x >= target.tf$min_ge))])
      pred.df <- pred.df %>%
        mutate(huf_mn=colMeans(pred.huf),
               huf_mn_exp=colMeans(pred.huf.exp),
               huf_md=apply(pred.huf, 2, median),
               huf_md_exp=apply(pred.huf.exp, 2, median),
               huf_mnpr0=apply(pred.huf.exp, 2,
                           function(x) mean(x>=target.tf$min_ge[1] &
                                              x<=target.tf$min_ge[2])),
               huf_mnpr1=apply(pred.huf.exp, 2,
                           function(x) mean(x>=target.tf$min_ge[2] &
                                              x<=target.tf$min_ge[3])),
               huf_mnpr2=apply(pred.huf.exp, 2,
                           function(x) mean(x>=target.tf$min_ge[3] &
                                              x<=target.tf$min_ge[4])),
               huf_mnpr3=apply(pred.huf.exp, 2,
                           function(x) mean(x>=target.tf$min_ge[4])),
               huf_prmax0=colMeans(cat.iter==1),
               huf_prmax1=colMeans(cat.iter==2),
               huf_prmax2=colMeans(cat.iter==3),
               huf_prmax3=colMeans(cat.iter==4),
               huf_mnCat=colMeans(cat.iter),
               huf_wtmnCat=huf_mnpr0 + 2*huf_mnpr1 + 3*huf_mnpr2 + 4*huf_mnpr3)
    }
    if(ord & huf) {
      write_csv(pred.df, glue("out{sep}weekFit{sep}pred_{target}_{week_pred}.csv"))
    } else if(ord) {
      write_csv(pred.df, glue("out{sep}weekFit{sep}pred_ord_{target}_{week_pred}.csv"))
    } else {
      write_csv(pred.df, glue("out{sep}weekFit{sep}pred_huf_{target}_{week_pred}.csv"))
    }
    
    cat("  Fitted model \n")
  }
}



