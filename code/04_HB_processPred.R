# HABReports Bayesian modelling
# Bayesian HAB model weekly forecast 
# Tim Szewczyk


# This script fits models with weekly data from 2018-2019 to assess forecasting performance



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "LaplacesDemon", "brms", "bayesplot", "WeStCOMS")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
walk(dir("code", "000_fn", full.names=T), source)

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

tl.col <- c("0_green"="green3", "1_yellow"="gold1", "2_orange"="orange", "3_red"="red")



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





sampling.df <- sampling.df %>% 
  mutate(wk=floor(as.numeric(date - min(date))/7),
         month=month(date))
for(sp in 1:length(species)) {
  target <- species[sp]
  target.tf <- thresh.df %>% filter(hab_parameter==target)
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
  
  ordPred.df <- dir(glue("out{sep}weekFitFull"), glue("pred_ord_{target}"), full.names=T) %>%
    map_dfr(~read_csv(.x, show_col_types=F))
  hufPred.df <- dir(glue("out{sep}weekFitFull"), glue("pred_huf_{target}"), full.names=T) %>%
    map_dfr(~read_csv(.x, show_col_types=F))
  
  if(nrow(hufPred.df)>0) {
    pred.df <- target.df %>% select(obs.id, site.id, date, N.catF, N.ln) %>%
      right_join(., full_join(ordPred.df, hufPred.df)) %>%
      mutate(sp=target)
    write_csv(pred.df, glue("out{sep}pred_{target}.csv"))
  } else if(nrow(ordPred.df)>0) {
    pred.df <- target.df %>% select(obs.id, site.id, date, N.catF, N.ln) %>%
      right_join(., ordPred.df) %>%
      mutate(sp=target)
    write_csv(pred.df, glue("out{sep}predFull_{target}.csv"))
  }
  
  

  
}



