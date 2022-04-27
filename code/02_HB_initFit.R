# HABReports Bayesian modelling
# Bayesian HAB model initial fits
# Tim Szewczyk


# This script fits models with data from 2013-2017 using horseshoe priors



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
  map(read_csv) %>%
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

predictors_main <- "ydayCos * ydaySin * short_wave_L_wk * temp_L_wk * 
  precip_L_wk * salinity_L_wk * water_L_wk * water_R_wk * wind_L_wk * 
  windDir_L_wk * waterDir_L_wk * waterDir_R_wk * fetch *
  N.PA_1 * N.PA_2 * N.lnWt_1 * N.lnWt_2 + (1|site.id)"

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







# initial fit -------------------------------------------------------------

out.huf <- out.ord <- vector("list", length(species))

for(sp in 1:length(species)) {
  target <- species[sp]
  target.tf <- thresh.df %>% filter(hab_parameter==target)
  
  train.df <- sampling.df %>% 
    rename(N=!!target) %>%
    select(obs.id, site.id, date, hour, grid, lon, lat, fetch, bearing, N) %>%
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
    full_join(hydro.df) %>%
    mutate(across(contains("Dir_"), ~cos(.x-bearing))) %>%
    mutate(across(one_of(grep("Dir", covars, invert=T, value=T)), CenterScale)) %>%
    arrange(site.id, date) %>%
    filter(complete.cases(.)) %>%
    select(site.id, lon, lat, date, year, obs.id, fetch, bearing,
           starts_with("N"), starts_with("date_"), starts_with("yday"),
           one_of(covars)) %>%
    filter(year <= 2017)
  
  # stanvars <- stanvar(scode=readr::read_file(glue("models{sep}hurdle_frechet_fn.stan")),
  #                     block="functions")
  # 
  # out.huf[[sp]] <- brm(form_huf, data=train.df,
  #                      family=hurdle_frechet(), stanvars=stanvars,
  #                      chains=4, cores=4, init="0",
  #                      iter=2000, warmup=1000, refresh=100,
  #                      control=list(adapt_delta=0.95, max_treedepth=20),
  #                      prior=c(prior(normal(2, 1), "Intercept", dpar="hu"),
  #                              prior(horseshoe(3, par_ratio=0.2), "b")),
  #                      save_model=glue("out{sep}initFit{sep}huf_{target}.stan"),
  #                      file=glue("out{sep}initFit{sep}huf_{target}"))
  
  out.ord[[sp]] <- brm(form_ordinal, data=train.df,
                       family=cumulative("probit"), 
                       chains=4, cores=4, 
                       iter=2000, warmup=1000, refresh=100, inits="0",
                       control=list(adapt_delta=0.95, max_treedepth=20),
                       prior=prior(horseshoe(3, par_ratio=0.2), class="b"),
                       save_model=glue("out{sep}initFitFull{sep}ordinal_{target}.stan"),
                       file=glue("out{sep}initFitFull{sep}ordinal_{target}"))
  
  cat("Finished", target, "\n")
}




