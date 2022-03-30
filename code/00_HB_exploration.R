# HABReports Bayesian modelling
# Bayesian HAB model exploration
# Tim Szewczyk


# This script is a translation of parts of docs/02_HB_exploration.Rmd for running remotely



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "ncdf4", "sf", "LaplacesDemon",
          "cmdstanr", "bayesplot", "posterior", "jsonlite")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
walk(dir("code", "00[0-9]_fn", full.names=T), source)

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


fsa.df <- fromJSON(glue("data{sep}copy_fsa.txt")) %>% 
  as_tibble %>% 
  select(-geom) %>%
  filter(easting < 4e10) %>% # entry error: Camb - Mid Yell Voe - 2013-04-16 
  filter(easting > 0 & northing > 0) %>%
  filter(!is.na(date_collected)) %>%
  filter(karenia_mikimotoi >= 0) %>% # -99 in karenia...?
  mutate(datetime_collected=as_datetime(date_collected),
         date_collected=date(datetime_collected),
         year=year(date_collected),
         month=month(date_collected),
         hour=pmin(20, pmax(5, hour(datetime_collected))),
         minutes=minute(datetime_collected),
         grid=if_else(date_collected < "2019-05-01", "minch2", "WeStCOMS2")) %>%
  filter(datetime_collected >= "2013-07-20") %>% 
  group_by(sin, area, site) %>% 
  mutate(lon=mean(easting), lat=mean(northing)) %>%
  ungroup %>% 
  mutate(site.id=as.numeric(factor(paste(sin, area, site))))







# munge datasets ----------------------------------------------------------

par.ls <- setSimParams()
# generate sampling details and lnlambdas
sampling.df <- fsa.df %>% 
  mutate(date=date_collected,
         dateChar=str_remove_all(date, "-"),
         time=hour + minutes/60,
         mode=sample(par.ls$modes, n(), T, c(0.7, 0.2, 0.1)),
         depth=sample(par.ls$depths, n(), T),
         tides=sample(par.ls$tides, n(), T, c(0.2, 0.05, 0.5, 0.05, 0.3)),
         lnlambda=rnorm(n(), 7, 2),
         lnlambdap1=log(exp(lnlambda)+1))

# extract element IDs, trinodes for each observation
mesh.sf <- map(mesh.f, st_read)
mesh <- map(mesh.f, loadMesh)
site.xy <- fsa.df %>% select("grid", "site.id", "lon", "lat") %>%
  group_by(grid, site.id) %>% slice_head(n=1) %>%
  group_by(grid) %>%
  group_split() %>%
  map2(.x=., .y=mesh.sf, 
       ~st_as_sf(.x, coords=c("lon", "lat"), remove=F, crs=27700) %>% 
         st_join(., .y %>% 
                   select(-area) %>% 
                   rename(depth.elem=depth,
                          site.elem=i)) %>% 
         filter(!is.na(site.elem)) %>%
         st_drop_geometry) %>%
  bind_rows
sampling.df <- sampling.df %>%
  right_join(., site.xy, by=c("grid", "site.id", "lon", "lat")) %>%
  mutate(short_wave=NA_real_, zeta=NA_real_, temp=NA_real_) %>%
  group_by(grid) %>% group_split
site_trinode <- map2(.x=mesh, .y=sampling.df,
                     ~ncvar_get(.x, "trinodes")[.y$site.elem,])
walk(mesh, nc_close)

# load relevant data
for(grid in 1:2) {
  sampling_dates <- unique(sampling.df[[grid]]$dateChar)
  hydrovars.ls <- vector("list", length(sampling_dates))
  for(i in 1:length(sampling_dates)) {
    date_i <- sampling_dates[i]
    rows_i <-  which(sampling.df[[grid]]$dateChar == date_i)
    hydroVars.ls[[i]] <- loadHydroVars(date_i, 
                                       site_trinode[[grid]][rows_i,], 
                                       sampling.df[[grid]]$hour[rows_i],
                                       sampling.df[[grid]]$depth[rows_i],
                                       westcoms.dir[[grid]], 
                                       sep, 
                                       vars=c("temp", "short_wave", "zeta"), 
                                       lag=2) %>%
      mutate(rows=rows_i)
  }
  sampling.df[[grid]] <- sampling.df[[grid]] %>%
    bind_cols(do.call(rbind, hydrovars.ls) %>% arrange(rows))
}
sampling.df <- sampling.df %>%
  bind_rows %>%
  mutate(across(one_of("time", "depth", "tides", "short_wave", "zeta", "temp"), 
                CenterScale, 
                .names="{.col}_sc"))

write_csv(sampling.df, glue("data{sep}obs_df.csv"))




# prepare model -----------------------------------------------------------
# 
# sampling.df <- read_csv(glue("data{sep}obs_df.csv")) %>% 
#   mutate(target_sp=round(alexandrium_sp/20)) %>%
#   # filter(target_sp > 20) %>%
#   mutate(#short_wave=log(short_wave+1),
#          yday=yday(date_collected),
#          ydayCos=cos(2*pi*yday/365),
#          ydaySin=sin(2*pi*yday/365)) %>%
#   mutate(across(one_of("time", "depth", "tides", "short_wave", "zeta", "temp"), 
#                 CenterScale, 
#                 .names="{.col}_sc"))
#   
# samp.mx <- model.matrix(~ 0 + time_sc + I(time_sc^2) + 
#                           short_wave_sc + zeta_sc + temp_sc, 
#                         data=sampling.df)
# 
# mod <- cmdstan_model(glue("models{sep}00_sampling_simple.stan"))
# sampling_hydro.data <- list(
#   N=nrow(sampling.df),
#   X=samp.mx,
#   nCov=ncol(samp.mx),
#   y=sampling.df$target_sp,
#   prior_ln_lambda=c(mean(log(sampling.df$target_sp)), 
#                     sd(log(sampling.df$target_sp)))
# )






# run model ---------------------------------------------------------------

# fit <- mod$sample(data=sampling_hydro.data,
#                   chains=4, parallel_chains=4,
#                   max_treedepth=20,
#                   refresh=500, iter_sampling=2000)
# fit$cmdstan_diagnose()
# out.sum <- fit$summary()
# 
# mcmc_areas(fit$draws("beta")) + scale_y_discrete(labels=colnames(samp.mx))
# ppc_dens_overlay(y=log(sampling_hydro.data$y+1),
#                  yrep=log(fit$draws("y_pred", format="matrix")[1:100,]+1))
# 
# out.df <- sampling.df %>% 
#   mutate(lnlambda_mean=filter(out.sum, grepl("ln_lambda", variable))$mean,
#          lnlambda_q5=filter(out.sum, grepl("ln_lambda", variable))$q5,
#          lnlambda_q95=filter(out.sum, grepl("ln_lambda", variable))$q95)
# ggplot(out.df, aes(log(target_sp), lnlambda_mean)) + geom_point() +
#   geom_linerange(aes(ymin=lnlambda_q5, ymax=lnlambda_q95)) +
#   geom_abline()
# ggplot(out.df, aes(target_sp, exp(lnlambda_mean))) + geom_point() +
#   geom_linerange(aes(ymin=exp(lnlambda_q5), ymax=exp(lnlambda_q95))) +
#   geom_abline()
# 
# ggplot(out.df, aes(temp_sc, log(target_sp)-lnlambda_mean)) + 
#   geom_point()
# ggplot(out.df, aes(time_sc, log(target_sp)-lnlambda_mean)) + 
#   geom_point()
# ggplot(out.df, aes(short_wave_sc, log(target_sp)-lnlambda_mean)) + 
#   geom_point()
# ggplot(out.df, aes(zeta_sc, log(target_sp)-lnlambda_mean)) + 
#   geom_point()
# 
# 
# 
# # Is this doing *anything* besides shrinking values toward the mean?
# 
# 
# library(brms)
# out <- brm(bf(target_sp ~ short_wave_sc + zeta_sc + temp_sc + ydayCos + ydaySin +
#              temp_sc:short_wave_sc +  short_wave_sc:ydayCos + temp_sc:ydayCos +
#              (1|site.id),
#              zi ~ temp_sc + ydayCos + ydaySin +
#                (1|site.id)), 
#            data=sampling.df, family=zero_inflated_negbinomial(), cores=4)
# pp_check(out)
# plot(conditional_effects(out), points=T, point_args=list(shape=1, alpha=0.5, size=0.75))
# # I think these values are too large for a poisson distribution and that I
# # should consider alternatives. For example, a hurdle model with a normal
# # distribution using log-transformed values might be an ok solution.
# 
