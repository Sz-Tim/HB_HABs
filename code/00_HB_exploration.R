# HABReports Bayesian modelling
# Bayesian HAB model exploration
# Tim Szewczyk


# This script is a translation of parts of docs/02_HB_exploration.Rmd for running remotely



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "ncdf4", "sf",
          "cmdstanr", "bayesplot", "posterior", "jsonlite")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
walk(dir("code", "00[0-9]_fn", full.names=T), source)
sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")

westcoms.dir <- ifelse(.Platform$OS.type=="unix",
                       "/media/archiver/common/sa01da-work/WeStCOMS2/Archive/",
                       "D:\\hydroOut\\WestCOMS2\\Archive\\")
mesh.f <- ifelse(.Platform$OS.type=="unix",
                 "/home/sa04ts/FVCOM_meshes/WeStCOMS2_mesh.gpkg",
                 "..\\..\\01_FVCOM\\data\\WeStCOMS2_mesh.gpkg")


fsa.df <- fromJSON(glue("data{sep}copy_fsa.txt")) %>% 
  as_tibble %>% 
  select(-geom) %>%
  filter(easting < 4e10) %>% # entry error: Camb - Mid Yell Voe - 2013-04-16 
  filter(easting > 0 & northing > 0) %>%
  filter(!is.na(date_collected)) %>%
  filter(karenia_mikimotoi >= 0) %>% # -99 in karenia...?
  mutate(date_collected=date(date_collected),
         year=year(date_collected),
         month=month(date_collected),
         yday=yday(date_collected),
         ydayCos=cos(yday/365),
         ydaySin=sin(yday/365),
         x_y=paste0(easting, "_", northing)) %>%
  filter(year > 2019 & year < 2022)






# simulation settings -----------------------------------------------------

hydro.par <- setSimParams(nDays=29)
mesh.sf <- st_read(mesh.f)
mesh <- loadMesh(mesh.f)
sampling.df <- fsa.df %>%
  mutate(dateChar=str_remove_all(date_collected, "-"),
         time=rbeta(n(), 2, 5)*12+7,
         hour=floor(time),
         minutes=round(time %% 1 * 60),
         mode=sample(hydro.par$modes, n(), T, c(0.7, 0.2, 0.1)),
         depth=sample(hydro.par$depths, n(), T),
         tides=sample(hydro.par$tides, n(), T, c(0.2, 0.05, 0.5, 0.05, 0.3))) %>%
  st_as_sf(coords=c("easting", "northing"), remove=F, crs=27700) %>% 
  st_join(., mesh.sf %>% select(i)) %>% 
  filter(!is.na(i)) %>%
  rename(site.elem=i) %>%
  mutate(site.id=as.numeric(as.factor(site.elem)),
         short_wave=NA, zeta=NA, temp=NA)
site_trinode <- ncvar_get(mesh, "trinodes")[sampling.df$site.elem,]
nc_close(mesh)

# load relevant data
sampling_dates <- unique(sampling.df$dateChar)
for(i in 1:length(sampling_dates)) {
  sampling.df <- sampling.df %>%
    loadHydroVars(sampling_dates[i], site_trinode, westcoms.dir, sep)
}
sampling.df <- sampling.df %>%
  mutate(across(one_of("time", "depth", "tides", "short_wave", "zeta", "temp"), 
                CenterScale, 
                .names="{.col}_sc"))
