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
                 "/home/sa04ts/FVCOM_meshes/WeStCOMS2_Mesh.gpkg",
                 "..\\..\\01_FVCOM\\data\\WeStCOMS2_Mesh.gpkg")


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
fsa.sf <- st_as_sf(fsa.df, coords=c("easting", "northing"), crs=st_crs(27700))






# simulation settings -----------------------------------------------------

hydro.par <- setSimParams(nDays=29)
site.xy <- fsa.df %>% group_by(sin, area, site) %>% 
  summarise(lon=mean(easting), lat=mean(northing))
sampling.df <- makeSimSamples(hydro.par, westcoms.dir, mesh.f, site.xy, sep)


