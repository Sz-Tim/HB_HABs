## -----------------------------------------------------------------------------
## HAB forecasts
## Initial setup
## Tim Szewczyk
## -----------------------------------------------------------------------------


# set up ------------------------------------------------------------------
library(tidyverse)
library(glue)
library(jsonlite)
library(lubridate)
library(WeStCOMS)
source("code/00_fn.R")

# directories
proj.dir <- getwd()
part.dir <- "/home/sa04ts/HAB_operational/init/"
wc_base.dir <- "/media/archiver/common/sa01da-work/"
mesh_base.dir <- "/home/sa04ts/FVCOM_meshes/"

cores <- 24

# WeStCOMS dates
v1_start <- "2014-10-28"
v1_end <- "2019-04-01"
v2_end <- "2022-06-01"

# site locations
fsa_v1 <- read_delim(paste0(part.dir, "fsa_sites_v1.dat"), 
                     delim="\t", col_names=c("site.id", "x", "y")) %>%
  filter(!site.id %in% c(70, 74, 75, 80, 88))
fsa_v2 <- read_delim(paste0(part.dir, "fsa_sites_v2.dat"), 
                     delim="\t", col_names=c("site.id", "x", "y"))



# particle tracking -------------------------------------------------------
setwd(part.dir)

# WeStCOMS v1
wc1.prop <- setPartTrackProperties(parallelThreads=cores, 
                                   destinationDirectory="out_v1/",
                                   datadir=paste0(wc_base.dir, "minch2/Archive/"),
                                   mesh1=paste0(mesh_base.dir, "WeStCOMS_mesh.nc"),
                                   location="minch",
                                   minchVersion="2",
                                   sitefile="large_elem_centroids_v1.dat",
                                   sitefileEnd="fsa_sites_v1.dat",
                                   start_ymd=str_remove_all(v1_start, "-"),
                                   end_ymd=str_remove_all(v1_end, "-"))
cat(wc1.prop, "\n", file="HAB_WeStCOMS1.properties")
system2("bash", c("runPTrack_smn.sh", "HAB_WeStCOMS1.properties"))

# WeStCOMS v2
wc2.prop <- setPartTrackProperties(parallelThreads=cores,
                                   destinationDirectory="out_v2/",
                                   restartParticles="out_v1/locationsEnd.dat",
                                   start_ymd=str_remove_all(v1_end, "-"),
                                   end_ymd=str_remove_all(v2_end, "-"))
cat(wc2.prop, "\n", file="HAB_WeStCOMS2.properties")
system2("bash", c("runPTrack_smn.sh", "HAB_WeStCOMS2.properties"))

# compile output
connect.f <- dir("out_v1/connectivity")
connect.v1 <- connect.f %>%
  map_dfr(~read_delim(glue("out_v1/connectivity/{.x}"), 
                      delim=" ", col_names=F, show_col_types=F) %>% 
            mutate(date=str_sub(.x, 14, 21)) %>%
            group_by(date) %>% 
            summarise(across(everything(), sum)) %>% 
            ungroup) %>%
  pivot_longer(-1, names_to="site.id", values_to="influx") %>%
  mutate(site.id=rep(fsa_v1$site.id, length(connect.f)),
         date=ymd(date))
connect.f <- dir("out_v2/connectivity")
connect.v2 <- connect.f %>%
  map_dfr(~read_delim(glue("out_v2/connectivity/{.x}"), 
                      delim=" ", col_names=F, show_col_types=F) %>% 
            mutate(date=str_sub(.x, 14, 21)) %>%
            group_by(date) %>% 
            summarise(across(everything(), sum)) %>% 
            ungroup) %>%
  pivot_longer(-1, names_to="site.id", values_to="influx") %>%
  mutate(site.id=rep(fsa_v2$site.id, length(connect.f)),
         date=ymd(date))

bind_rows(connect.v1 %>% filter(! date %in% connect.v2$date), 
          connect.v2) %>%
  write_csv(paste0(proj.dir, "/data/influx_init.csv"))

setwd(proj.dir)




# WeStCOMS ----------------------------------------------------------------

nLags <- 7
buffer_radius <- 5e3
westcoms.dir <- paste0(wc_base.dir,
                       c("minch2/Archive/", "WeStCOMS2/Archive/"))
mesh.f <- paste0(mesh_base.dir,
                 c("WeStCOMS_mesh.gpkg", "WeStCOMS2_mesh.gpkg"))

# load files
mesh.sf <- map(mesh.f, loadMesh) 
fsa.df <- glue("data/copy_fsa.txt") %>% fromJSON() %>% clean_fsa(v1_end)
sampling.local <- fsa.df %>%
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
  bind_rows %>%
  mutate(depth=pmin(10, depth.elem),
         obs.id=row_number())
sampling.regional <- sampling.local %>%
  select(-site.elem, -depth.elem, -starts_with("trinode")) %>%
  group_by(grid) %>%
  group_split() %>%
  map2(.x=., .y=mesh.sf, 
       ~st_as_sf(.x, coords=c("lon", "lat"), remove=F, crs=27700) %>% 
         st_buffer(dist=buffer_radius) %>%
         st_join(., .y %>% 
                   select(-area) %>% 
                   rename(depth.elem=depth,
                          site.elem=i)) %>% 
         filter(!is.na(site.elem)) %>%
         st_drop_geometry) %>%
  bind_rows %>%
  mutate(depth=pmin(10, depth.elem))
write_csv(sampling.local, glue("data{sep}sampling_local.csv"))
write_csv(sampling.regional, glue("data{sep}sampling_regional.csv"))

# extract hydro data
var.df <- tibble(var=c("temp", "salinity", "short_wave",
                       "u", "v", "ww", "km",
                       "uwind_speed", "vwind_speed", "precip"),
                 dayFn=c(mean, mean, integrateShortWave,
                         q90, q90, q90, q90,
                         q90, q90, sum),
                 depthFn=c(mean, mean, NA,
                           mean, mean, mean, mean,
                           NA, NA, NA))
hydro_L <- hydro_R <- vector("list", nLags)
for(i in 0:nLags) {
  cat("Starting lag", i, "\n")
  hydro_L[[i+1]] <- sampling.local %>% 
    select(obs.id, site.id, date, depth, grid, site.elem, starts_with("trinode")) %>%
    mutate(date=str_remove_all(as.character(date-i), "-")) %>%
    extractHydroVars(., westcoms.dir, var.df$var, var.df$dayFn, var.df$depthFn,
                     cores=cores, progress=T) %>%
    rename_with(~paste0(.x, "_L", i), .cols=any_of(var.df$var))
  hydro_R[[i+1]] <- sampling.regional %>% 
    select(obs.id, site.id, date, depth, grid, site.elem, starts_with("trinode")) %>%
    mutate(date=str_remove_all(as.character(date-i), "-")) %>%
    extractHydroVars(., westcoms.dir, var.df$var, var.df$dayFn, var.df$depthFn,
                     regional=T, cores=cores, progress=T) %>%
    rename_with(~paste0(.x, "_R", i), .cols=any_of(var.df$var))
}
hydro.df <- c(hydro_L, hydro_R) %>%
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
write_csv(glue("data{sep}hydro_init.csv"))



# Copernicus --------------------------------------------------------------






# algal densities ---------------------------------------------------------







# toxin concentrations ----------------------------------------------------








# compile -----------------------------------------------------------------


