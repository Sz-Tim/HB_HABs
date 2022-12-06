# HABReports Bayesian modelling
# WeStCOMS data extraction
# Tim Szewczyk

# This script extracts and summarises hydrodynamic data from WeStCOMS datasets.



# set up ------------------------------------------------------------------

# devtools::install_github("https://github.com/Sz-Tim/WeStCOMS")
pkgs <- c("raster", "tidyverse", "glue", "lubridate", "sf", "jsonlite", "WeStCOMS")
lapply(pkgs, library, character.only=T)

cores <- 60
nLags <- 7
buffer_radius <- 5e3

if(.Platform$OS.type=="unix") {
  sep <- "/"
  westcoms.dir <- paste0("/media/archiver/common/sa01da-work/",
                         c("minch2/Archive/", "WeStCOMS2/Archive/"))
  mesh.f <- paste0("/home/sa04ts/FVCOM_meshes/",
                   c("WeStCOMS_mesh.gpkg", "WeStCOMS2_mesh.gpkg"))
  fetch <- "../../gis/fetch/log10_eu200m1a.tif"
} else {
  sep <- "\\"
  westcoms.dir <- paste0("D:\\hydroOut\\", 
                         c("minch2\\Archive\\", "WeStCOMS2\\Archive\\"))
  mesh.f <- paste0("C:\\Users\\sa04ts\\OneDrive - SAMS\\Projects\\WeStCOMS\\data\\",
                   c("WeStCOMS_mesh.gpkg", "WeStCOMS2_mesh.gpkg"))
  fetch <- "../../00_gis/fetch/log10_eu200m1a.tif"
}





# load files --------------------------------------------------------------

mesh.sf <- map(mesh.f, loadMesh) 
fetch.tif <- raster(fetch)

fsa.df <- fromJSON(glue("data{sep}copy_fsa.txt")) %>% 
  as_tibble %>% 
  filter(!is.na(date_collected)) %>%
  filter(karenia_mikimotoi >= 0) %>% # -99 in karenia...?
  mutate(datetime_collected=as_datetime(date_collected),
         date=date(datetime_collected),
         hour=pmin(20, pmax(5, hour(datetime_collected))),
         grid=if_else(date_collected < "2019-07-02", 1, 2),
         site.id=as.numeric(factor(paste(sin, area, site)))) %>%
  filter(datetime_collected >= "2013-07-20") %>% 
  rename(obs.id=oid) %>%
  group_by(sin, area, site) %>% 
  mutate(lon=median(easting), lat=median(northing)) %>%
  ungroup %>%
  mutate(depth_recorded=depth %>% 
           str_to_lower() %>% 
           str_remove_all("m| |\\r|\\n|j|<|>|\\+") %>%
           str_replace("su?.face", "0") %>%
           str_remove("\\(0\\)") %>%
           str_replace("bucket", "NA") %>%
           str_replace("30c", ".3") %>%
           str_replace("55", "5.5") %>%
           str_remove("[0-9]-") %>%
           str_replace("-", "NA") %>%
           as.numeric) %>%
  select(-geom, -easting, -northing, -tide, -datetime_collected)

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

fetch.df <- sampling.local %>% 
  group_by(site.id) %>%
  summarise(lon=mean(lon), lat=mean(lat)) %>%
  st_as_sf(coords=c("lon", "lat"), crs=27700) %>%
  mutate(fetch=raster::extract(fetch.tif, .))

sampling.local <- full_join(sampling.local, st_drop_geometry(fetch.df))

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





# extract hydro data ------------------------------------------------------

var.df <- tibble(var=c("temp", "salinity", "short_wave",
                       "u", "v", "ww", "km",
                       "uwind_speed", "vwind_speed", "precip"),
                 dayFn=c(mean, mean, integrateShortWave,
                         q90, q90, q90, q90,
                         q90, q90, sum),
                 depthFn=c(mean, mean, NA,
                           mean, mean, mean, mean,
                           NA, NA, NA))

# * local -----------------------------------------------------------------

sampling.local_i <- sampling.local %>% 
  select(obs.id, site.id, date, depth, grid, site.elem, starts_with("trinode"))

for(i in 0:nLags) {
  cat("Starting lag", i, "local\n")
  sampling.local_i %>%
    mutate(date=str_remove_all(as.character(date-i), "-")) %>%
    extractHydroVars(., westcoms.dir, var.df$var, var.df$dayFn, var.df$depthFn,
                     cores=cores, progress=T, returnFullDf=T) %>%
    rename_with(~paste0(.x, "_L", i), .cols=any_of(var.df$var)) %>%
    mutate(depth=0) %>%
    extractHydroVars(., westcoms.dir, "temp", list(mean), NULL, 
                     cores=cores, progress=T, returnFullDf=T) %>%
    rename(temp_0m=temp) %>%
    mutate(depth=20) %>%
    extractHydroVars(., westcoms.dir, "temp", list(mean), NULL, 
                     cores=cores, progress=T, returnFullDf=T) %>%
    rename(temp_20m=temp) %>%
    mutate(tempStrat20m=temp_20m-temp_0m) %>%
    rename_with(~paste0(.x, "_L", i), .cols="tempStrat20m") %>%
    select(obs.id, contains("_L")) %>%
    write_csv(glue("data{sep}hydro_L{i}.csv"))
}




# * regional --------------------------------------------------------------

sampling.regional_i <- sampling.regional %>% 
  select(obs.id, site.id, date, depth, grid, site.elem, starts_with("trinode"))

for(i in 0:nLags) {
  cat("Starting lag", i, "regional\n")
  sampling.regional_i %>%
    mutate(date=str_remove_all(as.character(date-i), "-")) %>%
    extractHydroVars(., westcoms.dir, var.df$var, var.df$dayFn, var.df$depthFn,
                     regional=T, cores=cores, progress=T, returnFullDf=T) %>%
    rename_with(~paste0(.x, "_R", i), .cols=any_of(var.df$var)) %>%
    mutate(depth=0) %>%
    extractHydroVars(., westcoms.dir, "temp", list(mean), NULL, 
                     regional=T, cores=cores, progress=T, returnFullDf=T) %>%
    rename(temp_0m=temp) %>%
    mutate(depth=20) %>%
    extractHydroVars(., westcoms.dir, "temp", list(mean), NULL, 
                     regional=T, cores=cores, progress=T, returnFullDf=T) %>%
    rename(temp_20m=temp) %>%
    mutate(tempStrat20m=temp_20m-temp_0m) %>%
    rename_with(~paste0(.x, "_R", i), .cols="tempStrat20m") %>%
    group_by(obs.id) %>%
    summarise(across(contains("_R"), mean)) %>%
    ungroup %>%
    write_csv(glue("data{sep}hydro_R{i}.csv"))
}

