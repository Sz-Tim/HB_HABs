# HABReports Bayesian modelling
# Simulation functions
# Tim Szewczyk



setSimParams <- function(nSites=45, nDays=30,
                         modes=c("10m Tube Sampler", "Bucket", "Other"),
                         depths=c(0, 1, 2, 4, 5, 7, 10),
                         tides=c(-1, -0.5, 0, 0.5, 1)) {
  par.ls <- list(nSites=nSites, nDays=nDays, N=nSites*nDays,
                 modes=modes, depths=depths, tides=tides)
  return(par.ls)
}



makeSimSamples <- function(par.ls, westcoms.dir=NULL, mesh.f=NULL,
                           obs.df, sep=sep, lnlambda_prior=c(7,2)) {
  library(tidyverse); library(lubridate); library(glue); library(LaplacesDemon)
  
  # generate sampling details and lnlambdas
  sampling.df <- obs.df %>% 
    mutate(date=date_collected,
           dateChar=str_remove_all(date, "-"),
           time=hour + minutes/60,
           mode=sample(par.ls$modes, n(), T, c(0.7, 0.2, 0.1)),
           depth=sample(par.ls$depths, n(), T),
           tides=sample(par.ls$tides, n(), T, c(0.2, 0.05, 0.5, 0.05, 0.3)),
           lnlambda=rnorm(n(), lnlambda_prior[1], lnlambda_prior[2]),
           lnlambdap1=log(exp(lnlambda)+1))

  # extract element IDs, trinodes for each observation
  mesh.sf <- map(mesh.f, st_read)
  mesh <- map(mesh.f, loadMesh)
  site.xy <- obs.df %>% select("grid", "site.id", "lon", "lat") %>%
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
    group_by(grid) %>% group_split
  site_trinode <- map2(.x=mesh, .y=sampling.df,
                       ~ncvar_get(.x, "trinodes")[.y$site.elem,])
  walk(mesh, nc_close)
    
  # load relevant data
  for(grid in 1:2) {
    sampling_dates <- unique(sampling.df[[grid]]$dateChar)
    for(i in 1:length(sampling_dates)) {
      sampling.df[[grid]] <- sampling.df[[grid]] %>%
        loadHydroVars(sampling_dates[i], site_trinode[[grid]], westcoms.dir[[grid]], sep)
    }
  }
  sampling.df <- sampling.df %>% bind_rows
  
  return(sampling.df)
}




