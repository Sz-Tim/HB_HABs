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



makeSimSamples <- function(par.ls, westcoms.dir=NULL, mesh.f=NULL, site.xy=NULL, 
                           sep=sep, startDate="20190401", lnlambda_prior=c(7,2)) {
  library(tidyverse); library(lubridate); library(glue); library(LaplacesDemon)
  
  # generate sampling details and lnlambdas
  sampling.df <- tibble(sampleID=1:(par.ls$N),
                        site=rep(1:par.ls$nSites, times=par.ls$nDays),
                        day=rep(0:(par.ls$nDays-1), each=par.ls$nSites)) %>%
    mutate(date=ymd(startDate)+day,
           dateChar=str_remove_all(date, "-"),
           time=rbeta(n(), 2, 5)*12+7,
           hour=floor(time),
           minutes=round(time %% 1 * 60),
           mode=sample(par.ls$modes, n(), T, c(0.7, 0.2, 0.1)),
           depth=sample(par.ls$depths, n(), T),
           tides=sample(par.ls$tides, n(), T, c(0.2, 0.05, 0.5, 0.05, 0.3)),
           lnlambda=rnorm(n(), lnlambda_prior[1], lnlambda_prior[2]),
           lnlambdap1=log(exp(lnlambda)+1))

  # extract element IDs, trinodes for each observation
  mesh.sf <- st_read(mesh.f)
  mesh <- loadMesh(mesh.f)
  site.xy <- site.xy %>% 
    st_as_sf(coords=c("lon", "lat"), remove=F, crs=27700) %>% 
    st_join(., mesh.sf) %>% 
    filter(!is.na(i))
  sampling.df <- sampling.df %>%
    mutate(site.id=sample(1:nrow(site.xy), n(), replace=T),
           site.lat=site.xy$lat[site.id],
           site.lon=site.xy$lon[site.id],
           site.elem=site.xy$i[site.id],
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
  
  return(sampling.df)
}




