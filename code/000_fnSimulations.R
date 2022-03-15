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



makeSimSamples <- function(par.ls, westcoms.dir=NULL, mesh_sf.f=NULL, site.xy=NULL, 
                           startDate="20190401", lnlambda_prior=c(7,2)) {
  
  sampling.df <- tibble(sampleID=1:(par.ls$N),
                        site=rep(1:par.ls$nSites, times=par.ls$nDays),
                        day=rep(0:(par.ls$nDays-1), each=par.ls$nSites)) %>%
    mutate(date=ymd(startDate)+day,
           time=rbeta(n(), 2, 5)*12+7,
           hour=floor(time),
           minutes=round(time %% 1 * 60),
           mode=sample(par.ls$modes, n(), T, c(0.7, 0.2, 0.1)),
           depth=sample(par.ls$depths, n(), T),
           tides=sample(par.ls$tides, n(), T, c(0.2, 0.05, 0.5, 0.05, 0.3)),
           lnlambda=rnorm(n(), lnlambda_prior[1], lnlambda_prior[2]),
           lnlambdap1=log(exp(lnlambda)+1))
  
  if(!is.null(westcoms.dir)) {
    # TODO: extract most or all into separate functions
    # TODO: actual site locations
    
    # random site locations for now
    if(is.null(site.xy)){
      mesh <- loadMesh(westcoms.dir)
      site_elem <- sample(mesh$dim$elems$vals, par.ls$nSites)
      site_trinode <- ncvar_get(mesh, "trinodes")[site_elem,]
      sampling.df <- sampling.df %>%
        addCoords(mesh, site_elem, site_trinode) %>%
        mutate(short_wave=NA, zeta=NA, temp=NA)
      nc_close(mesh)
    } else {
      mesh.sf <- st_read(mesh_sf.f)
      mesh <- loadMesh(westcoms.dir)
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
    }
    
    
    
    # load relevant data
    for(i in 1:n_distinct(sampling.df$date)) {
      i_date <- unique(sampling.df$date)[i]
      i_rows <- which(sampling.df$date == i_date)
      
      hydro.f <- dir(westcoms.dir, str_remove_all(i_date, "-"))
      i_nc <- nc_open(paste0(westcoms.dir, hydro.f))
      short_wave <- meanOfNodes(ncvar_get(i_nc, "short_wave"), site_trinode)
      temp <- meanOfNodes(ncvar_get(i_nc, "temp"), site_trinode)
      zeta <- meanOfNodes(ncvar_get(i_nc, "zeta"), site_trinode)
      waterDepth <- zeta + c(meanOfNodes(ncvar_get(i_nc, "h"), site_trinode))
      siglay <- abs(ncvar_get(i_nc, "siglay")[1,])
      nc_close(i_nc)
      
      for(j in i_rows) {
        sampling.df$short_wave[j] <- integrateShortWave(short_wave, 
                                                        sampling.df$site[j], 
                                                        sampling.df$time[j])
        sampling.df$zeta[j] <- zeta[sampling.df$site[j],sampling.df$hour[j]]
        j_siglay <- which.min(abs(waterDepth[sampling.df$site[j]] * siglay - 
                                    sampling.df$depth[j]))
        sampling.df$temp[j] <- temp[sampling.df$site[j],j_siglay, sampling.df$hour[j]]
      }
    }
    
    sampling.df <- sampling.df %>%
      mutate(across(one_of("time", "depth", "tides", "short_wave", "zeta", "temp"), 
                    CenterScale, 
                    .names="{.col}_sc"))
  }
  
  return(sampling.df)
}




