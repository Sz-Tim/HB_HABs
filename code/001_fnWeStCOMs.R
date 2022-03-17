# HABReports Bayesian modelling
# WeStCOMS processing functions
# Tim Szewczyk



loadMesh <- function(mesh.f) {
  mesh <- nc_open(str_replace(mesh.f, "gpkg", "nc"))
}




loadMesh_sf <- function(path) {
  library(ncdf4); library(sf); library(tidyverse)
  mesh <- loadMesh(path)
  trinodes <- ncvar_get(mesh, "trinodes")
  nodexy_os <- ncvar_get(mesh, "nodexy_os")
  elemCoord <- map(1:nrow(trinodes), 
                   ~rbind(nodexy_os[trinodes[.x,1],], 
                          nodexy_os[trinodes[.x,2],],
                          nodexy_os[trinodes[.x,3],],
                          nodexy_os[trinodes[.x,1],]))
  mesh.sf <- map_df(seq_along(elemCoord), 
                    ~st_polygon(list(elemCoord[[.x]])) %>% 
                      st_sfc %>%
                      st_sf(i=.x, 
                            geometry=.)) %>%
    st_set_crs(27700) 
  nc_close(mesh)
  
  return(mesh.sf)
}



addCoords <- function(sampling.df, mesh, elems, trinodes) {
  sampling.df %>%
    mutate(lonc=ncvar_get(mesh, "uvnode_os")[elems,1][site],
           latc=ncvar_get(mesh, "uvnode_os")[elems,2][site],
           lon1=ncvar_get(mesh, "nodexy_os")[trinodes[,1],1][site],
           lon2=ncvar_get(mesh, "nodexy_os")[trinodes[,2],1][site],
           lon3=ncvar_get(mesh, "nodexy_os")[trinodes[,3],1][site],
           lat1=ncvar_get(mesh, "nodexy_os")[trinodes[,1],2][site],
           lat2=ncvar_get(mesh, "nodexy_os")[trinodes[,2],2][site],
           lat3=ncvar_get(mesh, "nodexy_os")[trinodes[,3],2][site])
}



loadHydroVars <- function(sampling.df, i_date, site_trinode, westcoms.dir, sep) {
  library(ncdf4); library(tidyverse); library(glue)
  i_rows <- which(sampling.df$dateChar == i_date)
  i_dir <- glue("{westcoms.dir}{sep}netcdf_{str_sub(i_date,1,4)}")
  hydro.f <- dir(i_dir, i_date)
  i_nc <- nc_open(glue("{i_dir}{sep}{hydro.f}"))
  short_wave <- meanOfNodes(ncvar_get(i_nc, "short_wave"), site_trinode)
  temp <- meanOfNodes(ncvar_get(i_nc, "temp"), site_trinode)
  zeta <- meanOfNodes(ncvar_get(i_nc, "zeta"), site_trinode)
  waterDepth <- zeta + c(meanOfNodes(ncvar_get(i_nc, "h"), site_trinode))
  siglay <- abs(ncvar_get(i_nc, "siglay")[1,])
  nc_close(i_nc)
  
  for(j in i_rows) {
    sampling.df$short_wave[j] <- integrateShortWave(short_wave, 
                                                    j, 
                                                    sampling.df$time[j])
    sampling.df$zeta[j] <- zeta[j,sampling.df$hour[j]]
    j_siglay <- which.min(abs(waterDepth[j] * siglay - sampling.df$depth[j]))
    sampling.df$temp[j] <- temp[j,j_siglay, sampling.df$hour[j]]
  }
  
  return(sampling.df)
}





meanOfNodes <- function(nc.ar, node.mx) {
  if (length(dim(nc.ar))==1) {
    node1 <- nc.ar[node.mx[,1]]
    node2 <- nc.ar[node.mx[,2]]
    node3 <- nc.ar[node.mx[,3]]
  } else if(length(dim(nc.ar))==2) {
    node1 <- nc.ar[node.mx[,1],]
    node2 <- nc.ar[node.mx[,2],]
    node3 <- nc.ar[node.mx[,3],]
  } else if(length(dim(nc.ar))==3) {
    node1 <- nc.ar[node.mx[,1],,]
    node2 <- nc.ar[node.mx[,2],,]
    node3 <- nc.ar[node.mx[,3],,]
  }
  return((node1 + node2 + node3)/3)
}




integrateShortWave <- function(nc.var, rowNum, endTime, startTime=0) {
  # linear interpolation between hours
  # short_wave units = joules/s
  hourlyJoules <- rep(0, ceiling(endTime)-floor(startTime))
  for(hour in seq_along(hourlyJoules)) {
    varStart <- nc.var[rowNum,hour]
    varEnd <- nc.var[rowNum,hour+1]
    dt <- 3600
    if(hour == ceiling(endTime)) {
      propHour <- endTime %% 1
      varEnd <- dVar*propHour + varStart
      dt <- 3600*propHour
    } 
    dVar <- varEnd - varStart
    hourlyJoules[hour] <- dt*min(varStart, varEnd) + (0.5*dt*abs(dVar))
  }
  return(sum(hourlyJoules))
}