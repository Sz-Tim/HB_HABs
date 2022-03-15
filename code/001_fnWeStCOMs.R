# HABReports Bayesian modelling
# WeStCOMS processing functions
# Tim Szewczyk



loadMesh <- function(path) {
  if(grepl("WestCOMS2", path)) {
    mesh <- nc_open("D:\\hydroOut\\WestCOMS2_Mesh.nc") 
  } else if (grepl("minch", path) || grepl("WestCOMS_")) {
    mesh <- nc_open("D:\\hydroOut\\WestCOMS_meshOS.nc")
  } else {
    stop(simpleError("Hydrodynamic data must be WestCOMS2 or minch"))
  }
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




findElement <- function(site.xy, mesh.sf) {
  
  return(site.elem)
}



findTrinodes <- function(site.elem=NULL, mesh.sf, site.xy=NULL) {
  if(is.null(site.elem)) {
    if(is.null(site.xy)) {
      stop("Error: must provide site coordinates or mesh element ids")
    }
    
    
  }
  return(site.trinodes)
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




integrateShortWave <- function(nc.var, siteNum, endTime, startTime=0) {
  # linear interpolation between hours
  # short_wave units = joules/s
  hourlyJoules <- rep(0, ceiling(endTime)-floor(startTime))
  for(hour in seq_along(hourlyJoules)) {
    varStart <- nc.var[siteNum,hour]
    varEnd <- nc.var[siteNum,hour+1]
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