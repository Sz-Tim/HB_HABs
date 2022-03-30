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





loadHydroVars <- function(date, trinodes, hours, depths, westcoms.dir, sep, 
                          vars=c("temp", "short_wave", "zeta"), 
                          lags=NULL, dayAvg=FALSE) {
  library(ncdf4); library(tidyverse); library(glue)
  elems <- 1:nrow(trinodes)
  
  # date = lag_0; date - n = lag_n 
  if(is.null(lags)) {
    vars_all <- paste0(vars, "_lag_0")
    dates <- date
  } else {
    vars_all <- unlist(map(0:lags, ~paste0(vars, "_lag_", .x)))
    dates <- map_chr(0:lags, ~str_remove_all(ymd(date)-.x, "-"))
  }
  
  out <- vector("list", length(vars_all)) %>% setNames(vars_all)
  hydro_temp <- vector("list", length(vars)) %>% setNames(vars)
  
  for(i in seq_along(dates)) {
    
    # load data for focal elements for each date from 0:lags
    dir_i <- glue("{westcoms.dir}{sep}netcdf_{str_sub(dates[i],1,4)}")
    file_i <- dir(dir_i, dates[i])
    nc_i <- nc_open(glue("{dir_i}{sep}{file_i}"))
    for(v in seq_along(vars)) {
      hydro_temp[[vars[v]]] <- meanOfNodes(ncvar_get(nc_i, vars[v]), trinodes)
    }
    waterDepth <- c(meanOfNodes(ncvar_get(nc_i, "h"), trinodes))
    siglay <- abs(ncvar_get(nc_i, "siglay")[1,])
    if("zeta" %in% vars) {
      waterDepth <- hydro_temp[["zeta"]] + waterDepth
    }
    nc_close(nc_i)
    
    # calculate and extract appropriate values for each element
    for(v in seq_along(vars)) {
      name_vi <- glue("{vars[v]}_lag_{i-1}")
      # indexes = 1:24, hours = 0:23
      
      # short_wave: integrate
      if(vars[v]=="short_wave") {
        if(dayAvg | i > 1) {
          out[[name_vi]] <- map_dbl(elems,
                                    ~integrateShortWave(hydro_temp[[vars[v]]], .x, 23)) 
        } else {
          out[[name_vi]] <- map_dbl(elems,
                                    ~integrateShortWave(hydro_temp[[vars[v]]], .x, hours[.x])) 
        }
      }
      
      # [elem, layer, hour]
      if(vars[v] %in% c("temp", "u", "v")) {
        if(dayAvg | i > 1) {
          # average temperature at a given depth
          siglays <- apply(depths/waterDepth, 1:2, function(x) which.min(abs(siglay - x)))
          indexes <- cbind(rep(elems, 24), c(siglays), rep(1:24, each=length(elems)))
          out[[name_vi]] <- rowMeans(matrix(hydro_temp[[vars[v]]][indexes], nrow=length(elems)))
        } else {
          indexes <- cbind(elems, hours+1)
          siglays <- apply(cbind(waterDepth[indexes]) %*% rbind(siglay) - depths,
                           1, function(x) which.min(abs(x)))
          indexes <- cbind(elems, siglays, hours+1)
          out[[name_vi]] <- hydro_temp[[vars[v]]][indexes] 
        }
      }
      
      # [elem, hour]
      if(vars[v] %in% c("zeta", "uwind_speed", "vwind_speed")) {
        if(dayAvg | i > 1) {
          out[[name_vi]] <- rowMeans(hydro_temp[[vars[v]]])
        } else {
          indexes <- cbind(elems, hours+1)
          out[[name_vi]] <- hydro_temp[[vars[v]]][indexes]
        }
      }
    }
  }
  
  return(as_tibble(out))
}







meanOfNodes <- function(nc.ar, node.mx) {
  if (length(dim(nc.ar))==1) {
    node1 <- nc.ar[node.mx[,1]]
    node2 <- nc.ar[node.mx[,2]]
    node3 <- nc.ar[node.mx[,3]]
  } else if(length(dim(nc.ar))==2) {
    node1 <- nc.ar[node.mx[,1],,drop=F]
    node2 <- nc.ar[node.mx[,2],,drop=F]
    node3 <- nc.ar[node.mx[,3],,drop=F]
  } else if(length(dim(nc.ar))==3) {
    node1 <- nc.ar[node.mx[,1],,,drop=F]
    node2 <- nc.ar[node.mx[,2],,,drop=F]
    node3 <- nc.ar[node.mx[,3],,,drop=F]
  }
  return((node1 + node2 + node3)/3)
}




integrateShortWave <- function(nc.var, rowNum, endTime, startTime=0) {
  # linear interpolation between hours
  # short_wave units = joules/s
  hourlyJoules <- rep(0, ceiling(endTime)-floor(startTime))
  for(hour in seq_along(hourlyJoules)) {
    varStart <- nc.var[rowNum,hour] # indexes = 1:24, hours = 0:23
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
