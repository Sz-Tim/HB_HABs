## -----------------------------------------------------------------------------
## HAB forecasts
## Functions
## Tim Szewczyk
## -----------------------------------------------------------------------------


clean_fsa <- function(x, v1_end) {
  x %>%
    as_tibble %>% 
    filter(!is.na(date_collected)) %>%
    filter(karenia_mikimotoi >= 0) %>% # -99 in karenia...?
    mutate(datetime_collected=as_datetime(date_collected)) %>%
    filter(datetime_collected >= "2013-07-20") %>% 
    mutate(date=date(datetime_collected),
           hour=pmin(20, pmax(5, hour(datetime_collected))),
           grid=if_else(date_collected < v1_end, 1, 2),
           site.id=as.numeric(factor(paste(sin, area, site)))) %>%
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
    select(-geom, -easting, -northing, -tide, -datetime_collected, -depth)
}



makeFormula <- function(data.df, covs, resp, covs_date=NULL, flist=NULL, sTerms=NULL) {
  library(tidyverse); library(brms); library(glue)
  
  # Exclude NcatF1, NcatF2 if only one level
  exclude_F1 <- n_distinct(data.df$NcatF1)==1
  exclude_F2 <- n_distinct(data.df$NcatF2)==1
  
  if(exclude_F1 & exclude_F2) {
    excl_term <- "catF"
  } else if(exclude_F1 & !exclude_F2) {
    excl_term <- "catF1"
  } else if(!exclude_F1 & exclude_F2) {
    excl_term <- "catF2"
  } else {
    excl_term <- "NO_EXCLUSIONS"
  }
  
  cov_incl <- grep(excl_term, covs, value=T, invert=T)
  
  # interactions between date:env
  if(is.null(flist)) {
    form <- bf(glue("{resp} ~ 1 + {paste(cov_incl, collapse='+')}", 
                    "{ifelse(length(cov_incl) > 0, '+', '')}",
                    "(1 {ifelse(length(cov_incl) > 0, '+', '')}",
                    "{paste(cov_incl, collapse='+')} | siteid)")) 
  } else {
    form <- bf(glue("{resp} ~ 1*bIntercept",
                    "{ifelse(length(sTerms$yday)>0, '+', '')}",
                    "{paste(sTerms$yday, covs_date, sep='*', collapse='+')}",
                    "{ifelse(length(sTerms$b)>0, '+', '')}",
                    "{ifelse(length(sTerms$p)>0, 
                             paste(sTerms$p, sTerms$b, covs, sep='*', collapse='+'),
                             paste(sTerms$b, covs, sep='*', collapse='+'))}"),
               flist=flist, nl=T)
  }
  return(form)
}





calc_ord_mnpr <- function(pred.ord, bloomThresh) {
  tibble(pr0=c(colMeans(pred.ord[,,1,drop=F])),
         pr1=c(colMeans(pred.ord[,,2,drop=F])),
         pr2=c(colMeans(pred.ord[,,3,drop=F])),
         pr3=c(colMeans(pred.ord[,,4,drop=F])),
         obsid=1:dim(pred.ord)[2]) %>%
    pivot_longer(1:4, names_to="ord_cat", values_to="ord_pr") %>%
    mutate(ord_cat=as.numeric(str_sub(ord_cat, -1L, -1L))) %>%
    filter(ord_cat >= bloomThresh) %>%
    group_by(obsid) %>%
    summarise(ord_mnpr=sum(ord_pr),
              across(everything(), ~last(.x))) %>%
    select(ord_mnpr) %>% 
    as_vector
}







best_xgb <- function(df_train, max_depth=2:5, eta=seq(0.01, 2, by=0.01),
                     nthread=1, folds=NULL) {
  library(tidyverse); library(xgboost)
  xg_grid <- expand_grid(max_depth=max_depth, eta=eta) %>%
    mutate(LL_test=NA)
  for (j in 1:nrow(xg_grid)) {
    set.seed(123)
    m_xgb_untuned <- xgb.cv(
      data=as.matrix(df_train[,-1]),
      label=as.numeric(df_train[,1])-1,
      nrounds=1000,
      objective="binary:logistic",
      early_stopping_rounds=4,
      nfold=5,
      folds=folds,
      max_depth=xg_grid$max_depth[j],
      eta=xg_grid$eta[j],
      verbose=F,
      nthread=nthread
    )
    xg_grid$LL_test[j] <- m_xgb_untuned$evaluation_log$test_logloss_mean[m_xgb_untuned$best_iteration]
  }  
  best.par <- filter(xg_grid, LL_test==min(LL_test))
  out <- xgboost(data=as.matrix(df_train[,-1]), label=as.numeric(df_train[,1])-1, 
                 max.depth=best.par$max_depth, eta=best.par$eta, 
                 nthread=nthread, nrounds=10, objective="binary:logistic")
  return(out)
}





# modified from astsa::trend
detrend_loess <- function (x, y, span=0.75, robust=TRUE) {
  if(length(y) < 30 | diff(range(y))==0) {
    return(y - mean(y))
  }
  if(span=="adapt") {
    span <- 0.1
    if(length(y) < 100) span <- 0.2
    if(length(y) < 40) span <- 0.4
  }
  fam = ifelse(robust, "symmetric", "gaussian")
  lo = stats::predict(stats::loess(y ~ x, span=span, family=fam), se = F)
  return(c(y - lo))
}





find_open_bearing <- function(site.df, mesh.fp, buffer=10e3, nDir=24) {
  library(tidyverse); library(sf)
  
  hub.df <- site.df %>% 
    select(site.id, lon, lat) %>%
    st_as_sf(coords=c("lon", "lat"), crs=27700)
  spoke.df <- hub.df %>%
    st_buffer(dist=buffer, nQuadSegs=nDir/4) %>%
    st_cast("POINT") %>%
    group_by(site.id) %>%
    mutate(spoke.id=row_number()) %>%
    filter(spoke.id <= nDir) 
  hubRep.df <- full_join(hub.df, spoke.df %>% st_drop_geometry())
  coords <- cbind(st_coordinates(hubRep.df), st_coordinates(spoke.df))
  spoke.lines <- lapply(1:nrow(coords),
                        function(i){
                          st_linestring(matrix(coords[i,], ncol=2, byrow=TRUE))
                        }) %>%
    st_sfc() %>% st_as_sf() %>% st_set_crs(27700) %>%
    rename(geometry=x) %>%
    mutate(site.id=spoke.df$site.id,
           spoke.id=spoke.df$spoke.id)
  spoke.mesh <- st_intersection(spoke.lines, mesh.fp) %>%
    st_cast("LINESTRING") %>%
    mutate(len=round(as.numeric(st_length(.)))) %>%
    group_by(site.id) %>%
    filter(len==max(len)) %>%
    ungroup %>%
    st_cast("POINT")
  bearings <- as.numeric(lwgeom::st_geod_azimuth(st_transform(spoke.mesh, 4326)))
  bearing.df <- spoke.mesh[1:nrow(spoke.mesh) %% 2 == 0,] %>%
    st_drop_geometry() %>%
    mutate(bearing=bearings[1:length(bearings) %% 2 == 1]) %>%
    group_by(site.id) %>%
    summarise(openBearing=median(bearing)) %>%
    ungroup
  return(full_join(site.df, bearing.df, by="site.id"))
}
