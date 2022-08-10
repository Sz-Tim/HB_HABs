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


calc_ord_mnpr <- function(pred.ord, bloomThresh) {
  tibble(pr0=c(colMeans(pred.ord[,,1,drop=F])),
         pr1=c(colMeans(pred.ord[,,2,drop=F])),
         pr2=c(colMeans(pred.ord[,,3,drop=F])),
         pr3=c(colMeans(pred.ord[,,4,drop=F])),
         obsid=pred.df$obsid) %>%
    pivot_longer(1:4, names_to="ord_cat", values_to="ord_pr") %>%
    mutate(ord_cat=as.numeric(str_sub(ord_cat, -1L, -1L))) %>%
    filter(ord_cat >= bloomThresh) %>%
    group_by(obsid) %>%
    summarise(ord_mnpr=sum(ord_pr),
              across(everything(), ~last(.x))) %>%
    select(ord_mnpr) %>% 
    as_vector
}


