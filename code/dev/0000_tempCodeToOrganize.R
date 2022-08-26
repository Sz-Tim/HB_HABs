# CODE THAT NEEDS TO BE ORGANIZED




# wind/water directions ---------------------------------------------------

expand_grid(x=seq(-1,1,length.out=20), 
            y=seq(-1,1,length.out=20)) %>% 
  mutate(dir=atan2(x,y)) %>% 
  mutate(Nshore=pi, 
         Sshore=0) %>% 
  ggplot(aes(x,y,colour=cos(dir-Sshore))) + 
  geom_point() + scale_colour_distiller(type="div")
# With this parameterisation, NEG = TOWARD shore, POS = AWAY from shore

# get fetch ---------------------------------------------------------------

library(raster); library(sf)
fetch <- raster("..\\..\\00_gis\\fetch\\log10_UK200m_depth40_0na.tif")
sites.sf <- read_csv("data\\sampling_local.csv") %>% group_by(site.id) %>% 
  summarise(lon=first(lon), lat=first(lat)) %>%
  st_as_sf(coords=c("lon", "lat"), crs=27700) %>%
  mutate(fetch=raster::extract(fetch, ., fun=mean, small=T))
read_csv("data\\sampling_local.csv") %>% left_join(st_drop_geometry(sites.sf)) %>%
  write_csv("data\\sampling_local.csv")

sites.sf <- read_csv("data\\sampling_regional.csv") %>% group_by(site.id) %>% 
  summarise(lon=first(lon), lat=first(lat)) %>%
  st_as_sf(coords=c("lon", "lat"), crs=27700) %>%
  mutate(fetch=raster::extract(fetch, ., fun=mean, small=T))
read_csv("data\\sampling_regional.csv") %>% left_join(st_drop_geometry(sites.sf)) %>%
  write_csv("data\\sampling_regional.csv")





library(tidyverse); library(sf)

mesh <- st_read("..\\..\\WeStCOMS\\data\\WeStCOMS2_meshOutline.gpkg")
fsa.pt <- st_read("data\\fsa_WeStCOMS2_lineID.gpkg")

mesh %>% filter(line_id %in% unique(fsa.pt$join_line_id)) %>%
  st_write("data\\WeStCOMS2_fsaLines.gpkg")

# st_read("data\\WeStCOMS2_fsaPts.gpkg") %>%
#   group_by(line_id) %>% arrange(distance) %>%
#   slice_head(n=1) %>%
#   ungroup %>%
#   st_write("data\\WeStCOMS2_fsaPts.gpkg", append=F)

mesh.pt <- st_read("data\\WeStCOMS2_fsaPts.gpkg")

site.pts <- fsa.pt %>% group_by(site.id) %>%
  slice_head(n=1) %>% ungroup %>%
  select(site.id, join_line_id)



bind_rows(mesh.pt %>% 
            select(-distance, -join_fid, -distance) %>%
            mutate(source="mesh"),
          site.pts %>% mutate(source="fsa") %>% rename(line_id=join_line_id)) %>%
  arrange(line_id, source) %>%
  group_by(line_id) %>%
  summarise(site.id=first(site.id)) %>%
  st_cast("LINESTRING") %>%
  st_transform(4326) %>%
  mutate(bearing=stplanr::line_bearing(.)) %>%
  st_drop_geometry %>%
  select(site.id, bearing) %>%
  write_csv("data\\site_coast_bearing.csv")
  # st_write("data\\fsa_to_mesh_lines.gpkg", append=F)
  

left_join(read_csv("data\\sampling_local.csv"),
          read_csv("data\\site_coast_bearing.csv"), 
          by="site.id") %>%
  write_csv("data\\sampling_local.csv")

left_join(read_csv("data\\sampling_regional.csv"),
          read_csv("data\\site_coast_bearing.csv"), 
          by="site.id") %>%
  write_csv("data\\sampling_regional.csv")





# https://github.com/AndriSignorell/DescTools

ConDisPairs <-function(x){
  
  # tab is a matrix of counts
  # Based on code of Michael Friendly and Laura Thompson
  
  # slooooow because of 2 nested for clauses O(n^2)
  # this is NOT faster when implemented with a mapply(...)
  
  # Lookin for alternatives in C
  # http://en.verysource.com/code/1169955_1/kendl2.cpp.html
  # cor(..., "kendall") is for dimensions better
  
  n <- nrow(x)
  m <- ncol(x)
  pi.c <- pi.d <- matrix(0, nrow = n, ncol = m)
  
  row.x <- row(x)
  col.x <- col(x)
  
  for(i in 1:n){
    for(j in 1:m){
      pi.c[i, j] <- sum(x[row.x<i & col.x<j]) + sum(x[row.x>i & col.x>j])
      pi.d[i, j] <- sum(x[row.x<i & col.x>j]) + sum(x[row.x>i & col.x<j])
    }
  }
  C <- sum(pi.c * x)/2
  D <- sum(pi.d * x)/2
  
  return(list(pi.c = pi.c, pi.d = pi.d, C = C, D = D))
  
}

SomersDelta <- function(x,  y = NULL, direction=c("row","column"), conf.level = NA, ...) {
  
  if(!is.null(y)) tab <- table(x, y, ...)
  else tab <- as.table(x)
  
  # tab is a matrix of counts
  x <- ConDisPairs(tab)
  
  # use .DoCount
  #   if(is.na(conf.level)) {
  #     d.tab <- as.data.frame.table(tab)
  #     x <- .DoCount(d.tab[,1], d.tab[,2], d.tab[,3])
  #   } else {
  #     x <- ConDisPairs(tab)
  #   }
  
  m <- min(dim(tab))
  n <- sum(tab)
  switch( match.arg( arg = direction, choices = c("row","column") )
          , "row" = { ni. <- colSums(tab) }
          , "column" = { ni. <- rowSums(tab) }
  )
  wt <- n^2 - sum(ni.^2)
  # Asymptotic standard error: sqrt(sigma2)
  sigma2 <- 4/wt^4 * (sum(tab * (wt*(x$pi.c - x$pi.d) - 2*(x$C-x$D)*(n-ni.))^2))
  # debug: print(sqrt(sigma2))
  
  somers <- (x$C - x$D) / (n * (n-1) /2 - sum(ni. * (ni. - 1) /2 ))
  
  pr2 <- 1 - (1 - conf.level)/2
  ci <- qnorm(pr2) * sqrt(sigma2) * c(-1, 1) + somers
  
  if(is.na(conf.level)){
    result <- somers
  } else {
    result <- c(somers = somers,  lwr.ci=max(ci[1], -1), upr.ci=min(ci[2], 1))
  }
  
  return(result)
  
}










# ML ----------------------------------------------------------------------

# https://www.r-bloggers.com/2021/04/random-forest-in-r/
library(randomForest); library(caret)

for(sp in 1:length(species)) {
  target <- species[sp]
  target.df <- read_csv(glue("out{sep}dataset_{target}.csv"))
  
  train.rf <- target.df %>% filter(year <= 2017) %>%
    select(N.bloom, 
           site.id, lon, lat, 
           ydaySin, ydayCos, fetch,
           N.ln_1, N.cat_1, 
           N.ln_2, N.cat_2, 
           N.lnWt_1, N.lnWt_2, 
           one_of(covars)) %>%
    mutate(N.bloom=factor(N.bloom)) %>%
    as.data.frame()
  test.rf <- target.df %>% filter(year > 2017) %>%
    select(N.bloom, 
           site.id, lon, lat, 
           ydaySin, ydayCos, fetch,
           N.ln_1, N.cat_1, 
           N.ln_2, N.cat_2, 
           N.lnWt_1, N.lnWt_2, 
           one_of(covars)) %>%
    mutate(N.bloom=factor(N.bloom)) %>%
    as.data.frame()
  rf <- randomForest(N.bloom ~ ., data=train.rf, 
                     xtest=test.rf[,-1], ytest=test.rf[,1], 
                     proximity=T, importance=T, keep.forest=T)
  p.rf <- predict(rf, test.rf, type="prob")
  
}



hist(treesize(rf),
     main = "No. of Nodes for the Trees",
     col = "green")
varImpPlot(rf,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")
importance(rf)

partialPlot(rf, train.rf, short_wave_L_wk, "1")
partialPlot(rf, train.rf, o2_wk, "1")
partialPlot(rf, train.rf, wind_L_wk, "1")
partialPlot(rf, train.rf, attn_wk, "1")
partialPlot(rf, train.rf, po4_wk, "1")
