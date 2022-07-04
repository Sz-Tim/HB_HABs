


pkgs <- c("tidyverse", "brms", "glue", "lubridate", "sf", "WeStCOMS")
invisible(lapply(pkgs, library, character.only=T))
theme_set(theme_classic() + theme(panel.grid.minor=element_blank()))
source("code/000_fnMisc.R")
if(.Platform$OS.type=="unix") {
  sep <- "/"
  westcoms.dir <- paste0("/media/archiver/common/sa01da-work/",
                         c("minch2/Archive/", "WeStCOMS2/Archive/"))
  mesh.f <- paste0("/home/sa04ts/FVCOM_meshes/",
                   c("WeStCOMS_mesh.gpkg", "WeStCOMS2_mesh.gpkg"))
} else {
  sep <- "\\"
  westcoms.dir <- paste0("D:\\hydroOut\\", 
                         c("minch2\\Archive\\", "WeStCOMS2\\Archive\\"))
  mesh.f <- paste0("C:\\Users\\sa04ts\\OneDrive - SAMS\\Projects\\WeStCOMS\\data\\",
                   c("WeStCOMS_mesh.gpkg", "WeStCOMS2_mesh.gpkg"))
}
mesh.sf <- map(mesh.f, loadMesh)



# make site files ---------------------------------------------------------

# Calculate shortest distance between sites:
# https://gis.stackexchange.com/questions/305070/finding-shortest-path-without-barriers-using-the-sf-package-in-r

mesh.sf[[1]] %>% 
  filter(area > exp(14)) %>% 
  st_centroid() %>% 
  mutate(lon=st_coordinates(.)[,1], lat=st_coordinates(.)[,2]) %>% 
  st_drop_geometry %>% 
  select(i, lon, lat) %>% 
  sample_n(size=5e3) %>% 
  write_tsv("particle_tracking/out/large_elem_centroids_v1.dat", col_names=F)
mesh.sf[[2]] %>% 
  filter(area > exp(14)) %>% 
  st_centroid() %>% 
  mutate(lon=st_coordinates(.)[,1], lat=st_coordinates(.)[,2]) %>% 
  st_drop_geometry %>% 
  select(i, lon, lat) %>% 
  sample_n(size=5e3) %>% 
  write_tsv("particle_tracking/out/large_elem_centroids_v2.dat", col_names=F)

mesh.sf[[1]] %>% 
  st_centroid() %>% 
  mutate(lon=st_coordinates(.)[,1], lat=st_coordinates(.)[,2]) %>%
  st_drop_geometry %>% 
  select(i, lon, lat) %>% 
  sample_n(size=5e3) %>% 
  write_tsv("particle_tracking/out/elem_centroids_v1.dat", col_names=F)
mesh.sf[[2]] %>% 
  st_centroid() %>% 
  mutate(lon=st_coordinates(.)[,1], lat=st_coordinates(.)[,2]) %>% 
  st_drop_geometry %>% 
  select(i, lon, lat) %>% 
  sample_n(size=5e3) %>% 
  write_tsv("particle_tracking/out/elem_centroids_v2.dat", col_names=F)







# plots -------------------------------------------------------------------

gis.dir <- "..\\..\\00_gis\\"
coast <- read_sf(paste0(gis.dir, "coastline\\ne_10m_land.shp")) %>%
  add_row(read_sf(paste0(gis.dir, "coastline\\ne_10m_minor_islands.shp"))) %>%
  st_crop(st_bbox(c(xmin=-7.75, xmax=-4.25, ymin=54.5, ymax=58.5), crs=st_crs(4326)))

out.dir <- "particle_tracking\\out\\"

fsa.df <- read_delim(paste0(out.dir, "fsa_sites_v1.dat"), 
                     delim="\t", col_names=c("site.id", "x", "y")) %>%
  filter(!site.id %in% c(70, 74, 75, 80, 88))




# shortest paths ----------------------------------------------------------

library(raster); library(gdistance); library(rSDM)
# need to turn the mesh into a raster into a transition layer 
# https://agrdatasci.github.io/gdistance/reference/index.html
mesh.full <- raster(str_replace(mesh.f[1], ".gpkg", ".tif"))
mesh.r <- raster(str_replace(mesh.f[1], ".gpkg", "Footprint.tif"))
mesh.tr <- transition(mesh.r, mean, 16)
mesh.tr <- geoCorrection(mesh.tr)
saveRDS(mesh.tr, glue("data{sep}WeStCOMS_transitionMx.rds"))

set_ll_warn(TRUE)
fsa.moved <- SpatialPointsDataFrame(fsa.df[,2:3], data=fsa.df[,1], 
                       proj4string=CRS("+init=epsg:27700")) %>%
  points2nearestcell(., mesh.r) %>%
  as.data.frame


paths.ls <- list()
for(i in 1:nrow(fsa.moved)) {
    paths.ls[[i]] <- shortestPath(mesh.tr, 
                                  as.matrix(fsa.moved[i,2:3]),
                                  as.matrix(fsa.moved[-i,2:3]),
                                  output="SpatialLines") %>%
      st_as_sf %>%
      mutate(origin=fsa.moved$site.id[i],
             destination=fsa.moved$site.id[-i],
             distance=st_length(.))
}
dist.df <- do.call(rbind, paths.ls) %>%
  st_drop_geometry()
write_csv(dist.df, glue("data{sep}fsa_site_pairwise_distances.csv"))


par(mfrow=c(5,5))
for(i in 1:nrow(fsa.moved)) {
  plot(mesh.full)
  points(fsa.moved[i,2:3], col="red3", pch=16)
  lines(shortestPath(mesh.tr, 
                     as.matrix(fsa.moved[i,2:3]),
                     as.matrix(fsa.moved[-i,2:3]),
                     output="SpatialLines"))
}


dist.df %>% group_by(origin) %>%
  summarise(proximity=log10(sum(1/as.numeric(distance)))) %>%
  ungroup %>%
  full_join(fsa.df, ., by=c("site.id"="origin")) %>%
  st_as_sf(coords=c("x", "y"), crs=27700) %>%
  ggplot() + 
  geom_sf(data=coast) +
  geom_sf(aes(colour=proximity)) +
  scale_colour_viridis_c()





# connectivity ------------------------------------------------------------

connect.f <- dir(glue("{out.dir}connectivity"))

connect.df <- map_dfr(connect.f, 
                      ~read_delim(glue("{out.dir}connectivity{sep}{.x}"), 
                                  delim=" ", col_names=F, show_col_types=F) %>% 
                        mutate(date=str_sub(.x, 14, 21)) %>%
                        group_by(date) %>% 
                        summarise(across(everything(), sum)) %>% 
                        ungroup) %>%
  pivot_longer(-1, names_to="site.id", values_to="influx") %>%
  mutate(site.id=rep(fsa.df$site.id, length(connect.f)),
         date=ymd(date))
write_csv(connect.df, glue("data{sep}connectivity_5e3parts.csv"))

connect.df

ggplot(connect.df, aes(date, log(influx+1), group=site.id)) + 
  stat_smooth(method="loess", se=F, span=0.1) + facet_wrap(~site.id)

connect.df %>%
  group_by(site.id) %>%
  mutate(rel_influx=influx/max(1, influx)) %>%
  ggplot(aes(date, rel_influx, group=site.id)) + 
  stat_smooth(method="loess", se=F, span=0.1) + facet_wrap(~site.id)

connect.df %>%
  group_by(site.id) %>%
  mutate(rel_influx=log(influx+1)/max(log(2), log(influx+1))) %>%
  ggplot(aes(date, rel_influx, group=site.id)) + 
  stat_smooth(method="loess", se=F, span=0.1) + facet_wrap(~site.id)



loc.df <- read_delim(dir(glue("{out.dir}locs"), "locationsEnd", full.names=T)) %>%
  arrange(ID, age)
loc.sf <- loc.df %>% 
  st_as_sf(coords=c("x", "y"), crs=27700) %>%
  group_by(ID) %>% 
  summarise(ID=first(ID), 
            startLocation=first(startLocation),
            do_union=F) %>%
  st_cast("LINESTRING")


ggplot() + 
  geom_sf(data=coast) +
  geom_sf(data=loc.sf, alpha=0.5, colour="blue3") +
  theme_classic() +
  scale_colour_viridis_c()

