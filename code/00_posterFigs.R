# HABReports Bayesian modelling
# GlobalHAB Workshop 2022 May
# Poster figures
# Tim Szewczyk







# set up ------------------------------------------------------------------
pkgs <- c("tidyverse", "glue", "lubridate", "sf", "WeStCOMS")
invisible(lapply(pkgs, library, character.only=T))
theme_set(theme_classic() + theme(panel.grid.minor=element_blank()))

mesh.f <- paste0("C:\\Users\\sa04ts\\OneDrive - SAMS\\Projects\\WeStCOMS\\data\\",
                 c("WeStCOMS_mesh.gpkg", "WeStCOMS2_mesh.gpkg"))
mesh.sf <- map(mesh.f, loadMesh)

gis.dir <- "..\\..\\00_gis\\"
coast <- read_sf(paste0(gis.dir, "coastline\\ne_10m_land.shp")) %>%
  add_row(read_sf(paste0(gis.dir, "coastline\\ne_10m_minor_islands.shp")))

sampling.df <- read_csv("data\\obs_df.csv") %>% filter(grid=="minch2") %>% filter(!is.na(site.elem))
sites.sf <- sampling.df %>% 
  group_by(site.id) %>%
  summarise(N=n(), nYr=n_distinct(year), 
            lon=median(lon), lat=median(lat)) %>%
  ungroup %>%
  st_as_sf(coords=c("lon", "lat"), crs=27700)

poster_theme <- theme_classic() + 
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16),
        plot.title=element_text(size=20))

tl.col <- c("green3", "gold1", "orange", "red")





# map of data -------------------------------------------------------------
p <- ggplot() +
  geom_sf(data=coast, fill="#1E525E", colour=NA) +
  geom_sf(data=mesh.sf[[1]], fill=NA, colour="#46BFDE", size=0.1) +
  geom_sf(data=sites.sf, colour="#E77334", size=3) +
  scale_x_continuous(limits=c(-7.75, -4.25), breaks=c(-7, -5)) +
  scale_y_continuous(limits=c(54.5, 58.5), breaks=c(55, 57)) +
  coord_sf(label_axes="-NE-") + poster_theme
ggsave("figs\\poster\\map.png", p, height=8, width=4.25, dpi=300)  






# distributions -----------------------------------------------------------

tl.examp <- tibble(tl=unique(thresh.df$tl),
                   val_min=c(-0.25, 0.25, 5, 7),
                   val=c(0, 0.25, 5, 7),
                   val_max=c(0.25, 5, 7, 21),
                   ymin=0, 
                   ymax=650,
                   z=c(0.5, 0.9, 1.5, NA),
                   zdens=dnorm(z))
bind_rows(tibble(distr="Frechet: >0",
                 value=rfrechet(3e3, scale=2.5, shape=2)),
          tibble(distr="Hurdle: Pr(0)",
                 value=rep(0, 6e2))) %>%
  mutate(distr=factor(distr, levels=rev(unique(distr)))) %>%
  ggplot() + 
  geom_segment(data=tl.examp, aes(colour=tl, x=val, xend=val, y=ymin, yend=ymax),
               size=2) +
  geom_histogram(aes(value, fill=distr),
                 breaks=seq(-0.25,20.25,by=0.5), colour="black", size=0.25) +
  xlim(-1,21) + 
  scale_fill_manual("Distribution", values=c("#D21E4C", "#2A5883")) +
  scale_colour_manual(values=tl.col, guide="none") +
  labs(x="ln(Phytoplankton)", y="count") +
  theme_classic() + poster_theme +
  theme(legend.position=c(0.8, 0.8))
ggsave("figs\\poster\\distr_frechet.png", height=5, width=7, dpi=300)  




tibble(x=seq(-3,3,length.out=1e4),
       y=dnorm(x),
       tl=case_when(x<tl.examp$z[1] ~ "0_green",
                    x>=tl.examp$z[1] & x<tl.examp$z[2] ~ "1_yellow",
                    x>=tl.examp$z[2] & x<tl.examp$z[3] ~ "2_orange",
                    x>=tl.examp$z[3] ~ "3_red")) %>%
  ggplot() + 
  geom_ribbon(aes(x, ymin=0, ymax=y, fill=tl, group=tl),
              colour="black", size=0.25) +
  geom_segment(data=tl.examp, aes(x=z, xend=z, y=0, yend=zdens),
               colour="black", size=0.25) +
  scale_fill_manual(values=tl.col, guide="none") +
  labs(x="Latent variable", y="density") +
  theme_classic() + poster_theme
ggsave("figs\\poster\\distr_ordinal.png", height=5, width=7, dpi=300)  

