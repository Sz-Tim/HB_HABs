# HABReports Bayesian modelling
# GlobalHAB Workshop 2022 May
# Poster figures
# Tim Szewczyk







# set up ------------------------------------------------------------------
pkgs <- c("tidyverse", "glue", "lubridate", "sf", "WeStCOMS")
invisible(lapply(pkgs, library, character.only=T))
theme_set(theme_classic() + theme(panel.grid.minor=element_blank()))
source("code/000_fnMisc.R")

mesh.f <- paste0("C:\\Users\\sa04ts\\OneDrive - SAMS\\Projects\\WeStCOMS\\data\\",
                 c("WeStCOMS_mesh.gpkg", "WeStCOMS2_mesh.gpkg"))
mesh.sf <- map(mesh.f, loadMesh)

gis.dir <- "..\\..\\00_gis\\"
coast <- read_sf(paste0(gis.dir, "coastline\\ne_10m_land.shp")) %>%
  add_row(read_sf(paste0(gis.dir, "coastline\\ne_10m_minor_islands.shp")))

sampling.df <- read_csv("data\\sampling_local.csv") 
sites.sf <- sampling.df %>% 
  group_by(site.id) %>%
  summarise(lon=median(lon), lat=median(lat)) %>%
  ungroup %>%
  st_as_sf(coords=c("lon", "lat"), crs=27700)

thresh.df <- read_csv("data\\hab_tf_thresholds.csv") %>%
  filter(!is.na(tl)) %>%
  group_by(hab_parameter, tl) %>%
  slice_head(n=1) %>%
  ungroup

pred.df <- dir("out", "pred_.*csv", full.names=T) %>% 
  map(~read_csv(.x, show_col_types=F)) %>% bind_rows

poster_theme <- theme_classic() + 
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=16),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16),
        plot.title=element_text(size=20),
        strip.text=element_text(size=18))

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
                   z=c(0.5, 0.7, 1.2, NA),
                   zdens=dnorm(z)) %>%
  mutate(tl2=lead(tl))
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
  labs(x=expression(italic(tilde(y))~~plain("(latent)")), y="density") +
  theme_classic() + poster_theme
ggsave("figs\\poster\\distr_ordinal.png", height=5, width=7, dpi=300)  



y_i <- tibble(x_min=c(-2, 1),
              x_max=c(-1, 2),
              x=c(-1.5, 1.5),
              y=c(-0.04, -0.04))

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
  geom_point(data=y_i, aes(x, y), size=3) +
  geom_errorbarh(data=y_i, aes(xmin=x_min, xmax=x_max, y=y), height=0.01) +
  scale_fill_manual(values=tl.col, guide="none") +
  labs(x=expression(italic(tilde(y))~~plain("(latent)")), y="density") +
  theme_classic() + poster_theme
ggsave("figs\\poster\\distr_ordinal_pts.png", height=5, width=7, dpi=300) 

y_i <- tribble(
  ~x_min, ~x_25, ~x, ~x_75, ~x_max, ~y,
  -1, -0.95, -0.8, -0.65, 0.6, -0.08,
  1, 1.2, 1.5, 1.7, 2, -0.08,
  -2.1, -1.9, -1.75, -1.6, -1.5, -0.08
)
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
  geom_point(data=y_i, aes(x, y), size=3) +
  geom_linerange(data=y_i, aes(xmin=x_25, xmax=x_75, y=y), size=1.5) +
  geom_errorbarh(data=y_i, aes(xmin=x_min, xmax=x_max, y=y), height=0.01) +
  scale_fill_manual(values=tl.col, guide="none") +
  labs(x=expression(italic(tilde(y))~~plain("(latent)")), y="density") +
  theme_classic() + poster_theme + coord_flip()
ggsave("figs\\poster\\distr_ordinal_pts_flip.png", height=7*.8, width=5*.8, dpi=300) 




y_tl <- bind_rows(
  tibble(tl=rep(unique(tl.examp$tl), each=4e3),
         pr=c(rnorm(4e3, 0.05, 0.02), rnorm(4e3, 0.2, 0.03), rnorm(4e3, 0.3, 0.05), rnorm(4e3, 0.6, 0.05)),
         pt="italic(tilde(y)[1][t])"),
  tibble(tl=rep(unique(tl.examp$tl), each=4e3),
         pr=c(rnorm(4e3, 0.8, 0.05), rnorm(4e3, 0.1, 0.04), rnorm(4e3, 0.03, 0.02), rnorm(4e3, 0.02, 0.02)),
         pt="italic(tilde(y)[2][t])"),
  tibble(tl=rep(unique(tl.examp$tl), each=4e3),
         pr=c(rnorm(4e3, 0.95, 0.02), rnorm(4e3, 0.025, 0.01), rnorm(4e3, 0.02, 0.01), rnorm(4e3, 0.015, 0.01)),
         pt="italic(tilde(y)[3][t])")
) %>% 
  mutate(pr=pmax(pmin(pr, 1), 0))

ggplot(y_tl, aes(pr, fill=tl)) + geom_density(alpha=0.5) +
  scale_fill_manual(values=tl.col, guide="none") +
  facet_wrap(~pt, labeller="label_parsed", scales="free_y", ncol=1) +
  theme_classic() + poster_theme +
  scale_x_continuous(breaks=c(0, 0.5, 1)) +
  labs(x="Forecast probability")
ggsave("figs\\poster\\sim_ordinal_pred.png", height=7, width=3, dpi=300)







# forecast barplots -------------------------------------------------------

predPrmax.df <- pred.df %>% 
  select(obs.id, sp, contains("ord_prmax"), N.catF) %>%
  pivot_longer(ends_with(c("0", "1", "2", "3")), 
               names_to="predCat", values_to="prob") %>%
  mutate(predCat=str_sub(predCat, -1, -1)) %>%
  mutate(predCat=factor(predCat, labels=sort(unique(N.catF)), ordered=T),
         N.catF=factor(N.catF, labels=sort(unique(N.catF)), ordered=T)) %>%
  group_by(sp, obs.id) %>%
  arrange(desc(prob)) %>%
  slice_head(n=1) %>%
  ungroup 

sp.rho <- predPrmax.df %>%
  group_by(sp) %>%
  summarise(r=cor(as.numeric(N.catF), as.numeric(predCat), method="spearman"))

pred.ls <- predPrmax.df %>% 
  group_by(sp, N.catF, predCat) %>% 
  summarise(N=n()) %>% ungroup %>%
  pivot_wider(names_from=predCat, values_from=N, names_prefix="pred_", values_fill=0) %>%
  group_by(sp) %>%
  group_split() 
sp.somersDelta <- map_dbl(pred.ls, ~.x %>% select(3:6) %>% as.matrix() %>% as.table() %>%
           SomersDelta(., direction="row"))

sp.df <- tibble(sp=c("alexandrium_sp", "dinophysis_sp", "karenia_mikimotoi",
                     "prorocentrum_lima", "pseudo_nitzschia_sp"),
                sp_clean=c("Alexandrium sp.", "Dinophysis sp.",
                           "Karenia mikimotoi", "Prorocentrum lima",
                           "Pseudo-nitzschia sp."),
                sp_expr=c("italic(Alexandrium~sp.)", "italic(Dinophysis~sp.)",
                           "italic(Karenia~mikimotoi)", "italic(Prorocentrum~lima)",
                           "italic(Pseudonitzschia~sp.)"),
                sp_rho=paste0(sp_clean, ", r: ", str_pad(round(sp.rho$r, 2), 4, "r", "0")),
                sp_Dxy=paste0(sp_clean, ", Dxy: ", str_pad(round(sp.somersDelta, 2), 4, "r", "0")),
                sp_rho_expr=paste0("list(", sp_expr, ",rho==", str_pad(round(sp.rho$r, 2), 4, "r", "0"), ")"),
                sp_Dxy_expr=paste0("list(", sp_expr, ",D==", str_pad(round(sp.somersDelta, 2), 4, "r", "0"), ")"))

predPrmax.df %>%
  full_join(., sp.df) %>%
  ggplot(aes(predCat, fill=N.catF)) + geom_bar(position="fill", colour="grey30") +
  scale_fill_manual("Observed", values=tl.col, labels=c("Absent", "Low", "Med.", "High")) +
  scale_x_discrete(labels=c("Absent", "Low", "Med.", "High")) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) +
  labs(title="Forecasting performance",
       x="Forecasted category", y="Proportion of observations") + 
  facet_wrap(~sp_Dxy_expr, ncol=2, labeller=label_parsed) + 
  poster_theme + 
  theme(legend.position=c(0.775, 0.15), 
        strip.text=element_text(face="italic"),
        axis.text=element_text(size=14),
        strip.text.x=element_text(size=12))
ggsave("figs\\poster\\forecast_proportions.png", height=10, width=6, dpi=300)

predPrmax.df %>%
  full_join(., sp.df) %>%
  ggplot(aes(predCat, fill=N.catF)) + geom_bar(position="fill", colour="grey30") +
  scale_fill_manual("Observed", values=tl.col, labels=c("Absent", "Low", "Med.", "High")) +
  scale_x_discrete(labels=c("Absent", "Low", "Med.", "High")) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) +
  labs(title="Forecasting performance",
       x="Forecasted category", y="Proportion of observations") + 
  facet_wrap(~sp_Dxy_expr, ncol=3, labeller=label_parsed) + 
  poster_theme + 
  theme(legend.position=c(0.775, 0.15), 
        strip.text=element_text(face="italic"),
        axis.text=element_text(size=14),
        strip.text.x=element_text(size=12))
ggsave("figs\\poster\\forecast_proportions_wide.png", height=6, width=9, dpi=300)

predPrmax.df %>%
  full_join(., sp.df) %>%
  ggplot(aes(N.catF, fill=predCat)) + geom_bar(position="fill", colour="grey30") +
  scale_fill_manual("Forecasted", values=tl.col, labels=c("Absent", "Low", "Med.", "High")) +
  scale_x_discrete(labels=c("Absent", "Low", "Med.", "High")) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) +
  labs(title="Forecasting performance",
       x="Observed category", y="Proportion of observations") + 
  facet_wrap(~sp_clean, ncol=2) + 
  poster_theme + 
  theme(legend.position=c(0.775, 0.15), 
        strip.text=element_text(face="italic"),
        axis.text=element_text(size=14),
        strip.text.x=element_text(size=16))
ggsave("figs\\poster\\observed_proportions.png", height=10, width=5.5, dpi=300)


predPrmax.df %>% 
  group_by(sp, predCat) %>%
  summarise(prCorrect=mean(N.catF==predCat)) %>%
  arrange(desc(predCat), sp)





# observed distribution ---------------------------------------------------

obs.df <- sampling.df %>% 
  select(site.id, date, any_of(sp.df$sp)) %>%
  filter(year(date) <= 2019) %>%
  pivot_longer(any_of(sp.df$sp), names_to="sp", values_to="N") %>%
  mutate(N.ln=log(N+1)) %>%
  rowwise() %>%
  mutate(N.cat=filter(thresh.df, hab_parameter==sp)$tl[max(which(N >= filter(thresh.df, hab_parameter==sp)$min_ge))])

obs.df %>%
  full_join(., sp.df) %>%
  ggplot(aes(N.ln, fill=N.cat)) + 
  geom_histogram(colour="grey30", size=0.5, bins=15) +
  scale_fill_manual(values=tl.col, guide="none") +
  poster_theme + 
  theme(strip.text=element_text(face="italic"),
        axis.text=element_text(size=14),
        strip.text.y=element_text(size=12)) +
  labs(title="Observed distributions",
       x="ln(cells/L)") +
  facet_wrap(~sp_clean, ncol=1, scales="free", strip.position="right")
ggsave("figs\\poster\\observed_data_distr.png", height=11, width=4, dpi=300)

obs.df %>%
  group_by(sp, N.cat) %>%
  summarise(nObs=n()) %>%
  group_by(sp) %>%
  mutate(prObs=nObs/sum(nObs)) %>%
  full_join(., sp.df) %>%
  ggplot(aes(x="", y=prObs, fill=N.cat)) + 
  geom_bar(stat="identity", colour="grey30") +
  scale_fill_manual(values=tl.col, guide="none") +
  theme(strip.text=element_text(face="italic"),
        axis.text=element_blank(),
        strip.text.y=element_text(size=12),
        axis.title=element_blank()) +
  facet_wrap(~sp_clean) +
  coord_polar(theta="y")
ggsave("figs\\poster\\observed_data_pie.png", height=4, width=5, dpi=300)
