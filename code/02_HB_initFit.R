# HABReports Bayesian modelling
# Bayesian HAB model initial fits
# Tim Szewczyk


# This script fits models with data from 2013-2017 using horseshoe priors



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "LaplacesDemon", "brms", "bayesplot", "WeStCOMS")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
walk(dir("code", "000_fn", full.names=T), source)

# minch2:    2013-06-20 to 2019-07-02
# WeStCOMS2: 2019-04-01 to 2022-01-26
if(.Platform$OS.type=="unix") {
  sep <- "/"
  westcoms.dir <- c("/media/archiver/common/sa01da-work/minch2/Archive/",
                    "/media/archiver/common/sa01da-work/WeStCOMS2/Archive/")
  mesh.f <- c("/home/sa04ts/FVCOM_meshes/WeStCOMS_mesh.gpkg",
              "/home/sa04ts/FVCOM_meshes/WeStCOMS2_mesh.gpkg")
} else {
  sep <- "\\"
  westcoms.dir <- c("D:\\hydroOut\\minch2\\Archive\\",
                    "D:\\hydroOut\\WestCOMS2\\Archive\\")
  mesh.f <- c("..\\..\\01_FVCOM\\data\\WeStCOMS_Mesh.gpkg",
              "..\\..\\01_FVCOM\\data\\WeStCOMS2_Mesh.gpkg")
}





# covariates --------------------------------------------------------------

species <- c("alexandrium_sp", "dinophysis_sp", "karenia_mikimotoi",
             "prorocentrum_lima", "pseudo_nitzschia_sp")

sampling.df <- read_csv(glue("data{sep}sampling_local.csv"))

thresh.df <- read_csv(glue("data{sep}hab_tf_thresholds.csv")) %>%
  filter(!is.na(tl)) %>%
  group_by(hab_parameter, tl) %>%
  slice_head(n=1) %>%
  ungroup

hydro.df <- dir("data", "hydro_", full.names=T) %>% 
  map(~read_csv(.x, show_col_types=F)) %>%
  reduce(., full_join) %>% 
  pivot_longer(c(contains("_L"), contains("_R")), names_to="param", values_to="value") %>%
  mutate(parType=str_sub(param, 1, -4),
         res=str_sub(param, -2, -2),
         lag=str_sub(param, -1, -1)) %>%
  select(-param) %>%
  pivot_wider(names_from="parType", values_from="value") %>%
  mutate(water=sqrt(u^2 + v^2 + ww^2),
         wind=sqrt(uwind_speed^2 + vwind_speed^2),
         waterDir=atan2(v, u),
         windDir=atan2(vwind_speed, uwind_speed),
         timespan=if_else(lag=="0", "0", "wk")) %>% 
  select(-u, -v, -ww, -vwind_speed, -uwind_speed) %>%
  group_by(obs.id, res, timespan) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm=T))) %>%
  ungroup %>%
  pivot_wider(names_from=c(res, timespan), values_from=4:ncol(.),
              names_glue="{.value}_{res}_{timespan}")
connect.df <- read_csv(glue("data{sep}connectivity.csv")) %>%
  group_by(site.id) %>%
  mutate(influx_wk=zoo::rollsum(influx, 7, align='right', fill="extend")) %>%
  mutate(across(starts_with("influx"), ~log(.x+1)))
# >0.85 cor(res,time) = temp, salinity
# >0.85 cor(res) = short_wave, km, precip, wind, windDir
# >0.85 cor(time) = water
# waterDir is more variable across all 4
covars <- c("temp_L_wk", "salinity_L_wk",
            "short_wave_L_0", "short_wave_L_wk",
            "km_L_0", "km_L_wk",
            "precip_L_0", "precip_L_wk",
            "wind_L_0", "wind_L_wk",
            "windDir_L_0", "windDir_L_wk",
            "water_L_wk", "water_R_wk",
            "waterDir_L_0", "waterDir_R_0", "waterDir_L_wk", "waterDir_R_wk",
            "influx", "influx_wk")

predictors_main <- "
  ydayCos + ydaySin + 
                 
  short_wave_L_wk:ydayCos:ydaySin +
  temp_L_wk:ydayCos:ydaySin +
  precip_L_wk:ydayCos:ydaySin +
  salinity_L_wk:ydayCos:ydaySin +
  water_L_wk:ydayCos:ydaySin +
  water_R_wk:ydayCos:ydaySin +
  wind_L_wk:ydayCos:ydaySin +
  windDir_L_wk:ydayCos:ydaySin +
  waterDir_L_wk:ydayCos:ydaySin +
  waterDir_R_wk:ydayCos:ydaySin +
  fetch + influx_wk +
                 
  mo(N.catF_1):ydayCos:ydaySin +
  mo(N.catF_2):ydayCos:ydaySin +

  N.lnWt_1:ydayCos:ydaySin + 
  N.lnWt_2:ydayCos:ydaySin + 

  influx_wk:ydayCos:ydaySin +

  mo(N.catF_1):mo(N.catF_2) +

  wind_L_wk:fetch +
  water_L_wk:fetch +
  windDir_L_wk:fetch +
  waterDir_L_wk:fetch +
  waterDir_R_wk:fetch +
                 
  (1|site.id)
"

# predictors_main <- "
#   temp_L_wk * precip_L_wk * short_wave_L_wk * salinity_L_wk * ydayCos * ydaySin +
#   water_L_wk * waterDir_L_wk * ydayCos * ydaySin * fetch +
#   water_R_wk * waterDir_R_wk * ydayCos * ydaySin * fetch +
#   wind_L_wk * windDir_L_wk * ydayCos * ydaySin * fetch +
#   mo(N.catF_1) * mo(N.catF_2) * ydayCos * ydaySin +
#   N.lnWt_1 * N.lnWt_2 * ydayCos * ydaySin +
#   influx_wk * ydayCos * ydaySin +
#   (1|site.id)"

form_ordinal <- bf(glue("N.catNum | thres(3) ~ {predictors_main}"))
form_bern <- bf(glue("N.bloom ~ {predictors_main}"))





# initial fit -------------------------------------------------------------

out.noBloom <- out.bloom <- out.ord <- vector("list", length(species))

for(sp in 1:length(species)) {
  target <- species[sp]
  target.tf <- thresh.df %>% filter(hab_parameter==target)
  
  target.df <- sampling.df %>%
    filter(grid==1) %>% 
    rename(N=!!target) %>%
    select(obs.id, site.id, date, hour, grid, lon, lat, fetch, bearing, N) %>%
    mutate(yday=yday(date),
           ydayCos=cos(2*pi*yday/365),
           ydaySin=sin(2*pi*yday/365),
           year=year(date),
           bearing=bearing*pi/180,
           N=round(N),
           N.ln=log(N+1),
           N.PA=as.numeric(N>0)) %>%
    rowwise() %>%
    mutate(N.cat=target.tf$tl[max(which(N >= target.tf$min_ge))]) %>%
    mutate(N.catF=factor(N.cat, levels=unique(target.tf$tl), ordered=T),
           N.catNum=as.numeric(N.catF),
           N.bloom=as.numeric(N.catNum > 2)) %>%
    arrange(site.id, date) %>%
    group_by(site.id) %>%
    multijetlag(N.ln, N.PA, N.cat, N.catF, N.bloom, date, n=2) %>%
    ungroup %>%
    mutate(across(starts_with("date_"), ~as.numeric(date-.x)),
           # I don't love this since small if N.ln_x is small OR date_x is large
           N.lnWt_1=N.ln_1/date_1,
           N.lnWt_2=N.ln_2/date_2) %>%
    full_join(hydro.df) %>%
    full_join(connect.df) %>%
    mutate(across(contains("Dir_"), ~cos(.x-bearing))) %>%
    mutate(across(one_of(grep("Dir", covars, invert=T, value=T)), CenterScale)) %>%
    arrange(site.id, date) %>%
    filter(complete.cases(.)) %>%
    select(site.id, lon, lat, date, year, obs.id, fetch, bearing,
           starts_with("N"), starts_with("date_"), starts_with("yday"),
           one_of(covars))
  
  write_csv(target.df, glue("out{sep}dataset_{target}.csv"))
  
  train.df <- target.df %>% filter(year <= 2017)
  
  out.ord[[sp]] <- brm(form_ordinal, data=train.df,
                       family=cumulative("probit"),
                       chains=4, cores=4,
                       iter=2000, warmup=1000, refresh=100, inits="0",
                       control=list(adapt_delta=0.95, max_treedepth=20),
                       prior=prior(horseshoe(3, par_ratio=0.2), class="b"),
                       save_model=glue("out{sep}initFit{sep}ord_{target}.stan"),
                       file=glue("out{sep}initFit{sep}ord_{target}"))
  out.noBloom[[sp]] <- brm(form_bern, data=train.df %>% filter(N.bloom_1==0),
                           family=bernoulli("probit"), 
                           chains=4, cores=4, 
                           iter=2000, warmup=1000, refresh=500, inits="0",
                           control=list(adapt_delta=0.95, max_treedepth=20),
                           prior=prior(horseshoe(3, par_ratio=0.2), class="b"),
                           save_model=glue("out{sep}initFit{sep}bern01_{target}.stan"),
                           file=glue("out{sep}initFit{sep}bern01_{target}"))
  out.bloom[[sp]] <- brm(form_bern, data=train.df %>% filter(N.bloom_1==1),
                         family=bernoulli("probit"), 
                         chains=4, cores=4, 
                         iter=2000, warmup=1000, refresh=500, inits="0",
                         control=list(adapt_delta=0.95, max_treedepth=20),
                         prior=prior(horseshoe(3, par_ratio=0.2), class="b"),
                         save_model=glue("out{sep}initFit{sep}bern11_{target}.stan"),
                         file=glue("out{sep}initFit{sep}bern11_{target}"))
  

  
  cat("Finished", target, "\n")
}



# 
# 
# # NOT RUN -----------------------------------------------------------------
# 
# for(sp in 1:length(species)) {
#   target <- species[sp]
#   condEff.01 <- conditional_effects(out.noBloom[[sp]], surface=T)
#   condEff.11 <- conditional_effects(out.bloom[[sp]], surface=T)
#   plot.01 <- plot(condEff.01, stype="raster", plot=F)
#   plot.11 <- plot(condEff.11, stype="raster", plot=F)
#   p <- map(plot.01[-(57:61)], ~.x +
#              scale_colour_viridis_c("pr(bloom)", option="B") +
#              scale_fill_viridis_c("pr(bloom)", option="B"))
#   p.01 <- ggpubr::ggarrange(plotlist=p, ncol=8, nrow=8, common.legend=T, legend="bottom")
#   ggsave(glue("figs{sep}temp{sep}condEff_01_{target}.png"), p.01, width=15, height=15, dpi=200)
#   p <- map(plot.11[-(57:61)], ~.x +
#              scale_colour_viridis_c("pr(bloom)", option="B") +
#              scale_fill_viridis_c("pr(bloom)", option="B"))
#   p.11 <- ggpubr::ggarrange(plotlist=p, ncol=5, nrow=5, common.legend=T, legend="bottom")
#   ggsave(glue("figs{sep}temp{sep}condEff_11_{target}.png"), p.11, width=15, height=15, dpi=200)
# }

# p <- train.df %>% 
#   filter(N.bloom_1==0) %>% 
#   mutate(pred_mn=colMeans(posterior_epred(out.noBloom[[sp]]))) %>%
#   ggplot(aes(pred_mn, N.bloom)) + 
#   geom_jitter(alpha=0.25, width=0, height=0.01)  + 
#   stat_smooth(method="glm", method.args=list(family="binomial"), fullrange=T) + 
#   ggtitle(glue("{target}: Bloom initiation")) +
#   xlim(0, 1)
# ggsave(glue("figs{sep}temp{sep}pred_y_01_{target}.png"), p, height=6, width=7)
# p <- train.df %>% 
#   filter(N.bloom_1==1) %>% 
#   mutate(pred_mn=colMeans(posterior_epred(out.bloom[[sp]]))) %>%
#   ggplot(aes(pred_mn, N.bloom)) + 
#   geom_jitter(alpha=0.25, width=0, height=0.01)  + 
#   stat_smooth(method="glm", method.args=list(family="binomial"), fullrange=T) + 
#   ggtitle(glue("{target}: Bloom continuation")) +
#   xlim(0, 1)
# ggsave(glue("figs{sep}temp{sep}pred_y_11_{target}.png"), p, height=6, width=7)
# 
# p <- train.df %>% 
#   filter(N.bloom_1==0) %>% 
#   mutate(pred_mn=colMeans(posterior_epred(out.noBloom[[sp]])),
#          pred_round=((pred_mn*100) %/% 5 * 5)/100) %>%
#   group_by(pred_round) %>%
#   summarise(prop_Bloom=mean(N.bloom)) %>%
#   ggplot(aes(pred_round, prop_Bloom)) + 
#   geom_abline(linetype=2) + 
#   geom_point() + stat_smooth(method="lm", se=F) +
#   xlim(0,1) + ylim(0,1) +
#   labs(title=glue("{target}: Bloom initiation"),
#        x="Mean prediction", y="Proportion of observed blooms")
# ggsave(glue("figs{sep}temp{sep}conditProb_01_{target}.png"), p, height=6, width=7)
# p <- train.df %>% 
#   filter(N.bloom_1==1) %>% 
#   mutate(pred_mn=colMeans(posterior_epred(out.bloom[[sp]])),
#          pred_round=((pred_mn*100) %/% 5 * 5)/100) %>%
#   group_by(pred_round) %>%
#   summarise(prop_Bloom=mean(N.bloom)) %>%
#   ggplot(aes(pred_round, prop_Bloom)) + 
#   geom_abline(linetype=2) + 
#   geom_point() + stat_smooth(method="lm", se=F) +
#   xlim(0,1) + ylim(0,1) +
#   labs(title=glue("{target}: Bloom continuation"),
#        x="Mean prediction", y="Proportion of observed blooms")
# ggsave(glue("figs{sep}temp{sep}conditProb_11_{target}.png"), p, height=6, width=7)

