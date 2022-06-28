# HABReports Bayesian modelling
# Bayesian HAB model initial fits
# Tim Szewczyk


# This script fits models with data from 2013-2017 using horseshoe priors



# set up ------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "LaplacesDemon", "brms", 
          "bayesplot", "WeStCOMS", "randomForest", "caret")
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
  ungroup %>%
  mutate(bloom=as.numeric(!is.na(alert)))

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
connect.df <- read_csv(glue("data{sep}connectivity_5e3parts.csv")) %>%
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

predictors_main <- c(
  "ydayCos", "ydaySin", 
  
  "short_wave_L_wk:ydayCos:ydaySin",
  "temp_L_wk:ydayCos:ydaySin",
  "precip_L_wk:ydayCos:ydaySin",
  "salinity_L_wk:ydayCos:ydaySin",
  "water_L_wk:ydayCos:ydaySin",
  "water_R_wk:ydayCos:ydaySin",
  "wind_L_wk:ydayCos:ydaySin",
  "windDir_L_wk:ydayCos:ydaySin",
  "waterDir_L_wk:ydayCos:ydaySin",
  "waterDir_R_wk:ydayCos:ydaySin",
  
  "fetch", 
  "influx_wk",
  
  "mo(N.catF_1):ydayCos:ydaySin",
  "mo(N.catF_2):ydayCos:ydaySin",
  
  "N.lnWt_1:ydayCos:ydaySin",
  "N.lnWt_2:ydayCos:ydaySin", 
  
  "influx_wk:ydayCos:ydaySin",
  
  "mo(N.catF_1):mo(N.catF_2)",
  
  "wind_L_wk:fetch",
  "water_L_wk:fetch",
  "water_R_wk:fetch",
  "windDir_L_wk:fetch",
  "waterDir_L_wk:fetch",
  "waterDir_R_wk:fetch",
  
  "(1|site.id)"
)


# predictors_main <- "
#   temp_L_wk * precip_L_wk * short_wave_L_wk * salinity_L_wk * ydayCos * ydaySin +
#   water_L_wk * waterDir_L_wk * ydayCos * ydaySin * fetch +
#   water_R_wk * waterDir_R_wk * ydayCos * ydaySin * fetch +
#   wind_L_wk * windDir_L_wk * ydayCos * ydaySin * fetch +
#   mo(N.catF_1) * mo(N.catF_2) * ydayCos * ydaySin +
#   N.lnWt_1 * N.lnWt_2 * ydayCos * ydaySin +
#   influx_wk * ydayCos * ydaySin +
#   (1|site.id)"

form_ordinal <- bf(glue("N.catNum | thres(3) ~ {paste(predictors_main, collapse=' + ')}"))
form_bern <- bf(glue("N.bloom ~ {paste(predictors_main, collapse=' + ')}"))
form_bern_noCatF <- bf(glue("N.bloom ~ {paste(grep('catF_1', predictors_main, value=T, invert=T), collapse=' + ')}"))





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
    ungroup %>%
    mutate(N.catF=factor(N.cat, levels=unique(target.tf$tl), ordered=T),
           N.catNum=as.numeric(N.catF),
           N.bloom=target.tf$bloom[match(N.cat, target.tf$tl)]) %>%
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
  test.df <- target.df %>% filter(year > 2017)
  
  out.ord[[sp]] <- brm(form_ordinal, data=train.df,
                       family=cumulative("probit"),
                       chains=4, cores=4,
                       iter=2000, warmup=1000, refresh=100, inits="0",
                       control=list(adapt_delta=0.95, max_treedepth=20),
                       prior=prior(horseshoe(3, par_ratio=0.2), class="b"),
                       save_model=glue("out{sep}initFit{sep}ord_{target}.stan"),
                       file=glue("out{sep}initFit{sep}ord_{target}"))
  if(n_distinct(filter(train.df, N.bloom_1==0)$N.catF_1)==1) {
    form_01 <- form_bern_noCatF
  } else {
    form_01 <- form_bern
  }
  out.noBloom[[sp]] <- brm(form_01, data=train.df %>% filter(N.bloom_1==0),
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
  train.rf <- train.df %>%
    select(N.bloom, 
           site.id, lon, lat, 
           ydaySin, ydayCos, fetch,
           N.ln_1, N.cat_1, 
           N.ln_2, N.cat_2, 
           N.lnWt_1, N.lnWt_2, 
           one_of(covars)) %>%
    mutate(N.bloom=factor(N.bloom)) %>%
    as.data.frame()
  test.rf <- test.df %>%
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
                     proximity=T, keep.forest=T)
  p.rf <- predict(rf, test.rf, type="prob")
  
  pred.df <- test.df
  pred.ord <- posterior_epred(out.ord[[sp]], 
                              newdata=test.df, 
                              allow_new_levels=T)
  cat.iter <- apply(pred.ord, 1:2, which.max)
  pred.df <- pred.df %>%
    mutate(ord_mnpr0=c(colMeans(pred.ord[,,1,drop=F])),
           ord_mnpr1=c(colMeans(pred.ord[,,2,drop=F])),
           ord_mnpr2=c(colMeans(pred.ord[,,3,drop=F])),
           ord_mnpr3=c(colMeans(pred.ord[,,4,drop=F])),
           ord_prmax0=colMeans(cat.iter==1),
           ord_prmax1=colMeans(cat.iter==2),
           ord_prmax2=colMeans(cat.iter==3),
           ord_prmax3=colMeans(cat.iter==4),
           ord_mnCat=colMeans(cat.iter),
           ord_wtmnCat=ord_mnpr0 + 2*ord_mnpr1 + 3*ord_mnpr2 + 4*ord_mnpr3)
  pred.bern01 <- posterior_epred(out.noBloom[[sp]], 
                                   newdata=test.df %>% filter(N.bloom_1==0) %>% droplevels, 
                                   allow_new_levels=T) 
  pred.01 <- test.df %>% filter(N.bloom_1==0) %>% select(obs.id) %>%
    mutate(bern_mnpr=colMeans(pred.bern01))
  pred.bern11 <- posterior_epred(out.bloom[[sp]], 
                                 newdata=test.df %>% filter(N.bloom_1==1) %>% droplevels, 
                                 allow_new_levels=T) 
  pred.11 <- test.df %>% filter(N.bloom_1==1) %>% select(obs.id) %>%
    mutate(bern_mnpr=colMeans(pred.bern11))
  pred.bern <- bind_rows(pred.01, pred.11)
  pred.df <- full_join(pred.df, pred.bern) %>%
    mutate(rf_mnpr=p.rf[,2])
  write_csv(pred.df, glue("out{sep}initFit{sep}pred_{target}.csv"))
  
  cat("Finished", target, "\n")
}



# 
# 
# # NOT RUN -----------------------------------------------------------------
# 

# tl.col <- c("green3", "gold1", "orange", "red")
# sams.cols <- c("alexandrium_sp"="#2A5883",
#                "dinophysis_sp"="#46BFDE", 
#                "karenia_mikimotoi"="#A9DAE0",
#                "prorocentrum_lima"="#E77334", 
#                "pseudo_nitzschia_sp"="#D21E4C")
# 
# pred.f <- dir(glue("out{sep}initFit"), "pred_")
# pred.df <- map_dfr(pred.f, ~read_csv(glue("out{sep}initFit{sep}{.x}")) %>%
#                      select(obs.id, N.bloom, N.catF, N.catNum, N.bloom_1,
#                             contains("_mnpr")) %>%
#                      mutate(species=str_remove(str_remove(.x, "pred_"), ".csv"),
#                             bloomThresh=max((!N.bloom)*N.catNum))) %>% 
#   pivot_longer(starts_with("ord"), names_to="ord_cat", values_to="ord_pr") %>%
#   mutate(ord_cat=as.numeric(str_remove(ord_cat, "ord_mnpr"))) %>% 
#   filter(ord_cat >= bloomThresh) %>% 
#   group_by(species, obs.id) %>% 
#   summarise(ord_pr=sum(ord_pr), 
#             across(everything(), ~first(.x))) %>%
#   mutate(avg_pr=(bern_mnpr + ord_pr + rf_mnpr)/3)
# 
# ggplot(pred.df, aes(bern_mnpr, N.bloom)) + xlim(0, 1) +
#   geom_jitter(alpha=0.25, width=0, height=0.01) +
#   stat_smooth(method="glm", method.args=list(family="binomial"), fullrange=T, se=F) +
#   facet_wrap(~species)
# ggplot(pred.df, aes(ord_pr, N.bloom)) + xlim(0, 1) +
#   geom_jitter(alpha=0.25, width=0, height=0.01) +
#   stat_smooth(method="glm", method.args=list(family="binomial"), fullrange=T, se=F) +
#   facet_wrap(~species)
# ggplot(pred.df, aes(rf_mnpr, N.bloom)) + xlim(0, 1) +
#   geom_jitter(alpha=0.25, width=0, height=0.01) +
#   stat_smooth(method="glm", method.args=list(family="binomial"), fullrange=T, se=F) +
#   facet_wrap(~species)
# 
# pred.df %>%
#   filter(N.bloom_1==0) %>%
#   mutate(pred_round=((bern_mnpr*100) %/% 5 * 5)/100) %>%
#   group_by(species, pred_round) %>%
#   summarise(prop_Bloom=mean(N.bloom)) %>%
#   ggplot(aes(pred_round, prop_Bloom)) +
#   geom_abline(linetype=2) +
#   geom_point() + stat_smooth(method="lm", se=F) +
#   xlim(0,1) + ylim(0,1) +
#   facet_wrap(~species) +
#   labs(title="Bern: Bloom initiation",
#        x="Mean prediction", y="Proportion of observed blooms")
# 
# pred.df %>%
#   filter(N.bloom_1==0) %>%
#   mutate(pred_round=(((ord_pr)*100) %/% 5 * 5)/100) %>%
#   group_by(species, pred_round) %>%
#   summarise(prop_Bloom=mean(N.bloom)) %>%
#   ggplot(aes(pred_round, prop_Bloom)) +
#   geom_abline(linetype=2) +
#   geom_point() + stat_smooth(method="lm", se=F) +
#   xlim(0,1) + ylim(0,1) +
#   facet_wrap(~species) +
#   labs(title="Ord: Bloom initiation",
#        x="Mean prediction", y="Proportion of observed blooms")
# 
# pred.df %>%
#   filter(N.bloom_1==0) %>%
#   mutate(pred_round=(((rf_mnpr)*100) %/% 5 * 5)/100) %>%
#   group_by(species, pred_round) %>%
#   summarise(prop_Bloom=mean(N.bloom)) %>%
#   ggplot(aes(pred_round, prop_Bloom)) +
#   geom_abline(linetype=2) +
#   geom_point() + stat_smooth(method="lm", se=F) +
#   xlim(0,1) + ylim(0,1) +
#   facet_wrap(~species) +
#   labs(title="RF: Bloom initiation",
#        x="Mean prediction", y="Proportion of observed blooms")
# 
# ggplot(pred.df, aes(as.factor(N.bloom_1), bern_mnpr, fill=N.catF)) + 
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   scale_fill_manual(values=tl.col) +
#   ggtitle("Dual bernoulli")
# 
# ggplot(pred.df, aes(as.factor(N.bloom_1), ord_pr, fill=N.catF)) + 
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   scale_fill_manual(values=tl.col) +
#   ggtitle("Ordinal")
# 
# ggplot(pred.df, aes(as.factor(N.bloom_1), rf_mnpr, fill=N.catF)) + 
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   scale_fill_manual(values=tl.col) +
#   ggtitle("RF")
# 
# ggplot(pred.df, aes(as.factor(N.bloom_1), avg_pr, fill=N.catF)) + 
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   scale_fill_manual(values=tl.col) +
#   ggtitle("Avg")
# 
# 
# ggplot(pred.df, aes(as.factor(N.bloom_1), bern_mnpr, fill=factor(N.bloom))) + 
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   ggtitle("Dual bernoulli")
# 
# ggplot(pred.df, aes(as.factor(N.bloom_1), ord_pr, fill=factor(N.bloom))) + 
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   ggtitle("Ordinal")
# 
# ggplot(pred.df, aes(as.factor(N.bloom_1), rf_mnpr, fill=factor(N.bloom))) + 
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   ggtitle("RF")
# 
# ggplot(pred.df, aes(as.factor(N.bloom_1), avg_pr, fill=factor(N.bloom))) + 
#   geom_boxplot() + facet_wrap(~species) + ylim(0, 1) +
#   ggtitle("Avg")
# 
# 
# binary.df <- pred.df %>% #filter(N.bloom_1==0) %>%
#   select(obs.id, species, ord_pr, bern_mnpr, rf_mnpr, avg_pr, N.bloom) %>%
#   group_by(species) %>%
#   group_split()
# 
# roc.ord <- map(binary.df, ~pROC::roc(.x$N.bloom ~ .x$ord_pr))
# roc.bern <- map(binary.df, ~pROC::roc(.x$N.bloom ~ .x$bern_mnpr))
# roc.rf <- map(binary.df, ~pROC::roc(.x$N.bloom ~ .x$rf_mnpr))
# roc.avg <- map(binary.df, ~pROC::roc(.x$N.bloom ~ .x$avg_pr))
# 
# par(mfrow=c(1,4))
# plot(NA, NA, xlim=c(0, 1), ylim=c(0, 1), xlab="", ylab="", axes=F, main="Ordinal")
# axis(side=1, at=c(0, 0.5, 1))
# axis(side=2, at=c(0, 0.5, 1))
# abline(a=0, b=1, col="grey30", lwd=0.5)
# map(1:5, ~lines(1-roc.ord[[.x]]$specificities, roc.ord[[.x]]$sensitivities, 
#                 col=sams.cols[binary.df[[.x]]$species[1]], lwd=2))
# 
# plot(NA, NA, xlim=c(0, 1), ylim=c(0, 1), xlab="", ylab="", axes=F, main="Bernoulli")
# axis(side=1, at=c(0, 0.5, 1))
# axis(side=2, at=c(0, 0.5, 1))
# abline(a=0, b=1, col="grey30", lwd=0.5)
# map(1:5, ~lines(1-roc.bern[[.x]]$specificities, roc.bern[[.x]]$sensitivities, 
#                 col=sams.cols[binary.df[[.x]]$species[1]], lwd=2))
# 
# plot(NA, NA, xlim=c(0, 1), ylim=c(0, 1), xlab="", ylab="", axes=F, main="RF")
# axis(side=1, at=c(0, 0.5, 1))
# axis(side=2, at=c(0, 0.5, 1))
# abline(a=0, b=1, col="grey30", lwd=0.5)
# map(1:5, ~lines(1-roc.rf[[.x]]$specificities, roc.rf[[.x]]$sensitivities, 
#                 col=sams.cols[binary.df[[.x]]$species[1]], lwd=2))
# 
# plot(NA, NA, xlim=c(0, 1), ylim=c(0, 1), xlab="", ylab="", axes=F, main="Avg")
# axis(side=1, at=c(0, 0.5, 1))
# axis(side=2, at=c(0, 0.5, 1))
# abline(a=0, b=1, col="grey30", lwd=0.5)
# map(1:5, ~lines(1-roc.avg[[.x]]$specificities, roc.avg[[.x]]$sensitivities, 
#                 col=sams.cols[binary.df[[.x]]$species[1]], lwd=2))
# 
# legend("bottomright", names(sams.cols), col=sams.cols, lty=1, lwd=2, bty="n")
# 
# tibble(
#   bern=map(roc.bern, ~.x$auc) %>% do.call(c, .) %>% round(3),
#   ord=map(roc.ord, ~.x$auc) %>% do.call(c, .) %>% round(3),
#   rf=map(roc.rf, ~.x$auc) %>% do.call(c, .) %>% round(3),
#   avg=map(roc.avg, ~.x$auc) %>% do.call(c, .) %>% round(3)
# ) %>%
#   mutate(species=species) %>%
#   pivot_longer(1:4, names_to="model", values_to="AUC") %>%
#   ggplot(aes(model, AUC)) + 
#   geom_point(stat="identity", position="dodge") + 
#   facet_wrap(~species, nrow=1) + ylim(0.5, 1)
# 
# 
# par(mfrow=c(1,1))
# 
# 
# 
# 
# out.b01 <- map(dir(glue("out{sep}initFit{sep}"), "bern01.*rds", full.names=T),
#                readRDS)
# p.ls <- map(out.b01, ~plot(conditional_effects(.x, surface=T), stype="raster", plot=F))
# walk2(p.ls, species, ~ggpubr::ggarrange(plotlist=.x, ncol=4, nrow=5) %>% 
#         ggsave(filename=glue("figs{sep}temp{sep}bern01_{.y}.png"), plot=., width=16, height=16))
# 
# 
# out.b11 <- map(dir(glue("out{sep}initFit{sep}"), "bern11.*rds", full.names=T),
#                readRDS)
# p.ls <- map(out.b11, ~plot(conditional_effects(.x, surface=T), stype="raster", plot=F))
# walk2(p.ls, species, ~ggpubr::ggarrange(plotlist=.x, ncol=4, nrow=5) %>% 
#         ggsave(filename=glue("figs{sep}temp{sep}bern11_{.y}.png"), plot=., width=16, height=16))
# 
# 
# out.ord <- map(dir(glue("out{sep}initFit{sep}"), "ord.*rds", full.names=T),
#                readRDS)
# p.ls <- map(out.ord, ~plot(conditional_effects(.x, surface=T), stype="raster", plot=F))
# walk2(p.ls, species, ~ggpubr::ggarrange(plotlist=.x, ncol=4, nrow=5) %>% 
#         ggsave(filename=glue("figs{sep}temp{sep}ord_{.y}.png"), plot=., width=16, height=16))
