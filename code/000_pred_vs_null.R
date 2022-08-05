## -----------------------------------------------------------------------------
## project
## scriptDescription
## Tim Szewczyk
## -----------------------------------------------------------------------------

library(tidyverse)
library(glue)
library(lubridate)
library(brms)
theme_set(theme_bw())

# testing simple idea of spatially varying effects
species <- str_sub(dir("out", "dataset"), 9, -5)
spPred.df <- vector("list", length(species))

for(i in seq_along(species)) {
  obs.df <- read_csv(glue("out/dataset_{species[i]}.csv")) %>%
    select(site.id, date, N.bloom_1, N.bloom, lon, lat, ydaySin, ydayCos, temp_L_wk, influx_wk) %>%
    mutate(tempDate=ydaySin*ydayCos*temp_L_wk,
           Date=ydaySin*ydayCos,
           influxDate=ydaySin*ydayCos*influx_wk) %>%
    mutate(lon=lon/1e5,
           lat=lat/1e5) # scale for the spline
  fit.0x <- obs.df %>% filter(year(date)<2018) %>% filter(N.bloom_1==0)
  fit.1x <- obs.df %>% filter(year(date)<2018) %>% filter(N.bloom_1==1)
  pred.0x <- obs.df %>% filter(year(date)>2017) %>% filter(N.bloom_1==0)
  pred.1x <- obs.df %>% filter(year(date)>2017) %>% filter(N.bloom_1==1)
  
  temp.f <- bf(N.bloom ~ 
                 bDateS*ydaySin + 
                 bDateC*ydayCos +
                 bDate*Date +
                 bInfluxDate*influx_wk + 
                 bTempDate*temp_L_wk,
               bDateS ~ 1 + (1|site.id),
               bDateC ~ 1 + (1|site.id),
               bDate ~ 1 + (1|site.id),
               bInfluxDate ~ s(ydaySin, ydayCos) + (1|site.id),
               bTempDate ~ s(ydaySin, ydayCos) + (1|site.id),
               nl=T)
  
  temp.priors <- c(prior(normal(0, 1), class="b", nlpar="bDateS"),
                   prior(normal(0, 1), class="b", nlpar="bDateC"),
                   prior(normal(0, 1), class="b", nlpar="bDate"),
                   prior(normal(0, 1), class="b", nlpar="bTempDate"),
                   prior(normal(0, 1), class="b", nlpar="bInfluxDate"))
  
  out.0x <- brm(temp.f, family=bernoulli("probit"), data=fit.0x, 
                iter=2000, init=0, cores=4, chains=4, 
                file=glue("out/spatTest_{species[i]}_0x"))
  out.1x <- brm(temp.f, family=bernoulli("probit"), data=fit.1x, 
                iter=2000, init=0, cores=4, chains=4,
                file=glue("out/spatTest_{species[i]}_1x"))
  
  spPred.df[[i]] <- bind_rows(
    pred.0x %>%
      mutate(sp_mnpr=colMeans(posterior_epred(out.0x, newdata=., allow_new_levels=T))),
    pred.1x %>%
      mutate(sp_mnpr=colMeans(posterior_epred(out.1x, newdata=., allow_new_levels=T)))) %>%
    mutate(species=species[i])
}

spPred.df <- do.call('rbind', spPred.df)










data.df <- dir("out", "dataset_.*csv", full.names=T) %>% 
  map_dfr(~read_csv(.x) %>%
            mutate(month=lubridate::month(date),
                   week=lubridate::week(date)) %>%
            select(site.id, date, obs.id, month, week, date, N.bloom_1) %>%
            mutate(species=str_remove(str_remove(.x, "out/dataset_"), ".csv"))) %>%
  arrange(species, site.id, date) %>%
  group_by(species, month) %>%
  mutate(mo_mnpr=cummean(N.bloom_1)) %>%
  group_by(species, site.id, month) %>%
  mutate(mo_xy_mnpr=cummean(N.bloom_1)) %>%
  group_by(species, site.id) %>%
  mutate(xy_mnpr=cummean(N.bloom_1)) %>%
  ungroup


pred.df <- map_dfr(dir("out", "pred_.*csv"), ~read_csv(glue("out/{.x}")) %>%
                     select(site.id, date, obs.id, N.bloom, N.catF, N.catNum, N.bloom_1,
                            contains("_mnpr")) %>%
                     mutate(species=str_remove(str_remove(.x, "pred_"), ".csv"),
                            bloomThresh=max((!N.bloom)*N.catNum))) %>%
  pivot_longer(starts_with("ord"), names_to="ord_cat", values_to="ord_pr") %>%
  mutate(ord_cat=as.numeric(str_remove(ord_cat, "ord_mnpr"))) %>%
  filter(ord_cat >= bloomThresh) %>%
  group_by(species, obs.id) %>%
  summarise(ord_mnpr=sum(ord_pr),
            across(everything(), ~first(.x))) %>%
  mutate(avg_mnpr=(bern_mnpr + ord_mnpr + rf_mnpr + rf_split_mnpr)/4,
         bayes_mnpr=(bern_mnpr + ord_mnpr)/2,
         ML_mnpr=(rf_mnpr+rf_split_mnpr)/2) %>%
  ungroup %>%
  left_join(., spPred.df %>% select(species, site.id, date, sp_mnpr)) %>%
  mutate(month=lubridate::month(date)) %>%
  select(species, obs.id, site.id, month, date, N.bloom, N.bloom_1, ends_with("_mnpr")) %>%
  left_join(., data.df %>% select(species, obs.id, ends_with("_mnpr")), 
            by=c("species", "obs.id")) %>%
  pivot_longer(ends_with("_mnpr"), names_to="model", values_to="pred") %>%
  mutate(model=str_remove(model, "_mnpr"),
         siteNum=as.numeric(as.factor(site.id))) 


mod_cols <- c(mo_xy="black", mo="grey40", xy="grey70",
              avg="red2", sp="orange",
              bayes="steelblue4", bern="dodgerblue3", ord="dodgerblue", 
              ML="green4", rf="green3", rf_split="green")
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model) %>%
  summarise(Brier=mean((pred-N.bloom)^2, na.rm=T)) %>%
  ggplot(aes(species, Brier, colour=model)) + 
  geom_point(position=position_jitter(width=0.01)) + 
  scale_colour_manual(values=mod_cols) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model) %>%
  summarise(Brier=mean((pred-N.bloom)^2, na.rm=T)) %>%
  group_by(species) %>%
  arrange(species, model) %>%
  mutate(Improve=first(Brier)-Brier) %>%
  ggplot(aes(species, Improve, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_point(position=position_jitter(width=0.05)) + 
  scale_colour_manual(values=mod_cols) +
  labs(x="", y="Improvement in Brier score\nrelative to site-monthly avg") +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model) %>%
  summarise(Brier=mean((pred-N.bloom)^2, na.rm=T)) %>%
  group_by(species) %>%
  arrange(species, model) %>%
  mutate(pctImprove=(first(Brier)-Brier)/first(Brier)*100) %>%
  ggplot(aes(species, pctImprove, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_point(position=position_jitter(width=0.05)) + 
  scale_colour_manual(values=mod_cols) +
  labs(x="", y="Improvement in Brier score\nrelative to site-monthly avg (%)") +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, N.bloom_1, N.bloom) %>%
  summarise(Brier=mean((pred-N.bloom)^2, na.rm=T)) %>%
  mutate(Bloom=if_else(N.bloom==0, "to No Bloom", "to Bloom"),
         Bloom_1=if_else(N.bloom_1==0, "from No Bloom", "from Bloom")) %>%
  ggplot(aes(species, Brier, colour=model)) + 
  geom_point(position=position_jitter(width=0.01)) + 
  scale_colour_manual(values=mod_cols) +
  facet_grid(Bloom_1~Bloom) + 
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, N.bloom_1, N.bloom) %>%
  summarise(Brier=mean((pred-N.bloom)^2, na.rm=T)) %>%
  mutate(Bloom=if_else(N.bloom==0, "to No Bloom", "to Bloom"),
         Bloom_1=if_else(N.bloom_1==0, "from No Bloom", "from Bloom")) %>%
  group_by(species, Bloom, Bloom_1) %>%
  arrange(species, model) %>%
  mutate(Improve=first(Brier)-Brier,
         pctImprove=(first(Brier)-Brier)/first(Brier)*100) %>%
  ggplot(aes(species, Improve, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_point(position=position_jitter(width=0.05)) + 
  scale_colour_manual(values=mod_cols) +
  facet_grid(Bloom_1~Bloom) + 
  labs(x="", y="Improvement in Brier score\nrelative to site-monthly avg") +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, N.bloom_1, N.bloom) %>%
  summarise(Brier=mean((pred-N.bloom)^2, na.rm=T)) %>%
  mutate(Bloom=if_else(N.bloom==0, "to No Bloom", "to Bloom"),
         Bloom_1=if_else(N.bloom_1==0, "from No Bloom", "from Bloom")) %>%
  group_by(species, Bloom, Bloom_1) %>%
  arrange(species, model) %>%
  mutate(Improve=first(Brier)-Brier,
         pctImprove=(first(Brier)-Brier)/first(Brier)*100) %>%
  ggplot(aes(species, pctImprove, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_point(position=position_jitter(width=0.05)) + 
  scale_colour_manual(values=mod_cols) +
  facet_grid(Bloom_1~Bloom, scales="free") + 
  labs(x="", y="Improvement in Brier score\nrelative to site-monthly avg (%)") +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))



# By site
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, siteNum) %>%
  summarise(Brier=mean((pred-N.bloom)^2, na.rm=T)) %>%
  ggplot(aes(siteNum, Brier, colour=model)) + 
  geom_point(position=position_jitter(width=0.01)) + 
  scale_colour_manual(values=mod_cols) +
  facet_wrap(~species, scales="free") 
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, siteNum) %>%
  summarise(Brier=mean((pred-N.bloom)^2, na.rm=T)) %>%
  group_by(species, model) %>%
  summarise(mn=mean(Brier, na.rm=T), 
            lo=quantile(Brier, 0.05, na.rm=T),
            hi=quantile(Brier, 0.95, na.rm=T)) %>%
  ggplot(aes(model, mn, ymin=lo, ymax=hi, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_pointrange() + 
  scale_colour_manual(values=mod_cols) +
  facet_wrap(~species, scales="free")

pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, siteNum) %>%
  summarise(Brier=mean((pred-N.bloom)^2, na.rm=T)) %>%
  group_by(species, siteNum) %>%
  arrange(species, model) %>%
  mutate(Improve=first(Brier)-Brier) %>%
  ggplot(aes(siteNum, Improve, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_point(position=position_jitter(width=0.05)) + 
  scale_colour_manual(values=mod_cols) +
  facet_wrap(~species, scales="free") +
  labs(x="", y="Improvement in Brier score\nrelative to site-monthly avg")
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, siteNum) %>%
  summarise(Brier=mean((pred-N.bloom)^2, na.rm=T)) %>%
  group_by(species, siteNum) %>%
  arrange(species, model) %>%
  mutate(Improve=first(Brier)-Brier) %>%
  group_by(species, model) %>%
  summarise(mn=mean(Improve, na.rm=T), 
            lo=quantile(Improve, 0.05, na.rm=T),
            hi=quantile(Improve, 0.95, na.rm=T)) %>%
  ggplot(aes(model, mn, ymin=lo, ymax=hi, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_pointrange() + 
  scale_colour_manual(values=mod_cols) +
  facet_wrap(~species, scales="free") +
  labs(x="", y="Improvement in Brier score\nrelative to site-monthly avg")

pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, siteNum) %>%
  summarise(Brier=mean((pred-N.bloom)^2, na.rm=T)) %>%
  group_by(species, siteNum) %>%
  arrange(species, model) %>%
  mutate(pctImprove=(first(Brier)-Brier)/first(Brier)*100) %>%
  ggplot(aes(siteNum, pctImprove, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_point(position=position_jitter(width=0.05)) + 
  scale_colour_manual(values=mod_cols) +
  facet_wrap(~species, scales="free")
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, siteNum) %>%
  summarise(Brier=mean((pred-N.bloom)^2, na.rm=T)) %>%
  group_by(species, siteNum) %>%
  arrange(species, model) %>%
  mutate(Brier=Brier+0.01) %>%
  mutate(pctImprove=(first(Brier)-Brier)/first(Brier)*100) %>%
  group_by(species, model) %>%
  summarise(mn=mean(pctImprove, na.rm=T), 
            lo=quantile(pctImprove, 0.05, na.rm=T),
            hi=quantile(pctImprove, 0.95, na.rm=T)) %>%
  ggplot(aes(model, mn, ymin=lo, ymax=hi, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_pointrange() + 
  scale_colour_manual(values=mod_cols) +
  facet_wrap(~species, scales="free") +
  labs(x="", y="Improvement in Brier score\nrelative to site-monthly avg (%)")





# By month
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, month) %>%
  summarise(Brier=mean((pred-N.bloom)^2, na.rm=T)) %>%
  ggplot(aes(month, Brier, colour=model)) + 
  geom_line(position=position_jitter(width=0.01)) + 
  scale_colour_manual(values=mod_cols) +
  scale_x_continuous(breaks=seq(1,12,by=2), labels=month.abb[seq(1,12,by=2)]) +
  facet_wrap(~species, scales="free") 
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, month) %>%
  summarise(Brier=mean((pred-N.bloom)^2, na.rm=T)) %>%
  group_by(species, month) %>%
  arrange(species, model) %>%
  mutate(Improve=first(Brier)-Brier) %>%
  ggplot(aes(month, Improve, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_line(position=position_jitter(width=0.05)) + 
  scale_colour_manual(values=mod_cols) +
  scale_x_continuous(breaks=seq(1,12,by=2), labels=month.abb[seq(1,12,by=2)]) +
  facet_wrap(~species, scales="free") +
  labs(x="", y="Improvement in Brier score\nrelative to site-monthly avg")
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, month) %>%
  summarise(Brier=mean((pred-N.bloom)^2, na.rm=T)) %>%
  group_by(species, month) %>%
  arrange(species, model) %>%
  mutate(Brier=Brier+0.01) %>%
  mutate(pctImprove=(first(Brier)-Brier)/first(Brier)*100) %>%
  ggplot(aes(month, pctImprove, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_line(position=position_jitter(width=0.05)) + 
  scale_colour_manual(values=mod_cols) +
  scale_x_continuous(breaks=seq(1,12,by=2), labels=month.abb[seq(1,12,by=2)]) +
  ylim(-50,NA) +
  facet_wrap(~species, scales="free") +
  labs(x="", y="Improvement in Brier score\nrelative to site-monthly avg (%)")
