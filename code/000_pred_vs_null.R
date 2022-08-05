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
species <- str_sub(dir("out/test_full", "dataset"), 9, -5)
spPred.df <- vector("list", length(species))

for(i in seq_along(species)) {
  obs.df <- read_csv(glue("out/test_full/dataset_{species[i]}.csv")) %>%
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










data.df <- dir("out/test_full", "dataset_.*csv", full.names=T) %>% 
  map_dfr(~read_csv(.x) %>%
            mutate(month=lubridate::month(date),
                   week=lubridate::week(date)) %>%
            select(site.id, date, obs.id, month, week, date, N.bloom_1) %>%
            mutate(species=str_remove(str_remove(.x, "out/test_full/dataset_"), ".csv"))) %>%
  arrange(species, site.id, date) %>%
  group_by(species, month) %>%
  mutate(mo_mnpr=cummean(N.bloom_1)) %>%
  group_by(species, site.id, month) %>%
  mutate(mo_xy_mnpr=cummean(N.bloom_1)) %>%
  group_by(species, site.id) %>%
  mutate(xy_mnpr=cummean(N.bloom_1)) %>%
  ungroup


pred.df <- map_dfr(dir("out/test_full", "pred_.*csv"), ~read_csv(glue("out/test_full/{.x}")) %>%
                     select(site.id, date, obs.id, N.bloom, N.catF, N.catNum, N.bloom_1,
                            contains("_mnpr")) %>%
                     mutate(species=str_remove(str_remove(.x, "pred_int_"), ".csv"),
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


mod_cols <- c(mo="black", xy="grey40", mo_xy="grey70",
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

# McFaddens R2
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  mutate(pred=pmax(pmin(pred, 1-1e-5), 1e-5)) %>%
  group_by(species, model) %>%
  summarise(LL=sum(dbinom(N.bloom, 1, pred, log=T))) %>%
  group_by(species) %>%
  arrange(species, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  ggplot(aes(species, R2, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_point(position=position_jitter(width=0.05)) + 
  scale_colour_manual(values=mod_cols) +
  labs(x="", y="McFadden's R2\nrelative to site-monthly avg") +
  ylim(-0.1, NA) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

# By Bloom_1
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

# McFaddens R2
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  mutate(pred=pmax(pmin(pred, 1-1e-5), 1e-5)) %>%
  group_by(species, model, N.bloom_1, N.bloom) %>%
  summarise(LL=sum(dbinom(N.bloom, 1, pred, log=T))) %>%
  mutate(Bloom=if_else(N.bloom==0, "to No Bloom", "to Bloom"),
         Bloom_1=if_else(N.bloom_1==0, "from No Bloom", "from Bloom")) %>%
  group_by(species, Bloom, Bloom_1) %>%
  arrange(species, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  ggplot(aes(species, R2, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_point(position=position_jitter(width=0.05)) + 
  scale_colour_manual(values=mod_cols) +
  facet_grid(Bloom_1~Bloom, scales="free") + 
  labs(x="", y="McFadden's R2 relative to site-monthly avg") +
  ylim(-0.1, NA) +
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








sams.cols <- c("alexandrium_sp"="#2A5883",
               "dinophysis_sp"="#46BFDE",
               "karenia_mikimotoi"="#A9DAE0",
               "prorocentrum_lima"="#E77334",
               "pseudo_nitzschia_sp"="#D21E4C")

pred.df %>%
  ggplot(aes(pred, N.bloom, colour=model)) + xlim(0,1) +
  geom_jitter(alpha=0.25, width=0, height=0.05) +
  stat_smooth(method="glm", method.args=list(family="binomial"), fullrange=T, se=F) +
  facet_wrap(~species) +
  labs(x="Predicted pr(bloom)", y="Bloom") +
  scale_colour_manual("Model", values=mod_cols) +
  theme_classic() +
  theme(legend.position=c(0.85, 0.2))
ggsave(glue("figs{sep}performance{sep}pr_v_obs_bloom.png"),
       width=8, height=5, dpi=300)

pred.df %>%
  ggplot(aes(date, pred, colour=as.factor(N.bloom), shape=as.factor(N.bloom))) +
  geom_hline(yintercept=0, colour="grey", size=0.5) +
  geom_point(alpha=0.5, size=0.9) +
  scale_colour_manual("Bloom", values=c("grey50", "red3"),
                      guide=guide_legend(override.aes=list(size=1, alpha=1))) +
  scale_shape_manual("Bloom", values=c(1, 19)) +
  facet_grid(species~model) +
  theme_bw() +
  theme(axis.title.x=element_blank(), panel.grid=element_blank()) +
  scale_y_continuous(breaks=c(0, 0.5, 1), limits=c(0, 1)) +
  labs(y="Pr(bloom)", title="All observations")
ggsave(glue("figs{sep}performance{sep}pr_v_obs_bloom_dates.png"),
       width=10, height=7.25, dpi=300)

pred.df %>%
  filter(N.bloom_1==0) %>%
  ggplot(aes(date, pred, colour=as.factor(N.bloom), shape=as.factor(N.bloom))) +
  geom_hline(yintercept=0, colour="grey", size=0.5) +
  geom_point(alpha=0.5, size=0.9) +
  scale_colour_manual("Bloom", values=c("grey50", "red3"),
                      guide=guide_legend(override.aes=list(size=1, alpha=1))) +
  scale_shape_manual("Bloom", values=c(1, 19)) +
  facet_grid(species~model) +
  theme_bw() +
  theme(axis.title.x=element_blank(), panel.grid=element_blank()) +
  scale_y_continuous(breaks=c(0, 0.5, 1), limits=c(0, 1)) +
  labs(y="Pr(bloom)", title="Bloom initiation")
ggsave(glue("figs{sep}performance{sep}pr_v_obs_bloom_dates_init.png"),
       width=10, height=7.25, dpi=300)

pred.df %>%
  group_by(species, date, N.bloom, model) %>%
  summarise(pr_mn=mean(pred),
            pr_lo=quantile(pred, 0.05),
            pr_hi=quantile(pred, 0.95)) %>%
  ggplot(aes(date, pr_mn, colour=as.factor(N.bloom), shape=as.factor(N.bloom))) +
  geom_hline(yintercept=0, colour="grey", size=0.5) +
  geom_point(alpha=0.5, size=0.9) +
  scale_colour_manual("Bloom", values=c("grey50", "red3"),
                      guide=guide_legend(override.aes=list(size=1, alpha=1))) +
  scale_shape_manual("Bloom", values=c(1, 19)) +
  facet_grid(species~model) +
  theme_bw() +
  theme(axis.title.x=element_blank(), panel.grid=element_blank()) +
  scale_y_continuous(breaks=c(0, 0.5, 1), limits=c(0, 1)) +
  labs(y="Pr(bloom)", title="All observations")
ggsave(glue("figs{sep}performance{sep}pr_v_obs_bloom_dates_dayMeans.png"),
       width=10, height=7.25, dpi=300)


pred.df %>%
  filter(N.bloom_1==0) %>%
  mutate(pred_round=((pred*100) %/% 5 * 5)/100) %>%
  group_by(species, model, pred_round) %>%
  summarise(obs_bloom=mean(N.bloom), N=n()) %>%
  ggplot(aes(pred_round, obs_bloom, colour=model)) +
  geom_abline(linetype=2) +
  geom_point(alpha=0.5, aes(size=N)) +
  stat_smooth(method="lm", se=F) +
  xlim(0,1) + ylim(0,1) +
  scale_colour_manual("Model", values=mod_cols) +
  scale_size_continuous(breaks=c(1, 100, 500)) +
  facet_wrap(~species) +
  theme_classic() +
  theme(legend.position="bottom") +
  labs(title="Bloom initiation", x="Mean prediction", y="Proportion observed blooms")
ggsave(glue("figs{sep}performance{sep}pred_v_obs_bloom_init.png"),
       width=8, height=6, dpi=300)

pred.df %>%
  filter(N.bloom_1==0) %>%
  group_by(species, date, N.bloom, model) %>%
  summarise(pr_mn=mean(pred),
            pr_lo=quantile(pred, 0.05),
            pr_hi=quantile(pred, 0.95)) %>%
  ggplot(aes(date, pr_mn, colour=as.factor(N.bloom), shape=as.factor(N.bloom))) +
  geom_hline(yintercept=0, colour="grey", size=0.5) +
  geom_point(alpha=0.5, size=0.9) +
  scale_colour_manual("Bloom", values=c("grey50", "red3"),
                      guide=guide_legend(override.aes=list(size=1, alpha=1))) +
  scale_shape_manual("Bloom", values=c(1, 19)) +
  facet_grid(species~model) +
  theme_bw() +
  theme(axis.title.x=element_blank(), panel.grid=element_blank()) +
  scale_y_continuous(breaks=c(0, 0.5, 1), limits=c(0, 1)) +
  labs(y="Pr(bloom)", title="Bloom initiation")
ggsave(glue("figs{sep}performance{sep}pr_v_obs_bloom_dates_dayMeans.png"),
       width=10, height=7.25, dpi=300)

pred.df %>%
  mutate(pred_round=((pred*100) %/% 5 * 5)/100) %>%
  group_by(species, model, pred_round) %>%
  summarise(obs_bloom=mean(N.bloom), N=n()) %>%
  ggplot(aes(pred_round, obs_bloom, colour=model)) +
  geom_abline(linetype=2) +
  geom_point(alpha=0.5, aes(size=N)) +
  stat_smooth(method="lm", se=F) +
  xlim(0,1) + ylim(0,1) +
  scale_colour_manual("Model", values=mod_cols) +
  scale_size_continuous(breaks=c(1, 100, 500)) +
  facet_wrap(~species) +
  theme_classic() +
  theme(legend.position="bottom") +
  labs(title="All predictions", x="Mean prediction", y="Proportion observed blooms")
ggsave(glue("figs{sep}performance{sep}pred_v_obs.png"),
       width=8, height=6, dpi=300)

pred.df %>%
  ggplot(aes(model, pred, fill=factor(N.bloom))) +
  geom_boxplot() + facet_grid(N.bloom_1~species)



binary.ls <- pred.df %>% #filter(N.bloom_1==0) %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  select(obs.id, species, N.bloom, model, pred) %>%
  group_by(model) %>%
  group_split() %>%
  map(~.x %>% group_by(species) %>% group_split)
names(binary.ls) <- names(mod_cols)

roc.ls <- map_depth(binary.ls, 2,
                    ~pROC::roc(.x$N.bloom ~ .x$pred))


for(i in 1:5) {
  par(mfrow=c(1,1))
  png(glue("figs/performance/ROC_{species[i]}.png"), width=5, height=5, res=300, units="in")
  plot(NA, NA, xlim=c(0, 1), ylim=c(0, 1),
       xlab="1 - Specificity", ylab="Sensitivity", axes=F, main=species[i])
  axis(side=1, at=c(0, 0.5, 1))
  axis(side=2, at=c(0, 0.5, 1))
  abline(a=0, b=1, col="grey30", lwd=0.5)
  map(seq_along(roc.ls),
      ~lines(1-roc.ls[[.x]][[i]]$specificities, roc.ls[[.x]][[i]]$sensitivities,
             col=mod_cols[.x], lwd=2.5))
  legend("bottomright",
         paste(names(mod_cols), map_dbl(roc.ls, ~round(.x[[i]]$auc, 2)),  sep=": "),
         col=mod_cols, lty=1, lwd=2, bty="n", cex=0.95)
  dev.off()
}

auc.df <- map_depth(roc.ls, 2, ~.x$auc) %>% 
  map_dfc(~do.call(c, .x)) %>%
  mutate(species=species) %>%
  pivot_longer(-species, names_to="model", values_to="AUC")
write_csv(auc.df, glue("figs/performance/auc_all.csv"))

auc.df <- tibble(
  bern=map(roc.ls[["bern"]], ~.x$auc) %>% do.call(c, .) %>% round(3),
  ord=map(roc.ls[["ord"]], ~.x$auc) %>% do.call(c, .) %>% round(3),
  rf=map(roc.ls[["rf"]], ~.x$auc) %>% do.call(c, .) %>% round(3),
  rf_s=map(roc.ls[["rf_s"]], ~.x$auc) %>% do.call(c, .) %>% round(3),
  avg=map(roc.ls[["avg"]], ~.x$auc) %>% do.call(c, .) %>% round(3)
) %>%
  mutate(species=species) %>%
  pivot_longer(1:5, names_to="model", values_to="AUC")
write_csv(auc.df, glue("figs{sep}performance{sep}auc_all.csv"))

par(mfrow=c(1,1))





auc.f <- dir(glue("figs{sep}performance{sep}"), "auc_.*csv", full.names=T)
auc.df <- map_dfr(auc.f, ~read_csv(.x, show_col_types=F) %>%
                    mutate(type=str_split_fixed(.x, "auc_", 2)[,2])) %>%
  mutate(data_subset=if_else(grepl("01", type), "01", "all"),
         covars=if_else(grepl("dayOnly", type), "dayOnly", "full"))
ggplot(auc.df, aes(model, AUC, colour=covars)) +
  geom_point() +
  facet_grid(species~data_subset)

auc.df %>% select(-type) %>%
  pivot_wider(names_from="covars", values_from="AUC") %>%
  mutate(AUC_diff=full-dayOnly) %>%
  ggplot(aes(model, AUC_diff, colour=data_subset)) +
  geom_hline(yintercept=0) +
  geom_point() +
  facet_wrap(~species) + theme_classic()






dateConditions <- tibble(yday=seq(1, 365, length.out=12)) %>%
  mutate(ydayCos=cos(yday*2*pi/365),
         ydaySin=sin(yday*2*pi/365))

ord.f <- dir("out/test_full", "ord_", full.names=T)
bern01.f <- dir("out/test_full", "bern01_", full.names=T)
ord.eff <- map(ord.f, 
               ~fixef(readRDS(.x), probs=c(0.1, 0.9)) %>% 
                 as_tibble(rownames="var") %>% 
                 filter(sign(Q10)==sign(Q90)))
ce.f <- map2(ord.f, ord.eff, 
             ~readRDS(.x) %>%
               conditional_effects(.,
                                   effects=grep("Intercept", .y$var, invert=T, value=T),
                                   conditions=dateConditions))

