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
species <- unique(str_sub(str_split_fixed(dir("out/test_full", "dataset"), "_", 3)[,3], 1, -5))


data.df <- dir("out/test_full", "dataset.*csv", full.names=T) %>% 
  map_dfr(~read_csv(.x) %>%
            mutate(month=lubridate::month(date),
                   week=lubridate::week(date)) %>%
            select(siteid, date, obsid, month, week, date, Nbloom1) %>%
            mutate(species=str_sub(str_split_fixed(.x, "_", 4)[,4], 1, -5),
                   covarSet=str_split_fixed(.x, "_", 4)[,3])) %>%
  arrange(species, covarSet, date, siteid) %>%
  group_by(species, covarSet) %>%
  mutate(mean_mnpr=cummean(Nbloom1)) %>%
  group_by(species, covarSet, month) %>%
  mutate(mo_mnpr=cummean(Nbloom1)) %>%
  group_by(species, covarSet, siteid, month) %>%
  mutate(mo_xy_mnpr=cummean(Nbloom1)) %>%
  group_by(species, covarSet, siteid) %>%
  mutate(xy_mnpr=cummean(Nbloom1)) %>%
  ungroup

# TODO: Calculate LL using cross-validated runs, k=nYears[fit

fit.df <- map_dfr(dir("out/test_full/", "^fit.*csv"), 
              ~read_csv(glue("out/test_full/{.x}")) %>%
                select(covarSet, siteid, date, obsid, Nbloom, Nbloom1, contains("_mnpr")) %>%
                mutate(month=lubridate::month(date), 
                       species=str_sub(str_split_fixed(.x, "_", 3)[,3], 1, -5)))
LL.tot <- fit.df %>% 
  group_by(covarSet, species) %>% 
  summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
  mutate(across(ends_with("wt"), ~.x/rowSums(across(ends_with("wt"))))) %>% 
  ungroup
LL.mo <- fit.df %>% 
  group_by(covarSet, species, month) %>% 
  summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
  mutate(across(ends_with("wt"), ~.x/rowSums(across(ends_with("wt"))))) %>%
  ungroup
LL.st <- fit.df %>% 
  group_by(covarSet, species, siteid) %>% 
  summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
  mutate(across(ends_with("wt"), ~.x/rowSums(across(ends_with("wt"))))) %>%
  ungroup
LL.B1 <- fit.df %>% 
  group_by(covarSet, species, Nbloom1) %>% 
  summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
  mutate(across(ends_with("wt"), ~.x/rowSums(across(ends_with("wt"))))) %>%
  ungroup
LL.moB1 <- fit.df %>% 
  group_by(covarSet, species, month, Nbloom1) %>% 
  summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
  mutate(across(ends_with("wt"), ~.x/rowSums(across(ends_with("wt"))))) %>%
  ungroup
LL.stB1 <- fit.df %>% 
  group_by(covarSet, species, siteid, Nbloom1) %>% 
  summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
  mutate(across(ends_with("wt"), ~.x/rowSums(across(ends_with("wt"))))) %>%
  ungroup
LL.moSt <- fit.df %>% 
  group_by(covarSet, species, month, siteid) %>% 
  summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
  mutate(across(ends_with("wt"), ~.x/rowSums(across(ends_with("wt"))))) %>%
  ungroup
LL.moStB1 <- fit.df %>% 
  group_by(covarSet, species, month, siteid, Nbloom1) %>% 
  summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
  mutate(across(ends_with("wt"), ~.x/rowSums(across(ends_with("wt"))))) %>%
  ungroup

pred.df <- map_dfr(dir("out/test_full/", "^pred.*csv"), 
                  ~read_csv(glue("out/test_full/{.x}")) %>%
                    select(covarSet, siteid, date, obsid, Nbloom, Nbloom1, contains("_mnpr")) %>%
                    mutate(month=lubridate::month(date), 
                           species=str_sub(str_split_fixed(.x, "_", 3)[,3], 1, -5)))

pred.ens <- pred.df %>% 
  left_join(LL.tot) %>% rowwise() %>%
  mutate(avg=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  ungroup() %>% select(-ends_with("wt")) %>% 
  left_join(LL.mo) %>% rowwise() %>%
  mutate(avgMo=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  ungroup() %>% select(-ends_with("wt")) %>% 
  left_join(LL.st) %>% rowwise() %>%
  mutate(avgSt=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  ungroup() %>% select(-ends_with("wt")) %>% 
  left_join(LL.B1) %>% rowwise() %>%
  mutate(avgB1=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  ungroup() %>% select(-ends_with("wt")) %>% 
  left_join(LL.moB1) %>% rowwise() %>%
  mutate(avgMoB1=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  ungroup() %>% select(-ends_with("wt")) %>% 
  left_join(LL.stB1) %>% rowwise() %>%
  mutate(avgStB1=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  ungroup() %>% select(-ends_with("wt")) %>% 
  left_join(LL.moSt) %>% rowwise() %>%
  mutate(avgMoSt=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  ungroup() %>% select(-ends_with("wt")) %>%
  left_join(LL.moStB1) %>% rowwise() %>%
  mutate(avgMoStB1=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  ungroup() %>% select(-ends_with("wt"))
pred.df <- pred.ens %>% 
  mutate(avgSt=if_else(is.na(avgSt), avg, avgSt),
         avgMo=if_else(is.na(avgMo), avg, avgMo),
         avgStB1=if_else(is.na(avgStB1), avgB1, avgStB1),
         avgMoSt=if_else(is.na(avgMoSt), avgMo, avgMoSt),
         avgMoB1=if_else(is.na(avgMoB1), avgB1, avgMoB1),
         avgMoStB1=if_else(is.na(avgMoStB1), avgMoB1, avgMoStB1)) %>%
  rename(avg_mnpr=avg, 
         avgMo_mnpr=avgMo,
         avgSt_mnpr=avgSt,
         avgB1_mnpr=avgB1,
         avgMoB1_mnpr=avgMoB1,
         avgStB1_mnpr=avgStB1,
         avgMoSt_mnpr=avgMoSt,
         avgMoStB1_mnpr=avgMoStB1) %>%
  group_by(species, covarSet) %>%
  left_join(., data.df %>% select(species, covarSet, obsid, ends_with("_mnpr")),
            by=c("species", "covarSet", "obsid")) %>%
  pivot_longer(ends_with("_mnpr"), names_to="model", values_to="pred") %>%
  mutate(model=str_remove(model, "_mnpr"),
         siteNum=as.numeric(as.factor(siteid))) %>%
  mutate(modType=case_when(grepl("mo|xy|mean", model) ~ "Null",
                           grepl("LP", model) ~ "Laplace",
                           grepl("rf", model) ~ "RF",
                           model %in% c("bern", "ord") ~ "intrct",
                           model %in% c("bernP", "ordP") ~ "indVar",
                           grepl("avg", model) ~ "ens"))
  

# LDA
lda.train.df <- fit.df %>% ungroup %>%
  filter(covarSet=="date") %>%
  filter(species=="alexandrium_sp") %>%
  dplyr::select(obsid, Nbloom, ends_with("mnpr")) %>% 
  rename_with(~str_sub(.x, 1, -6), ends_with("mnpr"))
lda.test.df <- pred.df %>% ungroup %>%
  filter(covarSet=="date") %>%
  filter(!grepl("^avg", model)) %>%
  filter(! model %in% c("mean", "mo", "mo_xy", "xy")) %>%
  filter(species=="alexandrium_sp") %>%
  dplyr::select(obsid, Nbloom, model, pred) %>% 
  pivot_wider(names_from="model", values_from="pred")
library(MASS)
fit.lda <- lda(formula=Nbloom ~ ord + ordP + ordLP + rf + rf_split + bern + bernP + bernLP,
               data=lda.train.df)
fitCV.lda <- lda(formula=Nbloom ~ ord + ordP + ordLP + rf + rf_split + bern + bernP + bernLP,
               data=lda.train.df, CV=T)
pred.lda <- predict(fit.lda, newdata=lda.test.df)

library(e1071)
fit.svm <- svm(formula=Nbloom ~ ord + ordP + ordLP + rf + rf_split + bern + bernP + bernLP,
               data=lda.train.df %>% mutate(Nbloom=factor(Nbloom)), probability=T)
tune_svm <- tune(svm, Nbloom ~ ord + ordP + ordLP + rf + rf_split + bern + bernP + bernLP, 
                 data=lda.train.df %>% mutate(Nbloom=factor(Nbloom)),
                 ranges=list(epsilon=seq(0,1,0.5), cost=2^(seq(5,7,0.5))), probability=T)
pred.svm <- predict(tune_svm$best.model, newdata=lda.test.df %>% mutate(Nbloom=factor(Nbloom)), probability=T)
lda.test.df %>%
  mutate(svm=attr(pred.svm, "probabilities")[,2],
         lda=pred.lda$posterior[,2],
         avg=filter(pred.df, covarSet=="date" & model=="avg" & species=="alexandrium_sp")$pred) %>%
  pivot_longer(3:13, names_to="model", values_to="pred") %>%
  group_by(model) %>%
  summarise(LL=sum(dbinom(Nbloom, 1, pred, log=T))) %>%
  ggplot(aes(model, LL)) + geom_point()

mod_cols <- c(mean="black", mo="grey40", xy="grey60", mo_xy="grey80",
              bern="dodgerblue", bernLP="dodgerblue", bernP="dodgerblue",
              ord="cadetblue", ordLP="cadetblue", ordP="cadetblue",
              rf="green4", rf_split="green3",
              avg="red", avgMoStB1="red3",
              avgMo="purple", avgSt="purple3", avgB1="purple4",
              avgMoB1="yellow", avgStB1="yellow3", avgMoSt="yellow4")
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  mutate(pred=pmax(pmin(pred, 1-1e-7), 1e-7)) %>%
  group_by(species, covarSet, model, modType) %>%
  summarise(LL=sum(dbinom(Nbloom, 1, pred, log=T))) %>%
  group_by(species, covarSet) %>%
  arrange(species, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  ggplot(aes(model, R2, colour=model, group=model)) +
  geom_point() +
  facet_grid(species~covarSet) + 
  scale_y_continuous(limits=c(-0.1, 0.5), breaks=c(0, 0.2, 0.4)) +
  scale_colour_manual(values=mod_cols) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), panel.grid.minor.y=element_blank())
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  mutate(pred=pmax(pmin(pred, 1-1e-7), 1e-7)) %>%
  mutate(Bloom_1=if_else(Nbloom1==0, "from No Bloom", "from Bloom")) %>%
  group_by(species, covarSet, model, modType, Bloom_1) %>%
  summarise(LL=sum(dbinom(Nbloom, 1, pred, log=T))) %>%
  group_by(species, covarSet, Bloom_1) %>%
  arrange(species, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  ggplot(aes(model, R2, colour=model, group=model, shape=Bloom_1)) +
  geom_point() +
  facet_grid(species~covarSet) + 
  ylim(-0.1, 0.5) +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19,1)) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  mutate(pred=pmax(pmin(pred, 1-1e-7), 1e-7)) %>%
  group_by(species, covarSet, model, modType) %>%
  summarise(LL=sum(dbinom(Nbloom, 1, pred, log=T))) %>%
  group_by(species, covarSet) %>%
  arrange(species, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  ggplot(aes(model, R2, colour=covarSet, group=covarSet)) +
  geom_point() + geom_line() +
  facet_grid(species~.) + 
  scale_y_continuous(limits=c(-0.1, 0.5), breaks=c(0, 0.2, 0.4)) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), panel.grid.minor.y=element_blank())

pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, modType) %>%
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
  ggplot(aes(modType, Brier, colour=model, shape=modType)) + 
  geom_point(position=position_jitter(width=0.1), size=2) + 
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(1,15,16,17,18,4,2)) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5)) +
  facet_wrap(~species, scales="free_y")

pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model) %>%
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
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
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
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
  mutate(pred=pmax(pmin(pred, 1-1e-7), 1e-7)) %>%
  group_by(species, model, modType) %>%
  summarise(LL=sum(dbinom(Nbloom, 1, pred, log=T))) %>%
  group_by(species) %>%
  arrange(species, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  ggplot(aes(modType, R2, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_point(position=position_jitter(width=0.05)) + 
  scale_colour_manual(values=mod_cols) +
  labs(x="", y="McFadden's R2\nrelative to site-monthly avg") +
  ylim(-0.1, NA) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5)) +
  facet_wrap(~species)
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  mutate(pred=pmax(pmin(pred, 1-1e-7), 1e-7)) %>%
  group_by(species, model, modType) %>%
  summarise(LL=sum(dbinom(Nbloom, 1, pred, log=T))) %>%
  group_by(species) %>%
  arrange(species, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  ggplot(aes(model, R2, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_point() + 
  scale_colour_manual(values=mod_cols) +
  labs(x="", y="McFadden's R2 relative to grand mean") +
  ylim(-0.1, NA) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5)) +
  facet_grid(species~.)

# By Bloom_1
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, modType, Nbloom1, Nbloom) %>%
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
  mutate(Bloom=if_else(Nbloom==0, "to No Bloom", "to Bloom"),
         Bloom_1=if_else(Nbloom1==0, "from No Bloom", "from Bloom")) %>%
  ggplot(aes(species, Brier, colour=model, shape=modType)) + 
  geom_point(position=position_jitter(width=0.1), size=2) + 
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(1,15,16,17,18,4,2)) +
  facet_grid(Bloom_1~Bloom) + 
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, modType, Nbloom1) %>%
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
  mutate(Bloom_1=if_else(Nbloom1==0, "from No Bloom", "from Bloom")) %>%
  ggplot(aes(species, Brier, colour=model, shape=modType)) + 
  geom_point(position=position_jitter(width=0.1), size=2) + 
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(1,15,16,17,18,4,2)) +
  facet_wrap(~Bloom_1) + 
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, Nbloom1, Nbloom) %>%
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
  mutate(Bloom=if_else(Nbloom==0, "to No Bloom", "to Bloom"),
         Bloom_1=if_else(Nbloom1==0, "from No Bloom", "from Bloom")) %>%
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
  group_by(species, model, modType, Nbloom1, Nbloom) %>%
  summarise(LL=sum(dbinom(Nbloom, 1, pred, log=T))) %>%
  mutate(Bloom=if_else(Nbloom==0, "to No Bloom", "to Bloom"),
         Bloom_1=if_else(Nbloom1==0, "from No Bloom", "from Bloom")) %>%
  group_by(species, Bloom, Bloom_1) %>%
  arrange(species, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  ggplot(aes(species, R2, colour=model, shape=modType)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_point(position=position_jitter(width=0.05)) + 
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(1,15,16,17,18,4,2)) +
  facet_grid(Bloom_1~Bloom, scales="free") + 
  labs(x="", y="McFadden's R2 relative to site-monthly avg") +
  ylim(-0.1, NA) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  mutate(pred=pmax(pmin(pred, 1-1e-5), 1e-5)) %>%
  group_by(species, model, modType, Nbloom1) %>%
  summarise(LL=sum(dbinom(Nbloom, 1, pred, log=T))) %>%
  mutate(Bloom_1=if_else(Nbloom1==0, "from No Bloom", "from Bloom")) %>%
  group_by(species, Bloom_1) %>%
  arrange(species, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  ggplot(aes(species, R2, colour=model, shape=modType)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_point(position=position_jitter(width=0.1), size=2) + 
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(1,15,16,17,18,4,2)) +
  facet_grid(.~Bloom_1, scales="free") + 
  labs(x="", y="McFadden's R2 relative to site-monthly avg") +
  ylim(-0.1, NA) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  mutate(pred=pmax(pmin(pred, 1-1e-5), 1e-5)) %>%
  group_by(species, model, modType, Nbloom) %>%
  summarise(LL=sum(dbinom(Nbloom, 1, pred, log=T))) %>%
  mutate(Bloom=if_else(Nbloom==0, "to No Bloom", "to Bloom")) %>%
  group_by(species, Bloom) %>%
  arrange(species, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  ggplot(aes(species, R2, colour=model, shape=modType)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_point(position=position_jitter(width=0.05)) + 
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(15,16,17,18,4,2)) +
  facet_grid(.~Bloom, scales="free") + 
  labs(x="", y="McFadden's R2 relative to site-monthly avg") +
  ylim(-0.1, NA) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))


# By site
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, siteNum) %>%
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
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
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
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
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
  ggplot(aes(month, Brier, colour=model)) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) +
  scale_x_continuous(breaks=seq(1,12,by=2), labels=month.abb[seq(1,12,by=2)]) +
  facet_wrap(~species, scales="free") 
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, month) %>%
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
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
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
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



binary.ls <- pred.df %>% filter(Nbloom1==0) %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  select(obsid, species, Nbloom, model, pred) %>%
  group_by(model) %>%
  group_split() %>%
  map(~.x %>% group_by(species) %>% group_split)
names(binary.ls) <- names(mod_cols)

roc.ls <- map_depth(binary.ls, 2,
                    ~pROC::roc(.x$Nbloom ~ .x$pred))


for(i in 1:5) {
  par(mfrow=c(1,1))
  png(glue("figs/performance/0x_ROC_{species[i]}.png"), width=5, height=5, res=300, units="in")
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
         col=mod_cols, lty=1, lwd=2, bty="n", cex=0.6)
  dev.off()
}

auc.df <- map_depth(roc.ls, 2, ~.x$auc) %>% 
  map_dfc(~do.call(c, .x)) %>%
  mutate(species=species) %>%
  pivot_longer(-species, names_to="model", values_to="AUC") %>%
  mutate(model=factor(model, levels=names(mod_cols)))
write_csv(auc.df, glue("figs/performance/auc_0x.csv"))

ggplot(auc.df, aes(model, AUC, colour=species, group=species)) + geom_line() + 
  scale_colour_manual(values=sams.cols)
ggplot(auc.df, aes(species, AUC, colour=model, group=model)) + geom_line() + 
  scale_colour_manual(values=mod_cols)
par(mfrow=c(1,1))





auc.f <- dir(glue("figs{sep}performance{sep}"), "auc_.*csv", full.names=T)
auc.df <- map_dfr(auc.f, ~read_csv(.x, show_col_types=F) %>%
                    mutate(type=str_split_fixed(.x, "auc_", 2)[,2])) %>%
  mutate(data_subset=if_else(grepl("0x", type), "0x", "all")) %>%
  mutate(model=factor(model, levels=names(mod_cols)))
ggplot(auc.df, aes(model, AUC, colour=model, group=model)) +
  geom_point() +
  facet_grid(species~data_subset) + 
  scale_colour_manual(values=mod_cols) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

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





dateConditions <- tibble(yday=seq(1, 365, length.out=365)) %>%
  mutate(ydayCos=cos(yday*2*pi/365),
         ydaySin=sin(yday*2*pi/365))



sams.cols <- c("alexandrium_sp"="#2A5883",
               "dinophysis_sp"="#46BFDE",
               "karenia_mikimotoi"="#A9DAE0",
               "prorocentrum_lima"="#E77334",
               "pseudo_nitzschia_sp"="#D21E4C")
nlpars <- c('bydayC', 'bydayS', 'bydaySC', 'btempLwk', 'bshortwaveLwk', 
            'bwindVel', 'binfluxwk', 'bfetch', 'bdinowk', 'bo2wk', 'bpo4wk', 
            'bNbloom1', 'bNbloom2')
out.alex <- readRDS("out/test_full/sord_laplace_alexandrium_sp.rds")
data.alex <- read_csv("out/test_full/dataset_alexandrium_sp.csv")
pred.alex <- map_dfr(nlpars, ~data.alex %>%
                       mutate(slope=colMeans(posterior_epred(out.alex, newdata=., allow_new_levels=T, nlpar=.x)),
                              var=.x))

ggplot(pred.alex, aes(yday, slope)) + 
  geom_hline(yintercept=0, colour="grey30") +
  geom_point(alpha=0.5) + facet_wrap(~var)


out.dino <- readRDS("out/test_full/sord_laplace_dinophysis_sp.rds")
data.dino <- read_csv("out/test_full/dataset_dinophysis_sp.csv")
pred.dino <- map_dfr(nlpars, ~data.dino %>%
                       mutate(slope=colMeans(posterior_epred(out.dino, newdata=., allow_new_levels=T, nlpar=.x)),
                              var=.x))

ggplot(pred.dino, aes(yday, slope)) + 
  geom_hline(yintercept=0, colour="grey30") +
  geom_point(alpha=0.5) + facet_wrap(~var)

pred.da <- bind_rows(pred.alex %>% mutate(species="alexandrium_sp"),
                     pred.dino %>% mutate(species="dinophysis_sp"))
ggplot(pred.da, aes(yday, slope, colour=species)) + 
  geom_hline(yintercept=0, colour="grey30") +
  geom_point(alpha=0.5, shape=1, size=0.7) + facet_wrap(~var) +
  scale_colour_manual(values=sams.cols) +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1)))



out.p <- dir("out/test_full", "ord_p_", full.names=T) %>%
  map(readRDS) %>% setNames(species)
p.effects <- c('tempLwk', 'salinityLwk', 
               'shortwaveLwk', 'windVel', 'waterVelL', 'waterVelR', 'influxwk', 
               'fetch', 'attnwk', 'chlwk', 'dinowk', 'o2wk', 'phwk', 'po4wk', 
               'Nbloom1', 'Nbloom2')
p.nonZero <- map(out.p, ~fixef(.x) %>% 
                   as_tibble(rownames="var") %>%
                   filter(grepl("^p", var)) %>%
                   filter(Q2.5 > 0.01))
library(bayesplot)
imap(out.p, ~mcmc_intervals(.x, regex_pars="^b_p.*Intercept") + ggtitle(.y))

data.alex <- read_csv("out/test_full/dataset_alexandrium_sp.csv")

smooths.df <- bind_cols(expand_grid(yday=seq(1, 365, length.out=52),
                                    siteid=unique(data.alex$siteid)) %>%
                          mutate(ydayCos=cos(yday*2*pi/365),
                                 ydaySin=sin(yday*2*pi/365),
                                 ydaySC=ydayCos*ydaySin),
                        map_dfr(p.effects, ~tibble(var=.x, val=0)) %>%
                          pivot_wider(names_from="var", values_from="val"))
pred.ls <- vector("list", length(species)) %>% setNames(species)
for(i in seq_along(pred.ls)) {
  pred.ls[[i]] <- map_dfr(p.effects, 
                          ~smooths.df %>%
                            mutate(slope=colMeans(posterior_epred(out.p[[i]], newdata=., 
                                                                  allow_new_levels=T,
                                                                  nlpar=paste0("b", .x))),
                                   p=colMeans(posterior_epred(out.p[[i]], newdata=., 
                                                              allow_new_levels=T,
                                                              nlpar=paste0("p", .x))),
                                   var=.x,
                                   species=species[i]))
}
pred.smooths <- do.call(rbind, pred.ls)
ggplot(pred.smooths, aes(yday, slope*p, colour=species, group=paste(species, siteid))) + 
  geom_hline(yintercept=0, colour="grey30", linetype=2) +
  geom_line(alpha=0.25) + facet_wrap(~var) +
  # scale_colour_manual(values=sams.cols) +
  scale_colour_brewer(type="qual", palette=2) +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1))) +
  theme(panel.grid=element_blank())






out.p <- dir("out/test_full", "ordP_", full.names=T) %>%
  map(readRDS) %>% setNames(c("autoreg", "cprn", "external", "local"))
p.effects <- c('tempLwk', 'salinityLwk', 
               'shortwaveLwk', 'windVel', 'waterVelL', 'waterVelR', 'influxwk', 
               'fetch', 'attnwk', 'chlwk', 'dinowk', 'o2wk', 'phwk', 'po4wk', 
               'Nbloom1', 'Nbloom2', 'NlnWt1', 'NlnWt2')
p.all <- map_dfr(out.p, ~fixef(.x) %>% 
                       as_tibble(rownames="var") %>%
                       filter(grepl("^p", var)))
library(bayesplot)
imap(out.p, ~mcmc_intervals(.x, regex_pars="^b_p.*Intercept") + ggtitle(.y))

smooths.df <- bind_cols(expand_grid(yday=seq(1, 365, length.out=52),
                                    siteid=unique(data.alex$siteid)) %>%
                          mutate(ydayCos=cos(yday*2*pi/365),
                                 ydaySin=sin(yday*2*pi/365),
                                 ydaySC=ydayCos*ydaySin),
                        map_dfr(p.effects, ~tibble(var=.x, val=0)) %>%
                          pivot_wider(names_from="var", values_from="val"))
pred.ls <- vector("list", length(out.p)) %>% setNames(names(out.p))
p.effects <- list(autoreg=c("Nbloom1", "Nbloom2", "NlnWt1", "NlnWt2"),
                  cprn=c("attnwk", "chlwk", "dinowk", "o2wk", "phwk", "po4wk"),
                  external=c("windVel", "waterVelL", "waterVelR", "influxwk", "fetch"),
                  local=c("tempLwk", "salinityLwk", "shortwaveLwk"))
for(i in seq_along(pred.ls)) {
  pred.ls[[i]] <- map_dfr(p.effects[[i]], 
                          ~smooths.df %>%
                            mutate(slope=colMeans(posterior_epred(out.p[[i]], newdata=., 
                                                                  allow_new_levels=T,
                                                                  nlpar=paste0("b", .x))),
                                   p=colMeans(posterior_epred(out.p[[i]], newdata=., 
                                                              allow_new_levels=T,
                                                              nlpar=paste0("p", .x))),
                                   var=.x,
                                   covarSet=names(out.p)[i]))
}
pred.smooths <- do.call(rbind, pred.ls)
ggplot(pred.smooths, aes(yday, slope*p, colour=covarSet, group=paste(covarSet, siteid))) + 
  geom_hline(yintercept=0, colour="grey30", linetype=2) +
  geom_line(alpha=0.25) + facet_wrap(~var) +
  # scale_colour_manual(values=sams.cols) +
  scale_colour_brewer(type="qual", palette=2) +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1))) +
  theme(panel.grid=element_blank())
