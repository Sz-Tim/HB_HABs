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

# species <- str_sub(dir("out/test_full", "dataset_all"), 9, -5)
species <- str_sub(dir("out/test_full", "dataset_all"), 13, -5)


# load --------------------------------------------------------------------

data.df <- dir("out/test_full", "dataset_all.*csv", full.names=T) %>% 
  map_dfr(~read_csv(.x, show_col_type=F) %>%
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
  # group_by(species, covarSet, siteid, month) %>%
  # mutate(mo_xy_mnpr=cummean(Nbloom1)) %>%
  # group_by(species, covarSet, siteid) %>%
  # mutate(xy_mnpr=cummean(Nbloom1)) %>%
  ungroup


cv.df <- map_dfr(dir("out/test_full/", "^CV.*csv"), 
                 ~read_csv(glue("out/test_full/{.x}"), show_col_types=F) %>%
                   select(covarSet, siteid, date, obsid, Nbloom, Nbloom1, contains("_mnpr")) %>%
                   mutate(month=lubridate::month(date), 
                          species=str_sub(str_split_fixed(.x, "_", 3)[,3], 1, -5),
                          # some RF pr == 0 | 1
                          across(ends_with("mnpr"), ~pmin(pmax(.x, 1e-6), 1-1e-6))))
fit.df <- map_dfr(dir("out/test_full/", "^fit_all.*csv"), 
              ~read_csv(glue("out/test_full/{.x}"), show_col_types=F) %>%
                select(covarSet, siteid, date, obsid, Nbloom, Nbloom1, contains("_mnpr")) %>%
                mutate(month=lubridate::month(date), 
                       species=str_sub(str_split_fixed(.x, "_", 3)[,3], 1, -5),
                       across(ends_with("mnpr"), ~pmin(pmax(.x, 1e-6), 1-1e-6))))
pred.df <- map_dfr(dir("out/test_full/", "^pred_all.*csv"), 
                   ~read_csv(glue("out/test_full/{.x}"), show_col_types=F) %>%
                     select(covarSet, siteid, date, obsid, Nbloom, Nbloom1, contains("_mnpr")) %>%
                     mutate(month=lubridate::month(date), 
                            species=str_sub(str_split_fixed(.x, "_", 3)[,3], 1, -5),
                            across(ends_with("mnpr"), ~pmin(pmax(.x, 1e-6), 1-1e-6))))


# calc weights ------------------------------------------------------------


LL.tot <- cv.df %>% 
  group_by(covarSet, species) %>% 
  summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
  mutate(across(ends_with("wt"), ~.x/rowSums(across(ends_with("wt"))))) %>% 
  ungroup
LL.tot <- bind_rows(LL.tot, tibble(covarSet="all", species=species)) %>% 
  mutate(across(ends_with("wt"), ~replace_na(.x, 1/8)))

LL.mo <- cv.df %>% 
  group_by(covarSet, species, month) %>% 
  summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
  mutate(across(ends_with("wt"), ~.x/rowSums(across(ends_with("wt"))))) %>%
  ungroup
LL.st <- cv.df %>% 
  group_by(covarSet, species, siteid) %>% 
  summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
  mutate(across(ends_with("wt"), ~.x/rowSums(across(ends_with("wt"))))) %>%
  ungroup
LL.B1 <- cv.df %>% 
  group_by(covarSet, species, Nbloom1) %>% 
  summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
  mutate(across(ends_with("wt"), ~.x/rowSums(across(ends_with("wt"))))) %>%
  ungroup
# LL.moB1 <- cv.df %>% 
#   group_by(covarSet, species, month, Nbloom1) %>% 
#   summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
#   mutate(across(ends_with("wt"), ~.x/rowSums(across(ends_with("wt"))))) %>%
#   ungroup
# LL.stB1 <- cv.df %>% 
#   group_by(covarSet, species, siteid, Nbloom1) %>% 
#   summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
#   mutate(across(ends_with("wt"), ~.x/rowSums(across(ends_with("wt"))))) %>%
#   ungroup
# LL.moSt <- cv.df %>% 
#   group_by(covarSet, species, month, siteid) %>% 
#   summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
#   mutate(across(ends_with("wt"), ~.x/rowSums(across(ends_with("wt"))))) %>%
#   ungroup
# LL.moStB1 <- cv.df %>% 
#   group_by(covarSet, species, month, siteid, Nbloom1) %>% 
#   summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
#   mutate(across(ends_with("wt"), ~.x/rowSums(across(ends_with("wt"))))) %>%
#   ungroup


pred.ens <- pred.df %>% 
  left_join(LL.tot) %>% rowwise() %>%
  mutate(avg=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  ungroup() %>% select(-ends_with("wt"))# %>% 
  # left_join(LL.mo) %>% rowwise() %>%
  # mutate(avgMo=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  # ungroup() %>% select(-ends_with("wt")) %>% 
  # left_join(LL.st) %>% rowwise() %>%
  # mutate(avgSt=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  # ungroup() %>% select(-ends_with("wt")) %>% 
  # left_join(LL.B1) %>% rowwise() %>%
  # mutate(avgB1=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  # ungroup() %>% select(-ends_with("wt")) %>% 
  # left_join(LL.moB1) %>% rowwise() %>%
  # mutate(avgMoB1=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  # ungroup() %>% select(-ends_with("wt")) %>% 
  # left_join(LL.stB1) %>% rowwise() %>%
  # mutate(avgStB1=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  # ungroup() %>% select(-ends_with("wt")) %>% 
  # left_join(LL.moSt) %>% rowwise() %>%
  # mutate(avgMoSt=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  # ungroup() %>% select(-ends_with("wt")) %>%
  # left_join(LL.moStB1) %>% rowwise() %>%
  # mutate(avgMoStB1=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  # ungroup() %>% select(-ends_with("wt"))
pred.df <- pred.ens %>% 
  # mutate(avgSt=if_else(is.na(avgSt), avg, avgSt),
  #        avgMo=if_else(is.na(avgMo), avg, avgMo)#,
         # avgStB1=if_else(is.na(avgStB1), avgB1, avgStB1),
         # avgMoSt=if_else(is.na(avgMoSt), avgMo, avgMoSt),
         # avgMoB1=if_else(is.na(avgMoB1), avgB1, avgMoB1),
         # avgMoStB1=if_else(is.na(avgMoStB1), avgMoB1, avgMoStB1)
         # ) %>%
  rename(avg_mnpr=avg#, 
         # avgMo_mnpr=avgMo,
         # avgSt_mnpr=avgSt,
         # avgB1_mnpr=avgB1,
         # avgMoB1_mnpr=avgMoB1,
         # avgStB1_mnpr=avgStB1,
         # avgMoSt_mnpr=avgMoSt,
         # avgMoStB1_mnpr=avgMoStB1
         ) %>%
  group_by(species, covarSet) %>%
  left_join(., data.df %>% select(species, covarSet, obsid, ends_with("_mnpr")),
            by=c("species", "covarSet", "obsid")) %>%
  pivot_longer(ends_with("_mnpr"), names_to="model", values_to="pred") %>%
  mutate(model=str_remove(model, "_mnpr"),
         siteNum=as.numeric(as.factor(siteid))) %>%
  mutate(modType=case_when(grepl("mo|xy|mean", model) ~ "Null",
                           grepl("rf", model) ~ "Random Forest",
                           grepl("ord", model) ~ "Bayes: Ordinal",
                           grepl("bern", model) ~ "Bayes: Logistic",
                           grepl("avg", model) ~ "Ensemble"))
  # mutate(modType=case_when(grepl("mo|xy|mean", model) ~ "Null",
  #                          grepl("LP", model) ~ "Laplace",
  #                          grepl("rf", model) ~ "RF",
  #                          model %in% c("bern", "ord") ~ "intrct",
  #                          model %in% c("bernP", "ordP") ~ "indVar",
  #                          grepl("avg", model) ~ "ens"))

fit.ens <- fit.df %>% 
  left_join(LL.tot) %>% rowwise() %>%
  mutate(avg=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  ungroup() %>% select(-ends_with("wt")) #%>% 
  # left_join(LL.mo) %>% rowwise() %>%
  # mutate(avgMo=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  # ungroup() %>% select(-ends_with("wt")) %>% 
  # left_join(LL.st) %>% rowwise() %>%
  # mutate(avgSt=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  # ungroup() %>% select(-ends_with("wt")) %>% 
  # left_join(LL.B1) %>% rowwise() %>%
  # mutate(avgB1=sum(c_across(ends_with("wt")) * c_across(ends_with("mnpr")))) %>%
  # ungroup() %>% select(-ends_with("wt")) 
fit.df_ <- fit.ens %>% 
  # mutate(avgSt=if_else(is.na(avgSt), avg, avgSt),
  #        avgMo=if_else(is.na(avgMo), avg, avgMo)
  # ) %>%
  rename(avg_mnpr=avg#, 
         # avgMo_mnpr=avgMo,
         # avgSt_mnpr=avgSt,
         # avgB1_mnpr=avgB1
  ) %>%
  group_by(species, covarSet) %>%
  left_join(., data.df %>% select(species, covarSet, obsid, ends_with("_mnpr")),
            by=c("species", "covarSet", "obsid")) %>%
  pivot_longer(ends_with("_mnpr"), names_to="model", values_to="pred") %>%
  mutate(model=str_remove(model, "_mnpr"),
         siteNum=as.numeric(as.factor(siteid))) %>%
  mutate(modType=case_when(grepl("mo|xy|mean", model) ~ "Null",
                           grepl("rf", model) ~ "Random Forest",
                           grepl("ord", model) ~ "Bayes: Ordinal",
                           grepl("bern", model) ~ "Bayes: Logistic",
                           grepl("avg", model) ~ "Ensemble"))
  # mutate(modType=case_when(grepl("mo|xy|mean", model) ~ "Null",
  #                          grepl("LP", model) ~ "Laplace",
  #                          grepl("rf", model) ~ "RF",
  #                          model %in% c("bern", "ord") ~ "intrct",
  #                          model %in% c("bernP", "ordP") ~ "indVar",
  #                          grepl("avg", model) ~ "ens"))
  

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



# R2 plots ----------------------------------------------------------------

mod_cols <- c(mean="black", mo="grey40", xy="grey60", mo_xy="grey80",
              bern="dodgerblue", bernLP="dodgerblue", bernP="dodgerblue",
              ord="cadetblue", ordLP="cadetblue", ordP="cadetblue",
              rf="green4", rf_split="green3",
              avg="red", #avgMoStB1="red3",
              avgMo="purple", avgSt="purple3", avgB1="purple4")#,
              # avgMoB1="yellow", avgStB1="yellow3", avgMoSt="yellow4")
mod.i <- tibble(mod_col=names(mod_cols),
                mod_clean=c("Grand mean", "Monthly mean", 
                            "Site mean", "Monthly-site mean",
                            paste0("Logistic-", 1:3),
                            paste0("Ordinal-", 1:3),
                            "RF-1", "RF-2",
                            paste0("Ensemble-", 1:4)),
                mod_short=c("Mean", "Month", "Site", "Month-Site",
                            "L-1", "L-2", "L-3",
                            "O-1", "O-2", "O-3",
                            "RF-1", "RF-2", 
                            "Ens-1", "Ens-2", "Ens-3", "Ens-4"))
modType_cols <- c("Null"="grey40", "Random Forest"="#d95f02", 
                  "Bayes: Ordinal"="#1b9e77", "Bayes: Logistic"="#7570b3",
                  "Ensemble"="#e7298a")

bind_rows(fit.df_ %>% mutate(set="Fitted"), 
          pred.df %>% mutate(set="Out-of-sample")) %>%
  ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols), labels=mod.i$mod_clean),
         pred=pmax(pmin(pred, 1-1e-7), 1e-7)) %>%
  group_by(species, set, model, modType) %>%
  summarise(LL=sum(dbinom(Nbloom, 1, pred, log=T))) %>%
  group_by(species, set) %>%
  arrange(species, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  ggplot(aes(model, R2, colour=modType, group=model)) +
  geom_point() + geom_rug(sides="l") +
  facet_grid(set~species) + 
  scale_y_continuous(limits=c(0, 0.6), breaks=c(0, 0.2, 0.4, 0.6)) +
  scale_colour_manual("Model type", values=modType_cols) +
  labs(x="Model", y=expression("McFadden's pseudo-R"^2)) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        panel.grid.minor.y=element_blank(),
        legend.position="bottom")
ggsave("figs/R2.jpg", width=10, height=5, units="in")

bind_rows(fit.df_ %>% mutate(set="Fitted"), 
          pred.df %>% mutate(set="Out-of-sample")) %>%
  ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols), labels=mod.i$mod_clean),
         pred=pmax(pmin(pred, 1-1e-7), 1e-7),
         Bloom=if_else(Nbloom==0, "to No Bloom", "to Bloom")) %>%
  group_by(species, set, model, modType, Bloom) %>%
  summarise(LL=sum(dbinom(Nbloom, 1, pred, log=T))) %>%
  group_by(species, set, Bloom) %>%
  arrange(species, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  ggplot(aes(model, R2, colour=modType, group=model, shape=Bloom)) +
  geom_point() + geom_rug(sides="l") +
  facet_grid(set~species) + 
  scale_y_continuous(limits=c(0, 0.7), breaks=c(0, 0.2, 0.4, 0.6)) +
  scale_colour_manual("Model type", values=modType_cols) +
  scale_shape_manual(values=c(19,1)) +
  labs(x="Model", y=expression("McFadden's pseudo-R"^2)) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        panel.grid.minor.y=element_blank(),
        legend.position="bottom")
ggsave("figs/R2_bloom.jpg", width=10, height=5, units="in")

bind_rows(fit.df_ %>% mutate(set="Fitted"), 
          pred.df %>% mutate(set="Out-of-sample")) %>%
  ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols), labels=mod.i$mod_clean),
         pred=pmax(pmin(pred, 1-1e-7), 1e-7),
         Bloom_1=if_else(Nbloom1==0, "from No Bloom", "from Bloom")) %>%
  group_by(species, set, model, modType, Bloom_1) %>%
  summarise(LL=sum(dbinom(Nbloom, 1, pred, log=T))) %>%
  group_by(species, set, Bloom_1) %>%
  arrange(species, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  ggplot(aes(model, R2, colour=modType, group=model, shape=Bloom_1)) +
  geom_point() + geom_rug(sides="l") +
  facet_grid(set~species) + 
  scale_y_continuous(limits=c(0, 0.7), breaks=c(0, 0.2, 0.4, 0.6)) +
  scale_colour_manual("Model type", values=modType_cols) +
  scale_shape_manual(values=c(19,1)) +
  labs(x="Model", y=expression("McFadden's pseudo-R"^2)) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        panel.grid.minor.y=element_blank(),
        legend.position="bottom")
ggsave("figs/R2_bloom1.jpg", width=10, height=5, units="in")


bind_rows(fit.df_ %>% mutate(set="Fitted"), 
          pred.df %>% mutate(set="Out-of-sample")) %>%
  ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols), labels=mod.i$mod_clean),
         pred=pmax(pmin(pred, 1-1e-7), 1e-7)) %>%
  group_by(species, set, model, modType) %>%
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
  ggplot(aes(model, Brier, colour=modType, group=model)) +
  geom_point() + geom_rug(sides="l") +
  facet_grid(set~species) + 
  scale_y_continuous(limits=c(0, 0.25), breaks=c(0, 0.1, 0.2)) +
  scale_colour_manual("Model type", values=modType_cols) +
  labs(x="Model", y="Brier score") +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        panel.grid.minor.y=element_blank(),
        legend.position="bottom")
ggsave("figs/Brier.jpg", width=10, height=5, units="in")

bind_rows(fit.df_ %>% mutate(set="Fitted"), 
          pred.df %>% mutate(set="Out-of-sample")) %>%
  ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols), labels=mod.i$mod_clean),
         pred=pmax(pmin(pred, 1-1e-7), 1e-7),
         Bloom_1=if_else(Nbloom1==0, "from No Bloom", "from Bloom")) %>%
  group_by(species, set, model, modType, Bloom_1) %>%
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
  ggplot(aes(model, Brier, colour=modType, group=model, shape=Bloom_1)) +
  geom_point() + geom_rug(sides="l") +
  facet_grid(set~species) + 
  scale_y_continuous(limits=c(0, 0.5), breaks=c(0, 0.2, 0.4)) +
  scale_colour_manual("Model type", values=modType_cols) +
  scale_shape_manual(values=c(19,1)) +
  labs(x="Model", y="Brier score") +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        panel.grid.minor.y=element_blank(),
        legend.position="bottom")
ggsave("figs/Brier_bloom1.jpg", width=10, height=5, units="in")

bind_rows(fit.df_ %>% mutate(set="Fitted"), 
          pred.df %>% mutate(set="Out-of-sample")) %>%
  ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols), labels=mod.i$mod_clean),
         pred=pmax(pmin(pred, 1-1e-7), 1e-7),
         Bloom=if_else(Nbloom==0, "to No Bloom", "to Bloom")) %>%
  group_by(species, set, model, modType, Bloom) %>%
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
  ggplot(aes(model, Brier, colour=modType, group=model, shape=Bloom)) +
  geom_point() + geom_rug(sides="l") +
  facet_grid(set~species) + 
  scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1)) +
  scale_colour_manual("Model type", values=modType_cols) +
  scale_shape_manual(values=c(19,1)) +
  labs(x="Model", y="Brier score") +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        panel.grid.minor.y=element_blank(),
        legend.position="bottom")
ggsave("figs/Brier_bloom.jpg", width=10, height=5, units="in")



# R2 site -----------------------------------------------------------------

pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, siteNum, covarSet) %>%
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
  ggplot(aes(siteNum, Brier, colour=model)) + 
  geom_point(position=position_jitter(width=0.01)) + 
  scale_colour_manual(values=mod_cols) +
  facet_grid(species~covarSet, scales="free") 
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






# R2 month ----------------------------------------------------------------

pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, month, covarSet) %>%
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
  ggplot(aes(month, Brier, colour=model)) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) +
  scale_x_continuous(breaks=seq(1,12,by=2), labels=month.abb[seq(1,12,by=2)]) +
  facet_grid(species~covarSet, scales="free") 
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, month, covarSet) %>%
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
  group_by(species, month, covarSet) %>%
  arrange(species, model) %>%
  mutate(Improve=first(Brier)-Brier) %>%
  ggplot(aes(month, Improve, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_line(position=position_jitter(width=0.05)) + 
  scale_colour_manual(values=mod_cols) +
  scale_x_continuous(breaks=seq(1,12,by=2), labels=month.abb[seq(1,12,by=2)]) +
  facet_grid(species~covarSet, scales="free")  +
  labs(x="", y="Improvement in Brier score\nrelative to site-monthly avg")
pred.df %>% ungroup %>%
  mutate(model=factor(model, levels=names(mod_cols))) %>%
  group_by(species, model, month, covarSet) %>%
  summarise(Brier=mean((pred-Nbloom)^2, na.rm=T)) %>%
  group_by(species, month, covarSet) %>%
  arrange(species, model) %>%
  mutate(Brier=Brier+0.01) %>%
  mutate(pctImprove=(first(Brier)-Brier)/first(Brier)*100) %>%
  ggplot(aes(month, pctImprove, colour=model)) + 
  geom_hline(yintercept=0, colour="grey", size=0.25) +
  geom_line(position=position_jitter(width=0.05)) + 
  scale_colour_manual(values=mod_cols) +
  scale_x_continuous(breaks=seq(1,12,by=2), labels=month.abb[seq(1,12,by=2)]) +
  ylim(-50,NA) +
  facet_grid(species~covarSet, scales="free")  +
  labs(x="", y="Improvement in Brier score\nrelative to site-monthly avg (%)")








# pred-obs ----------------------------------------------------------------

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




# ROC ---------------------------------------------------------------------

ROC_cols <- mod_cols[names(mod_cols) %in% unique(pred.df$model)]
binary.ls <- fit.df_ %>% ungroup %>% filter(covarSet=="all") %>% filter(Nbloom1==1) %>%
  mutate(model=factor(model, levels=names(ROC_cols))) %>%
  select(obsid, species, Nbloom, model, pred) %>%
  group_by(model) %>%
  group_split() %>%
  map(~.x %>% group_by(species) %>% group_split)
names(binary.ls) <- names(ROC_cols)

roc.ls <- map_depth(binary.ls, 2,
                    ~pROC::roc(.x$Nbloom ~ .x$pred))

# species <- c("alexandrium_sp", "karenia_mikimotoi", "prorocentrum_lima")
for(i in 1:5) {
  par(mfrow=c(1,1))
  png(glue("figs/performance/ROC_fit_1x_{species[i]}.png"), width=5, height=5, res=300, units="in")
  plot(NA, NA, xlim=c(0, 1), ylim=c(0, 1),
       xlab="1 - Specificity", ylab="Sensitivity", axes=F, main=species[i])
  axis(side=1, at=c(0, 0.5, 1))
  axis(side=2, at=c(0, 0.5, 1))
  abline(a=0, b=1, col="grey30", lwd=0.5)
  map(seq_along(roc.ls),
      ~lines(1-roc.ls[[.x]][[i]]$specificities, roc.ls[[.x]][[i]]$sensitivities,
             col=ROC_cols[.x], lwd=2.5))
  legend("bottomright",
         paste(names(ROC_cols), map_dbl(roc.ls, ~round(.x[[i]]$auc, 2)),  sep=": "),
         col=ROC_cols, lty=1, lwd=2, bty="n", cex=0.6)
  dev.off()
}

auc.df <- map_depth(roc.ls, 2, ~.x$auc) %>% 
  map_dfc(~do.call(c, .x)) %>%
  mutate(species=str_remove(species, "all_")) %>%
  pivot_longer(-species, names_to="model", values_to="AUC") %>%
  mutate(model=factor(model, levels=names(ROC_cols)))  %>%
  mutate(modType=case_when(grepl("mo|xy|mean", model) ~ "Null",
                           grepl("LP", model) ~ "Laplace",
                           grepl("rf", model) ~ "RF",
                           model %in% c("bern", "ord") ~ "intrct",
                           model %in% c("bernP", "ordP") ~ "indVar",
                           grepl("avg", model) ~ "ens")) %>%
  mutate(modType=factor(modType, levels=c("Null", "RF",  "intrct", "Laplace", "indVar", "ens")))
write_csv(auc.df, glue("figs/performance/auc_fit_1x_all.csv"))

ggplot(auc.df, aes(model, AUC, colour=species, group=species)) + 
  geom_point() + geom_line() + 
  # scale_colour_manual(values=sams.cols) + 
  facet_grid(.~modType, scales="free_x", space="free_x") + 
  theme(legend.position="bottom") + ylim(0.5, 1)
ggplot(auc.df, aes(species, AUC, colour=model, group=model)) + geom_line() + 
  scale_colour_manual(values=ROC_cols) + ylim(0.5, 1)
par(mfrow=c(1,1))





auc.f <- dir(glue("figs/performance/"), "auc_.*xx.*csv", full.names=T)
auc.df <- map_dfr(auc.f, ~read_csv(.x, show_col_types=F) %>%
                    mutate(type=str_split_fixed(.x, "auc_", 2)[,2])) %>%
  mutate(data_subset=if_else(grepl("fit", type), "fit", "forecast")) %>%
  mutate(modType=case_when(grepl("mo|xy|mean", model) ~ "Null",
                           grepl("rf", model) ~ "Random Forest",
                           grepl("ord", model) ~ "Bayes: Ordinal",
                           grepl("bern", model) ~ "Bayes: Logistic",
                           grepl("avg", model) ~ "Ensemble"),
         model=factor(model, levels=mod.i$mod_col, labels=mod.i$mod_short))
ggplot(auc.df, aes(model, AUC, colour=modType)) +
  geom_point() +
  facet_grid(data_subset~species) + 
  scale_colour_manual(values=modType_cols) + 
  scale_y_continuous(limits=c(0.45, 1), breaks=c(0.5, 0.75, 1)) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        axis.title.x=element_blank(),
        legend.position="bottom",
        legend.title=element_blank())
ggsave("figs/AUC.jpg", width=7.5, height=3.25)








# smooths -----------------------------------------------------------------

sams.cols <- c("alexandrium_sp"="#2A5883",
               "dinophysis_sp"="#46BFDE",
               "karenia_mikimotoi"="#A9DAE0",
               "prorocentrum_lima"="#E77334",
               "pseudo_nitzschia_sp"="#D21E4C")
nlpars <- c('bydayCos', 'bydaySin', 
            'btempLwk', 'bsalinityLwk', 'bshortwaveLwk', 
            'bwindVel', 'bwaterVelL', 'bwaterVelR', 'binfluxwk', 'bfetch', 
            'battnwk', 'bchlwk', 'bdinowk', 'bo2wk', 'bphwk', 'bpo4wk', 
            'bNbloom1', 'bNbloom2', 'bNlnWt1', 'bNlnWt2')



# * ordP ------------------------------------------------------------------
library(tidyverse)
library(glue)
library(lubridate)
library(brms)
theme_set(theme_bw())
data.alex <- read_csv("out/test_full/dataset_all_alexandrium_sp.csv")
species <- names(sams.cols)
out.p <- dir("out/test_full", "ordP_all_", full.names=T) %>%
  map(readRDS) %>% setNames(species)
p.effects <- c('tempLwk', 'salinityLwk', 'shortwaveLwk', 'kmLwk', 'precipLwk',
               'windVel', 'waterVelL', 'waterVelR', 'windLwk', 'waterLwk', 'waterRwk', 
               'fetch',
               'attnwk', 'chlwk', 'dinowk', 'o2wk', 'phwk', 'po4wk',
               'Nbloom1', 'Nbloom2', 'NlnWt1', 'NlnWt2', 'NlnRAvg1', 'NlnRAvg2')
p.nonZero <- map(out.p, ~fixef(.x) %>% 
                   as_tibble(rownames="var") %>%
                   filter(grepl("^p", var)) %>%
                   filter(Q2.5 > 0.01))
# library(bayesplot)
# imap(out.p, ~mcmc_intervals(.x, regex_pars="^b_p.*Intercept") + ggtitle(.y))

smooths.df <- bind_cols(expand_grid(yday=seq(0, 364, length.out=52),
                                    siteid=unique(data.alex$siteid)) %>%
                          mutate(ydayCos=cos(yday*2*pi/365),
                                 ydaySin=sin(yday*2*pi/365),
                                 ydaySC=ydayCos*ydaySin,
                                 date=as_date(yday, origin="2019-01-01")),
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
pred.ordP <- do.call(rbind, pred.ls)
ggplot(pred.ordP, aes(date, slope*p, colour=species, group=paste(species, siteid))) + 
  geom_hline(yintercept=0, colour="grey30", linetype=2) +
  geom_line(alpha=0.25) + facet_wrap(~var, ncol=5) +
  scale_colour_manual(values=sams.cols) +
  scale_x_date(date_labels="%b", date_breaks="2 months") +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1))) +
  theme(panel.grid.minor=element_blank(), legend.position="bottom")
ggsave("figs/smooths_ordP.jpg", width=9, height=9)

vars_nonZero <- pred.ordP %>% 
  group_by(var) %>% 
  summarise(maxB=max(abs(slope*p))) %>%
  filter(maxB > 0.05)

pred.ordP %>% 
  filter(var %in% vars_nonZero$var) %>%
  ggplot(aes(date, slope*p, colour=species, group=paste(species, siteid))) + 
  geom_hline(yintercept=0, colour="grey30", linetype=2) +
  geom_line(alpha=0.25) + facet_wrap(~var) +
  scale_colour_manual(values=sams.cols) +
  scale_x_date(date_labels="%b", date_breaks="2 months") +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1))) +
  theme(panel.grid=element_blank(), legend.position="bottom")




# * bernP01 ----------------------------------------------------------------
library(tidyverse)
library(glue)
library(lubridate)
library(brms)
theme_set(theme_bw())
data.alex <- read_csv("out/test_full/dataset_all_alexandrium_sp.csv")
species <- names(sams.cols)
out.p <- dir("out/test_full", "bernP01_all_", full.names=T) %>%
  map(readRDS) %>% setNames(species)
p.effects <- c('tempLwk', 'salinityLwk', 'shortwaveLwk', 'kmLwk', 'precipLwk',
               'windVel', 'waterVelL', 'waterVelR', 'windLwk', 'waterLwk', 'waterRwk', 
               'fetch',
               'attnwk', 'chlwk', 'dinowk', 'o2wk', 'phwk', 'po4wk',
               'Nbloom1', 'Nbloom2', 'NlnWt1', 'NlnWt2', 'NlnRAvg1', 'NlnRAvg2')
p.nonZero <- map(out.p, ~fixef(.x) %>% 
                   as_tibble(rownames="var") %>%
                   filter(grepl("^p", var)) %>%
                   filter(Q2.5 > 0.01))
# library(bayesplot)
# imap(out.p, ~mcmc_intervals(.x, regex_pars="^b_p.*Intercept") + ggtitle(.y))

smooths.df <- bind_cols(expand_grid(yday=seq(0, 364, length.out=52),
                                    siteid=unique(data.alex$siteid)) %>%
                          mutate(ydayCos=cos(yday*2*pi/365),
                                 ydaySin=sin(yday*2*pi/365),
                                 ydaySC=ydayCos*ydaySin,
                                 date=as_date(yday, origin="2019-01-01")),
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
pred.bernP01 <- do.call(rbind, pred.ls)
ggplot(pred.bernP01, aes(date, slope*p, colour=species, group=paste(species, siteid))) + 
  geom_hline(yintercept=0, colour="grey30", linetype=2) +
  geom_line(alpha=0.25) + facet_wrap(~var, ncol=5) +
  scale_colour_manual(values=sams.cols) +
  scale_x_date(date_labels="%b", date_breaks="2 months") +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1))) +
  theme(panel.grid=element_blank(), legend.position="bottom")
ggsave("figs/smooths_bernP01.jpg", width=9, height=9.5)






# * bernP11 ----------------------------------------------------------------
library(tidyverse)
library(glue)
library(lubridate)
library(brms)
theme_set(theme_bw())
data.alex <- read_csv("out/test_full/dataset_all_alexandrium_sp.csv")
species <- names(sams.cols)
out.p <- dir("out/test_full", "bernP11_all_", full.names=T) %>%
  map(readRDS) %>% setNames(species)
p.effects <- c('tempLwk', 'salinityLwk', 'shortwaveLwk', 'kmLwk', 'precipLwk',
               'windVel', 'waterVelL', 'waterVelR', 'windLwk', 'waterLwk', 'waterRwk', 
               'fetch', 
               'attnwk', 'chlwk', 'dinowk', 'o2wk', 'phwk', 'po4wk',
               'Nbloom1', 'Nbloom2', 'NlnWt1', 'NlnWt2', 'NlnRAvg1', 'NlnRAvg2')
p.nonZero <- map(out.p, ~fixef(.x) %>% 
                   as_tibble(rownames="var") %>%
                   filter(grepl("^p", var)) %>%
                   filter(Q2.5 > 0.01))
# library(bayesplot)
# imap(out.p, ~mcmc_intervals(.x, regex_pars="^b_p.*Intercept") + ggtitle(.y))

smooths.df <- bind_cols(expand_grid(yday=seq(0, 364, length.out=52),
                                    siteid=unique(data.alex$siteid)) %>%
                          mutate(ydayCos=cos(yday*2*pi/365),
                                 ydaySin=sin(yday*2*pi/365),
                                 ydaySC=ydayCos*ydaySin,
                                 date=as_date(yday, origin="2019-01-01")),
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
pred.bernP11 <- do.call(rbind, pred.ls)
ggplot(pred.bernP11, aes(date, slope*p, colour=species, group=paste(species, siteid))) + 
  geom_hline(yintercept=0, colour="grey30", linetype=2) +
  geom_line(alpha=0.25) + facet_wrap(~var, ncol=5) +
  scale_colour_manual(values=sams.cols) +
  scale_x_date(date_labels="%b", date_breaks="2 months") +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1))) +
  theme(panel.grid.minor=element_blank(), legend.position="bottom")
ggsave("figs/smooths_bernP11.jpg", width=9, height=9)






#####
# Compare ordLP vs ordP for one species
data.dino <- read_csv("out/test_full/dataset_autoreg_dinophysis_sp.csv")
out.p <- dir("out/test_full", "P.*autoreg_dinophysis_sp.rds", full.names=T) %>%
  map(readRDS) %>% setNames(c("ordLP", "ordP"))
p.effects <- c('Nbloom1', 'Nbloom2', 'NlnWt1', 'NlnWt2', 'NlnRAvg1', 'NlnRAvg2')
p.all <- map_dfr(out.p, ~fixef(.x) %>% 
                   as_tibble(rownames="var") %>%
                   filter(grepl("^p", var)))
library(bayesplot)
imap(out.p, ~mcmc_intervals(.x, regex_pars="^b_p.*Intercept") + ggtitle(.y))

smooths.df <- bind_cols(expand_grid(yday=seq(0, 364, length.out=52),
                                    siteid=unique(data.alex$siteid)) %>%
                          mutate(ydayCos=cos(yday*2*pi/365),
                                 ydaySin=sin(yday*2*pi/365),
                                 ydaySC=ydayCos*ydaySin,
                                 date=as_date(yday, origin="2019-01-01")),
                        map_dfr(p.effects, ~tibble(var=.x, val=0)) %>%
                          pivot_wider(names_from="var", values_from="val"))
pred.ls <- vector("list", length(out.p)) %>% setNames(names(out.p))
for(i in seq_along(pred.ls)) {
  if(grepl("LP", names(pred.ls)[i])) {
    pred.ls[[i]] <- map_dfr(p.effects, 
                            ~smooths.df %>%
                              mutate(slope=colMeans(posterior_epred(out.p[[i]], newdata=., 
                                                                    allow_new_levels=T,
                                                                    nlpar=paste0("b", .x))),
                                     var=.x,
                                     model=names(out.p)[i]))
  } else {
    pred.ls[[i]] <- map_dfr(p.effects, 
                            ~smooths.df %>%
                              mutate(b=colMeans(posterior_epred(out.p[[i]], newdata=., 
                                                                    allow_new_levels=T,
                                                                    nlpar=paste0("b", .x))),
                                     p=colMeans(posterior_epred(out.p[[i]], newdata=., 
                                                                allow_new_levels=T,
                                                                nlpar=paste0("p", .x))),
                                     var=.x,
                                     model=names(out.p)[i])  %>%
                              mutate(slope=b*p) %>%
                              select(-b, -p)) 
  }
}
pred.smooths <- do.call(rbind, pred.ls)
ggplot(pred.smooths, aes(date, slope, colour=model, group=paste(model, siteid))) + 
  geom_hline(yintercept=0, colour="grey30", linetype=2) +
  geom_line(alpha=0.25) + facet_wrap(~var, ncol=2) +
  scale_colour_brewer(type="qual", palette=2) +
  scale_x_date(date_labels="%b", date_breaks="2 months") +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1))) +
  theme(panel.grid.minor=element_blank())
