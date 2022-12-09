# HABReports Bayesian modelling
# Bayesian HAB output processing
# Tim Szewczyk


# This script reads and processes the output from 02*_fit_*.R


# setup -------------------------------------------------------------------

library(tidyverse)
library(glue)
library(lubridate)
library(brms)
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
walk(dir("code", "*00_fn", full.names=T), source)


out.dir <- "out/full/"

sp.i <- read_csv("data/sp_i.csv")



# load output -------------------------------------------------------------

fit.df <- map_dfr(dir(out.dir, "^fit.*csv"), 
                   ~read_csv(glue("{out.dir}/{.x}"), show_col_types=F) %>%
                     select(siteid, date, obsid, Nbloom, Nbloom1, contains("_mnpr")) %>%
                     mutate(month=month(date), 
                            species=str_sub(str_split_fixed(.x, "_", 4)[,4], 1, -5),
                            across(ends_with("mnpr"), ~pmin(pmax(.x, 1e-6), 1-1e-6))) %>%
                     pivot_longer(ends_with("mnpr"), names_to="model", values_to="pred")) %>%
  pivot_wider(names_from="model", values_from="pred")

pred.df <- map_dfr(dir(out.dir, "^pred.*csv"), 
                   ~read_csv(glue("{out.dir}/{.x}"), show_col_types=F) %>%
                     select(siteid, date, obsid, Nbloom, Nbloom1, contains("_mnpr")) %>%
                     mutate(month=month(date), 
                            species=str_sub(str_split_fixed(.x, "_", 4)[,4], 1, -5),
                            across(ends_with("mnpr"), ~pmin(pmax(.x, 1e-6), 1-1e-6))) %>%
                     pivot_longer(ends_with("mnpr"), names_to="model", values_to="pred")) %>%
  pivot_wider(names_from="model", values_from="pred")

cv.df <- map_dfr(dir(out.dir, "^CV.*csv"), 
                  ~read_csv(glue("{out.dir}/{.x}"), show_col_types=F) %>%
                    select(siteid, date, obsid, Nbloom, Nbloom1, contains("_mnpr")) %>%
                    mutate(month=month(date), 
                           species=str_sub(str_split_fixed(.x, "_", 4)[,4], 1, -5),
                           across(ends_with("mnpr"), ~pmin(pmax(.x, 1e-6), 1-1e-6))) %>%
                    pivot_longer(ends_with("mnpr"), names_to="model", values_to="pred")) %>%
  pivot_wider(names_from="model", values_from="pred")



# null models -------------------------------------------------------------

null_grand.df <- fit.df %>% 
  group_by(species) %>%
  summarise(mn=mean(Nbloom)) %>%
  ungroup
null_4wk.df <- expand_grid(species=sp.i$full, yday=1:365) %>% 
  mutate(mn=NA)
for(i in 1:nrow(null_4wk.df)) {
  null_4wk.df$mn[i] <- mean(
    filter(fit.df, 
           species==null_4wk.df$species[i],
           (abs(yday(date) - null_4wk.df$yday[i]) < (7*2) |
              abs(yday(date) - 365 - null_4wk.df$yday[i]) < (7*2) |
              abs(yday(date) + 365 - null_4wk.df$yday[i]) < (7*2)))$Nbloom)
}
null_4wk.df <- null_4wk.df %>% ungroup



# ensemble models ---------------------------------------------------------

# calculate LL-based weights
wt.tot <- cv.df %>% 
  group_by(species) %>% 
  summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
  mutate(across(ends_with("wt"), ~.x^2/rowSums(across(ends_with("wt"))^2))) %>% 
  ungroup
wt.month <- cv.df %>% 
  group_by(species, month) %>% 
  summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
  mutate(across(ends_with("wt"), ~.x^2/rowSums(across(ends_with("wt"))^2))) %>% 
  ungroup
wt.B1 <- cv.df %>% 
  group_by(species, Nbloom1) %>% 
  summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
  mutate(across(ends_with("wt"), ~.x^2/rowSums(across(ends_with("wt"))^2))) %>% 
  ungroup
wt.HB <- cv.df %>% 
  group_by(species, Nbloom1) %>% 
  summarise(across(ends_with("mnpr"), ~1/sum(-dbinom(Nbloom, 1, .x, log=T)), .names="{.col}_wt")) %>%
  mutate(across(contains("svm"), ~0), across(contains("xgb"), ~0), across(contains("rf"), ~0),
         across(contains("ord_"), ~0), across(contains("bern_"), ~0)) %>%
  mutate(across(ends_with("wt"), ~.x^2/rowSums(across(ends_with("wt"))^2))) %>% 
  ungroup

fit.df_ <- full_join(
  fit.df, 
  left_join(
    fit.df %>% 
      pivot_longer(ends_with("mnpr"), names_to="model", values_to="pred") %>%
      mutate(model=str_remove(model, "_mnpr")), 
    wt.tot %>% 
      pivot_longer(ends_with("wt"), names_to="model", values_to="wt") %>%
      mutate(model=str_remove(model, "_mnpr_wt"))) %>%
    group_by(obsid, species) %>%
    summarise(avg_mnpr=sum(pred*wt, na.rm=T))) %>%
  left_join(.,
            left_join(
              fit.df %>% 
                pivot_longer(ends_with("mnpr"), names_to="model", values_to="pred") %>%
                mutate(model=str_remove(model, "_mnpr")), 
              wt.month %>% 
                pivot_longer(ends_with("wt"), names_to="model", values_to="wt") %>%
                mutate(model=str_remove(model, "_mnpr_wt"))) %>%
              group_by(obsid, species) %>%
              summarise(avgMo_mnpr=sum(pred*wt, na.rm=T))) %>%
  left_join(.,
            left_join(
              fit.df %>% 
                pivot_longer(ends_with("mnpr"), names_to="model", values_to="pred") %>%
                mutate(model=str_remove(model, "_mnpr")), 
              wt.B1 %>% 
                pivot_longer(ends_with("wt"), names_to="model", values_to="wt") %>%
                mutate(model=str_remove(model, "_mnpr_wt"))) %>%
              group_by(obsid, species) %>%
              summarise(avgB1_mnpr=sum(pred*wt, na.rm=T))) %>%
  left_join(.,
            left_join(
              fit.df %>% 
                pivot_longer(ends_with("mnpr"), names_to="model", values_to="pred") %>%
                mutate(model=str_remove(model, "_mnpr")), 
              wt.HB %>% 
                pivot_longer(ends_with("wt"), names_to="model", values_to="wt") %>%
                mutate(model=str_remove(model, "_mnpr_wt"))) %>%
              group_by(obsid, species) %>%
              summarise(avgHB_mnpr=sum(pred*wt, na.rm=T))) %>%
  mutate(yday=yday(date)) %>%
  left_join(., null_grand.df %>% rename(grand_mnpr=mn)) %>%
  left_join(., null_4wk.df %>% rename(fourWk_mnpr=mn)) %>%
  pivot_longer(ends_with("_mnpr"), names_to="model", values_to="pred") %>%
  mutate(model=str_remove(model, "_mnpr"),
         siteNum=as.numeric(as.factor(siteid))) %>%
  mutate(modType=case_when(grepl("fourWk|grand", model) ~ "Null",
                           grepl("glm", model) ~ "GLM",
                           grepl("rf", model) ~ "ML: RF",
                           grepl("^svm", model) ~ "ML: SVM",
                           grepl("^xgb", model) ~ "ML: XGB",
                           grepl("ord", model) ~ "Bayes: Ordinal",
                           grepl("bern", model) ~ "Bayes: Logistic",
                           grepl("avg", model) ~ "Ensemble"))

pred.df_ <- full_join(
  pred.df, 
  left_join(
    pred.df %>% 
      pivot_longer(ends_with("mnpr"), names_to="model", values_to="pred") %>%
      mutate(model=str_remove(model, "_mnpr")), 
    wt.tot %>% 
      pivot_longer(ends_with("wt"), names_to="model", values_to="wt") %>%
      mutate(model=str_remove(model, "_mnpr_wt"))) %>%
    group_by(obsid, species) %>%
    summarise(avg_mnpr=sum(pred*wt, na.rm=T))) %>%
  left_join(.,
            left_join(
              pred.df %>% 
                pivot_longer(ends_with("mnpr"), names_to="model", values_to="pred") %>%
                mutate(model=str_remove(model, "_mnpr")), 
              wt.month %>% 
                pivot_longer(ends_with("wt"), names_to="model", values_to="wt") %>%
                mutate(model=str_remove(model, "_mnpr_wt"))) %>%
              group_by(obsid, species) %>%
              summarise(avgMo_mnpr=sum(pred*wt, na.rm=T))) %>%
  left_join(.,
            left_join(
              pred.df %>% 
                pivot_longer(ends_with("mnpr"), names_to="model", values_to="pred") %>%
                mutate(model=str_remove(model, "_mnpr")), 
              wt.B1 %>% 
                pivot_longer(ends_with("wt"), names_to="model", values_to="wt") %>%
                mutate(model=str_remove(model, "_mnpr_wt"))) %>%
              group_by(obsid, species) %>%
              summarise(avgB1_mnpr=sum(pred*wt, na.rm=T))) %>%
  left_join(.,
            left_join(
              pred.df %>% 
                pivot_longer(ends_with("mnpr"), names_to="model", values_to="pred") %>%
                mutate(model=str_remove(model, "_mnpr")), 
              wt.HB %>% 
                pivot_longer(ends_with("wt"), names_to="model", values_to="wt") %>%
                mutate(model=str_remove(model, "_mnpr_wt"))) %>%
              group_by(obsid, species) %>%
              summarise(avgHB_mnpr=sum(pred*wt, na.rm=T))) %>%
  mutate(yday=yday(date)) %>%
  left_join(., null_grand.df %>% rename(grand_mnpr=mn)) %>%
  left_join(., null_4wk.df %>% rename(fourWk_mnpr=mn)) %>%
  pivot_longer(ends_with("_mnpr"), names_to="model", values_to="pred") %>%
  mutate(model=str_remove(model, "_mnpr"),
         siteNum=as.numeric(as.factor(siteid))) %>%
  mutate(modType=case_when(grepl("fourWk|grand", model) ~ "Null",
                           grepl("glm", model) ~ "GLM",
                           grepl("rf", model) ~ "ML: RF",
                           grepl("^svm", model) ~ "ML: SVM",
                           grepl("^xgb", model) ~ "ML: XGB",
                           grepl("ord", model) ~ "Bayes: Ordinal",
                           grepl("bern", model) ~ "Bayes: Logistic",
                           grepl("avg", model) ~ "Ensemble"))

bind_rows(fit.df_ %>% mutate(dataSubset="fit"),
          pred.df_ %>% mutate(dataSubset="oos")) %>%
  write_csv(glue("{out.dir}/full_mean_out.csv"))



