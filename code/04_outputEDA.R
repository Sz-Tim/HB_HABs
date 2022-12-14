# HABReports Bayesian modelling
# Bayesian HAB output processing
# Tim Szewczyk


# This script explores the output from 03_processOut.R


# setup -------------------------------------------------------------------

library(tidyverse)
library(glue)
library(lubridate)
library(brms)
# library(gganimate)
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
walk(dir("code", "*00_fn", full.names=T), source)


prior_type <- c("1-loose", "2-medium", "3-tight")[2]
out.dir <- glue("out/{prior_type}/")

sp.i <- read_csv("data/sp_i.csv")

bloomThresh <- dir(out.dir, "dataset.*csv", full.names=T) %>% 
  map(read_csv) %>%
  map(~max((!.x$Nbloom)*.x$NcatNum))




# colours, etc ------------------------------------------------------------

mod_cols <- c(grand="black", fourWk="grey40", 
              glmRidge="grey60", glmRidge_split="grey60", glmLasso="grey80", glmLasso_split="grey80",
              bern="dodgerblue", bernP="dodgerblue", 
              ord="cadetblue", ordP="cadetblue", 
              rf="green4", rf_split="green3",
              svm="olivedrab", svm_split="olivedrab4", 
              xgb="yellow3", xgb_split="yellow4",
              avg="purple", avgHB="purple2", avgMo="purple3", avgB1="purple4")
mod.i <- tibble(mod_col=names(mod_cols),
                mod_clean=c("Grand mean", "4-week mean", 
                            paste0("GLM-", 1:4),
                            paste0("Logistic-", 1:2),
                            paste0("Ordinal-", 1:2),
                            "RF-1", "RF-2",
                            "SVM-1", "SVM-2",
                            "XGB-1", "XGB-2",
                            paste0("Ensemble-", c("Tot", "Bayes", "Month", "Bloom"))),
                mod_short=c("Mean", "4-wk", 
                            "GLM-1", "GLM-2", "GLM-3", "GLM-4",
                            "L-1", "L-2", 
                            "O-1", "O-2", 
                            "RF-1", "RF-2", 
                            "SVM-1", "SVM-2", 
                            "XGB-1", "XGB-2",
                            "Ens-Avg", "Ens-HB", "Ens-Mo", "Ens-B1"))
modType_cols <- c("Null"="black", "GLM"="grey60", "Ensemble"="#e7298a",
                  "ML: RF"="#d95f02", "ML: SVM"="#e6ab02", "ML: XGB"="#7570b3", 
                  "Bayes: Ordinal"="#1b9e77", "Bayes: Logistic"="#66a61e")


# load output -------------------------------------------------------------

outMn.df <- read_csv(glue("{out.dir}/full_mean_out.csv")) %>%
  mutate(model=factor(model, levels=mod.i$mod_col, labels=mod.i$mod_col)) 





# Pseudo-R2 ---------------------------------------------------------------

p <- outMn.df %>%
  group_by(species, dataSubset, model, modType) %>%
  summarise(LL=sum(dbinom(Nbloom, 1, pred, log=T))) %>%
  group_by(species, dataSubset) %>%
  arrange(species, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  ggplot(aes(model, R2, colour=modType, group=model)) +
  geom_hline(yintercept=0, colour="grey30", size=0.25) +
  geom_point() + geom_rug(sides="l") +
  facet_grid(dataSubset~species) + 
  scale_y_continuous(limits=c(0, 0.6), breaks=c(0, 0.2, 0.4, 0.6)) +
  scale_colour_manual("Model type", values=modType_cols) +
  labs(x="Model", y=expression("McFadden's pseudo-R"^2)) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        panel.grid.minor.y=element_blank(),
        legend.position="bottom")
ggsave(glue("figs/R2_McFadden_{prior_type}.png"), p, width=10, height=5, units="in")

p <- outMn.df %>%
  group_by(species, dataSubset, model, modType) %>%
  summarise(LL=sum(dbinom(Nbloom, 1, pred, log=T))) %>%
  group_by(species, dataSubset) %>%
  arrange(species, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  arrange(desc(R2)) %>% 
  filter(modType != "Ensemble") %>%
  mutate(R2_rank=row_number()) %>%
  mutate(modCat=case_when(modType=="Null" ~ "Null",
                          modType=="GLM" ~ "Reg. GLM",
                          grepl("Bayes", modType) ~ "Bayes",
                          grepl("ML", modType) ~ "ML")) %>%
  filter(!is.na(R2)) %>%
  ggplot(aes(species, R2_rank, fill=modCat)) +
  geom_hline(yintercept=0, colour="grey30", size=0.25) +
  geom_point(shape=22, size=3, colour="grey30") + 
  facet_grid(dataSubset~.) + 
  scale_fill_manual("Model type", values=c("#1b9e77", "#d95f02", "grey30", "#e6ab02")) +
  labs(x="", y=expression("Ranked McFadden's pseudo-R"^2)) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        panel.grid.minor.y=element_blank())
ggsave(glue("figs/R2_ranks_McFadden_{prior_type}.png"), p, width=3, height=6, units="in")

p <- outMn.df %>%
  group_by(species, dataSubset, model, modType) %>%
  summarise(LL=sum(dbinom(Nbloom, 1, pred, log=T)),
            N=n()) %>%
  group_by(species, dataSubset) %>%
  arrange(species, model) %>%
  mutate(R2=(1 - exp((LL - first(LL))/N))/(1 - exp(-first(LL)/N))) %>%
  ggplot(aes(model, R2, colour=modType, group=model)) +
  geom_hline(yintercept=0, colour="grey30", size=0.25) +
  geom_point() + geom_rug(sides="l") +
  facet_grid(dataSubset~species) + 
  scale_y_continuous(limits=c(0, 0.6), breaks=c(0, 0.2, 0.4, 0.6)) +
  scale_colour_manual("Model type", values=modType_cols) +
  labs(x="Model", y=expression("Nagelkerke's pseudo-R"^2)) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        panel.grid.minor.y=element_blank(),
        legend.position="bottom")
ggsave(glue("figs/R2_Nagelkerke_{prior_type}.png"), p, width=10, height=5, units="in")

p <- outMn.df %>%
  group_by(species, dataSubset, model, modType) %>%
  summarise(LL=sum(dbinom(Nbloom, 1, pred, log=T)),
            N=n()) %>%
  group_by(species, dataSubset) %>%
  arrange(species, model) %>%
  mutate(R2=(1 - exp((LL - first(LL))/N))/(1 - exp(-first(LL)/N))) %>%
  arrange(desc(R2)) %>% 
  filter(modType != "Ensemble") %>%
  mutate(R2_rank=row_number()) %>%
  mutate(modCat=case_when(modType=="Null" ~ "Null",
                          modType=="GLM" ~ "Reg. GLM",
                          grepl("Bayes", modType) ~ "Bayes",
                          grepl("ML", modType) ~ "ML")) %>%
  filter(!is.na(R2)) %>%
  ggplot(aes(species, R2_rank, fill=modCat)) +
  geom_hline(yintercept=0, colour="grey30", size=0.25) +
  geom_point(shape=22, size=3, colour="grey30") + 
  facet_grid(dataSubset~.) + 
  scale_fill_manual("Model type", values=c("#1b9e77", "#d95f02", "grey30", "#e6ab02")) +
  labs(x="", y=expression("Ranked Nagelkerke's pseudo-R"^2)) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        panel.grid.minor.y=element_blank())
ggsave(glue("figs/R2_ranks_Nagelkerke_{prior_type}.png"), p, width=3, height=6, units="in")





# TSS ---------------------------------------------------------------------

thresh <- seq(0, 1, by=0.01)
thresh.TSS <- map_dfr(thresh, 
                  ~outMn.df %>%
                    full_join(outMn.df %>% filter(model=="grand") %>%
                                group_by(species) %>% summarise(grandMean=mean(pred))) %>%
                    mutate(thresh=.x,
                           predBloom=as.numeric(pred >= thresh)) %>%
                    group_by(species, model, modType, dataSubset, thresh) %>%
                    summarise(TP=sum(Nbloom==1 & predBloom==1),
                              TN=sum(Nbloom==0 & predBloom==0),
                              FP=sum(Nbloom==0 & predBloom==1),
                              FN=sum(Nbloom==1 & predBloom==0),
                              N=n()) %>%
                    mutate(TPR=TP/(TP+FN),
                           TNR=TN/(FP+TN),
                           TSS=TPR+TNR-1) %>%
                    filter(model!="Grand mean"))
# anim <- ggplot(thresh.TSS, aes(model, TSS, shape=dataSubset, colour=modType)) + 
#   geom_point() + 
#   facet_grid(.~species) +
#   scale_y_continuous(limits=c(-0.1, 1), breaks=c(0, 0.25, 0.5, 0.75, 1)) +
#   scale_colour_manual("Model type", values=modType_cols) +
#   scale_shape_manual(values=c(1,19)) +
#   transition_states(thresh, transition_length=1, state_length=0) +
#   ease_aes('cubic-in-out') +
#   labs(x="Model", y="TSS", title="p(bloom) threshold:{closest_state}") +
#   theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
#         panel.grid.minor.y=element_blank(),
#         legend.position="bottom")
# anim_save("figs/TSS_pBloomThresh.gif", anim, nframes=3*length(thresh), fps=32,
#           width=10, height=4.5, res=300, units="in")


p <- thresh.TSS %>%
  filter(dataSubset=="oos") %>%
  filter(model %in% c("fourWk", "glm", "bernP", "ordP", "rf", "svm", "xgb", "avg")) %>%
  ggplot(aes(thresh, TSS, group=model, colour=model)) +
  geom_vline(data=outMn.df %>% filter(model=="grand") %>%
               group_by(species) %>% summarise(grandMean=mean(pred)),
             aes(xintercept=grandMean), linetype=2) +
  geom_line() +
  scale_colour_manual(values=mod_cols[c("fourWk", "glm", "bernP", "ordP", "rf", "svm", "xgb", "avg")]) +
  facet_grid(.~species) +
  scale_y_continuous(limits=c(-0.1, 1), breaks=c(0, 0.5, 1)) +
  labs(x="p(bloom) threshold", y="TSS", title="Out-of-sample TSS") +
  theme(panel.grid.minor.y=element_blank(),
        legend.position="bottom")
ggsave(glue("figs/TSS_pBloomThresh_range_{prior_type}.png"), p, width=10, height=4.5, dpi=300)

opt.TSS <- thresh.TSS %>%
  filter(dataSubset=="fit") %>%
  group_by(species, model) %>%
  arrange(desc(TSS)) %>%
  slice_head(n=1) %>%
  ungroup %>%
  select(species, model, thresh) %>%
  rename(opt_thresh=thresh)
p <- thresh.TSS %>%
  filter(dataSubset=="oos") %>%
  full_join(., opt.TSS) %>%
  filter(thresh==opt_thresh) %>%
  ggplot(aes(TSS, model, colour=modType)) + 
  geom_point() + 
  scale_colour_manual(values=modType_cols) +
  facet_grid(.~species) +
  scale_x_continuous(limits=c(-0.1, 1), breaks=seq(0,1,0.25)) +
  labs(x="TSS", y="", title="Out-of-sample TSS at fitted optimal p(bloom)") +
  theme(panel.grid.minor.y=element_blank(),
        legend.position="bottom")
ggsave(glue("figs/TSS_pBloomThresh_opt_{prior_type}.png"), p, width=10, height=4.5, dpi=300)






# HSS ---------------------------------------------------------------------

thresh <- seq(0, 1, by=0.01)
thresh.HSS <- map_dfr(thresh, 
                      ~outMn.df %>%
                        full_join(outMn.df %>% filter(model=="grand") %>%
                                    group_by(species) %>% summarise(grandMean=mean(pred))) %>%
                        mutate(thresh=.x,
                               predBloom=as.numeric(pred >= thresh)) %>%
                        group_by(species, model, modType, dataSubset, thresh) %>%
                        summarise(TP=sum(Nbloom==1 & predBloom==1),
                                  TN=sum(Nbloom==0 & predBloom==0),
                                  FP=sum(Nbloom==0 & predBloom==1),
                                  FN=sum(Nbloom==1 & predBloom==0),
                                  N=n()) %>%
                        mutate(expCorrRand=1/N * ((TP+FN)*(TP+FP) + (TN+FN)*(TN+FP)),
                               HSS=((TP+TN)-expCorrRand)/(N-expCorrRand)) %>%
                        filter(model!="Grand mean"))
p <- thresh.HSS %>%
  filter(dataSubset=="oos") %>%
  filter(model %in% c("fourWk", "glm", "bernP", "ordP", "rf", "svm", "xgb", "avg")) %>%
  ggplot(aes(thresh, HSS, group=model, colour=model)) +
  geom_vline(data=outMn.df %>% filter(model=="grand") %>%
               group_by(species) %>% summarise(grandMean=mean(pred)),
             aes(xintercept=grandMean), linetype=2) +
  geom_line() +
  scale_colour_manual(values=mod_cols[c("fourWk", "glm", "bernP", "ordP", "rf", "svm", "xgb", "avg")]) +
  facet_grid(.~species) +
  scale_y_continuous(limits=c(-0.1, 1), breaks=c(0, 0.5, 1)) +
  labs(x="p(bloom) threshold", y="TSS", title="Out-of-sample HSS") +
  theme(panel.grid.minor.y=element_blank(),
        legend.position="bottom")
ggsave(glue("figs/HSS_pBloomThresh_range_{prior_type}.png"), p, width=10, height=4.5, dpi=300)

opt.HSS <- thresh.HSS %>%
  filter(dataSubset=="fit") %>%
  group_by(species, model) %>%
  arrange(desc(HSS)) %>%
  slice_head(n=1) %>%
  ungroup %>%
  select(species, model, thresh) %>%
  rename(opt_thresh=thresh)
p <- thresh.HSS %>%
  filter(dataSubset=="oos") %>%
  full_join(., opt.HSS) %>%
  filter(thresh==opt_thresh) %>%
  ggplot(aes(HSS, model, colour=modType)) + 
  geom_point() + 
  scale_colour_manual(values=modType_cols) +
  facet_grid(.~species) +
  scale_x_continuous(limits=c(-0.1, 1), breaks=seq(0,1,0.25)) +
  labs(x="HSS", y="", title="Out-of-sample HSS at fitted optimal p(bloom)") +
  theme(panel.grid.minor.y=element_blank(),
        legend.position="bottom")
ggsave(glue("figs/HSS_pBloomThresh_opt_{prior_type}.png"), p, width=10, height=4.5, dpi=300)




# ROC ---------------------------------------------------------------------


ROC_cols <- mod_cols[names(mod_cols) %in% unique(outMn.df$model)]
for(j in c("fit", "oos")) {
  binary.ls <- outMn.df %>% filter(dataSubset==j) %>% 
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
    png(glue("figs/performance/ROC_{j}_{sp.i$full[i]}_{prior_type}.png"), width=5, height=5, res=300, units="in")
    plot(NA, NA, xlim=c(0, 1), ylim=c(0, 1),
         xlab="1 - Specificity", ylab="Sensitivity", axes=F, main=sp.i$full[i])
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
    mutate(species=sp.i$full) %>%
    pivot_longer(-species, names_to="model", values_to="AUC") %>%
    mutate(model=factor(model, levels=names(ROC_cols)))  %>%
    mutate(modType=case_when(grepl("fourWk|grand", model) ~ "Null",
                             grepl("glm", model) ~ "GLM",
                             grepl("rf", model) ~ "ML: RF",
                             grepl("^svm", model) ~ "ML: SVM",
                             grepl("^xgb", model) ~ "ML: XGB",
                             grepl("ord", model) ~ "Bayes: Ordinal",
                             grepl("bern", model) ~ "Bayes: Logistic",
                             grepl("avg", model) ~ "Ensemble"))
  write_csv(auc.df, glue("figs/performance/auc_{j}_{prior_type}_all.csv"))
}


# ggplot(auc.df, aes(model, AUC, colour=species, group=species)) + 
#   geom_point() + geom_line() + 
#   # scale_colour_manual(values=sams.cols) + 
#   facet_grid(.~modType, scales="free_x", space="free_x") + 
#   theme(legend.position="bottom")# + ylim(0.5, 1)
# ggplot(auc.df, aes(species, AUC, colour=model, group=model)) + geom_line() + 
#   scale_colour_manual(values=ROC_cols) #+ ylim(0.5, 1)
# par(mfrow=c(1,1))

auc.f <- dir(glue("figs/performance/"), "auc_.*csv", full.names=T)
auc.df <- map_dfr(auc.f, ~read_csv(.x, show_col_types=F) %>%
                    mutate(type=str_split_fixed(.x, "_", 4)[,2],
                           prior=str_split_fixed(.x, "_", 4)[,3])) %>%
  mutate(type=if_else(grepl("fit", type), "fit", "forecast")) %>%
  mutate(modType=case_when(grepl("fourWk|grand", model) ~ "Null",
                           grepl("glm", model) ~ "GLM",
                           grepl("rf", model) ~ "ML: RF",
                           grepl("^svm", model) ~ "ML: SVM",
                           grepl("^xgb", model) ~ "ML: XGB",
                           grepl("ord", model) ~ "Bayes: Ordinal",
                           grepl("bern", model) ~ "Bayes: Logistic",
                           grepl("avg", model) ~ "Ensemble"),
         model=factor(model, levels=mod.i$mod_col, labels=mod.i$mod_short))
p <- ggplot(auc.df, aes(model, AUC, colour=modType, shape=prior)) +
  geom_point() +
  facet_grid(type~species) + 
  scale_colour_manual(values=modType_cols) + 
  scale_shape_manual(values=c("1", "2", "3")) +
  scale_y_continuous(limits=c(0.45, 1), breaks=c(0.5, 0.7, 0.8, 0.9, 1)) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        axis.title.x=element_blank(),
        legend.position="bottom",
        legend.title=element_blank(),
        panel.grid.minor.y=element_blank())
ggsave(glue("figs/AUC_{prior_type}.png"), p, width=8, height=7)



# smooths -----------------------------------------------------------------

sams.cols <- c("alexandrium_sp"="#2A5883",
               "dinophysis_sp"="#46BFDE",
               "karenia_mikimotoi"="#A9DAE0",
               "prorocentrum_lima"="#E77334",
               "pseudo_nitzschia_sp"="#D21E4C")
site_sample <- sample(unique(outMn.df$siteid), 20)
p.effects <- c('tempLwk', 'salinityLwk', 'shortwaveLwk', 'kmLwk', 'precipLwk', 
               'tempStrat20mLwk', 'tempStrat20mRwk', 
               'windVel', 'waterVelL', 'waterVelR', 
               'windLwk', 'waterLwk', 'waterRwk', 
               'fetch', 'influxwk',
               'attnwk', 'chlwk', 'dinowk', 'o2wk', 'phwk', 'po4wk', 
               'tempLwkdt', 'salinityLwkdt', 'shortwaveLwkdt', 
               'precipLwkdt', 'tempStrat20mLwkdt', 'tempStrat20mRwkdt',
               'windLwkdt', 'waterLwkdt', 'waterRwkdt', 
               'chlwkdt', 'o2wkdt', 'phwkdt', 'po4wkdt',
               'NlnWt1', 'NlnWt2', 'NlnRAvg1', 'NlnRAvg2')
if(prior_type=="full") p.effects <- c(p.effects, 'Nbloom1', "Nbloom2")



# * ordP ------------------------------------------------------------------

out.p <- dir(out.dir, "ordP_all_", full.names=T) %>%
  map(readRDS) %>% setNames(sp.i$full)

smooths.df <- bind_cols(expand_grid(yday=seq(0, 364, length.out=52),
                                    siteid=site_sample) %>%
                          mutate(ydayCos=cos(yday*2*pi/365),
                                 ydaySin=sin(yday*2*pi/365),
                                 ydaySC=ydayCos*ydaySin,
                                 date=as_date(yday, origin="2019-01-01")),
                        map_dfr(p.effects, ~tibble(var=.x, val=0)) %>%
                          pivot_wider(names_from="var", values_from="val"))
pred.ls <- vector("list", nrow(sp.i)) %>% setNames(sp.i$full)
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
                                   species=sp.i$full[i])) %>%
    bind_rows(smooths.df %>%
                mutate(slope=colMeans(posterior_epred(out.p[[i]], newdata=.,
                                                      allow_new_levels=T,
                                                      nlpar="bIntercept")),
                       p=1,
                       var="Intercept",
                       species=sp.i$full[i]))
}
pred.ordP <- do.call(rbind, pred.ls) %>%
  rename_all(~str_remove(.x, "wk")) %>%
  rename_all(~str_replace(.x, "L", "_L")) %>%
  rename_all(~str_replace(.x, "R", "_R"))
saveRDS(pred.ordP, glue("out/effects_ordP_{prior_type}.rds"))
p <- ggplot(pred.ordP, aes(date, slope*p, colour=species, group=paste(species, siteid))) + 
  geom_hline(yintercept=0, colour="grey30", linetype=2) +
  geom_line(alpha=0.25) + facet_wrap(~var, ncol=5) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_x_date(date_labels="%b", date_breaks="2 months") +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1))) +
  theme(panel.grid.minor=element_blank(), legend.position="bottom")
ggsave(glue("figs/effects_ordP_{prior_type}.png"), p, width=9, height=9)



# * bernP01 ---------------------------------------------------------------

out.p <- dir(out.dir, "bernP01_all_", full.names=T) %>%
  map(readRDS) %>% setNames(sp.i$full)

smooths.df <- bind_cols(expand_grid(yday=seq(0, 364, length.out=52),
                                    siteid=site_sample) %>%
                          mutate(ydayCos=cos(yday*2*pi/365),
                                 ydaySin=sin(yday*2*pi/365),
                                 ydaySC=ydayCos*ydaySin,
                                 date=as_date(yday, origin="2019-01-01")),
                        map_dfr(p.effects, ~tibble(var=.x, val=0)) %>%
                          pivot_wider(names_from="var", values_from="val"))
pred.ls <- vector("list", nrow(sp.i)) %>% setNames(sp.i$full)
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
                                   species=sp.i$full[i])) %>%
    bind_rows(smooths.df %>%
                mutate(slope=colMeans(posterior_epred(out.p[[i]], newdata=.,
                                                      allow_new_levels=T,
                                                      nlpar="bIntercept")),
                       p=1,
                       var="Intercept",
                       species=sp.i$full[i]))
}
pred.bernP01 <- do.call(rbind, pred.ls) %>%
  rename_all(~str_remove(.x, "wk")) %>%
  rename_all(~str_replace(.x, "L", "_L")) %>%
  rename_all(~str_replace(.x, "R", "_R"))
saveRDS(pred.bernP01, glue("out/effects_bernP01_{prior_type}.rds"))
p <- ggplot(pred.bernP01, aes(date, slope*p, colour=species, group=paste(species, siteid))) + 
  geom_hline(yintercept=0, colour="grey30", linetype=2) +
  geom_line(alpha=0.25) + facet_wrap(~var, ncol=5) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_x_date(date_labels="%b", date_breaks="2 months") +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1))) +
  theme(panel.grid.minor=element_blank(), legend.position="bottom")
ggsave(glue("figs/effects_bernP01_{prior_type}.png"), p, width=9, height=9)



# * bernP11 ---------------------------------------------------------------

out.p <- dir(out.dir, "bernP11_all_", full.names=T) %>%
  map(readRDS) %>% setNames(sp.i$full)

smooths.df <- bind_cols(expand_grid(yday=seq(0, 364, length.out=52),
                                    siteid=site_sample) %>%
                          mutate(ydayCos=cos(yday*2*pi/365),
                                 ydaySin=sin(yday*2*pi/365),
                                 ydaySC=ydayCos*ydaySin,
                                 date=as_date(yday, origin="2019-01-01")),
                        map_dfr(p.effects, ~tibble(var=.x, val=0)) %>%
                          pivot_wider(names_from="var", values_from="val"))
pred.ls <- vector("list", nrow(sp.i)) %>% setNames(sp.i$full)
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
                                   species=sp.i$full[i])) %>%
    bind_rows(smooths.df %>%
                mutate(slope=colMeans(posterior_epred(out.p[[i]], newdata=.,
                                                      allow_new_levels=T,
                                                      nlpar="bIntercept")),
                       p=1,
                       var="Intercept",
                       species=sp.i$full[i]))
}
pred.bernP11 <- do.call(rbind, pred.ls) %>%
  rename_all(~str_remove(.x, "wk")) %>%
  rename_all(~str_replace(.x, "L", "_L")) %>%
  rename_all(~str_replace(.x, "R", "_R"))
saveRDS(pred.bernP11, glue("out/effects_bernP11_{prior_type}.rds"))
p <- ggplot(pred.bernP11, aes(date, slope*p, colour=species, group=paste(species, siteid))) + 
  geom_hline(yintercept=0, colour="grey30", linetype=2) +
  geom_line(alpha=0.25) + facet_wrap(~var, ncol=5) +
  scale_colour_brewer(type="qual", palette="Dark2") +
  scale_x_date(date_labels="%b", date_breaks="2 months") +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1))) +
  theme(panel.grid.minor=element_blank(), legend.position="bottom")
ggsave(glue("figs/effects_bernP11_{prior_type}.png"), p, width=9, height=9)






# * ord ---------------------------------------------------------------

out.p <- dir(out.dir, "ord_all_", full.names=T) %>%
  map(readRDS) %>% setNames(sp.i$full)
p.nonZero <- out.p %>% map(~fixef(.x, probs=c(0.25, 0.75)) %>% 
                             as_tibble(rownames="par") %>%
                             filter(!grepl("Intercept", par) & sign(Q25) == sign(Q75)))

pred.ls <- vector("list", nrow(sp.i)) %>% setNames(sp.i$full)
for(i in seq_along(pred.ls)) {
  pred.ls[[i]] <- map_dfr(p.nonZero[[i]]$par, 
                          ~bind_cols(expand_grid(yday=seq(0, 364, length.out=52),
                                                 siteid=site_sample,
                                                 x_var_temp=seq(-2, 2, 1)) %>%
                                       mutate(x_var=x_var_temp) %>%
                                       setNames(c("yday", "siteid", .x, "x_var")) %>%
                                       mutate(ydayCos=cos(yday*2*pi/365),
                                              ydaySin=sin(yday*2*pi/365),
                                              ydaySC=ydayCos*ydaySin,
                                              date=as_date(yday, origin="2019-01-01")),
                                     map_dfr(grep(glue("^{.x}$"), p.effects, value=T, invert=T),
                                             ~tibble(var=.x, val=0)) %>%
                                       pivot_wider(names_from="var", values_from="val")) %>%
                            mutate(pred=calc_ord_mnpr(posterior_epred(out.p[[i]], newdata=.,
                                                                      allow_new_levels=T),
                                                      bloomThresh[[i]]),
                                   var=.x,
                                   species=sp.i$full[i]))
}
pred.ord <- do.call(rbind, pred.ls) %>%
  rename_all(~str_remove(.x, "wk")) %>%
  rename_all(~str_replace(.x, "L", "_L")) %>%
  rename_all(~str_replace(.x, "R", "_R"))
saveRDS(pred.ord, glue("out/effects_ord_{prior_type}.rds"))
p <- ggplot(pred.ord, aes(date, pred, colour=x_var, group=paste(x_var, siteid))) + 
  geom_hline(yintercept=0, colour="grey30", linetype=2) +
  geom_line(alpha=0.25) + facet_grid(species~var) +
  scale_x_date(date_labels="%b", date_breaks="2 months") +
  scale_colour_viridis_c() +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1))) +
  theme(panel.grid.minor=element_blank(), legend.position="bottom")
ggsave(glue("figs/effects_ord_{prior_type}.png"), p, width=28, height=8)
p <- ggplot(pred.ord, aes(date, x_var, colour=pred)) + 
  geom_jitter(alpha=0.25, width=0, height=0.5) + facet_grid(species~var) +
  scale_x_date(date_labels="%b", date_breaks="2 months") +
  scale_colour_viridis_c() +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1))) +
  theme(panel.grid.minor=element_blank(), legend.position="bottom")
ggsave(glue("figs/effects_pts_ord_{prior_type}.png"), p, width=28, height=8)



# * bern01 ---------------------------------------------------------------

out.p <- dir(out.dir, "bern01_all_", full.names=T) %>%
  map(readRDS) %>% setNames(sp.i$full)
p.nonZero <- out.p %>% map(~fixef(.x, probs=c(0.25, 0.75)) %>% 
                             as_tibble(rownames="par") %>%
                             filter(!grepl("Intercept", par) & sign(Q25) == sign(Q75)))

pred.ls <- vector("list", nrow(sp.i)) %>% setNames(sp.i$full)
for(i in seq_along(pred.ls)) {
  pred.ls[[i]] <- map_dfr(p.nonZero[[i]]$par, 
                          ~bind_cols(expand_grid(yday=seq(0, 364, length.out=52),
                                                 siteid=site_sample,
                                                 x_var_temp=seq(-2, 2, 1)) %>%
                                       mutate(x_var=x_var_temp) %>%
                                       setNames(c("yday", "siteid", .x, "x_var")) %>%
                                       mutate(ydayCos=cos(yday*2*pi/365),
                                              ydaySin=sin(yday*2*pi/365),
                                              ydaySC=ydayCos*ydaySin,
                                              date=as_date(yday, origin="2019-01-01")),
                                     map_dfr(grep(glue("^{.x}$"), p.effects, value=T, invert=T),
                                             ~tibble(var=.x, val=0)) %>%
                                       pivot_wider(names_from="var", values_from="val")) %>%
                            mutate(pred=colMeans(posterior_epred(out.p[[i]], newdata=.,
                                                                 allow_new_levels=T)),
                                   var=.x,
                                   species=sp.i$full[i]))
}
pred.bern01 <- do.call(rbind, pred.ls) %>%
  rename_all(~str_remove(.x, "wk")) %>%
  rename_all(~str_replace(.x, "L", "_L")) %>%
  rename_all(~str_replace(.x, "R", "_R"))
saveRDS(pred.bern01, glue("out/effects_bern01_{prior_type}.rds"))
p <- ggplot(pred.bern01, aes(date, pred, colour=x_var, group=paste(x_var, siteid))) + 
  geom_hline(yintercept=0, colour="grey30", linetype=2) +
  geom_line(alpha=0.25) + facet_grid(species~var) +
  scale_x_date(date_labels="%b", date_breaks="2 months") +
  scale_colour_viridis_c() +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1))) +
  theme(panel.grid.minor=element_blank(), legend.position="bottom")
ggsave(glue("figs/effects_bern01_{prior_type}.png"), p, width=28, height=8)
p <- ggplot(pred.ord, aes(date, x_var, colour=pred)) + 
  geom_jitter(alpha=0.25, width=0, height=0.5) + facet_grid(species~var) +
  scale_x_date(date_labels="%b", date_breaks="2 months") +
  scale_colour_viridis_c() +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1))) +
  theme(panel.grid.minor=element_blank(), legend.position="bottom")
ggsave(glue("figs/effects_pts_bern01_{prior_type}.png"), p, width=28, height=8)




# * bern11 ---------------------------------------------------------------

out.p <- dir(out.dir, "bern11_all_", full.names=T) %>%
  map(readRDS) %>% setNames(sp.i$full)
p.nonZero <- out.p %>% map(~fixef(.x, probs=c(0.25, 0.75)) %>% 
                             as_tibble(rownames="par") %>%
                             filter(!grepl("Intercept", par) & sign(Q25) == sign(Q75)))

pred.ls <- vector("list", nrow(sp.i)) %>% setNames(sp.i$full)
for(i in seq_along(pred.ls)) {
  pred.ls[[i]] <- map_dfr(p.nonZero[[i]]$par, 
                          ~bind_cols(expand_grid(yday=seq(0, 364, length.out=52),
                                                 siteid=site_sample,
                                                 x_var_temp=seq(-2, 2, 1)) %>%
                                       mutate(x_var=x_var_temp) %>%
                                       setNames(c("yday", "siteid", .x, "x_var")) %>%
                                       mutate(ydayCos=cos(yday*2*pi/365),
                                              ydaySin=sin(yday*2*pi/365),
                                              ydaySC=ydayCos*ydaySin,
                                              date=as_date(yday, origin="2019-01-01")),
                                     map_dfr(grep(glue("^{.x}$"), p.effects, value=T, invert=T),
                                             ~tibble(var=.x, val=0)) %>%
                                       pivot_wider(names_from="var", values_from="val")) %>%
                            mutate(pred=colMeans(posterior_epred(out.p[[i]], newdata=.,
                                                                 allow_new_levels=T)),
                                   var=.x,
                                   species=sp.i$full[i]))
}
pred.bern11 <- do.call(rbind, pred.ls) %>%
  rename_all(~str_remove(.x, "wk")) %>%
  rename_all(~str_replace(.x, "L", "_L")) %>%
  rename_all(~str_replace(.x, "R", "_R"))
saveRDS(pred.bern11, glue("out/effects_bern11_{prior_type}.rds"))
p <- ggplot(pred.bern11, aes(date, pred, colour=x_var, group=paste(x_var, siteid))) + 
  geom_hline(yintercept=0, colour="grey30", linetype=2) +
  geom_line(alpha=0.25) + facet_grid(species~var) +
  scale_x_date(date_labels="%b", date_breaks="2 months") +
  scale_colour_viridis_c() +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1))) +
  theme(panel.grid.minor=element_blank(), legend.position="bottom")
ggsave(glue("figs/effects_bern11_{prior_type}.png"), p, width=28, height=8)
p <- ggplot(pred.ord, aes(date, x_var, colour=pred)) + 
  geom_jitter(alpha=0.25, width=0, height=0.5) + facet_grid(species~var) +
  scale_x_date(date_labels="%b", date_breaks="2 months") +
  scale_colour_viridis_c() +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=1))) +
  theme(panel.grid.minor=element_blank(), legend.position="bottom")
ggsave(glue("figs/effects_pts_bern11_{prior_type}.png"), p, width=28, height=8)




