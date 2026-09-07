#!/usr/bin/env Rscript
# analysis_update 03: the Stan longitudinal meta-regression, refit on v1.
# Same model file as submitted (multilevel_model.stan), new data; site gains
# a third level ('mixed') where trials never split jaws - the factor absorbs it.
suppressPackageStartupMessages({library(dplyr); library(readr); library(rstan); library(tidyr)})
options(mc.cores = parallel::detectCores()); set.seed(1234)
out <- "."
d <- read_csv(list.files("data", pattern="^MBL_data_extracted",
                         full.names=TRUE)[1], show_col_types=FALSE) %>%
  separate(group, c("length_class","graft_class","site_class"), sep="_", remove=FALSE) %>%
  filter(sd_mbl > 0, n > 0) %>%
  mutate(study_factor = factor(Study_ID),
         site_factor = factor(site_class), length_factor = factor(length_class),
         graft_factor = factor(graft_class),
         se_mbl = sd_mbl/sqrt(n),
         time_cen_yr = (time_of_followup_from_loading - 12)/12,
         time_sqrt = sqrt(time_of_followup_from_loading))
sd_ <- list(N=nrow(d), y=d$mean_mbl, se=d$se_mbl,
            time_cen_yr=d$time_cen_yr, time_sqrt=d$time_sqrt,
            n_site=nlevels(d$site_factor), n_length=nlevels(d$length_factor),
            n_graft=nlevels(d$graft_factor), n_study=nlevels(d$study_factor),
            site_id=as.integer(d$site_factor), length_id=as.integer(d$length_factor),
            graft_id=as.integer(d$graft_factor), study_id=as.integer(d$study_factor))
m <- stan_model("model/multilevel_model.stan")
fit <- sampling(m, data=sd_, chains=4, iter=4000, warmup=2000, seed=1234,
                control=list(adapt_delta=0.99, max_treedepth=15), refresh=500)
saveRDS(fit, file.path(out,"fits/stan_v1.rds"))
s <- summary(fit)$summary
write.csv(s, file.path(out,"output/stan_v1_summary.csv"))
cat(sprintf("Stan v1 done: N=%d rows, %d studies, max Rhat=%.4f, min n_eff=%.0f\n",
            sd_$N, sd_$n_study, max(s[,"Rhat"], na.rm=TRUE), min(s[,"n_eff"], na.rm=TRUE)))
