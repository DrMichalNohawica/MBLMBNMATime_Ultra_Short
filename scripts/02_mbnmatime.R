#!/usr/bin/env Rscript
# analysis_update 02: MBNMAtime on v1 - the contrast-based PRIMARY engine.
# Mirrors the submitted crosscheck's calls, pointed at the new extraction.
suppressPackageStartupMessages({library(dplyr); library(readr); library(MBNMAtime)})
set.seed(1234)
out <- "."
d <- read_csv(list.files("data", pattern="^MBL_data_extracted",
                         full.names=TRUE)[1], show_col_types=FALSE) %>%
  transmute(studyID = as.integer(factor(Study_ID)), study_label = Study_ID,
            time = as.numeric(time_of_followup_from_loading),
            treatment = group, y = mean_mbl, se = sd_mbl/sqrt(n), n = n) %>%
  filter(se > 0) %>% arrange(studyID, treatment, time)
# Amato 2020 (W3024572837) is the sole study with a grafted-short arm: one
# contrast at one timepoint cannot identify the node's two relative
# time-course parameters (d.1/d.2 sampled from prior, Rhat ~3.9). Excluded
# from this engine only; it stays in netmeta and the regression.
d <- filter(d, study_label != "W3024572837")
# The time-course model needs a within-study contrast: >=2 arms sharing a
# timepoint. Retrospective cohorts reporting only arm-specific mean follow-up
# (e.g. Schiegnitz 2022: 42.1mo TG vs 32.9mo CG) have none - they stay in
# netmeta bins and the continuous-time regression, but must be dropped here.
shared <- d %>% group_by(studyID, time) %>%
  summarise(k = n_distinct(treatment), .groups = "drop") %>%
  group_by(studyID) %>% summarise(any_shared = any(k >= 2), .groups = "drop")
dropped <- d %>% semi_join(filter(shared, !any_shared), by = "studyID") %>%
  distinct(study_label)
if (nrow(dropped)) cat("Dropped (no shared timepoint across arms):",
                       paste(dropped$study_label, collapse = ", "), "\n")
d <- d %>% semi_join(filter(shared, any_shared), by = "studyID") %>%
  mutate(studyID = as.integer(factor(studyID)))
write_csv(distinct(d, studyID, study_label), file.path(out,"output/mbnmatime_study_ids.csv"))
ref <- "long_grafted_mandible"
d$treatment <- factor(d$treatment, levels = c(ref, sort(setdiff(unique(d$treatment), ref))))
net <- mb.network(select(d, -study_label), reference = ref)
sink(file.path(out,"output/mbnmatime_network_summary.txt")); print(summary(net)); sink()
fp <- tfpoly(degree = 2,
             pool.1 = "rel", method.1 = "random", method.power1 = 0.5,
             pool.2 = "rel", method.2 = "random", method.power2 = 1)
fit <- mb.run(net, fun = fp, n.chain = 4, n.iter = 4000, n.burnin = 2000)
saveRDS(fit, file.path(out,"fits/mbnmatime_v1.rds"))
# summary(fit) errors on FP models in this MBNMAtime version - use BUGSoutput
conv <- round(fit$BUGSoutput$summary[, c("mean","2.5%","97.5%","Rhat","n.eff")], 3)
write.csv(conv, file.path(out,"output/mbnmatime_convergence_v1.csv"))
bad <- conv[!is.na(conv[,"Rhat"]) & conv[,"Rhat"] > 1.05, , drop = FALSE]
cat(sprintf("MBNMAtime: %d params, %d with Rhat>1.05\n", nrow(conv), nrow(bad)))
pred <- predict(fit, times = c(0, 12, 36, 60))
pdf(file.path(out,"output/mbnmatime_predictions_v1.pdf"), width = 10, height = 7)
plot(pred); dev.off()
sink(file.path(out,"output/mbnmatime_predictions.txt")); print(pred); sink()
cat("MBNMAtime v1 done\n")
