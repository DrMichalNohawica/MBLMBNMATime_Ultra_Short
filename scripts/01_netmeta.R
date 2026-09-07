#!/usr/bin/env Rscript
# analysis_update/scripts/01_netmeta.R - frequentist NMA per timepoint on the
# v1 extraction. First numbers + the STUDY-LEVEL funnel/Egger (R2.9 answer).
suppressPackageStartupMessages({library(dplyr); library(readr); library(netmeta); library(meta)})
d <- read_csv(list.files("data",
                         pattern="^MBL_data_extracted", full.names=TRUE)[1],
              show_col_types = FALSE)
dir.create("output", showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)
for (tp in c(0, 12, 36, 60)) {
  dt <- d %>% filter(time_of_followup_from_loading == tp,
                     sd_mbl > 0, n > 0) %>%
    group_by(Study_ID, group) %>%
    summarise(mean_mbl = weighted.mean(mean_mbl, n),
              sd_mbl = sqrt(weighted.mean(sd_mbl^2, n)),
              n = sum(n), .groups="drop")
  multi <- dt %>% count(Study_ID) %>% filter(n >= 2) %>% pull(Study_ID)
  dt <- dt %>% filter(Study_ID %in% multi)
  if (n_distinct(dt$group) < 2 || nrow(dt) < 4) { cat(tp, "m: too sparse\n"); next }
  pw <- pairwise(treat = group, n = n, mean = mean_mbl, sd = sd_mbl,
                 studlab = Study_ID, data = dt, sm = "MD")
  nt <- tryCatch(netmeta(pw, reference.group = grep("^long", unique(dt$group), value=TRUE)[1],
                         common = FALSE), error = function(e) {cat(tp,"m netmeta error:",conditionMessage(e),"\n"); NULL})
  if (is.null(nt)) next
  lg <- netleague(nt, digits = 2)$random
  write.csv(lg, sprintf("output/netmeta_league_%dm.csv", tp))
  ps <- netrank(nt, small.values = "desirable")$ranking.random
  write.csv(data.frame(group = names(ps), pscore = ps),
            sprintf("output/netmeta_pscore_%dm.csv", tp), row.names = FALSE)
  cat(sprintf("%dm: %d studies, %d comparisons, %d nodes | tau=%.3f\n",
              tp, nt$k, nrow(pw), nt$n, nt$tau))
  # headline short-vs-long direct MDs
  sl <- pw %>% filter(grepl("^short", treat1) & grepl("^long", treat2) |
                      grepl("^long", treat1) & grepl("^short", treat2))
  if (nrow(sl) > 0) {
    mg <- metagen(TE, seTE, data = sl, sm = "MD", common = FALSE)
    cat(sprintf("   pooled DIRECT short-vs-long MD = %+.3f mm [%.3f; %.3f], I2=%.0f%%, k=%d\n",
                mg$TE.random, mg$lower.random, mg$upper.random, 100*mg$I2, mg$k))
    if (tp == 12) {
      png("output/funnel_study_level_12m.png", 900, 700, res = 120)
      funnel(mg, studlab = TRUE, cex.studlab = .5)
      title("Study-level funnel - direct short vs long MD, 12 months")
      dev.off()
      eg <- tryCatch(metabias(mg, method.bias = "Egger", k.min = 5), error=function(e) NULL)
      if (!is.null(eg)) cat(sprintf("   Egger (study-level, 12m): intercept p = %.3f\n", eg$p.value))
    }
  }
}
