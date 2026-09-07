#!/usr/bin/env Rscript
# analysis_update 06: PROSPERO-registered length criterion as a sensitivity.
# The registration named ~4 mm "extra-short" implants; the review as run
# uses <=6 mm. Here the short arm is restricted to studies whose short
# implants are <=4 mm (nominal/endosseous), and short-vs-long is pooled
# per timepoint as a pairwise random-effects MA (REML + Hartung-Knapp -
# k is small, so exact per-study contrasts, no network machinery).
#
# VERDICT CRITERIA - fixed before the numbers were seen: the sensitivity
# is "consistent with the primary" if the pooled MD keeps the primary's
# direction (short not worse) and its 95% CI overlaps the primary NMA's
# 12m/36m short-vs-long interval. The trap named in advance: with ~6
# studies the CI will be wide - a wide straddling-zero CI is "no evidence
# of difference", not "evidence of no difference"; only the direction and
# overlap claims may be made.
suppressPackageStartupMessages({library(dplyr); library(readr); library(meta)})
out <- "."
leg <- read_csv(list.files("data",
                           pattern = "^MBL_data_extracted", full.names = TRUE)[1],
                show_col_types = FALSE)
ext <- read_csv(list.files("data",
                           pattern = "^extended_regression", full.names = TRUE)[1],
                show_col_types = FALSE)
s4 <- ext %>%
  filter(arm == "short",
         pmin(suppressWarnings(as.numeric(nominal_length_mm)),
              suppressWarnings(as.numeric(endosseous_length_mm)),
              na.rm = TRUE) <= 4.05) %>%
  distinct(study_id) %>% pull(study_id)
cat("<=4mm studies:", paste(s4, collapse = ", "), "\n")

# collapse each study to one short and one long cell per timepoint
# (multi-pair trials, e.g. per-jaw arms, pooled with combined-group formulas)
pool_arms <- function(g) {
  N <- sum(g$n); M <- weighted.mean(g$mean_mbl, g$n)
  V <- sum((g$n - 1) * g$sd_mbl^2 + g$n * (g$mean_mbl - M)^2) / (N - 1)
  tibble(n = N, mean = M, sd = sqrt(V))
}
res <- list()
for (tp in c(12, 36)) {
  dd <- leg %>%
    filter(Study_ID %in% s4, time_of_followup_from_loading == tp) %>%
    mutate(side = ifelse(implant_length == "short", "short", "long")) %>%
    group_by(Study_ID, side) %>% group_modify(~pool_arms(.x)) %>% ungroup() %>%
    tidyr::pivot_wider(names_from = side, values_from = c(n, mean, sd)) %>%
    filter(!is.na(n_short), !is.na(n_long))
  if (nrow(dd) < 2) { cat(sprintf("%dm: only %d studies, skipped\n", tp, nrow(dd))); next }
  ma <- metacont(n_short, mean_short, sd_short, n_long, mean_long, sd_long,
                 data = dd, studlab = Study_ID, sm = "MD",
                 random = TRUE, common = FALSE,
                 method.tau = "REML", method.random.ci = "HK")
  cat(sprintf("\n== %dm, k=%d ==\n", tp, nrow(dd)))
  print(summary(ma))
  res[[length(res) + 1]] <- data.frame(
    time_m = tp, k = nrow(dd),
    md = ma$TE.random, lo = ma$lower.random, hi = ma$upper.random,
    i2 = ma$I2, tau = sqrt(ma$tau2),
    pred_lo = ma$lower.predict, pred_hi = ma$upper.predict)
  png(file.path(out, sprintf("figures/sensitivity_4mm_forest_%dm.png", tp)),
      width = 2400, height = 240 + 160 * nrow(dd), res = 300)
  forest(ma, leftlabs = c("Study", "n", "mean", "SD", "n", "mean", "SD"),
         label.left = "favours short", label.right = "favours long")
  dev.off()
}
res <- bind_rows(res) %>% mutate(across(where(is.numeric), ~round(.x, 3)))
write.csv(res, file.path(out, "output/sensitivity_4mm_v1.csv"), row.names = FALSE)
print(res, row.names = FALSE)
