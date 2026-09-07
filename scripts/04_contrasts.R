#!/usr/bin/env Rscript
# analysis_update 04: derived short-vs-long contrasts + equivalence.
#
# VERDICT CRITERIA - fixed before the numbers were seen:
#   ROPE margins: |diff| < 0.2 mm (strict, ~measurement error of periapical
#   radiographs) and |diff| < 0.5 mm (the between-implant-system difference
#   the field treats as clinically negligible). Report P(in ROPE), P(diff<0)
#   at 12/36/60 months from loading. The trap named in advance: a contrast
#   can sit "in ROPE" simply because its CI is enormous at extrapolated
#   times - so the CI width is printed next to every probability; a P(ROPE)
#   claimed from a CI wider than the ROPE itself is not evidence of
#   equivalence and must not be reported as such.
#
# MBNMAtime: contrast = short_X - long_X within stratum (ungrafted
# maxilla/mandible/mixed), from posterior draws of d.1 (sqrt-t) and d.2
# (linear t) vs the network reference. Stan: contrast = length effect +
# time interactions, from the refitted submitted regression.
suppressPackageStartupMessages({library(dplyr)})
out <- "."
times <- c(12, 36, 60)

fit <- readRDS(file.path(out, "fits/mbnmatime_v1.rds"))
sims <- fit$BUGSoutput$sims.list
trt <- fit$network$treatments
idx <- function(nm) which(trt == nm)
pairs <- list(
  c("short_ungrafted_maxilla",  "long_ungrafted_maxilla"),
  c("short_ungrafted_mandible", "long_ungrafted_mandible"),
  c("short_ungrafted_mixed",    "long_ungrafted_mixed"))
rows <- list()
for (p in pairs) {
  i <- idx(p[1]); j <- idx(p[2])
  d1 <- sims$d.1[, i] - sims$d.1[, j]
  d2 <- sims$d.2[, i] - sims$d.2[, j]
  for (t in times) {
    diff <- d1 * sqrt(t) + d2 * t
    rows[[length(rows) + 1]] <- data.frame(
      engine = "MBNMAtime", stratum = sub("short_", "", p[1]), time_m = t,
      median = median(diff),
      lo = quantile(diff, .025), hi = quantile(diff, .975),
      ci_width = diff(quantile(diff, c(.025, .975))),
      p_lt0 = mean(diff < 0),
      p_rope_02 = mean(abs(diff) < 0.2), p_rope_05 = mean(abs(diff) < 0.5))
  }
}

# Stan refit: short - long trajectory from the scaled effects the model
# exposes. Factor levels alphabetical in 03 (long=1, short=2); time enters
# as sqrt(t) and time_cen_yr = (t-12)/12, exactly as 03 built the design.
sfit <- readRDS(file.path(out, "fits/stan_v1.rds"))
ex <- rstan::extract(sfit)
for (t in times) {
  diff <- (ex$length_effect[, 2] - ex$length_effect[, 1]) +
    (ex$length_time_sqrt[, 2] - ex$length_time_sqrt[, 1]) * sqrt(t) +
    (ex$length_time_cen[, 2] - ex$length_time_cen[, 1]) * (t - 12) / 12
  rows[[length(rows) + 1]] <- data.frame(
    engine = "Stan regression", stratum = "all (length effect)", time_m = t,
    median = median(diff),
    lo = quantile(diff, .025), hi = quantile(diff, .975),
    ci_width = diff(quantile(diff, c(.025, .975))),
    p_lt0 = mean(diff < 0),
    p_rope_02 = mean(abs(diff) < 0.2), p_rope_05 = mean(abs(diff) < 0.5))
}
res <- bind_rows(rows) %>% mutate(across(where(is.numeric), ~round(.x, 3)))
write.csv(res, file.path(out, "output/contrasts_short_vs_long_v1.csv"),
          row.names = FALSE)
print(res, row.names = FALSE)
