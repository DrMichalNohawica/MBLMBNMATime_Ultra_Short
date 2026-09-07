#!/usr/bin/env Rscript
# analysis_update 07: the PROSPERO question answered in its own terms.
# Minimal Bayesian random-effects MA of short(<=4mm)-vs-long MD at one
# timepoint (12m primary; 36m reported with k=3 caveat). Same study set
# and within-study arm pooling as 06; hierarchical normal-normal model,
# weakly informative priors on the mm scale: mu ~ normal(0,1),
# tau ~ half-normal(0,0.5).
#
# VERDICT CRITERIA - fixed before the numbers were seen: report posterior
# median MD [95% CrI], P(MD<0), P(|MD|<0.2), P(|MD|<0.5). Same trap as 04:
# a P(ROPE) from a CrI wider than the ROPE is not evidence of equivalence.
suppressPackageStartupMessages({library(dplyr); library(readr); library(rstan)})
set.seed(1234)
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
pool_arms <- function(g) {
  N <- sum(g$n); M <- weighted.mean(g$mean_mbl, g$n)
  V <- sum((g$n - 1) * g$sd_mbl^2 + g$n * (g$mean_mbl - M)^2) / (N - 1)
  tibble(n = N, mean = M, sd = sqrt(V))
}
code <- "
data { int<lower=1> K; vector[K] md; vector<lower=0>[K] se; }
parameters { real mu; real<lower=0> tau; vector[K] theta; }
model {
  mu ~ normal(0, 1);
  tau ~ normal(0, 0.5);
  theta ~ normal(mu, tau);
  md ~ normal(theta, se);
}"
sm <- stan_model(model_code = code)
res <- list()
for (tp in c(12, 36)) {
  dd <- leg %>%
    filter(Study_ID %in% s4, time_of_followup_from_loading == tp) %>%
    mutate(side = ifelse(implant_length == "short", "short", "long")) %>%
    group_by(Study_ID, side) %>% group_modify(~pool_arms(.x)) %>% ungroup() %>%
    tidyr::pivot_wider(names_from = side, values_from = c(n, mean, sd)) %>%
    filter(!is.na(n_short), !is.na(n_long)) %>%
    mutate(md = mean_short - mean_long,
           se = sqrt(sd_short^2 / n_short + sd_long^2 / n_long))
  fit <- sampling(sm, data = list(K = nrow(dd), md = dd$md, se = dd$se),
                  chains = 4, iter = 4000, warmup = 2000, refresh = 0,
                  control = list(adapt_delta = 0.99))
  mu <- rstan::extract(fit)$mu
  rh <- max(summary(fit)$summary[, "Rhat"], na.rm = TRUE)
  res[[length(res) + 1]] <- data.frame(
    time_m = tp, k = nrow(dd),
    md = median(mu), lo = quantile(mu, .025), hi = quantile(mu, .975),
    p_lt0 = mean(mu < 0),
    p_rope_02 = mean(abs(mu) < 0.2), p_rope_05 = mean(abs(mu) < 0.5),
    tau_med = median(rstan::extract(fit)$tau), max_rhat = round(rh, 4))
}
res <- bind_rows(res) %>% mutate(across(where(is.numeric), ~round(.x, 3)))
write.csv(res, file.path(out, "output/bayes_4mm_v1.csv"), row.names = FALSE)
print(res, row.names = FALSE)
