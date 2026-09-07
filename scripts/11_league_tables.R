#!/usr/bin/env Rscript
# analysis_update 11: full enumeration of pairwise comparisons - the league
# tables (PRISMA-NMA standard). Every node vs every node from the MBNMAtime
# posterior at 12/36/60 months. Fig 4 shows the pre-specified primary
# contrasts; these tables are the complete matrix it summarises.
suppressPackageStartupMessages({library(dplyr)})
out <- "."
fit <- readRDS(file.path(out, "fits/mbnmatime_v1.rds"))
sims <- fit$BUGSoutput$sims.list
trt <- fit$network$treatments
lab <- function(x) {
  x <- gsub("^short|_short", "ultra-short", x)
  x <- sub("^long", "longer", x)
  gsub("_", " ", x)
}
for (t in c(12, 36, 60)) {
  m <- sims$d.1 * sqrt(t) + sims$d.2 * t
  colnames(m) <- trt
  k <- length(trt)
  cell <- matrix("", k, k, dimnames = list(lab(trt), lab(trt)))
  for (i in seq_len(k)) for (j in seq_len(k)) {
    if (i == j) { cell[i, j] <- "-"; next }
    dd <- m[, i] - m[, j]
    cell[i, j] <- sprintf("%+.2f [%+.2f; %+.2f]", median(dd),
                          quantile(dd, .025), quantile(dd, .975))
  }
  write.csv(cell, file.path(out, sprintf("tables/league_mbnmatime_%dm.csv", t)))
}
cat("League tables written: row minus column, MD in mm, 12/36/60m\n")
