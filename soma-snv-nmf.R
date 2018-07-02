#!/usr/bin/env Rscript
library(doParallel)
suppressPackageStartupMessages(library(SomaticSignatures))

params = commandArgs(trailingOnly = TRUE)

burden = params[[1]]
outprefix = params[[2]]
seed = params[[3]]
kmin = params[[4]]
kmax = params[[5]]
B = params[[6]]
cores = params[[7]]

zeroes = apply(burden == 0, 2, all)
burden = burden[,!zeroes]

# Factorise
for (i in kmin:kmax) {
    outfile = sprintf("%s.%d.rds", outprefix, i)
    message(sprintf("Fitting cardinality %d", i))
    fit = identifySignatures(burden, i, nmfDecomposition, nrun = B, includeFit = TRUE, seed = seed, .options = sprintf("p%d", cores))
    saveRDS(fit, outfile)
}
