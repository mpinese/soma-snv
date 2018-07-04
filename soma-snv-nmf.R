#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(SomaticSignatures))

params = commandArgs(trailingOnly = TRUE)

burden = readRDS(params[[1]])
outprefix = params[[2]]
seed = as.integer(params[[3]])
kmin = as.integer(params[[4]])
kmax = as.integer(params[[5]])
B = as.integer(params[[6]])
cores = as.integer(params[[7]])

# Remove samples without any somatic muts at all
zero_samples = apply(burden == 0, 2, all)
burden = burden[,!zero_samples]

# Check for completely missing variant classes
zero_motifs = names(which(apply(burden == 0, 1, all)))
if (length(zero_motifs) > 0)
{
    constant = min(burden[burden != 0]) / 2
    if (length(zero_motifs) > 10)
        warning(sprintf("A very large number of motifs (%d: %s) have no somatic mutations.  Setting zero cells to %g as workaround.  Fits are likely to be unreliable; exercise great care with results.", length(zero_motifs), paste(zero_motifs, collapse = ","), constant))
    else
        warning(sprintf("Some motifs (%d: %s) have no somatic mutations.  Setting zero cells to %g as workaround.  Fits may be unreliable; exercise care with results.", length(zero_motifs), paste(zero_motifs, collapse = ","), constant))

    # Add a very small constant mut rate to address these zeroes.
    burden[burden == 0] = constant
}

# Factorise
for (i in kmin:kmax) {
    outfile = sprintf("%s.%d.rds", outprefix, i)
    message(sprintf("Fitting cardinality %d", i))
    fit = identifySignatures(burden, i, nmfDecomposition, nrun = B, includeFit = TRUE, seed = seed, .options = sprintf("p%d", cores))
    saveRDS(fit, outfile)
}
