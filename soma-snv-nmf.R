#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(SomaticSignatures))

params = commandArgs(trailingOnly = TRUE)

burden = readRDS(params[[1]])
outprefix = params[[2]]
seed = as.integer(params[[3]])
kmin = as.integer(params[[4]])
kmax = as.integer(params[[5]])
B = 1
cores = 1

# Remove samples without any somatic muts at all
zero_samples = apply(burden == 0, 2, all)
burden = burden[,!zero_samples]
message(sprintf("Dropped %d samples with no detected mutations", sum(zero_samples)))

if (ncol(burden) < 30)
    warning(sprintf("Very few (%d) samples: results are likely to be unreliable.", ncol(burden)))

# Check for completely missing variant classes
zero_motifs = names(which(apply(burden == 0, 1, all)))
if (length(zero_motifs) > 0)
{
    constant = min(burden[burden != 0]) / 2
    if (length(zero_motifs) > 10)
        warning(sprintf("A very large number of signatures (%d: %s) have no mutations.  Setting zero cells to runif(%g,%g) as workaround.  Fits will be uninformative for these signatures and are likely to be unreliable; exercise great care with results.", length(zero_motifs), paste(zero_motifs, collapse = ","), constant/10, constant))
    else
        warning(sprintf("Some signatures (%d: %s) have no mutations.  Setting zero cells to runif(%g,%g) as workaround.  Fits will be uninformative for these signatures; exercise care with results.", length(zero_motifs), paste(zero_motifs, collapse = ","), constant/10, constant))

    # Add a very small mut rate to address these zeroes.
    for (zero_motif_i in zero_motifs)
        burden[zero_motif_i,] = runif(ncol(burden), min = constant/10, max = constant)
}

# Factorise
for (i in kmin:kmax) {
    outfile = sprintf("%s.%d.%d.rds", outprefix, seed, i)
    message(sprintf("Fitting cardinality %d", i))
    fit = identifySignatures(burden, i, nmfDecomposition, nrun = B, includeFit = TRUE, seed = seed, .options = sprintf("p%d", cores))
    saveRDS(fit, outfile)
}

message("Done")
