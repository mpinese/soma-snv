args = commandArgs(trailingOnly = TRUE)

targetk = as.integer(args[[1]])
outfile = args[[2]]
infiles = args[c(-1,-2)]

suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(SomaticSignatures))


all_fits_for_k = llply(infiles, function(infile) {
    k = as.integer(gsub(".*\\.", "", gsub("\\.rds$", "", infile)))
    if (k != targetk)
        return(list(fit = NULL, evar = -1))
    fit = readRDS(infile)
    seed = as.integer(gsub(".*\\.", "", gsub("\\.[0-9]+\\.rds$", "", infile)))
    expvar = evar(object = fitted(fit), target = observed(fit))
    list(fit = fit, evar = expvar)
})

best_fit_for_k = all_fits_for_k[[which.max(sapply(all_fits_for_k, function(item) item$evar))]]$fit

saveRDS(best_fit_for_k, file = outfile)
