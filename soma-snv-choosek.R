args = commandArgs(trailingOnly = TRUE)

outpdf = args[[1]]
infiles = args[-1]

suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(SomaticSignatures))

pdf(outpdf, height = 6, width = 6)



cardinality_search = ldply(infiles, function(infile) {
    fit = readRDS(infile)
    k = as.integer(gsub(".*\\.", "", gsub("\\.rds$", "", infile)))
    seed = as.integer(gsub(".*\\.", "", gsub("\\.[0-9]+\\.rds$", "", infile)))
    expvar = evar(object = fitted(fit), target = observed(fit))
    data.frame(infile = infile, k = k, seed = seed, evar = expvar)
})


# Manually add a value for k=1, which is just the row-wise min
temp.data = observed(readRDS(infiles[[1]]))
temp.k1.sig = apply(temp.data, 1, min)
if (max(temp.k1.sig) == 0) {
    temp.k1.evar = 0
} else {
    temp.k1.score = temp.data / temp.k1.sig
    temp.k1.score[is.nan(temp.k1.score)] = Inf
    temp.k1.score = apply(temp.k1.score, 2, min)
    temp.k1.evar = evar(tcrossprod(temp.k1.sig, temp.k1.score), temp.data)
}
cardinality_search = rbind(data.frame(infile = NA, k = 1, seed = NA, evar = temp.k1.evar), cardinality_search)


cardinality_search_best = ddply(cardinality_search, .(k), function(d) d[which.max(d$evar),])

ggplot(cardinality_search, aes(x = k, y = evar)) + geom_jitter(width = 0.1, height = 0) + geom_line(data = cardinality_search_best) + ylim(0, 1) + xlab("Cardinality (k)") + ylab("Explained variance") + ggtitle("Cardinality search")


for (i in 2:nrow(cardinality_search_best))
{
    fit = readRDS(as.character(cardinality_search_best$infile[i]))
    print(plotSignatures(fit) + ggtitle(sprintf("k=%d seed=%d evar=%f", cardinality_search_best$k[i], cardinality_search_best$seed[i], cardinality_search_best$evar[i])))
}


dev.off()
