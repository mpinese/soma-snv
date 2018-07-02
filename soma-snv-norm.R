#!/usr/bin/env Rscript
library(plyr)

args = commandArgs(trailingOnly = TRUE)
# Arguments: soma-snv-combine.R <sampleid> <infile_variants> <infile_background> <outfile_burden>

sample = args[[1]]
infile_variants = args[[2]]
infile_background = args[[3]]
outfile_burden = args[[4]]

# Load variants
variants = read.table(infile_variants,
    colClasses = c(character(0), integer(0), character(0), character(0), integer(0), integer(0), integer(0), integer(0), numeric(0), numeric(0)),
    col.names = c("chrom", "pos", "context", "alt", "strand", "dp", "rd", "ad", "sens", "snr"),
    stringsAsFactors = FALSE)
variants$alt = as.character(variants$alt)
variants$alt[variants$alt == "TRUE"] = "T"
if (nrow(variants) > 0) {
    variants$sample = sample
} else {
    variants$sample = character(0)
}

# Load and merge backgrounds
backgrounds = read.table(infile_background,
    colClasses = c(character(0), numeric(0)),
    col.names = c("context", "background"),
    stringsAsFactors = FALSE
)
backgrounds = ddply(backgrounds, .(context), function(d) sum(d$background))
colnames(backgrounds)[2] = "background"
backgrounds$sample = sample
backgrounds = backgrounds[,c("sample", "context", "background")]

# Compute normalised signature burden
signatures = expand.grid(ref_upstream = c("A", "C", "G", "T"), ref = c("C", "T"), ref_downstream = c("A", "C", "G", "T"), alt = c("A", "C", "G", "T"))
signatures = signatures[as.character(signatures$ref) != signatures$alt,]
signatures$signature = sprintf("%s%s %s.%s", signatures$ref, signatures$alt, signatures$ref_upstream, signatures$ref_downstream)
signatures$context = paste(signatures$ref_upstream, signatures$ref, signatures$ref_downstream, sep = "")
signatures$count = unlist(alply(signatures, 1, function(d) sum(variants$context == d$context & variants$alt == d$alt)))
signatures$background = unlist(alply(signatures, 1, function(d) sum(backgrounds$background[backgrounds$context == d$context])))
signatures$burden = signatures$count / signatures$background
signatures$sample = sample
signatures = signatures[,c("sample", "context", "alt", "signature", "count", "background", "burden")]

saveRDS(signatures, outfile_burden)


