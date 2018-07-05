#!/usr/bin/env Rscript

# Define the variant signatures
signatures = expand.grid(ref_upstream = c("A", "C", "G", "T"), ref = c("C", "T"), ref_downstream = c("A", "C", "G", "T"), alt = c("A", "C", "G", "T"))
signatures = signatures[as.character(signatures$ref) != signatures$alt,]
signatures = signatures[order(signatures$ref, signatures$alt, signatures$ref_upstream, signatures$ref_downstream),]
signatures$signature = sprintf("%s%s %s.%s", signatures$ref, signatures$alt, signatures$ref_upstream, signatures$ref_downstream)
stopifnot(nrow(signatures) == 96)

# Get the output file
outfile = commandArgs(trailingOnly = TRUE)[[1]]

# Get the input files
infiles = commandArgs(trailingOnly = TRUE)[-1]
if (any(duplicated(infiles)))
    stop("Error: duplicate input files.  Include each input file only once on the command line.")

# Read the input files to create a file -> sampleID mapping
samples = sapply(infiles, function(infile) readRDS(infile)$sample[[1]])

if (any(duplicated(samples)))
    stop("Error: some input files share sample IDs.  Sample IDs must be unique.")

# Create signature x sample count and background matrices to populate
count = matrix(NA, nrow = nrow(signatures), ncol = length(infiles), dimnames = list(signature = signatures$signature, sample = samples))
background = matrix(NA, nrow = nrow(signatures), ncol = length(infiles), dimnames = list(signature = signatures$signature, sample = samples))

# Populate matrices
message("Processing:")
for (i in 1:length(infiles))
{
    message(sprintf("  %s", samples[[i]]))
    i.data = readRDS(infiles[[i]])
    stopifnot(i.data$sample == samples[[i]])
    stopifnot(i.data$sample == colnames(count)[i])
    stopifnot(i.data$sample == colnames(background)[i])
    count[,i] = i.data$count[match(rownames(count), i.data$signature)]
    background[,i] = i.data$background[match(rownames(background), i.data$signature)]
}

saveRDS(count / background * 1e6, outfile)
message("Done")
