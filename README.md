# soma-snv
Detect and quantify somatic variants in low-depth sequencing data


## Overview

Somatic variants are difficult to detect using error-prone massively-parallel sequencing technology.  This problem is particularly challenging for low-depth data, for which low-VAF (variant allele frequency) variants and machine errors are difficult to distinguish.  soma-snv is a tool to address this problem, by giving estimates of somatic variant burden using 30X whole-genome sequencing data.

Conceptually, soma-snv combines a permissive variant calling front-end with a spectral denoising stage.  The front-end stage tracks expected somatic variant detection senstivity per sample, enabling the normalisation of somatic burden by each sample's detection cross-section.


## Requirements

* samtools (tested with v1.5)
* Python 3 (tested with pypy 5.8)
  * intervaltree library https://pypi.python.org/pypi/intervaltree (tested with 2.1.0)
* R (tested with 3.4.2)
  * SomaticSignatures library http://bioconductor.org/packages/release/bioc/html/SomaticSignatures.html (tested with 2.14.0)
  * doParallel library (tested with 1.0.11)

Note that the use of pypy is strongly recommended for speed.


## Workflow

soma-snv currently only works on datasets containing many samples, although a single sample mode is planned.  The workflow proceeds in two stages:
1. Per-sample variant detection
2. Per-sample variant normalisation
3. Cohort merging
4. Cohort spectral denoising and burden estimation


## Example

### 1. Per-sample: variant detection and cross-section estimation
The following code detects putative somatic variants and estimates the background cross-section for each sample:
```
samtools mpileup -d 101 -a -C 50 -E -I -q 30 -Q 30 -s -r "${region}" -f "${REFERENCE}" "${INFILE}" | \
pypy soma-snv.py --blacklist "${BLACKLIST}" --background "${outfile_background}" --snr 5.0 --vaf 0.2 --error 2.0e-3 --het 1.0e-4 --scv 5.0e-7 --maxdp 100 > "${outfile_variants}"
```

${INFILE} is the input BAM
${BLACKLIST} is the path to a bed file of blacklist regions.  The included data/blacklist_truseq_nano_hiseqX_hs37d5x.bed is appropriate for Illumina TruSeq Nano libraries sequenced on the HiSeq X system and mapped to hs37d5 with additional PhiX decoy.
${region} is an optional genomic region to analyse, for parallelisation purposes

### 1b. Per-sample: combine sharded results from step 1
If parallelisation was used in step 1 (via the ${region} parameter), the individual shards must be combined into two files for each sample:

```
cat ${sampleid}.shard*.variants.tsv > ${sampleid}.variants.tsv
cat ${sampleid}.shard*.background.tsv > ${sampleid}.background.tsv
```


### 2. Per-sample: compute normalised variant burden
For each sample the following script is run to prepare for spectral denoising:
```
Rscript soma-snv-norm.R "${sampleid}" "${infile_variants}" "${infile_background}" "${outfile_burden}"
```

"${sampleid}" is the sample identifier
"${infile_variants}" is the path to a variant file produced in step 1.
"${infile_background}" is the path to the matching background file produced in step 1.
"${outfile_burden}" is a path to an output rds.


### 3. Cohort merging
```
Rscript soma-snv-merge.R "${outfile_merged}" infile_burden_1.rds infile_burden_2.rds ...
```

"${outfile_merged}" is the merged cohort data file.
infile_burden_i.rds is a normalised variant burden file from step 2.


### 4. NMF decomposition
```
Rscript soma-snv-nmf.R "${outfile_merged}" "${prefix}" <seed> <kmin> <kmax> <B> <cores>
```

## Generation of a custom blacklist bed

Section under construction.


## Future

Planned extensions to soma-snv:
* A single-sample mode that uses pre-trained spectral decompositions to denoise results.
* A third workflow stage that annotates variants identified in stage 1 with likelihoods of membership in each of the spectral components identified in stage 2.
