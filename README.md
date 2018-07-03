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

soma-snv currently only works on datasets containing many samples, although a single sample mode is planned.  The workflow proceeds as follows:
1. Per-sample variant detection
2. Per-sample variant normalisation
3. Cohort merging
4. Cohort spectral denoising and burden estimation


## Example

### 1. Per-sample: variant detection and cross-section estimation
The following code detects putative somatic variants and estimates the background cross-section for each sample:
```
samtools mpileup \
  -d 101 -a -C 50 -E -I -q 30 -Q 30 -s -f reference.fa sampleID.bam | \
pypy soma-snv.py \
  --blacklist data/blacklist_truseq_nano_hiseqX_hs37d5x.bed --snr 5.0 --vaf 0.2 \
  --error 2.0e-3 --het 1.0e-4 --scv 5.0e-7 --maxdp 100 \
  --background sampleID.background.tsv \
> sampleID.variants.tsv
```

`sampleID.bam` is the input BAM
`data/blacklist_truseq_nano_hiseqX_hs37d5x.bed` is the path to a bed file of blacklist regions.  The included data/blacklist_truseq_nano_hiseqX_hs37d5x.bed is appropriate for Illumina TruSeq Nano libraries sequenced on the HiSeq X system and mapped to hs37d5 with additional PhiX decoy.
`sampleID.background.tsv` and `sampleID.variants.tsv` are the output files from this stage.

This stage may be parallelised for speed, in which case add a `-r <region>` flag to each `samtools mpileup` command to define a subset of the genome to examine for somatic variation.  The output IDs must be changed to match (eg `sampleID.variants.tsv` might be changed to `sampleID.shard001.variants.tsv`, likewise for the background file).  These subsets are then combined in step 1b.

### 1b. Per-sample: combine sharded results from step 1
If parallelisation was used in step 1, the individual shards must be combined into two files for each sample:

```
cat sampleID.shard*.variants.tsv > sampleID.variants.tsv
cat sampleID.shard*.background.tsv > sampleID.background.tsv
```


### 2. Per-sample: compute normalised variant burden
For each sample the following script is run to prepare for spectral denoising:
```
Rscript soma-snv-norm.R sampleID sampleID.variants.tsv sampleID.background.tsv sampleID.burden.rds
```

`sampleID` is the sample identifier
`sampleID.variants.tsv` and `sampleID.background.tsv` are the output files from step 1.
`sampleID.burden.rds` is the output from this stage, containing normalised variant class burden estimates.


### 3. All samples: cohort merging
In this step, the per-sample output files from step 2 are merged into a single file, in preparation for spectral denoising.
```
Rscript soma-snv-merge.R merged.burden.rds sampleID1.burden.rds sampleID2.burden.rds ...
```

`merged.burden.rds` is the merged cohort data file.
`sampleID1.burden.rds`, `sampleID2.burden.rds`, etc are the normalised variant burden files for each sample from step 2.


### 4. All samples: spectral denoising
```
Rscript soma-snv-nmf.R merged.burden.rds "${prefix}" <seed> <kmin> <kmax> <B> <cores>
```

### 5. Extraction of scores and burden.


## Generation of a custom blacklist bed

Section under construction.


## Future

Planned extensions to soma-snv:
* A single-sample mode that uses pre-trained spectral decompositions to denoise results.
* A final workflow stage that annotates variants with likelihoods of membership in each of the spectral components identified by NMF.
