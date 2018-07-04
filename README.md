# soma-snv
Detect and quantify somatic variants in low-depth sequencing data


## Overview

Somatic variants are difficult to detect using error-prone massively-parallel sequencing technology.  This problem is particularly challenging for low-depth data, for which low-VAF (variant allele frequency) variants and machine errors are difficult to distinguish.  soma-snv is a tool to address this problem, by giving estimates of somatic variant burden using 30X whole-genome sequencing data.

Conceptually, soma-snv combines a permissive variant detector with a spectral denoising stage.  The detector tracks expected somatic variant senstivity per sample, enabling the normalisation of somatic burden by each sample's detection cross-section.  More details on the detector are available in [docs/detector.pdf](docs/detector.pdf).


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


### 1. Per-sample: variant detection and cross-section estimation
The following code detects putative somatic variants and estimates the background cross-section for each sample:
```bash
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

```bash
cat sampleID.shard*.variants.tsv > sampleID.variants.tsv
cat sampleID.shard*.background.tsv > sampleID.background.tsv
```


### 2. Per-sample: compute normalised variant burden
For each sample the following script is run to prepare for spectral denoising:
```bash
Rscript soma-snv-norm.R sampleID sampleID.variants.tsv sampleID.background.tsv sampleID.burden.rds
 ```

`sampleID` is the sample identifier
`sampleID.variants.tsv` and `sampleID.background.tsv` are the output files from step 1.
`sampleID.burden.rds` is the output from this stage, containing normalised variant class burden estimates.


### 3. All samples: cohort merging
In this step, the per-sample output files from step 2 are merged into a single file:
```bash
Rscript soma-snv-merge.R merged.burden.rds sampleID1.burden.rds sampleID2.burden.rds etc...
```

`merged.burden.rds` is the merged cohort data file.
`sampleID1.burden.rds`, `sampleID2.burden.rds`, etc are the normalised variant burden files for each sample from step 2.


### 4. All samples: spectral denoising and cardinality search
Spectral denoising is performed via non-negative matrix factorisation, by the script `soma-snv-nmf.R`:
```bash
Rscript soma-snv-nmf.R merged.burden.rds <prefix> <seed> <kmin> <kmax> <B> <cores>
```

`prefix` is a prefix path to be used for the output files of this script.
`seed` is a PRNG seed for reproducibility
`kmin` and `kmax` denote bounds for the cardinality search.
`B` is the number of random restarts for the NMF factorization algorithm.
`cores` is the number of cores to use.

`soma-snv-nmf.R` produces a file `<prefix>.<i>.rds` for every `<i>` in [`<kmin>`, `<kmax>`].

As part of the spectral denoising procedure it is necessary to perform a search over the factorization cardinality, to find the best fit to the data.  This search space is defined by the parameters `<kmin>` and `<kmax>`.  Suggested values are `<kmin>` = 2, `<kmax>` = 8, although larger values of kmax may be required for very complex or large datasets; this is discussed further in a later section.

NMF is a non-deterministic algorithm that often requires multiple random restarts to achieve a good fit.  The number of restarts per cardinality is controlled by the `<B>` parameter; in our use `<B>` = 100 was suitable.  An appropriate value of B will give stable results (in terms of evar) for varying seeds, which suggests a principled approach to the selection of B for a given dataset.

### 5. Extraction of scores and burden


## Example

Some files have been provided in the example/ directory to illustrate the use of soma-snv.


### 1. Per-sample: variant detection and cross-section estimation

Due to the size of the BAM input files the input data for this stage is not supplied.  However, output data are present for sequenced DNA from the NA12878 and NA24385 cell lines.  These files were generated as (for NA12878):
```bash
samtools mpileup \
  -d 101 -a -C 50 -E -I -q 30 -Q 30 -s -f hs37d5x.fa NA12878.bam | \
pypy soma-snv.py \
  --blacklist data/blacklist_truseq_nano_hiseqX_hs37d5x.bed --snr 5.0 --vaf 0.2 \
  --error 2.0e-3 --het 1.0e-4 --scv 5.0e-7 --maxdp 100 \
  --background example/NA12878.background.tsv \
> example/NA12878.variants.tsv
```
Likewise for NA24385.


### 2. Per-sample: compute normalised variant burden
```bash
Rscript soma-snv-norm.R NA12878 example/NA12878.variants.tsv example/NA12878.background.tsv example/NA12878.burden.rds
Rscript soma-snv-norm.R NA24385 example/NA24385.variants.tsv example/NA24385.background.tsv example/NA24385.burden.rds
 ```


### 3. All samples: cohort merging
```bash
Rscript soma-snv-merge.R example/merged.burden.rds example/NA12878.burden.rds example/NA24385.burden.rds example/deb905.burden.rds example/1c97ab.burden.rds ...
```
Note that many of the files used in this command are not available, for space and privacy reasons.  However, a full example/merged.burden.rds for 500 samples (including NA12878 and NA24385) is included in this repository, to enable the demonstration of later steps.


### 4. All samples: spectral denoising and cardinality search
This example illustrates the use of multiple restarts of the algorithm with varying seeds to evaluate fit stability.  In this case, 20 different runs are used, testing k from 2 to 6, each time using 40 random restarts of NMF, and 28 cores total.  The output will be a number of files with path `example/nmf.<seed>.<k>.rds`.
```bash
seq 314159 314178 | while read seed; do
    echo "Fitting seed ${seed}..." > /dev/stderr
    Rscript soma-snv-nmf.R example/merged.burden.rds example/nmf.${seed} ${seed} 2 6 40 28
done
```


## Generation of a custom blacklist bed

Section under construction.


## Future

Planned extensions to soma-snv:
* A single-sample mode that uses pre-trained spectral decompositions to denoise results.
* A final workflow stage that annotates variants with likelihoods of membership in each of the spectral components identified by NMF.
