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
  * doParallel (tested with 1.0.11)
  * ggplot2 (tested with 2.2.1)
  * plyr (tested with 1.8.4)
  * SomaticSignatures http://bioconductor.org/packages/release/bioc/html/SomaticSignatures.html (tested with 2.14.0)

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


### 4. All samples: spectral denoising

#### 4a. Cardinality search

soma-snv spectral denoising uses non-negative matrix factorisation, which requires as input an estimate of the number of mutational signatures present in the data (the "cardinality" of the factorisation).  Prior to performing the final spectral denoising (4c), we perform a series of trial factorisations in order to estimate a reasonable value for the cardinality.

```bash
seq <seed_start> <seed_end> | while read seed; do
  Rscript soma-snv-nmf.R merged.burden.rds <prefix> ${seed} <kmin> <kmax>
done
```

Or equivalently, in parallel:
```bash
parallel Rscript soma-snv-nmf.R merged.burden.rds <prefix> {} <kmin> <kmax> ::: $(seq <seed_start> <seed_end>)
```

`prefix` is a prefix path to be used for the output files of this script.
`seed` is a PRNG seed for reproducibility
`kmin` and `kmax` denote bounds for the cardinality search.
`seed_start` and `seed_end` denote bounds for the random repeats of the fit, see [Fit stability and multiple restarts](#####fit-stability-and-multiple-restarts).

`soma-snv-nmf.R` produces a file `<prefix>.<seed>.<i>.rds` for every `<i>` in [`<kmin>`, `<kmax>`].

This search space is defined by the parameters `<kmin>` and `<kmax>`.  Suggested values are `<kmin>` = 2, `<kmax>` = 8, although larger values of kmax may be required for very complex or large datasets -- this can be detected from the output of the `soma-snv-choosek.R` script.


##### Fit stability and multiple restarts

NMF relies on local optimisation, and often requires multiple random restarts to achieve a good fit.  Too few restarts leads to a high likelihood of poor fit and incorrect results, but too many is wasteful.  The number of restarts is controlled by the `<seed_start>` and `<seed_end>` parameters: seed_end - seed_start + 1 restarts will be used.  The fit achieved by these restarts is plotted by `soma-snv-choosek.R` (see [4b. Cardinality selection](####-4b.-cardinality-selection)), and can be used to evaluate whether sufficient restarts have been used.  A good starting point is seed_end = seed_start + 100.

##### Increasing search speed

For large cohorts and k, the cardinality search can be slow.  One approach to increase speed is to reduce the number of samples in the burden matrix for the cardinality search only, by dropping columns.  Ensure that the cardinality search is always performed on the same subset of the data.

In special circumstances in which there is a-priori known biological structure within the cohort, merging samples into homogeneous groups dictated by this structure is an effective strategy.  For example, if the samples consist of a number of different cancer types, then consider merging samples into blocks within each cancer type.  As a second example, if the samples vary across a continuous parameter (for example, age), then merge samples into age groups.  The appropriate merging in all cases is a row-wise mean for selected columns of the burden matrix.  If this strategy is used, the chosen k will favour the number of latent components that vary across the variable used to stratify the samples.  This merging strategy is not appropriate if the samples cannot be stratified in some way that reflects a property of interest.

Once an appropriate k is chosen using subset data, `soma-snv-nmf.R` can be re-run on the full data set, using only the final selected value of `k` for both `<kmin>` and `<kmax>`.  Random restarts with different seeds will still be required, and possibly more restarts will be needed than were suggested by the results on the subset data.


#### 4b. Cardinality selection

`soma-nmf-choosek.R` takes the results from 4a and produces a series of plots to aid in the selection of an appropriate cardinality.
```bash
Rscript soma-nmf-choosek.R <output_pdf> <infile1> <infile2> ...
```
The input files are `.rds` files produced by step 4a.  `soma-nmf-choosek.R` extracts some information from the names of these files, so ensure that they conform to the `<prefix>.<seed>.<k>.rds` format produced by `soma-snv-nmf.R`.  The first page of the output pdf is a summary figure of the explained variance of each cardinality, and is the primary tool used to verify an appropriate number of restarts were used, and select `k`.  Later pages supply information which may assist in the selection of a `k` when there is no clear optimum.

##### Verification of global fit

The first plot shows the change in explained variance of the burden matrix as a function of k:

![Selecting the number of restarts](/docs/evar_restarts.png?raw=true "Selecting the number of restarts")

Points (jittered in x) show the explained variance for individual runs of the algorithm with different seeds and k, and the line connects the best fits for each cardinality across all restarts.  A suitably high number of random restarts is indicated if the maximum value of the explained variance for each k appears to be well-sampled: many points are clustered near the maximum value for each k, as seen for example in the left panel.  Conversely, too few restarts are seen in the right panel, where the maximum value for a given k is only supported by one or two points.

Note that if subsetting was performed in 4a (see [Increasing search speed](#####Increasing-search-speed)), more restarts may be required for the full data than is suggested by the subsetted burden results.

##### Selection of k

Assuming sufficient random restarts were used and the global fit is good, the appropriate value of k can be estimated by identifying an inflection point in the optimum explained variance line.  In the example above, this inflection point is at cardinality 3.

There is some subjectivity to the selection of k, in which case the additional plots produced by `soma-nmf-choosek.R` can be useful.  The additional plots show the mutational signatures resulting from the fits, one plot per value of k.  Comparison of these spectra with canonical cancer signatures at [COSMIC](https://cancer.sanger.ac.uk/cosmic/signatures) can sometimes help resolve an equivocal explained variance plot.


### 5. Extraction of scores and burden
The convenience script `soma-snv-extractbest.R` can be used to extract the best fit for a given k, for further analysis:
```bash
Rscript soma-snv-extractbest.R <k> <outfile> <infile1> <infile2> ...
```
`<k>` is the selected best `k` from section 4b.
`<outfile>` is an output `.rds`
`<infile1>`, etc are input files, as were supplied to `soma-nmf-choosek.R`

An example of extracting data from the resultant `<outfile>` is given in [5. Extraction and plotting of results](###-5.extraction-and-plotting-of-results).


## Example

Some files have been provided in the example/ directory to illustrate the use of soma-snv.


### 1. Per-sample: variant detection and cross-section estimation

Due to the size of the BAM input files the input data for this stage is not supplied.  However, output data are given in `example/persample/*` for sequenced DNA from the NA12878 and NA24385 cell lines.  These files were generated as (for NA12878):
```bash
mkdir -p example/persample/NA12878
samtools mpileup \
  -d 101 -a -C 50 -E -I -q 30 -Q 30 -s -f hs37d5x.fa NA12878.bam | \
pypy soma-snv.py \
  --blacklist data/blacklist_truseq_nano_hiseqX_hs37d5x.bed --snr 5.0 --vaf 0.2 \
  --error 2.0e-3 --het 1.0e-4 --scv 5.0e-7 --maxdp 100 \
  --background example/persample/NA12878/NA12878.background.tsv \
> example/persample/NA12878/NA12878.variants.tsv
```
Likewise for NA24385.


### 2. Per-sample: compute normalised variant burden
```bash
Rscript soma-snv-norm.R NA12878 \
  example/persample/NA12878/NA12878.variants.tsv \
  example/persample/NA12878/NA12878.background.tsv \
  example/persample/NA12878/NA12878.burden.rds

Rscript soma-snv-norm.R NA24385 \
  example/persample/NA24385/NA24385.variants.tsv \
  example/persample/NA24385/NA24385.background.tsv \
  example/persample/NA24385/NA24385.burden.rds
 ```


### 3. All samples: cohort merging
```bash
Rscript soma-snv-merge.R example/merged.burden.rds \
  example/NA12878/NA12878.burden.rds \
  example/NA24385/NA24385.burden.rds \
  example/zrre/zrre.burden.rds \
  example/tusi/tusi.burden.rds \
  ...
```
Note that many of the files used in this command are not available, for space and privacy reasons.  However, a full `example/merged.burden.rds` for 500 anonymised samples (including NA12878 and NA24385) is included in this repository, to enable the demonstration of later steps.


### 4. All samples: spectral denoising

#### 4a. Cardinality search

This example illustrates the use of multiple restarts of the algorithm with varying seeds to evaluate fit stability.  In this case, 50 different seed runs are used ([314159, 314208]), testing k from 2 to 6, with GNU parallel for speed.  The output will be 50x5=250 files with path `example/nmf.<seed>.<k>.rds`.
```bash
mkdir -p example/search
parallel Rscript soma-snv-nmf.R example/merged.burden.rds example/search/nmf {} 2 6 ::: $(seq 314159 314208)
```


#### 4b. Cardinality selection
Generate plots to decide an appropriate k:
```bash
Rscript soma-snv-choosek.R example/choosek.pdf example/search/nmf.*.rds
```

The resulting explained variance plot suggests k=3 is an appropriate choice:

![Example: explained variance plot](/docs/evar.png?raw=true "Example: explained variance plot")


### 5. Extraction and plotting of results
`soma-snv-extractbest.R` is a convenience script to extract the best fit for a given k:
```bash
Rscript soma-snv-extractbest.R 3 example/finalfit.k3.rds example/search/nmf.*.rds
```

Signatures, scores, and mutational burden can then be extracted from this file:
```R
library(SomaticSignatures)

fit = readRDS("example/finalfit.k3.rds")

# Get signatures
head(signatures(fit))
#                S1           S2           S3
# CA A.A   1.366513 2.220446e-16 3.271910e+01
# CA A.C  22.920074 2.220446e-16 2.937612e+01
# CA A.G 125.339755 2.220446e-16 2.220446e-16
# CA A.T  29.627604 2.220446e-16 2.341390e+01
# CA C.A   7.266072 2.220446e-16 3.659436e+01
# CA C.C  29.998369 2.220446e-16 2.070437e+01

# Get scores of signatures for each sample
head(samples(fit))
#                   S1           S2           S3
# NA12878 4.848781e-04 0.0016471674 0.0078611438
# NA24385 2.220446e-16 0.0001955538 0.0014927038
# bihjwn  2.220446e-16 0.0002038287 0.0017219708
# jabtjk  2.220446e-16 0.0002768122 0.0006197808
# rhacxe  3.675279e-03 0.0001375060 0.0004815627
# tgxleg  2.220446e-16 0.0006764890 0.0075259051

# Plot signatures
plotSignatures(fit)

# Get the mutation rate for each signature in each sample,
# expressed in somatic variants per megabase.
sample_signature_burden = t(t(samples(fit)) * colSums(signatures(fit)))

head(sample_signature_burden)
#                   S1        S2        S3
# NA12878 2.025931e+00 1.5041075 25.795531
# NA24385 9.277530e-13 0.1785695  4.898153
# bihjwn  9.277530e-13 0.1861258  5.650469
# jabtjk  9.277530e-13 0.2527705  2.033747
# rhacxe  1.535615e+01 0.1255634  1.580198
# tgxleg  9.277530e-13 0.6177345 24.695480
```

## Generation of a custom blacklist bed

Section under construction.


## Future

Planned extensions to soma-snv:
* A single-sample mode that uses pre-trained spectral decompositions to denoise results.
* A final workflow stage that annotates variants with likelihoods of membership in each of the spectral components identified by NMF.
