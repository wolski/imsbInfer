[![Build Status](https://travis-ci.org/wolski/imsbInfer.svg?branch=master)](https://travis-ci.org/wolski/imsbInfer)
[![Project Stats](https://www.ohloh.net/p/imsbInfer/widgets/project_thin_badge.gif)](https://www.ohloh.net/p/imsbInfer)

## R-package for the quantitative analysis of SWATH-MS data.

Overwiev

- data import, filtering, transfromations, normalization,  scaling
- visualization
- QC on peptide and protein level

## How to install:
```sh
install.packages("devtools")
library(devtools)
install_github("wolski/imsbInfer")
```

## Description

Protein quantification experiments can be performed using data independent aquistion DIA LC-MS/MS experiments.
imsbInfer provides methods to import data from openSwath - feature-alinger or spectronaut DIA LC-MS/MS software.

Objects and functions (S3) to import, export, summarize, visualize and analyze quantification results of MS proteomic LFQ experiments. This package is used to generate QC reports, and to suggest data scaling and transformation to improve data quality and reduce batch effects.

* Programmed with S3. 
  * inspecting function code from within R easy (see internals, learn, modify from within your R session).

* Documentation is brief but there is plenty of examples for each function
  * Hopefullys help to understand how to use it and what the functions are doing
  * examples are also used to test the function - so there are some _stopifnot_ statements in the example code sections
  * There are vignettes (R jargon for documentation)

* Filtering and aggregation functions:
  * Filter top transitions and peptides and aggregate them.
  * remove decoys, retention time ranges (i.e. begin or end of gradient), filter for m_score.

* Methods to scale/transform intensities.
  - robust scaling of individual samples (robust __z__ transform)
  - robust scaling of samples taking RT into account (running __z__ transform)
  - median scaling (same median intensity given RT as in reference sample)

* Evaluation of normalization methods - Assessing correlation between peptides of same protein or different proteins can guide the choice of data scaling method:
  - for randomly chosen peptides one would expect on average **no** corration.
  - for peptides from same protein **postive** correlation is expected
  - for same peptide and a different charge **high positive** correlation is expected

* Export to MSstats [http://www.msstats.org/](http://www.msstats.org/)

A lot of cool stuff to do visualization of the data.

* various flawors of pairplots : 
  - scatter plot
  - foldchanges against RT
  - qqplots

* imageplot to visualize 
  - correlations among samples
  - similarity of distributions among samples (ks.test).

* nice looking volcanoplot

* Availability of prebuild packages:
  - Example code is run with every package check and build. The time to run all examples exceeds time limits at CRAN and BIOCONDUCTOR. Therefore,
  prebuild packages can be downloaded from [http://fgcz-data.uzh.ch/~cpanse/imsbInfer](http://fgcz-data.uzh.ch/~cpanse/imsbInfer)


