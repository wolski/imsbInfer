[![Build Status](https://travis-ci.org/wolski/imsbInfer.svg?branch=master)](https://travis-ci.org/wolski/imsbInfer)
[![Project Stats](https://www.ohloh.net/p/imsbInfer/widgets/project_thin_badge.gif)](https://www.ohloh.net/p/imsbInfer)

## R-package for the quantitative analysis of SWATH-MS data.

This project is __DISCONTINUED__. 

New and related work can be found here:

https://github.com/protviz/

Overwiev

- data import, filtering, transfromations, normalization,  scaling
- visualization (large part of visualization has moved to package [quantable](https://github.com/protViz/quantable)
- QC on peptide and protein level

## How to install:
```sh
#setup dependencies
install.packages("devtools")
packs <- c("roxygen2","data.table","gplots","reshape2","scales","RColorBrewer")
install.packages(packs)
library(devtools)
install_github("wolski/quantable")

#install the package from github
install_github("wolski/imsbInfer")
```

## Description

Protein quantification experiments can be performed using data independent aquistion DIA LC-MS/MS experiments.
imsbInfer provides methods to import data from openSwath - feature-alinger or spectronaut DIA LC-MS/MS software.

Objects and functions (S3) to import, export, summarize, visualize and analyze quantification results of MS proteomic LFQ experiments. This package is used to generate QC reports, and to suggest data scaling and transformation to improve data quality and reduce batch effects.

* Programmed with S3. 
  * inspecting function code from within R easy (just type function name without arguments -> see internals, learn, modify from within your R session).

* Documentation is brief but there is plenty of examples for each function
  * Should help to understand how to use function and what the functions are doing
  * examples are also used to test the function - so there are some _stopifnot_ statements in the example code sections
  * There are vignettes (R jargon for documentation)

* Filtering and aggregation functions:
  * Filter top transitions and peptides and aggregate them
  * remove decoys, retention time ranges (i.e. begin or end of gradient), filter for m_score.

* Methods to scale/transform intensities.
  - robust scaling of individual samples (robust __z__ transform)
  - robust scaling of samples taking RT into account (running __z__ transform) - experimental
  - median scaling (same median intensity given RT as in reference sample) - experimental

* Evaluation of normalization methods - Assessing correlation between peptides of same protein or different proteins can guide the choice of data scaling and normalization method:
  - for randomly chosen peptides one would expect on average **no** corration.
  - for peptides from same protein **postive** correlation is expected
  - for same peptide and a different charge **high positive** correlation is expected

* Export to MSstats [http://www.msstats.org/](http://www.msstats.org/)

A lot of cool stuff to do visualization of the data (mostly moved to package [quantable](https://github.com/protViz/quantable).

* various flawors of pairplots : 
  - scatter plot (quantable)
  - foldchanges against RT
  - qqplots

* imageplot to visualize 
  - correlations among samples
  - similarity of distributions among samples (ks.test).

* nice looking volcanoplot (moved to package [quantable](https://github.com/protviz/quantable))

* Availability of prebuild packages:
  - Extensive package checks and builds are run after each change to the package.
  - Therefore the most recent version of the package is the recommended one and the best way to get the package is by executing reading section how to install above.
  - The time to run all examples exceeds time limits at CRAN and BIOCONDUCTOR.
  - Still we provide regular releases for backward compatibility which can be downloaded from : [https://github.com/wolski/imsbInfer/releases](https://github.com/wolski/imsbInfer/releases)

# Related work

* [SWATH2stats](http://www.bioconductor.org/packages/release/bioc/html/SWATH2stats.html)
* [aLFQ](https://cran.r-project.org/web/packages/aLFQ/index.html)
