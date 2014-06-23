[![Build Status](https://travis-ci.org/wolski/imsbInfer.svg?branch=master)](https://travis-ci.org/wolski/imsbInfer)
[![Project Stats](https://www.ohloh.net/p/imsbInfer/widgets/project_thin_badge.gif)](https://www.ohloh.net/p/imsbInfer)

## R-package for the quantitative analysis of SWATH-MS data.

Overwiev

- data import, filtering, transfromations, normalization,  scaling
- visualization
- QC on peptide and protein level

## How to build:

- download the repository as a zip archive (see button __Download Zip__ to the right)
- then run the sequence of commands below:

```sh
unzip imsbInfer-master.zip
mv imsbInfer-master imsbInfer
# to create the package documentation 
cp imsbInfer/runrox2.R .
Rscript runrox2.R 
# install the package
R CMD INSTALL imsbInfer
```

## Description

Protein quantification experiments can be performed using data independent aquistion DIA LC-MS/MS experiments.
imsbInfer provides methods to import data from openSwath - feature-alinger or spectronaut DIA LC-MS/MS software.


* Programmed with S3. 
  * inspecting function code from within R easy (see internals, learn, modify from within your R session).

* Documentation is brief but there is plenty of example code for each function 
  * helps to understand what the functions are doing
  * examples are also used to test the function - so there are is a lot of stopifnot statements in the example code sections

* filter :
  * Filter top transitions and peptides and to aggregate them.
  * remove decoys, retention time ranges (i.e. begin or end of gradient), filter for m_score.

* Methods to scale/transform intensities.
  - robust scaling of individual samples
  - robust scaling of samples taking RT into account
  - median scaling (same median intensity given RT as in reference sample)

* Assessing correlation between peptides of same protein or different proteins can guide the choice of data scaling method:
  - for randomly chosen peptides one would expect on average no corration.
  - for peptides from same protein postive correlation is expected
  - for same peptide and a different charge high correlation is expected

* Export to MSstats [http://www.msstats.org/](http://www.msstats.org/)

A lot of cool stuff to do visualization of the data.

* various flawors of pairplots : 
  - scatter plot
  - foldchanges against RT
  - qqplots

* imageplot to visualize 
  - correlations among samples
  - similarity of distributions among samples (ks.test).
