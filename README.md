[![Build Status](https://travis-ci.org/wolski/imsbInfer.svg?branch=master)](https://travis-ci.org/wolski/imsbInfer)
[![Project Stats](https://www.ohloh.net/p/imsbInfer/widgets/project_thin_badge.gif)](https://www.ohloh.net/p/imsbInfer)

R-package for the quantitative analysis of SWATH-MS data.

- data import, filtering, transfromations, normalization,  scaling
- visualization
- QC on peptide and protein level

How to build:

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

Description:

Protein quantification experiments can be performed using data independent aquistion experiments DIA LC-MS experiments \cite{ludovic}.

schwatz provides methods to import data from openSwath/featurealinger or spectronaut analysis software.

Methods to filter top transtions and peptides,  and to aggregate transitions and peptides.

Easy to remove decoys, retention time ranges (begin or end of gradient), filter for mscore.

Programmed with S3. 
    - Clean functional programming style
    - insepecting function code from within R easy (makes seeing internals easy).


Documentation is brief but there is plenty of example code for each function 
    - helps to understand what the functions are doing
    - examples are also used to test the function - so there are is a lot of stopifnot statements in the example code sections


Methods to scale/transform intensities.
    - robust scaling of individual samples
    - robust scaling of samples taking RT into account
    - median scaling (same median intensity given RT as in reference sample)

A lot of cool stuff to do visualization of the data.

different kind of pairplots : 
    - scatter plot
    - foldchanges/differences between samples against RT
    - qqplots

imageplot to visualize 
    - correlations among samples
    - similarity of distributions among (ks.test) samples:


Assessing correlation between transitions of same peptide, peptides of same protein
can guide the choice of data scaling method.
    - for randomly chosen transtions/peptides one would expect on average no corration.
    - for peptides from same protein one would expect postive correlation
    - for same peptide and a different charge expect high correlation

