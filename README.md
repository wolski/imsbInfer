[![Build Status](https://travis-ci.org/wolski/imsbInfer.svg?branch=master)](https://travis-ci.org/wolski/imsbInfer)
[![Project Stats](https://www.ohloh.net/p/imsbInfer/widgets/project_thin_badge.gif)](https://www.ohloh.net/p/imsbInfer)

R-package for the quantitative analysis of SWATH-MS data.

- data import, filtering, transfromations
- visualization
- QC on peptide and protein level
- data normalization,  scaling

How to build:

- download the repository as a zip archive (see link on the left)
- extract the zip file into the folder imsbInfer


cp imsbInfer/runrox2.R .
Rscript runrox2.R
R CMD build imsbInfer
R CMD install imsbInfer