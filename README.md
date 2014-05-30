[![Build Status](https://travis-ci.org/wolski/imsbInfer.svg?branch=master)](https://travis-ci.org/wolski/imsbInfer)
[![Project Stats](https://www.ohloh.net/p/imsbInfer/widgets/project_thin_badge.gif)](https://www.ohloh.net/p/imsbInfer)

R-package for the quantitative analysis of SWATH-MS data.

- data import, filtering, transfromations
- visualization
- QC on peptide and protein level
- data normalization,  scaling

How to build:

- download the repository as a zip archive (see button _Download Zip_ to the right)
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