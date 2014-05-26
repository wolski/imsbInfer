# prepare the package for release
#http://yihui.name/en/2013/04/travis-ci-general-purpose/

PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: check clean

deps:
	Rscript -e 'if (!require("Rd2roxygen")) install.packages("Rd2roxygen", repos="http://cran.rstudio.com")'

docs:
	R -q -e 'library(Rd2roxygen); rab(".", build = FALSE)'

build:
	cd ..;\
	R CMD BATCH imsbInfer/runrox2.R;\
	R CMD build --no-manual --no-rebuild-vignettes $(PKGSRC)

install: build
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --as-cran

travis: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --no-manual --no-rebuild-vignettes

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/

