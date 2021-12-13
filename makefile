default:
	@echo targets: all [roxy, install]  build install test

all:  roxy install

roxy:
	R -e "devtools::document()"

vig:
	R -e "devtools::build_vignettes()"

build:
	R CMD build --no-build-vignettes .

install:
	R CMD INSTALL --no-test-load .


test:
	(cd inst/unitTests; make)


