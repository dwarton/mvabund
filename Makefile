PACKAGE := $(shell grep '^Package:' DESCRIPTION | sed -E 's/^Package:[[:space:]]+//')
RSCRIPT = Rscript --no-init-file

install:
	R CMD INSTALL .

build: compile
	R CMD build .

check: build
	R CMD check --no-manual `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -f `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -rf ${PACKAGE}.Rcheck


compile: attributes
	Rscript -e 'devtools::compile_dll()'

vignettes:
	Rscript -e "library('methods'); devtools::build_vignettes()"

test:
	Rscript -e 'library(methods); devtools::test()'

attributes:
	Rscript -e "Rcpp::compileAttributes()"

.PHONY: test attributes install build check vignettes
