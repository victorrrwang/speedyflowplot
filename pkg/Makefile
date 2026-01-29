all: build

build:
	R CMD build .
	mv *.tar.gz ../releases/

check: build
	R CMD check --as-cran ../releases/*.tar.gz

install: build
	R CMD INSTALL ../releases/*.tar.gz

clean:
	rm -rf ../releases/speedyflowplot_*.tar.gz
	rm -rf *.Rcheck
	rm -rf src/*.o
	rm -rf src/*.so
