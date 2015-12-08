build: makedirs roxygen
	R CMD build Bluejay
	mv Bluejay_*.tar.gz build/

makedirs:
	mkdir -p build

clean:
	rm -rf build

roxygen:
	R CMD BATCH --no-save  pkgsrc/roxygenize.R /dev/stdout

vignette: makedirs
	cp -rf Bluejay/vignettes build
	cp pkgsrc/vignette/Makefile build/vignettes
	make -C build/vignettes build

manual: makedirs
	R CMD Rd2pdf --no-preview Bluejay
	mv -f Bluejay.pdf build/

