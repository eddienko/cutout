# FITS cutouts

This software produces cutouts from FITS images.

## Installation

### Required libraries

  * [CFITSIO](http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html) - library to read and write FITS files.
  * [AST](http://starlink.jach.hawaii.edu/starlink/AST) -- library to handle worl coordinate systems and image transformation
  * OpenMP -- a gcc compiler with openmp support

## Usage

The basic usage to extract a cutout from an image is

	cutout image.fit cutout.hdr cutout.fit

where `image.fit` is the image to extract the cutout from and `cutout.hdr` contains the header of the desired extracted cutout provided by the user, e.g.

	BITPIX  =                  -32 / number of bits per data pixel
	NAXIS   =                    2 / number of data axes
	NAXIS1  =                  175 / length of data axis 1
	NAXIS2  =                  175 / length of data axis 2
	CRVAL1  =        157.010000 / Reference value on 1st axis in primary WCS
	CRVAL2  =        -19.271000 / Reference value on 2nd axis in primary WCS
	CRPIX1  =               87.500000 / Reference pixel on 1st axis in primary WCS
	CRPIX2  =               87.500000 / Reference pixel on 2nd axis in primary WCS
	CD1_1   =           -9.479372e-05 / Transformation matrix for primary WCS
	CD2_2   =           9.479372e-05 / Transformation matrix for primary WCS
	CTYPE1  = 'RA---TAN'           / &
	CTYPE2  = 'DEC--TAN'           / &
	RADESYS = 'ICRS    '           / Reference frame for RA/DEC values

To extract multiple cutouts one can do

	cutout driver

where `driver` is a text file containing one line per desired cutout, e.g.

	image.fit[1] cutout1.hdr cutout1.fit
	image.fit[1] cutout2.hdr cutout2.fit
	image.fit[2] cutout3.hdr cutout3.fit

Note that image extensions have to be specified (unless the image contains a single extension). In this case `cutout` makes use of the multiprocessor architectore and makes the cutouts in parallel using the available number of CPUs or the ones set by the environment variable `OMP_NUM_THREADS`,

	env OMP_NUM_THREADS=4 cutout driver

with example output:

	Number of processors available = 12
	Number of threads =              4
	
	Processor 0  :  o20131211_00044.fit[8] -> !faaf25c3751c930362542af6535aa074.fit
	Processor 2  :  o20131211_00045.fit[8] -> !bfa9e66fac18cea336f5c7b3c13cf9de.fit
	Processor 1  :  o20131211_00052.fit[8] -> !51eadf5cd3359f94d8ac1787c5cebf4e.fit
	Processor 3  :  o20131211_00058.fit[8] -> !a103afa15ddc1f4cb182c7fa63e026c2.fit
	Processor 2  :  o20140221_00048.fit[8] -> !bdde9d7b61a9170a97a8cb56e9c29238.fit
	Processor 0  :  o20140221_00049.fit[21] -> !62374224036691df973143bafa2b2091.fit
	Processor 1  :  o20140221_00055.fit[8] -> !511f7afec4d6d9c0be73a51cac068098.fit
	Processor 3  :  o20140221_00056.fit[8] -> !f80541da80714a93ad32a14cf0c10971.fit
	Processor 0  :  o20140221_00057.fit[21] -> !20ebdf286498406a3bb47df6ef5cbf40.fit
	Processor 2  :  o20140221_00063.fit[8] -> !98b9613779bef8255013a0a56e9c358f.fit
	Processor 1  :  o20140221_00064.fit[21] -> !56cda30fdc0fb57269470d918574d88a.fit


## License

* see [LICENSE](https://github.com/eddienko/cutout/blob/master/LICENSE) file.

## Version

## Contact
