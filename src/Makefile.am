
lib_LTLIBRARIES = libsimple_cosfitter_twoscriptm.la \
	libsimple_light_twoscriptm.la

libsimple_cosfitter_twoscriptm_la_SOURCES = \
	cosgrids.cc fitter.cc paramfile.cc snedata.cc param_results.cc \
	lumdist.cc cosfitterexcept.cc utility.cc covmatrix.cc \
	lumdist_array.cc covmatrix_woodbury.cc \
	cosgrids.h fitter.h paramfile.h snedata.h param_results.h lumdist.h \
	cosfitterexcept.h utility.h mainpage.h param_tags.h covmatrix.h \
	lumdist_array.h covmatrix_woodbury.h

libsimple_light_twoscriptm_la_SOURCES = \
	paramfile.cc cosgrids.cc utility.cc cosfitterexcept.cc snedata.cc \
	covmatrix.cc covmatrix_woodbury.cc \
	paramfile.h cosgrids.h utility.h cosfitterexcept.h snedata.h \
	param_tags.h covmatrix.h covmatrix_woodbury.h 

if USECFITSIO 
bin_PROGRAMS = runfitter_twoscriptm indivchi_twoscriptm \
	convert_to_fits convert_prob_to_fits
else
bin_PROGRAMS = runfitter_twoscriptm indivchi_twoscriptm
endif

runfitter_twoscriptm_SOURCES = runfitter.cc
runfitter_twoscriptm_LDFLAGS = -lsimple_cosfitter_twoscriptm -lgsl -lgslcblas

indivchi_twoscriptm_SOURCES = indivchi_twoscriptm.cc
indivchi_twoscriptm_LDFLAGS = -lsimple_cosfitter_twoscriptm

convert_prob_to_fits_SOURCES = convert_prob_to_fits.cc
convert_prob_to_fits_LDFLAGS = -lsimple_light_twoscriptm

convert_to_fits_SOURCES = convert_to_fits.cc
convert_to_fits_LDFLAGS = -lsimple_light_twoscriptm

docs :
	doxygen Doxyfile
