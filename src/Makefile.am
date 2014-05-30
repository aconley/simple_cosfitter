
lib_LTLIBRARIES = libsimple_cosfitter.la libsimple_light.la

libsimple_cosfitter_la_SOURCES = \
	cosgrids.cc fitter.cc paramfile.cc snedata.cc param_results.cc \
	lumdist.cc cosfitterexcept.cc utility.cc covmatrix.cc \
	lumdist_array.cc covmatrix_woodbury.cc \
	cosgrids.h fitter.h paramfile.h snedata.h param_results.h lumdist.h \
	cosfitterexcept.h utility.h mainpage.h param_tags.h covmatrix.h \
	lumdist_array.h covmatrix_woodbury.h

libsimple_light_la_SOURCES = \
	paramfile.cc cosgrids.cc utility.cc cosfitterexcept.cc snedata.cc \
	covmatrix.cc covmatrix_woodbury.cc \
	paramfile.h cosgrids.h utility.h cosfitterexcept.h snedata.h \
	param_tags.h covmatrix.h covmatrix_woodbury.h 

if USECFITSIO 
bin_PROGRAMS = runfitter combine_prob get_lumdist indivchi \
		baryonpeak_grid wmap_grid w2df_grid convert_prob_to_fits \
		convert_to_fits get_wmap_shift omprior_grid wmap9_grid \
		percival_grid get_acoustic
else
bin_PROGRAMS = runfitter combine_prob get_lumdist indivchi omprior_grid \
		baryonpeak_grid wmap_grid w2df_grid get_wmap_shift wmap9_grid \
		percival_grid
endif

check_PROGRAMS = lumdisttest snedatatest covmatrixtest cosgridtest

runfitter_SOURCES = runfitter.cc
runfitter_LDFLAGS = -lsimple_cosfitter -lgsl -lgslcblas

dnl getchi_SOURCES = getchi.cc
dnl getchi_LDFLAGS = -lsimple_cosfitter

combine_prob_SOURCES = combine_prob.cc
combine_prob_LDFLAGS = -lsimple_cosfitter

get_lumdist_SOURCES = get_lumdist.cc
get_lumdist_LDFLAGS = -lsimple_cosfitter

indivchi_SOURCES = indivchi.cc
indivchi_LDFLAGS = -lsimple_cosfitter

baryonpeak_grid_SOURCES = auxconstraint.cc baryonpeak_grid.cc auxconstraint.h
baryonpeak_grid_LDFLAGS = -lsimple_light

wmap_grid_SOURCES = auxconstraint.cc wmap_grid.cc auxconstraint.h
wmap_grid_LDFLAGS = -lsimple_light

wmap9_grid_SOURCES = auxconstraint.cc wmap9_grid.cc auxconstraint.h
wmap9_grid_LDFLAGS = -lsimple_light

percival_grid_SOURCES = auxconstraint.cc percival_grid.cc auxconstraint.h
percival_grid_LDFLAGS = -lsimple_light

omprior_grid_SOURCES = auxconstraint.cc omprior_grid.cc auxconstraint.h
omprior_grid_LDFLAGS = -lsimple_light

w2df_grid_SOURCES = w2df.cc w2df_grid.cc w2df.h
w2df_grid_LDFLAGS = -lsimple_light

get_acoustic_SOURCES = auxconstraint.cc get_acoustic.cc auxconstraint.h
get_acoustic_LDFLAGS = -lsimple_light

get_wmap_shift_SOURCES = auxconstraint.cc get_wmap_shift.cc auxconstraint.h
get_wmap_shift_LDFLAGS = -lsimple_light

cosgridtest_SOURCES = cosgridtest.cc
cosgridtest_LDFLAGS = -lsimple_cosfitter

covmatrixtest_SOURCES = covmatrixtest.cc
covmatrixtest_LDFLAGS = -lsimple_cosfitter

lumdisttest_SOURCES = lumdisttest.cc
lumdisttest_LDFLAGS = -lsimple_cosfitter

snedatatest_SOURCES = snedatatest.cc
snedatatest_LDFLAGS = -lsimple_cosfitter

convert_prob_to_fits_SOURCES = convert_prob_to_fits.cc
convert_prob_to_fits_LDFLAGS = -lsimple_light

convert_to_fits_SOURCES = convert_to_fits.cc
convert_to_fits_LDFLAGS = -lsimple_light

docs :
	doxygen Doxyfile
