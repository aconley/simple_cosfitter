ACLOCAL_AMFLAGS = -I m4

SUBDIRS=src 

#Extra files not part of build process that need to be included
#in distribution
EXTRA_DIST = sample/Astier_sne.dat sample/Astier_param.txt \
	sample/Astier_flat_param.txt sample/Astier_wom_param.txt \
        idl/colordex.pro idl/contourshow.pro \
	idl/deepcolor.pro sample/sample_covfile.txt Doxyfile
docs :
	doxygen Doxyfile
