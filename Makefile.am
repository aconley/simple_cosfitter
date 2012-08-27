ACLOCAL_AMFLAGS = -I m4

SUBDIRS=src 

#Extra files not part of build process that need to be included
#in distribution
EXTRA_DIST = sample/lowe_extcorr_param.txt sample/lowe_extcorr_wom_param.txt \
	sample/lowe_param.txt sample/knop03_extended_all.dat \
        sample/knop03_extended.dat sample/knop03_extended_lowe.dat \
        sample/Astier_sne.dat sample/Astier_param.txt \
	sample/Astier_flat_param.txt sample/Astier_wom_param.txt \
	sample/Riess98_sne.dat sample/Kowalski_sne.dat \
	sample/Kowalski_syscovmat.txt sample/Kowalski_wom_param.txt \
        idl/colordex.pro idl/contourshow.pro \
	idl/deepcolor.pro sample/sample_covfile.txt Doxyfile

docs :
	doxygen Doxyfile
