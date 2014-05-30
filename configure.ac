#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(simple_cosfitter_twoscriptm, 1.6.14, alexander.conley@colorado.edu)
AC_CONFIG_SRCDIR([src])
AM_INIT_AUTOMAKE( [1.9] )
AC_CONFIG_MACRO_DIR([m4])

AC_LANG_CPLUSPLUS

# Checks for programs.
AC_PROG_CXX
AC_PROG_LIBTOOL

dnl Initialize CXXFLAGS and CPPFLAGS
dnl CPPFLAGS are preprocessor flags, CXX are compiler flags
dnl Unless --enable-debug is set, then -O3 is used
AC_ARG_WITH(cppflags,
            [  --with-cppflags=VAL     CPPFLAGS value (def: )],
            [CPPFLAGS="$withval"],[CPPFLAGS=""] )
AC_ARG_WITH(cxxflags,
            [  --with-cxxflags=VAL     CXXFLAGS value (def: -Wall)],
            [ CXXFLAGS=$withval ], [CXXFLAGS="-Wall -D_LINUX"])
AC_ARG_ENABLE(debug,
              [  --enable-debug          Add -g to CXXFLAGS (else add -O3)],
              [ CXXFLAGS="-g $CXXFLAGS"], [ CXXFLAGS="$CXXFLAGS -O3"])
AC_ARG_ENABLE(profile,
              [  --enable-profile        Add gprof-friendly stuff to CXXFLAGS],
              [ CXXFLAGS="$CXXFLAGS -pg -fprofile-arcs" ],)

dnl Use rpath so you don't have to set LD_LIBRARY_PATH
AC_ARG_ENABLE( rpath,
  [  --enable-rpath          Use rpath so you don't have to set LD_LIBRARY_PATH ],
  rpathset=$enableval, rpathset="" )

dnl Timing test
AC_ARG_ENABLE( timing,
	       [  --enable-timing         Adds timing tests to some types of fits],
	       timingset=$enableval,timingset="" )
if test $timingset; then
   AC_MSG_RESULT(Compiling timing output)
   AM_CONDITIONAL( TIMING, true )
else
   AC_MSG_RESULT(Not compiling timing output)
   AM_CONDITIONAL( TIMING, false )
fi

dnl Handle cfitsio library stuff
AC_LANG_C
AC_ARG_ENABLE( cfitswrite,
	       [  --enable-cfitswrite     Adds ability to output to FITS files],
	       cfitswriteset=$enableval,cfitswriteset="" )
AC_ARG_ENABLE( covfits,
	       [  --enable-covfits        Adds ability to read FITS cov matrices],
	       covfits=$enableval,covfits="" )
if ( (test $cfitswriteset;) || (test $covfits;)) then
#if test $cfitswriteset; then
   if test $cfitswriteset; then
      AC_MSG_RESULT(Compiling CFITSIO based output to FITS files)
   fi
   if test $covfits; then
      AC_MSG_RESULT(Compiling CFITSIO based reading of FITS cov matrices)
      AC_DEFINE( COVFITS, 1 )
   fi
   dnl We'll need to possible handle a user specified CFITSIO location
   dnl We can't use --with-cfitsio because autoconf gets confused
   dnl with the library name
   AC_ARG_WITH( cfitsiolib, 
           [  --with-cfitsiolib=DIR   DIR=-LDIR, has libcfitsio.a],
           cfitsiolib="$withval/" , cfitsiolib= )
   if test $cfitsiolib; then
     if test $rpathset; then
      	  LDFLAGS="$LDFLAGS -L$cfitsiolib -Wl,-rpath -Wl,$cfitsiolib"
     else
          LDFLAGS="$LDFLAGS -L$cfitsiolib"
     fi
   fi
   AC_CHECK_LIB(cfitsio,main,,
            [AC_MSG_ERROR("Can't find libcfitsio; use --with-cfitsiolib")], 
	    [-lm "-L$cfitsiolib"])
   dnl And cfitsio headers
   AC_ARG_WITH(cfitsioinc,
        [  --with-cfitsioinc=VAL   VAL=-IVAL, has cfitsio include files ],
	[cfitsioinc="$withval/" CPPFLAGS="$CPPFLAGS -I$cfitsioinc"],
        [cfitsioinc=] )
   AC_CHECK_HEADER(fitsio.h,,
     [ if test $cfitsioinc; then
      AC_MSG_ERROR("cfitsio headers not in location specified by --with-cfitsioinc")
      elif test $cfitsiolib; then
       AC_CHECK_HEADER($cfitsiolib/../include/fitsio.h,
        [cfitsioinc="$cfitsiolib/../include/" CPPFLAGS="$CPPFLAGS -I$cfitsioinc"],
         AC_MSG_ERROR("Can't find cfitsio includes; use --with-cfitsioinc"))
       else
         AC_MSG_ERROR("Can't find cfitsio includes; use --with-cfitsioinc")
      fi ] )
    if test $cfitswriteset; then    
        AM_CONDITIONAL( USECFITSIO, true )
    else
        AM_CONDITIONAL( USECFITSIO, false )
    fi
    if test $covfits; then    
        AM_CONDITIONAL( USECOVFITS, true )
    else
        AM_CONDITIONAL( USECOVFITS, false )
    fi
else 
   AC_MSG_RESULT(Not compiling CFITSIO support)
   AM_CONDITIONAL( USECFITSIO, false )
   AM_CONDITIONAL( USECOVFITS, false )
fi  

dnl GSL
AC_LANG_C
AC_ARG_WITH(gsllib, 
            [  --with-gsllib=DIR       DIR=-LDIR, has libgsl and libgslcblas],
            gsllib="$withval/",gsllib= )
if test $gsllib; then
  AC_MSG_RESULT( Checking User specified GSL )
  if test $rpathset; then
     LDFLAGS="$LDFLAGS -L$gsllib -Wl,-rpath -Wl,$gsllib"
  else
     LDFLAGS="$LDFLAGS -L$gsllib"
  fi
  AC_CHECK_LIB( gsl, gsl_odeiv_step_alloc, ,
      AC_MSG_ERROR("libgsl not in location specified by --with-gsllib: $gsllib"),
      [ -lm -lgslcblas -L$gsllib ] )
  AC_CHECK_LIB( gslcblas, cblas_saxpy, ,
      AC_MSG_ERROR("libgslcblas not in location specified by --with-gsllib: $gsllib"),
      [ -lm -L$gsllib -lgslcblas ] )
else
  AC_MSG_RESULT( Checking GSL )
  AC_CHECK_LIB( gsl, gsl_odeiv_step_alloc, ,
      AC_MSG_ERROR("libgsl not found; use --with-gsllib"), [-lm -lgslcblas] )
  AC_CHECK_LIB( gslcblas, cblas_saxpy, ,
      AC_MSG_ERROR("libgslcblas not found; use --with-gsllib"), [-lm -lgslcblas lgsl ] )
fi

AC_ARG_WITH(gslinc,
        [  --with-gslinc=VAL       VAL=-IVAL, has gsl include files ],
	[gslinc="$withval/" CPPFLAGS="$CPPFLAGS -I$gslinc"],
        [gslinc=] )
AC_MSG_RESULT( Checking GSL Headers )
AC_CHECK_HEADERS([gsl/gsl_errno.h gsl/gsl_sf_ellint.h gsl/gsl_integration.h gsl/gsl_odeiv.h ],,
  [ if test $gslinc; then
   AC_MSG_ERROR("gsl headers not in location specified by --with-gslinc")
   elif test $gsllib; then
    AC_CHECK_HEADER($gsllib/../include/gsl_sf_ellint.h,
     [gslinc="$gsllib/../include/" CPPFLAGS="$CPPFLAGS -I$gslinc"],
      AC_MSG_ERROR("Can't find gsl includes in user specified location; use --with-gslinc"))
    else
      AC_MSG_ERROR("Can't find gsl includes; use --with-gslinc")
   fi ] )

dnl Handle Mac OS X Accelerate framework
AC_ARG_ENABLE( accelerate,    
   [  --enable-accelerate     Use OS X Accelerate to do matrix inversion ],
   useaccelerateset=$enableval, useaccelerateset="" )
if test $useaccelerateset; then
   AC_MSG_RESULT( Using OSX Accelerate framework for matrix inversion )
   AC_DEFINE( USEACCELERATE, 1 )
   LDFLAGS="$LDFLAGS -framework Accelerate"
fi

dnl Handle ATLAS 
AC_LANG_C
AC_ARG_ENABLE( atlas, 
	  [  --enable-atlas          Use ATLAS to do matrix inversion ],
	  useatlasset=$enableval,useatlasset="" )

if test $useatlasset; then
   AC_MSG_RESULT(Using ATLAS for matrix inversion)
   dnl Set so code understands
   AC_DEFINE( USEATLAS, 1 )  

   if test $useaccelerateset; then
     AC_MSG_ERROR( Can't use both Accelerate and ATLAS )
   fi

   dnl User specified path to ATLAS
   AC_ARG_WITH( atlaslib, 
           [  --with-atlaslib=DIR     DIR=-LDIR, has libatlas, etc.],
           atlaslib="$withval/", atlaslib= )
   if test $atlaslib; then
    if test $rpathset; then
       LDFLAGS="$LDFLAGS -L$atlaslib -Wl,-rpath -Wl,$atlaslib"
    else
     LDFLAGS="$LDFLAGS -L$atlaslib"
    fi
   fi

   dnl User specified headers
   AC_ARG_WITH(atlasinc,
     [  --with-atlasinc=VAL     VAL=-IVAL, has atlas include files (clapack.h) ],
      [ atlasinc="$withval/" CPPFLAGS="$CPPFLAGS -I$atlasinc"],
       [atlasinc=] )

	   	   
   dnl Now check for headers	
   AC_MSG_RESULT( Checking for ATLAS headers )  
   AC_CHECK_HEADER(clapack.h,,
    [ if test $atlasinc; then
       AC_MSG_ERROR("ATLAS headers not in location specified by --with-atlasinc: $atlasinc")
      elif test $atlaslib; then
       AC_CHECK_HEADER($atlaslib/../../include/clapack.h,
        [atlasinc="$atlaslib/../../include/" CPPFLAGS="$CPPFLAGS -I$atlasinc"],
         AC_MSG_ERROR("Can't find ATLAS includes; use --with-atlasinc"))
       else
         AC_MSG_ERROR("Can't find ATLAS includes; use --with-atlasinc")
      fi ] )

   dnl Apppend -l commands            
   LDFLAGS="$LDFLAGS -llapack -lcblas -lf77blas -latlas"
fi

dnl MKL
AC_LANG_C
AC_ARG_ENABLE( intel,
	[  --enable-intel          Use the Intel MKL to do matrix inversion],
	useintelset=$enableval, useintelset="")
AC_ARG_ENABLE( intel64,
	[  --enable-intel64        Use 64 bit Intel MKL for matrix inversion],
	useintel64set=$enableval, useintel64="")
AC_ARG_ENABLE( inteluniversal,
	[  --enable-inteluniversal Use universal binary Intel MKL for matrix inversion (on a mac)],
	useinteluniversalset=$enableval, useinteluniversal="")

if test $useintelset || test $useintel64set || test $useinteluniversalset; then
   if test $useatlas; then
      AC_MSG_ERROR(Can't use ATLAS and MKL both)
   fi
   if test $useaccelerateset; then
     AC_MSG_ERROR( Can't use both Accelerate and MKL )
   fi

   AC_MSG_RESULT(Using Intel MKL for matrix inversion)

   dnl Set so code understands
   AC_DEFINE( USEMKL, 1 )  

   dnl Get intel directory
   AC_ARG_WITH(inteldir,
	[  --with-inteldir=VAL     VAL=-IVAL, base directory for intel installation],
	[ inteldir="$withval" ], [inteldir= ] )

   dnl Now add directories to path
   if test $inteldir; then
        MKLINCDIR="$inteldir/mkl/include"
	INTELINCDIR="$inteldir/include"
	if test $useintelset; then
	   	AC_MSG_RESULT(Looking for 32 bit MKL)
	   	INTELLIBDIR="$inteldir/lib/ia32/"
		MKLLIBDIR="$inteldir/mkl/lib/32/"
	elif test $useintel64set; then
	   	AC_MSG_RESULT(Looking for 64 bit MKL)
		INTELLIBDIR="$inteldir/lib/ia64/"
		MKLLIBDIR="$inteldir/mkl/lib/64/"
	else 
	   	AC_MSG_RESULT(Looking for universal bit MKL)
		INTELLIBDIR="$inteldir/lib/"
		MKLLIBDIR="$inteldir/Frameworks/mkl/lib/universal/"
	fi

	CPPFLAGS="$CPPFLAGS -I$INTELINCDIR -I$MKLINCDIR"
	LDFLAGS="$LDFLAGS -L$INTELLIBDIR -L$MKLLIBDIR"	
   fi

   dnl Now check for headers	
   AC_MSG_RESULT( Checking for MKL headers )  
   AC_CHECK_HEADER(mkl_lapack.h,,
     AC_MSG_ERROR("mkl headers not in user specified location"))

   dnl Add -l values
   if test $rpathset; then
      LDFLAGS="$LDFLAGS -Wl,-rpath -Wl,-lmkl -Wl,-rpath -Wl,-lmkl_lapack"
      LDFLAGS="$LDFLAGS -Wl,-rpath -Wl,-lmkl_intel"
      LDFLAGS="$LDFLAGS -Wl,--start-group -lmkl_intel_thread"
      LDFLAGS="$LDFLAGS -lmkl_core -Wl,--end-group -Wl,-rpath -Wl,-lguide"
      LDFLAGS="$LDFLAGS -Wl,-rpath -Wl,-lpthread -Wl,-rpath -Wl,-liomp5"
   else
      if test $useinteluniversalset; then
        LDFLAGS="$LDFLAGS -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_lapack"
        LDFLAGS="$LDFLAGS -lmkl_core -lguide"
      else  
        LDFLAGS="$LDFLAGS -lmkl -lmkl_lapack"
        LDFLAGS="$LDFLAGS -lmkl_intel -Wl,--start-group -lmkl_intel_thread"
        LDFLAGS="$LDFLAGS -lmkl_core -Wl,--end-group -lguide -lpthread -liomp5"
      fi
   fi
fi

# Checks for standard header files.
AC_HEADER_STDC

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
