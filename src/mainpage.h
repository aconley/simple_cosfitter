//This isn't really a header file, but the main documentation
// for the Doxygen docs

/*!

\mainpage simple_cosfitter

This is the source documentation for simple_cosfitter.

\section intro Introduction 

This is an implementation of a fairly simple-minded luminosity
distance fitter, intended for use with supernova data.  The
calculational technique is based on evaluating the \f$\chi^2\f$
of the model fit on a grid and marginalization over various nuisance
parameters.  Of course, the nature of these things is that this code
has gotten steadily more complex, so perhaps the simple moniker is
no longer justified.

This code was written by Alex Conley (conley@astro.utoronto.ca).
Michael Wood-Vasey has contributed many improvements to contourshow.
Matrix related code was mostly taken from the 
<A HREF="http://math.nist.gov/tnt/index.html">Template Numerical Toolkit</A>
by Roldan Pozo.

\section installation Installation  

Get the software somehow.  If you have a SNLS account in Toronto,
you can get it from the CVS via something like
\code
cvs -d/data10/cvshome checkout -P simple_cosfitter
aclocal; autoconf; automake --add-missing 
\endcode 
Otherwise (and this is probably easier even for people in Toronto) you 
need a tar file or the like, which the author is happy
to provide.  In fact, there is one currently available 
at http://qold.astro.utoronto.ca/conley/simple_cosfitter/

Then simply run 
\code
./configure --prefix=<where to install>
make all
make install
\endcode
If you desire the ability to convert the probability surfaces to FITS,
and you have cfitsio installed, use --enable-fitswrite.  If your
cfitsio library is someplace unusual, use --with-cfitiolib=<path>
To create this documentation (assuming doxygen is installed on your system)
\code
make docs
\endcode
This will produce an html version of the documentation in the doc/html
subdirectory.  One can also produce a postscript version by
\code
cd doc/latex
make ps
\endcode
after the html documentation has been created.  Personally, I find the
ps version pretty useless.  The documentation assumes Doxygen 1.4 or above.

The code depends on the Gnu Scientific Library being installed and
in the include and library paths.  See http://www.gnu.org/software/gsl/ .
If installing via rpm, make sure that you install the gsl-devel package
as well as the gsl package so that the headers are available.  If you
want the ability to write the probability files as fits, you need
cfitsio.

There is a contour viewing program (contourshow) in the idl subdirectory,
along with some colour management routines.  You might want to add
that subdir to your idl path.

\section usage Usage
The way to run this fitter is via the command line program
runfitter.  This takes one or more arguments, which are the names of
parameter files which specifies how the fit is to be performed,
as well as the name of various input and output options.
\code
runfitter paramfile [paramfile2 ...]
\endcode
Of course, some data will also be necessary.  The location
of the data file is specified in the parameter file, and the
format of this file is detailed in SNeData::readData.  Each
paramfile is handled completely seperately, so providing a set
of several is equivalent to doing each one at a time in sequence
with separate calls.

Examples of parameter and data files can be found in the sample subdirectory.

The basic model for supernova magnitudes used by this code is that
there are two linear corrections allowed, labeled by \f$\alpha\f$ and
\f$\beta\f$.  Why the user is free to use these however they see fit,
the usage envisioned for them is that \f$\alpha\f$ represents the
width-luminosity relation for SNe Ia, while \f$\beta\f$ represents
some sort of colour correction, either for dust or for some sort
of intrinsic colour/magnitude relation.  That is, one can imagine
a corrected magnitude of the form
\f[
  m_{corr} = m_{raw} + \alpha \times \mbox{width} - \beta \times 
   \mbox{colour} .
\f]
So, if one is convinced that any colour variations are due to standard
Milky-Way dust, one could set \f$ \beta = 4.1 \f$ (assuming that
the magnitude was in the rest frame B).  The user has the freedom to
set either of these values (as would be appropriate when fitting to MLCS
data, where the light-curve model itself forces certain values of
\f$\alpha\f$ by training on the low redshift sample, and where all
residual colour terms are assumed due to dust), or fit for either or
both of them, as was done in the SNLS first year cosmology paper.

The contourshow IDL program included in the distribution (in the
idl subdirectory) can be used to view the confidence regions.  Currently 
contourshow can only view 1D-3D  confidence regions.  This code also
outputs error estimates based on the marginalized 1D distributions.
Contourshow works with text and fits based output files (text outputs
are the default from runfitter).  However, there is an option
to produce binary output files.  These are not human readable, but
do take up considerably less space and don't loose precision the way
the text files do (this is not a large effect, generally).

The probability surfaces can be converted to FITS format if the
convert routines were built.  There are two routines for this
purpose.  The first, convert_prob_to_fits, simply takes a
output probability surface file and converts it to a single FITS file.
The second, convert_to_fits, also includes information about the
fit parameters and the SN data used to do the fit.  In general, it
is more useful than the first routine.  However, if you've produced
probability surfaces outside of runfitter (by, say, combine_prob), then
you have to use convert_prob_to_fits.

In order to see the residuals of individual SNe from the best
fit, use VERBOSE YES in the parameter file.  If even more
information about the residuals is desired, use the
EXTENDEDOUTFILE mechanism.

The amount taken for this code to execute, and the amount of
memory used during the run, can literally vary by three orders
of magnitude with small changes to the inputs, so let the user
beware!  

A major complication is whether or not off-diagonal covariance
terms are included between supernovae in the calculation.  This is
a very non-trivial issue.  The model here is that the covariance
matrix is the sum of three terms, one independent of the nuisance
parameters, one proportional to \f$\alpha\f$, and one to \f$\beta\f$
(assuming one is using these parameters).  In addition, the diagonal
elements of the covariance matrix depend on the nuisance parameters
and the covariances between the light-curve fit parameters.  Since
matrix inversion is an expensive operation, the way the code loops
over the variables is completely different in this case.  Using the
covariance matricies makes execution much slower -- see MAGCOVFILE,
WIDTHCOVFILE, COLOURCOVFILE, etc. below for details on how to use this
information.

Note that, unlike previous versions of this code, there is absolutely
no support for analytically marginalizing over \f$\alpha\f$ or
\f$\beta\f$ as was done in Astier et al. (2006).  This is because
this procedure turns out to be intrinsically biased; the analytic
marginalization required that the nuisance parameters be held at fixed
values for the purposes of error propagation (as was done in Perlmutter
et al. (1999)), but this actually produces substantial bias in
\f$\alpha, \beta\f$.  This is somewhat ironic, given that the reason
for introducing this technique in P99 was to avoid biases due to effectively
fitting a line to data with errors in both dimensions.  In this case, the
cure was significantly worse than the disease, and so should not be done.

On the other hand, this version of the code always uses analytic
marginalization over \f$\mathcal{M}\f$, since this can be done without
fixing the errors and is much faster.  As a result of this, there is
no way to output the \f$\mathcal{M}\f$ probability surface.  Sorry.

\section other_constraints Other Constraints

There are also a handful of utility programs which can be used
to apply other cosmological constraints to the data.  In particular,
the Baryon Acoustic Peak measurements of 
<A HREF="http://dx.doi.org/10.1086/466512">Eisenstein '05</A>, the
2dFGRS growth of structure parameter determination, and the 
WMAP 3yr and 5yr shift parameter constraints are available.
There are command line utilities for calculating the probability
surfaces from these measurements: baryonpeak_grid, w2df_grid, wmap_grid,
and wmap5_grid.
The probabilities can then be combined using combine_prob, also
a command line utility.  To use these, you have to provide them
with the name of a parameter file as used for the supernova fits.
So, if sne_output.txt is the output from a SN fit that was run
using the parameter file sne_param.txt,
\code
 baryonpeak_grid sne_param.txt baryon_output.txt
 combine_prob baryon_output.txt sne_output.txt combined_output.txt
\endcode
will produce a combined confidence surface which can be viewed
using contourshow.  Generally speaking, the BAO and WMAP constraints
are much more useful than the 2DF one when combined with supernova
data. 

However, one should be a little careful using these, particularly the
WMAP shift parameters, which is an attempt to encapsulate a lot of
information in a simple way.  The derivation of the shift parameter and
it's errors assumes particular parameter spaces, so if one is considering
parameters not included in this derivation the constraint won't be correct.
The WMAP3 shift parameter doesn't allow for fitting for w in a non-flat
universe, but the WMAP5 one does.  However, this assumes certain
models for $w$ as a function of redshift -- see Komatsu et al. (2008)
for details.  So, caveat emptor!

For more information about these utilities, consult the
help:
\code
 baryonpeak_grid --help
 wmap_grid --help
 w2df_grid --help
\endcode
Again, convert_to_fits doesn't understand or know about these files,
since their names aren't specified in the parameter file.  Therefore,
the routine convert_prob_to_fits is available to make these into
fits files.

\section paramfile Parameter File

The parameter file should be a plain text file which consists
of TAG VALUE(S) sets.  Blank lines are ignored, and lines that
start with # are ignored.

The way to specify that you want to include a parameter in
your fit is to have a line like: PARAMNAME min max n, or,
if a fixed value is desired, PARAMNAME val.  The valid PARAMNAMES
are W, WA, OM, ODE, ALPHA, BETA.  Most of these have synonyms,
noted below.

The following tags are recognized:
<DL>

<DT>DATAFILE filename</DT>
<DD><P>The file containing the input data.  See SNeData::readData
 for a description of the format.</P></DD>

<DT>MAGCOVFILE filename</DT>
<DD><P>The file containing the covariance matrix between supernova
 magnitudes.  This should be a file with the first two entries the size
 of the array followed by the matrix values in C order.
 The diagonal terms in this file are added to the errors resulting
 from the information in DATAFILE, so be careful not to double count
 things -- that is, if you have uncertainties in the peak magnitude
 in DATAFILE, don't also add them here. The code doesn't check to see if this
 matrix is symmetric or diagonal.  If it isn't symmetric strange things
 will probably happen, and if it's diagonal the code will be a lot
 slower for no real reason.  So pay attention!  Further discussion can
 be found in \ref covmatrix .</P></DD>

<DT>WIDTHCOVFILE filename</DT>
<DD><P>The file containing the covariances in the stretch values for
 different SN.  This is multiplied by \f$\alpha^2\f$ when applied.
 Note that if you supply this the code expects you to actually
 be using \f$\alpha\f$ in some way (a fixed value is fine).  Note
 that you can specify a MAGCOVFILE, use \f$\alpha\f$, and not provide
 this with no problem.
 Further discussion can be found in \ref covmatrix .</P></DD>

<DT>COLOURCOVFILE filename</DT>
<DD><P>The file containing the covariances in the colour values for
 different SN, using the same format as for MAGCOVFILE.  This is
 multiplied by \f$\beta^2\f$ when applied.
 Note that if you supply this the code expects you to actually
 be using \f$\beta\f$ in some way (a fixed value is fine). Note
 that you can specify a MAGCOVFILE, use \f$\beta\f$, and not provide
 this with no problem.
 Further discussion can be found in \ref covmatrix .</P></DD>

<DT>MAGWIDTHCOVFILE filename</DT>
<DD><P>The covariances between the magnitudes and widths for different
 SN, using the same format as MAGCOVFILE.  This is multiplied by
 \f$2 \alpha\f$ when applied.  Further discussion can be found in
 \ref covmatrix .</P></DD>

<DT>MAGCOLOURCOVFILE filename</DT>
<DD><P>The covariances between the magnitudes and colours for different
 SN, using the same format as MAGCOVFILE.  This is multiplied by
 \f$- 2 \beta\f$ when applied.  Further discussion can be found in
 \ref covmatrix .</P></DD>

<DT>WIDTHCOLOURCOVFILE filename</DT>
<DD><P>The covariances between the widths and colours for different
 SN, using the same format as MAGCOVFILE.  This is multiplied by
 \f$- 2 \alpha \beta\f$ when applied.  Further discussion can be found in
 \ref covmatrix .</P></DD>

<DT>OUTFILE filename</DT>
<DD><P>File to output marginalized probability array to.  See
  cosgrid3D::writeFile for an example of the type of format.</P></DD>

<DT>W minw maxw nw OR W fixval</DT>
<DD><P>The \f$w_0\f$ specification. Unless this is set
 w is assumed to be -1 (which corresponds to the cosmological constant).
 W0 and W_0 are supported as synonyms.</P></DD>

<DT>WA minw maxw nw OR WA fixval</DT>
<DD><P>The \f$w_a\f$ specification. Unless this is set
 w is assumed to be 0.0 (which corresponds to constant \f$w\f$).
 The model is \f$w\left( a \right) = w_0 + w_a \left( 1 - a \right)\f$,
 where \f$a = \frac{1}{1+z}\f$ is the scale parameter.
 W_A is a synonym.</P></DD>

<DT>OM minom maxom nom OR OM fixval</DT>
<DD><P>The \f$\Omega_m\f$ specification. Supported synonyms
 are OMEGAM and OMEGA_M.</P></DD>

<DT>OL minol maxol nol</DT>
<DD><P>The \f$\Omega_{DE}\f$ specification.  Unless w is set,
  this is \f$\Omega_{\Lambda}\f$.  Supported synonyms are
  ODE, OMEGA_LAMBDA, OMEGALAMBDA, OMEGA_DE, OMEGADE.</P></DD>

<DT>ALPHA minal maxal nal OR ALPHA fixval</DT>
<DD><P>The \f$\alpha\f$ specification.  By including a line like
 this you are choosing to include the \f$\alpha\f$ correction.  This 
 parameter is the slope of the width-luminosity relationship -- that is, 
 it applies the stretch correction.  AL is a synonym.</P></DD>

<DT>BETA minal maxal nal OR BETA fixval</DT>
<DD><P>The \f$\beta\f$ specification.  By including a line like
 this you are choosing to include the \f$\beta\f$ correction.  This 
 parameter is the slope of the colour-luminosity relationship -- that is, 
 it applies the extinction/intrinsic colour correction.  
 BE is a synonym.</P></DD>

<DT>FLATONLYFIT YES/NO</DT>
<DD><P>If this is set then the Universe is assumed to be flat.  In this
 case the OL gridding is ignored, and \f$\Omega_{DE} = 1 - \Omega_m\f$.
 This is equivalent to setting FIXCURV 0.0
 </DD></P>

<DT>FIXCURV value</DT>
<DD><P>If this is set then the Universe is assumed to have a fixed
 value of \f$\Omega_{curv}\f$ set by the following value.  In this
 case, any OL gridding is ignored, and 
 \f$\Omega_{DE} = 1 - \Omega_{curv} - \Omega_m\f$. FIXCURVE is also
 accepted, in case you are attached to your e's.
 </DD></P>

<DT>USEKOMATSUFORM YES/NO</DT>
<DD><P>If this is set then the Komatsu et al. (2008) form is used
  for \f$w\left(a\right)\f$ rather than the usual Linder form.
  This is off by default, and only matters if you fit for \f$w_a\f$.
 </DD></P>

<DT>ATRANS value</DT>
<DD><P>The transition scale factor for the Komatsu form of the
 variable dark energy equation of state.  The default is 0.1.
 Has no effect unless USEKOMATSUFORM YES.
 </DD></P>

<DT>SNDISP dispersion OR SNDISP dataset dispersion</DT>
<DD><P>The amount of intrinsic scatter to be added in quadrature to the
 magnitudes of each SN.  The default value is zero.  If you use
 the second form (with dataset), then different dispersions are
 applied to each dataset as specified in the data file.  Note
 that you must provide an actual value for each dataset if you use
 this feature.</P></DD>

<DT>PECZ peculiar_v</DT>
<DD><P>The peculiar velocity error in redshift units.  Will be added as
 an error component to all SN, although it only really has an effect
 in the nearby sample.  The default value is 0.001, which is equivalent
 to 300 km/sec.</P></DD>

<DT>ALPHABETAOUTFILE filename</DT>
<DD><P>If this is set, then the confidence contours for
 \f$\alpha\f$ and \f$\beta\f$ are output.  If one of them is
 not used (or is fixed) then the 1D distribution over the other
 is output. ALBETAOUTFILE is a synonym.</P></DD>

<DT>EXTENDEDOUTFILE filename</DT>
<DD><P>If this is set, then extended information about each SN
 and its residual from the best fit is output to the specified
 filename.  This is similar to what happens when verbose is set,
 but is more detailed.  Note that if analytical \f$ {\mathcal M} \f$
 marginalization is used, then its value will be estimated for
 these purposes, even though it is not explicitly used in the fit.
 In particular, users should be careful when comparing the
 \f$ \chi^2 \f$ output to this file to that which the fit reports --
 the latter represents an equivalent \f$ \chi^2 \f$ after marginalizing
 over this parameter, the other represents it at a particular value.
 EXTENDEDOUT and EXTOUTFILE are supported as synonyms.</P></DD>

<DT>PARAMSUMMARYFILE filename</DT>
<DD><P>If this is set, then a summary of the fitted parameters with
approximate error estimates will be output to this file.  This
does not provide error estimates for \f$\mathcal{M}\f$.  Countourshow
is probably more accurate.</P></DD>

<DT>BINARYOUTFILE YES/NO</DT>
<DD><P>If this is set, then all of the probability surfaces are written
 as binary files instead of as human readable text.  These take up less
 space, but contourshow doesn't understand them.  The text files also
 have a potential loss of precision which this avoids.  Only really
 interesting if one later plans on converting the probability to some
 other format.</P></DD>

<DT>VERBOSE YES/NO</DT>
<DD><P>If this is set to yes, then more information is output as the
 fit runs.  The default is no.  This will turn on the progress
 bar.</P></DD>

<DT>SHOWPROGBAR YES/NO</DT>
<DD><P>If this is set to yes, then a progress bar is shown
 for some (but not all) types of fits.  The default is no.
 </P></DD>

<DT>DZINT dz_step</DT>
<DD><P>If you have a large number of SN, then the luminosity distances
 will be calculated using a non-adaptive trapezoidal rule.  This is the
 step size.  The default is 0.005.  See the documentation for lumdist
 for more information.</P></DD>

</DL>

\section examples Examples

There are some examples included in the sample subdirectory
demonstrating some of the capabilities of this code,
which focus on the data samples of Knop '03 and Astier '06.  
To keep things simple, the default of human readable text output 
is used for the examples. None of the sample files use a full covariance 
matrix.

The lowe_param.txt file is a parameter file for a fit to the low-extinction
subset of Knop'03 without extinction correction.  This uses
ALPHABETAOUTFILE to output the \f$\alpha\f$ distribution.
On my machine (3.2GHz PIV) this takes about 1 second to run.  To run the fit, 
\code
cd sample
../bin/runfitter lowe_param.txt
\endcode
The \f$\Omega_m\f$, \f$\Omega_{\Lambda}\f$
confidence regions can be seen by (from IDL, and assuming that
contourshow is in your IDL path)
\code
IDL> contourshow,'lowe_output.txt'
\endcode
and the \f$\alpha\f$ region by
\code
IDL> contourshow,'lowe_alpha.txt'
\endcode
You can see some more information about the results of the fit
in the lowe_results.txt file (which was specified by the
PARAMSUMMARY line in lowe_param.txt).  There is also some SN-by-SN
information in lowe_extout.txt (EXTENDEDOUT in the param file).

If the FITS conversion routines were built, then the output can be
converted into a handy FITS file via
\code
 ../bin/convert_to_fits lowe_param.txt lowe_output.fits
\endcode
Note that contourshow doesn't understand this format (yet).
This output file has three HDUs.  The first is the probability surface,
the second is a binary table of SN information, and the third is
the marginalized \f$\alpha\f$, \f${\mathcal M}\f$ surface.
They can be viewed in most fits readers -- i.e.,
\code
 ds9 lowe_output.fits
 ds9 lowe_output.fits[ALPHPROB]
\endcode
will show you the cosmology probability and the marginalized
\f$\alpha\f$ surface.  Note that the WCS
information has been loaded with information about the axes of the
probability surface.

The lowe_noscriptm_param.txt file is similar, but only analytically
marginalizes over \f$\mathcal{M}\f$, and loops over \f$\alpha\f$.
Therefore one can not output the confidence regions
for \f$\alpha\f$ and \f${\mathcal M}\f$.

The lowe_extcorr_param.txt is a fit to the same data set but
with an extinction correction based on the standard \f${\mathcal R}_B=4.1\f$
value.  This is set by the line
\code
BETA 4.1
\endcode
in the parameter file.

The lowe_extcorr_wom_param.txt file is another fit to the Knop'03 data
with an extinction correction using \f${\mathcal R}_B=4.1\f$,
but the cosmological variables are \f$w\f$ and \f$\Omega_m\f$, with
a flat universe assumed.  This is done by removing the OL line
from the parameter file and including
\code
W -5 -0.3 181
FLATONLYFIT YES
\endcode
The range on \f$\Omega_m\f$ has also been modified to get better sampling
over the region of interest. To view the results:
\code
IDL> contourshow,'lowe_extcorr_wom_output.txt'
\endcode
Note that the probability contours do not close off at the bottom
over the calculated range, and so the confidence limits are not
strictly accurate.  If this data is combined with other measurements,
this problem does not occur.  You can try fixing \f$\Omega_{curv}\f$ at
a slightly different value by something like
\code
FIXCURV 0.02
\endcode
instead of FLATONLYFIT.  The difference in the results is not large.

Here's how to add the baryon acoustic oscillations constraint to
this data.
\code
bin/baryonpeak_grid lowe_extcorr_wom_param.txt lowe_extcorr_wom_bap_output.txt
bin/combine_prob lowe_extcorr_wom_output.txt lowe_extcorr_wom_bap_output.txt lowe_extcorr_wom_combined.txt
\endcode
You can then view the last file in contourshow as above.  To overplot
the constraints nicely
\code
IDL> contourshow,'lowe_extcorr_wom_output.txt'
IDL> contourshow,'lowe_extcorr_wom_bap_output.txt',/OVERPLOT, $
     C_LINESTYLE=2,COLOR=colordex('blue')
IDL> contourshow,'lowe_extcorr_wom_combined.txt',/OVERPLOT,$
     COLOR=colordex('red')
\endcode

Some of the other example files included in this directory include
fits to the Astier et al. (2006) data to \f$\Omega_m\f$ in flat, LCDM
univers,  \f$\Omega_m, \Omega_{\Lambda}\f$, and to a flat
universe \f$\Omega_m, w\f$.  Note that the results you get aren't
quite the same as those given in the paper.  This is because A06
analytically marginalized \f$\alpha, \beta\f$ by fixing them for the
purposes of error propagation, which is now known to be a bad idea.
Fortunately, the effects were small.

\section algorithm Algorithms

The basic approach is to form a grid over the cosmological parameters
and calculate the \f$ \chi^2 \f$ at each point from the expression
\f[
 \chi^2 = \sum_{i=1}^{N} \frac{ \left( m_{i} - m_{theory} \left( 
   z_i , w, \Omega_m, \Omega_{\Lambda}, \alpha, \beta, {\mathcal M} \right) 
 \right)^2 }{ \sigma^2_i }
\f]
where the sum is over the SN, \f$ m_i \f$ is the observed rest frame
magnitude of each SN, and \f$ m_{theory} \f$ is given by the following
expression
\f[
 m_{theory} = 5 \log_{10} {\mathcal D}_L
   \left( z_i , w, \Omega_m, \Omega_{\Lambda} \right)
   - \alpha \left( s_i - 1 \right) + \beta {\mathcal C}_i + {\mathcal M},
\f]
where \f$ s_i \f$ is the stretch and \f$ {\mathcal C}_i \f$ is some
colour parameter.  \f$ {\mathcal D}_L \f$ is the \f$ c/H_0 \f$ reduced
luminosity distance (this factor is effectively absorbed into
\f$ {\mathcal M} \f$).  More explicitly,
\f[
 {\mathcal D}_L = \frac{1+z_{hel}}{
    \sqrt{\left| 1 - \Omega_m - \Omega_{\Lambda} \right| } }
 S_{k} \left[ \sqrt{ \left| 1 - \Omega_m - \Omega_{\Lambda} \right| }
   \int_0^{z_{cmb}} \, \frac{ dz } {
   \sqrt{ \left(1+z\right)^3 \Omega_m + \left(1+z\right)^{3\left(1+w\right)} 
   \Omega_{\Lambda}
   + \left( 1 + z\right)^2 \left(1 - \Omega_m - \Omega_{\Lambda}\right) } }
 \right],
\f]
where \f$ S_k \f$ is the identity for a spatially flat universe, sin for a
spatially closed one, and sinh for a spatially open universe.  Properly
speaking, \f$\Omega_{\Lambda}\f$ refers to the case \f$ w=-1 \f$, but 
I ignore this question of nomenclature.  This expression is only valid
for a constant dark energy equation of state; the code also supports
the \f$ w\left( a \right) = w_0 + w_a \left(1 - a\right)\f$
parameterization.
 
Not all of these terms are always used.  In fact, trying to fully 
grid over all six parameters would be a serious computing challenge. 

Once the \f$ \chi^2 \f$ is computed at every point, it is converted into
a probability via \f$ P \propto \exp \left( -\frac{1}{2} \chi^2 \right)\f$.
The proportionality is taken care of by normalizing over the grid, so the
errors become tricky to interpret if the region of high probability
is not enclosed by the grid considered.  Of course, the problem of
interpreting errors near a boundry is always difficult. 
The undesired (nuisance) parameters 
( \f$ \alpha, \beta, {\mathcal M} \f$ ) are then marginalized over  
and the resulting probability surface is output.  Note that this  
is therefore a somewhat Bayesian approach, as it requires priors 
over these parameters.  In this code flat priors are assumed. 
 
The definition of \f${\mathcal M}\f$ used in this code differs slighly
from that used by Perlmutter et al. 1999 and Knop et al. 2003 in that
the extra constant factor of c has also been absorbed.  Explicitly, 
\f[
 {\mathcal M} = M - 5 \log_{10} h + 42.38410 ,
\f] 
where \f$M\f$ is
the absolute magnitude of a fiducial SN Ia in whatever sense you
are using as a standard candle (the last term is \f$ 15.0 + 5 \log_{10} c \f$,
where c is in kilometers per second).

To convert the probability surface to confidence intervals, the 
cumulative probability distribution is constructed and limits are
found such that 68.3%, etc. of the probability is enclosed.  This 
is necessary because the fitting problem is sufficiently non-linear
in the cosmological parameters that the \f$\chi^2_{min} + 1\f$ (or
the multidimensional equivalent (\f$ \chi^2_{min}+2.3\f$ in 2D)) 
method of estimating the errors is no longer entirely trustworthy.

\section equationofstate The dark energy equation of state

There are two options for the dark energy equation of state.
The first is the standard Linder (2003) form:
\f[
w \left( a \right) = w_{0} + \left( 1 - a \right) w_a .
\f]
This has the disadvantage that it is not well behaved for large
redshifts, so can cause serious problems with CMB constraints.

The other form is that proposed by Komatsu et al. (2008), and
is like the Linder form but asymptotes to \f$w=-1\f$ at high redshift.
\f[
w \left( a \right) = \frac{ a \left( w_{0} + \left( 1-a \right) w_{a} \right)}
                          { a + a_{t} } - \frac{ a_t }{a + a_t}
\f]
which introduces another parameter, \f$a_t\f$.  This relates to
the transition between the asymptotic forms.  Note that this reduces
to the Linder form as \f$ a_{t} \to 0 \f$.  The code does not support
fitting \f$a_{t}\f$, as it will be fairly poorly constrained.  Instead,
it is held fixed at a value the user can specify using ATRANS.

\section covmatrix Off diagonal terms in the covariance matrix

Life is much simpler if there are no covariances between different
SNe.  The basic problem is that, generally speaking, the elements
of the covariance matrix depend on the values of the nuisance parameters
\f$\alpha\f$ and \f$\beta\f$.  This means that for each value of
these parameters we have to re-invert the covariance matrix to get the
\f$\chi^2\f$.  Inverting matricies is very, very expensive -- a 
\f$N^3\f$ process, where N is the number of SNe.  This is much more
expensive than computing luminosity distances. However, it turns out that
the inversion is not the dominant cost, but rather the cost of evaluating
the \f$\chi^2\f$ and including off-diagonal terms.  In any case, using
this option slows the code down considerably, and deciding whether or 
not to use a covariance matrix is a serious issue.  The various pieces are 
specified via the MAGCOVFILE, WIDTHCOVFILE, COLOURCOVFILE, MAGWIDTHCOVFILE,
MAGCOLOURCOVFILE, and WIDTHCOLOURCOVFILE lines in the parameter file.

The model for the combined covariance matrix is
\f[
\mathbf{V} = \mbox{diag} + \mathbf{V}_{mm} + \alpha^2 \mathbf{V}_{ww} + 
\beta^2 \mathbf{V}_{cc} + 2 \alpha \mathbf{V}_{mw} 
- 2 \beta \mathbf{V}_{mc} - 2 \alpha \beta \mathbf{V}_{wc}
\f]
where \f$V_{mm}\f$ is the part specified by MAGCOVFILE, \f$\mathbf{V}_{ww}\f$
by WIDTHCOVFILE, \f$\mathbf{V}_{cc}\f$ by COLOURCOVFILE, 
\f$\mathbf{V}_{mw}\f$ by MAGWIDTHCOVFILE, \f$\mathbf{V}_{mc}\f$ by 
MAGCOLOURCOVFILE, and \f$\mathbf{V}_{wc}\f$ by WIDTHCOLOURCOVFILE.  diag is
the diagonal elements corresponding to the \f$\sigma_i\f$ in our
formula for \f$\chi^2\f$ above.

Note that the diagonal entries in the input covaraince matrix are
considered additional information to that supplied in DATAFILE,
and are added to those values.  So the final diagonal term in the
covariance matrix is whatever was in the matrix files plus the magnitude
errors from DATAFILE, incorporating the stretch and colour correction,
the intrinsic dispersion, and the peculiar velocity errors.  So don't
put that information in the MAGCOVFILE.  Appropriate things to put there
are errors that you expect to be highly correlated.  A good example
are zero-point errors, which should be put in the covariance matrix
but not DATAFILE because they will affect almost all of the sample
identically.

The covariance matrix will by symmetric positive definite, so the
method used for inverting it is to calculate the Cholesky decomposition
\f$\left( \mathbf{V} = \mathbf{L} \cdot \mathbf{L}^T \right)\f$, 
where L is lower diagonal, invert L (which is relatively cheap for lower 
triangular matricies), and then compute \f$\mathbf{V}^{-1} = 
\left(\mathbf{L}^{-1}\right)^T \cdot \mathbf{L}^{-1}\f$.
This is a factor of two or so faster than other methods, unless your
matrix is very small (in that case, a direct application of Cramer's rule may
be faster) or it has a special form like Sherman-Morrison-Woodbury.  
The TNT implementation of Cholesky decomposition is
used, which is not nearly as fast as that of LAPACK.  However, support
is provided, which must be compiled in, for the ATLAS and Intel MKL 
libraries.  See configure --help for information.  If you are going to
be using the covariance matricies, compiling in this support is probably
well worth your while.  There may also be some numerical accuracy
benefits, since the by-hand code isn't really designed to worry about
roudning and all that.

Aficionados of matrix manipulations may point out that when one
is computing a \f$\chi^2\f$ using
\f[
\chi^2 = \vec{x}^T \cdot \mathbf{V}^{-1} \cdot \vec{x}
\f]
there is no need to actually compute the inverse.  Instead one can
decompose the matrix by LU decomposition, so that \f$V = \mathbf{L} 
\cdot \mathbf{U}\f$, where L is lower diagonal and U is upper diagonal, 
then solve the two equations \f$ \mathbf{L} \cdot \vec{y} = \vec{x} , 
\mathbf{U} \cdot \vec{z} = \vec{y} \f$
and finally calculate \f$ \chi^2 = \vec{x}^T \cdot \vec{z} \f$ (because
of memory usage patterns, the LU decomposition is typically faster
than the Cholesky one here).  This method is, in fact, far faster 
than finding the inverse for a single difference vector \f$\vec{x}\f$.  
However, it has the disadvantage that the
solver must be rerun for every different value of \f$\vec{x}\f$.  In
our case, we will be computing many different values of \f$\vec{x}\f$
for a covariance matrix which changes less frequently, or not at all,
and therefore computing the inverse and then using matrix multiply wins out.
At least, that was what I concluded after running a whole bunch of tests.
One also needs the sum of the entries of the inverse covariance matrix
to compute the analytic marginalization over \f$\mathcal{M}\f$ --
see below.

\section analyticmarg Analytically marginalizing over the nuisance parameters

One can analytically marginalize over any linear nuisance parameters
given a suitable simple prior.  This is always possible for \f${\mathcal M}\f$ 
because it is simply an additive parameter in \f$m_{theory}\f$ and because
the errors for each SN do not depend on its value.  The same trick can 
be employed to analytically marginalize over
\f$\alpha\f$ and \f$\beta\f$.  However, the situation is more complex
because the errors for each data point depend on their values.
In order to make the integral manageable, it is necessary to fix
the errors in the computation.  Unfortunately, this results in biased
estimates for these parameters, so is not supported in this version of
the code.

The basic method has been described in several places, but a useful
reference is <A HREF="http://www.edpsciences.org/articles/aa/full/2001/46/aah2780/aah2780.html?access=ok">Goliath et al., 2001 A&A 380, p. 6-18</A>.

Briefly, the \f$\chi^2\f$ can be analytically marginalized by converting
it to a probability, integrating, and then re-converting to \f$\chi^2\f$.
\f[
 \chi^2_{{\mathcal M} \, marg} = -2 \log \left[ 
   \int^{\infty}_{-\infty} \, d{\mathcal M} \, 
   \exp \left( -\frac{1}{2} \chi^2 \right)
      \pi \left( {\mathcal M} \right) 
 \right],
\f]
where \f$\pi \left( {\mathcal M} \right)\f$ is the prior over this parameter.
In this code flat priors are assumed.

If one defines \f$\Delta \vec{ m } = \vec{m} - \mathcal{D}_L \left(
\vec{z}; w, \Omega_{m}, \Omega_{DE}\right) - 
\alpha \left(\vec{s} -\vec{1} \right)+ \beta \vec{c}\f$
then one can write the \f$\chi^2\f$ as
\f[
  \chi^2 = \left( \Delta \vec{m} - \mathcal{M} \vec{1} \right)^T \cdot 
  \mathbf{V}^{-1} \cdot \left( \Delta \vec{m} - \mathcal{M} \vec{1} \right)
\f]
where \f$\mathbf{V}^{-1}\f$ is the inverse covariance matrix.
One then expands the products and does the integral, ending
up with the following terms:
\f[
 A  \equiv \Delta \vec{ m }^{T} \cdot \mathbf{V}^{-1} \cdot \Delta \vec{m} 
\f]
\f[
 B  \equiv \Delta \vec{ m }^{T} \cdot \mathbf{V}^{-1} \cdot \vec{1} 
\f]
\f[
 E  \equiv \vec{1}^T \cdot \mathbf{V}^{-1} \cdot \vec{1} 
\f]
It might be helpful to show what these parameters look like in
the case of diagonal errors:
\f[
A \equiv \sum_i \frac{ \Delta m^2_i }{ \sigma^2_i }
\f]
\f[
B \equiv \sum_i \frac{ \Delta m_i }{ \sigma^2_i }
\f]
\f[
E \equiv \sum_i \frac{ 1 }{ \sigma^2_i }
\f]



The marginalization over \f$\mathcal{M}\f$ is then
\f[
 \chi^2_{\mathcal{M}} = A + \log \frac{E}{2 \pi} - \frac{B^2}{E} .
\f]

The next question that comes to mind is how one gets estimates
for the value of \f$\mathcal{M}\f$.  The solution is actually
quite simple:  \f$\chi^2_{\mathcal{M}}\f$ is identical, to within a 
constant offset, to the value of \f$\chi^2\f$ one gets with
\f$\mathcal{M} = \frac{B}{E}\f$.  

The missing term is worth discussing in more detail.  It is common
practice in the cosmological fitting industry to compute the
\f$\chi^2\f$ by simply minimizing the nuisance parameters for fixed
values of the cosmological parameters and then evaluating the \f$\chi^2\f$
at that point - for example, this is how \f$\Omega_b h^2\f$ and \f$h\f$
are handled for the WMAP 5 yr shift parameters.
As discussed in Cash (1976) and Lampton et al. (1976), this procedure
usually does what you want, making it unecessary to explicitly grid
over your nuisance parameters.  Does this work here?  Unfortunately
the answer is no, and so we have to explicitly grid and marginalize,
unless we are lucky enough to be able to do it analytically as for 
\f$\mathcal{M}\f$.  

To see how this works, consider the above marginalization over 
\f$\mathcal{M}\f$ and imagine that we want to apply the minimize and
evaluate procedure.  Expanding the \f$\chi^2\f$ we get
\f[
  \chi^2 = A - 2 \mathcal{M} B + \mathcal{M}^2 E .
\f]
The minimum value occurs for \f$\mathcal{M} = B / E\f$.   If we insert
this, we find that this is the same as the analytic marginalization
except that the \f$\log E / 2 \pi\f$ term is missing.  The problem is this:
the missing term depends on \f$\alpha\f$ and \f$\beta\f$, and so if we
don't include this term, when we compare different values of
these parameters the normalization is different.  This means we can't 
compare the \f$\chi^2\f$ for different values of these nuisance parameters,
which makes it impossible to evaluate the fit.

Now, in this particular case, we are all right because we were able
to analytically figure out what the missing term was, and so we can
just stick it back in if we want.  But, if we try the same trick with
\f$\alpha\f$ and \f$\beta\f$, we don't know what the missing term is
because the integral is just too hard to do (this is particularly difficult
when the covariance matrix is included, because then one has to integrate
over the analytic inverse of a N by N matrix, which is hopeless for
any reasonable number of SN). Since the best fitting \f$\alpha\f$ and
\f$\beta\f$ inevitably depends on the cosmological parameters, if we
try this approach, we sacrifice the ability to compare the \f$\chi^2\f$
for different parameters, which is obviously disastrous.  

Looking deeper, the real issue is that our errors depend on \f$\alpha\f$
and \f$\beta\f$ -- if this were not the case, these issues would go away.
Also note that this doesn't apply to a purely frequentist, best fit
approach, but only comes in when one tries to marginalize over nuisance
parameters, which is fundamentally Bayesian.  Also note that if
\f$\alpha\f$ and \f$\beta\f$ are fixed, then this issue goes away again.
However, they truly have to be fixed for some external reason - simply
holding them fixed for the purposes of computing the errors while 
simultaneously fitting them as magnitude corrections, as was done
in Astier et al. (2006) or Perlmutter et al. (1999), is actually quite
biasing and should be avoided.

\section pecvel Handling peculiar velocities

The code as written is designed to properly handle the distinction
between the observer frame, heliocentric redshift, and the redshift
corrected for the velocity of the Earth relative to the cosmic
rest frame (the famous CMB dipole).  This is why we distinguish
between \f$z_{hel}\f$ and \f$z_{cmb}\f$ in the
input files.  However, if one also wants to correct for the peculiar
velocities of the SN hosts with respect to the CMB frame, then an
additional correction is necessary which can not simply be effected
by modifying these redshifts.  Here we follow the notation and
results of Hui and Greene (2006 Phys Rev D 73, id.\ 123526).

What we want to compare with observations is \f$d_L\f$, the luminosity
distance in the perturbed universe.  This is different than the
homogenous luminosity distance, \f$\bar{d}_L\f$.  It is this quantity
which is the familiar one calculated from the cosmological parameters.
If we restrict ourselves to the flat, cosmological constant case 
for simplicity, this is defined by
\f[
 \bar{d}_L \left( z \right) = \left(1 + z\right) \chi_e \left( z \right)
\f]
where
\f[
 \chi_e \left( z \right) = \int_0^z\, dz / \sqrt{ \Omega_m \left(1+z\right)^3
  + \Omega_k \left( 1 + z \right)^2 + \Omega_{\Lambda} } .
\f]
The effects of peculiar velocity are two-fold.  First, if we
denote the peculiar motion of the observer (i.e., the Earth)
with respect to the CMB by \f$\vec{v}_0\f$, that of
the SN by \f$\vec{v}_e\f$, and the unit vector between the observer
and the object by \f$\vec{n}\f$, then the relation between the observed
redshift \f$z_{hel}\f$ and the redshift in the absence of 
perturbations \f$\bar{z}\f$ (i.e., the redshift in the CMB frame) is
\f[
 1 + z_{hel} = \left( 1 + \bar{z} \right) \left( 1 +
  \vec{v}_e \cdot \vec{n} - \vec{v}_0 \cdot \vec{n} \right) 
\f]
where \f$c=1\f$.  So the first effect is to shift the redshift.
The second effect is more subtle.  As shown by Hui \& Greene,
the relation between the homogenous universe luminosity distance
evaluated at the redshift as it would be measured in the CMB
rest frame and the actual luminosity distance that we want to
compare with our magnitudes is
\f[
  d_L = \bar{d}_L \left(1 + 2 \vec{v}_e \cdot \vec{n} - 
   \vec{v}_0 \cdot \vec{n} \right) . 
\f]
It is the factor of 2 in front of \f$\vec{v}_e\f$ that is
perhaps surpising.

What simple_cosfitter actually calculates is \f$\tilde{d}_L\f$,
which is a function that takes two inputs and is defined by
\f$ \tilde{d}_L \left( z_1, z_2 \right) = \left(1 + z_1 \right)
 \chi_e \left( z_2 \right)\f$.  If we plug \f$z_{hel}\f$
and \f$\bar{z}\f$ into these two slots, then the relation between
\f$\tilde{d}_L\f$ and the desired \f$d_L\f$ is
\f[
  d_L = \tilde{d}_L \left( 1 + \vec{v}_e \cdot \vec{n} \right) .
\f]
Thus, in the absence of \f$\vec{v}_e\f$, this is the correct
thing to do.  However, if the SN peculiar velocity is not zero, 
a further adjustment is necessary.  The easiest approach is to
adjust \f$m_B\f$, since the only thing we care about is
\f$m_B - 5 \log_{10} d_L\f$.

To be more explicit, given \f$z_{hel}\f$, \f$\vec{v}_0 \cdot \vec{n}\f$,
and \f$\vec{v}_0 \cdot \vec{n}\f$, one should first determine
\f$\bar{z}\f$ as above.  \f$z_{hel}\f$ and \f$\bar{z}\f$
should then be used in the data file for simple_cosfitter as \f$z_{hel}\f$
and \f$z_{cmb}\f$.  Then the magnitude values should be
adjusted by \f$-5 \log_{10} \left( 1 + \vec{v}_e \cdot \vec{n}
\right)\f$.

*/

