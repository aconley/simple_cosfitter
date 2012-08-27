;+
;NAME
; contourshow
;PURPOSE
; Makes confidence intervals from the outputs of simple_cosfitter
;USAGE
; contourshow,infile
;INPUTS
; infile         Name of file to read probability from
;OPTIONAL INPUTS
; problevels     Cumulative probability levels. Default:[0.683,0.954,0.997]
; conflevel      The probability the errors are printed at. Def: 0.683
;KEYWORDS
; overplot       Overplot on top of previously existing plot
; transpose      Transpose the axes
; nosumdiff      Do not output sum and difference information for
;                 multi-dimensional probabilities
; ps             Set this if you are writing to a postscript device.
; noxlabel/noylabel Don't show x/y labels
;OPTIONAL OUTPUTS
; chilevels      In the two-2D case, this will return a vector
;                 giving the chisquare offsets above the minimum
;                 value that correspond the the desired confidence
;                 intervals.
; vals           The expectation value of the parameters.  Only works
;                 for the 1D and 2D cases
; lowlimits      The lower limits of the parameters corresponding to
;                 conflevel
; uplimits       The upper limits of the parameters
;RESTRICTIONS
; The 3D case is handled by showing 3 2D projections.
;NOTES
; This code prints error values.  These are the 1D marginalized
; errors.  The quoted estimate is the mean value weighted by the
; probability.  The +- errors are the 68.3% confidence interval
; centered around this best fit value in cumulative probability.
;
; Additional keywords for the plotting routines are handled via
; the _EXTRA mechanism.  Note that linestyles are set by C_LINESTYLE,
; not LINESTYLE, because the CONTOUR routine is used.
;MODIFICATION HISTORY
; Author: Alex Conley
; 2006/01/26 [Michael Wood-Vasey]: Added 3D grid capability
;-

;;Notes
;; Sep 20, 2006 Internal routines renamed to avoid global namespace
;;               issues.

;; Reads the header information from a simple_cosfitter grid file
;; Use this instead of 'readcol' to extract the header information
;;   because 'readcol' counts the number of lines in the file, which
;;   can take a non-negligible amount of time for large grids.
PRO contourshow_readfitshead, file, naxes, axisname, naxis,minaxis,$
                              maxaxis,stepsize
  
  COMPILE_OPT IDL2, STRICTARRSUBS, HIDDEN
  
  header = headfits(file)

  naxes = sxpar(header, 'NAXIS')
  axisname = STRARR(naxes)
  naxis    = INTARR(naxes) 
  minaxis  = FLTARR(naxes)
  maxaxis  = FLTARR(naxes)
  stepsize = FLTARR(naxes)

  FOR i=0, naxes-1 DO BEGIN
     intstr = STRTRIM(STRING(i+1),2) 
     axisname[i] = sxpar(header, 'CTYPE'+intstr)
     naxis[i]    = sxpar(header, 'NAXIS'+intstr)
     minaxis[i]  = sxpar(header, 'CRVAL'+intstr)
     stepsize[i] = sxpar(header, 'CRDELT'+intstr)
     maxaxis[i]  = minaxis[i] + naxis[i]*stepsize[i]
  ENDFOR

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Reads the header information from a simple_cosfitter grid file
;; Use this instead of 'readcol' to extract the header information
;;   because 'readcol' counts the number of lines in the file, which
;;   can take a non-negligible amount of time for large grids.
PRO contourshow_readtexthead, file, naxes, axisname, naxis, minaxis, $
                              maxaxis, stepsize
  
  COMPILE_OPT IDL2, STRICTARRSUBS, HIDDEN

  OPENR,unit,file,/GET_LUN
  tmp = ''
  naxes = 0l
  READF,unit,tmp,naxes,FORMAT='(A6,I)'
  axisname = STRARR(naxes)
  naxis    = LONARR(naxes) ;; I can't imagine truly needing longs here.
  minaxis  = FLTARR(naxes)
  maxaxis  = FLTARR(naxes)
  stepsize = FLTARR(naxes)

  temp = ' '
  delimiter = ' '
  FOR i=0, naxes-1 DO BEGIN   
     ;; Adapted these few lines from 'readcol.pro'
     READF, unit, temp
     temp = strtrim(temp,1)  ;; Remove leading spaces
     
     ;; Split into separate values by our delimiter (' ')
     var = strsplit(strcompress(temp),delimiter,/extract)
     
     axisname[i] = var[0]
     naxis[i]    =  LONG(var[1])
     minaxis[i]  = FLOAT(var[2])
     maxaxis[i]  = FLOAT(var[3])
     stepsize[i] = FLOAT(var[4])
  ENDFOR

  FREE_LUN,unit ;; Close the file

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO contourshow_readgridhead, infile, naxes, axisname, naxis, $
                              minaxis, maxaxis, _EXTRA=_e

  COMPILE_OPT IDL2, STRICTARRSUBS, HIDDEN

  ;; Very simplistic check to see if the fits file ends in '.fits'
  ;; Eventually should make better and check for compressed files, etc.
  IF STRPOS(infile,'.fits') GT -1 THEN BEGIN
     contourshow_readfitshead, infile, naxes, axisname, naxis, minaxis, maxaxis
  ENDIF ELSE BEGIN
     contourshow_readtexthead, infile, naxes, axisname, naxis, minaxis, maxaxis
  ENDELSE

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO contourshow_readfitsdata, infile, prob_ref, naxes, naxis, _EXTRA=_e

  COMPILE_OPT IDL2, STRICTARRSUBS, HIDDEN

  prob_ref = readfits(infile)
  prob_ref /= total(prob_ref) ;;normalize
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO contourshow_readtextdata, infile, prob_ref, naxes, naxis, _EXTRA=_e

  COMPILE_OPT IDL2, STRICTARRSUBS

  rdfloat,infile,prob,SKIPLINE=naxes+1,/SILENT,/DOUBLE

  ;; Impose structure and reform array
  CASE naxes OF
     1: prob_ref = TRANSPOSE( REFORM( prob,                     naxis[0] ) )
     2: prob_ref = TRANSPOSE( REFORM( prob,           naxis[1], naxis[0] ) )
     3: prob_ref = TRANSPOSE( REFORM( prob, naxis[2], naxis[1], naxis[0] ) )
  ENDCASE

  DELVARX,prob
  prob_ref /= total(prob_ref) ;;normalize
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO contourshow_readdata, infile, prob_ref, naxes, axisname, naxis, minaxis, $
                          maxaxis, _EXTRA=_e

  COMPILE_OPT IDL2, STRICTARRSUBS, HIDDEN

  contourshow_readgridhead, infile, naxes, axisname, naxis, minaxis, $
                            maxaxis, _EXTRA=_e
  IF STRPOS(infile,'.fits') GT -1 THEN BEGIN
     contourshow_readfitsdata, infile, prob_ref, naxes, naxis, _EXTRA=_e
  ENDIF ELSE BEGIN
     contourshow_readtextdata, infile, prob_ref, naxes, naxis, _EXTRA=_e
  ENDELSE

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Prints the mean and errors of a 1D distribution
PRO contourshow_error_estimate_1D, axis, probvalues, axisname,$
                                   bestval, lowlimit, uplimit, $
                                   CONFLEVEL=conflevel, PLOT=PLOT, $
                                   COLOR=COLOR, THICK=THICK, QUIET=quiet

  COMPILE_OPT IDL2, STRICTARRSUBS, HIDDEN

  IF N_ELEMENTS( conflevel ) EQ 0 THEN conflevel = 0.682690 ;;1sigma

  ;;bestval is our estimate for the parameter
  bestval = TOTAL(probvalues*axis)/TOTAL(probvalues)

  ;;cumbest is the cumulative probability at that value
  cumprob = TOTAL(probvalues,/CUMULATIVE)
  cumbest = INTERPOL( cumprob, axis, bestval,/SPLINE )

  IF cumbest GT (1.0 - conflevel/2.0) THEN BEGIN
     PRINT,"WARNING: ",axisname," estimate is too close to cumulative prob "
     PRINT," upper limit of 1."
     PRINT," Probability distribution may be too steep"
     PRINT," Perhaps try increasing the number of grid steps?"
     PRINT," Errors may not be accurate"
  ENDIF
  IF cumbest LT (conflevel/2.0) THEN BEGIN
     PRINT,"WARNING: ",axisname," estimate is too close to cumulative prob "
     PRINT," lower limit of 0."
     PRINT," Probability distribution may be too steep"
     PRINT," Perhaps try increasing the number of grid steps?"
     PRINT," Errors may not be accurate"
  ENDIF

  upcum = (cumbest + conflevel/2.0) 
  IF upcum GE 1 THEN BEGIN
     uplimit  = INTERPOL( axis, cumprob, 1 )    
  ENDIF ELSE uplimit  = INTERPOL( axis, cumprob, upcum, /SPLINE )
  botcum = (cumbest - conflevel/2.0) 
  IF upcum GE 1 THEN BEGIN
     lowlimit  = INTERPOL( axis, cumprob, 0 )    
  ENDIF ELSE lowlimit  = INTERPOL( axis, cumprob, botcum, /SPLINE )

  IF ~ KEYWORD_SET( quiet ) THEN $
     PRINT,axisname,": "+STRN(bestval,FORMAT='(F9.5)')+" +"+$
           STRN(uplimit-bestval,FORMAT='(F8.5)')+" -"+$
           STRN(bestval-lowlimit,FORMAT='(F8.5)')+" (mean: "+$
           STRN(0.5*(uplimit-lowlimit),FORMAT='(F8.5)')+ ")"

  IF KEYWORD_SET(PLOT) THEN BEGIN
     scale = 0.2           ; fraction of full axis that max probability extends
     normprob = probvalues/MAX(probvalues)
     normprob *= scale
     normaxis = AXIS/(MAX(axis)-MIN(axis))
     
     IF PLOT EQ 'X' THEN BEGIN
        PLOT, normaxis, normprob, COLOR=COLOR, THICK=THICK, /NORMAL, $
              /NOERASE, XSTYLE=4, YSTYLE=4
     ENDIF
     IF PLOT EQ 'Y' THEN BEGIN
        ;; Complete hack to put a w vs. O_M prob contour on right side
        OPLOT, 1-normprob, axis, COLOR=COLOR, THICK=THICK
     ENDIF
  ENDIF

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO contourshow_plotdoublecontour, prob_ref, axis0, axis1, axisname, $
                                   minaxis, maxaxis, problevels, PS=ps,$
                                   TRANSPOSE=transpose,OVERPLOT=overplot,$
                                   CHILEVELS=chilevels,COLOR=COLOR,$
                                   MARGINALIZE=marginalize, $
                                   XMARGINALIZE=xmarginalize,$
                                   YMARGINALIZE=ymarginalize,VALS=vals,$
                                   LOWLIMITS=lowlimits, UPLIMITS=uplimits, $
                                   CONFLEVEL=conflevel,$
                                   AREA=area,NOSUMDIFF=nosumdiff,$
                                   GIVEAREA=givearea,NOXLABEL=noxlabel,$
                                   NOYLABEL=noylabel,_EXTRA=_e, QUIET=quiet

  COMPILE_OPT IDL2, STRICTARRSUBS, HIDDEN

  IF KEYWORD_SET( ps ) THEN BEGIN
     thick = 4
     charsize = 1.5
  ENDIF ELSE thick = 1
  IF N_ELEMENTS( color ) EQ 0 THEN cprim = !P.COLOR ELSE cprim = color
  def_cprim = !P.COLOR   

;;Find integral levels
  sortprob = SORT( prob_ref )
  totprob = dblarr( N_ELEMENTS( prob_ref ) )
  totprob[0] = prob_ref[ sortprob[0] ]
  FOR i=1L,N_ELEMENTS( prob_ref )-1 DO $
     totprob[i] = totprob[i-1]+prob_ref[ sortprob[i] ]
  totprob = 1.0d0 - totprob
  TABINV, totprob, problevels, ieff
  contlevels = prob_ref[ sortprob[ reverse(ieff) ] ]

  nocontourplot = 0
  noplot = 0

  IF ~ KEYWORD_SET(noplot) THEN BEGIN

     ;; Check to make sure contlevels is increasing
     ;; We just save not plotting the contours because
     ;;   we do still want to print the plot frame
     ;;   to keep our plots in order since the caller
     ;;   is expecting something to be plotted.
     FOR i = 1, N_ELEMENTS(contlevels)-1 DO BEGIN
        IF contlevels[i-1] EQ contlevels[i] THEN BEGIN
           PRINT, "Couldn't find all contour levels within grid"
           PRINT, "Not plotting contours for ", axisname
           nocontourplot = 1
           BREAK
        ENDIF
     ENDFOR
     
     ;;Estimate corresponding chisquare levels
     chilevels = REVERSE( -2*alog( contlevels / max( prob_ref ) ) )

     IF ~ KEYWORD_SET(quiet) THEN PRINT, axisname[0], " vs ", axisname[1]

     ;;Figure out the plot ranges.  These may either be explicitly
     ;; calculated or passed via _EXTRA
     IF N_ELEMENTS(_e) NE 0 && TAG_EXIST(_e,'xrange',/TOP_LEVEL) THEN BEGIN
        xrange = _e.xrange
     ENDIF ELSE BEGIN
        IF KEYWORD_SET( transpose ) THEN $
           xrange = [ MIN(axis1), MAX(axis1) ] ELSE $
              xrange = [ MIN(axis0), MAX(axis0) ]
     ENDELSE

     IF N_ELEMENTS(_e) NE 0 && TAG_EXIST(_e,'yrange',/TOP_LEVEL) THEN BEGIN
        yrange = _e.yrange
     ENDIF ELSE BEGIN
        IF KEYWORD_SET( transpose ) THEN BEGIN
           yrange = [ MIN(axis0), MAX(axis0) ]
        ENDIF ELSE BEGIN
           yrange = [ MIN(axis1), MAX(axis1) ]
        ENDELSE
     ENDELSE
     
     ;;Plot the axes
     IF KEYWORD_SET(transpose) THEN BEGIN
        IF ~ KEYWORD_SET( overplot ) THEN BEGIN
           IF KEYWORD_SET(noxlabel) THEN xtitle='' ELSE $
              xtitle=textoidl(axisname[1])
           IF KEYWORD_SET(noylabel) THEN ytitle='' ELSE $
              ytitle=textoidl(axisname[0])
           PLOT,axis1,axis0,XTITLE=xtitle,YTITLE=ytitle,$
                /NODATA,XTHICK=thick,YTHICK=thick,CHARTHICK=thick,$
                XRANGE=xrange, YRANGE=yrange, $
                COLOR=def_cprim, XSTYLE=1,YSTYLE=1,_EXTRA=_e
        ENDIF
        
        IF ~ nocontourplot THEN BEGIN
           CONTOUR,TRANSPOSE(prob_ref),axis1,axis0,/OVERPLOT,LEVELS=contlevels,$
                   THICK=thick,COLOR=cprim,_EXTRA=_e
           IF KEYWORD_SET(givearea) OR ARG_PRESENT(area) THEN $
              CONTOUR,TRANSPOSE(prob_ref),axis1,axis0,LEVELS=contlevels,$
                      /PATH_DOUBLE, /PATH_DATA_COORDS, /OVERPLOT,$
                      PATH_XY=path_xy, PATH_INFO=path_info
        ENDIF ELSE IF KEYWORD_SET(givearea) OR ARG_PRESENT(area) THEN $
           CONTOUR,TRANSPOSE(prob_ref),axis1,axis0,LEVELS=contlevels,$
                   /PATH_DOUBLE, /PATH_DATA_COORDS,$
                   PATH_XY=path_xy, PATH_INFO=path_info

     ENDIF ELSE BEGIN
        IF ~ KEYWORD_SET( overplot ) THEN BEGIN
           IF KEYWORD_SET(noxlabel) THEN xtitle='' ELSE $
              xtitle=textoidl(axisname[0])
           IF KEYWORD_SET(noylabel) THEN ytitle='' ELSE $
              ytitle=textoidl(axisname[1])
           PLOT,axis0,axis1,XTITLE=xtitle,YTITLE=ytitle,$
                /NODATA,XTHICK=thick,YTHICK=thick,CHARTHICK=thick,$
                XRANGE=xrange, YRANGE=yrange, $
                COLOR=def_cprim, XSTYLE=1,YSTYLE=1,_EXTRA=_e
        ENDIF
        
        IF ~ nocontourplot THEN BEGIN
           CONTOUR,prob_ref,axis0,axis1,/OVERPLOT,LEVELS=contlevels,$
                   THICK=thick,COLOR=cprim,_EXTRA=_e
           IF KEYWORD_SET(givearea) OR ARG_PRESENT(area) THEN $
              CONTOUR,prob_ref,axis0,axis1,LEVELS=contlevels,$
                      /PATH_DOUBLE, /PATH_DATA_COORDS,/OVERPLOT,$
                      PATH_XY=path_xy, PATH_INFO=path_info
        ENDIF ELSE IF KEYWORD_SET(givearea) OR ARG_PRESENT(area) THEN $
           CONTOUR,prob_ref,axis0,axis1,LEVELS=contlevels,$
                   /PATH_DOUBLE, /PATH_DATA_COORDS,$
                   PATH_XY=path_xy, PATH_INFO=path_info
        
     ENDELSE
     
     IF KEYWORD_SET(TRANSPOSE) THEN BEGIN
        IF KEYWORD_SET(xmarginalize) THEN axisplot1="X"
        IF KEYWORD_SET(ymarginalize) THEN axisplot0="Y"
     ENDIF ELSE BEGIN
        IF KEYWORD_SET(xmarginalize) THEN axisplot0="X"
        IF KEYWORD_SET(ymarginalize) THEN axisplot1="Y"
     ENDELSE
     
  ENDIF

;;Estimate 1D errors
;;Plot 1D marginalizations if asked
  marg0 = TOTAL( prob_ref, 2 )
  contourshow_error_estimate_1D,axis0,marg0,axisname[0],PLOT=axisplot0,$
                                COLOR=cprim,THICK=THICK, bestval0, lowlimit0, uplimit0, CONFLEVEL=conflevel,$
                                QUIET=quiet
  marg1 = TOTAL( prob_ref, 1 )
  contourshow_error_estimate_1D,axis1,marg1,axisname[1],PLOT=axisplot1,$
                                COLOR=cprim,THICK=THICK, bestval1, lowlimit1, uplimit1, CONFLEVEL=conflevel,$
                                QUIET=quiet

;;Estimate covariance
  tprob = TOTAL(prob_ref)
  xyprod = (axis0 # axis1 ) * prob_ref / tprob
  cov_x0_x1 = TOTAL(xyprod) - bestval0*bestval1
  rho_x0_x1 = cov_x0_x1 / ( MEAN([bestval0-lowlimit0,uplimit0-bestval0]) * $
                            MEAN([bestval1-lowlimit1,uplimit1-bestval1]) )
  IF ~ KEYWORD_SET(quiet) THEN $
     PRINT,axisname[0]," ",axisname[1], " cov : ",$
           STRING(cov_x0_x1,FORMAT='(F9.6)')," rho: ",$
           STRING(rho_x0_x1,FORMAT='(F8.5)')

  axis0step = (MAX(axis0)-MIN(axis0))/(N_ELEMENTS(axis0)-1)
  axis1step = (MAX(axis1)-MIN(axis1))/(N_ELEMENTS(axis1)-1)

;;Print sum/difference statistics.  We want to generate the
;;sum/difference PDFs
  IF ~ KEYWORD_SET( nosumdiff ) THEN BEGIN
     xproj = axis0 # REPLICATE(1.0,N_ELEMENTS(axis1))
     yproj = REPLICATE(1.0,N_ELEMENTS(axis0)) # axis1
     sum = (xproj + yproj)*prob_ref
     sum /= TOTAL(sum)

     sumdiffbinsize=0.5*(axis0step+axis1step)
     hval = HISTOGRAM( xproj+yproj, BIN=sumdiffbinsize, REVERSE_INDICES=rev,$
                       LOCATIONS=sumloc)
     sumloc += 0.5* sumdiffbinsize ;;Locations at center of bins
     sumprob = FLTARR( N_ELEMENTS(sumloc) )
     FOR i=0, N_ELEMENTS(hval)-1 DO IF rev[i] NE rev[i+1] THEN $
        sumprob[i] = TOTAL( sum[ rev[ rev[i] : rev[i+1]-1 ] ] )
     contourshow_error_estimate_1D,sumloc,sumprob,axisname[0]+" + "+$
                                   axisname[1], sumval, lowlimitsum, $
                                   uplimitsum, CONFLEVEL=conflevel, QUIET=quiet
     
     diff = (xproj - yproj)*prob_ref
     diff /= TOTAL(diff)
     hval = HISTOGRAM( xproj-yproj, BIN=sumdiffbinsize, REVERSE_INDICES=rev,$
                       LOCATIONS=diffloc)
     diffloc += 0.5* sumdiffbinsize ;;Locations at center of bins
     diffprob = FLTARR( N_ELEMENTS(diffloc) )
     FOR i=0, N_ELEMENTS(hval)-1 DO IF rev[i] NE rev[i+1] THEN $
        diffprob[i] = TOTAL( sum[ rev[ rev[i] : rev[i+1]-1 ] ] )
     contourshow_error_estimate_1D,diffloc,diffprob,axisname[0]+$
                                   " - "+axisname[1],diffval, lowlimitdiff, $
                                   uplimitdiff, CONFLEVEL=conflevel,$
                                   QUIET=quiet
  ENDIF

;;Estimate area
;;We assume the grid is uniform
  IF KEYWORD_SET(givearea) OR ARG_PRESENT(area) THEN BEGIN
     cont_idx=N_ELEMENTS(contlevels)-1 ;;This is hte inner contour
     area = POLY_AREA( path_xy[ 0, path_info[cont_idx].offset:path_info[cont_idx].offset + $
                                path_info[cont_idx].n - 1],$
                       path_xy[1,path_info[cont_idx].offset: path_info[cont_idx].offset + $
                               path_info[cont_idx].n - 1])
     ;;Old method
     ;;boxsize = axis0step * axis1step
     ;;wabove = WHERE( prob_ref GE MAX(contlevels),count )
     ;;area = count * boxsize

     IF ~ KEYWORD_SET( quiet ) THEN $
        PRINT,"Area of inner contour: ",STRING( area, FORMAT='(F9.7)')
  ENDIF


  IF ARG_PRESENT(vals) THEN vals = [ bestval0, bestval1 ]
  IF ARG_PRESENT(lowlimits) THEN lowlimits = [ lowlimit0, lowlimit1 ]
  IF ARG_PRESENT(uplimits) THEN uplimits = [ uplimit0, uplimit1 ]

END ;; of plotdoublecontour

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO contourshow_singlecontour,infile,PS=ps,OVERPLOT=overplot,$
                              VALS=vals, LOWLIMITS=lowlimits, $
                              UPLIMITS=uplimits, CONFLEVEL=conflevel,$
                              NOXLABEL=noxlabel,NOYLABEL=noylabel,$
                              QUIET=quiet,_EXTRA=_e

  COMPILE_OPT IDL2, STRICTARRSUBS, HIDDEN

  ;;1D case -- not really a contour, per se.

  thick = 1
  IF KEYWORD_SET(color) THEN cprim = color ELSE cprim = !P.COLOR
  IF KEYWORD_SET( ps ) THEN thick = 3
  def_cprim = !P.COLOR

  contourshow_readdata, infile, prob_ref, naxes, axisname, naxis, $
                        minaxis, maxaxis

;;Make axis
  daxis0 = (maxaxis[0] - minaxis[0])/(naxis[0]-1)
  axis0 = dindgen( naxis[0] ) * daxis0 + minaxis[0]

  yrange = [0,MAX(prob_ref)*1.1]

  IF ~ KEYWORD_SET( overplot ) THEN BEGIN
     IF KEYWORD_SET(noxlabel) THEN xtitle='' ELSE xtitle=textoidl(axisname[0])
     PLOT,axis0,prob_ref,THICK=thick,COLOR=def_cprim,$
          XRANGE=[minaxis[0],maxaxis[0]],YRANGE=yrange,$
          XSTYLE=1,YSTYLE=1,PSYM=10,XTITLE=xtitle,$
          YTITLE='Relative probability',/NODATA,_EXTRA=_e
  ENDIF

  OPLOT,axis0,prob_ref,THICK=thick,COLOR=cprim,_EXTRA=_e

;;Error estimate
  contourshow_error_estimate_1D,axis0,prob_ref,axisname[0], bestval, lowlimit,$
                                uplimit, CONFLEVEL=conflevel,QUIET=quiet
  IF ARG_PRESENT(vals) THEN vals = bestval
  IF ARG_PRESENT(lowlimits) THEN lowlimits = lowlimit
  IF ARG_PRESENT(uplimits) THEN uplimits = uplimit

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO contourshow_doublecontour,infile,problevels,PS=ps,TRANSPOSE=transpose,$
                              OVERPLOT=overplot,CHILEVELS=chilevels,VALS=vals,$
                              LOWLIMITS=lowlimits, UPLIMITS=uplimits, $
                              CONFLEVEL=conflevel, AREA=area, $
                              GIVEAREA=givearea, NOSUMDIFF=nosumdiff, $
                              NOXLABEL=noxlabel,NOYLABEL=noylabel,QUIET=quiet,$
                              _EXTRA=_e

  COMPILE_OPT IDL2, STRICTARRSUBS, HIDDEN

;;2D case

  contourshow_readdata, infile, prob_ref, naxes, axisname, naxis, $
                        minaxis, maxaxis

;; Figure out which ones are non-unitary
  w = WHERE(naxis GT 1, nw)
  IF nw EQ 2 THEN BEGIN
     ;; Just save the axis with more than 1 element
    ;;; We will use this to reform the probability array below
     naxis    = naxis[w]
     axisname = axisname[w]
     minaxis  = minaxis[w]
     maxaxis  = maxaxis[w]
  ENDIF ELSE BEGIN
     MESSAGE,"More than two axis of length greater than 1.  Unable to continue."
     RETURN
  ENDELSE

;;Set up axes arrays
  daxis0 = (maxaxis[0] - minaxis[0])/(naxis[0]-1)
  axis0 = dindgen( naxis[0] ) * daxis0 + minaxis[0]
  daxis1 = (maxaxis[1] - minaxis[1])/(naxis[1]-1)
  axis1 = dindgen( naxis[1] ) * daxis1 + minaxis[1]

  contourshow_plotdoublecontour, prob_ref, axis0, axis1, axisname, $
                                 minaxis, maxaxis, problevels, $
                                 TRANSPOSE=transpose, OVERPLOT=overplot,$
                                 CHILEVELS=chilevels, MARGINALIZE=marginalize,$
                                 XMARGINALIZE=xmarginalize, $
                                 YMARGINALIZE=ymarginalize,VALS=vals, $
                                 LOWLIMITS=lowlimits, UPLIMITS=uplimits,$
                                 CONFLEVEL=conflevel, AREA=area,$
                                 NOSUMDIFF=nosumdiff,GIVEAREA=givearea,$
                                 NOXLABEL=noxlabel,NOYLABEL=noylabel,$
                                 QUIET=quiet,_EXTRA=_e
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; This uses doublecontour to print out 3 2D contours
PRO contourshow_triplecontour,infile,problevels,CLOBBER=CLOBBER,$
                              TRANSPOSE=transpose,OVERPLOT=overplot,$
                              CHILEVELS=chilevels, SHOWW0WA=showw0wa, $
                              CONFLEVEL=conflevel, NOSUMDIFF=nosumdiff,$
                              GIVEAREA=givearea, NOXLABEL=noxlabel, $
                              NOYLABEL=noylabel, QUIET=quiet, _EXTRA=_e

  COMPILE_OPT IDL2, STRICTARRSUBS, HIDDEN

;;3D case

  contourshow_readdata, infile, prob_ref, naxes, axisname, naxis, $
                        minaxis, maxaxis

;; Figure out which ones are non-unitary
  w = WHERE(naxis GT 1, nw)
  IF nw EQ 3 THEN BEGIN
     ;; Just save the axis with more than 1 element
    ;;; We will use this to reform the probability array below
     naxis    = naxis[w]
     axisname = axisname[w]
     minaxis  = minaxis[w]
     maxaxis  = maxaxis[w]
  ENDIF ELSE BEGIN
     errmsg ="More or less than three axis of length greater than 1.  "+$
             "Unable to continue." 
     MESSAGE,errmsg
  ENDELSE



  print, "Reading in: ", infile
  contourshow_readdata, infile, prob_ref, naxes, axisname, naxis, minaxis, $
                        maxaxis

;; Set up the various 2D marginalizations
;; Passing a second argument to TOTAL contracts along that dimension.
;; For reasons eluding comprehension and probably having to do
;;   with IDL's FORTRAN inspiration, the dimensional numbering starts at 1
;;   instead of 0.
  prob_ref_0_1 = TOTAL(prob_ref, 3)
  prob_ref_0_2 = TOTAL(prob_ref, 2)
  prob_ref_1_2 = TOTAL(prob_ref, 1)
  prob_ref_0_1 /= TOTAL(prob_ref_0_1) ;;normalize
  prob_ref_0_2 /= TOTAL(prob_ref_0_2) ;;normalize
  prob_ref_1_2 /= TOTAL(prob_ref_1_2) ;;normalize

;;Set up axes arrays
  daxis0 = (maxaxis[0] - minaxis[0])/(naxis[0]-1)
  axis0 = dindgen( naxis[0] ) * daxis0 + minaxis[0]
  daxis1 = (maxaxis[1] - minaxis[1])/(naxis[1]-1)
  axis1 = dindgen( naxis[1] ) * daxis1 + minaxis[1]
  daxis2 = (maxaxis[2] - minaxis[2])/(naxis[2]-1)
  axis2 = dindgen( naxis[2] ) * daxis2 + minaxis[2]

  IF KEYWORD_SET(showw0wa) THEN BEGIN
     ;; Assume this is Axis1 vs. Axis2
     contourshow_plotdoublecontour, prob_ref_1_2, axis1, axis2, $
                                    axisname[[1,2]], minaxis[[1,2]], $
                                    maxaxis[[1,2]], problevels, $
                                    TRANSPOSE=transpose,OVERPLOT=overplot, $
                                    CHILEVELS=chilevels,$
                                    MARGINALIZE=marginalize, $
                                    XMARGINALIZE=xmarginalize, $
                                    YMARGINALIZE=ymarginalize,$
                                    CONFLEVEL=conflevel,$
                                    NOSUMDIFF=nosumdiff,GIVEAREA=givearea,$
                                    QUIET=quiet,_EXTRA=_e

  ENDIF ELSE BEGIN
     old_pmulti = !P.MULTI
     !P.MULTI = [0,2,2]
     
     ;; 2006/01/26 - Michael Wood-Vasey:
     ;; Oh, this is horrible, but I have to use some COMMON block variables
     ;;   to keep track of the plotting positions so that we can overplot
     ;;   later from a separate call to 'contourshow' if desired.
     COMMON pmulticom, p00, p01, p10, p11, x00, x01, x10, x11, y00, y01, $
        y10, y11

     ;; This ordering of plot commands and axes is chosen so that the
     ;; Y-axes will line up for a row and the X-axes will line up for a column.

     IF KEYWORD_SET(overplot) THEN BEGIN
        ;; Basic check to make sure we actually have at least one of these
        ;;   variables
        IF N_ELEMENTS(p00) GT 0 THEN BEGIN
           !P=p00 & !X=x00 & !Y=y00
        ENDIF ELSE BEGIN
           ;; Remove the overplot since we don't have anything defined to
           ;;  usefully overplot
           OVERPLOT=0
        ENDELSE
     ENDIF
     ;; Axis0 vs. Axis2
     contourshow_plotdoublecontour, prob_ref_1_2, axis1, axis2, $
                                    axisname[[1,2]], minaxis[[1,2]], $
                                    maxaxis[[1,2]], problevels, $
                                    TRANSPOSE=transpose,OVERPLOT=overplot, $
                                    CHILEVELS=chilevels, CONFLEVEL=conflevel,$
                                    GIVEAREA=givearea,QUIET=quiet,_EXTRA=_e

     p00=!P & x00=!X & y00=!Y
     print, "---"
     
     IF KEYWORD_SET(overplot) THEN BEGIN
        !P=p01 & !X=x01 & !Y=y01
     ENDIF
     
     ;; Axis1 vs. Axis2
     contourshow_plotdoublecontour, prob_ref_0_2, axis0, axis2, $
                                    axisname[[0,2]], minaxis[[0,2]], $
                                    maxaxis[[0,2]], problevels, $
                                    TRANSPOSE=transpose,OVERPLOT=overplot,$
                                    CHILEVELS=chilevels, XRANGE=inyrange,$
                                    YRANGE=inzrange,GIVEAREA=givearea,$
                                    CONFLEVEL=conflevel,QUIET=quiet,_EXTRA=_e
     
     p01=!P & x01=!X & y01=!Y
     print, "---"
     
     IF KEYWORD_SET(OVERPLOT) THEN BEGIN
        !P=p10 & !X=x10 & !Y=y10
     ENDIF
     
     ;; Axis0 vs. Axis1
     contourshow_plotdoublecontour, prob_ref_0_1, axis0, axis1, $
                                    axisname[[0,1]], minaxis[[0,1]], $
                                    maxaxis[[0,1]], problevels, $
                                    TRANSPOSE=transpose,OVERPLOT=overplot, $
                                    CHILEVELS=chilevels, CONFLEVEL=conflevel,$
                                    GIVEAREA=givearea, QUIET=quiet, _EXTRA=_e
     
     p10=!P & x10=!X & y10=!Y
     print, "---"
     
     !P.MULTI = old_pmulti
  ENDELSE

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO contourshow,infile,PROBLEVELS=problevels,PS=ps,TRANSPOSE=transpose,$
                OVERPLOT=overplot,CHILEVELS=chilevels,VALS=vals,$
                LOWLIMITS=lowlimits, UPLIMITS=uplimits, CONFLEVEL=conflevel,$
                AREA=area,GIVEAREA=givearea,NOSUMDIFF=nosumdiff,$
                SHOWW0WA=showw0wa , NOXLABEL=noxlabel,$
                NOYLABEL=noylabel, QUIET=quiet,_EXTRA=_e

  COMPILE_OPT IDL2, STRICTARRSUBS

  IF ~ FILE_TEST(infile) THEN BEGIN
     print, "Couldn't find grid file: ", infile
     print, "Giving up."
     RETURN
  ENDIF

  old_pcolor = !P.COLOR
  old_pbackground = !P.BACKGROUND
  !P.COLOR = colordex('black')
  !P.BACKGROUND = colordex('white')

  IF N_ELEMENTS(problevels) EQ 0 THEN $
     problevels = [ 0.682690, 0.954500, 0.9973 ]
  IF N_ELEMENTS( conflevel ) EQ 0 THEN conflevel = 0.682690 ;;1sigma

;;Read in number of axes in file
  contourshow_readgridhead, infile, naxes, axisname, naxis, minaxis, maxaxis

;; Check the sizes of the axes.
;; Figure out which axes ones are non-unitary
  w = WHERE(naxis GT 1, nw)
  CASE nw OF
     1 : contourshow_singlecontour,infile,PS=ps,OVERPLOT=overplot,$
                                   VALS=vals, LOWLIMITS=lowlimits,$
                                   UPLIMITS=uplimits,QUIET=quiet,$
                                   CONFLEVEL=conflevel,NOXLABEL=noxlabel,$
                                   NOYLABEL=noylabel,_EXTRA=_e
     2 : contourshow_doublecontour,infile,problevels,PS=ps,TRANSPOSE=transpose,$
                                   OVERPLOT=overplot,CHILEVELS=chilevels,$
                                   VALS=vals,LOWLIMITS=lowlimits,$
                                   UPLIMITS=uplimits,CONFLEVEL=conflevel, $
                                   AREA=area, NOSUMDIFF=nosumdiff,$
                                   GIVEAREA=givearea,NOXLABEL=noxlabel,$
                                   NOYLABEL=noylabel,QUIET=quiet,$
                                   _EXTRA=_e                                 
     3 : contourshow_triplecontour,infile,problevels,PS=ps,TRANSPOSE=transpose,$
                                   OVERPLOT=overplot,CHILEVELS=chilevels,$
                                   CONFLEVEL=conflevel,NOSUMDIFF=nosumdiff,$
                                   GIVEAREA=givearea,SHOWW0WA=showw0wa,$
                                   NOXLABEL=noxlabel,NOYLABEL=noylabel,$
                                   QUIET=quiet,_EXTRA=_e  
     ELSE : BEGIN
        PRINT,"ERROR in contourshow: Unsupported number of axes: ",naxes
        PRINT,"ERROR in contourshow: Found ", nw, " non-unitary axes."
        RETURN
     END
  ENDCASE
;; Return to initial colours
  !P.COLOR = old_pcolor
  !P.BACKGROUND = old_pbackground

END
