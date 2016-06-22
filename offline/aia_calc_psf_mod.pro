FUNCTION aia_getmeshinfo, wavelength, $
                          use_preflightcore=use_preflightcore, $
                          dwavelength = dwavelength, $
                          qabort=qabort

; =========================================================================
;+
; PROJECT:
;
;       SDO / AIA
;
; NAME:
;
;       AIA_GETMESHINFO
;
; CATEGORY:
;
;       Calibration
;
; PURPOSE:
;
;       Return geometric parameters for meshes within the AIA entrance
;       and focal plane filters.
;
; CALLING SEQUENCE:
;
;       meshinfo = AIA_GETMESHINFO(wavelength, qabort=qabort)
;
; INPUTS:
;
;       WAVELENGTH - [Mandatory] (string scalar)
;                    Name of AIA UV/EUV passband for which the filter mesh
;                    info will be returned.
;                    Options: '94','131','171','193','211','304','335'
;
; KEYWORDS:
;
;       /USE_PREFLIGHTCORE - [Optional input] (Boolean scalar)
;                    Set this flag to select a model of the PSF core.
;                    0 (or default): Use a PSF core that is within the error
;                       budget limits and produces nice results for
;                       deconvolved images. This value has been optimized by
;                       visual evaluations.
;                    1 (or nonzero): Use a PSF core that is modeled directly
;                       on the component budget values, including the
;                       optical prescriptions.
;                    See Note #2 for further info.
;
;       QABORT     - [Optional output] (Boolean scalar)
;                    Status flag indicating whether the program ran to completion.
;                    0: Program ran to completion.
;                    1: Program aborted before completion.
;       DWAVELENGTH - dwavelength is the fractional change in the wavelength over the nominal wavelength
;                      so for 131, 131.5 dwavelength would be 0.5/131 = 0.00381679
             
;
; OUTPUTS:
;
;       Return - [Mandatory] (structure), with the following tags:
;         CHANNEL:    channel name
;         IMAGE:      image that as used to perform the measurements
;         REFIMAGE:   background image
;         SOLARNORTH: up or down, in the image
;         ANGLE1:     first angle
;         DELTA1:     error in first angle
;         ANGLE2:     etc.
;         DELTA2:
;         ANGLE3:
;         DELTA3:
;         ANGLE4:
;         DELTA4:
;         SPACING:    distance between diffraction spikes (from entrance filter)
;         DSPACING:   delta spacing
;         MESHPITCH:  pitch of the mesh in micrometers
;         MESHWIDTH:  width of the mesh in micrometers
;         FP_SPACING: distance between diffraction spikes (from focal plane filters)
;         GS_WIDTH:  the width applied to the Gaussian such that *after*
;           convolution we have the proper width
;
; COMMON BLOCKS:
;
;       none
;
; NOTES:
;
;       1) See AIA PSF documentation for detailed description of the
;          parameters and calculations performed herein.
;
;       2) Here is a sample of the pre-flight PSF core error budget, by
;          component, for the 171 A channel. (Excerpted from "Initial
;          Calibration of the Atmospheric Imaging Assembly (AIA) on the
;          Solar Dynamics Observatory (SDO)" by Boerner et al., and also
;          appearing in the AIA PSF documentation.)
;
;             ---------------------------------------------------------
;                                                      Contribution to
;                                                     RMS Spot diameter
;             Item                                        (arcsec)
;             ---------------------------------------------------------
;             Optical Prescription                          0.60
;             Fabrication, alignment and assembly effects   1.21
;             Launch shift effects                          0.10
;             On-orbit thermal effects                      0.21
;             Focus error (with on-orbit correction)        0.10
;             Jitter residual                               0.48
;             Detector Pixelization                         0.48
;             CCD Charge spreading                          0.80
;             ---------------------------------------------------------
;             On-orbit performance prediction               1.73
;             ---------------------------------------------------------
;
;          Setting /USE_PREFLIGHTCORE (in this example) will
;          use a Gaussian PSF core with a sigma of 1.73/2 arcsec.
;          This value will vary across the channels due to differences
;          in the optical prescriptions. NOT using this keyword
;          (i.e., the default functionality) will use some value
;          which is: (a) smaller than this; (b) constant across
;          channels; and (c) less pessimistic in the following sense.
;
;          The values represented in the table are conservative
;          estimates of the imaging performance based on pre-flight
;          measurements. They were originally derived in order to
;          verify that the instrument met its performance requirements,
;          and thus should be viewed as upper limits on the spot
;          diameter rather than as most-accurate predictions. For example, it
;          is known from the GT-to-science telescope alignment
;          shifts that the telescope mirrors did not move as much
;          as was allocated, and thus the contribution to the PSF
;          from optical misalignment is overstated. The AIA Team
;          is currently (January 2012) investigating the PSF core
;          empirically, using on-orbit data. This keyword's DEFAULT
;          mode (i.e., USE_PREFLIGHTCORE=0) will eventually be
;          updated with the results of that investigation.
;
; CONTACT:
;
;       Comments, feedback, and bug reports regarding this routine may be
;       directed to this email address:
;                boerner ~at~ lmsal.com
;
; MODIFICATION HISTORY:
;
progver = 'v2010-Jun-16' ;--- (P.Grigis (SAO)) Written.
progver = 'v2011-Jul-13' ;--- (P.Grigis (SAO)) Updated fp_spacing fields.
progver = 'v2011-Aug-15' ;--- (Y.Su (SAO)) Added gs_width information.
progver = 'v2012-Jan-03' ;--- (M.Weber (SAO)) Added abort case for wavelength.
;                             Also changed meshangle to meshinfo for consistency.
progver = 'v2012-Jan-06' ;--- (M.Weber (SAO)) Added use_preflightcore keyword,
;                             functionality, and documentation.
;
;-
; =========================================================================

  qpfcore = keyword_set(use_preflightcore)
  qabort = 0B

  ;; This defines the "wx" values for each channel, as will be used later in
  ;; the <aia_diffractionpattern> routine. The set used depends upon which
  ;; PSF core model is chosen with the /use_preflightcore keyword.
  wavenames = ['94','131','171','193','211','304','335']
  widths_default = fltarr(n_elements(wavenames)) + 4.5
  widths_preflight = [0.951, 1.033, 0.962, 1.512, 1.199, 1.247, 0.962]

  meshinfo     = {channel:'', $
                  image:'',$
                  refimage:'', $
                  solarnorth:'', $
                  angle1:0.0, $
                  delta1:0.0, $
                  angle2:0.0, $
                  delta2:0.0, $
                  angle3:0.0, $
                  delta3:0.0, $
                  angle4:0.0, $
                  delta4:0.0, $
                  spacing:0.0,$
                  dspacing:0.0, $
                  meshpitch:0.0, $
                  meshwidth:0.0, $
                  fp_spacing:0.0, $
                  gs_width:0.0}


  CASE wavelength OF

     '211': BEGIN
        ss_wave = where(wavenames EQ wavelength)
        meshinfo.channel='211'
        meshinfo.image='AIA20101016_191038_0211.fits'
        meshinfo.refimage='AIA20101016_190902_0211.fits'
        meshinfo.angle1=49.78
        meshinfo.delta1=0.02
        meshinfo.angle2=40.08
        meshinfo.delta2=0.02
        meshinfo.angle3=-40.34
        meshinfo.delta3=0.02
        meshinfo.angle4=-49.95
        meshinfo.delta4=0.02
        meshinfo.spacing=19.97
        meshinfo.dspacing=0.09
        meshinfo.solarnorth='UP'
        meshinfo.meshpitch=363.0
        meshinfo.meshwidth=34.0
        meshinfo.fp_spacing=0.465
        if (qpfcore EQ 1B) then meshinfo.gs_width=widths_preflight[ss_wave] $
                           else meshinfo.gs_width=widths_default[ss_wave]
     END

     '171': BEGIN
        ss_wave = where(wavenames EQ wavelength)
        meshinfo.channel='171'
        meshinfo.image='AIA20101016_191037_0171.fits'
        meshinfo.refimage='AIA20101016_190901_0171.fits'
        meshinfo.angle1=49.81
        meshinfo.delta1=0.02
        meshinfo.angle2=39.57
        meshinfo.delta2=0.02
        meshinfo.angle3=-40.13
        meshinfo.delta3=0.02
        meshinfo.angle4=-50.38
        meshinfo.delta4=0.02
        meshinfo.spacing=16.26
        meshinfo.dspacing=0.10
        meshinfo.solarnorth='UP'
        meshinfo.meshpitch=363.0
        meshinfo.meshwidth=34.0
        meshinfo.fp_spacing=0.377
        if (qpfcore EQ 1B) then meshinfo.gs_width=widths_preflight[ss_wave] $
                           else meshinfo.gs_width=widths_default[ss_wave]
     END

     '193' : BEGIN
        ss_wave = where(wavenames EQ wavelength)
        meshinfo.channel='193'
        meshinfo.image='AIA20101016_191056_0193.fits'
        meshinfo.refimage='AIA20101016_190844_0193.fits'
        meshinfo.angle1=49.82
        meshinfo.delta1=0.02
        meshinfo.angle2=39.57
        meshinfo.delta2=0.02
        meshinfo.angle3=-40.12
        meshinfo.delta3=0.03
        meshinfo.angle4=-50.37
        meshinfo.delta4=0.04
        meshinfo.spacing=18.39
        meshinfo.dspacing=0.20
        meshinfo.solarnorth='UP'
        meshinfo.meshpitch=363.0
        meshinfo.meshwidth=34.0
        meshinfo.fp_spacing=0.425
        if (qpfcore EQ 1B) then meshinfo.gs_width=widths_preflight[ss_wave] $
                           else meshinfo.gs_width=widths_default[ss_wave]
     END

     '335' : BEGIN
        ss_wave = where(wavenames EQ wavelength)
        meshinfo.channel='335'
        meshinfo.image='AIA20101016_191041_0335.fits'
        meshinfo.refimage='AIA20101016_190905_0335.fits'
        meshinfo.angle1=50.40
        meshinfo.delta1=0.02
        meshinfo.angle2=39.80
        meshinfo.delta2=0.02
        meshinfo.angle3=-39.64
        meshinfo.delta3=0.02
        meshinfo.angle4=-50.25
        meshinfo.delta4=0.02
        meshinfo.spacing=31.83
        meshinfo.dspacing=0.07
        meshinfo.solarnorth='UP'
        meshinfo.meshpitch=363.0
        meshinfo.meshwidth=34.0
        meshinfo.fp_spacing=0.738
        if (qpfcore EQ 1B) then meshinfo.gs_width=widths_preflight[ss_wave] $
                           else meshinfo.gs_width=widths_default[ss_wave]
     END

     '304' : BEGIN
        ss_wave = where(wavenames EQ wavelength)
        meshinfo.channel='304'
        meshinfo.image='AIA20101016_191021_0304.fits'
        meshinfo.refimage='AIA20101016_190845_0304.fits'
        meshinfo.angle1=49.76
        meshinfo.delta1=0.02
        meshinfo.angle2=40.18
        meshinfo.delta2=0.02
        meshinfo.angle3=-40.14
        meshinfo.delta3=0.02
        meshinfo.angle4=-49.90
        meshinfo.delta4=0.02
        meshinfo.spacing=28.87
        meshinfo.dspacing=0.05
        meshinfo.solarnorth='UP'
        meshinfo.meshpitch=363.0
        meshinfo.meshwidth=34.0
        meshinfo.fp_spacing=0.670
        if (qpfcore EQ 1B) then meshinfo.gs_width=widths_preflight[ss_wave] $
                           else meshinfo.gs_width=widths_default[ss_wave]
     END

     '131' : BEGIN
        ss_wave = where(wavenames EQ wavelength)
        meshinfo.channel='131'
        meshinfo.image='AIA20101016_191035_0131.fits'
        meshinfo.refimage='AIA20101016_190911_0131.fits'
        meshinfo.angle1=50.27
        meshinfo.delta1=0.02
        meshinfo.angle2=40.17
        meshinfo.delta2=0.02
        meshinfo.angle3=-39.70
        meshinfo.delta3=0.02
        meshinfo.angle4=-49.95
        meshinfo.delta4=0.02
        meshinfo.spacing=12.37
        meshinfo.dspacing=0.16
        meshinfo.solarnorth='UP'
        meshinfo.meshpitch=363.0
        meshinfo.meshwidth=34.0
        meshinfo.fp_spacing=0.289
        if (qpfcore EQ 1B) then meshinfo.gs_width=widths_preflight[ss_wave] $
                           else meshinfo.gs_width=widths_default[ss_wave]
     END

     '94' : BEGIN
        ss_wave = where(wavenames EQ wavelength)
        meshinfo.channel='94'
        meshinfo.image='AIA20101016_191039_0094.fits'
        meshinfo.refimage='AIA20101016_190903_0094.fits'
        meshinfo.angle1=49.81
        meshinfo.delta1=0.02
        meshinfo.angle2=40.16
        meshinfo.delta2=0.02
        meshinfo.angle3=-40.28
        meshinfo.delta3=0.02
        meshinfo.angle4=-49.92
        meshinfo.delta4=0.02
        meshinfo.spacing=8.99
        meshinfo.dspacing=0.13
        meshinfo.solarnorth='UP'
        meshinfo.meshpitch=363.0
        meshinfo.meshwidth=34.0
        meshinfo.fp_spacing=0.207
        if (qpfcore EQ 1B) then meshinfo.gs_width=widths_preflight[ss_wave] $
                           else meshinfo.gs_width=widths_default[ss_wave]
     END

     ELSE: BEGIN
        print, 'AIA_CALC_PSF_MOD: Input WAVELENGTH not a recognized value. ' $
        + 'See header documentation. Aborting.'
        qabort = 1B
        return, 0B
     END

  ENDCASE
  ;dwavelength is the fractional change in the wavelength over the nominal wavelength
  ;so for 131, 131.5 dwavelength would be 0.5/131 = 0.00381679

  If keyword_set( dwavelength ) then begin
     meshinfo.spacing *= (1.0 + dwavelength)
     print, meshinfo.spacing
     endif
  RETURN, meshinfo

END


; =========================================================================
; =========================================================================

FUNCTION aia_psfsinc, x

; =========================================================================
;+
; PROJECT:
;
;       SDO / AIA
;
; NAME:
;
;       AIA_PSFSINC
;
; CATEGORY:
;
;       Math utils
;
; PURPOSE:
;
;       Returns the sinc function y=sin(pi*x)/(pi*x).
;       Note that x=0 returns 1 (not NAN).
;
; CALLING SEQUENCE:
;
;       y = AIA_PSFSINC(x)
;
; INPUTS:
;
;       X - [Mandatory] (float or double precision real scalar)
;
; KEYWORDS:
;
;       n/a
;
; OUTPUTS:
;
;       Return - [Mandatory] (float or double precision real scalar)
;                The trigonometric sinc of input X:
;                y=sin(pi*x)/(pi*x) .
;
; COMMON BLOCKS:
;
;       none
;
; NOTES:
;
;       1) Y=1 (not NAN) when X=0.
;
; CONTACT:
;
;       Comments, feedback, and bug reports regarding this routine may be
;       directed to this email address:
;                boerner ~at~ lmsal.com
;
; MODIFICATION HISTORY:
;
progver = 'v2009-Sep-17' ;--- (P.Grigis (SAO)) Written.
;
;-
; =========================================================================

  ind=where( finite(x) EQ 1 AND x EQ x*0,count)

  arg=!Pi*x

  y=sin(arg)/arg
  IF count GT 0 THEN y[ind]=1.0

  RETURN,y


END


; =========================================================================
; =========================================================================

FUNCTION aia_diffractionpattern, meshinfo, $
                                 PSFEntranceFilter=PSFEntranceFilter, $
                                 PSFFocalPlane=PSFFocalPlane, $
                                 npix=npix


; =========================================================================
;+
; PROJECT:
;
;       SDO / AIA
;
; NAME:
;
;       AIA_DIFFRACTIONPATTERN
;
; CATEGORY:
;
;       Calibration
;
; PURPOSE:
;
;       Given the filter mesh geometric parameters for an AIA passband,
;       calculate the PSF for the diffraction and core effects, as a 2D array.
;
; CALLING SEQUENCE:
;
;       psf = AIA_DIFFRACTIONPATTERN(meshinfo)
;
; INPUTS:
;
;       MESHINFO - [Mandatory] (structure scalar)
;                  This is a structure that contains the geometric parameters
;                  for the filter meshes for a particular AIA passband. See the
;                  header documentation of AIA_GETMESHINFO for the definitions
;                  of the required tags.
;
; KEYWORDS:
;
;       PSFENTRANCEFILTER - [Optional output] (double float, 2D array, [4096,4096])
;                           A partial PSF from a midpoint in the calculations.
;                           Intended only for testing and debugging.
;
;       PSFFOCALPLANE     - [Optional output] (double float, 2D array, [4096,4096])
;                           A partial PSF from a midpoint in the calculations.
;                           Intended only for testing and debugging.
;
;
; OUTPUTS:
;
;       Return - [Mandatory] (double float, 2D array, [4096,4096])
;
; COMMON BLOCKS:
;
;       none
;
; NOTES:
;
;       1) See AIA PSF documentation for detailed description of the
;          parameters and calculations performed herein.
;
;       2) These PSFs account for the following components:
;          - the entrance filter mesh, producing a diffraction pattern;
;          - the mirror optics, producing a Gaussian spread to the PSF core;
;          - the focal plane filter mesh, producing a diffraction pattern; and
;          - the CCD charge-spreading, producing a Gaussian spread to the PSF core.
;
;       3) As of December 2011, the PSFs do **NOT** account for any
;          wide-angle spread by the optics. No analysis to-date has been able
;          to measure the PSF wings due to optical scatter, hence it is treated
;          as an insignificant effect.
;
;       4) Warning: This program can take several minutes or longer to run on
;          a workstation.
;
; CONTACT:
;
;       Comments, feedback, and bug reports regarding this routine may be
;       directed to this email address:
;                boerner ~at~ lmsal.com
;
; MODIFICATION HISTORY:
;
progver = 'v2011-Jun-16' ;--- (P.Grigis (SAO)) Written.
progver = 'v2012-Jan-03' ;--- (M.Weber (SAO)) Changed MESHINFO from keyword to
;                             parameter input. Modified formatting for consistency.
;
;-
; =========================================================================

  psf=dblarr(npix,npix)

  ;; This is the width applied to the Gaussian such that *after* convolution I
  ;; have the proper width (which is 4/3 at 1/e of max).
  wx=meshinfo.gs_width
  wy=meshinfo.gs_width

  print,'wx=',wx,'wy=',wy

  x=findgen(npix)+0.5
  y=x

  xx=x#(x*0+1)
  yy=(y*0+1)#y

  d=meshinfo.spacing
  dx1=meshinfo.spacing*cos(meshinfo.angle1/180*!pi)
  dy1=meshinfo.spacing*sin(meshinfo.angle1/180*!pi)

  dx2=meshinfo.spacing*cos(meshinfo.angle2/180*!pi)
  dy2=meshinfo.spacing*sin(meshinfo.angle2/180*!pi)

  dx3=meshinfo.spacing*cos(meshinfo.angle3/180*!pi)
  dy3=meshinfo.spacing*sin(meshinfo.angle3/180*!pi)

  dx4=meshinfo.spacing*cos(meshinfo.angle4/180*!pi)
  dy4=meshinfo.spacing*sin(meshinfo.angle4/180*!pi)

  meshratio=meshinfo.meshpitch/meshinfo.meshwidth
  k=1.0/(meshratio*d)

  dpx=0.5
  dpy=0.5

;==================================

  ;; This section computes the effect from the entrance filter mesh diffraction.
  FOR j=-100,100 DO BEGIN

     IF j NE 0 THEN BEGIN

        print,'First pass wide angle',j

        intensity=(aia_psfsinc(j*d*k))^2

        xc=npix/2.0+dx1*j+dpx
        yc=npix/2.0+dy1*j+dpy
        psf+=exp(-wx*(xx-xc)^2-wy*(yy-yc)^2)*intensity

        xc=npix/2.0+dx2*j+dpx
        yc=npix/2.0+dy2*j+dpy
        psf+=exp(-wx*(xx-xc)^2-wy*(yy-yc)^2)*intensity

        xc=npix/2.0+dx3*j+dpx
        yc=npix/2.0+dy3*j+dpy
        psf+=exp(-wx*(xx-xc)^2-wy*(yy-yc)^2)*intensity

        xc=npix/2.0+dx4*j+dpx
        yc=npix/2.0+dy4*j+dpy
        psf+=exp(-wx*(xx-xc)^2-wy*(yy-yc)^2)*intensity

     ENDIF

  ENDFOR

;==================================

  ;; This section applies modifications to the core of the entrance filter PSF
  ;; and normalizes.
  psf=psf/total(psf)*0.18
  psf2=exp(-wx*(xx-npix/2.0-dpx)^2-wy*(yy-npix/2.0-dpy)^2)
  psf2=psf2/total(psf2)*0.82
  PSFEntranceFilter=psf2+psf

;==================================

  ;; This section computes the effect from the focal plane filter mesh diffraction.
  psf=dblarr(npix,npix)

  d=meshinfo.fp_spacing

  meshratio=meshinfo.meshpitch/meshinfo.meshwidth
  k=1.0/(meshratio*d)


  dx1=d*cos(45.0/180*!pi)
  dy1=d*sin(45.0/180*!pi)

  dx2=d*cos(-45.0/180*!pi)
  dy2=d*sin(-45.0/180*!pi)

  dpx=0.5
  dpy=0.5

  FOR j=1,100 DO BEGIN

        print,j

        intensity=(aia_psfsinc(j*d*k))^2


        xc=npix/2.0+dx1*j+dpx
        yc=npix/2.0+dy1*j+dpy
        psf+=exp(-wx*(xx-xc)^2-wy*(yy-yc)^2)*intensity

        xc=npix/2.0+dx2*j+dpx
        yc=npix/2.0+dy2*j+dpy
        psf+=exp(-wx*(xx-xc)^2-wy*(yy-yc)^2)*intensity

        xc=npix/2.0-dx1*j+dpx
        yc=npix/2.0-dy1*j+dpy
        psf+=exp(-wx*(xx-xc)^2-wy*(yy-yc)^2)*intensity


        xc=npix/2.0-dx2*j+dpx
        yc=npix/2.0-dy2*j+dpy
        psf+=exp(-wx*(xx-xc)^2-wy*(yy-yc)^2)*intensity

  ENDFOR

;==================================

  ;; This section applies modifications to the core of the focal plane filter PSF
  ;; and normalizes.
  psf=psf/total(psf)*0.18
  psf2=exp(-wx*(xx-npix/2.0-dpx)^2-wy*(yy-npix/2.0-dpy)^2)
  psf2=psf2/total(psf2)*0.82
  PSFFocalPlane=psf2+psf

;==================================

  ;; This section generates the composite PSF and normalizes.

  psfnew=shift(abs(fft(fft(PSFFocalPlane)*fft(PSFEntranceFilter),-1)),npix/2,npix/2)

  psfnew=psfnew/total(psfnew)



  return, psfnew

END


; =========================================================================
; =========================================================================


FUNCTION aia_calc_psf_mod,	wavelength, $
                       		use_preflightcore=use_preflightcore, $
                       		qabort=qabort, npix=npix, $
                       		dwavelength = dwavelength

; =========================================================================
;+
; PROJECT:
;
;       SDO / AIA
;
; NAME:
;
;       AIA_CALC_PSF_MOD
;
; CATEGORY:
;
;       Calibration
;
; PURPOSE:
;
;       Return AIA PSF for a given passband, as a 2D array.
;
; CALLING SEQUENCE:
;
;       psf = AIA_CALC_PSF_MOD(wavelength)
;
; INPUTS:
;
;       /USE_PREFLIGHTCORE - [Optional input] (Boolean scalar)
;                    Set this flag to select a model of the PSF core.
;                    0 (or default): Use a PSF core that is within the error
;                       budget limits and produces nice results for
;                       deconvolved images. This value has been optimized by
;                       visual evaluations.
;                    1 (or nonzero): Use a PSF core that is modeled directly
;                       on the component budget values, including the
;                       optical prescriptions.
;                    See Note #2 in the header doc of <aia_getmeshinfo>
;                    for further info.
;
;       WAVELENGTH - [Mandatory] (string scalar)
;                    Name of AIA UV/EUV passband for which a PSF will be returned.
;                    Options: '94','131','171','193','211','304','335'
;
; KEYWORDS:
;
;       QABORT     - [Optional output] (Boolean scalar)
;                    Status flag indicating whether the program ran to completion.
;                    0: Program ran to completion.
;                    1: Program aborted before completion.
;       DWAVELENGTH - dwavelength is the fractional change in the wavelength over the nominal wavelength
;                      so for 131, 131.5 dwavelength would be 0.5/131 = 0.00381679
;       NPIX       - Number of pixels in psf, default is 4096L
;
;
; OUTPUTS:
;
;       Return - [Mandatory] (double float, 2D array, [npix,npix])
;
; COMMON BLOCKS:
;
; 	none
;
; NOTES:
;
;       1) See AIA PSF documentation for detailed description of the
;          parameters and calculations performed herein.
;
;       2) These PSFs account for the following components:
;          - the entrance filter mesh, producing a diffraction pattern;
;          - the mirror optics, producing a Gaussian spread to the PSF core;
;          - the focal plane filter mesh, producing a diffraction pattern; and
;          - the CCD charge-spreading, producing a Gaussian spread to the PSF core.
;
;       3) As of December 2011, the PSFs do **NOT** account for any
;          wide-angle spread by the optics. No analysis to-date has been able
;          to measure the PSF wings due to optical scatter, hence it is treated
;          as an insignificant effect.
;
; CONTACT:
;
;       Comments, feedback, and bug reports regarding this routine may be
;       directed to this email address:
;                boerner ~at~ lmsal.com
;
; MODIFICATION HISTORY:
;
progver = 'v2012-Jan-03' ;--- (M.Weber (SAO)) Written.
progver = 'v2012-Jan-06' ;--- (M.Weber (SAO)) Added use_preflightcore keyword,
;                             functionality, and documentation.
; added npix to allow for other than full 4096 psf
;
;-
; =========================================================================

  qabort = 0B
  default, npix, 4096L
  ;; Return geometric parameters for meshes within the AIA entrance
  ;; and focal plane filters.
  meshinfo = aia_getmeshinfo(wavelength, use_preflightcore=use_preflightcore, $
                             qabort=qabort, dwavelength = dwavelength )

  ;; Check whether aia_getmeshinfo subroutine has made a call to abort.
  if (qabort EQ 1B) then return, 0

  psf = aia_diffractionpattern(meshinfo, npix=npix)

  return, psf

END
